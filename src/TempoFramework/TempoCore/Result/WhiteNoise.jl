# src/TempoFramework/TempoCore/Result/WhiteNoise.jl
# -----------------------------------------------------------------------------
# White-noise estimation and per-backend/global aggregation.
#
# Dependencies imported by the parent module:
#   using Printf
#   using Statistics
#   using Optim
#   using Roots
#   using HypothesisTests
#   using Distributions
#
# Expected from ResidualStats.jl:
#   - NormalizedResidualStats
#   - build_normalized_residual_statistics(x::Vector{Float64})
#   - _group_by_backend(entries; in_fit::Union{Nothing,Bool}=nothing)
#   - _empty_norm_stats()
# -----------------------------------------------------------------------------
# Overview
# This module estimates white-noise calibration parameters (EFAC, EQUAD, constant
# offset) per backend by minimizing the Anderson–Darling A² of normalized residuals,
# and aggregates results into a global summary. The implementation is conservative
# and robust by design:
#   - Per-backend estimation runs independently under `try/catch`; failures are
#     recorded and excluded from the global summary.
#   - Working buffers are reused and many routines mutate their first argument
#     to avoid allocations on hot paths.
#
# Units
#   - Residual-related quantities **and** TOA uncertainties are expected in
#     microseconds (µs). Normalized residuals are unitless.
#   - EFAC is dimensionless; EQUAD has units of µs; offset is in µs.
#
# Error policy
#   - The low-level solver `estimate_white_noise_ad_with_offset` may throw on
#     severely degenerate inputs (e.g., zero variance, infeasible initialization).
#     The top-level `build_white_noise_fit` guards against this and marks such
#     backends as `success=false`.

# =============================================================================
# Grid defaults for initialization (tunable constants; no logic changes)
# These values were chosen empirically to provide stable seeds without slowing the solve
# (11×11 grid over EQUAD and target std of normalized residuals in [0.95, 1.05]).
# =============================================================================

const WHITENOISE_GRID_EQUADS       = 11
const WHITENOISE_GRID_STD_TARGETS  = 11
const WHITENOISE_STD_TARGET_RANGE  = (0.95, 1.05)

# =============================================================================
# utils — low-level helpers 
# =============================================================================

"""
    transform_uncertainties!(uncertainties_tr, uncertainties_orig, efac, equad)
        -> uncertainties_tr

In-place transformation: σ' = sqrt((efac * σ)^2 + equad^2). Mutates the first argument.
"""
function transform_uncertainties!(
    uncertainties_tr::Vector{Float64},
    uncertainties_orig::Vector{Float64},
    efac::Float64,
    equad::Float64,
)
    @. uncertainties_tr = sqrt(efac^2 * uncertainties_orig^2 + equad^2)
    return uncertainties_tr
end

"Out-of-place variant of `transform_uncertainties!`."
function transform_uncertainties(uncertainties_orig::Vector{Float64}, efac::Float64, equad::Float64)
    uncertainties_tr = similar(uncertainties_orig)
    transform_uncertainties!(uncertainties_tr, uncertainties_orig, efac, equad)
    return uncertainties_tr
end

"""
    weighted_mean_invvar(residuals, uncertainties_tr) -> (μ_w, se_μ)

Inverse-variance weighted mean with w = 1/σ′^2 *without allocations*.
Skips entries with non-finite residuals/σ′ or σ′ ≤ 0. Returns (NaN, NaN) if no valid points.
"""
function weighted_mean_invvar(residuals::Vector{Float64},
                              uncertainties_tr::Vector{Float64})
    @assert length(residuals) == length(uncertainties_tr)
    sw  = 0.0
    swr = 0.0
    @inbounds for i in eachindex(residuals, uncertainties_tr)
        r = residuals[i]
        s = uncertainties_tr[i]
        if isfinite(r) && isfinite(s) && s > 0.0
            w = 1.0 / (s * s)          # 1/σ′^2
            sw  += w
            swr += w * r
        end
    end
    if sw == 0.0
        return NaN, NaN
    else
        μ  = swr / sw
        se = sqrt(1.0 / sw)
        return μ, se
    end
end

# =============================================================================
# objective — AD-based objectives and normalization with offset
# =============================================================================

"AD objective from precomputed normalized residuals."
function ad_objective_function(residuals_norm::Vector{Float64})
    ad_test = OneSampleADTest(residuals_norm, Normal(0, 1))
    return ad_test.A²
end

"""
    ad_objective_function!(residuals_norm, residuals, uncertainties_tr) -> A²

Fill `residuals_norm .= residuals ./ uncertainties_tr` and return Anderson–Darling A².
Mutates `residuals_norm` as a working buffer.
"""
function ad_objective_function!(
    residuals_norm::Vector{Float64},
    residuals::Vector{Float64},
    uncertainties_tr::Vector{Float64},
)
    @. residuals_norm = residuals / uncertainties_tr
    return ad_objective_function(residuals_norm)
end


"""
    ad_objective_function_with_offset!(residuals_shifted_norm, residuals, uncertainties_tr, offset) -> A²

Fill `residuals_shifted_norm .= (residuals .- offset) ./ uncertainties_tr` and return A².
Mutates `residuals_shifted_norm` as a working buffer.
"""
function ad_objective_function_with_offset!(
    residuals_shifted_norm::Vector{Float64},
    residuals::Vector{Float64},
    uncertainties_tr::Vector{Float64},
    offset::Float64,
)
    @. residuals_shifted_norm = (residuals - offset) / uncertainties_tr
    return ad_objective_function(residuals_shifted_norm)
end

"""
    ad_objective_function_fit_offset!(residuals_shifted_norm, residuals, uncertainties_tr)
        -> (A²_min, offset*)

1D minimize A² over `offset`, given fixed `uncertainties_tr`. Uses `Optim.optimize`.
Mutates `residuals_shifted_norm` as a working buffer.
"""
function ad_objective_function_fit_offset!(
    residuals_shifted_norm::Vector{Float64},
    residuals::Vector{Float64},
    uncertainties_tr::Vector{Float64},
)
    ad_objective_scalar = offset ->
        ad_objective_function_with_offset!(residuals_shifted_norm, residuals, uncertainties_tr, offset)

    ad_objective_vector = args -> (@inbounds ad_objective_scalar(args[1]))

    offset0, _ = weighted_mean_invvar(residuals, uncertainties_tr)
    isfinite(offset0) || (offset0 = 0.0)

    optim_result = optimize(ad_objective_vector, [offset0], NelderMead())
    offset = optim_result.minimizer[1]
    ad_objective_val = ad_objective_scalar(offset)
    return ad_objective_val, offset
end

# =============================================================================
# solver — EFAC/EQUAD/offset search (initialization + optimization)
# =============================================================================

"""
    find_efac_for_fixed_equad_and_std!(residuals_norm, uncertainties_tr, residuals, uncertainties_orig,
                                       equad_fixed, std_res_norm_target) -> efac

Solve for `efac` such that std(residuals ./ σ′(efac, equad_fixed)) ≈ std_res_norm_target.

Notes:
- Monotone in `efac` (for fixed `equad`), so we first try to bracket a root and use `Bisection()`.
- Mutates `residuals_norm` and `uncertainties_tr` as work buffers; no allocations.
- Returns `NaN` if the target is infeasible or a robust root was not found.
"""
function find_efac_for_fixed_equad_and_std!(
    residuals_norm::Vector{Float64},
    uncertainties_tr::Vector{Float64},
    residuals::Vector{Float64},
    uncertainties_orig::Vector{Float64},
    equad_fixed::Float64,
    std_res_norm_target::Float64,
)
    # --- sanity checks
    n = length(residuals)
    @assert length(uncertainties_orig) == n
    @assert length(uncertainties_tr)   == n
    @assert length(residuals_norm)     == n

    # quick rejects / degenerate cases
    std_res = std(residuals)
    if !(isfinite(std_res) && std_res > 0.0) || !(isfinite(std_res_norm_target) && std_res_norm_target > 0.0)
        return NaN
    end

    # If EQUAD alone already enforces a std below the target, the solution lies at efac≈0
    if equad_fixed >= std_res / std_res_norm_target
        return NaN  # keep behavior: signal to skip this grid point
    end

    # Objective: f(efac) = target/std(residuals_norm(efac)) - 1
    # Monotone in efac ⇒ single root and easy bracketing.
    objective = function (efac::Float64)
        if !(isfinite(efac)) || efac < 0.0
            return NaN
        end
        transform_uncertainties!(uncertainties_tr, uncertainties_orig, efac, equad_fixed)
        @. residuals_norm = residuals / uncertainties_tr
        s = std(residuals_norm)
        return (isfinite(s) && s > 0.0) ? (std_res_norm_target / s - 1.0) : NaN
    end

    # Seed from homoscedastic approximation (matches the legacy initializer)
    mean_unc_sq = sum(abs2, uncertainties_orig) / n
    if !(isfinite(mean_unc_sq) && mean_unc_sq > 0.0)
        return NaN
    end
    efac0 = sqrt(((std_res / std_res_norm_target)^2 - equad_fixed^2) / mean_unc_sq)
    efac0 = isfinite(efac0) ? max(efac0, eps()) : 1.0

    # Try to bracket a sign change around efac0 and solve with Bisection()
    lo = efac0 * 0.5
    hi = efac0 * 2.0
    flo = objective(lo)
    fhi = objective(hi)

    it = 0
    while !(isfinite(flo) && isfinite(fhi) && flo * fhi < 0.0) && it < 30
        lo = max(lo * 0.5, eps())
        hi = hi * 2.0
        flo = objective(lo)
        fhi = objective(hi)
        it += 1
    end

    efac = NaN
    try
        if isfinite(flo) && isfinite(fhi) && flo * fhi < 0.0
            efac = find_zero(objective, (lo, hi), Bisection())
        else
            efac = find_zero(objective, efac0)
        end
    catch
        return NaN
    end

    return abs(efac)
end

"""
    estimate_white_noise_ad_with_offset(residuals, uncertainties_orig; kwargs...)
        -> (efac, equad, offset, ad_objective, converged)

Jointly estimate (efac, equad, offset) by minimizing AD A² of normalized residuals.
- If no initial guesses are provided, builds a small grid over `equad` and target `std(res_norm)`,
  solves `efac` for each pair, fits `offset`, and seeds `Optim.optimize`.
- Mutates working buffers internally for performance (documented below).
- `converged` is the boolean flag returned by `Optim.converged(...)` for the outer 3D solve.
Notes
- Inputs are expected in µs for residuals and TOA uncertainties; the objective
  uses unitless normalized residuals.
- On severely degenerate inputs (e.g., zero variance or infeasible seeds), the
  optimizer may throw; callers are expected to guard with `try/catch` (as done in
  `build_white_noise_fit`) and mark the backend as failed.
"""
function estimate_white_noise_ad_with_offset(
    residuals::Vector{Float64},
    uncertainties_orig::Vector{Float64};
    print_results::Bool = false,
    plot_results::Bool = false,  # plotting removed
    backend::String = "",
    efac_init_in::Union{Float64, Nothing} = nothing,
    equad_init_in::Union{Float64, Nothing} = nothing,
    offset_init_in::Union{Float64, Nothing} = nothing,
)
    @assert (!isempty(residuals) || !isempty(uncertainties_orig))  "empty input"
    @assert all(isfinite, residuals) "residuals contain non-finite values"
    @assert all(u -> isfinite(u) && u > 0.0, uncertainties_orig) "uncertainties_orig must be positive and finite"
    @assert length(residuals) == length(uncertainties_orig) "residuals and uncertainties_orig must match"

    uncertainties_tr       = similar(uncertainties_orig)  # σ' workspace
    residuals_norm         = similar(residuals)           # res/σ' workspace
    residuals_shifted_norm = similar(residuals)           # (res-offset)/σ' workspace

    equad_max = std(residuals)

    # Objective over [efac, equad, offset]; mutates `uncertainties_tr` & `residuals_shifted_norm`.
    function make_objective_function(
        residuals_shifted_norm::Vector{Float64},
        uncertainties_tr::Vector{Float64},
        residuals::Vector{Float64},
        uncertainties_orig::Vector{Float64},
    )
        objective_fun = (args) -> begin
            efac, equad, offset = args
            transform_uncertainties!(uncertainties_tr, uncertainties_orig, efac, equad)
            return ad_objective_function_with_offset!(
                residuals_shifted_norm, residuals, uncertainties_tr, offset
            )
        end
        return objective_fun
    end

    efac_equad_offset_ad_objective = make_objective_function(
        residuals_shifted_norm, uncertainties_tr, residuals, uncertainties_orig
    )

    if isnothing(efac_init_in) || isnothing(equad_init_in) || isnothing(offset_init_in)
        N_equads        = WHITENOISE_GRID_EQUADS
        N_std_res_norm  = WHITENOISE_GRID_STD_TARGETS
        std_lo, std_hi  = WHITENOISE_STD_TARGET_RANGE

        equad_init_arr        = collect(LinRange(0.0, equad_max, N_equads))
        std_res_norm_init_arr = collect(LinRange(std_lo, std_hi, N_std_res_norm))

        efac_init_best = equad_init_best = offset_init_best = NaN
        ad_objective_init_best = Inf

        for equad_init in equad_init_arr, std_res_norm_init in std_res_norm_init_arr
            efac_init = find_efac_for_fixed_equad_and_std!(
                residuals_norm, uncertainties_tr, residuals, uncertainties_orig,
                equad_init, std_res_norm_init
            )
            if isnan(efac_init); continue; end

            transform_uncertainties!(uncertainties_tr, uncertainties_orig, efac_init, equad_init)
            ad_objective_init, offset_init = ad_objective_function_fit_offset!(
                residuals_shifted_norm, residuals, uncertainties_tr
            )

            if ad_objective_init < ad_objective_init_best
                ad_objective_init_best = ad_objective_init
                efac_init_best   = efac_init
                equad_init_best  = equad_init
                offset_init_best = offset_init
            end
        end

        ad_result = optimize(efac_equad_offset_ad_objective, [efac_init_best, equad_init_best, offset_init_best], NelderMead())
    else
        ad_result = optimize(efac_equad_offset_ad_objective, [efac_init_in, equad_init_in, offset_init_in], NelderMead())
    end

    efac_best, equad_best, offset_best = Optim.minimizer(ad_result)
    efac_best  = max(abs(efac_best), eps())
    equad_best = abs(equad_best)
    ad_objective_best = ad_result.minimum
    converged_flag = Optim.converged(ad_result)

    if print_results
        print_white_noise_fit_report(
            String(backend),
            residuals, uncertainties_orig,
            efac_best, equad_best, offset_best, ad_objective_best;
            sections = Set([:header, :residuals, :uncertainties_orig, :uncertainties_tr, :normalized, :shifted_normalized, :chi2, :footer]),
        )
    end

    return efac_best, equad_best, offset_best, ad_objective_best, converged_flag
end

# =============================================================================
# aggregation — per-backend run + global summary (from TempoResult.jl)
# =============================================================================

"""
    WhiteNoiseBackendFitResult

Per-backend white-noise fit result.

Fields
- `backend`       : backend id as `Symbol`
- `efac`,`equad`  : white-noise parameters used to scale uncertainties
- `offset`        : constant offset applied before normalization
- `ad_objective`  : value of the optimized objective (e.g., AD A²)
- `stats`         : normalized residual stats with fitted efac/equad/offset
- `success`       : parameters are finite/usable; diagnostics/stats computed
- `converged`     : optimizer reported convergence (useful even if `success=true`)

Notes
- When `success == false`, `stats` contains a neutral placeholder (from `_empty_norm_stats()`),
  and the pretty-printer skips printing the stats block for this backend.
"""
struct WhiteNoiseBackendFitResult
    backend::Symbol
    efac::Float64
    equad::Float64
    offset::Float64
    ad_objective::Float64
    stats::NormalizedResidualStats
    success::Bool       # parameters are finite/usable; diagnostics/stats computed
    converged::Bool     # Optimizer reported convergence (useful even if success=true)
end

function Base.show(io::IO, ::MIME"text/plain", r::WhiteNoiseBackendFitResult)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    iop    = IOContext(io, :indent => indent + 4)

    println(io, pad, "WhiteNoiseBackendFitResult[$(r.backend)]")
    println(io, spad, "efac         = ", @sprintf("%.6g", r.efac))
    println(io, spad, "equad        = ", @sprintf("%.6g", r.equad))
    println(io, spad, "offset       = ", @sprintf("%.6g", r.offset))
    println(io, spad, "ad_objective = ", @sprintf("%.6g", r.ad_objective))
    println(io, spad, "success      = ", r.success)
    println(io, spad, "converged    = ", r.converged)
    println(io, spad, "stats:")
    if r.success
        show(iop, MIME"text/plain"(), r.stats)
    else
        println(iop, "(not available; success=false)")
    end
end


"""
    WhiteNoiseFitResult

White-noise fit for all backends + global summary over concatenated normalized residuals.

Fields
- `by_backend`      : map `backend => WhiteNoiseBackendFitResult`
- `global_stats`    : summary over *successful* backends only
- `failed_backends` : vector of backend ids that did not produce usable parameters (`success=false`)
"""
struct WhiteNoiseFitResult
    by_backend::Dict{Symbol, WhiteNoiseBackendFitResult}
    global_stats::NormalizedResidualStats
    failed_backends::Vector{Symbol}
end

function Base.show(io::IO, ::MIME"text/plain", r::WhiteNoiseFitResult)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    tpad   = repeat(" ", indent + 4)
    iop    = IOContext(io, :indent => indent + 4)
    ios    = IOContext(io, :indent => indent + 6)

    println(io, pad, "WhiteNoiseFitResult")
    if isempty(r.by_backend)
        println(io, spad, "by_backend: (none)")
    else
        println(io, spad, "by_backend:")
        println(io, tpad, "fit results (", length(r.by_backend), " backends):")
        print_aligned_table(
            ios, r.by_backend;
            namecol = :__key__,                     # Dict key = backend symbol
            cols    = (
                :efac, :equad, 
                "log_equad" => (x -> log10(x.equad) - 6.0),
                :offset, :ad_objective,
                "p_value" => (:stats, :ad_p_value),
                :success, :converged,
            ),
            sort_by = (:ad_objective,),             # sort by AD objective (ascending)
            order   = :asc,
            name_label = "backend",
            formats = Dict(
                :efac          => "%10.6f",
                :equad         => "%10.4f",
                "log_equad"    => "%10.6f",
                :offset        => "%10.6f",
                :ad_objective  => "%10.6f",
                "p_value"      => "%10.3e",
            ),
            renderers = Dict(
                :success   => x -> (x ? "✓" : "✗"),
                :converged => x -> (x ? "✓" : "✗"),
            ),
            aligns = Dict(
                :success   => :center,
                :converged => :center,
            ),
            header_align = :auto,
        )
        println(io, tpad, "fit stats (", length(r.by_backend), " backends):")
        print_aligned_table(
            ios, r.by_backend;
            namecol = :__key__,                     # Dict key = backend symbol
            cols    = (
                "n"          => (:stats, :n),
                "mean"       => (:stats, :mean),
                "median"     => (:stats, :median),
                "std"        => (:stats, :std),
                "kurtosis"   => (:stats, :kurtosis),
                "skewness"   => (:stats, :skewness),
                "chisqr"     => (:stats, :chisqr),
                "red_chisqr" => (:stats, :red_chisqr),
            ),
            sort_by = (:ad_objective,),             # sort by AD objective (ascending)
            order   = :asc,
            name_label = "backend",
            formats = Dict(
                "n"            => "%7d",
                "mean"         => "%10.6f",
                "median"       => "%10.6f",
                "std"          => "%10.6f",
                "kurtosis"     => "%10.6f",
                "skewness"     => "%10.6f",
                "chisqr"       => "%10.4f",
                "red_chisqr"   => "%10.6f",
            ),
            header_align = :auto,
        )
        if !isempty(r.failed_backends)
            println(io, tpad, "failed_backends (", length(r.failed_backends), "): ",
                    join(string.(r.failed_backends), ", "))
        end
    end
    println(io, spad, "global_stats:")
    show(iop, MIME"text/plain"(), r.global_stats)
end

"""
    _normalized_shifted_residuals_for_backend(residuals, uncertainties_orig, efac, equad, offset)
        -> Vector{Float64}

Transform uncertainties (σ' = sqrt((efac·σ)^2 + equad^2)) and compute shifted,
normalized residuals `(residuals - offset) ./ σ'` for a single backend.
Returns a vector aligned with the inputs (no filtering performed here).
Inputs are assumed to be in microseconds (µs); the result is unitless.
"""
function _normalized_shifted_residuals_for_backend(
    residuals::Vector{Float64},
    uncertainties_orig::Vector{Float64},
    efac::Float64,
    equad::Float64,
    offset::Float64,
)
    @assert length(residuals) == length(uncertainties_orig)

    uncertainties_tr = similar(uncertainties_orig)
    transform_uncertainties!(uncertainties_tr, uncertainties_orig, efac, equad)

    residuals_shifted_norm = similar(residuals)
    @. residuals_shifted_norm = (residuals - offset) / uncertainties_tr

    return residuals_shifted_norm
end

"""
    build_white_noise_fit(entries::Vector{CombinedTOAEntry}) -> WhiteNoiseFitResult

Per-backend estimation of (efac, equad, offset) on in-fit TOAs. Robust to empty groups and
solver failures. Produces per-backend normalized residual stats and a global summary.
Inputs are expected in microseconds (µs) for residuals and TOA uncertainties; normalized residuals are unitless.

Behavior
- If a backend fails (throws) or yields non-finite parameters, it is marked `success=false` and added
  to `failed_backends`. Its stats are set to a neutral placeholder and its contribution is excluded
  from `global_stats`.
- `global_stats` are computed by concatenating normalized residuals from successful backends only.
"""
function build_white_noise_fit(entries::Vector{CombinedTOAEntry})::WhiteNoiseFitResult
    results_by_backend            = Dict{Symbol, WhiteNoiseBackendFitResult}()
    residuals_shifted_norm_global = Float64[]
    failed_backends = Symbol[]

    entries_grouped = _group_by_backend(entries; in_fit=true)

    for (backend, es) in entries_grouped
        residuals          = Float64[e.residual_tn      for e in es]
        uncertainties_orig = Float64[e.uncertainty_orig for e in es]
        isempty(residuals) && continue

        efac = equad = offset = ad_obj = NaN
        success = true
        converged = false
        try
            efac, equad, offset, ad_obj, converged = estimate_white_noise_ad_with_offset(
                residuals, uncertainties_orig;
                backend = String(backend),
                print_results = false,
                plot_results  = false,
            )
        catch err
            @warn "White-noise fit failed" backend error=err
            success = false
            converged = false
        end

        if !isfinite(efac) || !isfinite(equad) || !isfinite(offset)
            @warn "White-noise fit returned non-finite parameters" backend efac=efac equad=equad offset=offset
            success = false
        end

        residuals_shifted_norm = success ? _normalized_shifted_residuals_for_backend(
            residuals, uncertainties_orig, efac, equad, offset
        ) : nothing

        stats = isnothing(residuals_shifted_norm) ? _empty_norm_stats() : build_normalized_residual_statistics(residuals_shifted_norm)

        results_by_backend[backend] = WhiteNoiseBackendFitResult(backend, efac, equad, offset, ad_obj, stats, success, converged)
        success || push!(failed_backends, backend)
        residuals_shifted_norm === nothing || append!(residuals_shifted_norm_global, residuals_shifted_norm)
    end

    global_stats = isempty(residuals_shifted_norm_global) ? _empty_norm_stats() : build_normalized_residual_statistics(residuals_shifted_norm_global)
    return WhiteNoiseFitResult(results_by_backend, global_stats, failed_backends)
end



