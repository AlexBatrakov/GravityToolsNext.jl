# src/TempoFramework/TempoCore/Result/GeneralResult.jl
# -----------------------------------------------------------------------------
# Top-level GeneralTempoResult container and constructors.
# -----------------------------------------------------------------------------
# defines:
#   ParamEstimate (alias)
#   GeneralTempoResult (struct)
#   GeneralTempoResult(iterations; ...) constructor
#   GeneralTempoResult(iter; ...) convenience constructor
# uses:
#   build_core_metrics, build_convergence_info

# -----------------------------------------------------------------------------
# Compact scalar metrics and helpers for GeneralTempoResult.
# -----------------------------------------------------------------------------
# defines:
#   build_core_metrics
#   result_metric   # safe lookup: get(res.metrics, key, NaN)
# uses:
#   InternalIterationResult, ConvergenceInfo

# --------------------------------------------------------------------------------------------------------------
# Final aggregated result of a TEMPO run / task
# --------------------------------------------------------------------------------------------------------------

# Small alias for parameter estimates we want to expose
const ParamEstimate = NamedTuple{(:value, :uncertainty), Tuple{Float64, Float64}}

"""
    build_core_metrics(final_it::InternalIterationResult, conv::ConvergenceInfo;
                       include_tim::Bool=true,
                       include_tim_even_if_same::Bool=false)
        -> Dict{Symbol,Float64}

Produce a compact, task-agnostic set of scalar metrics from the final iteration and
convergence info. Missing pieces are filled with `NaN`.

Included keys (always):
- :chi2_fit, :chi2r_fit, :wrms_fit, :wrms_tn_fit
- :pre_post_final
- :delta_wrms_tn, :delta_chi2
- :ad_white_fit   (global AD after white-noise fit; NaN if not available)

Additionally (only if `include_tim == true` and either the in-fit and in-tim sets differ
or `include_tim_even_if_same == true`):
- :chi2_tim, :chi2r_tim, :wrms_tim, :wrms_tn_tim
"""
function build_core_metrics(final_it::InternalIterationResult, conv::ConvergenceInfo;
                            include_tim::Bool=true,
                            include_tim_even_if_same::Bool=false)
    m = Dict{Symbol,Float64}()

    # Final TEMPO basic block
    basic = final_it.output.basic

    # Residual stats (may be missing)
    if final_it.stats !== nothing
        fit_all = final_it.stats.in_fit.all
        m[:wrms_fit]     = fit_all.raw.wrms
        m[:wrms_tn_fit]  = fit_all.tn.wrms
        m[:chi2_fit]     = fit_all.norm_global.chisqr
        m[:chi2r_fit]    = fit_all.norm_global.red_chisqr

        # Add _tim metrics only if requested and not a trivial alias (or forced)
        add_tim = include_tim && (final_it.stats.in_fit !== final_it.stats.in_tim || include_tim_even_if_same)
        if add_tim
            tim_all = final_it.stats.in_tim.all
            m[:wrms_tim]     = tim_all.raw.wrms
            m[:wrms_tn_tim]  = tim_all.tn.wrms
            m[:chi2_tim]     = tim_all.norm_global.chisqr
            m[:chi2r_tim]    = tim_all.norm_global.red_chisqr
        end
    else
        # No stats available
        m[:wrms_fit]     = NaN
        m[:wrms_tn_fit]  = NaN
        m[:chi2_fit]     = NaN
        m[:chi2r_fit]    = NaN
        # Note: _tim keys are omitted entirely when stats are missing
    end

    # Final pre/post from TEMPO (may be NaN)
    m[:pre_post_final] = (basic === nothing || !isfinite(basic.pre_post)) ? NaN : basic.pre_post

    # Convergence deltas between the last two points (NaN if <2 points)
    m[:delta_wrms_tn] = conv.wrms_tn.final_abs_delta
    m[:delta_chi2]    = conv.chisqr.final_abs_delta

    # Global AD after white-noise fit (if performed)
    if final_it.white_noise_fit !== nothing
        m[:ad_white_fit] = final_it.white_noise_fit.global_stats.ad_statistic
    else
        m[:ad_white_fit] = NaN
    end

    return m
end

"""
    GeneralTempoResult

Top-level container with everything you usually want after a run.

Fields
- `iterations`      : all internal iterations (may be empty if not saved)
- `final_index`     : index of the last successful iteration (falls back to last)
- `final`           : the final iteration result (by `final_index`)
- `convergence`     : convergence summary across iterations
- `par_file_final`  : written output par-file (if available)
- `param_estimates` : map `:NAME => (value, uncertainty)` extracted from the final iteration
- `residuals`       : residual statistics group for the final iteration (or `nothing`)
- `white_noise`     : per-backend white-noise fit for the final iteration (or `nothing`)
- `metrics`         : compact scalar metrics for quick ranking/comparison
- `subresults`      : optional nested results (e.g., per-epoch, per-band, grid cells)
- `subresult_type`  : tag describing what `subresults` mean (e.g., `:epoch`, `:band`, `:grid`)
- `metadata`        : extra info (paths, timings, seeds, etc.)
"""
struct GeneralTempoResult
    iterations::Vector{InternalIterationResult}
    final_index::Int
    final::InternalIterationResult

    convergence::ConvergenceInfo
    par_file_final::Union{TempoParFile, Nothing}

    param_estimates::Dict{Symbol, ParamEstimate}
    residuals::Union{ResidualStatisticsGroup, Nothing}
    white_noise::Union{WhiteNoiseFitResult, Nothing}

    metrics::Dict{Symbol, Float64}

    subresults::Vector{GeneralTempoResult}
    subresult_type::Union{Symbol, Nothing}

    metadata::Dict{Symbol, Any}
end

# ---------- pretty print ----------
function Base.show(io::IO, ::MIME"text/plain", r::GeneralTempoResult)
    indent = get(io, :indent, 0)
    pad  = repeat(" ", indent)
    spad = repeat(" ", indent + 2)
    iop  = IOContext(io, :indent => indent + 4)

    println(io, pad,  "GeneralTempoResult:")
    println(io, spad, "iterations:      ", length(r.iterations), " entries")
    println(io, spad, "final_index:     ", r.final_index)
    println(io, spad, "final:           present")
    println(io, spad, "convergence:     ", r.convergence.converged ? "yes" : "no")
    if r.par_file_final !== nothing
        println(io, spad, "par_file_final:  present")
    else
        println(io, spad, "par_file_final:  -")
    end
    println(io, spad, "param_estimates: ", length(r.param_estimates), " keys")
    println(io, spad, "residuals:       ", r.residuals === nothing ? "-" : "present")
    println(io, spad, "white_noise:     ", r.white_noise === nothing ? "-" : "present")
    println(io, spad, "metrics:         ", isempty(r.metrics) ? "-" : string(length(r.metrics), " keys"))
    if !isempty(r.subresults)
        println(io, spad, "subresults:     ", length(r.subresults), " (type=", r.subresult_type, ")")
    else
        println(io, spad, "subresults:      -")
    end
    isempty(r.metadata) || println(io, spad, "metadata keys:   ", isempty(r.metadata) ? "-" : string(length(r.metadata), " keys"))
end

# Base.show(io::IO, r::GeneralTempoResult) = show(io, MIME"text/plain"(), r)

# ---------- helpers ----------
# Extract (value, uncertainty) from the final fit parameters
function _extract_param_estimates(fit_params::Vector{FitParameter})
    out = Dict{Symbol, ParamEstimate}()
    @inbounds for p in fit_params
        if isfinite(p.post_fit) && isfinite(p.uncertainty)
            out[p.name_symbol] = (value = p.post_fit, uncertainty = p.uncertainty)
        elseif isfinite(p.post_fit) && !haskey(out, p.name_symbol)
            out[p.name_symbol] = (value = p.post_fit, uncertainty = NaN)
        end
    end
    return out
end

# Find last successful iteration (no TEMPO error); fallback to last index
function _last_successful_index(iters::Vector{InternalIterationResult})
    idx = findlast(it -> !iserror(it.output.error), iters)
    return idx === nothing ? length(iters) : idx
end

# ---------- main builder ----------
"""
    GeneralTempoResult(
        iterations;
        par_file_final=nothing,
        subresults=GeneralTempoResult[],
        subresult_type=nothing,
        metadata=Dict(),
        metrics_hook=nothing
    )

Build a `GeneralTempoResult` from a list of `InternalIterationResult`s.
- Picks the last successful iteration as `final`.
- Computes convergence across all iterations.
- Extracts parameter estimates from the final iteration.
- Carries over residual stats and white-noise fit from the final iteration.
- Computes `metrics` via `build_core_metrics(final, convergence)` and optionally merges
  a user-supplied `metrics_hook(final, convergence)` dictionary (values must be `Real`).
"""
function GeneralTempoResult(
    iterations::Vector{InternalIterationResult};
    par_file_final::Union{TempoParFile,Nothing}=nothing,
    subresults::Vector{GeneralTempoResult}=GeneralTempoResult[],
    subresult_type::Union{Symbol,Nothing}=nothing,
    metadata::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    metrics_hook::Union{Nothing,Function}=nothing,
)
    isempty(iterations) && error("GeneralTempoResult: iterations are empty")

    final_idx = _last_successful_index(iterations)
    final_it  = iterations[final_idx]

    conv   = build_convergence_info(iterations)
    params = _extract_param_estimates(final_it.output.fit_parameters)

    core = build_core_metrics(final_it, conv)

    if metrics_hook !== nothing
        extra_any = metrics_hook(final_it, conv)
        @assert extra_any isa AbstractDict "metrics_hook must return a Dict-like object"
        # coerce to Float64, allow Int/Real; drop non-finite to keep metrics clean
        extra = Dict{Symbol,Float64}()
        for (k, v) in extra_any
            ks = k isa Symbol ? k : Symbol(k)
            if v isa Real
                vf = Float64(v)
                if isfinite(vf)
                    extra[ks] = vf
                else
                    extra[ks] = NaN
                end
            end
        end
        merge!(core, extra)  # user values override core if keys collide
    end

    return GeneralTempoResult(
        iterations,
        final_idx,
        final_it,
        conv,
        par_file_final,
        params,
        final_it.stats,
        final_it.white_noise_fit,
        core,
        subresults,
        subresult_type,
        metadata,
    )
end

# Convenience constructor from a single iteration
function GeneralTempoResult(
    iter::InternalIterationResult;
    par_file_final::Union{TempoParFile,Nothing}=nothing,
    subresults::Vector{GeneralTempoResult}=GeneralTempoResult[],
    subresult_type::Union{Symbol,Nothing}=nothing,
    metadata::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    metrics_hook::Union{Nothing,Function}=nothing,
)
    return GeneralTempoResult([iter];
        par_file_final = par_file_final,
        subresults = subresults,
        subresult_type = subresult_type,
        metadata = metadata,
        metrics_hook = metrics_hook,
    )
end

"Safe lookup of a scalar metric from `GeneralTempoResult.metrics`."
result_metric(res::GeneralTempoResult, key::Symbol)::Float64 =
    get(res.metrics, key, NaN)