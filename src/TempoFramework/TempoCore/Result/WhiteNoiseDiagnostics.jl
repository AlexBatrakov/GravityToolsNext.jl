# src/TempoFramework/TempoCore/Result/WhiteNoiseDiagnostics.jl
# -----------------------------------------------------------------------------
# Diagnostics for WhiteNoise fit: weighted moments, χ² helpers, and a rich
# printable report. This file is self-contained (no dependency on WhiteNoise.jl)
# so it can be included in any order.
#
# Parent module should already have:
#   using Printf       # for @sprintf
#   using Statistics   # (optional; we don't call stats from here)
# -----------------------------------------------------------------------------

# =============================================================================
# diagnostics utils (internal) — only used by the report/printing
# =============================================================================

"""
    weighted_mean(res) -> (μ, se_μ)
    weighted_mean(res, unc) -> (μ_w, se_μ)

Unweighted fast path avoids allocating `ones(N)`. Weighted uses w = 1/unc^2.
Standard errors follow the same conventions as in the legacy printout.
"""
function weighted_mean(res::Vector{Float64})
    N = length(res)
    N == 0 && return (NaN, NaN)
    μ  = sum(res) / N
    se = 1 / sqrt(N)  # legacy convention for diagnostics
    return μ, se
end

function weighted_mean(res::Vector{Float64}, unc::Vector{Float64})
    @assert length(res) == length(unc)
    w  = @. inv(unc^2)
    sw = sum(w)
    μw = sum(w .* res) / sw
    se = sqrt(1 / sw)
    return μw, se
end

"""
    weighted_std(res) / weighted_std(res, unc) -> (σ, se_σ)

Centered at (weighted) mean. Unweighted uses population σ (divide by N),
to match the diagnostic behavior you had.
"""
function weighted_std(res::Vector{Float64})
    N = length(res)
    N == 0 && return (NaN, NaN)
    μ  = sum(res) / N
    σ  = sqrt(sum((res .- μ).^2) / N)
    N_eff = N
    se = σ / sqrt(2 * max(N_eff - 1, 1))
    return σ, se
end

function weighted_std(res::Vector{Float64}, unc::Vector{Float64})
    @assert length(res) == length(unc)
    w  = @. inv(unc^2)
    sw = sum(w)
    μw = sum(w .* res) / sw
    σ  = sqrt(sum(w .* (res .- μw).^2) / sw)
    N_eff = (sw^2) / sum(w.^2)
    se = σ / sqrt(2 * max(N_eff - 1, 1))
    return σ, se
end

"""
    weighted_skewness(res [, unc]) -> (g1, se_g1)
"""
function weighted_skewness(res::Vector{Float64})
    N = length(res)
    N == 0 && return (NaN, NaN)
    μ  = sum(res) / N
    σ  = sqrt(sum((res .- μ).^2) / N)
    num = sum((res .- μ).^3)
    den = N * σ^3
    g1  = num / den
    se  = sqrt(6.0 / N)
    return g1, se
end

function weighted_skewness(res::Vector{Float64}, unc::Vector{Float64})
    @assert length(res) == length(unc)
    w  = @. inv(unc^2)
    sw = sum(w)
    μw = sum(w .* res) / sw
    σ  = sqrt(sum(w .* (res .- μw).^2) / sw)
    num = sum(w .* (res .- μw).^3)
    den = sw * σ^3
    g1  = num / den
    N_eff = (sw^2) / sum(w.^2)
    se  = sqrt(6.0 / max(N_eff, 1))
    return g1, se
end

"""
    weighted_kurtosis(res [, unc]) -> (excess_g2, se_g2)
"""
function weighted_kurtosis(res::Vector{Float64})
    N = length(res)
    N == 0 && return (NaN, NaN)
    μ  = sum(res) / N
    σ  = sqrt(sum((res .- μ).^2) / N)
    num = sum((res .- μ).^4)
    den = N * σ^4
    g2  = num / den - 3.0
    se  = sqrt(24.0 / N)
    return g2, se
end

function weighted_kurtosis(res::Vector{Float64}, unc::Vector{Float64})
    @assert length(res) == length(unc)
    w  = @. inv(unc^2)
    sw = sum(w)
    μw = sum(w .* res) / sw
    σ  = sqrt(sum(w .* (res .- μw).^2) / sw)
    num = sum(w .* (res .- μw).^4)
    den = sw * σ^4
    g2  = num / den - 3.0
    N_eff = (sw^2) / sum(w.^2)
    se  = sqrt(24.0 / max(N_eff, 1))
    return g2, se
end

"""
    weighted_rms(res [, unc]) -> (rms, se_rms)
"""
function weighted_rms(res::Vector{Float64})
    N = length(res)
    N == 0 && return (NaN, NaN)
    rms = sqrt(sum(res.^2) / N)
    se  = rms / sqrt(2 * N)
    return rms, se
end

function weighted_rms(res::Vector{Float64}, unc::Vector{Float64})
    @assert length(res) == length(unc)
    w  = @. inv(unc^2)
    sw = sum(w)
    rms = sqrt(sum(w .* res.^2) / sw)
    N_eff = (sw^2) / sum(w.^2)
    se  = rms / sqrt(2 * max(N_eff, 1))
    return rms, se
end

"""
    chisq_stats(res_norm; dof=length(res_norm)) -> (χ², χ²_red)

For already normalized residuals (σ'=1 scale), χ² = ∑ r².
"""
function chisq_stats(res_norm::Vector{Float64}; dof::Int = length(res_norm))
    χ2 = sum(abs2, res_norm)
    return χ2, χ2 / dof
end

"""
    chisq_stats(res, unc; dof=length(res)) -> (χ², χ²_red)

For raw residuals with transformed uncertainties `unc`.
"""
function chisq_stats(res::Vector{Float64}, unc::Vector{Float64}; dof::Int = length(res))
    @assert length(res) == length(unc)
    χ2 = sum(abs2, res ./ unc)
    return χ2, χ2 / dof
end

# Local tiny helper to avoid include order issues with WhiteNoise.jl
# σ′ = sqrt((efac*σ)^2 + equad^2)
@inline function _diag_transform_uncertainties!(
    σ_tr::Vector{Float64}, σ_orig::Vector{Float64}, efac::Float64, equad::Float64
)
    @inbounds @simd for i in eachindex(σ_tr, σ_orig)
        σ_tr[i] = sqrt((efac * σ_orig[i])^2 + equad^2)
    end
    return σ_tr
end

# =============================================================================
# diagnostics — formatted console report (DRY, with optional expectations)
# =============================================================================

# Low-level: one-line formatter with optional expected value -> z-score
# Accepts either (val, err) tuple from weighted_* helpers or (val, err) scalars
function _fmt_stat(io::IO, label::AbstractString, val::Float64, err::Float64; expect::Union{Nothing,Float64}=nothing)
    if isfinite(err) && err > 0 && expect !== nothing
        z = (val - (expect::Float64)) / err
        println(io, @sprintf("    %-16s %12.6f ± %-10.6f   (z = %+6.2fσ)", label, val, err, z))
    else
        println(io, @sprintf("    %-16s %12.6f ± %-10.6f", label, val, err))
    end
end
_fmt_stat(io::IO, label::AbstractString, tup::Tuple{Float64,Float64}; expect=nothing) =
    _fmt_stat(io, label, tup[1], tup[2]; expect=expect)

# Print a four-line stats block (mean/std/skewness/kurtosis) for data `x`.
# If `w` is provided, uses weighted* with those weights; otherwise unit-weights.
# `normalized=true` adds expectations (0,1,0,0) and z-scores.
# `weighted_label=true` prefixes labels with "weighted ".
function _print_stats_block(
    io::IO,
    title::AbstractString,
    x::Vector{Float64};
    w::Union{Nothing,Vector{Float64}} = nothing,
    normalized::Bool = false,
    weighted_label::Bool = true,
)
    println(io, title)

    μ   = w === nothing ? weighted_mean(x)       : weighted_mean(x, w)
    s   = w === nothing ? weighted_std(x)        : weighted_std(x, w)
    g1  = w === nothing ? weighted_skewness(x)   : weighted_skewness(x, w)
    g2e = w === nothing ? weighted_kurtosis(x)   : weighted_kurtosis(x, w) # excess kurtosis

    exp_mean = normalized ? 0.0 : nothing
    exp_std  = normalized ? 1.0 : nothing
    exp_skew = normalized ? 0.0 : nothing
    exp_kurt = normalized ? 0.0 : nothing

    pref = (weighted_label ? "weighted " : "")

    _fmt_stat(io, pref * "mean:",     μ;  expect=exp_mean)
    _fmt_stat(io, pref * "std:",      s;  expect=exp_std)
    _fmt_stat(io, pref * "skewness:", g1; expect=exp_skew)
    _fmt_stat(io, pref * "kurtosis:", g2e; expect=exp_kurt)
end

"""
    print_white_noise_fit_report!(
        io::IO,
        backend::AbstractString,
        residuals::Vector{Float64},
        uncertainties_orig::Vector{Float64},
        uncertainties_tr::Vector{Float64},
        residuals_norm::Vector{Float64},
        residuals_shifted_norm::Vector{Float64},
        efac::Float64, equad::Float64, offset::Float64, ad_objective::Float64;
        n_toas::Int = length(residuals),
        sections::AbstractSet{Symbol} = Set([
            :header, :residuals, :uncertainties_orig, :uncertainties_tr, :normalized, :shifted_normalized, :chi2, :footer
        ]),
    )

Mutating workspace variant:
- Overwrites `uncertainties_tr`, `residuals_norm`, `residuals_shifted_norm`.
- Computes only what is required by `sections`.
"""
function print_white_noise_fit_report!(
    io::IO,
    backend::AbstractString,
    residuals::Vector{Float64},
    uncertainties_orig::Vector{Float64},
    uncertainties_tr::Vector{Float64},
    residuals_norm::Vector{Float64},
    residuals_shifted_norm::Vector{Float64},
    efac::Float64, equad::Float64, offset::Float64, ad_objective::Float64;
    n_toas::Int = length(residuals),
    sections::AbstractSet{Symbol} = Set([
        :header, :residuals, :uncertainties_orig, :uncertainties_tr, :normalized, :shifted_normalized, :chi2, :footer
    ]),
)
    has(s) = in(s, sections)

    # sanity
    n = length(residuals)
    @assert length(uncertainties_orig) == n
    @assert length(uncertainties_tr) == n
    @assert length(residuals_norm) == n
    @assert length(residuals_shifted_norm) == n

    # compute σ′ only if any section needs it
    need_tr = has(:residuals) || has(:uncertainties_tr) || has(:normalized) || has(:shifted_normalized) || has(:chi2)
    if need_tr
        _diag_transform_uncertainties!(uncertainties_tr, uncertainties_orig, efac, equad)
    end

    # precompute normalized variants only if needed
    if has(:normalized) || has(:shifted_normalized)
        @. residuals_norm = residuals / uncertainties_tr
        if has(:shifted_normalized)
            @. residuals_shifted_norm = (residuals - offset) / uncertainties_tr
        end
    end

    if has(:header)
        println(io, "backend: $backend")
        println(io, "    N_TOAs: $n_toas")
    end

    if has(:residuals)
        _print_stats_block(io, "residuals stats:", residuals; w=uncertainties_tr, normalized=false, weighted_label=true)
    end

    if has(:uncertainties_orig)
        _print_stats_block(io, "original uncertainties stats:", uncertainties_orig; w=nothing, normalized=false, weighted_label=false)
    end

    if has(:uncertainties_tr)
        _print_stats_block(io, "transformed uncertainties stats:", uncertainties_tr; w=nothing, normalized=false, weighted_label=false)
    end

    if has(:normalized)
        _print_stats_block(io, "normalized residuals stats:", residuals_norm; w=nothing, normalized=true, weighted_label=false)
    end

    if has(:shifted_normalized)
        _print_stats_block(io, "normalized (with offset) residuals stats:", residuals_shifted_norm; w=nothing, normalized=true, weighted_label=false)
    end

    if has(:chi2)
        println(io, "chi squared stats:")
        χ2, χ2r = chisq_stats(residuals, uncertainties_tr)
        println(io, "    chi2: $χ2")
        println(io, "    chi2r: $χ2r")
    end

    if has(:footer)
        println(io, "ad_objective_best:  $ad_objective")
        println(io, "efac_best:   $efac")
        println(io, "equad_best:  $equad")
        println(io, "offset_best: $offset")
        println(io)
    end
end

"""
    print_white_noise_fit_report(
        io::IO,
        backend, residuals, uncertainties_orig,
        efac, equad, offset, ad_objective;
        sections=Set([...])
    )

Non-mutating convenience:
allocates local work buffers and calls `print_white_noise_fit_report!`.
"""
function print_white_noise_fit_report(
    io::IO,
    backend::AbstractString,
    residuals::Vector{Float64},
    uncertainties_orig::Vector{Float64},
    efac::Float64, equad::Float64, offset::Float64, ad_objective::Float64;
    sections::AbstractSet{Symbol} = Set([
        :header, :residuals, :uncertainties_orig, :uncertainties_tr, :normalized, :shifted_normalized, :chi2, :footer
    ]),
)
    n = length(residuals)
    uncertainties_tr       = similar(uncertainties_orig)
    residuals_norm         = similar(residuals)
    residuals_shifted_norm = similar(residuals)

    print_white_noise_fit_report!(
        io, backend,
        residuals, uncertainties_orig,
        uncertainties_tr, residuals_norm, residuals_shifted_norm,
        efac, equad, offset, ad_objective;
        n_toas  = n,
        sections = sections,
    )
end

# Non-mutating convenience: allocate work buffers and print to STDOUT
function print_white_noise_fit_report(
    backend::Union{AbstractString,Symbol},
    residuals::Vector{Float64},
    uncertainties_orig::Vector{Float64},
    efac::Float64, equad::Float64, offset::Float64, ad_objective::Float64;
    sections::AbstractSet{Symbol} = Set([
        :header, :residuals, :uncertainties_orig, :uncertainties_tr,
        :normalized, :shifted_normalized, :chi2, :footer
    ]),
)
    n = length(residuals)
    uncertainties_tr       = similar(uncertainties_orig)
    residuals_norm         = similar(residuals)
    residuals_shifted_norm = similar(residuals)

    print_white_noise_fit_report!(
        stdout, String(backend),
        residuals, uncertainties_orig,
        uncertainties_tr, residuals_norm, residuals_shifted_norm,
        efac, equad, offset, ad_objective;
        n_toas  = n,
        sections = sections,
    )
end