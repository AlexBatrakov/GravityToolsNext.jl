using Printf                    # @sprintf
using Statistics                # mean, median, std, quantile
using StatsBase: skewness, kurtosis
using HypothesisTests: OneSampleADTest, ApproximateOneSampleKSTest,
                       JarqueBeraTest, ShapiroWilkTest, pvalue
using Distributions: Normal

#--------------------------------------------------------------------------------------------------------------
# Basic (non-normalized) residual stats
#--------------------------------------------------------------------------------------------------------------

"""
    BasicResidualStats

Summary for raw (non-normalized) residuals.

Fields
- `n`       : number of points
- `min,max` : range
- `mean`    : arithmetic mean
- `wmean`   : weighted mean (if weights given; else = mean)
- `median`  : sample median
- `rms`     : √(mean(x^2))             — uncentered
- `wrms`    : √(∑w x^2 / ∑w)           — uncentered, weighted
- `std`     : √(mean((x-mean)^2))      — centered
- `wstd`    : √(∑w (x-wmean)^2 / ∑w)   — centered, weighted
- `mad`     : median(|x - median(x)|)  — robust spread
"""
struct BasicResidualStats
    n::Int
    min::Float64
    max::Float64
    mean::Float64
    wmean::Float64
    median::Float64
    rms::Float64
    wrms::Float64
    std::Float64
    wstd::Float64
    mad::Float64
end

function Base.show(io::IO, s::BasicResidualStats)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)

    println(io, pad,  "Basic residual stats:")
    println(io, spad, "n      = ", s.n)
    println(io, spad, "min    = ", @sprintf("%.5f", s.min),    "   max   = ", @sprintf("%.5f", s.max))
    println(io, spad, "mean   = ", @sprintf("%.5f", s.mean),   "   wmean = ", @sprintf("%.5f", s.wmean))
    println(io, spad, "median = ", @sprintf("%.5f", s.median), "   mad   = ", @sprintf("%.5f", s.mad))
    println(io, spad, "rms    = ", @sprintf("%.5f", s.rms),    "   wrms  = ", @sprintf("%.5f", s.wrms))
    println(io, spad, "std    = ", @sprintf("%.5f", s.std),    "   wstd  = ", @sprintf("%.5f", s.wstd))
end

# Unweighted
function build_basic_residual_statistics(x::AbstractVector{<:Real})::BasicResidualStats
    n = length(x)
    @assert n > 0 "x must be non-empty"

    sx = 0.0
    sx2 = 0.0
    xmin = Inf
    xmax = -Inf

    @inbounds for xi0 in x
        xi = float(xi0)
        sx  += xi
        sx2 += xi*xi
        if xi < xmin; xmin = xi; end
        if xi > xmax; xmax = xi; end
    end

    μ    = sx / n
    med  = Statistics.median(x)
    rms  = sqrt(sx2 / n)
    std  = sqrt(max(rms*rms - μ*μ, 0.0))   # centered via second moment
    mad  = Statistics.median(abs.(x .- med))

    return BasicResidualStats(n, xmin, xmax, μ, μ, med, rms, rms, std, std, mad)
end

# Weighted
function build_basic_residual_statistics(x::AbstractVector{<:Real},
                                         w::AbstractVector{<:Real})::BasicResidualStats
    @assert length(x) == length(w) "x and w must have the same length"
    n = length(x)
    @assert n > 0 "x must be non-empty"

    sx = 0.0
    sx2 = 0.0
    sw = 0.0
    swx = 0.0
    swx2 = 0.0
    xmin = Inf
    xmax = -Inf

    @inbounds for i in eachindex(x, w)
        xi = float(x[i])
        wi = float(w[i])
        @assert wi >= 0 "weights must be non-negative"
        sx   += xi
        sx2  += xi*xi
        sw   += wi
        swx  += wi*xi
        swx2 += wi*xi*xi
        if xi < xmin; xmin = xi; end
        if xi > xmax; xmax = xi; end
    end
    @assert sw > 0 "sum(weights) must be > 0"

    μ    = sx / n
    μw   = swx / sw
    med  = Statistics.median(x)
    rms  = sqrt(sx2 / n)
    wrms = sqrt(swx2 / sw)
    std  = sqrt(max(rms*rms - μ*μ, 0.0))
    wstd = sqrt(max(wrms*wrms - μw*μw, 0.0))
    mad  = Statistics.median(abs.(x .- med))

    return BasicResidualStats(n, xmin, xmax, μ, μw, med, rms, wrms, std, wstd, mad)
end

#--------------------------------------------------------------------------------------------------------------
# Normalized residuals stats (for residuals ideally ~ N(0,1))
#--------------------------------------------------------------------------------------------------------------

const _P68  = 0.6826894921370859
const _P95  = 0.9544997361036416
const _P997 = 0.9973002039367398

struct NormalizedResidualStats
    n::Int
    min::Float64
    max::Float64
    mean::Float64
    median::Float64
    std::Float64
    skewness::Float64
    kurtosis::Float64
    chisqr::Float64
    rchisqr::Float64
    q01::Float64
    q05::Float64
    q95::Float64
    q99::Float64
    within1σ::Int
    within2σ::Int
    within3σ::Int
    ad_statistic::Float64
    ad_p_value::Float64
    ks_statistic::Float64
    ks_p_value::Float64
    jb_statistic::Float64
    jb_p_value::Float64
    sw_statistic::Float64
    sw_p_value::Float64
end

# pretty-print helper
_show_stat(io, label, val, expect, sigma) = begin
    z = (val - expect) / sigma
    println(io, @sprintf("    %-12s = %10.5f   (%+6.2fσ)", label, val, z))
end

function Base.show(io::IO, s::NormalizedResidualStats)
    n = s.n
    println(io, "Normalized residual statistics:")
    println(io, "  Summary:")
    println(io, @sprintf("    %-12s = %10d", "count", n))
    _show_stat(io, "minimum",  s.min,    0.0, s.std)
    _show_stat(io, "maximum",  s.max,    0.0, s.std)
    _show_stat(io, "mean",     s.mean,   0.0, 1 / sqrt(max(n, 1)))
    _show_stat(io, "median",   s.median, 0.0, 1.253 / sqrt(max(n, 1)))
    _show_stat(io, "std",      s.std,    1.0, 1 / sqrt(2 * max(n - 1, 1)))  # ≈ se(std) for Normal
    _show_stat(io, "skewness", s.skewness, 0.0, sqrt(6 / max(n, 1)))
    _show_stat(io, "kurtosis", s.kurtosis, 0.0, sqrt(24 / max(n, 1)))
    println(io, @sprintf("    %-12s = %10.5f   (red = %.5f)", "chisqr", s.chisqr, s.rchisqr))

    println(io, "\n  Quantiles:")
    println(io, @sprintf("    q01=%.3f  q05=%.3f  q95=%.3f  q99=%.3f", s.q01, s.q05, s.q95, s.q99))

    println(io, "\n  Empirical coverage:")
    p1 = s.within1σ / max(n, 1)
    p2 = s.within2σ / max(n, 1)
    p3 = s.within3σ / max(n, 1)
    println(io, @sprintf("    |x|≤1σ : %5d / %d  (%.3f, expected %.3f)", s.within1σ, n, p1, _P68))
    println(io, @sprintf("    |x|≤2σ : %5d / %d  (%.3f, expected %.3f)", s.within2σ, n, p2, _P95))
    println(io, @sprintf("    |x|≤3σ : %5d / %d  (%.3f, expected %.3f)", s.within3σ, n, p3, _P997))

    println(io, "\n  Normality tests:")
    println(io, @sprintf("    %-22s = %10.5f   (p = %6.5f)", "Anderson-Darling  A²", s.ad_statistic, s.ad_p_value))
    println(io, @sprintf("    %-22s = %10.5f   (p = %6.5f)", "Kolmogorov-Smirnov D", s.ks_statistic, s.ks_p_value))
    println(io, @sprintf("    %-22s = %10.5f   (p = %6.5f)", "Jarque-Bera        JB", s.jb_statistic, s.jb_p_value))
    println(io, @sprintf("    %-22s = %10.5f   (p = %6.5f)", "Shapiro-Wilk       W",  s.sw_statistic, s.sw_p_value))
end

# Safe p-values/tests
_safepv(t) = try pvalue(t) catch; NaN end

function _safe_sw(x)
    try
        t = ShapiroWilkTest(x); return t.W, _safepv(t)
    catch
        return NaN, NaN
    end
end

function _safe_ad(x)
    try
        t = OneSampleADTest(x, Normal(0,1)); return t.A², _safepv(t)
    catch
        return NaN, NaN
    end
end

function _safe_ks(x)
    try
        t = ApproximateOneSampleKSTest(x, Normal(0,1)); return t.δ, _safepv(t)
    catch
        return NaN, NaN
    end
end

function _safe_jb(x)
    try
        t = JarqueBeraTest(x); return t.JB, _safepv(t)
    catch
        return NaN, NaN
    end
end

function build_normalized_residual_statistics(x::Vector{Float64})::NormalizedResidualStats
    # filter non-finite just in case
    xf = Vector{Float64}(undef, 0)
    sizehint!(xf, length(x))
    @inbounds for v in x
        isfinite(v) && push!(xf, v)
    end
    n = length(xf)
    @assert n > 0 "normalized residuals are empty"

    # summary
    mn   = minimum(xf)
    mx   = maximum(xf)
    μ    = mean(xf)
    med  = median(xf)
    s    = std(xf)
    sk   = skewness(xf)
    ku   = kurtosis(xf)
    χ2   = sum(v->v*v, xf)
    rχ2  = χ2 / n

    # quantiles
    q01 = quantile(xf, 0.01)
    q05 = quantile(xf, 0.05)
    q95 = quantile(xf, 0.95)
    q99 = quantile(xf, 0.99)

    # empirical coverage
    absx = abs.(xf)
    w1 = count(<=(1.0), absx)
    w2 = count(<=(2.0), absx)
    w3 = count(<=(3.0), absx)

    # tests
    a2, ap  = _safe_ad(xf)
    ks, kp  = _safe_ks(xf)
    jb, jbp = _safe_jb(xf)
    sw, swp = _safe_sw(xf)

    return NormalizedResidualStats(
        n, mn, mx, μ, med, s, sk, ku,
        χ2, rχ2,
        q01, q05, q95, q99,
        w1, w2, w3,
        a2, ap, ks, kp, jb, jbp, sw, swp
    )
end

#--------------------------------------------------------------------------------------------------------------

struct ResidualStatistics
    raw::BasicResidualStats
    tn::BasicResidualStats
    norm_global::NormalizedResidualStats
    norm_local::NormalizedResidualStats
end

function Base.show(io::IO, stats::ResidualStatistics)
    println(io, "Residual Statistics:\n")

    print(io, " → Raw (residual):\n")
    show(IOContext(io, :compact => true), stats.raw)
    println(io)

    print(io, "\n → Whitened (residual_tn):\n")
    show(IOContext(io, :compact => true), stats.tn)
    println(io)

    print(io, "\n → Globally centered normalized whitened ((residual_tn - global_wmean) / uncertainty):\n")
    show(IOContext(io, :compact => true), stats.norm_global)
    println(io)

    print(io, "\n → Locally centered normalized whitened ((residual_tn - local_wmean) / uncertainty):\n")
    show(IOContext(io, :compact => true), stats.norm_local)
end

function build_residual_statistics(
    entries::Vector{CombinedTOAEntry},
    global_wmean::Float64
)::ResidualStatistics
    n = length(entries)

    residuals     = Vector{Float64}(undef, n)
    residuals_tn  = Vector{Float64}(undef, n)
    uncertainties = Vector{Float64}(undef, n)
    weights       = Vector{Float64}(undef, n)

    @inbounds for i in 1:n
        e = entries[i]
        residuals[i]     = e.residual
        residuals_tn[i]  = e.residual_tn
        uncertainties[i] = e.uncertainty
        weights[i]       = e.weight
    end

    raw = build_basic_residual_statistics(residuals, weights)
    tn  = build_basic_residual_statistics(residuals_tn, weights)

    # normalization: skip zero/invalid uncertainties
    norm_global = Float64[]
    norm_local  = Float64[]
    local_wmean = tn.wmean

    @inbounds for i in 1:n
        u = uncertainties[i]
        if isfinite(u) && u > 0.0
            push!(norm_global, (residuals_tn[i] - global_wmean) / u)
            push!(norm_local,  (residuals_tn[i] - local_wmean)  / u)
        end
    end

    norm_global_stats = isempty(norm_global) ? _empty_norm_stats() : build_normalized_residual_statistics(norm_global)
    norm_local_stats  = isempty(norm_local)  ? _empty_norm_stats() : build_normalized_residual_statistics(norm_local)

    return ResidualStatistics(raw, tn, norm_global_stats, norm_local_stats)
end

#--------------------------------------------------------------------------------------------------------------

struct ResidualStatisticsEntry
    all::ResidualStatistics
    by_backend::Dict{Symbol, ResidualStatistics}
end

"""
    _group_by_backend(entries; in_fit=true) -> Dict{Symbol, Vector{CombinedTOAEntry}}

Group entries by backend (Symbol). Optionally filter by `in_fit`.
"""
# разрешаем: nothing (=не фильтровать), true, false
function _group_by_backend(entries::Vector{CombinedTOAEntry}; in_fit::Union{Nothing,Bool}=nothing)
    groups = Dict{Symbol, Vector{CombinedTOAEntry}}()
    @inbounds for e in entries
        if isnothing(in_fit) || e.in_fit == in_fit
            be = Symbol(e.backend)
            push!(get!(groups, be, CombinedTOAEntry[]), e)
        end
    end
    return groups
end

function build_residual_statistics_entry(
    entries::Vector{CombinedTOAEntry},
    global_wmean::Float64
)::ResidualStatisticsEntry
    all_stats = build_residual_statistics(entries, global_wmean)

    grouped = _group_by_backend(entries; in_fit=nothing)
    backend_stats = Dict{Symbol, ResidualStatistics}()
    for (backend, group_entries) in grouped
        backend_stats[backend] = build_residual_statistics(group_entries, global_wmean)
    end

    return ResidualStatisticsEntry(all_stats, backend_stats)
end

function Base.show(io::IO, entry::ResidualStatisticsEntry)
    indent = get(io, :indent, 0)
    pad = repeat(" ", indent)

    println(io, pad, "All TOAs:")
    show(IOContext(io, :indent => indent + 2), entry.all)
    println(io)

    if isempty(entry.by_backend)
        println(io, pad, "Backends: (none)")
        return
    end

    println(io, pad, "Backends:")
    for be in sort!(collect(keys(entry.by_backend)))
        stats = entry.by_backend[be]
        n     = stats.norm_global.n
        wrms  = stats.tn.wrms
        println(io, pad, "  [", String(be), "] n=", n, "  wrms_tn=", @sprintf("%.3f", wrms))
    end
end

#--------------------------------------------------------------------------------------------------------------

struct ResidualStatisticsGroup
    in_fit::ResidualStatisticsEntry
    in_tim::ResidualStatisticsEntry
end

function Base.show(io::IO, group::ResidualStatisticsGroup)
    println(io, "Residual Statistics Group:")

    n_fit = group.in_fit.all.norm_global.n
    println(io, "\n→ In-Fit (", n_fit, " TOAs):")
    show(IOContext(io, :indent => 4), group.in_fit)

    n_tim = group.in_tim.all.norm_global.n
    println(io, "\n→ In-TIM (", n_tim, " TOAs):")
    show(IOContext(io, :indent => 4), group.in_tim)
end

# robust weighted mean for global centering; ignores non-finite/invalid weights
function compute_global_wmean(entries::Vector{CombinedTOAEntry})
    isempty(entries) && return 0.0
    wsum  = 0.0
    wrsum = 0.0
    @inbounds for e in entries
        w = e.weight
        r = e.residual_tn
        if isfinite(w) && isfinite(r) && w > 0.0
            wsum  += w
            wrsum += w * r
        end
    end
    if wsum > 0.0
        return wrsum / wsum
    else
        vals = [e.residual_tn for e in entries if isfinite(e.residual_tn)]
        return isempty(vals) ? 0.0 : mean(vals)
    end
end

"""
    build_residual_statistics_group(entries) -> ResidualStatisticsGroup

Build statistics for in-fit and in-tim sets. If there is no time window
(i.e. all entries are in-fit), the same `ResidualStatisticsEntry` object
is reused (aliased) for both fields to avoid duplicate computation/memory.
You can detect this via `group.in_fit === group.in_tim`.
"""
function build_residual_statistics_group(entries::Vector{CombinedTOAEntry})::ResidualStatisticsGroup
    # If all entries are marked in_fit, avoid creating a filtered copy
    all_in_fit = all(e -> e.in_fit, entries)

    # Choose the vector to use for the "in-fit" group
    entries_in_fit = all_in_fit ? entries : [e for e in entries if e.in_fit]
    if isempty(entries_in_fit)
        # No in-fit points; fall back to all entries to keep pipeline alive
        entries_in_fit = entries
    end

    # Global centering is computed from the in-fit set
    global_wmean = compute_global_wmean(entries_in_fit)

    if all_in_fit
        # Compute once and alias for both fields
        entry = build_residual_statistics_entry(entries, global_wmean)
        return ResidualStatisticsGroup(entry, entry)
    else
        # Compute separately for in-fit vs in-tim (all)
        return ResidualStatisticsGroup(
            build_residual_statistics_entry(entries_in_fit, global_wmean),
            build_residual_statistics_entry(entries,        global_wmean),
        )
    end
end

"""
    in_fit_equals_in_tim(group::ResidualStatisticsGroup) -> Bool

Returns `true` if the group used a single aliased entry for both fields
(no time window), i.e. `group.in_fit === group.in_tim`.
"""
in_fit_equals_in_tim(group::ResidualStatisticsGroup) = (group.in_fit === group.in_tim)

# --------------------------------------------------------------------------------------------------------------
# White-noise fit (per backend) + global aggregation
# --------------------------------------------------------------------------------------------------------------

"""
    WhiteNoiseBackendFitResult

Per-backend white-noise fit result.

Fields
- `backend`       : backend id as `Symbol`
- `efac`, `equad` : white-noise parameters used to scale uncertainties
- `offset`        : constant offset applied before normalization
- `ad_objective`  : value of the optimized objective (e.g., AD A²)
- `stats`         : normalized residual stats with fitted efac/equad/offset
"""
struct WhiteNoiseBackendFitResult
    backend::Symbol
    efac::Float64
    equad::Float64
    offset::Float64
    ad_objective::Float64
    stats::NormalizedResidualStats
end

"""
    WhiteNoiseFitResult

White-noise fit for all backends + global summary over concatenated normalized residuals.
"""
struct WhiteNoiseFitResult
    by_backend::Dict{Symbol, WhiteNoiseBackendFitResult}
    global_stats::NormalizedResidualStats
end

# -- helpers --------------------------------------------------------------------

"""
    _normalized_residuals_for_backend(residuals, sigma_orig, efac, equad, offset)
        -> Union{Vector{Float64}, Nothing}

Safely transform uncertainties and compute normalized residuals.
Returns `nothing` if no valid points remain.
"""
function _normalized_residuals_for_backend(
    residuals::Vector{Float64},
    sigma_orig::Vector{Float64},
    efac::Float64,
    equad::Float64,
    offset::Float64,
)
    @assert length(residuals) == length(sigma_orig)
    n = length(residuals)

    σt = similar(sigma_orig)
    transform_uncertainties!(σt, sigma_orig, efac, equad) # σ' = sqrt((efac*σ)^2 + equad^2)

    norm = Float64[]
    sizehint!(norm, n)
    @inbounds for i in 1:n
        r = residuals[i]
        s = σt[i]
        if isfinite(r) && isfinite(s) && s > 0.0
            push!(norm, (r - offset) / s)
        end
    end
    return isempty(norm) ? nothing : norm
end

"""
    _empty_norm_stats() -> NormalizedResidualStats

Return a neutral/empty `NormalizedResidualStats` so downstream code doesn't crash on empty groups.
"""
function _empty_norm_stats()::NormalizedResidualStats
    return NormalizedResidualStats(
        0, NaN, NaN, NaN, NaN, NaN, NaN, NaN,     # n, min, max, mean, median, std, skew, kurt
        NaN, NaN,                                 # chisqr, rchisqr
        NaN, NaN, NaN, NaN,                       # q01, q05, q95, q99
        0, 0, 0,                                  # within1σ, within2σ, within3σ
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN    # AD, KS, JB, SW stats+pvals
    )
end

# -- main -----------------------------------------------------------------------

"""
    build_white_noise_fit(entries::Vector{CombinedTOAEntry}) -> WhiteNoiseFitResult

Per-backend estimation of (efac, equad, offset) on in-fit TOAs. Robust to empty groups and
fitter failures. Produces per-backend normalized residual stats and a global summary.
"""
function build_white_noise_fit(entries::Vector{CombinedTOAEntry})::WhiteNoiseFitResult
    by_backend  = Dict{Symbol, WhiteNoiseBackendFitResult}()
    global_norm = Float64[]

    grouped = _group_by_backend(entries; in_fit=true)

    for (backend, es) in grouped
        # collect raw vectors
        residuals  = Float64[e.residual_tn      for e in es]
        sigma_orig = Float64[e.uncertainty_orig for e in es]
        isempty(residuals) && continue

        # fit efac/equad/offset; guard exceptions
        efac = equad = offset = ad_obj = NaN
        ok = true
        try
            efac, equad, offset, ad_obj = estimate_WhiteNoise_AD_with_offset(
                residuals, sigma_orig;
                backend = String(backend),
                print_results = false,
                plot_results  = false,
            )
        catch err
            @warn "White-noise fit failed" backend error=err
            ok = false
        end

        # sane fallbacks
        if !ok || !isfinite(efac) || !isfinite(equad) || !isfinite(offset)
            efac   = isfinite(efac)   ? efac   : 1.0
            equad  = isfinite(equad)  ? equad  : 0.0
            offset = isfinite(offset) ? offset : 0.0
            ad_obj = isfinite(ad_obj) ? ad_obj : NaN
        end

        # normalized residuals (may be nothing if all invalid)
        norm = _normalized_residuals_for_backend(residuals, sigma_orig, efac, equad, offset)

        # per-backend stats
        stats = isnothing(norm) ? _empty_norm_stats() : build_normalized_residual_statistics(norm)

        by_backend[backend] = WhiteNoiseBackendFitResult(
            backend, efac, equad, offset, ad_obj, stats
        )

        if norm !== nothing
            append!(global_norm, norm)
        end
    end

    # global stats across all backends
    global_stats = isempty(global_norm) ? _empty_norm_stats() : build_normalized_residual_statistics(global_norm)

    return WhiteNoiseFitResult(by_backend, global_stats)
end

# -- optional pretty printers ---------------------------------------------------

function Base.show(io::IO, r::WhiteNoiseBackendFitResult)
    println(io, "WhiteNoiseBackendFitResult[$(r.backend)]")
    println(io, "  efac        = ", @sprintf("%.6g", r.efac))
    println(io, "  equad       = ", @sprintf("%.6g", r.equad))
    println(io, "  offset      = ", @sprintf("%.6g", r.offset))
    println(io, "  ad_objective= ", @sprintf("%.6g", r.ad_objective))
    println(io, "  stats:")
    show(IOContext(io, :indent => 4), r.stats)
end

function Base.show(io::IO, R::WhiteNoiseFitResult)
    println(io, "WhiteNoiseFitResult")
    if isempty(R.by_backend)
        println(io, "  by_backend: (none)")
    else
        println(io, "  by_backend:")
        for be in sort!(collect(keys(R.by_backend)))
            res = R.by_backend[be]
            println(io, "    • ", be,
                    "  efac=",   @sprintf("%.4f", res.efac),
                    "  equad=",  @sprintf("%.4g", res.equad),
                    "  offset=", @sprintf("%.4g", res.offset),
                    "  AD=",     @sprintf("%.4g", res.ad_objective))
        end
    end
    println(io, "\n  global_stats:")
    show(IOContext(io, :indent => 2), R.global_stats)
end


#--------------------------------------------------------------------------------------------------------------
struct InternalIterationResult
    output::InternalIterationOutput
    residuals::Union{Vector{TempoResidualEntry}, Nothing}
    stats::Union{ResidualStatisticsGroup, Nothing}
    white_noise_fit::Union{WhiteNoiseFitResult, Nothing}
    metadata::Dict{Symbol, Any}
end

function Base.show(io::IO, r::InternalIterationResult)
    println(io, "InternalIterationResult:")

    println(io, "  ▸ TEMPO output (basic):")
    show(IOContext(io, :indent => 4), r.output.basic)
    println(io)

    if iserror(r.output.error)
        println(io, "  ▸ TEMPO error:")
        show(IOContext(io, :indent => 4), r.output.error)
    else
        println(io, "  ▸ TEMPO error: none")
    end

    if isnothing(r.stats)
        println(io, "  ▸ Residual statistics: none")
    else
        println(io, "  ▸ Residual statistics (in-fit only):")
        show(IOContext(io, :indent => 4), r.stats.in_fit)
        println(io)
    end

    if isnothing(r.residuals)
        println(io, "  ▸ Residuals: not stored")
    else
        println(io, "  ▸ Residuals: ", length(r.residuals), " entries loaded")
    end

    if isnothing(r.white_noise_fit)
        println(io, "  ▸ White-noise fit: none")
    else
        println(io, "  ▸ White-noise fit:")
        show(IOContext(io, :indent => 4), r.white_noise_fit)
        println(io)
    end

    if isempty(r.metadata)
        println(io, "  ▸ Metadata: none")
    else
        println(io, "  ▸ Metadata keys: ", join(string.(keys(r.metadata)), ", "))
    end
end

"""
    build_internal_iteration_result(output, residual_path, tim_entries; kwargs...) -> InternalIterationResult

Build a single-iteration result from parsed TEMPO output and (optionally) residuals on disk.
- Reads residuals if `residual_path` is a valid file.
- Combines TIM entries with residuals to compute residual statistics.
- Optionally performs white-noise fitting (`analyze_white_noise=true`).
"""
function build_internal_iteration_result(
    output::InternalIterationOutput,
    residual_path::Union{Nothing, String},
    tim_entries::Vector{TimTOAEntry};
    save_residuals::Bool = false,
    analyze_white_noise::Bool = false,
    time_start::Union{Nothing, Float64} = nothing,
    time_finish::Union{Nothing, Float64} = nothing
)::InternalIterationResult

    # If TEMPO reported an error — return early with the error recorded
    if iserror(output.error)
        return InternalIterationResult(
            output,
            nothing,
            nothing,
            nothing,  # white_noise_fit
            Dict(:error => output.error),
        )
    end

    # Read residuals if a file path is provided and exists
    residuals::Union{Nothing, Vector{TempoResidualEntry}} = nothing
    if residual_path !== nothing && isfile(residual_path)
        try
            residuals = read_residual_file(residual_path)
        catch err
            @warn "Failed to read residual file" path=residual_path error=err
            residuals = nothing
        end
    end

    # Combine TIM TOAs with residuals (if we have any)
    combined_entries = nothing
    if residuals !== nothing
        try
            combined_entries = combine_tim_and_residuals(
                tim_entries, residuals;
                time_start = time_start,
                time_finish = time_finish,
            )
        catch err
            @warn "Failed to combine TIM and residuals" error=err
            combined_entries = nothing
        end
    end

    # Residual statistics (if we have combined data)
    stats = nothing
    if combined_entries !== nothing && !isempty(combined_entries)
        try
            stats = build_residual_statistics_group(combined_entries)
        catch err
            @warn "Failed to build residual statistics" error=err
            stats = nothing
        end
    end

    # White-noise fit (optional, requires combined data)
    white_noise_fit = nothing
    if analyze_white_noise && combined_entries !== nothing && !isempty(combined_entries)
        try
            white_noise_fit = build_white_noise_fit(combined_entries)
        catch err
            @warn "White-noise fit failed" error=err
            white_noise_fit = nothing
        end
    end

    return InternalIterationResult(
        output,
        save_residuals ? residuals : nothing,
        stats,
        white_noise_fit,
        Dict{Symbol, Any}(),
    )
end

#--------------------------------------------------------------------------------------------------------------
# Convergence tracking
#--------------------------------------------------------------------------------------------------------------

struct ConvergenceSeries
    values::Vector{Float64}     # metric per iteration
    abs_deltas::Vector{Float64} # |Δₙ| = |xₙ - xₙ₋₁|
    rel_deltas::Vector{Float64} # |Δₙ / xₙ₋₁| with safe denom
    final_abs_delta::Float64    # NaN if not defined (fewer than 2 points)
    final_rel_delta::Float64    # NaN if not defined
end

function build_convergence_series(values::Vector{Float64})::ConvergenceSeries
    N = length(values)
    if N < 2
        return ConvergenceSeries(values, Float64[], Float64[], NaN, NaN)
    end
    abs_deltas = Vector{Float64}(undef, N-1)
    rel_deltas = Vector{Float64}(undef, N-1)
    @inbounds for i in 1:N-1
        d = values[i+1] - values[i]
        abs_deltas[i] = abs(d)
        denom = max(abs(values[i]), eps())   # stable division
        rel_deltas[i] = abs(d) / denom
    end
    return ConvergenceSeries(values, abs_deltas, rel_deltas, abs_deltas[end], rel_deltas[end])
end

#--------------------------------------------------------------------------------------------------------------

function find_worst_parameter(fit_params::Vector{FitParameter})
    worst_ratio = 0.0
    worst_parameter = nothing
    for p in fit_params
        if p.fit_flag && isfinite(p.uncertainty) && p.uncertainty > 0 && isfinite(p.difference)
            ratio = p.difference / p.uncertainty
            if abs(ratio) > abs(worst_ratio)
                worst_ratio = ratio
                worst_parameter = (
                    name        = p.name_symbol,
                    delta       = p.difference,
                    uncertainty = p.uncertainty,
                    ratio       = ratio,
                )
            end
        end
    end
    return worst_parameter
end

#--------------------------------------------------------------------------------------------------------------

struct ConvergenceInfo
    wrms::ConvergenceSeries
    wrms_tn::ConvergenceSeries
    chisqr::ConvergenceSeries
    worst_parameter::Union{NamedTuple, Nothing}

    final_pre_post::Float64                  # NaN if unavailable
    threshold::Dict{Symbol, Float64}
    converged::Bool
end

# pretty print
function Base.show(io::IO, info::ConvergenceInfo)
    fmtv(x) = isfinite(x) ? @sprintf("%.6f", x) : "-"
    fmtδ(x) = isfinite(x) ? @sprintf("%.6e", x) : "-"

    println(io, "ConvergenceInfo:")
    println(io, "  ✓ Converged: ", info.converged ? "✅ Yes" : "❌ No")
    println(io)

    for (name, conv) in [(:wrms, info.wrms), (:wrms_tn, info.wrms_tn), (:chisqr, info.chisqr)]
        last_val = isempty(conv.values) ? NaN : conv.values[end]
        println(io, "  ── $name ───────────────────────────")
        println(io, "     last value        : ", fmtv(last_val))
        println(io, "     abs delta         : ", fmtδ(conv.final_abs_delta))
        println(io, "     rel delta         : ", fmtδ(conv.final_rel_delta))

        thresh_abs = get(info.threshold, Symbol("abs_$name"), nothing)
        thresh_rel = get(info.threshold, Symbol("rel_$name"), nothing)
        if thresh_abs !== nothing
            println(io, "     threshold abs     : ", @sprintf("%.1e", thresh_abs))
        end
        if thresh_rel !== nothing
            println(io, "     threshold rel     : ", @sprintf("%.1e", thresh_rel))
        end
        println(io)
    end

    println(io, "  ── final_pre_post ───────────────")
    println(io, "     value             : ", isfinite(info.final_pre_post) ? @sprintf("%.8f", info.final_pre_post) : "-")
    if haskey(info.threshold, :final_pre_post)
        println(io, "     threshold         : ", @sprintf("%.1e", info.threshold[:final_pre_post]))
    end
    println(io)

    if info.worst_parameter !== nothing
        wp = info.worst_parameter
        println(io, "  ── worst_parameter ─────────────────────")
        println(io, "     name              : ", wp.name)
        println(io, "     Δ (difference)    : ", @sprintf("%.6g", wp.delta))
        println(io, "     σ (uncertainty)   : ", @sprintf("%.6g", wp.uncertainty))
        println(io, "     Δ / σ             : ", @sprintf("%.3g", wp.ratio))
    end
end

# map a threshold key → actual last-delta; return NaN if not defined
function get_convergence_metric(info::ConvergenceInfo, key::Symbol)::Float64
    if key === :abs_wrms
        return info.wrms.final_abs_delta
    elseif key === :rel_wrms
        return info.wrms.final_rel_delta
    elseif key === :abs_wrms_tn
        return info.wrms_tn.final_abs_delta
    elseif key === :rel_wrms_tn
        return info.wrms_tn.final_rel_delta
    elseif key === :abs_chisqr
        return info.chisqr.final_abs_delta
    elseif key === :rel_chisqr
        return info.chisqr.final_rel_delta
    elseif key === :final_pre_post
        return isfinite(info.final_pre_post) ? abs(info.final_pre_post - 1.0) : NaN
    else
        throw(ArgumentError("Unknown convergence key: $key"))
    end
end

function default_convergence_thresholds()::Dict{Symbol, Float64}
    Dict(
        :abs_wrms_tn    => 1e-2,
        :abs_chisqr     => 1e-2,
        :final_pre_post => 1e-6
    )
end

function is_converged(info::ConvergenceInfo)::Bool
    for (key, limit) in info.threshold
        actual = get_convergence_metric(info, key)
        if isnan(actual) || actual > limit
            return false
        end
    end
    return true
end

function is_converged_by(info::ConvergenceInfo, keys::Symbol...)
    all(key -> begin
        limit  = get(info.threshold, key, Inf)
        actual = get_convergence_metric(info, key)
        isfinite(actual) && actual <= limit
    end, keys)
end

function build_convergence_info(
    iterations::Vector{InternalIterationResult};
    threshold::Dict{Symbol, Float64} = default_convergence_thresholds()
)::ConvergenceInfo

    wrms_vals    = Float64[]
    wrms_tn_vals = Float64[]
    chisqr_vals  = Float64[]

    for iter in iterations
        if isnothing(iter.stats); continue; end
        push!(wrms_vals,    iter.stats.in_fit.all.raw.wrms)
        push!(wrms_tn_vals, iter.stats.in_fit.all.tn.wrms)
        push!(chisqr_vals,  iter.stats.in_fit.all.norm_global.chisqr)
    end

    wrms   = build_convergence_series(wrms_vals)
    wrms_tn = build_convergence_series(wrms_tn_vals)
    chisqr = build_convergence_series(chisqr_vals)

    # last successful iteration (no TEMPO error)
    idx_ok = findlast(it -> !iserror(it.output.error), iterations)
    if idx_ok === nothing
        info = ConvergenceInfo(wrms, wrms_tn, chisqr, nothing, NaN, threshold, false)
        return ConvergenceInfo(wrms, wrms_tn, chisqr, nothing, NaN, threshold, is_converged(info))
    end

    last_ok = iterations[idx_ok]
    final_pre_post = last_ok.output.basic.pre_post
    final_pre_post = isfinite(final_pre_post) ? final_pre_post : NaN
    worst_parameter = find_worst_parameter(last_ok.output.fit_parameters)

    info = ConvergenceInfo(wrms, wrms_tn, chisqr, worst_parameter, final_pre_post, threshold, false)
    return ConvergenceInfo(wrms, wrms_tn, chisqr, worst_parameter, final_pre_post, threshold, is_converged(info))
end

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
        m[:chi2r_fit]    = fit_all.norm_global.rchisqr

        # Add _tim metrics only if requested and not a trivial alias (or forced)
        add_tim = include_tim && (final_it.stats.in_fit !== final_it.stats.in_tim || include_tim_even_if_same)
        if add_tim
            tim_all = final_it.stats.in_tim.all
            m[:wrms_tim]     = tim_all.raw.wrms
            m[:wrms_tn_tim]  = tim_all.tn.wrms
            m[:chi2_tim]     = tim_all.norm_global.chisqr
            m[:chi2r_tim]    = tim_all.norm_global.rchisqr
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
    println(io, "GeneralTempoResult")
    println(io, "  iterations    : ", length(r.iterations))
    println(io, "  final_index   : ", r.final_index)
    println(io, "  converged     : ", r.convergence.converged ? "yes" : "no")
    if r.par_file_final !== nothing
        println(io, "  par_file_final: ", r.par_file_final.path)
    else
        println(io, "  par_file_final: -")
    end
    println(io, "  params        : ", length(r.param_estimates))
    println(io, "  residuals     : ", r.residuals === nothing ? "-" : "present")
    println(io, "  white_noise   : ", r.white_noise === nothing ? "-" : "present")
    println(io, "  metrics       : ", isempty(r.metrics) ? "-" : string(length(r.metrics), " keys"))
    if !isempty(r.subresults)
        println(io, "  subresults    : ", length(r.subresults), " (type=", r.subresult_type, ")")
    else
        println(io, "  subresults    : -")
    end
    isempty(r.metadata) || println(io, "  metadata keys : ", join(keys(r.metadata), ", "))
end

Base.show(io::IO, r::GeneralTempoResult) = show(io, MIME"text/plain"(), r)

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