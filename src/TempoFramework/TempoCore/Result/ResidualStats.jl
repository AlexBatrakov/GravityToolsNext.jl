# src/TempoFramework/TempoCore/Result/ResidualStats.jl
# -----------------------------------------------------------------------------
# Residual statistics: basic, normalized, grouping by backend, and utilities.
# -----------------------------------------------------------------------------
# Expects the parent module to import:
#   using Printf
#   using Statistics
#   using StatsBase: skewness, kurtosis
#   using HypothesisTests: OneSampleADTest, ApproximateOneSampleKSTest,
#                          JarqueBeraTest, ShapiroWilkTest, pvalue
#   using Distributions: Normal
# -----------------------------------------------------------------------------

# ------------------------------------------------------------
# Basic (non-normalized) residual stats
# ------------------------------------------------------------

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

Base.summary(io::IO, s::BasicResidualStats) =
    print(io, "BasicResidualStats(n=", s.n,
              ", mean=", @sprintf("%.4f", s.mean),
              ", std=",  @sprintf("%.4f", s.std), ")")

function Base.show(io::IO, s::BasicResidualStats)
    # Compact one-liner for container contexts; fuller one-liner otherwise
    compact = get(io, :compact, true)
    fmt(x) = @sprintf("%.5f", x)

    if compact
        # very short summary for listings (Dict/Vector elements etc.)
        print(io, "BasicResidualStats(")
        print(io, "n=", s.n,
                  ", mean=",   fmt(s.mean),
                  ", std=",    fmt(s.std),
                  ", wrms=",   fmt(s.wrms))
        print(io, ")")
    else
        # complete one-liner with all fields, but still single-line
        print(io, "BasicResidualStats(")
        print(io, "n=", s.n, ", ")
        print(io, "min=",   fmt(s.min),   ", ")
        print(io, "max=",   fmt(s.max),   ", ")
        print(io, "mean=",  fmt(s.mean),  ", ")
        print(io, "wmean=", fmt(s.wmean), ", ")
        print(io, "median=",fmt(s.median),", ")
        print(io, "rms=",   fmt(s.rms),   ", ")
        print(io, "wrms=",  fmt(s.wrms),  ", ")
        print(io, "std=",   fmt(s.std),   ", ")
        print(io, "wstd=",  fmt(s.wstd),  ", ")
        print(io, "mad=",   fmt(s.mad))
        print(io, ")")
    end
end

function Base.show(io::IO, ::MIME"text/plain", s::BasicResidualStats)
    indent = get(io, :indent, 0)
    pad  = repeat(" ", indent)
    spad = repeat(" ", indent + 2)

    println(io, pad,  "Basic residual stats:")
    println(io, spad, "n      = ", s.n)
    println(io, spad, "min    = ", @sprintf("%.5f", s.min),    "   max   = ", @sprintf("%.5f", s.max))
    println(io, spad, "mean   = ", @sprintf("%.5f", s.mean),   "   wmean = ", @sprintf("%.5f", s.wmean))
    println(io, spad, "median = ", @sprintf("%.5f", s.median), "   mad   = ", @sprintf("%.5f", s.mad))
    println(io, spad, "rms    = ", @sprintf("%.5f", s.rms),    "   wrms  = ", @sprintf("%.5f", s.wrms))
    println(io, spad, "std    = ", @sprintf("%.5f", s.std),    "   wstd  = ", @sprintf("%.5f", s.wstd))
end

# Unweighted variant
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

# Weighted variant
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

# ------------------------------------------------------------
# Normalized residual stats (ideally ~ N(0,1))
# ------------------------------------------------------------

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
    red_chisqr::Float64
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


Base.summary(io::IO, s::NormalizedResidualStats) =
    print(io, "NormalizedResidualStats(n=", s.n,
              ", mean=",  @sprintf("%.4f", s.mean),
              ", std=",   @sprintf("%.4f", s.std),
              ", rχ2=",   @sprintf("%.4f", s.red_chisqr), ")")

function Base.show(io::IO, s::NormalizedResidualStats)
    # Compact one-liner for container contexts; fuller one-liner otherwise
    compact = get(io, :compact, true)
    fmt(x) = @sprintf("%.5f", x)

    if compact
        # very short summary for listings (Dict/Vector elements etc.)
        print(io, "NormalizedResidualStats(")
        print(io, "n=", s.n,
                  ", mean=",   fmt(s.mean),
                  ", std=",    fmt(s.std))
        print(io, ")")
    else
        # complete one-liner with all fields, but still single-line
        print(io, "NormalizedResidualStats(")
        print(io, "n=", s.n, ", ")
        print(io, "min=",   fmt(s.min),   ", ")
        print(io, "max=",   fmt(s.max),   ", ")
        print(io, "mean=",  fmt(s.mean),  ", ")
        print(io, ")")
    end
end

# pretty-print helper
_show_stat(io, pad, label, val, expect, sigma) = begin
    z = (val - expect) / sigma
    println(io, pad, @sprintf("%-12s = %10.5f   (%+6.2fσ)", label, val, z))
end

function Base.show(io::IO, ::MIME"text/plain", s::NormalizedResidualStats)
    indent = get(io, :indent, 0)
    pad  = repeat(" ", indent)
    spad = repeat(" ", indent + 2)
    tpad = repeat(" ", indent + 4)

    n = s.n
    println(io, pad,  "Normalized residual statistics:")
    println(io, spad, "Summary:")
    println(io, tpad, @sprintf("%-12s = %10d", "count", n))
    _show_stat(io, tpad, "minimum",  s.min,      0.0, s.std)
    _show_stat(io, tpad, "maximum",  s.max,      0.0, s.std)
    _show_stat(io, tpad, "mean",     s.mean,     0.0, 1 / sqrt(max(n, 1)))
    _show_stat(io, tpad, "median",   s.median,   0.0, 1.253 / sqrt(max(n, 1)))
    _show_stat(io, tpad, "std",      s.std,      1.0, 1 / sqrt(2 * max(n - 1, 1)))
    _show_stat(io, tpad, "skewness", s.skewness, 0.0, sqrt(6 / max(n, 1)))
    _show_stat(io, tpad, "kurtosis", s.kurtosis, 0.0, sqrt(24 / max(n, 1)))
    println(io, tpad, @sprintf("%-12s = %10.5f   (red = %.5f)", "chisqr", s.chisqr, s.red_chisqr))

    println(io, spad, "Quantiles:")
    println(io, tpad, @sprintf("q01=%.3f  q05=%.3f  q95=%.3f  q99=%.3f", s.q01, s.q05, s.q95, s.q99))

    println(io, spad, "Empirical coverage:")
    p1 = s.within1σ / max(n, 1)
    p2 = s.within2σ / max(n, 1)
    p3 = s.within3σ / max(n, 1)
    println(io, tpad, @sprintf("|x|≤1σ : %5d / %d  (%.3f, expected %.3f)", s.within1σ, n, p1, _P68))
    println(io, tpad, @sprintf("|x|≤2σ : %5d / %d  (%.3f, expected %.3f)", s.within2σ, n, p2, _P95))
    println(io, tpad, @sprintf("|x|≤3σ : %5d / %d  (%.3f, expected %.3f)", s.within3σ, n, p3, _P997))

    println(io, spad, "Normality tests:")
    println(io, tpad, @sprintf("%-22s = %10.5f   (p = %6.5f)", "Anderson-Darling  A²", s.ad_statistic, s.ad_p_value))
    println(io, tpad, @sprintf("%-22s = %10.5f   (p = %6.5f)", "Kolmogorov-Smirnov D", s.ks_statistic, s.ks_p_value))
    println(io, tpad, @sprintf("%-22s = %10.5f   (p = %6.5f)", "Jarque-Bera        JB", s.jb_statistic, s.jb_p_value))
    println(io, tpad, @sprintf("%-22s = %10.5f   (p = %6.5f)", "Shapiro-Wilk       W",  s.sw_statistic, s.sw_p_value))
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
    χ2   = sum(abs2, xf)
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

"""
    _empty_norm_stats() -> NormalizedResidualStats

Return a neutral/empty `NormalizedResidualStats` so downstream code doesn't crash on empty groups.
Internal helper.
"""
function _empty_norm_stats()::NormalizedResidualStats
    return NormalizedResidualStats(
        0, NaN, NaN, NaN, NaN, NaN, NaN, NaN,     # n, min, max, mean, median, std, skew, kurt
        NaN, NaN,                                 # chisqr, red_chisqr
        NaN, NaN, NaN, NaN,                       # q01, q05, q95, q99
        0, 0, 0,                                  # within1σ, within2σ, within3σ
        NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN    # AD, KS, JB, SW stats+pvals
    )
end

# ------------------------------------------------------------
# Aggregated residual statistics (per-iteration, per-backend)
# ------------------------------------------------------------

struct ResidualStatistics
    raw::BasicResidualStats
    tn::BasicResidualStats
    norm_global::NormalizedResidualStats
    norm_local::NormalizedResidualStats
end

function Base.show(io::IO, ::MIME"text/plain", stats::ResidualStatistics)
    indent = get(io, :indent, 0)
    pad  = repeat(" ", indent)
    spad = repeat(" ", indent + 2)
    iop  = IOContext(io, :indent => indent + 4)

    println(io, pad, "Residual Statistics:")

    println(io, spad, "raw → Raw (residual):")
    show(iop, MIME"text/plain"(), stats.raw)

    println(io, spad, "tn → Whitened (residual_tn):")
    show(iop, MIME"text/plain"(), stats.tn)

    println(io, spad, "norm_global → Globally centered normalized whitened ((residual_tn - global_wmean) / uncertainty):")
    show(iop, MIME"text/plain"(), stats.norm_global)

    println(io, spad, "norm_local → Locally centered normalized whitened ((residual_tn - local_wmean) / uncertainty):")
    show(iop, MIME"text/plain"(), stats.norm_local)
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
    sizehint!(norm_global, n)
    sizehint!(norm_local,  n)

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

# ------------------------------------------------------------
# Grouping by backend and by in-fit / in-tim
# ------------------------------------------------------------

struct ResidualStatisticsEntry
    all::ResidualStatistics
    by_backend::Dict{Symbol, ResidualStatistics}
end

"""
    _group_by_backend(entries; in_fit=nothing) -> Dict{Symbol, Vector{CombinedTOAEntry}}

Group entries by backend (Symbol). Optionally filter by `in_fit`:
- `nothing` : no filtering
- `true`    : only in-fit
- `false`   : only not in-fit
Internal helper.
"""
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

function Base.show(io::IO, ::MIME"text/plain", entry::ResidualStatisticsEntry)
    indent = get(io, :indent, 0)
    pad  = repeat(" ", indent)
    spad = repeat(" ", indent + 2)
    tpad = repeat(" ", indent + 4)
    iop  = IOContext(io, :indent => indent + 4)

    println(io, pad, "all → All TOAs:")
    show(iop, MIME"text/plain"(), entry.all)

    if isempty(entry.by_backend)
        println(io, spad, "by_backend → Backends: (none)")
        return
    end

    println(io, spad, "by_backend → Backends:")
    for be in sort!(collect(keys(entry.by_backend)))
        stats = entry.by_backend[be]
        n     = stats.norm_global.n
        wrms  = stats.tn.wrms
        println(io, tpad, "[", String(be), "] n=", n, "  wrms_tn=", @sprintf("%.3f", wrms))
    end
end

struct ResidualStatisticsGroup
    in_fit::ResidualStatisticsEntry
    in_tim::ResidualStatisticsEntry
end

function Base.show(io::IO, ::MIME"text/plain", group::ResidualStatisticsGroup)
    indent = get(io, :indent, 0)
    pad  = repeat(" ", indent)
    spad = repeat(" ", indent + 2)
    tpad = repeat(" ", indent + 4)
    iop  = IOContext(io, :indent => indent + 4)

    println(io, pad, "Residual Statistics Group:")

    n_fit = group.in_fit.all.norm_global.n
    println(io, spad, "in_fit → In-Fit (", n_fit, " TOAs):")
    show(iop, MIME"text/plain"(), group.in_fit)

    n_tim = group.in_tim.all.norm_global.n
    if n_fit != n_tim
        println(io, spad, "in_tim → In-TIM (", n_tim, " TOAs):")
        show(iop, MIME"text/plain"(), group.in_tim)
    end
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