# src/TempoFramework/TempoCore/Result/Convergence.jl
# -----------------------------------------------------------------------------
# Convergence time-series, deltas, thresholds and summary helpers.
#
# This module builds convergence time series from per-iteration results, computes
# last-step deltas and applies user-defined thresholds to decide whether a run
# "converged". It is intentionally conservative:
#  - time-series use only iterations that have valid residual statistics;
#  - `pre_post_final` is taken strictly from the last iteration (if unavailable → NaN);
#  - the overall flag `converged` requires all relevant thresholds to be satisfied;
#  - when only a single iteration is available, convergence falls back to the
#    `pre_post_final` check.
#
# The worst (least stable) fitted parameter is reported from the last successful
# iteration that provides a fit-parameter table.
#

"""
    ConvergenceSeries

Holds a metric's values across iterations together with absolute/relative
last-step deltas. If fewer than two points are available, `final_*_delta = NaN`.
"""
struct ConvergenceSeries
    values::Vector{Float64}     # metric per iteration
    abs_deltas::Vector{Float64} # |Δₙ| = |xₙ - xₙ₋₁|
    rel_deltas::Vector{Float64} # |Δₙ / xₙ₋₁| with safe denom
    final_abs_delta::Float64    # NaN if not defined (fewer than 2 points)
    final_rel_delta::Float64    # NaN if not defined
end

"""
    build_convergence_series(values) -> ConvergenceSeries

Constructs absolute and relative last-step deltas for a numeric sequence.
Relative delta uses a stable denominator `max(abs(prev), eps())`.
"""
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

"""
    WorstParam

Named tuple describing the least stable fitted parameter:
`(name::Symbol, delta::Float64, uncertainty::Float64, ratio::Float64)`
where `ratio = delta / uncertainty`.
"""
const WorstParam = NamedTuple{(:name,:delta,:uncertainty,:ratio),Tuple{Symbol,Float64,Float64,Float64}}

"""
    find_worst_parameter(fit_params) -> WorstParam | nothing

Finds the parameter with the largest absolute ratio `|Δ/σ|` among parameters
with `fit_flag = true` and finite uncertainties (`σ > 0`). Returns a named tuple
`WorstParam` or `nothing` if no suitable parameters are present.
"""
function find_worst_parameter(fit_params::Vector{FitParameter})::Union{WorstParam, Nothing}
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

"""
    ConvergenceInfo

Summary of convergence diagnostics.
- `wrms`, `wrms_tn`, `chisqr` — time series with last-step deltas.
- `pre_post_final` — value taken **strictly from the last iteration** (NaN if unavailable).
- `worst_parameter` — `WorstParam` from the last successful iteration providing a
  fit-parameter table; `nothing` otherwise.
- `thresholds` — dictionary of limits used by `is_converged`.
- `converged` — final decision under the given thresholds.
"""
struct ConvergenceInfo
    wrms::ConvergenceSeries
    wrms_tn::ConvergenceSeries
    chisqr::ConvergenceSeries

    pre_post_final::Float64                  # NaN if unavailable
    worst_parameter::Union{WorstParam, Nothing}
    thresholds::Dict{Symbol, Float64}
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

        thresh_abs = get(info.thresholds, Symbol("abs_$name"), nothing)
        thresh_rel = get(info.thresholds, Symbol("rel_$name"), nothing)
        if thresh_abs !== nothing
            println(io, "     threshold abs     : ", @sprintf("%.1e", thresh_abs))
        end
        if thresh_rel !== nothing
            println(io, "     threshold rel     : ", @sprintf("%.1e", thresh_rel))
        end
        println(io)
    end

    println(io, "  ── pre_post_final ───────────────")
    println(io, "     value             : ", isfinite(info.pre_post_final) ? @sprintf("%.8f", info.pre_post_final) : "-")
    if haskey(info.thresholds, :pre_post_final)
        println(io, "     threshold         : ", @sprintf("%.1e", info.thresholds[:pre_post_final]))
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

"""
    get_convergence_metric(info, key) -> Float64

Maps a threshold key to the corresponding last-step delta (or |pre_post_final−1|).
Returns `NaN` if the metric is not available.
"""
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
    elseif key === :pre_post_final
        return isfinite(info.pre_post_final) ? abs(info.pre_post_final - 1.0) : NaN
    else
        throw(ArgumentError("Unknown convergence key: $key"))
    end
end

"""
    default_convergence_thresholds() -> Dict{Symbol,Float64}

Default limits used by `is_converged`:
- `:abs_wrms_tn`  — absolute delta of TN-weighted WRMS
- `:abs_chisqr`   — absolute delta of χ²
- `:pre_post_final` — |pre_post_final − 1|
"""
function default_convergence_thresholds()::Dict{Symbol, Float64}
    Dict(
        :abs_wrms_tn    => 1e-2,
        :abs_chisqr     => 1e-2,
        :pre_post_final => 1e-6
    )
end

"""
    is_converged(info) -> Bool

Checks all keys present in `info.thresholds`. If a required metric is `NaN`,
convergence fails.
"""
function is_converged(info::ConvergenceInfo)::Bool
    for (key, limit) in info.thresholds
        actual = get_convergence_metric(info, key)
        if isnan(actual) || actual > limit
            return false
        end
    end
    return true
end

"""
    is_converged_by(info, keys...) -> Bool

Checks convergence only for the specified keys. Metrics that are not finite are
ignored for the decision.
"""
function is_converged_by(info::ConvergenceInfo, keys::Symbol...)
    all(key -> begin
        limit  = get(info.thresholds, key, Inf)
        actual = get_convergence_metric(info, key)
        isfinite(actual) && actual <= limit
    end, keys)
end

"""
    ConvergenceInfo(wrms, wrms_tn, chisqr, pre_post_final, worst_parameter, thresholds)

Builds a `ConvergenceInfo` and computes `converged`. For sequences with fewer
than two points, convergence falls back to the `:pre_post_final` criterion.
"""
function ConvergenceInfo(
    wrms::ConvergenceSeries,
    wrms_tn::ConvergenceSeries,
    chisqr::ConvergenceSeries,
    pre_post_final::Float64,
    worst_parameter::Union{WorstParam, Nothing},
    thresholds::Dict{Symbol, Float64}
)::ConvergenceInfo
    info = ConvergenceInfo(wrms, wrms_tn, chisqr, pre_post_final, worst_parameter, thresholds, false)
    converged = length(info.wrms.values) > 1 ? is_converged(info) : is_converged_by(info, :pre_post_final)
    return ConvergenceInfo(wrms, wrms_tn, chisqr, pre_post_final, worst_parameter, thresholds, converged)
end

"""
    build_convergence_info(iterations; thresholds=default_convergence_thresholds())
        -> ConvergenceInfo

Constructs convergence time series from iterations that have residual
statistics. `pre_post_final` is taken **strictly from the last iteration**; if
unavailable there, it is `NaN`. The `worst_parameter` is taken from the last
successful iteration that provides a fit-parameter table.

When only one valid point is available, the decision `converged` falls back to
`pre_post_final`.
"""
function build_convergence_info(
    iterations::Vector{InternalIterationResult};
    thresholds::Dict{Symbol, Float64} = default_convergence_thresholds()
)::ConvergenceInfo

    if isempty(iterations)
        empty_series = ConvergenceSeries(Float64[], Float64[], Float64[], NaN, NaN)
        return ConvergenceInfo(empty_series, empty_series, empty_series, NaN, nothing, thresholds)
    end

    wrms_vals    = Float64[]
    wrms_tn_vals = Float64[]
    chisqr_vals  = Float64[]

    for iter in iterations
        if isnothing(iter.stats)
            push!(wrms_vals,    NaN)
            push!(wrms_tn_vals, NaN)
            push!(chisqr_vals,  NaN)
            continue
        end
        push!(wrms_vals,    iter.stats.in_fit.all.raw.wrms)
        push!(wrms_tn_vals, iter.stats.in_fit.all.tn.wrms)
        push!(chisqr_vals,  iter.stats.in_fit.all.norm_global.chisqr)
    end

    wrms    = build_convergence_series(wrms_vals)
    wrms_tn = build_convergence_series(wrms_tn_vals)
    chisqr  = build_convergence_series(chisqr_vals)

    # pre_post_final: strictly from the last iteration; if it's unavailable there → NaN
    last_iter = iterations[end]
    pre_post_final = !isnothing(last_iter.output.basic) ? last_iter.output.basic.pre_post : NaN

    # worst_parameter: from the last iteration that has a fit-parameters table
    idx_wp = findlast(it -> !isnothing(it.output.fit_parameters), iterations)
    worst_parameter = idx_wp === nothing ? nothing : find_worst_parameter(iterations[idx_wp].output.fit_parameters)

    return ConvergenceInfo(wrms, wrms_tn, chisqr, pre_post_final, worst_parameter, thresholds)
end