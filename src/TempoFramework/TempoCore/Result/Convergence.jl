# src/TempoFramework/TempoCore/Result/Convergence.jl
# -----------------------------------------------------------------------------
# Convergence time-series, deltas, thresholds and summary helpers.
# -----------------------------------------------------------------------------
# defines:
#   ConvergenceSeries, ConvergenceInfo
#   build_convergence_series, build_convergence_info
#   get_convergence_metric, default_convergence_thresholds
#   is_converged, is_converged_by, find_worst_parameter

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