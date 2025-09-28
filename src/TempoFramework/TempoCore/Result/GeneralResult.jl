# src/TempoFramework/TempoCore/Result/GeneralResult.jl
# -----------------------------------------------------------------------------
# Top-level GeneralTempoResult container and constructors.
# -----------------------------------------------------------------------------
# defines:
#   ParamEstimate (alias)
#   GeneralTempoResult (struct)
#   GeneralTempoResult(iterations; ...) constructor
#   GeneralTempoResult(iter; ...) convenience constructor
#   (virtual props) res.final, res.last_successful, res.residual_stats, res.white_noise
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
function build_core_metrics(final_iter::InternalIterationResult, conv::ConvergenceInfo;
                            include_tim::Bool=true,
                            include_tim_even_if_same::Bool=false)
    metrics = Dict{Symbol,Float64}()

    # Final TEMPO basic block
    basic = final_iter.output.basic

    metrics[:chi2_fit_basic]  = isfinite(basic.chisqr) ? basic.chisqr : NaN

    # Residual stats (may be missing)
    if final_iter.stats !== nothing
        fit_all = final_iter.stats.in_fit.all
        metrics[:wrms_fit]     = fit_all.raw.wrms
        metrics[:wrms_tn_fit]  = fit_all.tn.wrms
        metrics[:chi2_fit]     = fit_all.norm_global.chisqr
        metrics[:chi2r_fit]    = fit_all.norm_global.red_chisqr

        # Add _tim metrics only if requested and not a trivial alias (or forced)
        add_tim = include_tim && (final_iter.stats.in_fit !== final_iter.stats.in_tim || include_tim_even_if_same)
        if add_tim
            tim_all = final_iter.stats.in_tim.all
            metrics[:wrms_tim]     = tim_all.raw.wrms
            metrics[:wrms_tn_tim]  = tim_all.tn.wrms
            metrics[:chi2_tim]     = tim_all.norm_global.chisqr
            metrics[:chi2r_tim]    = tim_all.norm_global.red_chisqr
        end
    else
        # No stats available
        metrics[:wrms_fit]     = NaN
        metrics[:wrms_tn_fit]  = NaN
        metrics[:chi2_fit]     = NaN
        metrics[:chi2r_fit]    = NaN
        # Note: _tim keys are omitted entirely when stats are missing
    end

    # Final pre/post from TEMPO (may be NaN)
    metrics[:pre_post_final] = conv.pre_post_final

    # Convergence deltas between the last two points (NaN if <2 points)
    metrics[:delta_wrms_tn] = conv.wrms_tn.final_abs_delta
    metrics[:delta_chi2]    = conv.chisqr.final_abs_delta

    # Global AD after white-noise fit (if performed)
    if final_iter.white_noise_fit !== nothing
        metrics[:ad_white_fit] = final_iter.white_noise_fit.global_stats.ad_statistic
    else
        metrics[:ad_white_fit] = NaN
    end

    return metrics
end

"""
    GeneralTempoResult

Top-level container with everything you usually want after a run.

Fields (stored)
- `iterations`            : all internal iterations
- `final_index`           : index of the *last* iteration (== `length(iterations)`)
- `last_successful_index` : 0 if there were no successful iterations
- `success`               : quick boolean flag; true iff the run is considered successful
- `status`                : low-level run status (`:ok | :engine_failed | :parse_failed | :files_missing | :unknown`)
- `convergence`           : convergence summary across iterations
- `metrics`               : compact scalar metrics for quick ranking/comparison
- `param_estimates`       : map `:NAME => (value, uncertainty)` extracted from the *final* iteration
- `par_file_final`        : written output par-file (if available)
- `subresults`            : optional nested results (e.g., per-epoch, per-band, grid cells)
- `subresult_type`        : tag describing what `subresults` mean (e.g., `:epoch`, `:band`, `:grid`)
- `metadata`              : extra info (paths, timings, seeds, etc.)

Virtual properties (accessible via `res.final` / `res.last_successful`):
- `final`                 : alias for `res.iterations[res.final_index]`
- `last_successful`       : `nothing` if `last_successful_index == 0`, otherwise
                            `res.iterations[res.last_successful_index]`
- `residual_stats`        : alias for `res.final.stats`
- `white_noise`           : alias for `res.final.white_noise_fit`
  (`white_noise_fit` is also available)
"""
struct GeneralTempoResult
    # core timeline
    iterations::Vector{InternalIterationResult}
    final_index::Int                  # == length(iterations)
    last_successful_index::Int        # 0 if there were no successful iterations

    # run verdict
    success::Bool
    status::Symbol                    # :ok | :engine_failed | :parse_failed | :files_missing | :unknown

    # analysis
    convergence::ConvergenceInfo
    metrics::Dict{Symbol, Float64}    # compact numeric metrics
    param_estimates::Dict{Symbol, ParamEstimate}

    # artifacts from the (true) final iteration
    par_file_final::Union{TempoParFile, Nothing}

    # optional nesting
    subresults::Vector{GeneralTempoResult}
    subresult_type::Union{Symbol, Nothing}

    # everything else (paths, timings, policies…)
    metadata::Dict{Symbol, Any}
end

# --- virtual properties for convenient access ---
function Base.getproperty(r::GeneralTempoResult, s::Symbol)
    if s === :final
        return getfield(r, :iterations)[getfield(r, :final_index)]
    elseif s === :last_successful
        lsi = getfield(r, :last_successful_index)
        return lsi == 0 ? nothing : getfield(r, :iterations)[lsi]
    elseif s === :residual_stats
        return getproperty(getproperty(r, :final), :stats)
    elseif s === :white_noise_fit
        return getproperty(getproperty(r, :final), :white_noise_fit)
    else
        return getfield(r, s)
    end
end

# --- property names for REPL/tab-completion ---
function Base.propertynames(r::GeneralTempoResult, private::Bool=false)
    # real struct fields
    names = fieldnames(typeof(r))
    # virtual properties we provide via getproperty
    return (names..., :final, :last_successful, :residual_stats, :white_noise_fit)
end

# # Keep hasproperty in sync with propertynames/getproperty
Base.hasproperty(r::GeneralTempoResult, s::Symbol) = (s in propertynames(r))

# ---------- pretty print ----------
function Base.show(io::IO, ::MIME"text/plain", r::GeneralTempoResult)
    indent = get(io, :indent, 0)
    pad  = repeat(" ", indent)
    spad = repeat(" ", indent + 2)
    iop  = IOContext(io, :indent => indent + 4)

    println(io, pad,  "GeneralTempoResult:")
    println(io, spad, "iterations:            ", length(r.iterations), " entries")
    println(io, spad, "final_index:           ", r.final_index)
    println(io, spad, "final:                 ", (r.final_index > 0 && length(r.iterations) >= r.final_index) ? "present" : "-")
    println(io, spad, "last_successful_index: ", r.last_successful_index == 0 ? "-" : string(r.last_successful_index, r.last_successful_index == r.final_index ? " (==final)" : ""))
    println(io, spad, "last_successful:       ", (r.last_successful_index == 0 ? "-" : "present"))

    println(io, spad, "success:               ", r.success ? "yes" : "no")
    println(io, spad, "status:                ", r.status)

    println(io, spad, "convergence:           ", r.convergence.converged ? "yes" : "no")
    println(io, spad, "metrics:               ", isempty(r.metrics) ? "-" : string(length(r.metrics), " keys"))
    println(io, spad, "param_estimates:       ", isempty(r.param_estimates) ? "-" : string(length(r.param_estimates), " keys"))

    println(io, spad, "par_file_final:        ", r.par_file_final === nothing ? "-" : "present")

    if !isempty(r.subresults)
        println(io, spad, "subresults:            ", length(r.subresults), " (type=", r.subresult_type, ")")
    else
        println(io, spad, "subresults:            -")
    end

    println(io, spad, "metadata keys:         ", isempty(r.metadata) ? "-" : string(length(r.metadata), " keys"))
end


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

function _extract_param_estimates(fit_params::Nothing)
    return Dict{Symbol, ParamEstimate}()
end

# Find last successful iteration (no TEMPO error and has stats); return 0 if none
function _last_successful_index(iters::Vector{InternalIterationResult})::Int
    idx = findlast(ir -> (!iserror(ir.output.error) && ir.stats !== nothing), iters)
    return idx === nothing ? 0 : idx
end

# ---------- main builder ----------
"""
    build_general_tempo_result(
        iterations;
        par_file_final=nothing,
        subresults=GeneralTempoResult[],
        subresult_type=nothing,
        metadata=Dict(),
        metrics_hook=nothing
    )

Build a `GeneralTempoResult` from a list of `InternalIterationResult`s.
- Uses the *last* iteration as `final` (index `final_index = length(iterations)`).
- Computes `last_successful_index` (last iteration with no error and with stats); 0 if none.
- Computes convergence across all iterations with available stats.
- Extracts parameter estimates from the *final* iteration.
- Computes `metrics` via `build_core_metrics(final, convergence)` and optionally merges
  a user-supplied `metrics_hook(final, convergence)` dictionary (values must be `Real`).
"""
function build_general_tempo_result(
    iterations::Vector{InternalIterationResult};
    par_file_final::Union{TempoParFile,Nothing}=nothing,
    subresults::Vector{GeneralTempoResult}=GeneralTempoResult[],
    subresult_type::Union{Symbol,Nothing}=nothing,
    metadata::Dict{Symbol,Any}=Dict{Symbol,Any}(),
    metrics_hook::Union{Nothing,Function}=nothing,
)
    isempty(iterations) && error("build_general_tempo_result: iterations are empty")

    final_index = length(iterations)
    final_iteration  = iterations[final_index]

    last_success_index = _last_successful_index(iterations)  # Int, 0 if none

    # Success policy: prefer explicit metadata keys; fall back to iteration status
    success = haskey(metadata, :success) ? Bool(metadata[:success]) :
          haskey(metadata, :status)  ? (metadata[:status] == :ok) :
          !iserror(final_iteration.output.error)

    status  = get(metadata, :status, :unknown)

    conv   = build_convergence_info(iterations)
    param_estimates = _extract_param_estimates(final_iteration.output.fit_parameters)
    metrics = build_core_metrics(final_iteration, conv)

    if metrics_hook !== nothing
        extra_any = metrics_hook(final_iteration, conv)
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
        merge!(metrics, extra)  # user values override core if keys collide
    end

    return GeneralTempoResult(
        iterations,
        final_index,
        last_success_index,
        success,
        status,
        conv,
        metrics,
        param_estimates,
        par_file_final,
        subresults,
        subresult_type,
        metadata,
    )
end

"Safe lookup of a scalar metric from `GeneralTempoResult.metrics`."

result_metric(res::GeneralTempoResult, key::Symbol)::Float64 =
    get(res.metrics, key, NaN)





