# src/TempoFramework/TempoCore/Result/InternalIteration.jl
# -----------------------------------------------------------------------------
# Per-internal-iteration container and builder.
# -----------------------------------------------------------------------------
# defines:
#   InternalIterationResult
#   build_internal_iteration_result
# uses:
#   ResidualStatisticsGroup, build_residual_statistics_group
#   WhiteNoiseFitResult, build_white_noise_fit
#   read_residual_file_safe, combine_tim_and_residuals_safe
# -----------------------------------------------------------------------------


"""
    InternalIterationResult

A single TEMPO iteration snapshot:

Fields
- `output`          : parsed TEMPO output blocks (incl. basic/error/parameters)
- `residuals`       : raw residual rows as read from disk (or `nothing` if not kept)
- `stats`           : residual statistics group (in-fit / in-tim), if available
- `white_noise_fit` : per-backend white-noise fit result (optional)
- `metadata`        : extra info (paths, counts, flags, time window, etc.)
"""
struct InternalIterationResult
    output::InternalIterationOutput
    residuals::Union{Vector{TempoResidualEntry}, Nothing}
    stats::Union{ResidualStatisticsGroup, Nothing}
    white_noise_fit::Union{WhiteNoiseFitResult, Nothing}
    metadata::Dict{Symbol, Any}
end

function Base.show(io::IO, ::MIME"text/plain", r::InternalIterationResult)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    ipad   = indent + 4
    iop    = IOContext(io, :indent => indent + 4)

    println(io, pad, "InternalIterationResult:")

    # TEMPO basic block
    println(io, spad, "▸ output:")
    show(iop, MIME"text/plain"(), r.output.basic)
    println(io)

    # TEMPO error block
    if iserror(r.output.error)
        println(io, spad, "▸ error:")
        show(iop, MIME"text/plain"(), r.output.error)
        println(io)
    else
        println(io, spad, "▸ error: none")
    end

    # Stage summary
    println(io, spad, "▸ stages: ",
                 "output=", get(r.metadata, :stage_output_ok, false) ? "ok" : "fail",
                 ", residual_read=", get(r.metadata, :stage_residual_read, :skipped),
                 ", combine=",       get(r.metadata, :stage_combine, :skipped),
                 ", stats=",         get(r.metadata, :stage_stats, :skipped),
                 ", white_noise=",   get(r.metadata, :stage_white_noise, :skipped))
    println(io, spad, "▸ iter_ok: ", get(r.metadata, :iter_ok, false),
                 get(r.metadata, :iter_ok, false) ? "" : " (fail_stage=" * string(get(r.metadata, :iter_fail_stage, :unknown)) * ")")

    # Residual statistics (in-fit only)
    if r.stats === nothing
        println(io, spad, "▸ Residual statistics: none")
    else
        println(io, spad, "▸ Residual statistics (in-fit only):")
        show(iop, MIME"text/plain"(), r.stats.in_fit)
        println(io)
    end

    # # Residual rows presence
    # if r.residuals === nothing
    #     println(io, spad, "▸ Residuals: not stored")
    # else
    #     println(io, spad, "▸ Residuals: ", length(r.residuals), " entries loaded")
    # end

    # White-noise fit
    if r.white_noise_fit === nothing
        println(io, spad, "▸ white_noise_fit: none")
    else
        println(io, spad, "▸ white_noise_fit:")
        show(iop, MIME"text/plain"(), r.white_noise_fit)
        println(io)
    end

    # Metadata
    if isempty(r.metadata)
        println(io, spad, "▸ metadata: none")
    else
        println(io, spad, "▸ metadata: ", join(string.(keys(r.metadata)), ", "))
    end
end

# Delegate generic show to the MIME-specialized version
Base.show(io::IO, r::InternalIterationResult) = show(io, MIME"text/plain"(), r)

"""
    build_internal_iteration_result(output, residual_path, tim_entries; kwargs...) -> InternalIterationResult

Build a single-iteration result from parsed TEMPO output and (optionally) residuals on disk.

Behavior
- Does not return early on TEMPO errors; instead records the error and proceeds
  to attempt residual loading/combination guardedly.
- Reads residuals via a safe wrapper; if the file is missing or unreadable,
  residuals are `nothing`.
- Combines TIM entries with residuals using a safe combiner; on length mismatch
  or missing inputs, combination yields `nothing`.
- Optionally performs white-noise fitting (`analyze_white_noise=true`) on
  combined in-fit data when available.
- Records per-stage flags in `metadata` for downstream diagnostics.

Keywords
- `save_residuals::Bool=false` — keep raw residual rows inside the result
- `analyze_white_noise::Bool=false` — run per-backend white-noise fit
- `time_start::Union{Nothing,Float64}=nothing` — lower MJD bound for combining
- `time_finish::Union{Nothing,Float64}=nothing` — upper MJD bound for combining
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

    # Stage flags and quick verdicts
    stage_output_ok = !iserror(output.error)

    # Residuals: safe read
    residuals::Union{Nothing, Vector{TempoResidualEntry}} = nothing
    residual_path_exists = residual_path !== nothing && isfile(residual_path)
    if residual_path_exists
        residuals = read_residual_file_safe(residual_path)
    end
    stage_residual_read = residual_path === nothing ? :skipped : (residuals === nothing ? :failed : :ok)

    # Combine TIM + residuals: safe combine (non-throwing)
    combined_entries = combine_tim_and_residuals_safe(
        tim_entries,
        residuals;
        time_start = time_start,
        time_finish = time_finish,
    )
    stage_combine = residuals === nothing ? :skipped : (combined_entries === nothing ? :failed : :ok)

    # Residual statistics
    stats = nothing
    stage_stats = :skipped
    stats_error_msg = nothing
    if combined_entries !== nothing && !isempty(combined_entries) && any(isfinite, getfield.(combined_entries, :residual))
        try
            stats = build_residual_statistics_group(combined_entries)
            stage_stats = :ok
        catch err
            @warn "Failed to build residual statistics" error=err
            stats = nothing
            stage_stats = :failed
            stats_error_msg = sprint(showerror, err)
        end
    end

    # White-noise fit (optional, requires combined data)
    white_noise_fit = nothing
    stage_white_noise = :skipped
    wn_error_msg = nothing
    if analyze_white_noise && combined_entries !== nothing && !isempty(combined_entries) && any(isfinite, getfield.(combined_entries, :residual))
        try
            white_noise_fit = build_white_noise_fit(combined_entries)
            stage_white_noise = :ok
        catch err
            @warn "White-noise fit failed" error=err
            white_noise_fit = nothing
            stage_white_noise = :failed
            wn_error_msg = sprint(showerror, err)
        end
    end

    # Time window applied?
    time_window_applied = !(time_start === nothing && time_finish === nothing)

    # Combined MJD range (if any)
    combined_mjd_min = nothing
    combined_mjd_max = nothing
    if combined_entries !== nothing && !isempty(combined_entries)
        combined_mjd_min = minimum(getfield.(combined_entries, :toa))
        combined_mjd_max = maximum(getfield.(combined_entries, :toa))
    end

    # Was white-noise actually performed?
    wn_performed = analyze_white_noise && stage_white_noise == :ok

    # Iteration-level verdict
    iter_ok = stage_output_ok && (stats !== nothing)
    iter_fail_stage = nothing
    if !iter_ok
        # pick the first non-ok stage in a meaningful order
        if !stage_output_ok
            iter_fail_stage = :output
        elseif stage_residual_read == :failed
            iter_fail_stage = :residual_read
        elseif stage_combine == :failed
            iter_fail_stage = :combine
        elseif stage_stats == :failed || (stage_stats == :skipped && stats === nothing)
            iter_fail_stage = :stats
        elseif analyze_white_noise && stage_white_noise == :failed
            iter_fail_stage = :white_noise
        else
            iter_fail_stage = :unknown
        end
    end

    # Metadata (downstream diagnostics)
    meta = Dict{Symbol, Any}(
        :residual_path         => residual_path,
        :residual_path_exists  => residual_path_exists,
        :n_tim_entries         => length(tim_entries),
        :n_residual_rows       => residuals === nothing ? 0 : length(residuals),
        :n_combined_entries    => combined_entries === nothing ? 0 : length(combined_entries),
        :analyze_white_noise   => analyze_white_noise,
        :time_start            => time_start,
        :time_finish           => time_finish,
        :stage_output_ok       => stage_output_ok,
        :stage_residual_read   => stage_residual_read,
        :stage_combine         => stage_combine,
        :stage_stats           => stage_stats,
        :stage_white_noise     => stage_white_noise,
        :iter_ok               => iter_ok,
        :iter_fail_stage       => iter_fail_stage,
        :residuals_stored      => save_residuals,
        :time_window_applied   => time_window_applied,
        :combined_mjd_min      => combined_mjd_min,
        :combined_mjd_max      => combined_mjd_max,
        :wn_performed          => wn_performed,
        :stats_error_msg       => stats_error_msg,
        :white_noise_error_msg => wn_error_msg,
    )

    return InternalIterationResult(
        output,
        save_residuals ? residuals : nothing,
        stats,
        white_noise_fit,
        meta,
    )
end