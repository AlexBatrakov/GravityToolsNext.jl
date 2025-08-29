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
#   read_residual_file, combine_tim_and_residuals
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
- Returns early with the TEMPO error recorded if `output.error` is set.
- Reads residuals iff `residual_path` is a valid file.
- Combines TIM entries with residuals to compute residual statistics (guarded).
- Optionally performs white-noise fitting (`analyze_white_noise=true`) on combined in-fit data.
- Packs useful counts/flags/paths into `metadata`.

Keywords
- `save_residuals::Bool=false` — keep raw residual rows inside the result
- `analyze_white_noise::Bool=false` — run per-backend white-noise fit
- `time_start::Union{Nothing,Float64}=nothing` — lower MJD bound for combining
- `time_finish::Union{Nothing,Float64}=nothing` — upper MJD bound for combining
- `white_noise_kwargs...` — forwarded as keywords to `build_white_noise_fit` (optional future use)
"""
function build_internal_iteration_result(
    output::InternalIterationOutput,
    residual_path::Union{Nothing, String},
    tim_entries::Vector{TimTOAEntry};
    save_residuals::Bool = false,
    analyze_white_noise::Bool = false,
    time_start::Union{Nothing, Float64} = nothing,
    time_finish::Union{Nothing, Float64} = nothing,
    white_noise_kwargs...,
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
            # pass-through hook for future knobs (keeps API stable)
            white_noise_fit = build_white_noise_fit(combined_entries; white_noise_kwargs...)
        catch err
            @warn "White-noise fit failed" error=err
            white_noise_fit = nothing
        end
    end

    # Build metadata (useful for downstream logging/diagnostics)
    meta = Dict{Symbol, Any}(
        :residual_path        => residual_path,
        :n_tim_entries        => length(tim_entries),
        :n_residual_rows      => residuals === nothing ? 0 : length(residuals),
        :n_combined_entries   => combined_entries === nothing ? 0 : length(combined_entries),
        :analyze_white_noise  => analyze_white_noise,
        :time_start           => time_start,
        :time_finish          => time_finish,
    )

    return InternalIterationResult(
        output,
        save_residuals ? residuals : nothing,
        stats,
        white_noise_fit,
        meta,
    )
end