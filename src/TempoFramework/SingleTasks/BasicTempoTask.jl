# TempoFramework/SingleTasks/BasicTempoTask.jl

"""
    BasicTempoTask(settings::BasicTempoSettings)

Single TEMPO/TEMPO2 task with parsed iterations, optional residual statistics,
and optional white-noise analysis controlled by `settings.analysis`.
"""
struct BasicTempoTask <: SingleTempoTask
    settings::BasicTempoSettings
end

function Base.show(io::IO, ::MIME"text/plain", task::BasicTempoTask)
    indent = get(io, :indent, 0)
    pad    = repeat(" ", indent)
    spad   = repeat(" ", indent + 2)
    iop    = IOContext(io, :indent => indent + 4)
    println(io, pad,  "BasicTempoTask")
    println(io, spad, "settings:")
    show(iop, MIME"text/plain"(), task.settings)
end


# --- helpers ----------------------------------------------------------------------------------------------------

# Resolve `.tim` path against work_dir if the path is relative
_resolve_tim_path(s::BasicTempoSettings) =
    isabspath(s.files.tim_file) ? s.files.tim_file : joinpath(s.files.work_dir, s.files.tim_file)

_resolve_par_input_path(s::BasicTempoSettings) =
    isabspath(s.files.par_file_input) ? s.files.par_file_input :
                                        joinpath(s.files.work_dir, s.files.par_file_input)

# Pick residual file path for iteration `i`; prefer per-iteration files if present.
# Falls back to `residuals.dat` only for the final iteration.
function _pick_residual_path(s::BasicTempoSettings, i::Int, niter::Int)
    if s.behavior.write_residuals
        cand = joinpath(s.files.work_dir, "residuals_$(i).dat")
        if isfile(cand)
            return cand
        end
    end
    # fallback: final residuals only for the last iteration
    cand_final = joinpath(s.files.work_dir, "residuals.dat")
    if isfile(cand_final) && i == niter
        return cand_final
    end
    return nothing
end

# Try to read the final out-par into a TempoParFile (optional)
function _load_final_par_file(s::BasicTempoSettings)
    out_name = basename(s.files.par_file_output)
    out_path = joinpath(s.files.work_dir, out_name)
    if isfile(out_path)
        pf = TempoParFile(out_name, s.files.work_dir)
        try
            read_par_file!(pf)
            return pf
        catch err
            @warn "Failed to read final par file" path=out_path error=err
        end
    end
    return nothing
end

# Build a boolean mask selecting iterations for white-noise analysis
#   scope = :all   -> all iterations
#   scope = :final -> final iteration only
function _wn_mask(niter::Int, scope::Symbol)
    if scope === :all
        return trues(niter)
    elseif scope === :final
        m = falses(niter)
        niter >= 1 && (m[end] = true)
        return m
    else
        # be defensive: treat unknown scopes as :final
        m = falses(niter)
        niter >= 1 && (m[end] = true)
        return m
    end
end

"""
    run_task(task::BasicTempoTask) -> GeneralTempoResult

Run TEMPO/TEMPO2 using `task.settings`, parse iteration-wise outputs, attach
residual/TIM-derived statistics, optionally perform white-noise fitting, and
return an aggregated `GeneralTempoResult`.
"""
function run_task(task::BasicTempoTask)::GeneralTempoResult
    s = task.settings

    # 1) Run TEMPO and parse iteration-wise outputs
    parsed_iter_outputs = run_tempo_parsed(s)  # Vector{InternalIterationOutput}
    niter = length(parsed_iter_outputs)
    niter == 0 && error("run_tempo_parsed returned no iterations")

    # 2) Read TIM entries
    tim_path = _resolve_tim_path(s)
    tim_entries = TimTOAEntry[]
    try
        tim_entries = read_tim_file(tim_path)
    catch err
        @warn "Failed to read TIM file; statistics will be limited" path=tim_path error=err
    end

    # 3) Decide which iterations get white-noise analysis
    wn_enabled = s.analysis.enabled
    wn_mask_vec = wn_enabled ? _wn_mask(niter, s.analysis.scope) : falses(niter)

    # 4) Build InternalIterationResult per iteration
    results = InternalIterationResult[]
    sizehint!(results, niter)

    for (i, out) in enumerate(parsed_iter_outputs)
        residual_path = _pick_residual_path(s, i, niter)

        iter_res = build_internal_iteration_result(
            out, residual_path, tim_entries;
            time_start          = s.modifiers.time_start,
            time_finish         = s.modifiers.time_finish,
            save_residuals      = s.behavior.save_residuals,
            analyze_white_noise = wn_mask_vec[i]
        )
        push!(results, iter_res)
    end

    # 5) Optionally load the final par-file written by TEMPO
    par_file_final = _load_final_par_file(s)

    # 6) Assemble the top-level result (convergence/params are computed inside)
    meta = Dict{Symbol,Any}(
        :work_dir       => s.files.work_dir,
        :tim_path       => tim_path,
        :par_out        => basename(s.files.par_file_output),
        :tempo_ver      => typeof(s.options.tempo_version),
        :flags          => s.options.flags,
        :nits           => s.options.nits,
        :gain           => s.options.gain,
        :time_start     => s.modifiers.time_start,
        :time_finish    => s.modifiers.time_finish,
        :wn_enabled     => wn_enabled,
        :wn_scope       => s.analysis.scope,
        :save_residuals => s.behavior.save_residuals,
    )

    return GeneralTempoResult(results; par_file_final=par_file_final, metadata=meta)
end

# ——— interface hook: where this task runs
task_workdir(task::BasicTempoTask)::AbstractString = task.settings.files.work_dir

# ——— interface hook: clone task with overrides and a new work_dir
function task_with_overrides(task::BasicTempoTask,
                             overrides::AbstractVector{TempoParameter};
                             work_dir::AbstractString)::BasicTempoTask
    s = task.settings
    s2 = copy_with(s;
        work_dir               = work_dir,
        override_params_upsert = collect(overrides),
    )
    return BasicTempoTask(s2)
end

"""
    task_stage_inputs!(task::BasicTempoTask, node_dir::AbstractString) -> Nothing

Copy the input `.par` and `.tim` files into `node_dir` so the task can run
with name-only file references.
"""
function task_stage_inputs!(task::BasicTempoTask, node_dir::AbstractString)
    s = task.settings
    mkpath(node_dir)

    # par input
    src_par = _resolve_par_input_path(s)
    dst_par = joinpath(node_dir, basename(s.files.par_file_input))
    cp(src_par, dst_par; force=true)

    # tim file
    src_tim = _resolve_tim_path(s)
    dst_tim = joinpath(node_dir, basename(s.files.tim_file))
    cp(src_tim, dst_tim; force=true)

    return nothing
end
