# TempoFramework/SingleTasks/BasicTempoTask.jl

"""
    BasicTempoTask(settings::TempoRunSettings)

Single TEMPO/TEMPO2 task with parsed iterations, optional residual statistics,
and optional white-noise analysis controlled by `settings.analysis`.
"""
struct BasicTempoTask <: SingleTempoTask
    settings::TempoRunSettings
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

"""
    _load_final_par_file(par_out_path::AbstractString) :: Union{TempoParFile,Nothing}

Attempt to read the final output `.par` file produced by TEMPO/TEMPO2.
Returns a `TempoParFile` on success or `nothing` when the file is absent or
cannot be parsed.
"""
function _load_final_par_file(par_out_path::AbstractString)::Union{TempoParFile,Nothing}
    if isfile(par_out_path)
        try
            pf = read_par_file(par_out_path)
            return pf
        catch err
            @warn "Failed to read final par file" path=par_out_path error=err
        end
    end
    return nothing
end

"""
    _wn_mask(niter::Int, scope::Symbol) :: BitVector

Construct a boolean mask that selects which iterations will undergo
white-noise analysis:
  * `:all`   — all iterations are selected
  * `:final` — only the last iteration is selected

Unknown scopes are treated defensively as `:final`.
"""
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
    settings = task.settings

    # 1) Run TEMPO and parse iteration-wise outputs (low-level)
    run_out = run_tempo_parsed(settings)             # ::TempoRunOutput
    artifacts = run_out.artifacts
    parsed_iter_outputs = run_out.parsed
    niter = length(parsed_iter_outputs)

    # 2) Read TIM entries
    tim_entries = read_tim_file_safe(artifacts.tim_path)
    tim_count = tim_entries === nothing ? missing : length(tim_entries)

    # 3) Decide which iterations get white-noise analysis
    wn_enabled = settings.analysis.enabled
    wn_mask_vec = wn_enabled ? _wn_mask(niter, settings.analysis.scope) : falses(niter)

    # 4) Build InternalIterationResult per iteration
    iter_results = InternalIterationResult[]
    sizehint!(iter_results, niter)

    for (i, out) in enumerate(parsed_iter_outputs)
        residual_path = artifacts.residual_paths[i]

        iter_res = build_internal_iteration_result(
            out, residual_path, tim_entries;
            time_start          = settings.modifiers.time_start,
            time_finish         = settings.modifiers.time_finish,
            save_residuals      = settings.retention.save_residuals,
            analyze_white_noise = wn_mask_vec[i]
        )
        push!(iter_results, iter_res)
    end

    # Indices for transparency: last iteration vs last successful with stats
    last_index = niter

    # 5) Optionally load the final par-file written by TEMPO
    par_file_final = _load_final_par_file(artifacts.par_out_path)

    # 6) Assemble the top-level result (convergence/params are computed inside)
    metadata = Dict{Symbol,Any}(
        :work_dir       => settings.paths.work_dir,
        :job_root       => artifacts.job_root,
        :run_cwd        => artifacts.run_cwd,
        :input_dir      => artifacts.input_dir,
        :output_dir     => artifacts.output_dir,
        :tim_path       => artifacts.tim_path,
        :par_out_path   => artifacts.par_out_path,
        :out_path       => artifacts.out_path,
        :tempo_ver      => typeof(settings.engine.tempo_version),
        :wn_enabled     => wn_enabled,
        :wn_scope       => settings.analysis.scope,
        :save_residuals => settings.retention.save_residuals,
        :flags          => settings.engine.flags,
        :nits           => settings.engine.nits,
        :gain           => settings.engine.gain,
        :time_start     => settings.modifiers.time_start,
        :time_finish    => settings.modifiers.time_finish,
        :io_mirror      => settings.workspace.io_mirror,
        # helpful run statistics
        :n_iter         => niter,
        :tim_entries_n  => tim_count,
        :run_success    => run_out.success,
        :run_status     => run_out.status,
        :exit_code      => run_out.exit_code,
        :started_at     => run_out.started_at,
        :finished_at    => run_out.finished_at,
        :duration_s     => run_out.duration_s,
        :stderr_tail    => run_out.stderr_tail,
        :status         => run_out.status,      # standardized duplicate of :run_status
        :success        => run_out.success,     # standardized duplicate of :run_success
        :last_index     => last_index,
    )

    result = build_general_tempo_result(iter_results; par_file_final, metadata)

    # Post-run hygiene: remove tmp/CWD per keep_tmp_on_* policies (safe after reading artifacts)
    cleanup_run!(run_out, settings)

    return result
end

"""
    task_workdir(task::BasicTempoTask) :: AbstractString

Return the working directory where this `BasicTempoTask` will execute.
"""
task_workdir(task::BasicTempoTask)::AbstractString = task.settings.paths.work_dir

"""
    task_with_overrides(task::BasicTempoTask, overrides; work_dir) :: BasicTempoTask

Clone the task with a new `work_dir` and a list of `TempoParameter` overrides
that will be upserted into the underlying `TempoRunSettings`.
"""
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


function task_copy_with(task::BasicTempoTask; kwargs...)
    s = task.settings
    s2 = copy_with(s; kwargs...)
    return BasicTempoTask(s2)
end

function task_stage_inputs!(task::BasicTempoTask, dest_dir::AbstractString)
    # Copy par/tim files into `dest_dir` if present. Accepts relative (to work_dir) or absolute paths.
    s = task.settings
    mkpath(dest_dir)
    for rel_path in (s.paths.par_input, s.paths.tim_input)
        rel_path === nothing && continue
        src_path = isabspath(rel_path) ? rel_path : joinpath(s.paths.work_dir, rel_path)
        if isfile(src_path)
            dest_path = joinpath(dest_dir, basename(rel_path))
            if abspath(src_path) != abspath(dest_path)
                cp(src_path, dest_path; force=true)
            end
        else
            @warn "Input file to stage not found" path=src_path
        end
    end
    return nothing
end

#
# BasicTempoTask-specific derivation:
# Take the stem of base `par_output`, drop one trailing `_out` if present,
# then append `_" * node_tag * "_out.par"`.
function task_derive_par_output(t::BasicTempoTask, node_tag::AbstractString)
    base = basename(t.settings.paths.par_output)
    stem, _ = splitext(base)
    # drop one _out suffix if present
    if endswith(stem, "_out")
        stem = first(stem, lastindex(stem) - 4)
    end
    return string(stem, "_", String(node_tag), "_out.par")
end