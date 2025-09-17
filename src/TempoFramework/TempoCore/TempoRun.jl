#
# TempoRun.jl — low-level TEMPO/TEMPO2 runner
#
# Responsibilities of this file:
#  • materialize a reproducible workspace (jobdir layout, staging inputs);
#  • run the TEMPO/TEMPO2 engine with relative paths from the chosen CWD;
#  • parse stdout into iteration records and assemble filesystem artifacts;
#  • expose a small public API:
#       run_tempo_parsed(settings)::TempoRunOutput
#       cleanup_run!(::RunArtifacts/::TempoRunOutput, settings)
#  • DO NOT perform scientific analysis or post-run deletion of artifacts — that is task-level code.
#  • JSON helpers below are placeholders for lightweight manifest writing; subject to change.

# Assumes:
# - AbstractTempoVersion / Tempo / Tempo2 (AbstractTempo.jl)
# - tempo_cmd_path(v), tempo_env(v)
# - TempoParFile / read_par_file! / copy_par_file / update_par_file! / write_par_file!
# - parse_tempo_output (TempoOutput.jl)

# --------------------------------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------------------------------

# strip .par (case-insensitive) from a filename (no dirs)
_strip_par_ext(s::AbstractString) = replace(basename(s), r"\.[Pp][Aa][Rr]$" => "")

# Resolve path against work_dir if it’s not absolute
_resolve(work_dir::AbstractString, p::AbstractString) = isabspath(p) ? String(p) : joinpath(work_dir, p)

# Compute absolute paths for a *materialized* settings bundle
# Assumes: `paths.work_dir` is absolute; other paths are relative to it.
function _abs_paths(s::TempoRunSettings)
    @assert isabspath(s.paths.work_dir) "_abs_paths requires materialized settings (absolute work_dir)"
    wd      = s.paths.work_dir
    par_in  = joinpath(wd, s.paths.par_input)
    tim_in  = joinpath(wd, s.paths.tim_input)
    # par_output is a file name; route to output/ for split layout
    par_out = s.workspace.layout === :split ? joinpath(wd, "output", s.paths.par_output) : joinpath(wd, s.paths.par_output)
    run_cwd = s.workspace.temp_dir === nothing ? wd : joinpath(wd, s.workspace.temp_dir)
    return (work_dir = wd, par_in = par_in, tim_in = tim_in, par_out = par_out, run_cwd = run_cwd)
end


# Make a safe, unique directory name with numeric suffixes
function _unique_dir(base::AbstractString)
    if !isdir(base) && !isfile(base)
        return String(base)
    end
    for i in 1:999
        cand = string(base, "-", lpad(string(i), 3, '0'))
        if !isdir(cand) && !isfile(cand)
            return cand
        end
    end
    error("Failed to create unique directory name based on $(base)")
end

# Compute subdirectory under input/output to mirror (part of) temp_dir according to workspace.io_mirror
function _io_mirror_subdir(ws)
    # No mirroring for flat layout
    if ws.layout === :flat
        return ""
    end
    td = ws.temp_dir
    if td === nothing || isempty(String(td))
        return ""
    end
    # Normalize temp_dir and split into parts (platform-agnostic)
    nd = normpath(String(td))
    # replace backslashes with forward slashes to unify separators, then split
    nd = replace(nd, '\\' => '/')
    parts = split(nd, '/')
    parts = filter(p -> !isempty(p), parts)
    if isempty(parts)
        return ""
    end
    im = ws.io_mirror
    if im === :none
        return ""
    elseif im === :full
        return joinpath(parts...)
    elseif im isa Integer
        n = max(min(Int(im), length(parts)), 0)
        return n == 0 ? "" : joinpath(parts[1:n]...)
    elseif im isa Tuple{Symbol,Int}
        mode, k = im
        if mode === :depth_minus
            n = max(length(parts) - k, 0)
            return n == 0 ? "" : joinpath(parts[1:n]...)
        else
            error("Unsupported io_mirror tuple mode: $(mode)")
        end
    else
        error("Invalid io_mirror setting: $(im)")
    end
end

# --------------------------------------------------------------------------------------------------------------
# Materialization (job layout)
# --------------------------------------------------------------------------------------------------------------

"""
    MaterializedJob

Concrete, self-contained description of a single run workspace. It keeps the
*original* `TempoRunSettings` (never mutated) and records the realized
filesystem layout (job_root, input/output/tmp) and resolved paths to staged
inputs/outputs. Internal to TempoRun; not exported. Do not construct manually.

Fields
- `settings`       : the original settings used to materialize the job
- `job_root`       : absolute job root directory
- `input_dir`      : directory where input files live for this job
- `output_dir`     : directory where final outputs are written
- `tmp_dir`        : working CWD for the engine (may equal `job_root`)
- `run_cwd`        : actual CWD for the engine process
- `par_in_path`    : resolved path to the *input* par file used to create `_init.par`
- `tim_in_path`    : resolved path to the staged TIM inside the job
- `par_out_path`   : absolute path where the final `.par` is placed
- `created_paths`  : a list of paths created during materialization
- `notes`          : assorted flags (e.g., `:snapshot_par`, `:linked_tim`)
"""
struct MaterializedJob
    settings::TempoRunSettings
    job_root::String
    input_dir::String
    output_dir::String
    tmp_dir::String
    run_cwd::String
    par_in_path::String
    tim_in_path::String
    par_out_path::String
    created_paths::Vector{String}
    notes::Dict{Symbol,Any}
end

function _ensure_empty_dir(path::AbstractString)
    if isdir(path)
        for f in readdir(path)
            fp = joinpath(path, f)
            try
                isdir(fp) ? rm(fp; recursive=true) : rm(fp)
            catch err
                @warn "Failed to clean path" path=fp error=err
            end
        end
    else
        mkpath(path)
    end
end

"""
    cleanup_old_tempo_files(job::MaterializedJob) -> Nothing

Remove common TEMPO/TEMPO2 byproducts in `job.job_root` (and the resolved
`par_out_path`). This is a *pre-run* hygiene step and never touches user inputs
or previous persisted results outside the known byproduct set.
"""
function cleanup_old_tempo_files(job::MaterializedJob)
    dir = job.job_root
    base = _strip_par_ext(basename(job.par_in_path))

    files_to_delete = Set([
        "base_derivatives.dat",
        "derivatives.dat",
        "model_parameters.dat",
        "residuals.dat",
        "constraints.matrix",
        "design.matrix",
        "cov.matrix",
        "param.labels",
        "param.vals",
        "prefit.res",
        "postfit.res",
        "tnred.meta",
        "new.par",
        string(base, "_init.par"),
        string(base, "_init.out"),
        basename(job.par_out_path),
    ])

    residual_re = r"^residuals_\d+\.dat$"

    # Also try to remove the resolved par_output path (abs)
    try
        if isfile(job.par_out_path)
            rm(job.par_out_path)
        end
    catch err
        @warn "Failed to delete output par file" error=err
    end

    for file in readdir(dir)
        if file in files_to_delete || occursin(residual_re, file)
            full_path = joinpath(dir, file)
            if isfile(full_path)
                try
                    rm(full_path)
                catch err
                    @warn "Failed to delete tempo file" path=full_path error=err
                end
            end
        end
    end
end

"""
    materialize_job(s::TempoRunSettings) -> MaterializedJob

Create a reproducible run workspace according to `s.workspace` and stage inputs.
The returned `MaterializedJob` stores the original settings and the realized
filesystem layout and paths.

Semantics
- `work_mode`:
  * `:inplace` — use `paths.work_dir` as `job_root`.
  * `:jobdir`  — create or reuse a subdirectory inside `paths.work_dir`.
- `overwrite` (when the target jobdir already exists):
  * `:error`  — throw;
  * `:reuse`  — use as-is (no cleanup here; see `cleanup_before_run`);
  * `:unique` — create a suffixed sibling (…-001, …-002);
  * `:clean`  — purge contents and use the same directory.
- `layout`:
  * `:flat`  — `input_dir = output_dir = job_root`;
  * `:split` — create `input/`, `output/`, and `tmp/` (or `temp_dir`).
- `cleanup_before_run`:
  * if `tmp_dir != job_root` — empty only `tmp_dir`;
  * otherwise — remove typical TEMPO byproducts from `job_root`.

Staging
- TIM: symlink when `link_tim=true` (fallback to copy); no-op if already in place.
- PAR: copy snapshot when `snapshot_par=true` or the source is outside `job_root`.

Note: engine will be invoked from the chosen CWD with **relative** file names.
"""
function materialize_job(s::TempoRunSettings)
    # 1) Decide job_root
    ws = s.workspace
    work_root = s.paths.work_dir
    isdir(work_root) || mkpath(work_root)

    job_root = if ws.work_mode === :inplace
        work_root
    else
        ts = Dates.format(Dates.now(), "yyyymmdd-HHMMSS")
        r  = Random.randstring(5)
        auto_name = string("job-", ts, "-", r)
        base = joinpath(work_root, something(ws.job_name, auto_name))
        if isdir(base) || isfile(base)
            if ws.overwrite === :error
                error("Job directory already exists: $(base)")
            elseif ws.overwrite === :reuse
                # Use as-is; do not clean here. Cleanup (if any) is handled later by `cleanup_before_run`.
                base
            elseif ws.overwrite === :unique
                base_u = _unique_dir(base)
                mkpath(base_u)
                base_u
            elseif ws.overwrite === :clean
                _ensure_empty_dir(base)
                base
            else
                error("Unknown overwrite policy: $(ws.overwrite)")
            end
        else
            mkpath(base)
            base
        end
    end

    # 2) Layout
    if ws.layout === :flat
        input_dir = job_root
        output_dir = job_root
        tmp_dir = ws.temp_dir === nothing ? job_root : joinpath(job_root, ws.temp_dir)
        tmp_dir != job_root && mkpath(tmp_dir)
    elseif ws.layout === :split
        input_dir  = joinpath(job_root, "input")
        output_dir = joinpath(job_root, "output")
        tmp_dir    = ws.temp_dir === nothing ? joinpath(job_root, "tmp") : joinpath(job_root, ws.temp_dir)
        mkpath(input_dir)
        mkpath(output_dir)
        mkpath(tmp_dir)
    else
        error("Unknown layout: $(ws.layout)")
    end

    # Decide actual process CWD for the engine
    run_cwd = if ws.temp_dir === nothing
        tmp_dir != job_root ? tmp_dir : job_root
    else
        tmp_dir
    end
    run_cwd != job_root && mkpath(run_cwd)

    created = String[job_root]
    push!(created, input_dir)
    output_dir != job_root && push!(created, output_dir)
    tmp_dir    != job_root && push!(created, tmp_dir)

    # 3) Bring inputs
    par_src = _resolve(s.paths.work_dir, s.paths.par_input)
    tim_src = _resolve(s.paths.work_dir, s.paths.tim_input)
    isfile(par_src) || error("Input par file not found: $(par_src)")
    isfile(tim_src) || error("Input tim file not found: $(tim_src)")

    par_name = basename(s.paths.par_input)
    tim_name = basename(s.paths.tim_input)

    # Derive naming base from par_output: strip extension and a single trailing _out if present
    po_base = replace(_strip_par_ext(basename(s.paths.par_output)), r"_out$" => "")
    base_in_name = string(po_base, "_in.par")

    # Stage destinations under input_dir
    par_dst_default = joinpath(input_dir, par_name)  # only used when we don't snapshot/rename
    io_sub_for_in   = _io_mirror_subdir(ws)
    par_dst_custom  = isempty(io_sub_for_in) ? joinpath(input_dir, base_in_name) : joinpath(input_dir, io_sub_for_in, base_in_name)  # used when we snapshot
    tim_dst  = joinpath(input_dir, tim_name)

    notes = Dict{Symbol,Any}()

    # One-time archival of the original external par file at the job level.
    # We keep the original filename and store it flat in input/ (for :split) or job_root (for :flat).
    # Only do this if the source is outside job_root and the destination doesn't already exist.
    begin
        src_outside_job = !startswith(realpath(par_src), realpath(job_root))
        root_input_dir  = (ws.layout === :split) ? input_dir : job_root
        orig_dst        = joinpath(root_input_dir, par_name)
        if src_outside_job && !isfile(orig_dst)
            try
                cp(par_src, orig_dst; force=false)
                push!(created, orig_dst)
                notes[:original_par_archived] = true
            catch err
                @warn "Failed to archive original par into job root" src=par_src dst=orig_dst error=err
                notes[:original_par_archived] = false
            end
        end
    end

    # TIM: symlink or copy (skip if already in place)
    if abspath(tim_src) == abspath(tim_dst)
        # Nothing to do: TIM is already in the input_dir (inplace/flat case)
        notes[:linked_tim] = false
        notes[:tim_in_place] = true
    else
        if ws.link_tim
            try
                isfile(tim_dst) && rm(tim_dst)
                symlink(tim_src, tim_dst)
                notes[:linked_tim] = true
            catch err
                @warn "Symlink failed, fallback to copy" src=tim_src dst=tim_dst error=err
                cp(tim_src, tim_dst; force=true)
                notes[:linked_tim] = false
                push!(created, tim_dst)
            end
        else
            cp(tim_src, tim_dst; force=true)
            notes[:linked_tim] = false
            push!(created, tim_dst)
        end
    end

    # Early cleanup for tmp_dir to avoid wiping files we stage into run_cwd later
    if ws.cleanup_before_run && (tmp_dir != job_root)
        _ensure_empty_dir(tmp_dir)
    end

    # Ensure TIM is also available in run_cwd under the same basename
    tim_local = joinpath(run_cwd, tim_name)
    if !isfile(tim_local)
        try
            symlink(tim_dst, tim_local)
            notes[:tim_linked_in_cwd] = true
        catch err
            @warn "Symlink into run_cwd failed, fallback to copy" src=tim_dst dst=tim_local error=err
            cp(tim_dst, tim_local; force=true)
            notes[:tim_linked_in_cwd] = false
            push!(created, tim_local)
        end
    else
        notes[:tim_linked_in_cwd] = get(notes, :tim_linked_in_cwd, false)
    end

    # Ensure the raw input par is also available in run_cwd under <base>_in.par
    # This simplifies debugging (diff against <base>_init.par) regardless of snapshot_par policy.
    par_in_local = joinpath(run_cwd, base_in_name)
    if !isfile(par_in_local)
        try
            cp(par_src, par_in_local; force=true)
            push!(created, par_in_local)
            notes[:par_in_local] = true
        catch err
            @warn "Failed to stage raw par into run_cwd" src=par_src dst=par_in_local error=err
            notes[:par_in_local] = false
        end
    else
        notes[:par_in_local] = get(notes, :par_in_local, true)
    end

    # PAR: snapshot logic
    must_copy_par = ws.snapshot_par || !startswith(realpath(par_src), realpath(job_root))
    if must_copy_par
        # Archive the raw input par under a name tied to this run: <base>_in.par
        par_dst = par_dst_custom
        mkpath(dirname(par_dst))
        cp(par_src, par_dst; force=true)
        push!(created, par_dst)
        notes[:snapshot_par] = true
    else
        # Already inside job_root; use original file in place
        par_dst = par_src
        notes[:snapshot_par] = false
    end

    # 4) Resolve output par path within the job
    out_par_name = basename(s.paths.par_output)
    io_sub = _io_mirror_subdir(ws)
    par_file_out_path = isempty(io_sub) ? joinpath(output_dir, out_par_name) : joinpath(output_dir, io_sub, out_par_name)
    # Ensure parent exists for mirrored output
    mkpath(dirname(par_file_out_path))

    # 5) Optional cleanup before run (job-scoped)
    #    We create a provisional MaterializedJob just to reuse the cleanup helper,
    #    then replace it with the final one below.
    provisional_job = MaterializedJob(s, job_root, input_dir, output_dir, tmp_dir, run_cwd,
                                      par_dst, tim_dst, par_file_out_path, created, notes)
    if ws.cleanup_before_run
        if tmp_dir == job_root
            cleanup_old_tempo_files(provisional_job)
        end
    end

    job = provisional_job
    return job
end

# No-op overload for already materialized jobs

materialize_job(job::MaterializedJob) = job

# Internal invariant checks for a materialized job
# Throws an AssertionError if something is off. Users should never see these unless
# they circumvent the public API or mutate the job between materialization and run.
function _validate_job(job::MaterializedJob)
    @assert isabspath(job.job_root) "job_root must be absolute"
    @assert isabspath(job.input_dir) && isabspath(job.output_dir) && isabspath(job.tmp_dir) "job dirs must be absolute"
    @assert isabspath(job.run_cwd) "run_cwd must be absolute"
    @assert job.run_cwd == job.job_root || job.run_cwd == job.tmp_dir "run_cwd must be either job_root or tmp_dir"

    @assert isdir(job.job_root) "job_root does not exist: $(job.job_root)"
    @assert isdir(job.input_dir) "input_dir does not exist: $(job.input_dir)"
    @assert isdir(job.output_dir) "output_dir does not exist: $(job.output_dir)"
    @assert isdir(job.tmp_dir) || job.tmp_dir == job.job_root "tmp_dir does not exist: $(job.tmp_dir)"
    @assert isdir(job.run_cwd) "run_cwd does not exist: $(job.run_cwd)"

    @assert isfile(job.par_in_path) "input par not found: $(job.par_in_path)"
    @assert isfile(job.tim_in_path) "tim file not found: $(job.tim_in_path)"
    @assert isfile(joinpath(job.run_cwd, basename(job.tim_in_path))) "tim file not staged in run_cwd: $(joinpath(job.run_cwd, basename(job.tim_in_path)))"

    @assert isabspath(job.par_out_path) "par_out_path must be absolute"
    @assert isdir(dirname(job.par_out_path)) "parent of par_out_path must exist: $(dirname(job.par_out_path))"

    return nothing
end

# --- minimal JSON writer for Dict{String,Any}/Vector/String/Bool/Number ---
function _json_escape(s::AbstractString)
    buf = IOBuffer()
    for c in codeunits(s)
        if c == 0x22
            write(buf, "\\\"")
        elseif c == 0x5C
            write(buf, "\\\\")
        elseif c == 0x0A
            write(buf, "\\n")
        elseif c == 0x0D
            write(buf, "\\r")
        elseif c == 0x09
            write(buf, "\\t")
        else
            write(buf, Char(c))
        end
    end
    return String(take!(buf))
end

function _to_json(x)
    if x === nothing
        return "null"
    elseif x isa Bool
        return x ? "true" : "false"
    elseif x isa Integer || x isa AbstractFloat
        return string(x)
    elseif x isa AbstractString
        return string('"', _json_escape(x), '"')
    elseif x isa AbstractDict
        parts = String[]
        for (k,v) in x
            push!(parts, string('"', _json_escape(String(k)), '"', ":", _to_json(v)))
        end
        return string("{", join(parts, ","), "}")
    elseif x isa AbstractVector
        return string("[", join(_to_json.(x), ","), "]")
    else
        return string('"', _json_escape(string(x)), '"')
    end
end

function _write_json(path::AbstractString, obj)
    open(path, "w") do io
        print(io, _to_json(obj))
    end
end

# NOTE: legacy helper. Prefer calling `cleanup_run!(...)` from task-level code
# after you have finished reading artifacts. This function is no longer invoked
# automatically by `run_tempo_parsed`.
function finalize_job(job::MaterializedJob; success::Bool)
    options = job.settings.workspace
    # cleanup tmp based on policy, but never delete job_root implicitly
    if job.tmp_dir != job.job_root
        if success
            options.keep_tmp_on_success || (isdir(job.tmp_dir) && try rm(job.tmp_dir; recursive=true); catch; end)
        else
            options.keep_tmp_on_error || (isdir(job.tmp_dir) && try rm(job.tmp_dir; recursive=true); catch; end)
        end
    end
    return nothing
end


# --------------------------------------------------------------------------------------------------------------
# Run (raw)
# --------------------------------------------------------------------------------------------------------------

"""
    run_tempo_raw(job::MaterializedJob) -> (stdout::String, stderr::String, exit_code::Int)

Invoke TEMPO/TEMPO2 from the job's working directory (CWD) using **relative**
paths (`-f <init.par> <tim> -outpar <name.par>`). Returns captured stdout/stderr
and a coarse exit code: returns 0 on normal completion, 1 on Julia-side exception.
"""
function run_tempo_raw(job::MaterializedJob)
    settings = job.settings

    # Validate job invariants and existence of required files/dirs
    _validate_job(job)

    # Read input par; create an `_init.par` inside run_cwd; apply overrides
    run_cwd = job.run_cwd
    par_in = TempoParFile(basename(job.par_in_path), dirname(job.par_in_path))
    read_par_file!(par_in)

    # Derive base from par_output: strip extension and a single trailing _out
    po_base = replace(_strip_par_ext(basename(settings.paths.par_output)), r"_out$" => "")
    par_init_name = string(po_base, "_init.par")

    par_init = copy_par_file(par_in; new_name = par_init_name, new_dir = run_cwd)
    update_par_file!(par_init;
        override_params   = settings.modifiers.override_params,
        nits              = settings.engine.nits,
        gain              = settings.engine.gain,
        time_start        = settings.modifiers.time_start,
        time_finish       = settings.modifiers.time_finish,
        couple_f1_to_ddot = settings.modifiers.couple_f1_to_ddot,
    )
    write_par_file!(par_init)

    # Archive a copy of the exact engine input par into input_dir
    # This preserves what was actually fed to TEMPO even if tmp/ is later cleaned.
    begin
        par_init_src = joinpath(run_cwd, par_init.name)
        io_sub = _io_mirror_subdir(settings.workspace)
        par_init_dst = isempty(io_sub) ? joinpath(job.input_dir, par_init.name) : joinpath(job.input_dir, io_sub, par_init.name)
        mkpath(dirname(par_init_dst))
        if isfile(par_init_src)
            try
                cp(par_init_src, par_init_dst; force=true)
            catch e
                @warn "Failed to archive _init.par into input_dir" src=par_init_src dst=par_init_dst error=e
            end
        end
    end

    # Build argv using relative filenames (from CWD)
    tim_name = basename(job.tim_in_path)
    exe = tempo_cmd_path(settings.engine.tempo_version)
    # par_init_name already set above
    out_name = basename(job.par_out_path)
    args = String[
        exe,
        "-f", par_init_name,
        tim_name,
        "-outpar", out_name,
    ]

    # residuals flag(s)
    if settings.capture.write_residuals
        if settings.engine.tempo_version isa Tempo2
            push!(args, "-write_residuals")
        else
            push!(args, "-residuals")
        end
    end

    # user flags (split on whitespace)
    fl = strip(settings.engine.flags)
    if !isempty(fl)
        append!(args, split(fl))
    end

    cmd = Cmd(args)
    env = tempo_env(settings.engine.tempo_version)
    cmd = setenv(cmd, env; dir = run_cwd)

    output_io = IOBuffer(); stderr_io = IOBuffer()
    out_str = ""
    err_str = ""
    exit_code = 1
    try
        run(pipeline(cmd, stdout = output_io, stderr = stderr_io))
        out_str = String(take!(output_io))
        err_str = String(take!(stderr_io))
        exit_code = 0
    catch err
        # still flush buffers if any
        out_str = String(take!(output_io)) * "\nERROR: tempo process crashed with exception: $(err)"
        err_str = String(take!(stderr_io))
        exit_code = 1
    finally
        close(output_io); close(stderr_io)
    end

    # NOTE: we **copy** the engine-written outpar from CWD to `par_out_path` here.
    # The tmp/CWD remains a faithful snapshot of everything the engine produced.
    # Tasks and success/files_ok checks look at `artifacts.par_out_path`.
    outpar_path_cwd = joinpath(run_cwd, out_name)
    if isfile(outpar_path_cwd) && abspath(outpar_path_cwd) != abspath(job.par_out_path)
        try
            cp(outpar_path_cwd, job.par_out_path; force=true)
        catch e
            @warn "Failed to copy output .par" src=outpar_path_cwd dst=job.par_out_path error=e
        end
    end

    # Write stdout capture to CWD if requested; tasks can later prune/keep it per retention policy.
    if settings.capture.write_output
        out_log_name = string(po_base, ".out")
        out_path = joinpath(run_cwd, out_log_name)
        try
            open(out_path, "w") do io
                print(io, out_str)
            end
        catch e
            @warn "Failed to write .out file" path = out_path error = e
        end
        # Also copy .out into output_dir for persistence (tmp may be cleaned later)
        try
            io_sub = _io_mirror_subdir(settings.workspace)
            out_dst = isempty(io_sub) ? joinpath(job.output_dir, out_log_name) : joinpath(job.output_dir, io_sub, out_log_name)
            mkpath(dirname(out_dst))
            cp(out_path, out_dst; force=true)
        catch e
            @warn "Failed to copy .out into output_dir" src=out_path dst=out_dst error=e
        end
    end

    return out_str, err_str, exit_code
end

"""
    run_tempo_raw(settings::TempoRunSettings)

ERROR: This method is disabled. Materialize first: `job = materialize_job(settings)` and then call `run_tempo_raw(job)`. Or use `run_tempo_parsed(settings)`.
"""
function run_tempo_raw(::TempoRunSettings)
    throw(ArgumentError("run_tempo_raw expects a MaterializedJob. Call `job = materialize_job(settings)` first, then `run_tempo_raw(job)`, or use `run_tempo_parsed(settings)`."))
end


# --------------------------------------------------------------------------------------------------------------
# RunArtifacts
# --------------------------------------------------------------------------------------------------------------

"""
    RunArtifacts

Filesystem view of a concrete run. All paths are **absolute** and can be used
by task-level code to read, archive, or prune artifacts. Constructed by
`run_tempo_parsed`.

Fields
- `job_root`     : job root directory
- `run_cwd`      : actual engine CWD (`tmp/` for split layout by default)
- `input_dir`    : directory with staged inputs for the run
- `output_dir`   : directory where final outputs reside
- `tim_path`     : the staged TIM used by the run
- `par_out_path` : final `.par` path (may or may not exist on failure)
- `out_path`     : stdout capture file path (when enabled), otherwise `nothing`
- `residual_paths` : per-iteration residual file paths (`nothing` if absent)
"""
struct RunArtifacts
    job_root::String              # absolute root of this job
    run_cwd::String               # actual engine CWD (usually job_root or job_root/tmp)

    input_dir::String             # where staged inputs live (abs); for :flat == job_root
    output_dir::String            # where final outputs reside (abs); for :flat == job_root

    tim_path::String              # abs path to staged TIM (link or copy inside job)
    par_out_path::String          # abs path to final .par (may be missing on failure)

    out_path::Union{Nothing,String}         # abs path to .out if captured (write_output=true) else nothing
    residual_paths::Vector{Union{Nothing,String}}  # length == #iterations; abs path to residuals_i.dat or nothing
end

# --------------------------------------------------------------------------------------------------------------
# TempoRunOutput
# --------------------------------------------------------------------------------------------------------------

"""
    TempoRunOutput

Low-level result of invoking TEMPO/TEMPO2 for a single settings bundle. This is
*not* the scientific analysis; it only captures parsed iterations, file
artifacts and a coarse success/status for the engine run.

Fields
- `parsed`     : `Vector{InternalIterationOutput}` parsed from stdout
- `artifacts`  : `RunArtifacts` with absolute paths to all on-disk products
- `success`    : `true` when `engine_ok ∧ parse_ok ∧ files_ok`
- `status`     : `:ok | :engine_failed | :parse_failed | :files_missing`
- `n_iter`     : number of parsed iterations
- `exit_code`  : engine exit code (`0` ok, `1` exception)
- `stderr_tail`: last ~1000 characters of stderr, or `nothing`
- `started_at`, `finished_at`, `duration_s` : coarse timing info
"""
struct TempoRunOutput
    parsed::Vector{InternalIterationOutput}  # iterations parsed from TEMPO stdout
    artifacts::RunArtifacts                  # paths: job_root, run_cwd, input/output, tim, out.par, .out, residuals
    success::Bool                            # low-level run success (engine_ok ∧ parse_ok ∧ files_ok)

    # Diagnostics (lightweight but helpful)
    status::Symbol                           # :ok | :engine_failed | :parse_failed | :files_missing
    n_iter::Int                              # length(parsed)
    exit_code::Int                           # engine exit code
    stderr_tail::Union{Nothing,String}       # tail of stderr for quick debugging

    # Timings
    started_at::DateTime
    finished_at::DateTime
    duration_s::Float64
end

# --------------------------------------------------------------------------------------------------------------
# Run (parsed)
# --------------------------------------------------------------------------------------------------------------

"""
    run_tempo_parsed(settings::TempoRunSettings) :: TempoRunOutput

Materialize the workspace, run TEMPO/TEMPO2 with relative paths, parse stdout
into iteration records, and assemble `RunArtifacts`. No post-run deletion is
performed here; call `cleanup_run!` from task-level code **after** reading the
artifacts.
"""
function run_tempo_parsed(settings::TempoRunSettings)::TempoRunOutput
    job = materialize_job(settings)

    # timings
    started_at = Dates.now()
    out_str, err_str, exit_code = run_tempo_raw(job)
    finished_at = Dates.now()
    duration_s = Dates.value(finished_at - started_at) / 1e3

    # parse stdout
    parsed = parse_tempo_output(out_str, typeof(job.settings.engine.tempo_version))
    n_iter = length(parsed)

    # run_cwd (must match raw-run)
    run_cwd = job.run_cwd

    # residuals per iteration
    residual_paths = Vector{Union{Nothing,String}}(undef, n_iter)
    for i in 1:n_iter
        per_iter = joinpath(run_cwd, "residuals_$(i).dat")
        if isfile(per_iter)
            residual_paths[i] = per_iter
        elseif i == n_iter
            final_res = joinpath(run_cwd, "residuals.dat")
            residual_paths[i] = isfile(final_res) ? final_res : nothing
        else
            residual_paths[i] = nothing
        end
    end

    # .out (stdout capture), if enabled. Prefer persistent copy in output_dir if present.
    po_base  = replace(_strip_par_ext(basename(settings.paths.par_output)), r"_out$" => "")
    out_base = string(po_base, ".out")
    out_cwd  = joinpath(run_cwd, out_base)
    io_sub = _io_mirror_subdir(job.settings.workspace)
    out_copy = isempty(io_sub) ? joinpath(job.output_dir, out_base) : joinpath(job.output_dir, io_sub, out_base)
    if settings.capture.write_output
        if isfile(out_copy)
            out_path = out_copy
        elseif isfile(out_cwd)
            out_path = out_cwd
        else
            out_path = nothing
        end
    else
        out_path = nothing
    end

    artifacts = RunArtifacts(
        job.job_root,
        run_cwd,
        job.input_dir,
        job.output_dir,
        job.tim_in_path,
        job.par_out_path,
        out_path,
        residual_paths,
    )

    # success/status
    engine_ok = exit_code === 0
    parse_ok  = n_iter > 0 && !iserror(parsed[end].error)
    files_ok  = isfile(job.par_out_path) && (
        settings.capture.write_residuals ?
            (n_iter == 0 ? true : (residual_paths[end] !== nothing && isfile(something(residual_paths[end], "")))) :
            true
    )

    status = if !engine_ok
        :engine_failed
    elseif !parse_ok
        :parse_failed
    elseif !files_ok
        :files_missing
    else
        :ok
    end

    stderr_tail = if isempty(err_str)
        nothing
    else
        n = lastindex(err_str)
        firstidx = max(firstindex(err_str), n - 999)
        err_str[firstidx:n]
    end

    return TempoRunOutput(
        parsed,
        artifacts,
        engine_ok && parse_ok && files_ok,
        status,
        n_iter,
        exit_code,
        stderr_tail,
        started_at,
        finished_at,
        duration_s,
    )
end


# --------------------------------------------------------------------------------------------------------------
# Post-run cleanup API
# --------------------------------------------------------------------------------------------------------------

"""
    cleanup_run!(artifacts::RunArtifacts, settings::TempoRunSettings; success::Bool=true) -> Nothing

Remove the *working* directory (`run_cwd`) when it is a dedicated subdirectory
(e.g., `tmp/`). This is a **post-run** hygiene step controlled by
`keep_tmp_on_success` / `keep_tmp_on_error`. It never touches `input_dir` or
`output_dir` and is safe to call multiple times.
"""
function cleanup_run!(artifacts::RunArtifacts, settings::TempoRunSettings; success::Bool=true)
    if abspath(artifacts.run_cwd) != abspath(artifacts.job_root)
        if success
            settings.workspace.keep_tmp_on_success || (isdir(artifacts.run_cwd) && try rm(artifacts.run_cwd; recursive=true); catch; end)
        else
            settings.workspace.keep_tmp_on_error   || (isdir(artifacts.run_cwd) && try rm(artifacts.run_cwd; recursive=true); catch; end)
        end
    end
    return nothing
end

"""
    cleanup_run!(output::TempoRunOutput, settings::TempoRunSettings; success::Bool=output.success) -> Nothing

Ergonomic overload that accepts a `TempoRunOutput`. Delegates to the
`RunArtifacts`-based method. Override `success` if your task-level success
criterion differs from the engine-level one.
"""
function cleanup_run!(output::TempoRunOutput, settings::TempoRunSettings; success::Bool=output.success)
    return cleanup_run!(output.artifacts, settings; success=success)
end

