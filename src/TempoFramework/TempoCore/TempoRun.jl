# TempoRun.jl
# Running TEMPO/TEMPO2 and capturing its outputs.

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

# --------------------------------------------------------------------------------------------------------------
# Cleanup
# --------------------------------------------------------------------------------------------------------------

"""
    cleanup_old_tempo_files(settings::BasicTempoSettings)

Remove typical TEMPO/TEMPO2 byproducts from `settings.files.work_dir`.
Also removes `*_init.*` and the configured output `.par`.
"""
function cleanup_old_tempo_files(settings::BasicTempoSettings)
    dir  = settings.files.work_dir
    # base name derived from the *input par file name* (not path)
    base = _strip_par_ext(settings.files.par_file_input)

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
        basename(settings.files.par_file_output),
    ])

    residual_re = r"^residuals_\d+\.dat$"

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

# --------------------------------------------------------------------------------------------------------------
# Run (raw)
# --------------------------------------------------------------------------------------------------------------

"""
    run_tempo_raw(settings::BasicTempoSettings) -> (stdout::String, stderr::String)

Prepare working directory, create an `_init.par` with overrides, run TEMPO/TEMPO2,
capture stdout/stderr, and (optionally) write `.out`.
All non-absolute paths in settings are interpreted relative to `work_dir`.
"""
function run_tempo_raw(settings::BasicTempoSettings)
    work_dir = settings.files.work_dir
    isdir(work_dir) || error("Working directory not found: $work_dir")

    # resolve input/output paths w.r.t. work_dir
    par_in_path  = _resolve(work_dir, settings.files.par_file_input)
    tim_path     = _resolve(work_dir, settings.files.tim_file)
    out_par_name = basename(settings.files.par_file_output)  # we will write it *in* work_dir

    isfile(par_in_path) || error("Input par file not found: $par_in_path")
    isfile(tim_path)    || error("Tim file not found: $tim_path")

    # clean old artifacts *before* we start
    cleanup_old_tempo_files(settings)

    # read input par from its real location; we'll copy into work_dir
    par_in = TempoParFile(basename(par_in_path), dirname(par_in_path))
    read_par_file!(par_in)

    # make _init.par in work_dir, apply overrides
    par_init = copy_par_file(par_in; suffix="init", new_dir=work_dir)  # keep your copy semantics
    update_par_file!(par_init;
        override_params = settings.modifiers.override_params,
        nits  = settings.options.nits,
        gain  = settings.options.gain,
        time_start  = settings.modifiers.time_start,
        time_finish = settings.modifiers.time_finish,
    )
    write_par_file!(par_init)

    # build argv
    exe  = tempo_cmd_path(settings.options.tempo_version)
    args = String[
        exe,
        "-f", par_init.name,         # the _init.par in work_dir
        settings.files.tim_file,     # run with path relative to work_dir (we cd below)
        "-outpar", out_par_name,
    ]

    # residuals flag(s)
    if settings.behavior.write_residuals
        if settings.options.tempo_version isa Tempo2
            push!(args, "-write_residuals") # your Tempo2 feature: residuals_N.dat per internal iteration
        else
            push!(args, "-residuals")       # original Tempo: final residuals
        end
    end

    # user flags (split on whitespace)
    fl = strip(settings.options.flags)
    if !isempty(fl)
        append!(args, split(fl))
    end

    # NB: Cmd(::Vector{String}) — без kwargs
    cmd = Cmd(args)

    output_io = IOBuffer()
    stderr_io = IOBuffer()
    out    = ""
    errtxt = ""

    env = tempo_env(settings.options.tempo_version)  # Dict("TEMPO2" => "/...") или пустой Dict
    cmd = setenv(cmd, env; dir = settings.files.work_dir)

    try
        run(pipeline(cmd, stdout=output_io, stderr=stderr_io))
        out    = String(take!(output_io))
        errtxt = String(take!(stderr_io))
    catch err
        out    = "ERROR: tempo process crashed with exception: $(err)"
        errtxt = String(take!(stderr_io))
    finally
        close(output_io); close(stderr_io)
    end

    if settings.behavior.write_output
        out_name = string(_strip_par_ext(par_init.name), ".out")
        out_path = joinpath(settings.files.work_dir, out_name)
        try
            write(out_path, out)
        catch e
            @warn "Failed to write .out file" path=out_path error=e
        end
    end

    return out, errtxt
end

# --------------------------------------------------------------------------------------------------------------
# Run (parsed)
# --------------------------------------------------------------------------------------------------------------

"""
    run_tempo_parsed(settings::BasicTempoSettings) -> Vector{InternalIterationOutput}

Run TEMPO/TEMPO2 and parse the output into iteration-wise structures.
"""
function run_tempo_parsed(settings::BasicTempoSettings)
    output, stderr = run_tempo_raw(settings)
    all_parsed = parse_tempo_output(output, typeof(settings.options.tempo_version))
    return all_parsed
end