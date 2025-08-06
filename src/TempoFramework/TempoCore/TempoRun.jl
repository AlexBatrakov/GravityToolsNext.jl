function cleanup_old_tempo_files(settings::BasicTempoSettings)
    dir = settings.files.work_dir
    base = settings.files.par_file_input[1:end-4]  # без .par

    # Жёстко заданные файлы
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
        base * "_init.par",
        base * "_init.out",
        settings.files.par_file_output,
    ])

    # Удаляем по шаблону: residuals_*.dat
    residual_pattern = r"^residuals_\d+\.dat$"

    for file in readdir(dir)
        if file in files_to_delete || occursin(residual_pattern, file)
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


function run_tempo_raw(settings::BasicTempoSettings)
    if !isdir(settings.files.work_dir)
        error("Working directory not found: $(settings.files.work_dir)")
    end
    cd(settings.files.work_dir)

    cleanup_old_tempo_files(settings)

    if !isfile(settings.files.par_file_input)
        error("Input par file not found: $(settings.files.par_file_input)")
    end

    par_file_input = TempoParFile(settings.files.par_file_input, settings.files.work_dir)
    read_par_file!(par_file_input)

    par_file_init = copy_par_file(par_file_input, suffix="init")
    update_par_file!(par_file_init,
                     override_params = settings.modifiers.override_params, 
                     nits = settings.options.nits, 
                     gain = settings.options.gain,
                     time_start = settings.modifiers.time_start,
                     time_finish = settings.modifiers.time_finish)
    write_par_file!(par_file_init)

    all_flags = settings.options.flags
    all_flags = all_flags * (settings.behavior.write_residuals ? " -write_residuals" : "")

    tempo_command = `$(get_tempo_command(settings.options.tempo_version)) -f $(par_file_init.name) $(settings.files.tim_file) -outpar $(settings.files.par_file_output) $([split(all_flags)...])`

    # Инициализация IOBuffer
    output_io = IOBuffer()
    stderr_io = IOBuffer()

    # **Предварительная инициализация переменных для избежания UndefVarError**
    output = ""
    stderr = ""

    try
        process = run(pipeline(tempo_command, stdout=output_io, stderr=stderr_io), wait=false)

        wait(process)

        output = String(take!(output_io))
        stderr = String(take!(stderr_io))
    catch err
        output = output = "ERROR: tempo process crashed with exception: $(err)"
    finally
        close(output_io)
        close(stderr_io)
    end

    if settings.behavior.write_output
        output_file_name = par_file_init.name[1:end-4] * ".out"
        write(output_file_name, output)
    end

    return output, stderr
end

function run_tempo_parsed(settings::BasicTempoSettings)
    output, stderr = run_tempo_raw(settings)

    all_parsed_outputs = parse_tempo_output(output, typeof(settings.options.tempo_version))

    return all_parsed_outputs
end