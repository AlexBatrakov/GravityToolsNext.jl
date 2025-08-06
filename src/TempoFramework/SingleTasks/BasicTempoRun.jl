struct BasicTempoRun <: SingleTempoTask
    settings::BasicTempoSettings
end

function Base.show(io::IO, ::MIME"text/plain", run::BasicTempoRun)
    println(io, "BasicTempoRun")
    println(io, "  settings:")
    show(IOContext(io, :indent => 4), MIME"text/plain"(), run.settings)
end

function run_task(task::BasicTempoRun)::GeneralTempoResult
    settings = task.settings
    all_parsed_outputs = run_tempo_parsed(settings)

    selected_index = length(all_parsed_outputs)

    tim_file_entries = read_tim_file(joinpath(settings.files.work_dir, settings.files.tim_file))

    all_parsed_results = InternalIterationResult[]

    for (i, parsed_output) in enumerate(all_parsed_outputs)
        residual_path = settings.behavior.write_residuals ? joinpath(settings.files.work_dir, "residuals_$i.dat") : nothing
        iter_result = build_internal_iteration_result(parsed_output, residual_path, tim_file_entries; time_start = settings.modifiers.time_start, time_finish = settings.modifiers.time_finish, save_residuals = true)#settings.behavior.save_residuals)
        push!(all_parsed_results, iter_result)
    end

    convergence_info = build_convergence_info(all_parsed_results)



    final_iter = internal_iteration_results[selected_index]
    final_output = final_iter.output.output
    final_error  = final_iter.output.error

    return GeneralTempoResult(
        all_parsed_outputs,
        selected_index,
        final_output,
        final_error,
        Dict{Symbol, Any}(),  # convergence_info
        nothing,              # или settings.files.par_file_output, если хочешь читать обратно
        Dict{Symbol, NamedTuple{(:value, :uncertainty), Tuple{Float64, Float64}}}(),  # позже
        Dict(:all => final_iter.stats),  # можно будет развернуть
        nothing,  # white_noise_fit
        nothing, nothing,
        Dict{Symbol, Any}(),
    )

end