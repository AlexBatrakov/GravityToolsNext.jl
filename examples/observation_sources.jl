using GravityToolsNext

"""Prepare schedules and a TIM subset in an explicitly supplied output directory."""
function observation_sources_example(output_dir::AbstractString)
    mkpath(output_dir)
    schedule = regular_schedule("55825.123456789012345678901234";
        count=12, cadence_s=600, observatory="pks", frequency_mhz=1400,
        uncertainty_us=2, backend=repeat(["A", "B"], 6))
    schedule_path = write_tim_document(schedule.document, joinpath(output_dir, "schedule.tim"))
    source = tim_source(schedule_path)
    selected = select_toas(source; method=:stratified, count=6, seed=42)
    selected_path = write_tim_document(selected.document, joinpath(output_dir, "selected.tim"))
    println("Prepared ", length(selected.lineage), " observations at ", selected_path)
    return (; schedule, source, selected, schedule_path, selected_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Pass an output directory: observation_sources.jl OUTPUT_DIR")
    observation_sources_example(only(ARGS))
end
