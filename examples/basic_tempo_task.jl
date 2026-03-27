using GravityToolsNext

fixture_root = normpath(joinpath(@__DIR__, "..", "test", "data", "ddstg_smoke"))

settings = TempoRunSettings(
    work_dir = fixture_root,
    par_input = "j1141-6545_ddstg_gr_pdfb1_white_noise.par",
    tim_input = "j1141-6545_ddstg_gr_pdfb1_white_noise.tim",
    par_output = "j1141-6545_ddstg_gr_pdfb1_white_noise_out.par",
    tempo_version = Tempo2(fixture_root),  # placeholder local path for a construction-only example
    nits = 1,
    gain = 1.0,
    write_output = true,
    write_residuals = true,
    work_mode = :jobdir,
    job_name = "example_basic",
    overwrite = :reuse,
    layout = :split,
)

task = BasicTempoTask(settings)

println(task)

# Replace `Tempo2(fixture_root)` with your real TEMPO2 runtime directory, then uncomment:
# validate(task.settings)
# result = run_task(task)
# @show result.success result.status
