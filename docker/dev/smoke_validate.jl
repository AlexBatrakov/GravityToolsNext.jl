using GravityToolsNext

smoke_work_dir = get(ENV, "GTN_SMOKE_WORK_DIR", nothing)
smoke_par = get(ENV, "GTN_SMOKE_PAR", nothing)
smoke_tim = get(ENV, "GTN_SMOKE_TIM", nothing)
smoke_out = get(ENV, "GTN_SMOKE_OUT", nothing)
smoke_job_name = get(ENV, "GTN_SMOKE_JOB_NAME", nothing)

smoke_work_dir === nothing && error("GTN_SMOKE_WORK_DIR is not set")
smoke_par === nothing && error("GTN_SMOKE_PAR is not set")
smoke_tim === nothing && error("GTN_SMOKE_TIM is not set")
smoke_out === nothing && error("GTN_SMOKE_OUT is not set")
smoke_job_name === nothing && error("GTN_SMOKE_JOB_NAME is not set")

tempo = Tempo()
tempo2 = Tempo2()

validate(tempo)
validate(tempo2)

tempo_ddstg = joinpath(tempo_data_dir(tempo), "data_ddstg")
tempo2_ddstg = joinpath(tempo_data_dir(tempo2), "data_ddstg")

isdir(tempo_ddstg) || error("TEMPO DDSTG data directory is missing: $(tempo_ddstg)")
ispath(tempo2_ddstg) || error("TEMPO2 DDSTG link is missing: $(tempo2_ddstg)")
realpath(tempo2_ddstg) == realpath(tempo_ddstg) || error("TEMPO2 DDSTG link does not point at TEMPO data_ddstg")

println("Tempo validated:  ", tempo_cmd_path(tempo), " | ", tempo_data_dir(tempo))
println("Tempo2 validated: ", tempo_cmd_path(tempo2), " | ", tempo_data_dir(tempo2))
println("DDSTG link ok:    ", tempo2_ddstg, " -> ", realpath(tempo2_ddstg))

settings = TempoRunSettings(
    work_dir = smoke_work_dir,
    par_input = smoke_par,
    par_output = smoke_out,
    tim_input = smoke_tim,
    tempo_version = tempo2,
    flags = "",
    nits = 1,
    gain = 1.0,
    write_output = true,
    write_residuals = true,
    save_internal_iterations = false,
    save_residuals = false,
    white_noise_enabled = false,
    work_mode = :jobdir,
    job_name = smoke_job_name,
    overwrite = :clean,
    layout = :split,
    keep_tmp_on_success = false,
    keep_tmp_on_error = true,
    write_manifest = false,
)

validate(settings)

result = run_tempo_parsed(settings)

result.success || error("Smoke run failed with status=$(result.status), exit_code=$(result.exit_code), stderr_tail=$(repr(result.stderr_tail))")
isfile(result.artifacts.par_out_path) || error("Expected output par file is missing: $(result.artifacts.par_out_path)")

if result.artifacts.out_path === nothing || !isfile(result.artifacts.out_path)
    error("Expected captured stdout file is missing")
end

has_residual = any(path -> path !== nothing && isfile(path), result.artifacts.residual_paths)
has_residual || error("Expected at least one residuals file from the smoke run")

println("Smoke run ok:")
println("  status     = ", result.status)
println("  iterations = ", result.n_iter)
println("  duration_s = ", round(result.duration_s; digits = 3))
println("  job_root   = ", result.artifacts.job_root)
println("  par_out    = ", result.artifacts.par_out_path)
