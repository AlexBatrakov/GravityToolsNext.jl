using GravityToolsNext
using JLD2

basic_settings = settings = BasicTempoSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_GR",
    par_file_input = "DDSTG_GR_SCINT_RN84_test.par",
    par_file_output = "DDSTG_GR_SCINT_RN84_test_out.par",
    tim_file = "J1141-6545_pn_new.tim",
    tempo_version = Tempo2(),
    flags = "",
    nits = 3,
    gain = 1,
    override_params = [TP("DM", 116)],
    write_output = true,
    write_residuals = true,
    save_internal_iterations = true,
    save_residuals = false,
    white_noise_enabled = true,
    white_noise_scope = :all
)

# basic_settings_new = copy(basic_settings; gain = 0.5, nits = 3)

# out = run_tempo_parsed(basic_settings)

task = BasicTempoRun(basic_settings)

parsed_iter_outputs = run_tempo_parsed(settings)

# jldsave("parsed_iter_outputs.jld2"; parsed_iter_outputs = parsed_iter_outputs)

parsed_iter_outputs = load("parsed_iter_outputs.jld2", "parsed_iter_outputs")
