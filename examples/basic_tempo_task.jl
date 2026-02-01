using GravityToolsNext
using JLD2

basic_settings = TempoRunSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_GR",
    par_input = "DDSTG_GR_SCINT_RN84_test.par",
    par_output = "DDSTG_GR_SCINT_RN84_test_out.par",
    tim_input = "J1141-6545_pn_new.tim",
    tempo_version = Tempo2(),
    flags = "",
    nits = 1,
    gain = 1,
    override_params = [TP("DM", 116)],
    # time_start = nothing,
    # time_finish = nothing,
    write_output = true,
    write_residuals = true,
    save_internal_iterations = true,
    save_residuals = false,
    white_noise_enabled = true,
    white_noise_scope = :all,                     # :final, :all
    work_mode = :jobdir,                          # :inplace, :jobdir
    job_name = "TEST",                            # nothing, "TEST"
    overwrite = :reuse,                           # :error, :reuse, :unique, :clean
    layout = :split,                               # :flat, :split
    temp_dir = nothing,                            # nothing, "nodes/node_1"
    link_tim = false,                              # false
    snapshot_par = true,                           # true
    cleanup_before_run = true,                     # true
    keep_tmp_on_success = false,                   # false
    keep_tmp_on_error = true,                      # true
    # timeout_s = nothing,
    write_manifest = false,                         # false
    # manifest_style = :json,                      # :json, :toml
    # verbosity = :info,                           # :silent, :warn, :info, :debug
    # with_timestamps = true,
)

# basic_settings_new = copy(basic_settings; gain = 0.5, nits = 3)

# out = run_tempo_parsed(basic_settings)

basic_task = BasicTempoTask(basic_settings)

basic_task_result = run_task(basic_task)

# jldsave("basic_task_result.jld2"; basic_task_result = basic_task_result)

basic_task_result = load("basic_task_result.jld2", "basic_task_result")

# parsed_iter_outputs = run_tempo_parsed(settings)

# jldsave("parsed_iter_outputs.jld2"; parsed_iter_outputs = parsed_iter_outputs)

# parsed_iter_outputs = load("parsed_iter_outputs.jld2", "parsed_iter_outputs")

#------------------------------------------------------------------
# manual node rules with BasicTempoTask through TempoRunSettings
#------------------------------------------------------------------

basic_settings_1 = copy_with(
    basic_settings; 
    override_params = [TP("DDOT", 5e-19)],
    par_output = "DDSTG_GR_SCINT_RN84_test_node_1_out.par",
    overwrite = :reuse,
    layout = :split,
    temp_dir = "nodes/node_1",
    io_mirror = (:depth_minus, 2),
    keep_tmp_on_success = true
    )

basic_task_1 = BasicTempoTask(basic_settings_1)

basic_task_result_1 = run_task(basic_task_1)

#------------------------------------------------------------------

basic_settings_2 = copy_with(
    basic_settings; 
    override_params = [TP("DDOT", 5e-19)],
    par_output = "DDSTG_GR_SCINT_RN84_test_node_2_out.par",
    overwrite = :reuse,
    layout = :split,
    temp_dir = "nodes/node_2",
    io_mirror = (:depth_minus, 2),
    keep_tmp_on_success = true
    )

basic_task_2 = BasicTempoTask(basic_settings_2)

basic_task_result_2 = run_task(basic_task_2)





#------------------------------------------------------------------
# manual node rules with IterativeTempoTask through TempoRunSettings
#------------------------------------------------------------------


basic_settings_1_1 = copy_with(
    basic_settings; 
    override_params = [TP("DDOT", 4e-19)],
    par_output = "DDSTG_GR_SCINT_RN84_test_node_1_iter_1_out.par",
    overwrite = :reuse,
    layout = :split,
    temp_dir = "nodes/node_1/iters/iter_1",
    io_mirror = (:depth_minus, 2),
    keep_tmp_on_success = true
    )

basic_task_1_1 = BasicTempoTask(basic_settings_1_1)

basic_task_result_1_1 = run_task(basic_task_1_1)

#------------------------------------------------------------------

basic_settings_1_2 = copy_with(
    basic_settings; 
    override_params = [TP("DDOT", 4e-19)],
    par_input  = "TEST/output/nodes/node_1/DDSTG_GR_SCINT_RN84_test_node_1_iter_1_out.par",
    par_output = "DDSTG_GR_SCINT_RN84_test_node_1_iter_2_out.par",
    overwrite = :reuse,
    layout = :split,
    temp_dir = "nodes/node_1/iters/iter_2",
    io_mirror = (:depth_minus, 2),
    keep_tmp_on_success = true
    )

basic_task_1_2 = BasicTempoTask(basic_settings_1_2)

basic_task_result_1_2 = run_task(basic_task_1_2)

#------------------------------------------------------------------


