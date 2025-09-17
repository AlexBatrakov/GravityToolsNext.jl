using GravityToolsNext
using JLD2

basic_settings = TempoRunSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_GR",
    par_input = "DDSTG_GR_SCINT_RN84_test.par",
    par_output = "DDSTG_GR_SCINT_RN84_test_out.par",
    tim_input = "J1141-6545_pn_new.tim",
    tempo_version = Tempo2(),
    flags = "",
    nits = 4,
    gain = 1,
    override_params = [TP("DM", 116)],
    couple_f1_to_ddot = true,
    write_output = true,
    write_residuals = true,
    save_internal_iterations = true,
    save_residuals = false,
    white_noise_enabled = true,
    white_noise_scope = :all,                     # :final, :all
    work_mode = :jobdir,                           # :inplace, :jobdir
    job_name = "PRIOR_TEST",                       # nothing, "TEST"
    overwrite = :reuse,                            # :error, :reuse, :unique, :clean
    layout = :split,                               # :flat, :split
    temp_dir = nothing,                            # nothing, "nodes/node_1"
    link_tim = false,                              # false
    snapshot_par = true,                           # true
    cleanup_before_run = true,                     # true
    keep_tmp_on_success = true,                   # false
    keep_tmp_on_error = true,                      # true
    # timeout_s = nothing,
    write_manifest = false,                         # false
    # manifest_style = :json,                      # :json, :toml
    # verbosity = :info,                           # :silent, :warn, :info, :debug
    # with_timestamps = true,
)

base_task = BasicTempoTask(basic_settings)

prior_settings = PriorMarginalizationSettings(
        parameter = :DDOT,
        pin_mode = :fixed,
        prior = SampledPrior("J1141-6545_DDOT_prior.dat"),
        nodes = ClenshawCurtisNodes(6),
        likelihood_source = :chi2_fit,
        ref_strategy = :prior_median,
        representative = :prior_median,
        save_node_results = true,
        exec_options = PriorExecutionOptions(
        mode               = :chained,                   # :independent | :chained
        chain_direction    = :backward,                  # :forward | :backward
        chain_snapshot_par = true,
        # scheduler          = :serial,                   # :serial | :threads | :distributed
        # max_workers        = 0,
        # workdir_layout     = :per_node,                 # :per_node
        # node_dir_prefix    = "nodes/node_",
        # keep_node_dirs     = true,
        # on_error           = :collect,
        # dir_name_mode      = :index_only,               # :index_only | :with_value
        # index_pad          = 3,
        # value_sig          = 6
        )
)

task = PriorMarginalizedTempoTask(base_task, prior_settings)


# s      = task.settings
# thetas = GravityToolsNext._build_thetas(s)
# node_results = GravityToolsNext._run_nodes_serial!(task, thetas)

# jldsave("/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_GR/node_results.jld2"; node_results = node_results)

result = run_task(task)