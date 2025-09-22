using Distributed
using GravityToolsNext

addprocs(8)

@everywhere using GravityToolsNext


basic_settings = TempoRunSettings(
    work_dir = "/Users/abatrakov/Documents/Work/PhD/projects/J1141-6545/final_thesis/DDSTG_GR",
    par_input = "DDSTG_GR_SCINT_RN84_test.par",
    par_output = "DDSTG_GR_SCINT_RN84_test_out.par",
    tim_input = "J1141-6545_pn_new.tim",
    tempo_version = Tempo2(),
    flags = "",
    nits = 3,
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
    # work_mode = :jobdir,                          # :inplace, :jobdir
    # job_name = "TEST",                            # nothing, "TEST"
    # overwrite = :reuse,                           # :error, :reuse, :unique, :clean
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


basic_task = BasicTempoTask(basic_settings)

# result = run_task(basic_task)


grid_task = Adaptive2DGridTask(
    base_task = basic_task,
    x = GridAxis(:I,  min=74.0, max=79.0, N=2, rule = LinRule()),      # name, min, max, initial nodes
    y = GridAxis(:M2, min=0.99, max=1.01, N=2, rule = LinRule()),     # name, min, max, initial nodes
    ref_settings = RefinementSettings(
        params_to_save = (:chi2_fit, :wrms_fit, :wrms_tn_fit, :ad_white_fit, :pre_post_final, :F0, :F1, :F2, :PB, :T0, :A1, :OM, :ECC, :PBDOT, :XDOT, :OMDOT, :M2, :MTOT, :GAMMA, :I, :IDOT),
        desired_refinement_level = 1,
        parallel = true,
        LocalMinimaUnit(:chi2_fit, from_min=true, max=100.0)
    ),
    opts = GridWorkspaceOptions(
        grid_root = "GRID_TEST"
    ),
)

grid_result = run_task(grid_task)