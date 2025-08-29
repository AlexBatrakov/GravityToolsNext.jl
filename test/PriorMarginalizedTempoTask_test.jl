using GravityToolsNext

basic_settings = BasicTempoSettings(
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
    white_noise_enabled = false,
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
        save_node_results = true
)

task = PriorMarginalizedTempoTask(base_task, prior_settings)
