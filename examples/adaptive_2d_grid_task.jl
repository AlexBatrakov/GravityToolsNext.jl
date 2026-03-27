using GravityToolsNext
using Distributions: Normal

fixture_root = normpath(joinpath(@__DIR__, "..", "test", "data", "ddstg_smoke"))

base_settings = TempoRunSettings(
    work_dir = fixture_root,
    par_input = "j1141-6545_ddstg_gr_pdfb1_white_noise.par",
    tim_input = "j1141-6545_ddstg_gr_pdfb1_white_noise.tim",
    par_output = "j1141-6545_ddstg_gr_pdfb1_white_noise_out.par",
    tempo_version = Tempo2(fixture_root),  # placeholder local path for a construction-only example
    nits = 1,
    gain = 1.0,
    work_mode = :jobdir,
    job_name = "example_grid",
    overwrite = :reuse,
    layout = :split,
)

base_task = BasicTempoTask(base_settings)

prior_task = PriorMarginalizedTempoTask(
    base_task,
    PriorMarginalizationSettings(
        parameter = :DDOT,
        pin_mode = :fixed,
        prior = AnalyticPrior(Normal(0.0, 1.0)),
        nodes = ClenshawCurtisNodes(4),
        likelihood_source = :chi2_fit,
        representative = :prior_median,
    ),
)

x = GridAxis(:STG_BETA0; min = -6.0, max = 6.0, N = 5, rule = LinRule())
y = GridAxis(:STG_ALPHA0; min = -1e-4, max = -1e-1, N = 5, rule = LogRule())

ref_settings = RefinementSettings(
    DiffContourUnit(:chi2_marginalized; from_min = true, diffs = [1.0], contours = [lvl_2sigma]),
    LocalMinimaUnit(:chi2_marginalized; from_min = true, max = 20.0, max_diff = 0.1),
    params_to_save = (:chi2_marginalized, :wrms_fit, :wrms_tn_fit),
    desired_refinement_level = 0,
)

task = Adaptive2DGridTask(
    base_task = prior_task,
    x = x,
    y = y,
    ref_settings = ref_settings,
    opts = GridWorkspaceOptions(grid_root = "example_grid_scan"),
)

println(task)

# Replace `Tempo2(fixture_root)` with your real TEMPO2 runtime directory, then uncomment:
# result = run_task(task)
# @show result.min[:chi2_marginalized]
