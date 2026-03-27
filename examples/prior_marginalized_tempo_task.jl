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
    job_name = "example_prior",
    overwrite = :reuse,
    layout = :split,
)

base_task = BasicTempoTask(base_settings)

prior_settings = PriorMarginalizationSettings(
    parameter = :DDOT,
    pin_mode = :fixed,
    prior = AnalyticPrior(Normal(0.0, 1.0)),
    nodes = ClenshawCurtisNodes(4),
    likelihood_source = :chi2_fit,
    representative = :prior_median,
    save_node_results = true,
    exec_options = PriorExecutionOptions(
        mode = :chained,
        chain_direction = :backward,
        chain_snapshot_par = true,
        scheduler = :serial,
    ),
)

task = PriorMarginalizedTempoTask(base_task, prior_settings)

println(task)

# Replace `Tempo2(fixture_root)` with your real TEMPO2 runtime directory, then uncomment:
# result = run_task(task)
# @show result.success result.status
