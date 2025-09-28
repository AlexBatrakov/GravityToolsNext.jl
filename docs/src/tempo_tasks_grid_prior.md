# Adaptive Grid with Prior Marginalization

This advanced workflow runs `PriorMarginalizedTempoTask` over a 2D parameter space using `Adaptive2DGridTask`.

## Example
```julia
base = BasicTempoTask(TempoRunSettings(
  work_dir="/abs/work", par_input="base.par", tim_input="base.tim", par_output="base_out.par", tempo_version=Tempo2()))

ps = PriorMarginalizationSettings(
  parameter=:DDOT, pin_mode=:fixed, prior=SampledPrior("DDOT_prior.dat"), nodes=ClenshawCurtisNodes(4),
  likelihood_source=:chi2_fit, representative=:prior_median,
  save_node_results=true, exec_options=PriorExecutionOptions(mode=:chained, chain_direction=:backward, chain_snapshot_par=true)
)

prior_task = PriorMarginalizedTempoTask(base, ps)

x = GridAxis(:STG_BETA0, min=-6.0, max=6.0, N=4, rule=LinRule())
y = GridAxis(:STG_ALPHA0, min=-1e-4, max=-1e-1, N=4, rule=LogRule())

ref = RefinementSettings(
  params_to_save = (:chi2_marginalized, :wrms_fit, :ad_white_fit),
  desired_refinement_level = 0,
  DiffContourUnit(:chi2_marginalized, from_min=true, diffs=[1.0], contours=[lvl_2sigma]),
  LocalMinimaUnit(:chi2_marginalized, from_min=true, max=20.0, max_diff=0.1)
)

opts = GridWorkspaceOptions(grid_root="GRID_TEST")

gtask = Adaptive2DGridTask(base_task=prior_task, x=x, y=y, ref_settings=ref, opts=opts)

result = run_task(gtask)
```
