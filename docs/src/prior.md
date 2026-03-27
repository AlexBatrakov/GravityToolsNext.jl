# Prior Marginalization

Use `PriorMarginalizedTempoTask` to marginalize over a parameter using a prior and a node rule.
It chains or independently executes a series of runs with per-node overrides, aggregating results.

Current support boundary:
- serial node execution is implemented
- `scheduler = :distributed` is reserved for future work
- `keep_node_dirs` is currently an advisory setting rather than an enforced cleanup policy

## Settings
- `PriorMarginalizationSettings` controls the parameter (`parameter::Symbol`), the prior (`AnalyticPrior | GridPrior | SampledPrior`),
  node selection (e.g., `ClenshawCurtisNodes`), and execution options (`PriorExecutionOptions`).
- Likelihood source and representative can be customized (e.g., `:chi2_fit`, `:prior_median`).

## Example
```julia
using Distributions: Normal

base = BasicTempoTask(TempoRunSettings(
  work_dir="/abs/work",
  par_input="a.par",
  tim_input="a.tim",
  par_output="a_out.par",
  tempo_version=Tempo2("/path/to/TEMPO2"),
))

ps = PriorMarginalizationSettings(
  parameter=:DDOT,
  pin_mode=:fixed,
  prior=AnalyticPrior(Normal(0.0, 1.0)),
  nodes=ClenshawCurtisNodes(6),
  likelihood_source=:chi2_fit,
  representative=:prior_median,
  save_node_results=true,
  exec_options = PriorExecutionOptions(mode=:chained, chain_direction=:backward, chain_snapshot_par=true)
)

task = PriorMarginalizedTempoTask(base, ps)
# Uncomment when TEMPO2 and input files are available:
# res = run_task(task)
```
