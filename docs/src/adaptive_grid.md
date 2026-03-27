# Adaptive Grid — Overview

Use the Adaptive Grid Framework to evaluate a Tempo task over a 2D parameter space efficiently.
It refines the grid where signal changes rapidly and saves per-point results when requested.

- Grid axes: defined with `GridAxis(...; rule=...)`
- Refinement strategy: configured with `RefinementSettings` and units like `LocalMinimaUnit`, `FullUnit`, etc.
- Execution: orchestrated by `Adaptive2DGridTask`, which runs a supported `SingleTempoTask` at each grid point

## When to use
- Parameter scans and likelihood/chi² maps
- Contour extraction and coarse-to-fine searches

## Minimal pipeline

```julia
using GravityToolsNext
using Distributions: Normal

# 1) Define your base Tempo task
s = TempoRunSettings(
    work_dir="/abs/work",
    par_input="a.par",
    tim_input="a.tim",
    par_output="a_out.par",
    tempo_version=Tempo2("/path/to/TEMPO2"),
)
base = BasicTempoTask(s)

# 2) Wrap it in a prior-marginalized task if you want grid metrics such as
#    :chi2_marginalized
prior = AnalyticPrior(Normal(0.0, 1.0))
ps = PriorMarginalizationSettings(
    parameter = :DDOT,
    pin_mode = :fixed,
    prior = prior,
    nodes = ClenshawCurtisNodes(4),
    likelihood_source = :chi2_fit,
    representative = :prior_median,
)
prior_task = PriorMarginalizedTempoTask(base, ps)

# 3) Define axes
x = GridAxis(:PX; min=1.0, max=10.0, N=21, rule=LinRule())
y = GridAxis(:PY; min=1e-3, max=1.0, N=21, rule=LogRule(+1))

# 4) Define refinement
ref = RefinementSettings(
    LocalMinimaUnit(:chi2_marginalized),
    desired_refinement_level = 0,
    params_to_save = (:chi2_marginalized, :wrms_fit),
)

# 5) Build the grid task
opts = GridWorkspaceOptions(grid_root = "scan")
gtask = Adaptive2DGridTask(base_task=prior_task, x=x, y=y, ref_settings=ref, opts=opts)

# Uncomment when TEMPO2 and input files are available:
# result = run_task(gtask)
```

`Adaptive2DGridTask` relies on wrapper hooks implemented by the base task:
`task_copy_with`, `task_derive_par_output`, and `run_task`. Per-point JLD2
persistence additionally uses `save_result_jld2`.
