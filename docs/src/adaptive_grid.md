# Adaptive Grid — Overview

Use the Adaptive Grid Framework to evaluate a Tempo task over a 2D parameter space efficiently.
It refines the grid where signal changes rapidly and saves per-point results when requested.

- Grid axes: defined via `LinAxis`, `LogAxis`, or `ExplicitAxis` from `GridAxis` rules
- Refinement strategy: configured with `RefinementSettings` and units like `LocalMinimaUnit`, `FullUnit`, etc.
- Execution: orchestrated by `Adaptive2DGridTask`, which runs your `SingleTempoTask` at each grid point

## When to use
- Parameter scans and likelihood/chi² maps
- Contour extraction and coarse-to-fine searches

## Minimal pipeline

```julia
using GravityToolsNext

# 1) Define your base Tempo task
s = TempoRunSettings(
    work_dir="/abs/work", par_input="a.par", tim_input="a.tim", par_output="a_out.par", tempo_version=Tempo2())
base = BasicTempoTask(s)

# 2) Define axes
x = LinAxis(:PX, 1.0, 10.0, 21)
y = LogAxis(:PY, 1e-3, 1.0, 21)

# 3) Define refinement
ref = RefinementSettings(LocalMinimaUnit(:chi2_marginalized))

# 4) Build and run the grid task
opts = GridWorkspaceOptions(grid_root = "scan")
gtask = Adaptive2DGridTask(base_task=base, x=x, y=y, ref_settings=ref, opts=opts)

result = run_task(gtask)
```
