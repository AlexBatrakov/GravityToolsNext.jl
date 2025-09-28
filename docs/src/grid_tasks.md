# Adaptive 2D Grid Task

The `Adaptive2DGridTask` adapts a `SingleTempoTask` over a 2D parameter space.
It manages per-point workspaces, parameter overrides, and optional result persistence.

## Workspace options
`GridWorkspaceOptions` controls where grid artifacts live:
- `grid_root::String` — base directory for the grid (absolute or relative to base task work_dir)
- `point_job_prefix::String = "grid_points"` — subdir under grid_root for per-point job dirs
- `results_dirname::String   = "results"` — subdir for saved per-point results (.jld2)
- `input_dirname::String     = "input"` — optional staging location for base inputs

Tagging and staging:
- `tag_mode::Symbol = :with_value | :hash` (controls point tag format)
- `save_results_jld2::Bool = true`
- `stage_inputs::Bool = true`; `stage_inputs_mode = :root | :subdir`

## Example
```julia
base = BasicTempoTask(TempoRunSettings(
    work_dir="/abs/work", par_input="a.par", tim_input="a.tim", par_output="a_out.par", tempo_version=Tempo2()))

x = LinAxis(:PX, 0.0, 5.0, 51)
y = LogAxis(:PY, 1e-4, 1e0, 41)
ref = RefinementSettings(LocalMinimaUnit(:chi2_marginalized))

opts = GridWorkspaceOptions(grid_root="scan")

gtask = Adaptive2DGridTask(base_task=base, x=x, y=y, ref_settings=ref, opts=opts)
res = run_task(gtask)
```
