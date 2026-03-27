# Tasks

Higher-level orchestration built on top of the runner.

## Basic task
- Prepares settings
- Executes run and parses outputs
- Optionally stages inputs into a node directory for batch sweeps

Example skeleton:

```julia
# using GravityToolsNext
# t = BasicTempoTask(settings)
# res = run_task(t)
```

## Wrapper hooks
Higher-level wrappers currently rely on a small explicit extension contract for
single tasks:
- `task_copy_with`
- `task_derive_par_output`
- `task_stage_inputs!`

Wrappers that persist per-point results additionally use `save_result_jld2`.

## Parameter sweeps and adaptive grids
- `PriorMarginalizedTempoTask` is supported with the serial scheduler.
- `Adaptive2DGridTask` is supported for tasks implementing the wrapper hooks above.
- Distributed scheduling remains deferred.
