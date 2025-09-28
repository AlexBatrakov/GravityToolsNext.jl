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

## Parameter sweeps and adaptive grids
- MultiPointTasks enable structured exploration (WIP docs)
