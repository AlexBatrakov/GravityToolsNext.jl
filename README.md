# GravityToolsNext.jl

[![Docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alexbatrakov.github.io/GravityToolsNext.jl/dev/)
[![Docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://alexbatrakov.github.io/GravityToolsNext.jl/stable/)
[![CI](https://github.com/AlexBatrakov/GravityToolsNext.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/AlexBatrakov/GravityToolsNext.jl/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Modern Julia tooling for orchestrating TEMPO/TEMPO2 runs with reproducible job layouts, robust path handling, and convenient analysis hooks.

At a high level, this package aims to turn “a TEMPO2 run on my laptop” into a repeatable, inspectable, scriptable workflow:
- deterministic job directory layouts (easy to diff, archive, and share)
- explicit run settings (inputs, overrides, output policies)
- structured parsing of outputs into Julia-friendly objects
- higher-level task abstractions (single runs, priors/nodes, adaptive grids)

Typical use-cases:
- Reproducible TEMPO/TEMPO2 batch runs with clean job directories
- Parameter sweeps / adaptive grids around expensive likelihood surfaces
- Parsing and post-processing of outputs (including white-noise diagnostics)

## Features
- Clean settings API (RunPaths, EngineOptions, WorkspaceOptions, LoggingOptions)
- Materialized job execution with flat/split layouts
- Optional manifests and cleanup policies
- White-noise analysis toggles and scope control
- Extensible task layer for parameter sweeps and adaptive grids

## What’s inside

The codebase is split into two major parts:

- **TempoFramework** (external dependency at execution time)
	- `TempoRunSettings` and related option structs define *what* to run and *how* to lay out files.
	- `run_tempo_raw` / `run_tempo_parsed` execute the external binary and parse outputs.
	- Task layer (`BasicTempoTask`, `PriorMarginalizedTempoTask`, `Adaptive2DGridTask`) composes workflows.

- **AdaptiveGridFramework** (pure Julia)
	- `GridAxis` + grid rules (`LinRule`, `LogRule`, `ExplicitRule`) define parameter grids.
	- `RefinementSettings` + refinement units decide where refinement happens.
	- `AdaptiveRefinement2DGrid` performs solver-agnostic grid evaluation/refinement.

This split is intentional: **unit tests and CI cover the pure-Julia pieces**, while end-to-end task runs require a local TEMPO/TEMPO2 installation.

## Mental model

- If you want to **run TEMPO/TEMPO2 reproducibly**, start from `TempoRunSettings` → create a task (e.g. `BasicTempoTask`) → `run_task`.
- If you want to **build and refine grids without any external binaries**, start from `GridAxis` + `RefinementSettings` → `AdaptiveRefinement2DGrid`.
- If you want to **do parameter exploration**, use the task layer (priors/nodes, adaptive grids) to compose many runs.

## Installation

Supported Julia versions: **1.10+**.

```julia
using Pkg
Pkg.add(url = "https://github.com/AlexBatrakov/GravityToolsNext.jl")
```

## Quickstart

### Pure Julia (no TEMPO/TEMPO2 required)

The `AdaptiveGridFramework` is solver-agnostic and can be used without external binaries.

```julia
using GravityToolsNext

x = GridAxis(:x; min=-1.0, max=1.0, N=5, rule=LinRule())
y = GridAxis(:y; min=-1.0, max=1.0, N=5, rule=LinRule())

ref_sets = RefinementSettings(
	FullUnit(:z);
	desired_refinement_level=0,
	parallel=false,
	params_to_save=(:z,)
)

grid = AdaptiveRefinement2DGrid(x, y, ref_sets)

target(xv, yv) = (z = xv^2 + yv^2,)
precalculate_2DGrid!(grid, target, _ -> nothing)

@show grid.min[:z] grid.max[:z]
```

### TEMPO/TEMPO2 (external dependency)

End-to-end task execution requires installed TEMPO/TEMPO2 binaries and data directories.
For `Tempo2()`, make sure `TEMPO2` is set (data dir). Optionally set `TEMPO2_CMD` to the executable.

If you only plan to use the pure-Julia parts (like `AdaptiveGridFramework`), you don’t need any TEMPO/TEMPO2 setup.

```julia
using GravityToolsNext

s = TempoRunSettings(
	work_dir   = "/abs/workdir",
	par_input  = "example.par",
	tim_input  = "example.tim",
	par_output = "example_out.par",
	# Either rely on ENV["TEMPO2"], or pass a directory explicitly:
	tempo_version = Tempo2("/path/to/tempo2_data"),
)

# validate inputs and run
validate(s)
# result = run_tempo_parsed(s)
```

## Examples

End-to-end usage scripts (sandbox / manual runs) live in [examples/](examples/). Many of them require installed TEMPO/TEMPO2 and appropriate environment variables (e.g. `TEMPO2`, `TEMPO2_CMD`) pointing to the data directory / executable.

## Testing

Unit tests focus on pure-Julia components (for example `AdaptiveGridFramework`) and do **not** require TEMPO/TEMPO2.

CI runs this unit-test subset on clean machines.

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Documentation
- Dev: https://alexbatrakov.github.io/GravityToolsNext.jl/dev/
- Stable: https://alexbatrakov.github.io/GravityToolsNext.jl/stable/

Key pages:
- Quickstart • Settings • Runner • Tasks • API Reference

## Contributing
Issues and PRs are welcome. Please run tests and keep public APIs documented.

## License
MIT. See [LICENSE](LICENSE).

## Citation
If this package contributes to your work, please cite the repository. A CITATION.cff can be added upon request.