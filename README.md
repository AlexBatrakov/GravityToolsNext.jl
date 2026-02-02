# GravityToolsNext.jl

[![Docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://alexbatrakov.github.io/GravityToolsNext.jl/dev/)
[![Docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://alexbatrakov.github.io/GravityToolsNext.jl/stable/)
[![CI](https://github.com/AlexBatrakov/GravityToolsNext.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/AlexBatrakov/GravityToolsNext.jl/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Modern Julia tooling for orchestrating TEMPO/TEMPO2 runs with reproducible job layouts, robust path handling, and convenient analysis hooks.

## Features
- Clean settings API (RunPaths, EngineOptions, WorkspaceOptions, LoggingOptions)
- Materialized job execution with flat/split layouts
- Optional manifests and cleanup policies
- White-noise analysis toggles and scope control
- Extensible task layer for parameter sweeps and adaptive grids

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/AlexBatrakov/GravityToolsNext.jl")
```

## Quickstart

```julia
using GravityToolsNext

s = TempoRunSettings(
	work_dir   = "/abs/workdir",
	par_input  = "example.par",
	tim_input  = "example.tim",
	par_output = "example_out.par",
	tempo_version = Tempo2(),
)

# validate inputs and run
validate(s)
# result = run_tempo_parsed(s)
```

## Examples

End-to-end usage scripts (sandbox / manual runs) live in [examples/](examples/). Many of them require installed TEMPO/TEMPO2 and appropriate environment variables (e.g. `TEMPO2`, `TEMPO2_CMD`) pointing to the data directory / executable.

## Testing

Unit tests focus on pure-Julia components (for example `AdaptiveGridFramework`) and do **not** require TEMPO/TEMPO2.

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