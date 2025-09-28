# Quickstart

This guide shows how to run a simple TEMPO/TEMPO2 job with GravityToolsNext.jl.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/AlexBatrakov/GravityToolsNext.jl")
```

## Minimal run

```julia
using GravityToolsNext

s = TempoRunSettings(
    work_dir   = "/path/to/workdir",           # absolute path
    par_input  = "example.par",                # relative to work_dir, must end with .par
    tim_input  = "example.tim",                # relative to work_dir, must end with .tim
    par_output = "example_out.par",            # filename only
    tempo_version = Tempo2(),                   # or Tempo()
)

# Validate inputs (optional)
validate(s)

# Execute and parse
# result = run_tempo_parsed(s)
# show(result)
```

## Next steps
- See Settings for all options.
- See Runner for materialization and job layouts.
- See Tasks for higher-level orchestration.
