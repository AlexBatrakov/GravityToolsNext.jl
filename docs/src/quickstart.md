# Quickstart

This guide shows the currently supported minimal run path for a TEMPO2-backed
job in GravityToolsNext.jl.

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
    tempo_version = Tempo2("/path/to/TEMPO2"), # active parser path is Tempo2-oriented
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
- See Tasks for higher-level orchestration and wrapper hooks.
