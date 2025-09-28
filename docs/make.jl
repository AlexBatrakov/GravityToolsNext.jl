# docs/make.jl
using Documenter
using GravityToolsNext

# Make doctests know what to import
DocMeta.setdocmeta!(GravityToolsNext, :DocTestSetup, :(using GravityToolsNext); recursive=true)

makedocs(
    sitename = "GravityToolsNext.jl",
    modules  = [GravityToolsNext],
    format   = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        assets = String[],
        size_threshold      = 1024*1024,  # 1 MiB вместо 200 KiB
        size_threshold_warn = 512*1024,   # 512 KiB вместо 100 KiB
    ),
    checkdocs = :none,
    pages = [
        "Home"       => "index.md",
        "Quickstart" => "quickstart.md",
        "Tempo Framework" => Any[
            "tempo_overview.md",
            "tempo_tasks_basic.md",
            "results.md",
            "prior.md",
            "tempo_tasks_grid_prior.md",
        ],
        "Guides" => Any[
            "settings.md",
            "runner.md",
            "tasks.md",
        ],
        "Adaptive Grid" => Any[
            "adaptive_grid.md",
            "grid_axes.md",
            "refinement.md",
            "grid_tasks.md",
        ],
        "API"        => "api.md",
    ],
    doctest = true,
)

# Autodeploy (CI can be set up later via GitHub Actions)
deploydocs(
    repo = "github.com/AlexBatrakov/GravityToolsNext.jl.git",
    devbranch = "main",
)