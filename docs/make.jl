# docs/make.jl
using Documenter
using GravityToolsNext

# чтобы doctest-ы из докстрингов знали, что импортировать
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
    pages = [
        "Home" => "index.md",
        "API"  => "api.md",
    ],
    doctest = true,
)

# автодеплой (позже добавим GitHub Actions)
deploydocs(
    repo = "github.com/AlexBatrakov/GravityToolsNext.jl.git",  # замени на свой, например: github.com/AlexBatrakov/GravityToolsNext.jl
    devbranch = "main",                 # или "master"
)