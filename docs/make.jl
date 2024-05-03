using Documenter, LumenManufaktur

makedocs(;
    modules = [LumenManufaktur],
    sitename = "LumenManufaktur.jl",
    authors = "Tamas Gal",
    format = Documenter.HTML(;
        assets = ["assets/custom.css"],
        sidebar_sitename = true,
        collapselevel = 2,
        warn_outdated = true,
    ),
    warnonly = [:missing_docs],
    pages = [
        "Home" => "index.md",
        "Examples" => Any[
            "examples/an_example.md",
        ],
        "API" => "api.md"
    ],
    repo = Documenter.Remotes.URL(
        "https://git.km3net.de/tgal/LumenManufaktur.jl/blob/{commit}{path}#L{line}",
        "https://git.km3net.de/tgal/LumenManufaktur.jl"
    ),
)

deploydocs(;
  repo = "git.km3net.de/tgal/LumenManufaktur.jl",
  devbranch = "main",
  push_preview=true
)
