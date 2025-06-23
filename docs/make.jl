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
            "examples/directlightfrommuon.md",
            "examples/npes_and_zenith_angles.md",
        ],
        "API" => "api.md",
    ],
    repo = Documenter.Remotes.URL(
        "https://git.km3net.de/simulation/LumenManufaktur.jl/blob/{commit}{path}#L{line}",
        "https://git.km3net.de/simulation/LumenManufaktur.jl",
    ),
)

deploydocs(;
    repo = "git.km3net.de/simulation/LumenManufaktur.jl",
    devbranch = "main",
    push_preview = true,
)
