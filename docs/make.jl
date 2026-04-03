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
        "Light Properties" => Any[
            "light_properties/absorption.md",
            "light_properties/scattering.md",
            "light_properties/dispersion.md",
        ],
        "Light Yield" => Any[
            "light_yield/direct_light.md",
            "light_yield/scattered_light.md",
            "light_yield/combined_light.md",
        ],
        "PDF Tables" => Any[
            "pdf_tables/overview.md",
            "pdf_tables/muon_table.md",
            "pdf_tables/shower_brightpoint_tables.md",
        ],
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
