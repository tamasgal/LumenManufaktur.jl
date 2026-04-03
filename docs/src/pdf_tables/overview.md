# PDF Tables

Evaluating the analytical light-yield functions — especially the scattered-light
integrals — can be expensive: a single call to
[`scatteredlightfromEMshowers`](@ref) takes roughly **480 µs** on a modern CPU.
When fitting tracks or showers against many PMT hits, the same PDF must be
evaluated thousands of times at slightly different geometry parameters.

PDF tables replace the analytical evaluation with:

1. **Build once** — evaluate the analytical PDF on a parameter grid, convert
   each node's time axis to a CDF-quantile representation, and store the result.
2. **Look up fast** — at runtime, multilinear-interpolate the grid and do a
   short linear scan over 100 quantile knots.

Because the quantile knots are placed automatically where the PDF has high
density (near the Cherenkov peak), the representation is accurate without any
hand-tuned coordinate transforms.

!!! note "Build time vs. evaluation time"
    Building a table is a one-time cost, typically a few seconds to a few
    minutes depending on grid size and the number of physics integration
    points used.  The resulting table can be saved to an HDF5 file and
    loaded in subsequent runs.

## Available table types

| Type | Grid dimensions | Channels | Make function |
|---|---|---|---|
| [`DirectMuonPDFTable`](@ref) | `(R, θ, ϕ)` | 1 — direct muon | [`makepdf`](@ref) |
| [`MuonPDFTable`](@ref) | `(R, θ, ϕ)` | 6 — all muon-light contributions | [`makemuonpdf`](@ref) |
| [`EMShowerPDFTable`](@ref) | `(D, cd, θ, ϕ)` | 2 — direct + scattered EM shower (per GeV) | [`makeemshowerpdf`](@ref) |
| [`BrightPointPDFTable`](@ref) | `(D, ct)` | 2 — direct + scattered bright point | [`makebrightpointpdf`](@ref) |

## Performance

Benchmarked at `R = D = 10 m`, `θ = ϕ = π/2`, `Δt ∈ [−1, 15] ns`,
`n_quantiles = 100` (Apple M-series, single thread, Julia 1.12):

| Analytical call(s) | Analytical | Table | Speedup |
|---|---|---|---|
| `directlightfrommuon` | 667 ns | 451 ns (`DirectMuonPDFTable`) | 1.5× |
| `directlightfrommuon` | 667 ns | — | — |
| `scatteredlightfrommuon` | 22 200 ns | — | — |
| `directlightfromdeltarays` | 1 800 ns | — | — |
| `scatteredlightfromdeltarays` | 101 000 ns | — | — |
| `directlightfromEMshowers` | 3 800 ns | — | — |
| `scatteredlightfromEMshowers` | 477 000 ns | — | — |
| **all 6 muon channels** | **606 000 ns** | **1 261 ns** (`MuonPDFTable`) | **~495×** |
| `directlightfromEMshower` (5-param) | 400 ns | — | — |
| `scatteredlightfromEMshower` (5-param) | 271 000 ns | — | — |
| **direct + scattered EM shower** | **271 000 ns** | **1 370 ns** (`EMShowerPDFTable`) | **~190×** |
| `directlightfrombrightpoint` | 307 ns | — | — |
| `scatteredlightfrombrightpoint` | 9 200 ns | — | — |
| **direct + scattered bright point** | **9 500 ns** | **378 ns** (`BrightPointPDFTable`) | **~25×** |

The large speedups come from the expensive scattered-light integrals (Legendre
quadrature over emission positions and azimuth angles), which dominate runtime
in the analytical path but are fully amortised into the table build.

## Portability

All tables are saved and loaded as HDF5 files using
[`HDF5.jl`](https://github.com/JuliaIO/HDF5.jl).  The on-disk format is
readable from Python (`h5py`), C++, MATLAB, and any other HDF5-aware tool,
making it straightforward to share pre-built tables or consume them in a
reconstruction framework written in another language.

## Configuration

Each table type has a corresponding `*Config` struct with keyword arguments
for the grid spacing, number of quantiles, time window, and integration
resolution.  All fields have sensible defaults; override only what you need:

```julia
cfg = MuonPDFTableConfig(
    n_quantiles   = 150,    # more knots → higher accuracy
    t_max_scattered = 300.0, # wider window for large-R scattered light [ns]
    n_integration = 8000,   # finer trapezoidal CDF integration
)
tbl = makemuonpdf(params, pmt, cfg)
```

The scattered-light channels use a separate, wider time window
(`t_max_scattered`, default 200 ns) because their PDF tails extend much
further than the direct-light channels (`t_max`, default 25 ns).

See the sub-pages for detailed usage and examples for each table type.
