# Muon PDF Tables

Two table types cover muon-related light: [`DirectMuonPDFTable`](@ref) for the
direct channel only, and [`MuonPDFTable`](@ref) for all six contributions that
make up [`lightfrommuon`](@ref).

## DirectMuonPDFTable

A compact, single-channel table for [`directlightfrommuon`](@ref) on a
`(R, θ, ϕ)` grid.  Useful when only direct light is needed or when memory is
at a premium.

```@example muon_table
using CairoMakie, LumenManufaktur

params = LMParameters(dispersion_model = DispersionORCA)
pmt    = LumenManufaktur.KM3NeTPMT

# Compact grid for the docs example; production code uses the full defaults
cfg_direct = PDFTableConfig(
    R_values = [1.0, 5.0, 10.0, 20.0, 50.0, 100.0],
    n_quantiles = 100, n_integration = 200,
)
tbl_direct = makepdf(params, pmt, cfg_direct)
println(tbl_direct)
```

Calling the table returns the time PDF [npe/ns]:

```@example muon_table
R, θ, ϕ = 10.0, π/2, π/2

Δts = range(-2, 10, 400)
analytical = [directlightfrommuon(params, pmt, R, θ, ϕ, Δt) for Δt in Δts]
table_vals = [tbl_direct(R, θ, ϕ, Δt)                       for Δt in Δts]

fig = Figure(size = (700, 380))
ax  = Axis(fig[1,1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns",
    title  = "DirectMuonPDFTable — R = $(Int(R)) m",
    xgridstyle = :dash, ygridstyle = :dash,
)
lines!(ax, collect(Δts), analytical, label = "analytical", linewidth = 2)
lines!(ax, collect(Δts), table_vals, label = "table",      linewidth = 2, linestyle = :dash)
axislegend(ax; position = :rt)
fig
```

Call without `Δt` to get the integrated npe:

```@example muon_table
npe_analytical = directlightfrommuon(params, pmt, R, θ, ϕ)
npe_table      = tbl_direct(R, θ, ϕ)
println("npe — analytical: $(round(npe_analytical, digits=5)),  table: $(round(npe_table, digits=5))")
```

### Save and load

```julia
save_pdftable("direct_muon.h5", tbl_direct)
tbl_loaded = load_pdftable("direct_muon.h5")
```

---

## MuonPDFTable

Stores all six muon-light channels on the same `(R, θ, ϕ)` grid:

| Channel key | Source function | Units |
|---|---|---|
| `direct_muon` | [`directlightfrommuon`](@ref) | npe/ns |
| `scattered_muon` | [`scatteredlightfrommuon`](@ref) | npe/ns |
| `direct_dr` | [`directlightfromdeltarays`](@ref) | npe/ns per GeV/m |
| `scattered_dr` | [`scatteredlightfromdeltarays`](@ref) | npe/ns per GeV/m |
| `direct_em` | [`directlightfromEMshowers`](@ref) | npe/ns per GeV |
| `scattered_em` | [`scatteredlightfromEMshowers`](@ref) | npe/ns per GeV |

```@example muon_table
# Minimal grid for the docs example; production code uses the full defaults.
# Building with scattered-light channels is slow — see the overview page for
# production-recommended grid sizes and the typical build-time tradeoffs.
cfg_muon = MuonPDFTableConfig(
    R_values = [5.0, 10.0, 20.0],
    θ_values = [π/4, π/2, 3π/4],
    ϕ_values = [π/4, π/2, 3π/4],
    n_quantiles = 100, n_integration = 100,
    t_max_scattered = 100.0,
)
tbl_muon = makemuonpdf(params, pmt, cfg_muon)
println(tbl_muon)
```

Calling with `(R, θ, ϕ, Δt)` returns a `NamedTuple` of all six channels:

```@example muon_table
channels = tbl_muon(R, θ, ϕ, 2.0)
println(channels)
```

Calling with `(R, θ, ϕ, E, Δt)` applies the same energy weighting as
[`lightfrommuon`](@ref) and returns a scalar total [npe/ns]:

```@example muon_table
E   = 100.0   # muon energy [GeV]
Δts = range(-2, 80, 600)

analytical_total = [lightfrommuon(params, pmt, E, R, θ, ϕ, Δt) for Δt in Δts]
table_total      = [tbl_muon(R, θ, ϕ, E, Δt)                   for Δt in Δts]

fig2 = Figure(size = (700, 380))
ax2  = Axis(fig2[1,1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns",
    title  = "MuonPDFTable total — R = $(Int(R)) m, E = $(Int(E)) GeV",
    xgridstyle = :dash, ygridstyle = :dash,
)
lines!(ax2, collect(Δts), analytical_total, label = "analytical", linewidth = 2)
lines!(ax2, collect(Δts), table_total,      label = "table",      linewidth = 2, linestyle = :dash)
axislegend(ax2; position = :rt)
fig2
```

Calling without `Δt` returns a `NamedTuple` of integrated npe per channel:

```@example muon_table
npe = tbl_muon(R, θ, ϕ)
println("direct_muon npe:  $(round(npe.direct_muon, digits=5))")
println("scattered_muon npe: $(round(npe.scattered_muon, digits=5))")
```

### Save and load

```julia
save_muonpdftable("muon.h5", tbl_muon)
tbl_loaded = load_muonpdftable("muon.h5")
```
