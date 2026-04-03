# EM Shower and Bright Point PDF Tables

## EMShowerPDFTable

[`EMShowerPDFTable`](@ref) covers the 5-parameter (per-GeV, direction-resolved)
forms of [`directlightfromEMshower`](@ref) and
[`scatteredlightfromEMshower`](@ref) on a `(D, cd, θ, ϕ)` grid, where `cd` is
the cosine of the angle between the shower direction and the shower–PMT vector.

```@example shower_bp_table
using CairoMakie, LumenManufaktur

params = LMParameters(dispersion_model = DispersionORCA)
pmt    = LumenManufaktur.KM3NeTPMT

# Minimal grid for the docs example; production code uses the full defaults.
cfg_em = EMShowerPDFTableConfig(
    D_values  = [5.0, 10.0, 20.0],
    cd_values = [-0.5, 0.0, 0.5],
    θ_values  = [π/4, π/2, 3π/4],
    ϕ_values  = [π/4, π/2, 3π/4],
    n_quantiles = 100, n_integration = 100,
    t_max_scattered = 100.0,
)
tbl_em = makeemshowerpdf(params, pmt, cfg_em)
println(tbl_em)
```

Calling with `(D, cd, θ, ϕ, Δt)` returns `(direct, scattered)` as a
`NamedTuple` [npe/ns/GeV]:

```@example shower_bp_table
D, cd, θ, ϕ = 10.0, 0.0, π/2, π/2

Δts  = range(-2, 50, 500)
d_an = [directlightfromEMshower(params,  pmt, D, cd, θ, ϕ, Δt) for Δt in Δts]
s_an = [scatteredlightfromEMshower(params, pmt, D, cd, θ, ϕ, Δt) for Δt in Δts]
d_tb = [tbl_em(D, cd, θ, ϕ, Δt).direct    for Δt in Δts]
s_tb = [tbl_em(D, cd, θ, ϕ, Δt).scattered for Δt in Δts]

fig = Figure(size = (700, 380))
ax  = Axis(fig[1,1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns / GeV",
    title  = "EMShowerPDFTable — D = $(Int(D)) m",
    xgridstyle = :dash, ygridstyle = :dash,
)
lines!(ax, collect(Δts), d_an, label = "direct (analytical)",    linewidth = 2)
lines!(ax, collect(Δts), s_an, label = "scattered (analytical)", linewidth = 2)
lines!(ax, collect(Δts), d_tb, label = "direct (table)",    linewidth = 2, linestyle = :dash)
lines!(ax, collect(Δts), s_tb, label = "scattered (table)", linewidth = 2, linestyle = :dash)
axislegend(ax; position = :rt)
fig
```

Calling with `(D, cd, θ, ϕ, E, Δt)` scales by shower energy `E` [GeV] and
returns a total [npe/ns]:

```@example shower_bp_table
E = 50.0  # GeV

Δts2     = range(-2, 50, 500)
total_an = [lightfromEMshower(params, pmt, D, cd, θ, ϕ, Δt) * E for Δt in Δts2]
total_tb = [tbl_em(D, cd, θ, ϕ, E, Δt)                          for Δt in Δts2]

fig2 = Figure(size = (700, 380))
ax2  = Axis(fig2[1,1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns",
    title  = "EMShowerPDFTable total — D = $(Int(D)) m, E = $(Int(E)) GeV",
    xgridstyle = :dash, ygridstyle = :dash,
)
lines!(ax2, collect(Δts2), total_an, label = "analytical", linewidth = 2)
lines!(ax2, collect(Δts2), total_tb, label = "table",      linewidth = 2, linestyle = :dash)
axislegend(ax2; position = :rt)
fig2
```

Calling without `Δt` returns per-GeV npe:

```@example shower_bp_table
npe = tbl_em(D, cd, θ, ϕ)
println("npe per GeV — direct: $(round(npe.direct, digits=5)), scattered: $(round(npe.scattered, digits=5))")
```

### Save and load

```julia
save_emshowerpdftable("emshower.h5", tbl_em)
tbl_loaded = load_emshowerpdftable("emshower.h5")
```

---

## BrightPointPDFTable

[`BrightPointPDFTable`](@ref) covers [`directlightfrombrightpoint`](@ref) and
[`scatteredlightfrombrightpoint`](@ref) on a `(D, ct)` grid, where `ct` is the
cosine of the PMT orientation angle.

```@example shower_bp_table
# Compact grid for the docs example; production code uses the full defaults
cfg_bp = BrightPointPDFTableConfig(
    D_values  = [1.0, 5.0, 10.0, 20.0, 50.0],
    ct_values = collect(range(-1.0, 1.0; length = 7)),
    n_quantiles = 100, n_integration = 200,
    t_max_scattered = 150.0,
)
tbl_bp = makebrightpointpdf(params, pmt, cfg_bp)
println(tbl_bp)
```

Calling with `(D, ct, Δt)` returns `(direct, scattered)` [npe/ns]:

```@example shower_bp_table
D_bp, ct_bp = 10.0, 0.0

Δts3  = range(-2, 50, 500)
d_an3 = [directlightfrombrightpoint(params,  pmt, D_bp, ct_bp, Δt) for Δt in Δts3]
s_an3 = [scatteredlightfrombrightpoint(params, pmt, D_bp, ct_bp, Δt) for Δt in Δts3]
d_tb3 = [tbl_bp(D_bp, ct_bp, Δt).direct    for Δt in Δts3]
s_tb3 = [tbl_bp(D_bp, ct_bp, Δt).scattered for Δt in Δts3]

fig3 = Figure(size = (700, 380))
ax3  = Axis(fig3[1,1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns",
    title  = "BrightPointPDFTable — D = $(Int(D_bp)) m",
    xgridstyle = :dash, ygridstyle = :dash,
)
lines!(ax3, collect(Δts3), d_an3, label = "direct (analytical)",    linewidth = 2)
lines!(ax3, collect(Δts3), s_an3, label = "scattered (analytical)", linewidth = 2)
lines!(ax3, collect(Δts3), d_tb3, label = "direct (table)",    linewidth = 2, linestyle = :dash)
lines!(ax3, collect(Δts3), s_tb3, label = "scattered (table)", linewidth = 2, linestyle = :dash)
axislegend(ax3; position = :rt)
fig3
```

Calling without `Δt` returns integrated npe:

```@example shower_bp_table
npe_bp = tbl_bp(D_bp, ct_bp)
println("npe — direct: $(round(npe_bp.direct, digits=5)), scattered: $(round(npe_bp.scattered, digits=5))")
```

### Save and load

```julia
save_brightpointpdftable("brightpoint.h5", tbl_bp)
tbl_loaded = load_brightpointpdftable("brightpoint.h5")
```
