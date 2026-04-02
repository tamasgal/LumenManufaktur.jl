# Scattered Light

Scattered light arrives at the PMT after one or more scattering interactions
in the water. Because of the additional path length, scattered photons arrive
later than direct photons, forming a broad tail in the time PDF. All functions
accept an [`LMParameters`](@ref) and a [`PMTModel`](@ref), and return the
time PDF [npe/ns] (or [npe/ns/GeV] for shower functions normalised per unit
energy).

## Muon

[`scatteredlightfrommuon`](@ref) takes the closest-approach distance `D`, the
cosine of the track direction `cd`, the PMT orientation `(θ, ϕ)`, and the
time offset `Δt`.

```@example scattered_light
using CairoMakie
using LumenManufaktur

params = LMParameters(dispersion_model = DispersionORCA)
pmt    = LumenManufaktur.KM3NeTPMT

Δts = range(-5, 200, 1000)

fig = Figure(size = (650, 380))
ax = Axis(fig[1, 1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns",
    xgridstyle = :dash, ygridstyle = :dash,
)
for D in [5, 10, 20]
    pdf = [scatteredlightfrommuon(params, pmt, Float64(D), 0.0, π/2, π/2, Δt) for Δt in Δts]
    lines!(ax, collect(Δts), pdf, label = "D = $D m")
end
axislegend(ax; position = :rt)
fig
```

## EM Shower

See [`scatteredlightfromEMshower`](@ref) for the full function signature.

```@example scattered_light
Δts2 = range(-5, 200, 1000)

fig2 = Figure(size = (650, 380))
ax2 = Axis(fig2[1, 1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns / GeV",
    xgridstyle = :dash, ygridstyle = :dash,
)
for D in [5, 10, 20]
    pdf = [scatteredlightfromEMshower(params, pmt, Float64(D), 0.0, π/2, π/2, Δt) for Δt in Δts2]
    lines!(ax2, collect(Δts2), pdf, label = "D = $D m")
end
axislegend(ax2; position = :rt)
fig2
```

## Bright Point

See [`scatteredlightfrombrightpoint`](@ref) for the full function signature.

```@example scattered_light
Δts3 = range(-5, 200, 1000)

fig3 = Figure(size = (650, 380))
ax3 = Axis(fig3[1, 1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns",
    xgridstyle = :dash, ygridstyle = :dash,
)
for D in [5, 10, 20]
    pdf = [scatteredlightfrombrightpoint(params, pmt, Float64(D), 0.0, Δt) for Δt in Δts3]
    lines!(ax3, collect(Δts3), pdf, label = "D = $D m")
end
axislegend(ax3; position = :rt)
fig3
```
