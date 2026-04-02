# Combined Light

The convenience functions `lightfrommuon`, `lightfromEMshower`, and
`lightfrombrightpoint` sum the direct and scattered contributions (and, for
muons, the delta-ray and track-integrated EM-shower contributions) into a
single time PDF.

## Muon (direct + scattered + delta rays + EM showers)

`lightfrommuon(params, pmt, E, R, θ, ϕ, Δt)` sums six contributions:

- direct light from the muon track
- scattered light from the muon track
- direct + scattered light from track-integrated EM showers (scaled by energy `E`)
- direct + scattered light from delta rays (scaled by `deltarayenergyloss(E)`)

```@example combined_light
using CairoMakie
using LumenManufaktur

params = LMParameters(dispersion_model = DispersionORCA)
pmt    = LumenManufaktur.KM3NeTPMT

Δts = range(-5, 200, 1000)
E   = 100.0   # muon energy [GeV]
R   = 10.0    # closest-approach distance [m]

direct_muon   = [directlightfrommuon(params, pmt, R, π/2, π/2, Δt)            for Δt in Δts]
scattered_muon = [scatteredlightfrommuon(params, pmt, R, π/2, π/2, Δt)        for Δt in Δts]
total         = [lightfrommuon(params, pmt, E, R, π/2, π/2, Δt)               for Δt in Δts]

fig = Figure(size = (700, 400))
ax = Axis(fig[1, 1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns",
    title = "Muon light at R = $(Int(R)) m, E = $(Int(E)) GeV",
    xgridstyle = :dash, ygridstyle = :dash,
)
lines!(ax, collect(Δts), direct_muon,    label = "direct (muon)")
lines!(ax, collect(Δts), scattered_muon, label = "scattered (muon)")
lines!(ax, collect(Δts), total,          label = "total (lightfrommuon)", linestyle = :dash, linewidth = 2)
axislegend(ax; position = :rt)
fig
```

## EM Shower

`lightfromEMshower` is available in two forms:
- 5-parameter: returns `d²P/dt/dE` [npe/ns/GeV], independent of shower energy
- 6-parameter: includes shower energy `E` [GeV] and returns `dP/dt` [npe/ns]

```@example combined_light
Δts2 = range(-5, 100, 800)

fig2 = Figure(size = (700, 400))
ax2 = Axis(fig2[1, 1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns",
    title = "EM shower light at D = 10 m",
    xgridstyle = :dash, ygridstyle = :dash,
)
for E_sh in [10.0, 50.0, 100.0]
    pdf = [lightfromEMshower(params, pmt, E_sh, 10.0, 0.0, π/2, π/2, Δt) for Δt in Δts2]
    lines!(ax2, collect(Δts2), pdf, label = "E = $(Int(E_sh)) GeV")
end
axislegend(ax2; position = :rt)
fig2
```

## Bright Point

`lightfrombrightpoint(params, pmt, D, ct, Δt)` combines direct and scattered
light from an isotropic point source at distance `D`.

```@example combined_light
Δts3 = range(-5, 150, 800)

fig3 = Figure(size = (700, 400))
ax3 = Axis(fig3[1, 1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns",
    title = "Bright-point light",
    xgridstyle = :dash, ygridstyle = :dash,
)
for D in [5, 10, 20]
    pdf = [lightfrombrightpoint(params, pmt, Float64(D), 0.0, Δt) for Δt in Δts3]
    lines!(ax3, collect(Δts3), pdf, label = "D = $D m")
end
axislegend(ax3; position = :rt)
fig3
```
