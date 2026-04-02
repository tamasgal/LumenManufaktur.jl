# Direct Light

Direct light arrives at the PMT via a straight-line path from the Cherenkov
emission point without scattering.  All functions return the expected number
of photo-electrons (npe) when called without a time argument, or the time PDF
[npe/ns] when a time difference `Δt` is provided.

## Muon

```@example direct_light
using CairoMakie
using LumenManufaktur

params = LMParameters(dispersion_model = DispersionORCA)
pmt    = LumenManufaktur.KM3NeTPMT

# npe vs closest-approach distance for several zenith angles
Rs = range(1, 40, 500)
θs = [0.0, π/4, π/2, 3π/4, π]

fig = Figure(size = (650, 380))
ax = Axis(fig[1, 1],
    xlabel = "closest-approach distance R [m]",
    ylabel = "npe",
    xgridstyle = :dash, ygridstyle = :dash,
)
for θ in θs
    npe = [directlightfrommuon(params, pmt, R, θ, π/2) for R in Rs]
    lines!(ax, collect(Rs), npe, label = "θ = $(round(θ, digits=2)) rad")
end
axislegend(ax; position = :rt)
fig
```

The time PDF of direct muon light (npe per ns) as a function of the time
residual `Δt`:

```@example direct_light
Δts = range(-5, 30, 500)

fig2 = Figure(size = (650, 380))
ax2 = Axis(fig2[1, 1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns",
    xgridstyle = :dash, ygridstyle = :dash,
)
for R in [5, 10, 20]
    pdf = [directlightfrommuon(params, pmt, Float64(R), π/2, π/2, Δt) for Δt in Δts]
    lines!(ax2, collect(Δts), pdf, label = "R = $R m")
end
axislegend(ax2; position = :rt)
fig2
```

## EM Shower

The 5-parameter form takes the distance `D`, the cosine of the shower axis to
PMT direction `cd`, the PMT orientation `(θ, ϕ)`, and the time offset `Δt`.
The shower profile is normalised per GeV of shower energy.

```@example direct_light
Δts3 = range(-5, 60, 600)

fig3 = Figure(size = (650, 380))
ax3 = Axis(fig3[1, 1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns / GeV",
    xgridstyle = :dash, ygridstyle = :dash,
)
for D in [5, 10, 20]
    pdf = [directlightfromEMshower(params, pmt, Float64(D), 0.0, π/2, π/2, Δt) for Δt in Δts3]
    lines!(ax3, collect(Δts3), pdf, label = "D = $D m")
end
axislegend(ax3; position = :rt)
fig3
```

## Bright Point

An isotropic bright point source (e.g. a bioluminescent flash or calibration
LED) radiates uniformly in all directions.

```@example direct_light
Δts4 = range(-5, 40, 500)

fig4 = Figure(size = (650, 380))
ax4 = Axis(fig4[1, 1],
    xlabel = "time residual Δt [ns]",
    ylabel = "npe / ns",
    xgridstyle = :dash, ygridstyle = :dash,
)
for D in [5, 10, 20]
    pdf = [directlightfrombrightpoint(params, pmt, Float64(D), 0.0, Δt) for Δt in Δts4]
    lines!(ax4, collect(Δts4), pdf, label = "D = $D m")
end
axislegend(ax4; position = :rt)
fig4
```
