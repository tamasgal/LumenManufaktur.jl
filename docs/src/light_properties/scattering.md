# Scattering

Scattering in the medium is governed by two independent, swappable components
in `LMParameters`:

| Parameter | Abstract type | Purpose |
|---|---|---|
| `scattering_model` | `ScatteringModel` | Mean free path between scattering events |
| `scattering_probability_model` | `ScatteringProbabilityModel` | Angular distribution after a scatter |

## Scattering length

The scattering length is the mean free path between interactions.  The
built-in `Kopelevich` model separates contributions from pure sea water and
from small (<1 μm) and large (>1 μm) particulates, following the spectral
volume scattering functions of Mobley (1994).

```@example scattering
using CairoMakie
using LumenManufaktur

λs = range(300, 700, 500)
ls = [scatteringlength(Kopelevich(), λ) for λ in λs]

fig = Figure(size = (600, 350))
ax = Axis(
    fig[1, 1],
    xlabel = "wavelength [nm]",
    ylabel = "scattering length [m]",
    xgridstyle = :dash,
    ygridstyle = :dash,
)
lines!(ax, collect(λs), ls, color = :teal, label = "Kopelevich")
axislegend(ax; position = :rt)
fig
```

A custom scattering model only needs to subtype `ScatteringModel` and
implement `scatteringlength`:

```@example scattering
# Power-law scattering as a minimal custom example
struct PowerLawScattering <: ScatteringModel
    ls0::Float64  # scattering length at 550 nm [m]
    α::Float64    # spectral slope
end

LumenManufaktur.scatteringlength(m::PowerLawScattering, λ) = m.ls0 * (λ / 550.0)^m.α

ls_custom = [scatteringlength(PowerLawScattering(50.0, 1.0), λ) for λ in λs]

fig2 = Figure(size = (600, 350))
ax2 = Axis(fig2[1, 1], xlabel = "wavelength [nm]", ylabel = "scattering length [m]",
           xgridstyle = :dash, ygridstyle = :dash)
lines!(ax2, collect(λs), ls,        label = "Kopelevich (default)")
lines!(ax2, collect(λs), ls_custom, label = "PowerLawScattering", linestyle = :dash)
axislegend(ax2; position = :rb)
fig2
```

## Scattering probability

The scattering probability model describes the angular distribution of a
scattered photon. The built-in `Scatteringp00075` model is a mixture of
Rayleigh (17%) and Henyey-Greenstein (83%, *g* = 0.924):

```@example scattering
xs  = range(-1, 1, 1000)
hg  = [henyey_greenstein(x) for x in xs]
ray = [rayleigh(x) for x in xs]
mix = [scatteringprobability(Scatteringp00075(), x) for x in xs]

fig3 = Figure(size = (600, 350))
ax3 = Axis(fig3[1, 1], xlabel = "cos θ", ylabel = "probability density",
           xgridstyle = :dash, ygridstyle = :dash)
lines!(ax3, collect(xs), hg,  label = "Henyey-Greenstein (g = 0.924)")
lines!(ax3, collect(xs), ray, label = "Rayleigh")
lines!(ax3, collect(xs), mix, label = "p00075 mixture (default)", linestyle = :dash, linewidth = 2)
axislegend(ax3; position = :lt)
fig3
```

The strong forward peak means most scattered photons deviate only slightly
from their original direction. A custom probability model follows the same
pattern — subtype `ScatteringProbabilityModel` and implement
`scatteringprobability(::MyModel, cos_theta)`.
