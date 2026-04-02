# Dispersion

The refractive index of deep-sea water varies with wavelength (dispersion),
which directly determines both the Cherenkov angle and the photon travel time
used in timing PDFs.  The dispersion model is pluggable via the
`dispersion_model` field of `LMParameters`.

## Built-in models

`LumenManufaktur.jl` ships the **Bailey dispersion model**, parameterised by
ambient pressure (atm), with two pre-defined instances for the KM3NeT sites:

| Constant | Pressure | Depth |
|---|---|---|
| `DispersionORCA` | 240 atm | ~2450 m (KM3NeT/ORCA) |
| `DispersionARCA` | 350 atm | ~3500 m (KM3NeT/ARCA) |

```@example dispersion
using CairoMakie
using LumenManufaktur

λs = range(300, 700, 500)

models = [
    (DispersionORCA,      "ORCA (240 atm)",   Makie.wong_colors()[1]),
    (DispersionARCA,      "ARCA (350 atm)",   Makie.wong_colors()[2]),
    (BaileyDispersion(1), "Surface (1 atm)",  Makie.wong_colors()[3]),
]

fig = Figure(size = (700, 420))
ax_ph = Axis(fig[1, 1], xlabel = "wavelength [nm]", ylabel = "phase refractive index",
             title = "Phase index", xgridstyle = :dash, ygridstyle = :dash)
ax_gr = Axis(fig[1, 2], xlabel = "wavelength [nm]", ylabel = "group refractive index",
             title = "Group index", xgridstyle = :dash, ygridstyle = :dash)

for (model, label, color) in models
    nph = [refractionindexphase(model, λ) for λ in λs]
    ngr = [refractionindexgroup(model, λ) for λ in λs]
    lines!(ax_ph, collect(λs), nph, label = label, color = color)
    lines!(ax_gr, collect(λs), ngr, label = label, color = color)
end
axislegend(ax_ph; position = :rt)
fig
```

The phase index determines the Cherenkov angle; the group index governs photon
travel time and hence the timing PDFs used in track reconstruction. Higher
pressure at depth slightly increases both indices compared to surface water.

## Custom dispersion models

Any subtype of `DispersionModel` with the appropriate methods can be used.
The minimum interface is `refractionindexphase` and `refractionindexgroup`.
For example, a simple constant-index model:

```@example dispersion
struct ConstantDispersion <: DispersionModel
    n::Float64
end

LumenManufaktur.refractionindexphase(m::ConstantDispersion, λ) = m.n
LumenManufaktur.refractionindexgroup(m::ConstantDispersion, λ) = m.n

params = LMParameters(dispersion_model = ConstantDispersion(1.38))
directlightfrommuon(params, LumenManufaktur.KM3NeTPMT, 10.0, π/2, π/2)
```
