# Absorption

Deep-sea water absorbs light strongly outside the blue–green window. The
absorption length ``\lambda_\text{abs}(\lambda)`` is the mean free path of a
photon before absorption. Longer values mean a more transparent medium at that
wavelength.

## Default model

The built-in [`DefaultAbsorption`](@ref) model uses tabulated values from the
Jpp framework, interpolated over 290–715 nm.  It is the default when
constructing [`LMParameters`](@ref) without an explicit `absorption_model`.

```@example absorption
using CairoMakie
using LumenManufaktur

λs = range(300, 720, 500)
labs = [absorptionlength(DefaultAbsorption(), λ) for λ in λs]

fig = Figure(size = (600, 350))
ax = Axis(
    fig[1, 1],
    xlabel = "wavelength [nm]",
    ylabel = "absorption length [m]",
    xgridstyle = :dash,
    ygridstyle = :dash,
)
lines!(ax, collect(λs), labs, color = :royalblue)
fig
```

The peak transparency is around 412–440 nm (blue light), with the absorption
length exceeding 60 m. Beyond ~700 nm the medium becomes essentially opaque.

## Custom absorption models

Any struct that subtypes [`AbsorptionModel`](@ref) and implements
[`absorptionlength`](@ref) can be passed to [`LMParameters`](@ref) as a
drop-in replacement.  For example, a simple Gaussian-shaped transparency
window:

```@example absorption
struct GaussianAbsorption <: AbsorptionModel
    λ0::Float64   # peak wavelength [nm]
    σ::Float64    # width [nm]
    peak::Float64 # peak absorption length [m]
end

LumenManufaktur.absorptionlength(m::GaussianAbsorption, λ) =
    m.peak * exp(-0.5 * ((λ - m.λ0) / m.σ)^2)

my_abs = GaussianAbsorption(450.0, 60.0, 60.0)
params = LMParameters(dispersion_model = DispersionORCA, absorption_model = my_abs)

# Compare default vs custom at a few wavelengths
labs_custom = [absorptionlength(my_abs, λ) for λ in λs]

fig2 = Figure(size = (600, 350))
ax2 = Axis(fig2[1, 1], xlabel = "wavelength [nm]", ylabel = "absorption length [m]",
           xgridstyle = :dash, ygridstyle = :dash)
lines!(ax2, collect(λs), labs,        label = "DefaultAbsorption")
lines!(ax2, collect(λs), labs_custom, label = "GaussianAbsorption", linestyle = :dash)
axislegend(ax2; position = :rt)
fig2
```
