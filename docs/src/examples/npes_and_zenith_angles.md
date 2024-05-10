# npe vs zenith angle

```@example zenith2
using CairoMakie
using LumenManufaktur

fig = Figure(size = (600, 400));
ax = Axis(
	fig[1, 1],
	xlabel = "PMT zenith angle [rad]",
	ylabel = "npe",
	xgridstyle = :dash, 
	ygridstyle = :dash
)
θs = range(-π, π, 1000)
params = LMParameters(dispersion_model=BaileyDispersion(240))
pmt = LumenManufaktur.KM3NeTPMT
for R in 2:10
	lines!(ax, θs, [directlightfrommuon(params, pmt, R, θ, π/2) for θ in θs], label="R = $R m")
end
axislegend(; position = :ct)
fig
```

```@example zenith1
using CairoMakie
using LumenManufaktur

fig = Figure(size = (600, 400));
ax = Axis(
	fig[1, 1],
	xlabel = "distance between muon and PMT / m",
	ylabel = "npe",
	xgridstyle = :dash, 
	ygridstyle = :dash
)
Rs = range(1, 10, 1000)
params = LMParameters(dispersion_model=BaileyDispersion(240))
pmt = LumenManufaktur.KM3NeTPMT
for theta_deg in 0:30:180
    lines!(ax, Rs, [directlightfrommuon(params, pmt, R, deg2rad(theta_deg), π/2) for R in Rs], label="θ = $(theta_deg) deg")
end
axislegend(; position = :rt)
fig
```

