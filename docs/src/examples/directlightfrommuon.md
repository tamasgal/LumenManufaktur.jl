# Direct light from muon

The following example demonstrates how to calculate the number of photoelectrons
for a given environment and PMT using [`directlightfrommuon`](@ref).

```@example directlightfrommuon
using LumenManufaktur
```

Let's create the [`LMParameters`](@ref) with the defaults but override the
dispersion model with [`DispersionORCA`](@ref) for KM3NeT ORCA detectors:

```@example directlightfrommuon
params = LMParameters(dispersion_model=DispersionORCA)
```

To calculate the number of photoelectrons from direct light induced by a muon at
the closest distance of `R=23.5`, with a standard KM3NeT PMT pointing in the
direction of `θ=π` and `ϕ=π/2`, we call [`directlightfrommuon`](@ref):

```@example directlightfrommuon
directlightfrommuon(params, LumenManufaktur.KM3NeTPMT, 23.5, π, π/2)
```
