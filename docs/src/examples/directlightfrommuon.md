# Direct light from muon

The following example demonstrates how to calculate the number of photoelectrons
for a given environment and PMT.

```@example directlightfrommuon
using LumenManufaktur
```

Let's create the parameters with the defaults but override the dispersion model
with the one used in KM3NeT ORCA detectors in the Mediterranean:

```@example directlightfrommuon
params = LMParameters(dispersion_model=DispersionORCA)
```

To calculate the number of photoelectrons from direct light induced by a muon at
the closest distance of `R=23.5`, with a standard KM3NeT PMT pointing in the
direction of `θ=π` and `ϕ=π/2`, we call `directlightfrommuon()`:

```@example directlightfrommuon
directlight(params, LumenManufaktur.KM3NeTPMT, 23.5, π, π/2)
```
