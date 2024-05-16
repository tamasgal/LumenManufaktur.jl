# LumenManufaktur.jl

`LumenManufaktur.jl` is a software package written in Julia with the aim of
providing high-performance routines to calculate probability density functions
of the arrival time of light on PMTs induced by muons and showers including
production (Cherenkov radiation, photons from energy loss of a muon, shower
light), transmission (dispersion, absorption and scattering) and detection
(angular acceptance, quantum efficiency, shadowing) of light with a PMT.

The software is written in a way that allows the user to plug in their own models for dispersion, absorption, scattering and also for the quantum efficiency and angular acceptance of PMTs.

The implementation is based on the C++ framework `Jpp` developed by Maarten de
Jong for the KM3NeT neutrino telescope. More details about the underlying
physics and models can be found in **"The probability density function of the
arrival time of Čerenkov light", M. de Jong, E. van Campenhout** [arXiv:2305.19626
[astro-ph.IM]](https://arxiv.org/abs/2305.19626).

## Installation

`LumenManufaktur.jl` is **not an officially registered Julia package** but it's available on
the **[KM3NeT Julia registry](https://git.km3net.de/common/julia-registry)**. To add
the KM3NeT Julia registry to your local Julia registry list, follow the
instructions there. After that, you can install `LumenManufaktur.jl` just like any other Julia package:

    julia> import Pkg; Pkg.add("LumenManufaktur")
    
## Quickstart

```@example quickstart
using LumenManufaktur

params = LMParameters(dispersion_model=DispersionORCA)
directlightfrommuon(params, LumenManufaktur.KM3NeTPMT, 23.5, π, π/2)
```
