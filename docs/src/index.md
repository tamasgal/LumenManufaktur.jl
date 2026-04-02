# LumenManufaktur.jl

`LumenManufaktur.jl` is a Julia package for computing probability density
functions of Cherenkov-light arrival times on PMTs from muons and
electromagnetic showers, covering the full chain of light production
(Cherenkov radiation, delta rays, EM showers), transmission (dispersion,
absorption, scattering), and detection (quantum efficiency, angular
acceptance, module shadowing).

The implementation is based on the C++ framework `Jpp` developed by Maarten
de Jong for the KM3NeT neutrino telescope. More details about the underlying
physics can be found in **"The probability density function of the arrival
time of Čerenkov light", M. de Jong, E. van Campenhout**
[arXiv:2305.19626 [astro-ph.IM]](https://arxiv.org/abs/2305.19626).

## Modular Design

The key advantage of `LumenManufaktur.jl` over the original Jpp implementation
is its **modular, type-based architecture**: every physics ingredient is a
distinct type implementing a small interface, so any component can be replaced
with a drop-in custom implementation without touching the core code.

| Component | Abstract type | Built-in defaults |
|---|---|---|
| Dispersion | `DispersionModel` | `BaileyDispersion`, `DispersionORCA`, `DispersionARCA` |
| Scattering length | `ScatteringModel` | `Kopelevich` |
| Scattering probability | `ScatteringProbabilityModel` | `Scatteringp00075` |
| Absorption | `AbsorptionModel` | `DefaultAbsorption` |

All models are passed through the central `LMParameters` struct.  Swapping one
is a one-liner:

```@example modularity
using LumenManufaktur

# Use ORCA pressure dispersion with everything else at defaults
params_orca = LMParameters(dispersion_model = DispersionORCA)
```

Custom models require only a type declaration and a single method:

```@example modularity
# A flat (wavelength-independent) absorption length for illustration
struct FlatAbsorption <: AbsorptionModel end
LumenManufaktur.absorptionlength(::FlatAbsorption, λ) = 55.0   # [m], constant across wavelengths

params_flat = LMParameters(
    dispersion_model = DispersionORCA,
    absorption_model = FlatAbsorption(),
)
directlightfrommuon(params_flat, LumenManufaktur.KM3NeTPMT, 10.0, π/2, π/2)
```

The PMT is equally flexible: `PMTModel` accepts any callable for the quantum
efficiency and angular acceptance, so detector-specific response curves can be
injected directly without modifying library code.

## Installation

`LumenManufaktur.jl` is **not an officially registered Julia package** but it
is available on the
**[KM3NeT Julia registry](https://git.km3net.de/common/julia-registry)**.
Follow the instructions there to add the registry, then install as usual:

    julia> import Pkg; Pkg.add("LumenManufaktur")

## Quickstart

```@example quickstart
using LumenManufaktur

params = LMParameters(dispersion_model = DispersionORCA)
directlightfrommuon(params, LumenManufaktur.KM3NeTPMT, 23.5, π, π/2)
```
