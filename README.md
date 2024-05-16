# LumenManufaktur

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tgal.pages.km3net.de/LumenManufaktur.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tgal.pages.km3net.de/LumenManufaktur.jl/dev)
[![Build Status](https://git.km3net.de/tgal/LumenManufaktur.jl/badges/main/pipeline.svg)](https://git.km3net.de/tgal/LumenManufaktur.jl/pipelines)
[![Coverage](https://git.km3net.de/tgal/LumenManufaktur.jl/badges/main/coverage.svg)](https://git.km3net.de/tgal/LumenManufaktur.jl/commits/main)

`LumenManufaktur.jl` is a software package written in Julia with the aim of
providing high-performance routines to calculate probability density functions
of the arrival time of light on PMTs induced by muons and showers including
production (Cherenkov radiation, photons from energy loss of a muon, shower
light), transmission (dispersion, absorption and scattering) and detection
(angular acceptance, quantum efficiency, shadowing) of light with a PMT.

The software is written in a way that allows the user to plug in their own models for dispersion, absorption, scattering and also for the quantum efficiency and angular acceptance of PMTs.

The implementation is based on the C++ framework `Jpp` developed by Maarten de
Jong for the KM3NeT neutrino telescope. More details about the underlying
physics and models can be found in ["The probability density function of the
arrival time of Čerenkov light", M. de Jong, E. van Campenhout, arXiv:2305.19626
[astro-ph.IM]](https://arxiv.org/abs/2305.19626).


## Documentation

Check out the **[Latest Documention](https://tgal.pages.km3net.de/LumenManufaktur.jl/dev)**
which also includes tutorials and examples.


## Installation

`LumenManufaktur.jl` is not an officially registered Julia package but it's available via
the [KM3NeT Julia registry](https://git.km3net.de/common/julia-registry). To add
the KM3NeT Julia registry to your local Julia registry list, follow the
instructions there. After that, you can add `LumenManufaktur.jl` just like any other Julia package:

    julia> import Pkg; Pkg.add("LumenManufaktur")
    

## Quickstart

``` julia-repl
julia> using LumenManufaktur
```

## Benchmarks

Below are some benchmarks made on a MacBook Air M1 2020.

``` julia-repl
julia> using LumenManufaktur

julia> params = LMParameters(dispersion_model=BaileyDispersion(240))
LMParameters:
  minimum distance = 0.1 m
  module_radius = 0.25 m
  lambda min / max = 300.0 nm / 700.0 nm
  degree of Legendre polynomials = 5
  dispersion model = BaileyDispersion(240.0, 1.3201, 1.4e-5, 16.2566, -4383.0, 1.1455e6)
  scattering model = LumenManufaktur.Kopelevich()
  scattering probability model = LumenManufaktur.Scatteringp00075()
  absorption model = LumenManufaktur.DefaultAbsorption()
  
julia> @benchmark directlightfrommuon(params, LumenManufaktur.KM3NeTPMT, R, θ, ϕ) setup=begin; R=(rand()+1)*300; θ=rand()*2π; ϕ=rand()*2π; end
BenchmarkTools.Trial: 10000 samples with 200 evaluations.
 Range (min … max):  403.125 ns …  2.225 μs  ┊ GC (min … max): 0.00% … 77.28%
 Time  (median):     422.500 ns              ┊ GC (median):    0.00%
 Time  (mean ± σ):   428.128 ns ± 29.850 ns  ┊ GC (mean ± σ):  0.08% ±  1.08%

     ▁ ▁▄▆▆▇█▇▇▆▄▃▃▃▃▃▃▂▂▁▁▁▁▂▂▂▁▁▁                            ▂
  ▇▅▅█▆██████████████████████████████▇█▆▆▇▇▆▆▆▇▇▆▇▇▆▅▅▅▅▅▄▄▅▅▄ █
  403 ns        Histogram: log(frequency) by time       507 ns <

 Memory estimate: 64 bytes, allocs estimate: 4.

julia> @benchmark scatteredlightfrommuon(params, LumenManufaktur.KM3NeTPMT, D, cd, θ, ϕ, Δt) setup=begin; D=rand(1:100); cd=rand(); θ=rand()*2π; ϕ=rand()*2π; Δt=rand()*100; end
BenchmarkTools.Trial: 10000 samples with 7 evaluations.
 Range (min … max):  4.435 μs …  16.631 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     4.881 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   4.928 μs ± 256.698 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

                  ▅▂▇▅▅█▆▄▅▁▂
  ▂▁▂▂▂▂▂▃▃▃▄▄▅▇▇████████████▇▆▆▄▄▅▄▄▄▃▄▄▃▄▃▄▃▃▃▃▃▃▃▃▂▃▃▂▂▂▂▂ ▄
  4.43 μs         Histogram: frequency by time        5.62 μs <

 Memory estimate: 80 bytes, allocs estimate: 5.
```

Notice that the allocation estimate for a single call is not accurate, in fact,
the calculation itself is non-allocating. The example below where demonstrates
that calling the function 1000 times with random values in a tight for-loop only
requires 2 allocations in total and has a memory estimate of 192 bytes:

``` julia-repl
julia> function manycalls()
           params = LMParameters(dispersion_model=DispersionARCA)
           for i in 1:1000
               directlightfrommuon(params, LumenManufaktur.PMTKM3NeT, rand()*300+1, rand()*2π, rand()*2π)
           end
       end
manycalls (generic function with 1 method)

julia> @benchmark manycalls()
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  387.208 μs … 651.334 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     387.875 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   391.571 μs ±  11.914 μs  ┊ GC (mean ± σ):  0.00% ± 0.00%

  █▆▂▂▄▆  ▂▁            ▂▁  ▁                                   ▁
  ███████████▇▇█▇█▇▇▇▆▇▇██████▇█▇▇▆▆▆▅▆▅▅▅▆▆▆▆▅▇▇▆▆▆▅▆▅▅▃▅▅▄▃▃▄ █
  387 μs        Histogram: log(frequency) by time        425 μs <

 Memory estimate: 192 bytes, allocs estimate: 2.
```
