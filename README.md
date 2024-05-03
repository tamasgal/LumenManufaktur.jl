# MuonLight

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tgal.pages.km3net.de/MuonLight.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tgal.pages.km3net.de/MuonLight.jl/dev)
[![Build Status](https://git.km3net.de/tgal/MuonLight.jl/badges/main/pipeline.svg)](https://git.km3net.de/tgal/MuonLight.jl/pipelines)
[![Coverage](https://git.km3net.de/tgal/MuonLight.jl/badges/main/coverage.svg)](https://git.km3net.de/tgal/MuonLight.jl/commits/main)

Welcome to the `MuonLight.jl` repository!


## Documentation

Check out the **[Latest Documention](https://tgal.pages.km3net.de/MuonLight.jl/dev)**
which also includes tutorials and examples.


## Installation

`MuonLight.jl` is not an officially registered Julia package but it's available via
the [KM3NeT Julia registry](https://git.km3net.de/common/julia-registry). To add
the KM3NeT Julia registry to your local Julia registry list, follow the
instructions in its
[README](https://git.km3net.de/common/julia-registry#adding-the-registry) or simply do

    git clone https://git.km3net.de/common/julia-registry ~/.julia/registries/KM3NeT
    
After that, you can add `MuonLight.jl` just like any other Julia package:

    julia> import Pkg; Pkg.add("MuonLight")
    

## Quickstart

``` julia-repl
julia> using MuonLight
```

## Benchmarks

Below are some benchmarks made on a MacBook Air M1 2020.

``` julia-repl
julia> using MuonLight

julia> params = MuonLight.Parameters(dispersion_model=MuonLight.DispersionARCA)
Parameters:
  minimum distance = 0.1
  lambda min/max = 300.0/700.0
  degree of Legendre polynomials = 5
  dispersion model = BasicDispersion(350.0, 1.3201, 1.4e-5, 16.2566, -4383.0, 1.1455e6)
  scattering model = MuonLight.Kopelevich()
  absorption model = MuonLight.DefaultAbsorption()

julia> @benchmark MuonLight.directlight($params, $MuonLight.PMTKM3NeT, R, θ, ϕ) setup=begin; R=(rand()+1)*300; θ=rand()*2π; ϕ=rand()*2π; end
BenchmarkTools.Trial: 10000 samples with 198 evaluations.
 Range (min … max):  434.131 ns …  2.985 μs  ┊ GC (min … max): 0.00% … 81.73%
 Time  (median):     445.076 ns              ┊ GC (median):    0.00%
 Time  (mean ± σ):   451.110 ns ± 59.088 ns  ┊ GC (mean ± σ):  0.28% ±  1.93%

       ▆█▂                                                      
  ▂▃▄▄████▅▄▃▄▃▄▄▄▃▃▂▂▂▂▂▂▂▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▂▂▂▁▁▂▂▂▂▁▂▂▂▂▂▂▂ ▃
  434 ns          Histogram: frequency by time          532 ns <

 Memory estimate: 160 bytes, allocs estimate: 5.
```

Notice that the allocation estimate for a single call is not accurate, in fact,
the calculation itself is non-allocating. The example below where demonstrates
that calling the function 1000 times with random values in a tight for-loop only
requires 2 allocations in total and has a memory estimate of 192 bytes:

``` julia-repl
julia> function manycalls()
           params = MuonLight.Parameters(dispersion_model=MuonLight.DispersionARCA)
           for i in 1:1000
               MuonLight.directlight(params, MuonLight.PMTKM3NeT, rand()*300+1, rand()*2π, rand()*2π)
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
