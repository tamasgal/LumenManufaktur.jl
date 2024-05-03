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

Here are some benchmarks made on a MacBook Air M1 2020:

``` julia-repl
julia> using MuonLight

julia> dp = MuonLight.DispersionARCA
BasicDispersion(350.0, 1.3201, 1.4e-5, 16.2566, -4383.0, 1.1455e6)

julia> @benchmark MuonLight.directlight($dp, R, θ, ϕ) setup=begin; R=(rand()+1)*300; θ=rand()*2π; ϕ=rand()*2π; end
BenchmarkTools.Trial: 10000 samples with 11 evaluations.
 Range (min … max):  1.004 μs …  37.780 μs  ┊ GC (min … max): 0.00% … 92.37%
 Time  (median):     1.042 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   1.062 μs ± 521.187 ns  ┊ GC (mean ± σ):  0.66% ±  1.30%

         ▄██▇ ▄▁                                               
  ▂▂▂▃▅▁█████▁███▆▆▁▆▇▆▅▅▁▄▃▃▃▃▁▄▄▄▃▃▁▃▃▂▂▂▁▂▂▂▂▂▁▂▂▂▂▂▁▂▂▂▂▂ ▃
  1 μs            Histogram: frequency by time        1.19 μs <

 Memory estimate: 1.12 KiB, allocs estimate: 62.

```
