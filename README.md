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
