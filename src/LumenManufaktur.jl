module LumenManufaktur

using BasicInterpolators
using FastGaussQuadrature

const C = 0.299792458  # Speed of light in vacuum [m/ns]
const α = 1.0 / 137.036  # Fine-structure constant
@inline tanthetac(n) = √((n - 1.0) * (n + 1.0))  # Average tangent corresponding to the group velocity

include("dispersion.jl")
include("scattering.jl")
include("absorption.jl")
include("pmt.jl")
include("parameters.jl")
include("utils.jl")
include("muons.jl")
include("showers.jl")
include("deltarays.jl")
include("exports.jl")

end
