"""
The parameter set for light detection.

# Arguments
- `minimum_distance`: the minimum distance [m] between muon track and PMT
- `module_radius`: radius of the optical module [m] used to implement shadowing of the PMT by the optical module
- `lambda_min`: minimum wavelength [ns]
- `lambda_max`: maximum wavelength [ns]
- `n`: average index of refraction of medium (default: water) corresponding to the group velocity
- `legendre_coefficients`: a tuple of two vectors which contain the Legendre coefficients
- `dispersion_model`: the dispersion model (default: Bailey at 1 atm pressure)
- `scattering_model`: the scattering model (default: Kopelevich)
- `scattering_probability_model`: the scattering probability model (default: p00075)
- `absorption_model`: the absorption model
"""
Base.@kwdef struct LMParameters{
    D<:DispersionModel,
    S<:ScatteringModel,
    SP<:ScatteringProbabilityModel,
    A<:AbsorptionModel,
}
    minimum_distance::Float64 = 1.0e-1
    module_radius::Float64 = 0.25
    lambda_min::Float64 = 300.0
    lambda_max::Float64 = 700.0
    n::Float64 = 1.3800851282
    legendre_coefficients::Tuple{Vector{Float64},Vector{Float64}} = gausslegendre(5)
    dispersion_model::D = BaileyDispersion()
    scattering_model::S = Kopelevich()
    scattering_probability_model::SP = Scatteringp00075()
    absorption_model::A = DefaultAbsorption()
end
function Base.show(io::IO, p::LMParameters)
    println(io, "LMParameters:")
    println(io, "  minimum distance = $(p.minimum_distance) m")
    println(io, "  module_radius = $(p.module_radius) m")
    println(io, "  lambda min / max = $(p.lambda_min) nm / $(p.lambda_max) nm")
    println(
        io,
        "  degree of Legendre polynomials = $(length(first(p.legendre_coefficients)))",
    )
    println(io, "  dispersion model = $(p.dispersion_model)")
    println(io, "  scattering model = $(p.scattering_model)")
    println(io, "  scattering probability model = $(p.scattering_probability_model)")
    print(io, "  absorption model = $(p.absorption_model)")
end
