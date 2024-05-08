"""
The parameter set for light detection.

# Arguments
- `minimum_distance`: the minimum distance [m] between muon track and PMT
- `lambda_min`: minimum wavelength [ns]
- `lambda_max`: maximum wavelength [ns]
- `legendre_coefficients`: a tuple of two vectors which contain the Legendre coefficients
- `dispersion_model`: the dispersion model
- `scattering_model`: the scattering model
- `absorption_model`: the absorption model
"""
Base.@kwdef struct LMParameters{D<:DispersionModel,S<:ScatteringModel,A<:AbsorptionModel}
    minimum_distance::Float64 = 1.0e-1
    lambda_min::Float64 = 300.0
    lambda_max::Float64 = 700.0
    legendre_coefficients::Tuple{Vector{Float64},Vector{Float64}} = gausslegendre(5)
    dispersion_model::D = BaileyDispersion()
    scattering_model::S = Kopelevich()
    absorption_model::A = DefaultAbsorption()
end
function Base.show(io::IO, p::LMParameters)
    println(io, "LMParameters:")
    println(io, "  minimum distance = $(p.minimum_distance)")
    println(io, "  lambda min/max = $(p.lambda_min)/$(p.lambda_max)")
    println(
        io,
        "  degree of Legendre polynomials = $(length(first(p.legendre_coefficients)))",
    )
    println(io, "  dispersion model = $(p.dispersion_model)")
    println(io, "  scattering model = $(p.scattering_model)")
    print(io, "  absorption model = $(p.absorption_model)")
end
