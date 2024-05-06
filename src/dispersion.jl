abstract type DispersionModel end


"""

Light dispersion model for water in deep sea based on David J.L. Bailey's work:
"Monte Carlo tools and analysis methods for understanding the ANTARES experiment
and predicting its sensitivity to Dark Matter" - PhD thesis, University of
Oxford, United Kingdom, 2002.

"""
Base.@kwdef struct BaileyDispersion <: DispersionModel
    P::Float64 = 1                 #  ambient pressure [atm]
    a0::Float64 = 1.3201           #  offset
    a1::Float64 = 1.4e-5           #  dn/dP
    a2::Float64 = 16.2566          #  d^1n/(dx)^1
    a3::Float64 = -4383.0          #  d^2n/(dx)^2
    a4::Float64 = 1.1455e6         #  d^3n/(dx)^3
end

BaileyDispersion(P) = BaileyDispersion(P=P)
const DispersionORCA = BaileyDispersion(240)
const DispersionARCA = BaileyDispersion(350)


@inline function refractionindexphase(dp::BaileyDispersion, λ)
    x = 1.0 / λ
    dp.a0  +  dp.a1*dp.P  +  x*(dp.a2 + x*(dp.a3 + x*dp.a4))
end


@inline function refractionindexgroup(dp::BaileyDispersion, λ)
    error("Not implemented yet")
end
