abstract type DispersionModel end


"""

Light dispersion model for water in deep sea based on David J.L. Bailey's work:
"Monte Carlo tools and analysis methods for understanding the ANTARES experiment
and predicting its sensitivity to Dark Matter" - PhD thesis, University of
Oxford, United Kingdom, 2002.

The Bailey dispersion interface is based on the implementation in the Jpp
framework by Maarten de Jong.

"""
Base.@kwdef struct BaileyDispersion <: DispersionModel
    P::Float64 = 1                 #  ambient pressure [atm]
    a0::Float64 = 1.3201           #  offset
    a1::Float64 = 1.4e-5           #  dn/dP
    a2::Float64 = 16.2566          #  d^1n/(dx)^1
    a3::Float64 = -4383.0          #  d^2n/(dx)^2
    a4::Float64 = 1.1455e6         #  d^3n/(dx)^3
end

BaileyDispersion(P) = BaileyDispersion(P = P)
const DispersionORCA = BaileyDispersion(240)
const DispersionARCA = BaileyDispersion(350)


@inline function refractionindexphase(dp::BaileyDispersion, λ)
    x = 1.0 / λ
    dp.a0 + dp.a1 * dp.P + x * (dp.a2 + x * (dp.a3 + x * dp.a4))
end

@inline function refractionindexgroup(dp::BaileyDispersion, λ)
    n = refractionindexphase(dp, λ)
    y = dispersionphase(dp, λ)
    n / (1.0 + y*λ/n)
end

@inline function dispersionphase(dp::BaileyDispersion, λ)
    x = 1.0 / λ
    -(x^2)*(dp.a2 + x*(2.0*dp.a3 + x*3.0*dp.a4))
end

@inline function dispersiongroup(dp::BaileyDispersion, λ)
    x = 1.0 / λ

    n   = refractionindexphase(dp, λ)
    np  = dispersionphase(dp, λ)
    npp = x^3*(2.0*dp.a2 + x*(6.0*dp.a3 + x*12.0*dp.a4))
    ng  = n / (1.0 + np*λ/n)

    ng^2 * (2*np^2 - n*npp) * λ / (n^3);
end

"""
Determine wavelength [nm] for a given index of refraction corresponding to the
group velocity. The estimate of the wavelength is made by successive linear
extrapolations. The procedure starts from the given wavelength and
terminates if the index of refraction is equal to the target value within
the given precision.

# Arguments

- `n`: of refraction
- `w`: value wavelength [nm]
- `eps`: precision index of refraction

"""
function wavelength(dp::BaileyDispersion, n, w, eps)
    v = w

    while true
      y = refractionindexgroup(dp, v)
      abs(y - n) < eps && break
      v += (n - y) / dispersiongroup(dp, v)
    end

    v
end
