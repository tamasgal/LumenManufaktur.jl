# EM-shower photon emission profile parameters (from JGeanx in Jpp).
# P(cos θ) ∝ exp(b × |cos θ − 1/n|^a), normalised to unit solid angle.
#
# Reference: R. Mirani, "Parametrisation of EM-showers in the ANTARES detector volume",
#            Doctoral thesis in computational physics, University of Amsterdam.
const _GEANX_A = 0.35
const _GEANX_B = -5.40


# Integral of exp(b × |ct − 1/n|^a) dct from xmin to xmax.
#
# Derived by substituting t = (−b) × |ct − 1/n|^a, which maps the integrand to a
# lower incomplete gamma: ∫ exp(−t) × (1/a)(1/(−b))^(1/a) t^(1/a−1) dt = (1/a)(1/(−b))^(1/a) γ(1/a, t).
# Splitting at the Cherenkov cone ct = 1/n and applying sign conventions for the two
# half-integrals follows the implementation in JGeanx::evaluate (Jpp/JPhysics/JGeanx.hh).
function _geanx_integral(n::AbstractFloat, xmin::AbstractFloat, xmax::AbstractFloat)
    x  = 1.0 / n
    ai = 1.0 / _GEANX_A

    xl = (-_GEANX_B) * abs(x - xmin)^_GEANX_A
    xr = (-_GEANX_B) * abs(x - xmax)^_GEANX_A

    gp = exp(_lgamma(ai))           # Γ(1/a)
    yl = _gammainc(ai, xl) * gp    # lower incomplete gamma γ(1/a, xl)
    yr = _gammainc(ai, xr) * gp

    xmin > x && (yl = -yl)
    xmax < x && (yr = -yr)

    return ai * (1.0 / (-_GEANX_B))^ai * (yl + yr)
end


"""
    geant(n, ct)

Number of photons from an EM-shower per unit solid angle (d²N/dcos θ/dφ) as a function
of the phase index of refraction `n` and the cosine of the emission angle `ct`.

Parametrisation: c × exp(b × |ct − 1/n|^a), normalised so 2π ∫₋₁¹ P dcos = 1.

Reference: R. Mirani, "Parametrisation of EM-showers in the ANTARES detector volume",
           Doctoral thesis, University of Amsterdam.
"""
function geant(n::AbstractFloat, ct::AbstractFloat)
    norm = _geanx_integral(n, -1.0, 1.0)
    c = 1.0 / (2π * norm)
    return c * exp(_GEANX_B * abs(ct - 1.0 / n)^_GEANX_A)
end


"""
    geant(n, ct_min, ct_max)

Integral of the EM-shower photon emission profile over cos θ ∈ `[ct_min, ct_max]` (= dN/dφ).
Used to integrate over a small angular bin in scattered light calculations.
"""
function geant(n::AbstractFloat, ct_min::AbstractFloat, ct_max::AbstractFloat)
    norm = _geanx_integral(n, -1.0, 1.0)
    c = 1.0 / (2π * norm)
    return c * _geanx_integral(n, ct_min, ct_max)
end
