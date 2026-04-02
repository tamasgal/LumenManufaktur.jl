# Log-gamma function using the Lanczos approximation.
#
# Reference: Numerical Recipes in C++, W.H. Press, S.A. Teukolsky, W.T. Vetterling
#            and B.P. Flannery, Cambridge University Press, §6.1.
const _LGAMMA_COEFF = (
    76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 1.208650973866179e-3, -5.395239384953e-6,
)

function _lgamma(xx::Float64)::Float64
    x = xx
    y = xx
    tmp = x + 5.5
    tmp = (x + 0.5) * log(tmp) - tmp
    ser = 1.000000000190015
    for c in _LGAMMA_COEFF
        y += 1.0
        ser += c / y
    end
    return tmp + log(2.5066282746310005 * ser / x)  # 2.506... = sqrt(2π)
end


# Regularized lower incomplete gamma function P(a, x) = γ(a, x) / Γ(a).
#
# Computed via:
#  - Series expansion (§6.2, Numerical Recipes) for x < a + 1
#  - Continued fraction via modified Lentz's method (§6.2, Numerical Recipes) for x ≥ a + 1
#
# Reference: Numerical Recipes in C++, W.H. Press, S.A. Teukolsky, W.T. Vetterling
#            and B.P. Flannery, Cambridge University Press.
function _gammainc(a::Float64, x::Float64)::Float64
    x == 0.0 && return 0.0
    gln = _lgamma(a)
    if x < a + 1.0
        ap  = a
        del = 1.0 / a
        s   = del
        for _ in 1:1_000_000
            ap  += 1.0
            del *= x / ap
            s   += del
            abs(del) < abs(s) * eps(Float64) && return s * exp(-x + a * log(x) - gln)
        end
        error("_gammainc series did not converge for a=$a, x=$x")
    else
        FPMIN = floatmin(Float64) / eps(Float64)
        b = x + 1.0 - a
        c = 1.0 / FPMIN
        d = 1.0 / b
        h = d
        for i in 1:1_000_000
            an = -i * (i - a)
            b  += 2.0
            d   = an * d + b
            abs(d) < FPMIN && (d = FPMIN)
            c   = b + an / c
            abs(c) < FPMIN && (c = FPMIN)
            d   = 1.0 / d
            del = d * c
            h  *= del
            abs(del - 1.0) < eps(Float64) && return 1.0 - exp(-x + a * log(x) - gln) * h
        end
        error("_gammainc continued fraction did not converge for a=$a, x=$x")
    end
end


"""
Longitudinal emission profile of an EM-shower.

    P(z, E) ∝ z^(a-1) × exp(-z/b),   a = a₀ + a₁ × ln(E)

# Fields
- `a0`: power term constant
- `a1`: power term energy dependence
- `b`: exponential slope [m]
- `Emin`: energy below which the shower is treated as point-like

Reference: C. Kopper, "Performance Studies for the KM3NeT Neutrino Telescope",
           PhD thesis, University of Erlangen.
"""
struct Geanz
    a0::Float64
    a1::Float64
    b::Float64
    Emin::Float64
    Geanz(a0, a1, b) = new(a0, a1, b, exp(-a0 / a1))
end

const _MINIMAL_SHOWER_SIZE = 1.0e-6  # [m]

"""
    getmaximum(gz::Geanz, E)

Depth of shower maximum [m] for shower energy `E` [GeV].
"""
@inline function getmaximum(gz::Geanz, E)
    a = gz.a0 + gz.a1 * log(E)
    (a - 1.0) * gz.b
end


"""
    getprobability(gz::Geanz, E, z)

Probability density dP/dz of photon emission at depth `z` [m] for shower energy `E` [GeV].
"""
function getprobability(gz::Geanz, E, z)
    if E > gz.Emin
        a = gz.a0 + gz.a1 * log(E)
        return z^(a - 1.0) * exp(-z / gz.b) / (gz.b^a * exp(_lgamma(a)))
    end
    return z <= _MINIMAL_SHOWER_SIZE ? 1.0 / _MINIMAL_SHOWER_SIZE : 0.0
end


"""
    getintegral(gz::Geanz, E, z)

Integral of the shower profile from 0 to `z` [m] for energy `E` [GeV].
Returns the regularized lower incomplete gamma P(a, z/b).
"""
function getintegral(gz::Geanz, E, z)
    if E > gz.Emin
        a = gz.a0 + gz.a1 * log(E)
        return _gammainc(a, z / gz.b)
    end
    return z <= _MINIMAL_SHOWER_SIZE ? z / _MINIMAL_SHOWER_SIZE : 1.0
end


"""
    getlength(gz::Geanz, E, P, prec=1e-5)

Shower length [m] for integrated probability `P` at energy `E` [GeV].
Uses bisection / expansion search.
"""
function getlength(gz::Geanz, E, P, prec=1.0e-5)
    E <= gz.Emin && return 0.0
    z0 = 0.0
    z1 = getmaximum(gz, E)
    for _ in 1:1000
        p = getintegral(gz, E, z1)
        abs(p - P) < P * prec && return z1
        if p > P
            z1 = 0.5 * (z0 + z1)
        else
            z0 = z1
            z1 *= 1.5
        end
    end
    return z1
end


"""
Default longitudinal EM-shower profile for KM3NeT/ANTARES.
"""
const geanz = Geanz(1.85, 0.62, 0.54)
