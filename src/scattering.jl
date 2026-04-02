"""
Abstract base type for scattering-length models.  Subtypes must implement
`scatteringlength(::MyModel, λ)` returning the scattering length [m] for
wavelength `λ` [nm].
"""
abstract type ScatteringModel end

"""
Abstract base type for angular scattering-probability models.  Subtypes must
implement `scatteringprobability(::MyModel, cos_theta)` returning the
probability density for a scattered photon with direction cosine `cos_theta`.
"""
abstract type ScatteringProbabilityModel end

"""
Kopelevich volume scattering model for deep-sea water.

Separates contributions from pure sea water and small (<1 μm) and large
(>1 μm) particulates. Reference: C.D. Mobley, *Light and Water*, ISBN
0-12-502750-8, p. 119.
"""
struct Kopelevich <: ScatteringModel end

"""
p00075 angular scattering-probability model for deep-sea water.

A mixture of 17% Rayleigh and 83% Henyey-Greenstein (g = 0.924) scattering,
matching the KM3NeT/ANTARES p00075 parametrisation.
"""
struct Scatteringp00075 <: ScatteringProbabilityModel end

"""
    scatteringlength(λ)

Calculates the scattering length [m] for a given wavelength `λ` [ns] using the
Kopelevich model for spectral volume scattering functions. The model separates
the contributions by:

    - pure_sea:   pure sea water;
    - small_par: 'small' particles (size < 1 micro meter);
    - large_par: 'large' particles (size > 1 micro meter).

Values are taken from reference C.D. Mobley "Light and Water", ISBN 0-12-502750-8, pag. 119.

"""
function scatteringlength(::Kopelevich, λ)
    Vs = 0.0075
    Vl = 0.0075
    bw = 0.0017
    bs = 1.340
    bl = 0.312

    x = 550.0 / λ

    pure_sea = bw * x^4.3
    small_par = bs * Vs * x^1.7
    large_par = bl * Vl * x^0.3

    1.0 / (pure_sea + small_par + large_par)
end

"""
    scatteringprobability(model, cos_theta)

Probability density for a scattered photon with direction cosine `cos_theta`,
evaluated using `model`.
"""
scatteringprobability(::Scatteringp00075, λ) = p00075(λ)


"""
    p00075(x)

Model specific function to describe light scattering probability in water (p00075).

# Arguments
- `x`: cosine scattering angle
"""
@inline function p00075(x)
    g = 0.924
    f = 0.17

    f * rayleigh(x) + (1.0 - f) * henyey_greenstein(g, x)
end

"""
    henyey_greenstein(x)

Light scattering probability in water (Heneyey-Greenstein).

# Arguments
- `x`: cosine scattering angle
"""
@inline function henyey_greenstein(x)
    g = 0.924
    return henyey_greenstein(g, x)
end

"""
    henyey_greenstein(g, x)

Light scattering probability in water (Heneyey-Greenstein).

# Arguments
- `g`: angular dependence parameter
- `x`: cosine scattering angle
"""
@inline function henyey_greenstein(g, x)
    a0 = (1.0 - g^2) / (4π)
    y = 1.0 + g^2 - 2.0 * g * x

    a0 / (y * sqrt(y))
end

"""
    rayleigh(x)

Light scattering probability in water (Rayleigh).

# Arguments
- `x`: cosine scattering angle
"""
rayleigh(x) = rayleigh(0.835, x)

"""
    rayleigh(x)

Light scattering probability in water (Rayleigh).

# Arguments
- `g`: angular dependence parameter
- `x`: cosine scattering angle
"""
@inline function rayleigh(g, x)
    a0 = 1.0 / (1.0 + g / 3.0) / (4π)
    a0 * (1.0 + g * x^2)
end


"""
    inverseattenuationlength(::Scatteringp00075, l_abs, ls, cts)

Get the inverse of the attenuation length [m^-1].

# Arguments

- `l_abs`: absorption length [m]
- `ls`: scattering length [m]
- `cts`: cosine scattering angle

"""
function inverseattenuationlength(::Scatteringp00075, l_abs, ls, cts)
    1.0 / l_abs + inverseattenuationlengthinterpolator(cts) / ls
end

"""
Interpolator for the p00075 model based inverse attenutation calculation.
"""
const inverseattenuationlengthinterpolator = let
    xs = range(-1, 1; length = 100000)
    dx = xs.step.hi
    xs = collect(xs)
    ys = Float64[]
    W = 0.0
    for x in xs
        push!(ys, W)
        W += 2π * dx * scatteringprobability(Scatteringp00075(), x + 0.5 * dx)
    end
    # xs[1] = 0.0
    # xs[end] = 1.0

    LinearInterpolator(xs, ys, NoBoundaries())
end
