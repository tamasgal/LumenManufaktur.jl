abstract type ScatteringModel end

struct Kopelevich <: ScatteringModel end

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

    pure_sea  = bw *      x^4.3
    small_par = bs * Vs * x^1.7
    large_par = bl * Vl * x^0.3

    1.0 / (pure_sea + small_par + large_par)
end
