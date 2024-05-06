abstract type AbsorptionModel end

struct DefaultAbsorption <: AbsorptionModel end


const absorptionlengthinterpolator = LinearInterpolator(
    [0, 290, 310, 330, 350, 375, 412, 440, 475, 488, 510, 532, 555, 650, 676, 715, 720, 999999],
    [0.0, 0.0, 11.9, 16.4, 20.6, 29.5, 48.5, 67.5, 59.0, 55.1, 26.1, 19.9, 14.7, 2.8, 2.3, 1.0, 0.0, 0.0],
    NoBoundaries()
)
"""
    absorptionlength(λ)

Returns the absorption length [m] in deep sea water for a given wavelength [nm].
Reference missing! The interpolation values are taken from the Jpp framework,
written by Maarten de Jong.

"""
absorptionlength(::DefaultAbsorption, λ) = absorptionlengthinterpolator(λ)
