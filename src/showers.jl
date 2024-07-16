"""

Correction of energy determined by JShowerEnergy.cc in Jpp.

"""
struct EMShowerCorrection{D<:DispersionModel}
    lambda_min::Float64
    lambda_max::Float64
    nx::Int
    dispersion_model::D
    f1
    f2
    function EMShowerCorrection(dispersion_model::DispersionMdoel; lambda_min=300, lambda_max=700, nx=10000)
        xmin = 1.0 / lambda_max
        xmax = 1.0 / lambda_min
        value = 0.0

        for x in xmin:(xmax-xmin)/nx:xmax

            w  = 1.0 / x
            dw = dx * w*w

            n = refractionindexphase(dispersion_model, w)

            value += cherenkov(w, n) * dw
        end

        value *= geanc()      # number of photons per GeV

        f1 = LinearInterpolator(
            [
                0.510998946e-3,  # electron mass [GeV]
                1.0e-3,
                2.0e-3,
                3.0e-3,
                4.0e-3,
                5.0e-3,
                6.0e-3,
                7.0e-3,
                8.0e-3,
                9.0e-3,
                10.0e-3
            ],
            [
                0.00,
                90.96 /  1.0e-3,               # number of photons per GeV
                277.36 /  2.0e-3,              # from Geant4 simulation by Daniel Lopez
                485.82 /  3.0e-3,
                692.83 /  4.0e-3,
                890.01 /  5.0e-3,
                1098.53 /  6.0e-3,
                1285.47 /  7.0e-3,
                1502.86 /  8.0e-3,
                1687.15 /  9.0e-3,
                1891.00 / 10.0e-3
            ] / value,
            WeakBoundaries(),
        )


        f2 = LinearInterpolator(
            [-2.0, -1.0, 0.0, 1.0, 2.0],
            map(log, [1.891e3 / 1.0e-2, 1.905e4 / 1.0e-1, 1.889e5 / 1.0e+0, 1.875e6 / 1.0e+1, 1.881e7 / 1.0e+2] / value),
            WeakBoundaries(),
        )

        return new(lambda_min, lambda_max, nx, dispersion_model, f1, f2)
    end
end

"""
Get correction of number of photons from EM-shower as a function of energy [GeV].
"""
function (emsc::EMShowerCorrection)(energy)
    energy <= emsc.f1.r.xa && return 0.0
    energy <= emsc.f1.r.xb && return emsc.f1(energy)
    x = log10(energy)
    x <= emsc.f2.r.xa && return exp(first(emsc.f2.r.y))
    x <= emsc.f2.r.xb && return exp(emsc.f2(x))
    exp(last(emsc.f2.r.y))
end
