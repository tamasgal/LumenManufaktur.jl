module LumenManufaktur

using BasicInterpolators
using FastGaussQuadrature

include("exports.jl")

abstract type MediumProperties end
abstract type DispersionModel end
abstract type ScatteringModel end
abstract type AbsorptionModel end

struct Kopelevich <: ScatteringModel end
struct DefaultAbsorption <: AbsorptionModel end


struct PMTModel{T1, T2}
    photocathode_area::Float64
    quantum_efficiency::T1
    angular_acceptance::T2
end

# abstract type PMTModel end
# function quantumefficiency end
# photocathodearea(::Type{KM3NeTPMT}) = 45.4e-4
# quantumefficiency(::Type{KM3NeTPMT}) = LinearInterpolator(
#         [0, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335,
#         340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410,
#         415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485,
#         490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560,
#         565, 570, 575, 580, 585, 590, 595, 600, 605, 610, 615, 620, 625, 630, 635,
#         640, 645, 650, 655, 660, 665, 670, 675, 680, 685, 690, 695, 700, 705, 710,
#         999999],
#         # collection efficiency (with factor 0.9) correction
#         0.01 * 0.9 * [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.9, 1.8, 4.6, 7.4,
#         11., 15., 18., 21., 22., 23., 24., 25., 25., 26., 26., 26., 26., 27., 26.,
#         26., 26., 26., 26., 25., 25., 25., 25., 24., 24., 24., 23., 22., 22., 21.,
#         20., 19., 19., 19., 18., 18., 18., 17., 16., 15., 14., 12., 11., 10., 9.5,
#         8.7, 8.1, 7.6, 7.2, 6.8, 6.4, 6.0, 5.6, 5.1, 4.8, 4.4, 4.0, 3.6, 3.3, 3.0,
#         2.7, 2.3, 2.1, 1.9, 1.7, 1.5, 1.4, 1.2, 1.0, 0.9, 0.7, 0.6, 0.5, 0.4, 0.3,
#         0.2, 0.1, 0.0, 0.0],
#         NoBoundaries()
# )
#
const α = 1.0/137.036  # Fine-structure constant

const PMTKM3NeT = PMTModel(
    45.4e-4,
    # Yields slightly different values compared to the Jpp getQE() function,
    # very probably due to the interpolation implementation
    LinearInterpolator(
        [0, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335,
        340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410,
        415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485,
        490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560,
        565, 570, 575, 580, 585, 590, 595, 600, 605, 610, 615, 620, 625, 630, 635,
        640, 645, 650, 655, 660, 665, 670, 675, 680, 685, 690, 695, 700, 705, 710,
        999999],

        # collection efficiency (with factor 0.9) correction
        0.01 * 0.9 * [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.9, 1.8, 4.6, 7.4,
        11., 15., 18., 21., 22., 23., 24., 25., 25., 26., 26., 26., 26., 27., 26.,
        26., 26., 26., 26., 25., 25., 25., 25., 24., 24., 24., 23., 22., 22., 21.,
        20., 19., 19., 19., 18., 18., 18., 17., 16., 15., 14., 12., 11., 10., 9.5,
        8.7, 8.1, 7.6, 7.2, 6.8, 6.4, 6.0, 5.6, 5.1, 4.8, 4.4, 4.0, 3.6, 3.3, 3.0,
        2.7, 2.3, 2.1, 1.9, 1.7, 1.5, 1.4, 1.2, 1.0, 0.9, 0.7, 0.6, 0.5, 0.4, 0.3,
        0.2, 0.1, 0.0, 0.0],

        NoBoundaries()
    ),

    # using eps(), max angle and two zero values in the second vector to
    # get real zeros outside of the upper boundary defined in Jpp (i.e. >=0.4)
    LinearInterpolator(
        [-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55,
        -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00,
        0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.40+eps(), 1.0],

        [1.621, 1.346, 1.193, 1.073, 0.973, 0.877, 0.790, 0.711, 0.640, 0.575,
        0.517, 0.450, 0.396, 0.341, 0.295, 0.249, 0.207, 0.166, 0.128, 0.095, 0.065,
        0.038, 0.017, 0.006, 0.003, 0.002, 0.001, 0.001, 0.001, 0.000, 0.000],

        NoBoundaries()
    )
)


Base.@kwdef struct BasicDispersion <: DispersionModel
    P::Float64 = 1                 #  ambient pressure [atm]
    a0::Float64 = 1.3201           #  offset
    a1::Float64 = 1.4e-5           #  dn/dP
    a2::Float64 = 16.2566          #  d^1n/(dx)^1
    a3::Float64 = -4383.0          #  d^2n/(dx)^2
    a4::Float64 = 1.1455e6         #  d^3n/(dx)^3
end
BasicDispersion(P) = BasicDispersion(P=P)
const DispersionORCA = BasicDispersion(240)
const DispersionARCA = BasicDispersion(350)


Base.@kwdef struct Parameters{D<:DispersionModel, S<:ScatteringModel, A<:AbsorptionModel}
    minimum_distance::Float64 = 1.0e-1
    lambda_min::Float64 = 300.0
    lambda_max::Float64 = 700.0
    legendre_coefficients::Tuple{Vector{Float64}, Vector{Float64}} = gausslegendre(5)
    dispersion_model::D = BasicDispersion()
    scattering_model::S = Kopelevich()
    absorption_model::A = DefaultAbsorption()
end
function Base.show(io::IO, p::Parameters)
    println(io, "Parameters:")
    println(io, "  minimum distance = $(p.minimum_distance)")
    println(io, "  lambda min/max = $(p.lambda_min)/$(p.lambda_max)")
    println(io, "  degree of Legendre polynomials = $(length(first(p.legendre_coefficients)))")
    println(io, "  dispersion model = $(p.dispersion_model)")
    println(io, "  scattering model = $(p.scattering_model)")
    print(io, "  absorption model = $(p.absorption_model)")
end

@inline function refractionindexphase(λ, dp::DispersionModel=BasicDispersion())
    x = 1.0 / λ
    dp.a0  +  dp.a1*dp.P  +  x*(dp.a2 + x*(dp.a3 + x*dp.a4))
end
@inline function refractionindexgroup(λ, dp::DispersionModel=BasicDispersion())
    error("Not implemented yet")
end

const absorptionlengthinterpolator = LinearInterpolator(
    [0, 290, 310, 330, 350, 375, 412, 440, 475, 488, 510, 532, 555, 650, 676, 715, 720, 999999],
    [0.0, 0.0, 11.9, 16.4, 20.6, 29.5, 48.5, 67.5, 59.0, 55.1, 26.1, 19.9, 14.7, 2.8, 2.3, 1.0, 0.0, 0.0],
    NoBoundaries()
)
"""
    absorptionlength(λ)

Returns the absorption length [m] for a given wavelength [nm].
"""
absorptionlength(::DefaultAbsorption, λ) = absorptionlengthinterpolator(λ)

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


"""
    cherenkov(λ, n)

Compute the number of Cherenkov photons induced by a muon per unit track length
and per unit wavelength [m^-1 nm^-1].

# Arguments
- `λ`: the wavelength of light [nm]
- `n`: index of refraction

"""
function cherenkov(λ, n)
  x = n*λ
  1.0e9 * 2 * π * α * (n^2 - 1.0) / (x^2)
end


"""
    directlightfrommuon(params::Parameters, pmt::PMTModel, R, θ, ϕ)

Returns the number of photo-electrons from direct Cherenkov light from a muon
with a distance of `R` [m] to the PMT and the angles `θ` [rad] (zenith) and `ϕ`
[rad] (azimuth) with respect to the PMT axis.

# Arguments

- `ϕ`: the zenith angle [rad] which is 0 when the PMT points away from the muon
  track (in x-direction) and rotates counter clockwise to the y-axis when
  viewed from above, while the z-axis points upwards.
- `θ`: the azimuth angle [rad] which is 0 when the PMT points upwards (along the
  z-axis) and π/2 when pointing downwards. Between 0 and π/2, the PMT points to
  the z-axis


"""
function directlightfrommuon(params::Parameters, pmt::PMTModel, R, θ, ϕ)
    value = 0

    R  =  max(R, params.minimum_distance)
    A  =  pmt.photocathode_area

    px =  sin(θ)*cos(ϕ)
    pz =  cos(θ)

    # @show px
    # @show pz
    # @show wmin
    # @show wmax

    # for (m_x, m_y) in [(cos(x), sin(x)) for x in 0.5dx:dx:π]
    for (m_x, m_y) in zip(params.legendre_coefficients...)
    # gl = gausslegendre(5)
    # m_x = gl[1][1]
    # m_y = gl[2][1]

      # println("-----------")
      # @show m_x
      # @show m_y

      w  = 0.5 * (params.lambda_max + params.lambda_min)  +  m_x * 0.5 * (params.lambda_max - params.lambda_min)
      # @show w

      dw = m_y * 0.5 * (params.lambda_max - params.lambda_min)
      # @show dw

      n     = refractionindexphase(w, params.dispersion_model)
      # @show n

      l_abs = absorptionlength(params.absorption_model, w);  # 59 allocs, 500ns
      # @show l_abs
      ls    = scatteringlength(params.scattering_model, w);  # 0 allocs, 180ns
      # @show ls

      # @show cherenkov(w,n)
      # @show PMTKM3NeT.quantum_efficiency(w)
      npe   = cherenkov(w,n) * dw * pmt.quantum_efficiency(w);
      # @show npe

      ct0 = 1.0 / n;
      st0 = sqrt((1.0 + ct0)*(1.0 - ct0));

      d  = R / st0;                                 # distance traveled by photon
      ct = st0*px + ct0*pz;                         # cosine angle of incidence on PMT

      U  = pmt.angular_acceptance(ct);                # PMT angular acceptance
      V  = exp(-d/l_abs) * exp(-d/ls);              # absorption & scattering
      W  = A / (2.0*π*R*st0);                      # solid angle

      value += npe * U * V * W;
    end

    value
end


"""
    directlightfrommuon(params::Parameters, pmt::PMTModel, R, θ, ϕ, Δt)

Probability density function for direct Cherenkov light from a muon with a
distance of `R` [m] to the PMT, time difference `Δt` [ns] relative to direct
Cherenkov light and the angles `θ` [rad] (zenith) and `ϕ` [rad] (azimuth) with
respect to the PMT axis.

Returns dP/dt [npe/ns]

"""
# function directlight(R, θ, ϕ, Δt)
#     int N = 100;         # maximal number of iterations
#     double eps = 1.0e-6; # precision index of refraction

#     R = max(R_m, getRmin());
#     t = R * getTanThetaC() / C + t_ns; # time [ns]
#     a = C * t / R;                     # target value
#     A = getPhotocathodeArea();

#     px = sin(theta) * cos(phi);
#     pz = cos(theta);

#     # check validity range for index of refraction

#     for (double buffer[] = {wmin, wmax, 0.0}, *p = buffer; *p != 0.0; ++p) {

#       n = getIndexOfRefractionPhase(*p);
#       ng = getIndexOfRefractionGroup(*p);

#       ct0 = 1.0 / n;
#       st0 = sqrt((1.0 + ct0) * (1.0 - ct0));

#       b = (ng - ct0) / st0; // running value

#       if (*p == wmin && b < a) {
#         return 0.0;
#       }
#       if (*p == wmax && b > a) {
#         return 0.0;
#       }
#     }

#     umin = wmin;
#     umax = wmax;

#     for (int i = 0; i != N; ++i) { // binary search for wavelength

#       w = 0.5 * (umin + umax);

#       n = getIndexOfRefractionPhase(w);
#       ng = getIndexOfRefractionGroup(w);

#       ct0 = 1.0 / n;
#       st0 = sqrt((1.0 + ct0) * (1.0 - ct0));

#       b = (ng - ct0) / st0; // running value

#       if (fabs(b - a) < a * eps) {

#         np = getDispersionPhase(w);
#         ngp = getDispersionGroup(w);

#         l_abs = getAbsorptionLength(w);
#         ls = getScatteringLength(w);

#         d = R / st0; // distance traveled by photon
#         ct = st0 * px + ct0 * pz; // cosine angle of incidence on PMT

#         i3 = ct0 * ct0 * ct0 / (st0 * st0 * st0);

#         U = angularacceptance(ct); // PMT angular acceptance
#         V = exp(-d / l_abs - d / ls); // absorption & scattering
#         W = A / (2.0 * PI * R * st0); // solid angle

#         Ja = R * (ngp / st0 + np * (n - ng) * i3) / C; // dt/dlambda

#         return cherenkov(w, n) * quantumefficiency(w) * U * V * W / fabs(Ja);
#       }

#       if (b > a)
#         umin = w;
#       else
#         umax = w;
#     }

#     return 0.0;
#   }

end