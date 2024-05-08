module LumenManufaktur

using BasicInterpolators
using FastGaussQuadrature

const α = 1.0 / 137.036  # Fine-structure constant

include("dispersion.jl")
include("scattering.jl")
include("absorption.jl")
include("pmt.jl")
include("parameters.jl")
include("exports.jl")


"""
    cherenkov(λ, n)

Compute the number of Cherenkov photons induced by a muon per unit track length
and per unit wavelength [m^-1 nm^-1].

# Arguments
- `λ`: the wavelength of light [nm]
- `n`: index of refraction

"""
function cherenkov(λ, n)
    x = n * λ
    1.0e9 * 2 * π * α * (n^2 - 1.0) / (x^2)
end


"""
    directlightfrommuon(params::LMParameters, pmt::PMTModel, R, θ, ϕ)

Returns the number of photo-electrons from direct Cherenkov light from a muon
with a closest distance of `R` [m] to the PMT and the angles `θ` [rad] (zenith)
and `ϕ` [rad] (azimuth) with respect to the PMT axis.

# Arguments

- `R`: the (closest) distance [m] between muon and PMT
- `ϕ`: the zenith angle [rad] which is 0 when the PMT points away from the muon
  track (in x-direction) and rotates counter clockwise to the y-axis when
  viewed from above, while the z-axis points upwards.
- `θ`: the azimuth angle [rad] which is 0 when the PMT points upwards (along the
  z-axis) and π/2 when pointing downwards. Between 0 and π/2, the PMT points to
  the z-axis


"""
function directlightfrommuon(params::LMParameters, pmt::PMTModel, R, θ, ϕ)
    value = 0

    R = max(R, params.minimum_distance)
    A = pmt.photocathode_area
    px = sin(θ) * cos(ϕ)
    pz = cos(θ)

    for (m_x, m_y) in zip(params.legendre_coefficients...)
        w =
            0.5 * (params.lambda_max + params.lambda_min) +
            m_x * 0.5 * (params.lambda_max - params.lambda_min)
        dw = m_y * 0.5 * (params.lambda_max - params.lambda_min)
        n = refractionindexphase(params.dispersion_model, w)
        l_abs = absorptionlength(params.absorption_model, w)  # 59 allocs, 500ns
        ls = scatteringlength(params.scattering_model, w)  # 0 allocs, 180ns
        npe = cherenkov(w, n) * dw * pmt.quantum_efficiency(w)
        ct0 = 1.0 / n
        st0 = sqrt((1.0 + ct0) * (1.0 - ct0))
        d = R / st0                                 # distance traveled by photon
        ct = st0 * px + ct0 * pz                         # cosine angle of incidence on PMT
        U = pmt.angular_acceptance(ct)                # PMT angular acceptance
        V = exp(-d / l_abs) * exp(-d / ls)              # absorption & scattering
        W = A / (2.0 * π * R * st0)                      # solid angle
        value += npe * U * V * W
    end
    value
end


"""
    directlightfrommuon(params::LMParameters, pmt::PMTModel, R, θ, ϕ, Δt)

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
