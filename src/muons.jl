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
    deltayprobability(x)

Emission profile of photons from delta-rays.

Profile is taken from reference ANTARES-SOFT-2002-015, J. Brunner (fig. 3).

Returns the probability.

# Arguments
- `x`: cosine emission angle
"""
function deltarayprobability(x)
    0.188 * exp(-1.25 * abs(x - 0.90)^1.30)
end



"""
    geanc()

Equivalent muon track length per unit shower energy [m/GeV].

See ANTARES internal note ANTARES-SOFT-2002-015, J. Brunner.
"""
@inline geanc() = 4.7319  # dx/dE [m/GeV]


"""
    directlightfrommuon(params::LMParameters, pmt::PMTModel, R, θ, ϕ)

Returns the number of photo-electrons from direct Cherenkov light from a muon
with a closest distance of `R` [m] to the PMT and the angles `θ` \\[rad\\] (zenith)
and `ϕ` \\[rad\\] (azimuth) with respect to the PMT axis.

# Arguments

- `params`: parameters of the setup
- `pmt`: the PMT model
- `R`: (closest) distance [m] between muon and PMT
- `ϕ`: zenith angle \\[rad\\] which is 0 when the PMT points away from the muon
  track (in x-direction) and rotates counter clockwise to the y-axis when
  viewed from above, while the z-axis points upwards.
- `θ`: azimuth angle \\[rad\\] which is 0 when the PMT points upwards (along the
  z-axis) and π/2 when pointing downwards. Between 0 and π/2, the PMT points to
  the z-axis
"""
function directlightfrommuon(params::LMParameters, pmt::PMTModel, R, θ, ϕ)
    value = 0.0

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
closest distance of `R` [m] to the PMT, the angles `θ` \\[rad\\] (zenith) and `ϕ`
\\[rad\\] (azimuth) with respect to the PMT axis and a time difference `Δt` [ns]
relative to direct Cherenkov light. Returns dP/dt [npe/ns].

# Arguments

- `params`: parameters of the setup
- `pmt`: the PMT model
- `R`: (closest) distance [m] between muon and PMT
- `ϕ`: zenith angle \\[rad\\] which is 0 when the PMT points away from the muon
  track (in x-direction) and rotates counter clockwise to the y-axis when
  viewed from above, while the z-axis points upwards.
- `θ`: azimuth angle \\[rad\\] which is 0 when the PMT points upwards (along the
  z-axis) and π/2 when pointing downwards. Between 0 and π/2, the PMT points to
  the z-axis
- `Δt`: time difference [ns] relative to the Cherenkov light
"""
function directlightfrommuon(params::LMParameters, pmt::PMTModel, R, θ, ϕ, Δt)
    N = 100      # maximal number of iterations
    eps = 1.0e-6 # precision index of refraction

    R = max(R, params.minimum_distance)
    A = pmt.photocathode_area
    t = R * tanthetac(params.n) / C + Δt # time [ns]
    a = C * t / R                      # target value
    px = sin(θ) * cos(ϕ)
    pz = cos(θ)

    # check validity range for index of refraction
    for λ in (params.lambda_min, params.lambda_max)

        n = refractionindexphase(params.dispersion_model, λ)
        ng = refractionindexgroup(params.dispersion_model, λ)

        ct0 = 1.0 / n
        st0 = sqrt((1.0 + ct0) * (1.0 - ct0))

        b = (ng - ct0) / st0 # running value

        λ == params.lambda_min && b < a && return 0.0
        λ == params.lambda_max && b > a && return 0.0

    end

    umin = params.lambda_min
    umax = params.lambda_max

    for _ = 1:N  # binary search for wavelength

        w = 0.5 * (umin + umax)

        n = refractionindexphase(params.dispersion_model, w)
        ng = refractionindexgroup(params.dispersion_model, w)

        ct0 = 1.0 / n
        st0 = sqrt((1.0 + ct0) * (1.0 - ct0))

        b = (ng - ct0) / st0  # running value

        if abs(b - a) < a * eps

            np = dispersionphase(params.dispersion_model, w)
            ngp = dispersiongroup(params.dispersion_model, w)
            l_abs = absorptionlength(params.absorption_model, w)
            ls = scatteringlength(params.scattering_model, w)

            d = R / st0  # distance traveled by photon
            ct = st0 * px + ct0 * pz  # cosine angle of incidence on PMT

            i3 = ct0 * ct0 * ct0 / (st0 * st0 * st0)

            U = pmt.angular_acceptance(ct)
            V = exp(-d / l_abs - d / ls)  # absorption & scattering
            W = A / (2.0 * π * R * st0)  # solid angle

            Ja = R * (ngp / st0 + np * (n - ng) * i3) / C  # dt/dlambd

            return cherenkov(w, n) * pmt.quantum_efficiency(w) * U * V * W / abs(Ja)
        end

        if b > a
            umin = w
        else
            umax = w
        end
    end

    0.0
end

"""
    scatteredlightfrommuon(params::LMParameters, pmt::PMTModel, D, cd, θ, ϕ, Δt)

Probability density function for scattered light from muon. Returns [d^2P/dt/dx].

# Arguments
- `params`: parameters of the setup
- `pmt`: PMT model
- `D`: distance between track segment and PMT [m]
- `cd`: cosine angle of muon direction and track segment - PMT position
- `θ`: zenith  angle orientation PMT \\[rad\\]
- `ϕ`: azimuth angle orientation PMT \\[rad\\]
- `Δt`: time difference relative to direct Cherenkov light
"""
function scatteredlightfrommuon(params::LMParameters, pmt::PMTModel, D, cd, θ, ϕ, Δt)
    eps = 1.0e-10

    value = 0.0

    sd = sqrt((1.0 + cd) * (1.0 - cd))
    D = max(D, params.minimum_distance)
    R = sd * D  # minimal distance of approach [m]
    Z = -cd * D  # photon emission point
    L = D
    t = D * params.n / C + Δt  # time [ns]
    A = pmt.photocathode_area

    px = sin(θ) * cos(ϕ)
    py = sin(θ) * sin(ϕ)
    pz = cos(θ)


    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max)
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min)

    ni = C * t / L  # maximal index of refraction

    n0 >= ni && return value

    nj = min(ni, n1)

    w = params.lambda_max

    for (m_x, m_y) in zip(params.legendre_coefficients...)

        ng = 0.5 * (nj + n0) + m_x * 0.5 * (nj - n0)
        dn = m_y * 0.5 * (nj - n0)

        w = wavelength(params.dispersion_model, ng, w, 1.0e-5)

        dw = dn / abs(dispersiongroup(params.dispersion_model, w))

        n = refractionindexphase(params.dispersion_model, w)

        l_abs = absorptionlength(params.absorption_model, w)
        ls = scatteringlength(params.scattering_model, w)

        npe = cherenkov(w, n) * dw * pmt.quantum_efficiency(w)

        npe <= 0 && continue

        Jc = 1.0 / ls  # dN/dx

        ct0 = 1.0 / n  # photon direction before scattering
        st0 = sqrt((1.0 + ct0) * (1.0 - ct0))

        d = C * t / ng  # photon path

        cta = cd * ct0 + sd * st0
        dca = d - 0.5 * (d + L) * (d - L) / (d - L * cta)
        tip = -log(L * L / (dca * dca) + eps) / π

        ymin = exp(tip * π)
        ymax = 1.0

        for (q_x, q_y) in zip(params.legendre_coefficients...)

            y = 0.5 * (ymax + ymin) + q_x * 0.5 * (ymax - ymin)
            dy = q_y * 0.5 * (ymax - ymin)

            ϕ = log(y) / tip
            dp = -dy / (tip * y)

            cp0 = cos(ϕ)
            sp0 = sin(ϕ)

            u = (R * R + Z * Z - d * d) / (2 * R * st0 * cp0 - 2 * Z * ct0 - 2 * d)
            v = d - u

            u <= 0 && continue
            v <= 0 && continue

            vi = 1.0 / v
            cts = (R * st0 * cp0 - Z * ct0 - u) * vi  # cosine scattering angle

            V = exp(
                -d * inverseattenuationlength(
                    params.scattering_probability_model,
                    l_abs,
                    ls,
                    cts,
                ),
            )

            cts < 0.0 &&
                v * sqrt((1.0 + cts) * (1.0 - cts)) < params.module_radius &&
                continue

            vx = R - u * st0 * cp0
            vy = -u * st0 * sp0
            vz = -Z - u * ct0

            ct = (  # cosine angle of incidence on PMT
                (vx * px + vy * py + vz * pz) * vi,
                (vx * px - vy * py + vz * pz) * vi,
            )

            U = pmt.angular_acceptance(ct[1]) + pmt.angular_acceptance(ct[2])  # PMT angular acceptance
            W = min(A * vi * vi, 2π)  # solid angle

            Ja = scatteringprobability(params.scattering_probability_model, cts)  # d^2P/dcos/dϕ
            Jd = ng * (1.0 - cts) / C           # dt/du

            value += (npe * dp / (2 * π)) * U * V * W * Ja * Jc / abs(Jd)
        end
    end

    return value
end



"""
    directlightfromdeltarays(params::LMParameters, pmt::PMTModel, R, θ, ϕ, Δt)

Probability density function for direct light from delta-rays with a
closest distance of `R` [m] to the PMT, the angles `θ` \\[rad\\] (zenith) and `ϕ`
\\[rad\\] (azimuth) with respect to the PMT axis and a time difference `Δt` [ns]
relative to direct Cherenkov light. Returns d^2P/dt/dE [npe/ns ⋅ m/GeV].

# Arguments

- `params`: parameters of the setup
- `pmt`: the PMT model
- `R`: (closest) distance [m] between muon and PMT
- `ϕ`: zenith angle \\[rad\\] which is 0 when the PMT points away from the muon
  track (in x-direction) and rotates counter clockwise to the y-axis when
  viewed from above, while the z-axis points upwards.
- `θ`: azimuth angle \\[rad\\] which is 0 when the PMT points upwards (along the
  z-axis) and π/2 when pointing downwards. Between 0 and π/2, the PMT points to
  the z-axis
- `Δt`: time difference [ns] relative to the Cherenkov light
"""
function directlightfromdeltarays(params::LMParameters, pmt::PMTModel, R, θ, ϕ, Δt)
    value = 0.0

    R = max(R, params.minimum_distance)
    A = pmt.photocathode_area
    t = R * tanthetac(params.n) / C + Δt  # time [ns]
    D = 2.0 * sqrt(A / π)

    px = sin(θ) * cos(ϕ)
    pz = cos(θ)

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max)
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min)
    ni = sqrt(R * R + C * t * C * t) / R  # maximal index of refraction

    n0 >= ni && return value

    nj = min(n1, ni)

    w = params.lambda_max

    for (m_x, m_y) in zip(params.legendre_coefficients...)

        ng = 0.5 * (nj + n0) + m_x * 0.5 * (nj - n0)
        dn = m_y * 0.5 * (nj - n0)

        w = wavelength(params.dispersion_model, ng, w, 1.0e-5)

        dw = dn / abs(dispersiongroup(params.dispersion_model, w))

        n = refractionindexphase(params.dispersion_model, w)

        l_abs = absorptionlength(params.absorption_model, w)
        ls = scatteringlength(params.scattering_model, w)

        npe = cherenkov(w, n) * dw * pmt.quantum_efficiency(w)

        root = Root(R, ng, t)  # square root

        !root.isvalid && continue

        for z in (root.most_upstream, root.most_downstream)

            C * t <= z && continue

            d = sqrt(z * z + R * R)  # distance traveled by photon

            ct0 = -z / d;
            st0 = R / d;

            dz = D / st0  # average track length
            ct = st0 * px + ct0 * pz  # cosine angle of incidence on PMT

            U = pmt.angular_acceptance(ct)  # PMT angular acceptance
            V = exp(-d / l_abs - d / ls)  # absorption & scattering
            W = A / (d * d)  # solid angle

            Ja = deltarayprobability(ct0)  # d^2N/dcos/dϕ
            Jd = (1.0 - ng * ct0) / C  # d  t/ dz
            Je = ng * st0 * st0 * st0 / (R * C)  # d^2t/(dz)^2

            value += npe * geanc() * U * V * W * Ja / (abs(Jd) + 0.5 * Je * dz);
        end

    end
    return value
end

# /**
#  * Probability density function for scattered light from delta-rays.
#  *
#  * \param  R_m        distance between muon and PMT [m]
#  * \param  θ      zenith  angle orientation PMT [rad]
#  * \param  ϕ        azimuth angle orientation PMT [rad]
#  * \param  Δt       time difference relative to direct Cherenkov light [ns]
#  * \return            d^2P/dt/dx [npe/ns x m/GeV]
#  */
# double getScatteredLightFromDeltaRays(const double R_m, const double θ,
#                                       const double ϕ,
#                                       const double Δt) const {
#   double value = 0;

#   const double R = std::max(R_m, getRmin());
#   const double t = R * getTanΘC() / C + Δt; // time [ns]
#   const double A = pmt.photocathode_area;

#   const double px = sin(θ) * cos(ϕ);
#   const double py = sin(θ) * sin(ϕ);
#   const double pz = cos(θ);

#   const double n0 = refractionindexgroup(params.dispersion_model, params.lambda_max);
#   const double n1 = refractionindexgroup(params.dispersion_model, params.lambda_min);
#   const double ni =
#       sqrt(R * R + C * t * C * t) / R; // maximal index of refraction

#   if (n0 >= ni) {
#     return value;
#   }

#   const double nj = std::min(ni, n1);

#   double w = params.lambda_max;

#   for (const_iterator i = begin(); i != end(); ++i) {

#     const double ng = 0.5 * (nj + n0) + i->getX() * 0.5 * (nj - n0);
#     const double dn = i->getY() * 0.5 * (nj - n0);

#     w = getWavelength(ng, w, 1.0e-5);

#     const double dw = dn / fabs(getDispersionGroup(w));

#     const double n = getIndexOfRefractionPhase(w);

#     const double l_abs = getAbsorptionLength(w);
#     const double ls = getScatteringLength(w);

#     const double npe = cherenkov(w, n) * dw * getQE(w);

#     if (npe <= 0) {
#       continue;
#     }

#     const double Jc = 1.0 / ls; // dN/dx

#     JRoot rz(R, ng, t); // square root

#     if (!rz.is_valid) {
#       continue;
#     }

#     const double zmin = rz.first;
#     const double zmax = rz.second;

#     const double zap = 1.0 / l_abs;

#     const double xmin = exp(zap * zmax);
#     const double xmax = exp(zap * zmin);

#     for (const_iterator j = begin(); j != end(); ++j) {

#       const double x = 0.5 * (xmax + xmin) + j->getX() * 0.5 * (xmax - xmin);
#       const double dx = j->getY() * 0.5 * (xmax - xmin);

#       const double z = log(x) / zap;
#       const double dz = -dx / (zap * x);

#       const double D = sqrt(z * z + R * R);
#       const double cd = -z / D;
#       const double sd = R / D;

#       const double qx = cd * px + 0 - sd * pz;
#       const double qy = 1 * py;
#       const double qz = sd * px + 0 + cd * pz;

#       const double d = (C * t - z) / ng; // photon path

#       // const double V  = exp(-d/l_abs);                         //
#       // absorption

#       const double ds = 2.0 / (size() + 1);

#       for (double sb = 0.5 * ds; sb < 1.0 - 0.25 * ds; sb += ds) {

#         for (int buffer[] = {-1, +1, 0}, *k = buffer; *k != 0; ++k) {

#           const double cb = (*k) * sqrt((1.0 + sb) * (1.0 - sb));
#           const double dcb = (*k) * ds * sb / cb;

#           const double v = 0.5 * (d + D) * (d - D) / (d - D * cb);
#           const double u = d - v;

#           if (u <= 0) {
#             continue;
#           }
#           if (v <= 0) {
#             continue;
#           }

#           const double cts = (D * cb - v) / u; // cosine scattering angle

#           const double V =
#               exp(-d * getInverseAttenuationLength(l_abs, ls, cts));

#           if (cts < 0.0 &&
#               v * sqrt((1.0 + cts) * (1.0 - cts)) < MODULE_RADIUS_M) {
#             continue;
#           }

#           const double W = std::min(A / (v * v), 2.0 * PI); // solid angle
#           const double Ja = getScatteringProbability(cts);  // d^2P/dcos/dϕ
#           const double Jd = ng * (1.0 - cts) / C;           // dt/du

#           const double ca = (D - v * cb) / u;
#           const double sa = v * sb / u;

#           const double dp = PI / phd.size();
#           const double dom = dcb * dp * v * v / (u * u);

#           for (const_iterator l = phd.begin(); l != phd.end(); ++l) {

#             const double cp = l->getX();
#             const double sp = l->getY();

#             const double ct0 = cd * ca - sd * sa * cp;

#             const double vx = sb * cp * qx;
#             const double vy = sb * sp * qy;
#             const double vz = cb * qz;

#             const double U =
#                 pmt.angular_acceptance(vx + vy + vz) +
#                 pmt.angular_acceptance(vx - vy + vz); // PMT angular acceptance

#             const double Jb = deltarayprobability(ct0); // d^2N/dcos/dϕ

#             value += dom * npe * geanc() * dz * U * V * W * Ja * Jb * Jc /
#                      fabs(Jd);
#           }
#         }
#       }
#     }
#   }

#   return value;
# }

"""
    Root(R, n, t)

Helper struct to find solution(s) to ``z`` of the square root expression:

```math
\begin{aligned}
ct(z=0) & = & z + n \\sqrt{z^2 + R^2}
\end{aligned}
```

where ``n = 1/\\cos(\\theta_{c})`` is the index of refraction.

# Arguments

- `R`: minimal distance of approach [m]
- `n`: index of refraction
- `t`: time at z = 0 [ns]

# Fields

- `most_upstream`: most upstream solution
- `most_downstream`: most downstream solution
- `isvalid`: validity of the solution
```
"""
struct Root
    most_upstream::Float64
    most_downstream::Float64
    isvalid::Bool

    function Root(R, n, t)
        a = n^2 - 1.0
        b = 2 * C * t
        c = R^2 * n^2 - C^2 * t^2

        q = b * b - 4 * a * c

        isvalid = false
        if q >= 0.0

            first = (-b - sqrt(q)) / (2a)
            second = (-b + sqrt(q)) / (2a)

            isvalid = C * t > second
        end
        new(first, second, isvalid)
    end
end
