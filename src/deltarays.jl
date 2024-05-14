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

            ct0 = -z / d
            st0 = R / d

            dz = D / st0  # average track length
            ct = st0 * px + ct0 * pz  # cosine angle of incidence on PMT

            U = pmt.angular_acceptance(ct)  # PMT angular acceptance
            V = exp(-d / l_abs - d / ls)  # absorption & scattering
            W = A / (d * d)  # solid angle

            Ja = deltarayprobability(ct0)  # d^2N/dcos/dϕ
            Jd = (1.0 - ng * ct0) / C  # d  t/ dz
            Je = ng * st0 * st0 * st0 / (R * C)  # d^2t/(dz)^2

            value += npe * geanc() * U * V * W * Ja / (abs(Jd) + 0.5 * Je * dz)
        end

    end
    return value
end

"""
    scatteredlightfromdeltarays(params::LMParameters, pmt::PMTModel, R, θ, ϕ, Δt)

Probability density function for scattered light from delta-rays with a
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
function scatteredlightfromdeltarays(params::LMParameters, pmt::PMTModel, R, θ, ϕ, Δt)
    value = 0.0

    R = max(R, params.minimum_distance)
    A = pmt.photocathode_area
    t = R * tanthetac(params.n) / C + Δt  # time [ns]

    px = sin(θ) * cos(ϕ)
    py = sin(θ) * sin(ϕ)
    pz = cos(θ)

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max)
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min)
    ni = sqrt(R^2 + C^2 * t^2) / R  # maximal index of refraction

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

        root = Root(R, ng, t)  # square root

        !root.isvalid && continue

        zmin = root.most_upstream
        zmax = root.most_downstream

        zap = 1.0 / l_abs

        xmin = exp(zap * zmax)
        xmax = exp(zap * zmin)

        n_coefficients = length(params.legendre_coefficients[1])

        for (p_x, p_y) in zip(params.legendre_coefficients...)

            x = 0.5 * (xmax + xmin) + p_x * 0.5 * (xmax - xmin)
            dx = p_y * 0.5 * (xmax - xmin)

            z = log(x) / zap
            dz = -dx / (zap * x)

            D = sqrt(z^2 + R^2)
            cd = -z / D
            sd = R / D

            qx = cd * px + 0 - sd * pz
            qy = 1 * py
            qz = sd * px + 0 + cd * pz

            d = (C * t - z) / ng  # photon path

            ds = 2.0 / (n_coefficients + 1)

            sb = 0.5ds
            while sb < 1.0 - 0.25ds

                for k in (-1, 1)

                    cb = k * sqrt((1.0 + sb) * (1.0 - sb))
                    dcb = k * ds * sb / cb

                    v = 0.5 * (d + D) * (d - D) / (d - D * cb)
                    u = d - v

                    u <= 0 && continue
                    v <= 0 && continue

                    cts = (D * cb - v) / u  # cosine scattering angle

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

                    W = min(A / (v * v), 2.0 * π)  # solid angle
                    Ja = scatteringprobability(params.scattering_probability_model, cts)  # d^2P/dcos/dϕ
                    Jd = ng * (1.0 - cts) / C  # dt/du

                    ca = (D - v * cb) / u
                    sa = v * sb / u

                    dp = π / length(params.integration_points)
                    dom = dcb * dp * v^2 / (u^2)

                    for (cp, sp) in params.integration_points.xy
                        ct0 = cd * ca - sd * sa * cp

                        vx = sb * cp * qx
                        vy = sb * sp * qy
                        vz = cb * qz

                        U =
                            pmt.angular_acceptance(vx + vy + vz) +
                            pmt.angular_acceptance(vx - vy + vz)  # PMT angular acceptance

                        Jb = deltarayprobability(ct0)  # d^2N/dcos/dϕ

                        value +=
                            dom * npe * geanc() * dz * U * V * W * Ja * Jb * Jc / abs(Jd)
                    end
                end
                sb += ds

            end
        end
    end

    value
end
