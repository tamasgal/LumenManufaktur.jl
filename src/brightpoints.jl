"""
    directlightfrombrightpoint(params::LMParameters, pmt::PMTModel, D, ct, Δt)

Probability density function for direct light from an isotropic light source. Returns
``d^2P/dt/dE`` [npe/ns/GeV].

# Arguments
- `params`: parameters of the setup
- `pmt`: PMT model
- `D`: distance between light source and PMT [m]
- `ct`: cosine angle between light source direction and PMT direction
- `Δt`: time difference relative to direct light
"""

function directlightfrombrightpoint(params::LMParameters, pmt::PMTModel, D, ct, Δt)
    D = max(D, params.minimum_distance)
    t = D * params.n / C + Δt  # time [ns]
    A = pmt.photocathode_area

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max)
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min)
    ng = C * t / D  # index of refraction

    if (n0 >= ng)
        return 0.0
    end
    if (n1 <= ng)
        return 0.0
    end

    w = wavelength(params.dispersion_model, ng, params.lambda_max, 1.0e-5)
    n = refractionindexphase(params.dispersion_model, w)

    l_abs = absorptionlength(params.absorption_model, w)
    ls = scatteringlength(params.scattering_model, w)

    npe = cherenkov(w, n) * pmt.quantum_efficiency_(w)

    U = pmt.angular_acceptance(ct)  # PMT angular acceptance
    V = exp(-D / l_abs - D / ls)  # absorption & scattering
    W = A / (D * D)  # solid angle

    ngp = dispersiongroup(params.dispersion_model, w)

    Ja = D * ngp / C  # dt/dlambda
    Jb = 1.0/(4π)  # d^2N/dΩ

    return npe * geanc() * U * V * W * Jb / abs(Ja)
end

"""
    scatteredlightfrombrightpoint(params::LMParameters, pmt::PMTModel, D, ct, Δt)

Probability density function for scattered light from an isotropic light source. Returns
``d^2P/dt/dE`` [npe/ns/GeV].

# Arguments
- `params`: parameters of the setup
- `pmt`: PMT model
- `D`: distance between light source and PMT [m]
- `ct`: cosine angle between light source direction and PMT direction
- `Δt`: time difference relative to direct light
"""
function scatteredlightfrombrightpoint(params::LMParameters, pmt::PMTModel, D, ct, Δt)
    value = 0.0

    st = sqrt((1.0 + ct) * (1.0 - ct))
    D = max(D, params.minimum_distance)
    t = D * params.n / C + Δt  # time [ns]

    A = pmt.photocathode_area

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max)
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min)

    ni = C * t / D  # maximal index of refraction

    if (n0 >= ni)
        return value
    end

    nj = min(ni, n1)

    w = wmax

    n_coefficients = length(params.legendre_coefficients[1])

    for (m_x, m_y) in zip(params.legendre_coefficients...)

        ng = 0.5 * (nj + n0) + m_x * 0.5 * (nj - n0)
        dn = m_y * 0.5 * (nj - n0)

        w = getWavelength(ng, w, 1.0e-5)

        dw = dn / abs(dispersiongroup(params.dispersion_model, w))

        n = refractionindexphase(params.dispersion_model, w)

        npe = cherenkov(w, n) * dw * pmt.quantum_efficiency(w)
        
        if (npe <= 0)
            continue
        end
        
        l_abs = absorptionlength(params.absorption_model, w)
        ls = scatteringlength(params.scattering_model, w)

        Jc = 1.0 / ls  # dN/dx
        Jb = 1.0/(4π)  # d^2N/dΩ

        d = C * t / ng  # photon path

        ds = 2.0 / (n_coefficients + 1)

        sb = 0.5ds
        while sb < 1.0 - 0.25ds

            for k in (-1, 1)

                cb = k * sqrt((1.0 + sb) * (1.0 - sb))
                dcb = k * ds * sb / cb

                v = 0.5 * (d + D) * (d - D) / (d - D * cb)
                u = d - v

                if (u <= 0)
                    continue
                end
                if (v <= 0)
                    continue
                end

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

                W = min(A / (v * v), 2.0 * PI)  # solid angle
                Ja = scatteringprobability(params.scattering_probability_model, cts)  # d^2P/dcos/dϕ
                Jd = ng * (1.0 - cts) / C  # dt/du

                dp = π / length(params.integration_points)
                dom = dcb * dp * v * v / (u * u)

                for (cp, sp) in params.integration_points.xy

                    dot = cb * ct - sb * cp * st

                    U = 2 * pmt.angular_acceptance(dot)  # PMT angular acceptance

                    value += npe * geanc() * dom * U * V * W * Ja * Jb * Jc / abs(Jd)
                end
            end
            sb += ds
        end
    end

    return value
end
