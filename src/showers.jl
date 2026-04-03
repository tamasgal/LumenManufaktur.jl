"""
    directlightfromEMshower(params::LMParameters, pmt::PMTModel, D, cd, θ, ϕ, Δt)

Probability density function for direct light from EM-shower. Returns
``d^2P/dt/dE`` [npe/ns/GeV].

# Arguments
- `params`: parameters of the setup
- `pmt`: PMT model
- `D`: distance between EM-shower and PMT [m]
- `cd`: cosine angle of EM-shower direction and EM-shower - PMT
- `θ`: zenith  angle orientation PMT [rad]
- `ϕ`: azimuth angle orientation PMT [rad]
- `Δt`: time difference relative to direct Cherenkov light
"""
function directlightfromEMshower(params::LMParameters, pmt::PMTModel, D, cd, θ, ϕ, Δt)
    ct0 = cd
    st0 = sqrt((1.0 + ct0) * (1.0 - ct0))

    D = max(D, params.minimum_distance)
    t = D * params.n / C + Δt  # time [ns]
    A = pmt.photocathode_area

    px = sin(θ) * cos(ϕ)
    pz = cos(θ)

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max)
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min)
    ng = C * t / D  # index of refraction

    if (n0 >= ng)
        return 0.0
    end
    if (n1 <= ng)
        return 0.0
    end

    w_start = params.lambda_max + (ng - n0) / (n1 - n0) * (params.lambda_min - params.lambda_max)
    w = wavelength(params.dispersion_model, ng, w_start, 1.0e-5)
    n = refractionindexphase(params.dispersion_model, w)

    l_abs = absorptionlength(params.absorption_model, w)
    ls = scatteringlength(params.scattering_model, w)

    npe = cherenkov(w, n) * pmt.quantum_efficiency(w)

    ct = st0 * px + ct0 * pz  # cosine angle of incidence on PMT

    U = pmt.angular_acceptance(ct)  # PMT angular acceptance
    V = exp(-D / l_abs - D / ls)  # absorption & scattering
    W = A / (D * D)  # solid angle

    ngp = dispersiongroup(params.dispersion_model, w)

    Ja = D * ngp / C  # dt/dlambda
    Jb = geant(n, ct0)  # d^2N/dcos/dϕ

    return npe * geanc() * U * V * W * Jb / abs(Ja)
end

"""
    scatteredlightfromEMshower(params::LMParameters, pmt::PMTModel, D, cd, θ, ϕ, Δt)

Probability density function for scattered light from EM-shower. Returns
``d^2P/dt/dE`` [npe/ns/GeV].

# Arguments
- `params`: parameters of the setup
- `pmt`: PMT model
- `D`: distance between EM-shower and PMT [m]
- `cd`: cosine angle of EM-shower direction and EM-shower - PMT
- `θ`: zenith  angle orientation PMT [rad]
- `ϕ`: azimuth angle orientation PMT [rad]
- `Δt`: time difference relative to direct Cherenkov light
"""
function scatteredlightfromEMshower(params::LMParameters, pmt::PMTModel, D, cd, θ, ϕ, Δt)
    value = 0.0

    sd = sqrt((1.0 + cd) * (1.0 - cd))
    D = max(D, params.minimum_distance)
    L = D
    t = D * params.n / C + Δt  # time [ns]

    A = pmt.photocathode_area

    px = sin(θ) * cos(ϕ)
    py = sin(θ) * sin(ϕ)
    pz = cos(θ)

    qx = cd * px + 0 - sd * pz
    qy = 1 * py
    qz = sd * px + 0 + cd * pz

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max)
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min)

    ni = C * t / L  # maximal index of refraction

    if (n0 >= ni)
        return value
    end

    nj = min(ni, n1)

    w = params.lambda_max

    n_coefficients = length(params.legendre_coefficients[1])

    for (m_x, m_y) in zip(params.legendre_coefficients...)

        ng = 0.5 * (nj + n0) + m_x * 0.5 * (nj - n0)
        dn = m_y * 0.5 * (nj - n0)

        w = wavelength(params.dispersion_model, ng, w, 1.0e-5)

        dw = dn / abs(dispersiongroup(params.dispersion_model, w))

        n = refractionindexphase(params.dispersion_model, w)

        l_abs = absorptionlength(params.absorption_model, w)
        ls = scatteringlength(params.scattering_model, w)

        npe = cherenkov(w, n) * dw * pmt.quantum_efficiency(w)

        if (npe <= 0)
            continue
        end

        Jc = 1.0 / ls  # dN/dx

        d = C * t / ng  # photon path

        ds = 2.0 / (n_coefficients + 1)

        sb = 0.5ds
        while sb < 1.0 - 0.25ds

            for k in (-1, 1)

                cb = k * sqrt((1.0 + sb) * (1.0 - sb))
                dcb = k * ds * sb / cb

                v = 0.5 * (d + L) * (d - L) / (d - L * cb)
                u = d - v

                if (u <= 0)
                    continue
                end
                if (v <= 0)
                    continue
                end

                cts = (L * cb - v) / u  # cosine scattering angle

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

                ca = (L - v * cb) / u
                sa = v * sb / u

                dp = π / length(params.integration_points)
                dom = dcb * dp * v * v / (u * u)
                dos = sqrt(dom)

                for (cp, sp) in params.integration_points.xy

                    ct0 = cd * ca - sd * sa * cp

                    vx = -sb * cp * qx
                    vy = -sb * sp * qy
                    vz = cb * qz

                    U =
                        pmt.angular_acceptance(vx + vy + vz) +
                        pmt.angular_acceptance(vx - vy + vz)  # PMT angular acceptance

                    Jb = geant(n, ct0 - 0.5 * dos, ct0 + 0.5 * dos)  # dN/dϕ

                    value += npe * geanc() * dos * U * V * W * Ja * Jb * Jc / abs(Jd)
                end
            end
            sb += ds
        end
    end

    return value
end

"""
    directlightfromEMshower(params::LMParameters, pmt::PMTModel, E, D, cd, θ, ϕ, Δt)

Probability density function for direct light from EM-shower. Returns
``dP/dt`` [npe/ns].

# Arguments
- `params`: parameters of the setup
- `pmt`: PMT model
- `E`: EM-shower energy [GeV]
- `D`: distance between EM-shower and PMT [m]
- `cd`: cosine angle of EM-shower direction and EM-shower - PMT
- `θ`: zenith  angle orientation PMT [rad]
- `ϕ`: azimuth angle orientation PMT [rad]
- `Δt`: time difference relative to direct Cherenkov light
"""
function directlightfromEMshower(params::LMParameters, pmt::PMTModel, E, D, cd, θ, ϕ, Δt)
    value = 0.0

    sd = sqrt((1.0 + cd) * (1.0 - cd))
    D = max(D, params.minimum_distance)
    R = D * sd  # minimal distance of approach [m]
    Z = -D * cd
    t = D * params.n / C + Δt  # time [ns]

    n_coefficients = length(params.legendre_coefficients[1])


    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max)
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min)

    zmin = 0.0  # minimal shower length [m]
    zmax = getlength(geanz, E, 1.0)  # maximal shower length [m]

    d = sqrt((Z + zmax) * (Z + zmax) + R * R)

    if (C * t > max(n1 * D, zmax + n1 * d))
        return value
    end

    if (C * t < n0 * D)
        root = Root(R, n0, t + Z / C)  # square root
        !root.isvalid && return value
        if root.most_downstream > Z
            zmin = root.most_downstream - Z
        end
        if root.most_upstream > Z
            zmin = root.most_upstream - Z
        end
    end

    if (C * t > n1 * D)
        root = Root(R, n1, t + Z / C)  # square root

        !root.isvalid && return value

        if (root.most_downstream > Z)
            zmin = root.most_downstream - Z
        end
        if (root.most_upstream > Z)
            zmin = root.most_upstream - Z
        end
    end

    if (C * t < zmax + n0 * d)

        root = Root(R, n0, t + Z / C)  # square root

        !root.isvalid && return value

        if (root.most_upstream > Z)
            zmax = root.most_upstream - Z
        end
        if (root.most_downstream > Z)
            zmax = root.most_downstream - Z
        end
    end

    if (zmin < 0.0)
        zmin = 0.0
    end

    if (zmax > zmin)

        ymin = getintegral(geanz, E, zmin)
        ymax = getintegral(geanz, E, zmax)
        dy = (ymax - ymin) / n_coefficients

        if dy > 2 * eps()

            for y = ymin+0.5dy:dy:ymax

                z = Z + getlength(geanz, E, y)
                d = sqrt(R * R + z * z)
                t1 = t + (Z - z) / C - d * params.n / C

                value += dy * E * directlightfromEMshower(params, pmt, d, -z / d, θ, ϕ, t1)
            end
        end
    end

    return value
end

"""
    scatteredlightfromEMshower(params::LMParameters, pmt::PMTModel, E, D, cd, θ, ϕ, Δt)

Probability density function for scattered light from EM-shower. Returns
``dP/dt`` [npe/ns].

# Arguments
- `params`: parameters of the setup
- `pmt`: PMT model
- `E`: EM-shower energy [GeV]
- `D`: distance between EM-shower and PMT [m]
- `cd`: cosine angle of EM-shower direction and EM-shower - PMT
- `θ`: zenith  angle orientation PMT [rad]
- `ϕ`: azimuth angle orientation PMT [rad]
- `Δt`: time difference relative to direct Cherenkov light
"""
function scatteredlightfromEMshower(params::LMParameters, pmt::PMTModel, E, D, cd, θ, ϕ, Δt)
    value = 0.0

    sd = sqrt((1.0 + cd) * (1.0 - cd))
    D = max(D, params.minimum_distance)
    R = D * sd  # minimal distance of approach [m]
    Z = -D * cd
    t = D * params.n / C + Δt  # time [ns]

    n_coefficients = length(params.legendre_coefficients[1])

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max)

    zmin = 0.0  # minimal shower length [m]
    zmax = getlength(geanz, E, 1.0)  # maximal shower length [m]

    d = sqrt((Z + zmax) * (Z + zmax) + R * R)

    if (C * t < n0 * D)

        root = Root(R, n0, t + Z / C)  # square root

        !root.isvalid && return value

        if (root.most_downstream > Z)
            zmin = root.most_downstream - Z
        end
        if (root.most_upstream > Z)
            zmin = root.most_upstream - Z
        end
    end

    if (C * t < zmax + n0 * d)

        root = Root(R, n0, t + Z / C)  # square root

        !root.isvalid && return value

        if (root.most_upstream > Z)
            zmax = root.most_upstream - Z
        end
        if (root.most_downstream > Z)
            zmax = root.most_downstream - Z
        end
    end

    ymin = getintegral(geanz, E, zmin)
    ymax = getintegral(geanz, E, zmax)
    dy = (ymax - ymin) / n_coefficients

    if (dy > 2 * eps())

        for y = ymin+0.5dy:dy:ymax

            z = Z + getlength(geanz, E, y)
            d = sqrt(R * R + z * z)
            t1 = t + (Z - z) / C - d * params.n / C

            value += dy * E * scatteredlightfromEMshower(params, pmt, d, -z / d, θ, ϕ, t1)
        end
    end

    return value
end


"""
    directlightfromEMshowers(params::LMParameters, pmt::PMTModel, R, θ, ϕ, Δt)

Probability density function for direct light from EM-showers produced along a muon track.
Returns dP/dt [npe/ns/GeV].

# Arguments
- `params`: parameters of the setup
- `pmt`: PMT model
- `R`: (closest) distance [m] between muon and PMT
- `θ`: zenith  angle orientation PMT [rad]
- `ϕ`: azimuth angle orientation PMT [rad]
- `Δt`: time difference relative to direct Cherenkov light [ns]
"""
function directlightfromEMshowers(params::LMParameters, pmt::PMTModel, R, θ, ϕ, Δt)
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

            Ja = geant(n, ct0)  # d^2N/dcos/dϕ
            Jd = (1.0 - ng * ct0) / C  # d  t/ dz
            Je = ng * st0 * st0 * st0 / (R * C)  # d^2t/(dz)^2

            value += npe * gwater() * U * V * W * Ja / (abs(Jd) + 0.5 * Je * dz)
        end

    end
    return value
end


"""
    scatteredlightfromEMshowers(params::LMParameters, pmt::PMTModel, R, θ, ϕ, Δt)

Probability density function for scattered light from EM-showers produced along a muon track.
Returns dP/dt [npe/ns/GeV].

# Arguments
- `params`: parameters of the setup
- `pmt`: PMT model
- `R`: (closest) distance [m] between muon and PMT
- `θ`: zenith  angle orientation PMT [rad]
- `ϕ`: azimuth angle orientation PMT [rad]
- `Δt`: time difference relative to direct Cherenkov light [ns]
"""
function scatteredlightfromEMshowers(params::LMParameters, pmt::PMTModel, R, θ, ϕ, Δt)
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

                        Jb = geant(n, ct0)  # d^2N/dcos/dϕ

                        value +=
                            dom * npe * gwater() * dz * U * V * W * Ja * Jb * Jc / abs(Jd)
                    end
                end
                sb += ds

            end
        end
    end

    value
end


"""
    directlightfrombrightpoint(params::LMParameters, pmt::PMTModel, D, ct, Δt)

Probability density function for direct light from an isotropic (bright) point source.
Returns d²P/dt/dE [npe/ns/GeV].

# Arguments
- `params`: parameters of the setup
- `pmt`: PMT model
- `D`: distance between source and PMT [m]
- `ct`: cosine angle of PMT orientation
- `Δt`: time difference relative to direct Cherenkov light [ns]
"""
function directlightfrombrightpoint(params::LMParameters, pmt::PMTModel, D, ct, Δt)
    D = max(D, params.minimum_distance)
    t = D * params.n / C + Δt  # time [ns]
    A = pmt.photocathode_area

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max)
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min)
    ng = C * t / D  # index of refraction

    n0 >= ng && return 0.0
    n1 <= ng && return 0.0

    w_start = params.lambda_max + (ng - n0) / (n1 - n0) * (params.lambda_min - params.lambda_max)
    w = wavelength(params.dispersion_model, ng, w_start, 1.0e-5)
    n = refractionindexphase(params.dispersion_model, w)

    l_abs = absorptionlength(params.absorption_model, w)
    ls = scatteringlength(params.scattering_model, w)

    npe = cherenkov(w, n) * pmt.quantum_efficiency(w)

    U = pmt.angular_acceptance(ct)  # PMT angular acceptance
    V = exp(-D / l_abs - D / ls)    # absorption & scattering
    W = A / (D * D)                 # solid angle

    ngp = dispersiongroup(params.dispersion_model, w)
    Ja = D * ngp / C   # dt/dlambda
    Jb = 1.0 / (4π)    # isotropic: d^2N/dOmega

    return npe * geanc() * U * V * W * Jb / abs(Ja)
end


"""
    scatteredlightfrombrightpoint(params::LMParameters, pmt::PMTModel, D, ct, Δt)

Probability density function for scattered light from an isotropic (bright) point source.
Returns d²P/dt/dE [npe/ns/GeV].

# Arguments
- `params`: parameters of the setup
- `pmt`: PMT model
- `D`: distance between source and PMT [m]
- `ct`: cosine angle of PMT orientation
- `Δt`: time difference relative to direct Cherenkov light [ns]
"""
function scatteredlightfrombrightpoint(params::LMParameters, pmt::PMTModel, D, ct, Δt)
    value = 0.0

    D = max(D, params.minimum_distance)
    t = D * params.n / C + Δt  # time [ns]
    st = sqrt((1.0 + ct) * (1.0 - ct))
    A = pmt.photocathode_area

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max)
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min)
    ni = C * t / D  # maximal index of refraction

    n0 >= ni && return value

    nj = min(ni, n1)

    w = params.lambda_max

    n_coefficients = length(params.legendre_coefficients[1])

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

        Jc = 1.0 / ls    # dN/dx
        Jb = 1.0 / (4π)  # isotropic: d^2N/dcos/dϕ

        d = C * t / ng  # photon path

        dcb = 2.0 / (n_coefficients + 1)

        cb = -1.0 + 0.5dcb
        while cb < 1.0

            sb = sqrt((1.0 + cb) * (1.0 - cb))

            v = 0.5 * (d + D) * (d - D) / (d - D * cb)
            u = d - v

            if u > 0 && v > 0

                cts = (D * cb - v) / u  # cosine scattering angle

                V = exp(
                    -d * inverseattenuationlength(
                        params.scattering_probability_model,
                        l_abs,
                        ls,
                        cts,
                    ),
                )

                if !(cts < 0.0 && v * sqrt((1.0 + cts) * (1.0 - cts)) < params.module_radius)

                    W = min(A / (v * v), 2π)  # solid angle
                    Ja = scatteringprobability(params.scattering_probability_model, cts)  # d^2P/dcos/dϕ
                    Jd = ng * (1.0 - cts) / C  # dt/du

                    dp = π / length(params.integration_points)
                    dom = dcb * dp * v * v / (u * u)

                    for (cp, sp) in params.integration_points.xy
                        dot = cb * ct + sb * cp * st
                        U = 2 * pmt.angular_acceptance(dot)  # PMT angular acceptance (both ±φ)
                        value += npe * geanc() * dom * U * V * W * Ja * Jb * Jc / abs(Jd)
                    end
                end
            end

            cb += dcb
        end
    end

    return value
end


"""
    lightfrommuon(params, pmt, E, R, θ, ϕ, Δt)

Total light yield from a muon of energy `E` [GeV]: sum of direct and scattered
contributions from the muon track, EM showers, and delta rays. Returns dP/dt [npe/ns].
"""
function lightfrommuon(params::LMParameters, pmt::PMTModel, E, R, θ, ϕ, Δt)
    dE = deltarayenergyloss(E)
    return (directlightfrommuon(params, pmt, R, θ, ϕ, Δt)             +
            scatteredlightfrommuon(params, pmt, R, θ, ϕ, Δt)          +
            directlightfromEMshowers(params, pmt, R, θ, ϕ, Δt)    * E  +
            scatteredlightfromEMshowers(params, pmt, R, θ, ϕ, Δt) * E  +
            directlightfromdeltarays(params, pmt, R, θ, ϕ, Δt)    * dE +
            scatteredlightfromdeltarays(params, pmt, R, θ, ϕ, Δt) * dE)
end


"""
    lightfromEMshower(params, pmt, D, cd, θ, ϕ, Δt)
    lightfromEMshower(params, pmt, E, D, cd, θ, ϕ, Δt)

Total (direct + scattered) light from an EM shower. Returns d²P/dt/dE [npe/ns/GeV]
when called without energy, or dP/dt [npe/ns] when called with energy `E` [GeV].
"""
function lightfromEMshower(params::LMParameters, pmt::PMTModel, D, cd, θ, ϕ, Δt)
    return (directlightfromEMshower(params, pmt, D, cd, θ, ϕ, Δt) +
            scatteredlightfromEMshower(params, pmt, D, cd, θ, ϕ, Δt))
end

function lightfromEMshower(params::LMParameters, pmt::PMTModel, E, D, cd, θ, ϕ, Δt)
    return (directlightfromEMshower(params, pmt, E, D, cd, θ, ϕ, Δt) +
            scatteredlightfromEMshower(params, pmt, E, D, cd, θ, ϕ, Δt))
end


"""
    lightfrombrightpoint(params, pmt, D, ct, Δt)

Total (direct + scattered) light from an isotropic bright point source.
Returns d²P/dt/dE [npe/ns/GeV].
"""
function lightfrombrightpoint(params::LMParameters, pmt::PMTModel, D, ct, Δt)
    return (directlightfrombrightpoint(params, pmt, D, ct, Δt) +
            scatteredlightfrombrightpoint(params, pmt, D, ct, Δt))
end
