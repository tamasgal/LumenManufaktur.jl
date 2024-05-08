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

- `R`: (closest) distance [m] between muon and PMT
- `ϕ`: zenith angle [rad] which is 0 when the PMT points away from the muon
  track (in x-direction) and rotates counter clockwise to the y-axis when
  viewed from above, while the z-axis points upwards.
- `θ`: azimuth angle [rad] which is 0 when the PMT points upwards (along the
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
closest distance of `R` [m] to the PMT, the angles `θ` [rad] (zenith) and `ϕ`
[rad] (azimuth) with respect to the PMT axis and a time difference `Δt` [ns]
relative to direct Cherenkov light. Returns dP/dt [npe/ns].

# Arguments

- `R`: (closest) distance [m] between muon and PMT
- `ϕ`: zenith angle [rad] which is 0 when the PMT points away from the muon
  track (in x-direction) and rotates counter clockwise to the y-axis when
  viewed from above, while the z-axis points upwards.
- `θ`: azimuth angle [rad] which is 0 when the PMT points upwards (along the
  z-axis) and π/2 when pointing downwards. Between 0 and π/2, the PMT points to
  the z-axis
- `Δt`: time difference [ns] relative to the Cherenkov light
"""
function directlightfrommuon(params::LMParameters, pmt::PMTModel, R, θ, ϕ, Δt)
    N = 100      # maximal number of iterations
    eps = 1.0e-6 # precision index of refraction

    R = max(R, params.minimum_distance)
    A = pmt.photocathode_area
    t = R * tanthetac(params.n) / C + Δt; # time [ns]
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

    for _ in 1:N  # binary search for wavelength

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
