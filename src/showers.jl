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
    ct0 = cd;
    st0 = sqrt((1.0 + ct0) * (1.0 - ct0));

    D = max(D, params.minimum_distance);
    t = D * params.n / C + Δt  # time [ns]
    A = pmt.photocathode_area;

    px = sin(θ) * cos(ϕ);
    pz = cos(θ);

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max);
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min);
    ng = C * t / D  # index of refraction

    if (n0 >= ng)
      return 0.0;
    end
    if (n1 <= ng)
      return 0.0;
    end

    w = wavelength(params.dispersion_model, ng)
    n = refractionindexphase(params.dispersion_model, w);

    l_abs = absorptionlength(params.absorption_model, w);
    ls = scatteringlength(params.scattering_model, w);

    npe = cherenkov(w, n) * getQE(w);

    ct = st0 * px + ct0 * pz  # cosine angle of incidence on PMT

    U = pmt.angular_acceptance(ct)  # PMT angular acceptance
    V = exp(-D / l_abs - D / ls)  # absorption & scattering
    W = A / (D * D)  # solid angle

    ngp = getDispersionGroup(w);

    Ja = D * ngp / C  # dt/dlambda
    Jb = geant(n, ct0)  # d^2N/dcos/dϕ

    return npe * geanc() * U * V * W * Jb / abs(Ja);
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

    sd = sqrt((1.0 + cd) * (1.0 - cd));
    D = max(D, params.minimum_distance);
    L = D;
    t = D * params.n / C + Δt  # time [ns]

    A = pmt.photocathode_area;

    px = sin(θ) * cos(ϕ);
    py = sin(θ) * sin(ϕ);
    pz = cos(θ);

    qx = cd * px + 0 - sd * pz;
    qy = 1 * py;
    qz = sd * px + 0 + cd * pz;

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max);
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min);

    ni = C * t / L  # maximal index of refraction

    if (n0 >= ni)
      return value;
    end

    nj = min(ni, n1);

    w = wmax;

    n_coefficients = length(params.legendre_coefficients[1])

    for (m_x, m_y) in zip(params.legendre_coefficients...)

      ng = 0.5 * (nj + n0) + m_x * 0.5 * (nj - n0);
      dn = m_y * 0.5 * (nj - n0);

      w = getWavelength(ng, w, 1.0e-5);

      dw = dn / abs(getDispersionGroup(w));

      n = refractionindexphase(params.dispersion_model, w);

      l_abs = absorptionlength(params.absorption_model, w);
      ls = scatteringlength(params.scattering_model, w);

      npe = cherenkov(w, n) * dw * getQE(w);

      if (npe <= 0)
        continue;
      end

      Jc = 1.0 / ls  # dN/dx

      d = C * t / ng  # photon path

      ds = 2.0 / (n_coefficients + 1);

      sb = 0.5ds
      while sb < 1.0 - 0.25ds

        for k in (-1, 1)

          cb = k * sqrt((1.0 + sb) * (1.0 - sb));
          dcb = k * ds * sb / cb;

          v = 0.5 * (d + L) * (d - L) / (d - L * cb);
          u = d - v;

          if (u <= 0)
            continue;
          end
          if (v <= 0)
            continue;
          end

          cts = (L * cb - v) / u  # cosine scattering angle

          V = exp(-d * inverseattenuationlength(params.scattering_probability_model, l_abs, ls, cts));

          cts < 0.0 && v * sqrt((1.0 + cts) * (1.0 - cts)) < params.module_radius && continue

          W = min(A / (v * v), 2.0 * PI)  # solid angle
          Ja = scatteringprobability(params.scattering_probability_model, cts)  # d^2P/dcos/dϕ
          Jd = ng * (1.0 - cts) / C  # dt/du

          ca = (L - v * cb) / u;
          sa = v * sb / u;

          dp = π / length(params.integration_points)
          dom = dcb * dp * v * v / (u * u);
          dos = sqrt(dom);

          for (cp, sp) in params.integration_points.xy

            ct0 = cd * ca - sd * sa * cp;

            vx = -sb * cp * qx;
            vy = -sb * sp * qy;
            vz = cb * qz;

            U =
                pmt.angular_acceptance(vx + vy + vz) +
                pmt.angular_acceptance(vx - vy + vz)  # PMT angular acceptance

            Jb = geant(n, ct0 - 0.5 * dos, ct0 + 0.5 * dos)  # dN/dϕ

            value += npe * geanc() * dos * U * V * W * Ja * Jb * Jc / abs(Jd);
          end
        end
        sb += ds
      end
    end

    return value;
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

    sd = sqrt((1.0 + cd) * (1.0 - cd));
    D = max(D, params.minimum_distance);
    R = D * sd  # minimal distance of approach [m]
    Z = -D * cd;
    t = D * params.n / C + Δt  # time [ns]

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max);
    n1 = refractionindexgroup(params.dispersion_model, params.lambda_min);

    zmin = 0.0  # minimal shower length [m]
    zmax = geanz.getLength(E, 1.0)  # maximal shower length [m]

    d = sqrt((Z + zmax) * (Z + zmax) + R * R);

    if (C * t > max(n1 * D, zmax + n1 * d))
      return value;
    end

    if (C * t < n0 * D)
      root = Root(R, n0, t + Z / C)  # square root
      !root.isvalid && return value
      if root.most_downstream > Z
        zmin = root.most_downstream - Z;
      end
      if root.most_upstream > Z
        zmin = root.most_upstream - Z;
      end
    end

    if (C * t > n1 * D)
      root = Root(R, n1, t + Z / C)  # square root

      !root.isvalid && return value

      if (root.most_downstream > Z)
        zmin = root.most_downstream - Z;
      end
      if (root.most_upstream > Z)
        zmin = root.most_upwnstream - Z;
      end
    end

    if (C * t < zmax + n0 * d)

      root = Root(R, n0, t + Z / C)  # square root

      !root.isvalid && return value

      if (root.most_upstream > Z)
        zmax = root.most_upstream - Z;
      end
      if (root.most_downstream > Z)
        zmax = root.most_downstream - Z;
      end
    end

    if (zmin < 0.0)
      zmin = 0.0;
    end

    if (zmax > zmin)

      ymin = geanz.getIntegral(E, zmin);
      ymax = geanz.getIntegral(E, zmax);
      dy = (ymax - ymin) / size();

      if dy > 2 * eps()

        # for (double y = ymin + 0.5 * dy; y < ymax; y += dy)

        #   z = Z + geanz.getLength(E, y);
        #   d = sqrt(R * R + z * z);
        #   t1 = t + (Z - z) / C - d * params.n / C;

        #   value += dy * E * directlightfromEMshower(params, pmt, d, -z / d, θ, ϕ, t1);
        # end
      end
    end

    return value;
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

    sd = sqrt((1.0 + cd) * (1.0 - cd));
    D = max(D, params.minimum_distance);
    R = D * sd  # minimal distance of approach [m]
    Z = -D * cd;
    t = D * params.n / C + Δt  # time [ns]

    n0 = refractionindexgroup(params.dispersion_model, params.lambda_max);

    zmin = 0.0  # minimal shower length [m]
    zmax = geanz.getLength(E, 1.0)  # maximal shower length [m]

    d = sqrt((Z + zmax) * (Z + zmax) + R * R);

    if (C * t < n0 * D)

      root = Root(R, n0, t + Z / C)  # square root

      !root.isvalid && return value

      if (root.most_downstream > Z)
        zmin = root.most_downstream - Z;
      end
      if (root.most_upstream > Z)
        zmin = root.most_upstream - Z;
      end
    end

    if (C * t < zmax + n0 * d)

      root = Root(R, n0, t + Z / C)  # square root

      !root.isvalid && return value

      if (root.most_upstream > Z)
        zmax = root.most_upstream - Z;
      end
      if (root.most_downstream > Z)
        zmax = root.most_downstream - Z;
      end
    end

    ymin = geanz.getIntegral(E, zmin);
    ymax = geanz.getIntegral(E, zmax);
    dy = (ymax - ymin) / size();

    if (dy > 2 * eps())

      for (double y = ymin + 0.5 * dy; y < ymax; y += dy)

        z = Z + geanz.getLength(E, y);
        d = sqrt(R * R + z * z);
        t1 = t + (Z - z) / C - d * params.n / C;

        value += dy * E * scatteredlightfromEMshower(params, pmt, d, -z / d, θ, ϕ, t1);
      end
    end

    return value;
  end
