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
    const double ct0 = cd;
    const double st0 = sqrt((1.0 + ct0) * (1.0 - ct0));

    const double D = std::max(D, getRmin());
    const double t = D * getIndexOfRefraction() / C + t_ns; // time [ns]
    const double A = getPhotocathodeArea();

    const double px = sin(theta) * cos(phi);
    const double pz = cos(theta);

    const double n0 = getIndexOfRefractionGroup(wmax);
    const double n1 = getIndexOfRefractionGroup(wmin);
    const double ng = C * t / D; // index of refraction

    if (n0 >= ng) {
      return 0.0;
    }
    if (n1 <= ng) {
      return 0.0;
    }

    const double w = getWavelength(ng);
    const double n = getIndexOfRefractionPhase(w);

    const double l_abs = getAbsorptionLength(w);
    const double ls = getScatteringLength(w);

    const double npe = cherenkov(w, n) * getQE(w);

    const double ct = st0 * px + ct0 * pz; // cosine angle of incidence on PMT

    const double U = getAngularAcceptance(ct); // PMT angular acceptance
    const double V = exp(-D / l_abs - D / ls); // absorption & scattering
    const double W = A / (D * D);              // solid angle

    const double ngp = getDispersionGroup(w);

    const double Ja = D * ngp / C;   // dt/dlambda
    const double Jb = geant(n, ct0); // d^2N/dcos/dphi

    return npe * geanc() * U * V * W * Jb / fabs(Ja);
  }

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
    double value = 0;

    const double sd = sqrt((1.0 + cd) * (1.0 - cd));
    const double D = std::max(D_m, getRmin());
    const double L = D;
    const double t = D * getIndexOfRefraction() / C + t_ns; // time [ns]

    const double A = getPhotocathodeArea();

    const double px = sin(theta) * cos(phi);
    const double py = sin(theta) * sin(phi);
    const double pz = cos(theta);

    const double qx = cd * px + 0 - sd * pz;
    const double qy = 1 * py;
    const double qz = sd * px + 0 + cd * pz;

    const double n0 = getIndexOfRefractionGroup(wmax);
    const double n1 = getIndexOfRefractionGroup(wmin);

    const double ni = C * t / L; // maximal index of refraction

    if (n0 >= ni) {
      return value;
    }

    const double nj = std::min(ni, n1);

    double w = wmax;

    for (const_iterator i = begin(); i != end(); ++i) {

      const double ng = 0.5 * (nj + n0) + i->getX() * 0.5 * (nj - n0);
      const double dn = i->getY() * 0.5 * (nj - n0);

      w = getWavelength(ng, w, 1.0e-5);

      const double dw = dn / fabs(getDispersionGroup(w));

      const double n = getIndexOfRefractionPhase(w);

      const double l_abs = getAbsorptionLength(w);
      const double ls = getScatteringLength(w);

      const double npe = cherenkov(w, n) * dw * getQE(w);

      if (npe <= 0) {
        continue;
      }

      const double Jc = 1.0 / ls; // dN/dx

      const double d = C * t / ng; // photon path

      // const double V  = exp(-d/l_abs);                           //
      // absorption

      const double ds = 2.0 / (size() + 1);

      for (double sb = 0.5 * ds; sb < 1.0 - 0.25 * ds; sb += ds) {

        for (int buffer[] = {-1, +1, 0}, *k = buffer; *k != 0; ++k) {

          const double cb = (*k) * sqrt((1.0 + sb) * (1.0 - sb));
          const double dcb = (*k) * ds * sb / cb;

          const double v = 0.5 * (d + L) * (d - L) / (d - L * cb);
          const double u = d - v;

          if (u <= 0) {
            continue;
          }
          if (v <= 0) {
            continue;
          }

          const double cts = (L * cb - v) / u; // cosine scattering angle

          const double V =
              exp(-d * getInverseAttenuationLength(l_abs, ls, cts));

          if (cts < 0.0 &&
              v * sqrt((1.0 + cts) * (1.0 - cts)) < MODULE_RADIUS_M) {
            continue;
          }

          const double W = std::min(A / (v * v), 2.0 * PI); // solid angle
          const double Ja = getScatteringProbability(cts);  // d^2P/dcos/dphi
          const double Jd = ng * (1.0 - cts) / C;           // dt/du

          const double ca = (L - v * cb) / u;
          const double sa = v * sb / u;

          const double dp = PI / phd.size();
          const double dom = dcb * dp * v * v / (u * u);
          const double dos = sqrt(dom);

          for (const_iterator l = phd.begin(); l != phd.end(); ++l) {

            const double cp = l->getX();
            const double sp = l->getY();

            const double ct0 = cd * ca - sd * sa * cp;

            const double vx = -sb * cp * qx;
            const double vy = -sb * sp * qy;
            const double vz = cb * qz;

            const double U =
                getAngularAcceptance(vx + vy + vz) +
                getAngularAcceptance(vx - vy + vz); // PMT angular acceptance

            // const double Jb = geant(n,ct0);                      //
            // d^2N/dcos/dphi

            // value += npe * geanc() * dom * U * V * W * Ja * Jb * Jc /
            // fabs(Jd);

            const double Jb = geant(n, ct0 - 0.5 * dos,
                                    ct0 + 0.5 * dos); // dN/dphi

            value += npe * geanc() * dos * U * V * W * Ja * Jb * Jc / fabs(Jd);
          }
        }
      }
    }

    return value;
  }

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
    double value = 0;

    const double sd = sqrt((1.0 + cd) * (1.0 - cd));
    const double D = std::max(D_m, getRmin());
    const double R = D * sd; // minimal distance of approach [m]
    const double Z = -D * cd;
    const double t = D * getIndexOfRefraction() / C + t_ns; // time [ns]

    const double n0 = getIndexOfRefractionGroup(wmax);
    const double n1 = getIndexOfRefractionGroup(wmin);

    double zmin = 0.0;                     // minimal shower length [m]
    double zmax = geanz.getLength(E, 1.0); // maximal shower length [m]

    const double d = sqrt((Z + zmax) * (Z + zmax) + R * R);

    if (C * t > std::max(n1 * D, zmax + n1 * d)) {
      return value;
    }

    if (C * t < n0 * D) {

      JRoot rz(R, n0, t + Z / C); // square root

      if (!rz.is_valid) {
        return value;
      }

      if (rz.second > Z) {
        zmin = rz.second - Z;
      }
      if (rz.first > Z) {
        zmin = rz.first - Z;
      }
    }

    if (C * t > n1 * D) {

      JRoot rz(R, n1, t + Z / C); // square root

      if (!rz.is_valid) {
        return value;
      }

      if (rz.second > Z) {
        zmin = rz.second - Z;
      }
      if (rz.first > Z) {
        zmin = rz.first - Z;
      }
    }

    if (C * t < zmax + n0 * d) {

      JRoot rz(R, n0, t + Z / C); // square root

      if (!rz.is_valid) {
        return value;
      }

      if (rz.first > Z) {
        zmax = rz.first - Z;
      }
      if (rz.second > Z) {
        zmax = rz.second - Z;
      }
    }

    if (zmin < 0.0) {
      zmin = 0.0;
    }

    if (zmax > zmin) {

      const double ymin = geanz.getIntegral(E, zmin);
      const double ymax = geanz.getIntegral(E, zmax);
      const double dy = (ymax - ymin) / size();

      if (dy > 2 * std::numeric_limits<double>::epsilon()) {

        for (double y = ymin + 0.5 * dy; y < ymax; y += dy) {

          const double z = Z + geanz.getLength(E, y);
          const double d = sqrt(R * R + z * z);
          const double t1 = t + (Z - z) / C - d * getIndexOfRefraction() / C;

          value +=
              dy * E * getDirectLightFromEMshower(d, -z / d, theta, phi, t1);
        }
      }
    }

    return value;
  }

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
  double getScatteredLightFromEMshower(const double E, const double D_m,
                                       const double cd, const double theta,
                                       const double phi,
                                       const double t_ns) const {
    double value = 0;

    const double sd = sqrt((1.0 + cd) * (1.0 - cd));
    const double D = std::max(D_m, getRmin());
    const double R = D * sd; // minimal distance of approach [m]
    const double Z = -D * cd;
    const double t = D * getIndexOfRefraction() / C + t_ns; // time [ns]

    const double n0 = getIndexOfRefractionGroup(wmax);

    double zmin = 0.0;                     // minimal shower length [m]
    double zmax = geanz.getLength(E, 1.0); // maximal shower length [m]

    const double d = sqrt((Z + zmax) * (Z + zmax) + R * R);

    if (C * t < n0 * D) {

      JRoot rz(R, n0, t + Z / C); // square root

      if (!rz.is_valid) {
        return value;
      }

      if (rz.second > Z) {
        zmin = rz.second - Z;
      }
      if (rz.first > Z) {
        zmin = rz.first - Z;
      }
    }

    if (C * t < zmax + n0 * d) {

      JRoot rz(R, n0, t + Z / C); // square root

      if (!rz.is_valid) {
        return value;
      }

      if (rz.first > Z) {
        zmax = rz.first - Z;
      }
      if (rz.second > Z) {
        zmax = rz.second - Z;
      }
    }

    const double ymin = geanz.getIntegral(E, zmin);
    const double ymax = geanz.getIntegral(E, zmax);
    const double dy = (ymax - ymin) / size();

    if (dy > 2 * std::numeric_limits<double>::epsilon()) {

      for (double y = ymin + 0.5 * dy; y < ymax; y += dy) {

        const double z = Z + geanz.getLength(E, y);
        const double d = sqrt(R * R + z * z);
        const double t1 = t + (Z - z) / C - d * getIndexOfRefraction() / C;

        value +=
            dy * E * getScatteredLightFromEMshower(d, -z / d, theta, phi, t1);
      }
    }

    return value;
  }
