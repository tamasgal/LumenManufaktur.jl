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

    D = max(D, getRmin());
    t = D * getIndexOfRefraction() / C + t_ns  # time [ns]
    A = getPhotocathodeArea();

    px = sin(theta) * cos(phi);
    pz = cos(theta);

    n0 = getIndexOfRefractionGroup(wmax);
    n1 = getIndexOfRefractionGroup(wmin);
    ng = C * t / D  # index of refraction

    if (n0 >= ng) {
      return 0.0;
    end
    if (n1 <= ng) {
      return 0.0;
    end

    w = getWavelength(ng);
    n = getIndexOfRefractionPhase(w);

    l_abs = getAbsorptionLength(w);
    ls = getScatteringLength(w);

    npe = cherenkov(w, n) * getQE(w);

    ct = st0 * px + ct0 * pz  # cosine angle of incidence on PMT

    U = getAngularAcceptance(ct)  # PMT angular acceptance
    V = exp(-D / l_abs - D / ls)  # absorption & scattering
    W = A / (D * D)  # solid angle

    ngp = getDispersionGroup(w);

    Ja = D * ngp / C  # dt/dlambda
    Jb = geant(n, ct0)  # d^2N/dcos/dphi

    return npe * geanc() * U * V * W * Jb / fabs(Ja);
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
    double value = 0;

    sd = sqrt((1.0 + cd) * (1.0 - cd));
    D = max(D_m, getRmin());
    L = D;
    t = D * getIndexOfRefraction() / C + t_ns  # time [ns]

    A = getPhotocathodeArea();

    px = sin(theta) * cos(phi);
    py = sin(theta) * sin(phi);
    pz = cos(theta);

    qx = cd * px + 0 - sd * pz;
    qy = 1 * py;
    qz = sd * px + 0 + cd * pz;

    n0 = getIndexOfRefractionGroup(wmax);
    n1 = getIndexOfRefractionGroup(wmin);

    ni = C * t / L  # maximal index of refraction

    if (n0 >= ni) {
      return value;
    end

    nj = min(ni, n1);

    w = wmax;

    for (const_iterator i = begin(); i != end(); ++i) {

      ng = 0.5 * (nj + n0) + i->getX() * 0.5 * (nj - n0);
      dn = i->getY() * 0.5 * (nj - n0);

      w = getWavelength(ng, w, 1.0e-5);

      dw = dn / fabs(getDispersionGroup(w));

      n = getIndexOfRefractionPhase(w);

      l_abs = getAbsorptionLength(w);
      ls = getScatteringLength(w);

      npe = cherenkov(w, n) * dw * getQE(w);

      if (npe <= 0) {
        continue;
      end

      Jc = 1.0 / ls  # dN/dx

      d = C * t / ng  # photon path

      // V  = exp(-d/l_abs)  #
      // absorption

      ds = 2.0 / (size() + 1);

      for (double sb = 0.5 * ds; sb < 1.0 - 0.25 * ds; sb += ds) {

        for k in (-1, 1)

          cb = k * sqrt((1.0 + sb) * (1.0 - sb));
          dcb = k * ds * sb / cb;

          v = 0.5 * (d + L) * (d - L) / (d - L * cb);
          u = d - v;

          if (u <= 0) {
            continue;
          end
          if (v <= 0) {
            continue;
          end

          cts = (L * cb - v) / u  # cosine scattering angle

          V =
              exp(-d * getInverseAttenuationLength(l_abs, ls, cts));

          if (cts < 0.0 &&
              v * sqrt((1.0 + cts) * (1.0 - cts)) < MODULE_RADIUS_M) {
            continue;
          end

          W = min(A / (v * v), 2.0 * PI)  # solid angle
          Ja = getScatteringProbability(cts)  # d^2P/dcos/dphi
          Jd = ng * (1.0 - cts) / C  # dt/du

          ca = (L - v * cb) / u;
          sa = v * sb / u;

          dp = PI / phd.size();
          dom = dcb * dp * v * v / (u * u);
          dos = sqrt(dom);

          for (const_iterator l = phd.begin(); l != phd.end(); ++l) {

            cp = l->getX();
            sp = l->getY();

            ct0 = cd * ca - sd * sa * cp;

            vx = -sb * cp * qx;
            vy = -sb * sp * qy;
            vz = cb * qz;

            U =
                getAngularAcceptance(vx + vy + vz) +
                getAngularAcceptance(vx - vy + vz)  # PMT angular acceptance

            // Jb = geant(n,ct0)  #
            // d^2N/dcos/dphi

            // value += npe * geanc() * dom * U * V * W * Ja * Jb * Jc /
            // fabs(Jd);

            Jb = geant(n, ct0 - 0.5 * dos,
                                    ct0 + 0.5 * dos)  # dN/dphi

            value += npe * geanc() * dos * U * V * W * Ja * Jb * Jc / fabs(Jd);
          end
        end
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
    double value = 0;

    sd = sqrt((1.0 + cd) * (1.0 - cd));
    D = max(D_m, getRmin());
    R = D * sd  # minimal distance of approach [m]
    Z = -D * cd;
    t = D * getIndexOfRefraction() / C + t_ns  # time [ns]

    n0 = getIndexOfRefractionGroup(wmax);
    n1 = getIndexOfRefractionGroup(wmin);

    double zmin = 0.0  # minimal shower length [m]
    double zmax = geanz.getLength(E, 1.0)  # maximal shower length [m]

    d = sqrt((Z + zmax) * (Z + zmax) + R * R);

    if (C * t > max(n1 * D, zmax + n1 * d)) {
      return value;
    end

    if (C * t < n0 * D) {

      JRoot rz(R, n0, t + Z / C)  # square root

      if (!rz.is_valid) {
        return value;
      end

      if (rz.second > Z) {
        zmin = rz.second - Z;
      end
      if (rz.first > Z) {
        zmin = rz.first - Z;
      end
    end

    if (C * t > n1 * D) {

      JRoot rz(R, n1, t + Z / C)  # square root

      if (!rz.is_valid) {
        return value;
      end

      if (rz.second > Z) {
        zmin = rz.second - Z;
      end
      if (rz.first > Z) {
        zmin = rz.first - Z;
      end
    end

    if (C * t < zmax + n0 * d) {

      JRoot rz(R, n0, t + Z / C)  # square root

      if (!rz.is_valid) {
        return value;
      end

      if (rz.first > Z) {
        zmax = rz.first - Z;
      end
      if (rz.second > Z) {
        zmax = rz.second - Z;
      end
    end

    if (zmin < 0.0) {
      zmin = 0.0;
    end

    if (zmax > zmin) {

      ymin = geanz.getIntegral(E, zmin);
      ymax = geanz.getIntegral(E, zmax);
      dy = (ymax - ymin) / size();

      if dy > 2 * eps()

        for (double y = ymin + 0.5 * dy; y < ymax; y += dy) {

          z = Z + geanz.getLength(E, y);
          d = sqrt(R * R + z * z);
          t1 = t + (Z - z) / C - d * getIndexOfRefraction() / C;

          value +=
              dy * E * getDirectLightFromEMshower(d, -z / d, theta, phi, t1);
        end
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
    D = max(D_m, getRmin());
    R = D * sd  # minimal distance of approach [m]
    Z = -D * cd;
    t = D * getIndexOfRefraction() / C + t_ns  # time [ns]

    n0 = getIndexOfRefractionGroup(wmax);

    zmin = 0.0  # minimal shower length [m]
    zmax = geanz.getLength(E, 1.0)  # maximal shower length [m]

    d = sqrt((Z + zmax) * (Z + zmax) + R * R);

    if (C * t < n0 * D) {

      JRoot rz(R, n0, t + Z / C)  # square root

      if (!rz.is_valid) {
        return value;
      end

      if (rz.second > Z) {
        zmin = rz.second - Z;
      end
      if (rz.first > Z) {
        zmin = rz.first - Z;
      end
    end

    if (C * t < zmax + n0 * d) {

      JRoot rz(R, n0, t + Z / C)  # square root

      if (!rz.is_valid) {
        return value;
      end

      if (rz.first > Z) {
        zmax = rz.first - Z;
      end
      if (rz.second > Z) {
        zmax = rz.second - Z;
      end
    end

    ymin = geanz.getIntegral(E, zmin);
    ymax = geanz.getIntegral(E, zmax);
    dy = (ymax - ymin) / size();

    if (dy > 2 * eps()) {

      for (double y = ymin + 0.5 * dy; y < ymax; y += dy) {

        z = Z + geanz.getLength(E, y);
        d = sqrt(R * R + z * z);
        t1 = t + (Z - z) / C - d * getIndexOfRefraction() / C;

        value +=
            dy * E * getScatteredLightFromEMshower(d, -z / d, theta, phi, t1);
      end
    end

    return value;
  end
