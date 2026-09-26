// Magnetohydrodynamics: resistive decay against exp(-eta k^2 t), the Alfven period 2L / v_A,
// div B = 0 to rounding, the magnetopause at the Chapman-Ferraro distance, and the tokamak's
// fields and kink against the energy principle.
#include "TestRunner.h"
#include "Tests.h"

void testMagneticField() {
    // The walls are perfect conductors: the magnetic flux through them is frozen. The analytic
    // checks are therefore set up away from the walls the field crosses.
    auto faceVelocities = [](Field3& u, Field3& v, Field3& w, int nx, int ny, int nz) {
        u.init(nx + 1, ny, nz, {0, 0.5f, 0.5f});
        v.init(nx, ny + 1, nz, {0.5f, 0, 0.5f});
        w.init(nx, ny, nz + 1, {0.5f, 0.5f, 0});
    };
    Field3 u, v, w;

    // 1. Resistive diffusion at rest: Bz = B1 cos(pi x / L) decays as exp(-eta k^2 t) (k^2 of the
    //    discrete Laplacian, (2 - 2 cos(k dx)) / dx^2), measured mid-way between the z walls -
    //    0.5 m from them, far beyond the diffusion length sqrt(eta t) = 0.13 m.
    {
        const int nx = 32, ny = 4, nz = 32;
        const float dx = 1.0f / nx, kx = kPi;
        MagneticField m;
        m.conductivity = 1e5f; // eta = 7.96 m^2/s
        m.numericalDissipation = 0;
        m.reset(nx, ny, nz, dx, Vector3(0.0f));
        for (int k = 0; k <= nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) m.bz.at(i, j, k) = 1e-3f * std::cos(kx * (i + 0.5f) * dx);
        faceVelocities(u, v, w, nx, ny, nz);
        const float eta = m.resistivity(), t = 2e-3f;
        m.induce(u, v, w, 1.0f, t);
        double amp = 0, norm = 0;
        for (int i = 0; i < nx; ++i) {
            const double c = std::cos(kx * (i + 0.5f) * dx);
            amp += m.bz.at(i, 1, nz / 2) * c;
            norm += c * c;
        }
        amp /= norm * 1e-3;
        const double k2 = (2.0 - 2.0 * std::cos(kx * dx)) / (dx * dx), expected = std::exp(-eta * k2 * t);
        std::printf("  resistive decay: B/B0 %.5f after %.1f ms, exp(-eta k^2 t) = %.5f\n", amp, t * 1e3, expected);
        CHECK(std::fabs(amp / expected - 1) < 0.01, "resistive decay %f vs %f", amp, expected);
    }

    // 2. Standing torsional Alfven wave along B0 = x: a swirl (a vortex tube inside the cross-section,
    //    clear of the walls) twists the field lines, their tension untwists it - the Lorentz force and
    //    Faraday's law exchange kinetic and magnetic energy with period 2L / v_A,
    //    v_A = B0 / sqrt(mu0 rho), whatever the swirl's profile (torsional waves do not disperse).
    {
        const int nx = 32, ny = 16, nz = 16;
        const float L = 1.0f, dx = L / nx, kx = kPi / L, a = 0.2f, yc = 0.25f, zc = 0.25f;
        MagneticField m;
        m.applied = {0.01f, 0, 0};
        m.conductivity = 1e12f; // ideal
        m.numericalDissipation = 0;
        m.reset(nx, ny, nz, dx, Vector3(0.0f));
        faceVelocities(u, v, w, nx, ny, nz);
        const float rho = 1.0f, U = 0.05f;
        auto swirl = [&](float y, float z) { // angular velocity of the vortex tube
            const float r2 = (sqr(y - yc) + sqr(z - zc)) / (a * a);
            return r2 < 1 ? U / a * sqr(1 - r2) : 0.0f;
        };
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j <= ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    const float x = (i + 0.5f) * dx, y = j * dx, z = (k + 0.5f) * dx;
                    v.at(i, j, k) = -swirl(y, z) * (z - zc) * std::sin(kx * x);
                }
        for (int k = 0; k <= nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    const float x = (i + 0.5f) * dx, y = (j + 0.5f) * dx, z = k * dx;
                    w.at(i, j, k) = swirl(y, z) * (y - yc) * std::sin(kx * x);
                }
        auto kinetic = [&] {
            double e = 0;
            for (float s : v.d) e += s * s;
            for (float s : w.d) e += s * s;
            return 0.5 * rho * e * dx * dx * dx;
        };
        const float vA = 0.01f / std::sqrt(MagneticField::kMu0 * rho), period = 2 * L / vA;
        const float dt = 0.2f * dx / vA;
        const double k0 = kinetic(), m0 = m.energy();
        const int probeJ = int(yc / dx), probeK = int((zc + 0.5f * a) / dx); // a point in the swirl
        float t = 0, lastCross = -1, firstCross = -1, prev = v.at(nx / 2, probeJ, probeK);
        int crossings = 0;
        std::vector<uint8_t> none;
        while (t < 2.2f * period) {
            m.applyLorentzForce(u, v, w, none, rho, dt);
            m.induce(u, v, w, rho, dt);
            t += dt;
            const float now = v.at(nx / 2, probeJ, probeK);
            if ((prev > 0) != (now > 0)) {
                const float tc = t - dt * now / (now - prev);
                if (firstCross < 0) firstCross = tc;
                lastCross = tc;
                ++crossings;
            }
            prev = now;
        }
        const float measured = crossings > 1 ? 2 * (lastCross - firstCross) / float(crossings - 1) : 0;
        const double energyKept = (kinetic() + m.energy() - m0) / k0;
        std::printf("  torsional Alfven wave: v_A %.3f m/s, period %.4f s vs 2L/v_A %.4f s (%d zero crossings); energy kept %.1f %%\n",
                    vA, measured, period, crossings, 100 * energyKept);
        CHECK(std::fabs(measured / period - 1) < 0.02f, "Alfven period %f vs %f", measured, period);
        // The staggered interpolations of B and J make the scheme slightly dissipative (~10 % over two
        // periods here); what must never happen is energy growth (an instability).
        CHECK(energyKept > 0.8 && energyKept < 1.01, "ideal MHD energy: no growth, little loss (%f)", energyKept);
    }

    // 3. div B stays zero (to rounding) while a vortex winds up the field lines of an ideal
    //    conductor (at sigma = 1e6 S/m the field would just diffuse out of a 1 m box: eta = 0.8 m^2/s).
    {
        MagneticField m;
        m.conductivity = 1e9f;
        m.applied = {0.02f, 0.005f, 0.0f};
        m.reset(24, 24, 24, 1.0f / 24, Vector3(0.0f));
        Field3 uu, vv, ww;
        uu.init(25, 24, 24, {0, 0.5f, 0.5f});
        vv.init(24, 25, 24, {0.5f, 0, 0.5f});
        ww.init(24, 24, 25, {0.5f, 0.5f, 0});
        for (int k = 0; k < 24; ++k)
            for (int j = 0; j < 24; ++j)
                for (int i = 1; i < 24; ++i) uu.at(i, j, k) = -std::sin(kPi * i / 24.0f) * std::cos(kPi * (j + 0.5f) / 24.0f);
        for (int k = 0; k < 24; ++k)
            for (int j = 1; j < 24; ++j)
                for (int i = 0; i < 24; ++i) vv.at(i, j, k) = std::cos(kPi * (i + 0.5f) / 24.0f) * std::sin(kPi * j / 24.0f);
        for (int s = 0; s < 200; ++s) m.induce(uu, vv, ww, 1.0f, 2e-3f);
        std::printf("  vortex winding the field: max |B| %.4f T (from %.4f), max |div B| dx / |B| %.1e\n", m.maxField(),
                    std::sqrt(0.02f * 0.02f + 0.005f * 0.005f), m.maxDivergence());
        CHECK(m.maxDivergence() < 1e-5f, "div B not zero (%e)", m.maxDivergence());
        CHECK(m.maxField() > 1.1f * std::sqrt(0.02f * 0.02f + 0.005f * 0.005f), "the vortex must wind up (stretch) the field");
    }
}

void testMagnetosphere() {
    // A plasma wind against a magnetised sphere: it must stop at the Chapman-Ferraro distance,
    // where the dipole's magnetic pressure balances the ram pressure - 0.29 m (pressure balance)
    // to 0.37 m (field doubled by the magnetopause currents) - and flow around.
    Simulation sim;
    sim.loadPreset(Preset::Magnetosphere);
    GasSolver& g = sim.grid;
    while (g.time() < 1.0f) sim.stepFrame();
    float standoff = -1; // where, coming from upstream, the flow has slowed to half the wind speed
    for (float x = -0.16f; x > g.domain().lo.x + 0.05f; x -= 0.01f)
        if (length(g.velocityAt({x, 0.0f, 0.0f})) > 0.5f * g.params.inflowSpeed) { standoff = -x; break; }
    const Vector3 flank = g.velocityAt({-0.3f, 0.0f, 0.35f});
    std::printf("  magnetosphere after %.2f s: flow slowed to half at %.2f m upstream (theory 0.29-0.37 m), sideways flow at the flank "
                "%.2f m/s, div B %.1e\n",
                g.time(), standoff, flank.z, g.magnetic.maxDivergence());
    CHECK(standoff > 0.27f && standoff < 0.40f, "magnetopause distance %f", standoff);
    CHECK(flank.z > 0.2f, "the plasma must be deflected around the magnet (%f)", flank.z);
    CHECK(g.magnetic.maxDivergence() < 1e-5f && std::isfinite(g.magnetic.energy()), "field broken");
}

void testTokamak() {
    // 1) The fields on the grid follow the formulas: the coils' B_phi = B0 R0 / R; the plasma
    //    current by Ampere's law around the channel (the loop potential is exact, so the torus
    //    changes nothing there); B_theta above the axis against the straight-column value (the
    //    torus bends it by ~r / R0); the single loop's potential against its on-axis field.
    Simulation sim;
    sim.loadPreset(Preset::Tokamak);
    {
        const Tokamak& t = sim.tokamak;
        const MagneticField& m = sim.grid.magnetic;
        const float R0 = t.majorRadius, a = t.minorRadius, dx = sim.grid.dx();
        auto toroidal = [&](float R) { return m.fieldAt(t.centre + Vector3(R, 0.0f, 0.0f)).z; }; // phi_hat = +z at phi = 0
        auto poloidal = [&](float r) { return -m.fieldAt(t.centre + Vector3(R0, r, 0.0f)).x; };  // above the axis B_theta = -B_x
        const float Ip = t.measuredCurrent(m, dx);
        // One loop: B on its axis is mu0 I / (2 Rl); from the potential, B_y = (1/R) d(R A)/dR.
        const double Rl = 0.4, h = 1e-4, Ay = (Tokamak::loopPotential(2 * h, 0.1, Rl, 0, 100.0) * 2 * h - Tokamak::loopPotential(h, 0.1, Rl, 0, 100.0) * h) / (h * 1.5 * h);
        const double Bexact = MagneticField::kMu0 * 100.0 * Rl * Rl / (2 * std::pow(Rl * Rl + 0.01, 1.5));
        std::printf("  tokamak: B_phi on the axis %.3f mT (B0 %.3f), at R0 -/+ a: %.3f / %.3f mT (theory %.3f / %.3f); B_theta above the "
                    "axis at a/2 %.3f mT (column %.3f), at 1.5 a %.3f (column %.3f); I_p by Ampere %.1f A (set %.1f), B_v %.3f mT, "
                    "q_a %.2f; loop potential -> on-axis B %.4g vs exact %.4g T\n",
                    toroidal(R0) * 1000, t.toroidalField * 1000, toroidal(R0 - a) * 1000, toroidal(R0 + a) * 1000,
                    t.toroidalField * R0 / (R0 - a) * 1000, t.toroidalField * R0 / (R0 + a) * 1000, poloidal(0.5f * a) * 1000,
                    t.poloidalField(0.5f * a) * 1000, poloidal(1.5f * a) * 1000, t.poloidalField(1.5f * a) * 1000, Ip,
                    t.plasmaCurrent(), t.verticalFieldStrength() * 1000, t.safetyFactorEdge, Ay, Bexact);
        CHECK(std::fabs(toroidal(R0) / t.toroidalField - 1) < 0.03f, "B_phi on the axis %f", toroidal(R0));
        CHECK(std::fabs(toroidal(R0 - a) * (R0 - a) / (t.toroidalField * R0) - 1) < 0.03f, "B_phi is not B0 R0 / R inside");
        CHECK(std::fabs(toroidal(R0 + a) * (R0 + a) / (t.toroidalField * R0) - 1) < 0.03f, "B_phi is not B0 R0 / R outside");
        CHECK(std::fabs(Ay / Bexact - 1) < 1e-3, "loop potential %g vs %g", Ay, Bexact);
        CHECK(std::fabs(poloidal(0.5f * a) / t.poloidalField(0.5f * a) - 1) < 0.25f, "B_theta(a/2) %f", poloidal(0.5f * a));
        CHECK(std::fabs(poloidal(1.5f * a) / t.poloidalField(1.5f * a) - 1) < 0.25f, "B_theta(1.5a) %f", poloidal(1.5f * a));
        CHECK(std::fabs(Ip / t.plasmaCurrent() - 1) < 0.05f, "plasma current %f vs %f", Ip, t.plasmaCurrent());
        CHECK(m.maxDivergence() < 1e-5f, "div B");
    }

    // 2) The m = 1, n = 1 kink of the constant-current column with a conducting wall at b = 2a
    //    (Tokamak.h): unstable for 2a^2/(a^2+b^2) = 0.4 < q_a < 1, held by the line tension above
    //    (Kruskal-Shafranov) and by the wall below. The amplitude is the n = 1 harmonic of the
    //    current centroid's displacement around the torus; it starts at the 2 % seed. Its growth
    //    rate is compared with gamma of the energy principle. (Higher m are other modes with
    //    their own windows; the centroid does not see them.)
    sim.grid.params.resolutionX = 44; // dx = 3.2 cm: enough for the m = 1 mode, 2.5x faster
    struct Run { float qa; bool unstable; float seconds; };
    const Run runs[3] = {{0.7f, true, 1.6f}, {1.5f, false, 1.0f}, {0.25f, false, 1.0f}};
    for (const Run& r : runs) {
        sim.tokamak.safetyFactorEdge = r.qa;
        sim.reset();
        const Tokamak& t = sim.tokamak;
        const MagneticField& m = sim.grid.magnetic;
        const float dx = sim.grid.dx(), a = t.minorRadius, gap = t.vesselRadius - a;
        std::vector<std::pair<float, float>> series = {{0.0f, t.kinkAmplitude(m, dx)}};
        while (sim.time() < r.seconds) {
            sim.stepFrame();
            if (sim.time() >= series.back().first + 0.2f - 1e-4f) series.push_back({sim.time(), t.kinkAmplitude(m, dx)});
        }
        const float seed = series.front().second, last = series.back().second;
        // Growth rate: the log slope over the samples between 1.5 x seed and half way to the wall.
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        int n = 0;
        for (auto [time, amp] : series)
            if (amp > 1.5f * seed && amp < 0.5f * gap) { sx += time; sy += std::log(amp); sxx += time * time; sxy += time * std::log(amp); ++n; }
        const float gamma = n >= 2 ? float((n * sxy - sx * sy) / (n * sxx - sx * sx)) : 0.0f;
        std::printf("  q_a = %.2f (%s): kink amplitude", r.qa, r.unstable ? "unstable" : "stable");
        for (auto [time, amp] : series) std::printf(" %.1f", amp * 1000);
        std::printf(" mm at 0.2 s steps (a = %.0f mm); growth rate %.2f 1/s (theory %.2f), I_p %.0f of %.0f A, ring shift %.1f mm, "
                    "max speed %.2f m/s, div B %.1e\n",
                    a * 1000, gamma, t.kinkGrowthRate(sim.grid.params.fluidDensity), t.measuredCurrent(m, dx), t.plasmaCurrent(),
                    t.measuredShift(m, dx) * 1000, sim.grid.maxVelocity(), m.maxDivergence());
        CHECK(std::isfinite(last) && m.maxDivergence() < 1e-4f, "tokamak field broken at q_a = %f", r.qa);
        CHECK(std::fabs(t.measuredShift(m, dx)) < 0.1f * (t.vesselRadius - a), "position control lost the ring: shift %f m",
              t.measuredShift(m, dx));
        if (r.unstable) {
            CHECK(last > 0.2f * a && last > 5.0f * seed, "q_a = %.2f must kink: %f -> %f m", r.qa, seed, last);
            const float g0 = t.kinkGrowthRate(sim.grid.params.fluidDensity);
            CHECK(gamma > 0.5f * g0 && gamma < 2.0f * g0, "kink growth rate %f vs theory %f", gamma, g0);
        } else {
            CHECK(last < 2.5f * seed && last < 0.06f * a, "q_a = %.2f must keep the m = 1 shape: %f -> %f m", r.qa, seed, last);
        }
    }
}
