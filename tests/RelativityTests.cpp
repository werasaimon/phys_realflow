// Light and free particles in the Kerr and Schwarzschild spacetimes against the exact results of
// general relativity: the constants of motion along a geodesic (Noether: E, L; Carter: Q), the
// deflection of light, the photon sphere and the ISCO, the perihelion precession, the shadow of
// the hole, and the horizon crossing in horizon-penetrating coordinates.
//
// Registration (tests/Tests.h and tests/main.cpp):
//   run("geodesics: E, L and Carter's Q along a Kerr orbit (RK45 vs RK4)", testGeodesicInvariants);
//   run("light deflection by a mass: 4M/b + second order", testLightDeflection);
//   run("photon sphere at 3M and the ISCO at 6M", testPhotonSphere);
//   run("perihelion precession: 6 pi M / (a (1 - e^2))", testPerihelionPrecession);
//   run("shadow of a Schwarzschild hole: 3 sqrt(3) M (ray tracer)", testShadow);
//   run("horizon crossing: Eddington-Finkelstein vs Boyer-Lindquist", testHorizonPenetration);
#include "TestRunner.h"

#include "relativity/Geodesic.h"
#include "relativity/Metric.h"
#include "relativity/RayTracer.h"

#include <cmath>

using namespace rf;

static const double kPiD = 3.14159265358979323846;

void testGeodesicInvariants() {
    Metric m;
    m.M = 1.0;
    m.a = 0.9;
    Geodesic geo(m);
    // 1) A unit-mass particle on an inclined, eccentric bound orbit: r = 10 M, leaving the
    //    equatorial plane with a local speed of 0.34 c (the circular speed there is 0.35 c).
    auto boundParticle = [&] {
        GeodesicState s;
        s.x[1] = 10.0;
        s.x[2] = 0.5 * kPiD;
        double dir[3] = {0.0, 0.3, std::sqrt(1.0 - 0.09)};
        Geodesic::momentumFromLocal(m, s.x, dir, 0.34, 1.0, s.p);
        return s;
    };
    Geodesic::Options opt;
    opt.mu = 1;
    opt.tol = 1e-12;
    opt.lambdaMax = 20000.0; // ~ 100 orbits of the particle
    opt.hMax = 2.0;
    GeodesicState s = boundParticle();
    GeodesicResult r = geo.integrate(s, opt);
    // E and L are exact to rounding whatever the integrator: the metric does not depend on t and
    // phi, so their derivatives are exactly zero and p_t, p_phi are never touched. Q is the real
    // measure: it needs the r and theta motion to be right.
    const double dE = std::fabs(r.E1 / r.E0 - 1), dL = std::fabs(r.L1 / r.L0 - 1), dQ = std::fabs(r.Q1 - r.Q0) / std::max(std::fabs(r.Q0), 1e-30);
    std::printf("  massive, Kerr a = 0.9, RK45 tol 1e-12, lambda %.0f M (%d steps, r in the end %.2f M): dE %.1e, dL %.1e, dQ %.1e, H %.2e -> %.2e\n",
                s.lambda, r.steps, s.x[1], dE, dL, dQ, r.H0, r.H1);
    CHECK(r.end == GeodesicEnd::LambdaLimit, "the bound particle left the orbit (end %d)", int(r.end));
    CHECK(dE < 1e-8 && dL < 1e-8 && dQ < 1e-8, "constants of motion drift: E %e L %e Q %e", dE, dL, dQ);
    CHECK(std::fabs(r.H1 + 0.5) < 1e-8, "H of a unit-mass particle is %.10f, not -1/2", r.H1);
    // The same orbit with fixed-step RK4: the drift grows with the step - the reason a symplectic
    // integrator is the next step for long integrations.
    for (double h : {2.0, 0.5}) {
        Geodesic::Options rk4 = opt;
        rk4.adaptive = false;
        rk4.hFixed = h;
        GeodesicState s4 = boundParticle();
        GeodesicResult r4 = geo.integrate(s4, rk4);
        std::printf("  same orbit, RK4 h = %.1f: dE %.1e, dL %.1e, dQ %.1e\n", h, std::fabs(r4.E1 / r4.E0 - 1),
                    std::fabs(r4.L1 / r4.L0 - 1), std::fabs(r4.Q1 - r4.Q0) / std::max(std::fabs(r4.Q0), 1e-30));
    }
    // 2) A photon skimming the hole: from r = 10 M inwards and sideways, out of the plane.
    GeodesicState ph;
    ph.x[1] = 10.0;
    ph.x[2] = 0.4 * kPiD;
    double dirPh[3] = {-0.55, 0.2, std::sqrt(1.0 - 0.55 * 0.55 - 0.04)};
    Geodesic::momentumFromLocal(m, ph.x, dirPh, 1.0, 1.0, ph.p);
    Geodesic::Options optPh = opt;
    optPh.mu = 0;
    optPh.rMax = 1e4;
    optPh.hMax = 50.0;
    GeodesicResult rp = geo.integrate(ph, optPh);
    const double dEp = std::fabs(rp.E1 / rp.E0 - 1), dLp = std::fabs(rp.L1 / rp.L0 - 1), dQp = std::fabs(rp.Q1 - rp.Q0) / std::max(std::fabs(rp.Q0), 1e-30);
    std::printf("  photon, Kerr a = 0.9: %s after %d steps (r %.1f M, phi %.2f turns): dE %.1e, dL %.1e, dQ %.1e, H %.1e\n",
                rp.end == GeodesicEnd::Escaped ? "escaped" : rp.end == GeodesicEnd::Captured ? "captured" : "stopped", rp.steps, ph.x[1],
                ph.x[3] / (2 * kPiD), dEp, dLp, dQp, rp.H1);
    CHECK(dEp < 1e-8 && dLp < 1e-8 && dQp < 1e-8, "photon constants drift: E %e L %e Q %e", dEp, dLp, dQp);
    CHECK(std::fabs(rp.H1) < 1e-8, "H of a photon is %e, not 0", rp.H1);
}

void testLightDeflection() {
    // A photon from far away with impact parameter b = 200 M past a Schwarzschild mass: the
    // deflection 4 M / b + 15 pi M^2 / (4 b^2) (Einstein 1915 and the second order, e.g. Keeton &
    // Petters 2005). Started at r0 = 1e6 M on the straight line of flat space; the deflection
    // still to come beyond r0 (~4 M / r0 = 4e-6) is below the tolerance.
    Metric m;
    Geodesic geo(m);
    const double b = 200.0, r0 = 1e6;
    GeodesicState s;
    const double x0 = -std::sqrt(r0 * r0 - b * b), y0 = b;
    s.x[1] = r0;
    s.x[2] = 0.5 * kPiD;
    s.x[3] = std::atan2(y0, x0);
    // Velocity (1, 0) in the plane: radial part x0 / r0, azimuthal part -b / r0 (clockwise).
    double dir[3] = {x0 / r0, 0.0, -b / r0};
    Geodesic::momentumFromLocal(m, s.x, dir, 1.0, 1.0, s.p);
    Geodesic::Options opt;
    opt.tol = 1e-11;
    opt.rMax = r0 * (1.0 + 1e-9);
    opt.lambdaMax = 1e7;
    opt.hMax = 2000.0;
    GeodesicResult r = geo.integrate(s, opt);
    // Outgoing direction: the local velocity in the static frame back to Cartesian.
    double g[4][4], ginv[4][4];
    m.covariant(s.x, g);
    m.contravariant(s.x, ginv);
    double pu[4] = {0, 0, 0, 0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) pu[i] += ginv[i][j] * s.p[j];
    const double vr = std::sqrt(g[1][1]) * pu[1], vphi = std::sqrt(g[3][3]) * pu[3];
    const double norm = std::sqrt(vr * vr + vphi * vphi);
    const double cx = std::cos(s.x[3]), sx = std::sin(s.x[3]);
    const double vx = (vr * cx - vphi * sx) / norm, vy = (vr * sx + vphi * cx) / norm;
    const double deflection = std::fabs(std::atan2(vy, vx)); // from the incoming direction (1, 0), towards the mass (-y)
    const double theory = 4.0 / b + 15.0 * kPiD / (4.0 * b * b);
    std::printf("  deflection at b = %.0f M: %.7f rad, theory 4M/b + 15 pi M^2 / 4b^2 = %.7f (first order alone %.7f); %d steps, r_min %s\n",
                b, deflection, theory, 4.0 / b, r.steps, r.end == GeodesicEnd::Escaped ? "reached" : "not reached");
    CHECK(r.end == GeodesicEnd::Escaped, "the photon did not come back out");
    CHECK(std::fabs(deflection / theory - 1.0) < 1e-3, "deflection %.7f vs %.7f", deflection, theory);
}

void testPhotonSphere() {
    Metric m;
    Geodesic geo(m);
    // 1) A photon sent tangentially at r = 3 M circles the hole: the orbit is unstable (Lyapunov
    //    exponent 1 / (3 sqrt 3 M) in t: every revolution multiplies a perturbation by e^{2 pi}),
    //    so rounding lets it go after a few revolutions - a light ray from the back of your head
    //    reaches your eye along it.
    {
        GeodesicState s;
        s.x[1] = 3.0;
        s.x[2] = 0.5 * kPiD;
        double dir[3] = {0.0, 0.0, 1.0};
        Geodesic::momentumFromLocal(m, s.x, dir, 1.0, 1.0, s.p);
        double h = 0.01;
        double turns = 0;
        int steps = 0;
        while (std::fabs(s.x[1] / 3.0 - 1.0) < 0.01 && steps < 200000) {
            geo.stepRK45(s, h, 1e-13, 0.05);
            turns = s.x[3] / (2 * kPiD);
            ++steps;
        }
        std::printf("  photon sphere: tangential start at r = 3 M stays within 1%% for %.2f revolutions (b = L/E = %.6f, 3 sqrt 3 = %.6f)\n",
                    turns, s.p[3] / -s.p[0], 3.0 * std::sqrt(3.0));
        CHECK(turns >= 3.0, "the photon left the photon sphere after %.2f revolutions", turns);
    }
    // 2) The ISCO: the circular orbit at r = 6 M with its exact E and L stays circular; at
    //    5.9 M (inside the ISCO, an unstable circular orbit) the smallest inward push ends in the
    //    hole within a few orbits.
    {
        double E, L;
        m.circularOrbit(6.0, E, L);
        GeodesicState s;
        s.x[1] = 6.0;
        s.x[2] = 0.5 * kPiD;
        s.p[0] = -E;
        s.p[3] = L;
        Geodesic::Options opt;
        opt.mu = 1;
        opt.tol = 1e-11;
        opt.hMax = 1.0;
        double worst = 0;
        opt.lambdaMax = 10.0 * 2.0 * kPiD / m.keplerianOmega(6.0) / 1.2247; // 10 orbits of proper time (u^t = 1.2247 at 6 M)
        GeodesicResult r = geo.integrate(s, opt, [&](const GeodesicState&, const GeodesicState& after) {
            worst = std::max(worst, std::fabs(after.x[1] / 6.0 - 1.0));
            return false;
        });
        std::printf("  ISCO: circular orbit at 6 M (E %.6f, L %.6f) over %.1f turns: |r/6M - 1| <= %.1e\n", E, L, s.x[3] / (2 * kPiD), worst);
        CHECK(r.end == GeodesicEnd::LambdaLimit && worst < 1e-6, "the ISCO orbit drifted by %e", worst);
        m.circularOrbit(5.9, E, L);
        GeodesicState u;
        u.x[1] = 5.9;
        u.x[2] = 0.5 * kPiD;
        u.p[0] = -E;
        u.p[3] = L;
        u.p[1] = -1e-3; // an inward nudge
        Geodesic::Options opt2 = opt;
        opt2.lambdaMax = 5000.0;
        GeodesicResult r2 = geo.integrate(u, opt2);
        std::printf("  inside the ISCO: circular orbit at 5.9 M with an inward nudge %s after %.1f turns (lambda %.0f M)\n",
                    r2.end == GeodesicEnd::Captured ? "fell in" : "did not fall in", u.x[3] / (2 * kPiD), u.lambda);
        CHECK(r2.end == GeodesicEnd::Captured, "the orbit inside the ISCO did not fall in");
    }
}

void testPerihelionPrecession() {
    // A unit-mass particle on an orbit with perihelion 80 M and aphelion 120 M (a = 100 M,
    // e = 0.2) about a Schwarzschild mass: the perihelion advances 6 pi M / (a (1 - e^2)) per
    // orbit (Einstein 1915; for Mercury 43" per century).
    Metric m;
    Geodesic geo(m);
    const double rp = 80.0, ra = 120.0;
    // E^2 = (1 - 2M/r)(1 + L^2/r^2) at both turning points.
    const double fp = 1.0 - 2.0 / rp, fa = 1.0 - 2.0 / ra;
    const double L2 = (fp - fa) / (fa / (ra * ra) - fp / (rp * rp));
    const double E = std::sqrt(fp * (1.0 + L2 / (rp * rp)));
    GeodesicState s;
    s.x[1] = rp;
    s.x[2] = 0.5 * kPiD;
    s.p[0] = -E;
    s.p[3] = std::sqrt(L2);
    Geodesic::Options opt;
    opt.mu = 1;
    opt.tol = 1e-11;
    opt.hMax = 20.0;
    opt.lambdaMax = 1e7;
    std::vector<double> perihelia; // phi at successive minima of r
    double prevDr = 0;
    GeodesicState prev; // the state before `before`
    int orbits = 0;
    GeodesicResult r = geo.integrate(s, opt, [&](const GeodesicState& before, const GeodesicState& after) {
        const double dr = after.x[1] - before.x[1];
        if (prevDr < 0 && dr > 0 && before.lambda > 0) {
            // r turned upwards around `before`: the vertex of the parabola through the three
            // states gives the perihelion's lambda, and phi is interpolated there (a step is up
            // to 0.2 rad of phi here; the raw phi of `before` would jitter by that much).
            const double l0 = prev.lambda, l1 = before.lambda, l2 = after.lambda;
            const double r0 = prev.x[1], r1 = before.x[1], r2 = after.x[1];
            const double d1 = (r1 - r0) / (l1 - l0), d2 = (r2 - r1) / (l2 - l1);
            const double curvature = (d2 - d1) / (0.5 * (l2 - l0));
            const double lStar = 0.5 * (l0 + l1) - d1 / curvature;
            const double phi = lStar < l1 ? prev.x[3] + (before.x[3] - prev.x[3]) * (lStar - l0) / (l1 - l0)
                                          : before.x[3] + (after.x[3] - before.x[3]) * (lStar - l1) / (l2 - l1);
            perihelia.push_back(phi);
            ++orbits;
        }
        prevDr = dr;
        prev = before;
        return orbits >= 12;
    });
    (void)r;
    double sum = 0;
    for (size_t i = 1; i < perihelia.size(); ++i) sum += perihelia[i] - perihelia[i - 1] - 2.0 * kPiD;
    const double perOrbit = sum / double(perihelia.size() - 1);
    // Einstein's 6 pi M / p (p = a (1 - e^2)) is the first order in M / p; at p = 96 M the next
    // order, (3 pi / 2)(18 + e^2)(M / p)^2, is already 4.7 %. The exact value is the quadrature
    //   Delta phi = 2 int_{r_p}^{r_a} L dr / (r^2 sqrt(E^2 - (1 - 2M/r)(1 + L^2/r^2))) - 2 pi
    // (the radial equation of the same E and L), taken with r = c - d cos(chi) so that the
    // square-root endpoints are regular (Simpson, 40 000 intervals).
    const double first = 6.0 * kPiD / (100.0 * (1.0 - 0.04));
    const double L = std::sqrt(L2), c = 0.5 * (rp + ra), d = 0.5 * (ra - rp);
    auto f = [&](double chi) {
        const double r = c - d * std::cos(chi);
        const double V = E * E - (1.0 - 2.0 / r) * (1.0 + L2 / (r * r));
        const double drdchi = d * std::sin(chi); // vanishes like chi at the ends, as sqrt(V) does: the ratio is finite
        return L / (r * r) * drdchi / std::sqrt(std::max(V, 1e-300));
    };
    const int n = 40000;
    double integral = 0;
    for (int i = 1; i < n; ++i) { // the endpoints (0/0 in floating point) are left out: an O(1/n) part of a finite value
        const double chi = kPiD * i / n;
        integral += (i % 2 ? 4.0 : 2.0) * f(chi);
    }
    integral *= kPiD / n / 3.0;
    const double exact = 2.0 * integral - 2.0 * kPiD;
    std::printf("  perihelion precession, a = 100 M, e = 0.2: %.6f rad per orbit over %zu orbits; exact quadrature %.6f (%.3f%%), Einstein's first order 6 pi M / p = %.6f (%.1f%% below the exact)\n",
                perOrbit, perihelia.size() - 1, exact, 100.0 * (perOrbit / exact - 1.0), first, 100.0 * (1.0 - first / exact));
    CHECK(perihelia.size() >= 10, "only %zu perihelia found", perihelia.size());
    CHECK(std::fabs(perOrbit / exact - 1.0) < 2e-3, "precession %.6f vs exact %.6f", perOrbit, exact);
    CHECK(std::fabs(perOrbit / first - 1.0) < 0.1, "precession %.6f vs first order %.6f", perOrbit, first);
}

void testShadow() {
    // The shadow of a non-rotating hole seen from r0 = 1000 M: a black disk of angular radius
    // 3 sqrt(3) M / r0 (Bardeen 1973) - every ray inside it spirals into the hole. Measured from
    // the ray-traced image by counting black pixels.
    Metric m;
    RayTracer rt(m);
    const double r0 = 1000.0, thetaShadow = Metric::shadowRadiusSchwarzschild(1.0) / r0;
    rt.camera.r = r0;
    rt.camera.theta = 0.5 * kPiD;
    rt.camera.fovDeg = 4.0 * thetaShadow * 180.0 / kPiD;
    rt.camera.width = rt.camera.height = 96;
    rt.options.tol = 1e-6; // the edge is set by the pixel size, not by the integrator
    rt.options.hMax = 500.0; // the step control shrinks it near the hole; far out the field is flat
    rt.options.lambdaMax = 1e5;
    auto t0 = std::chrono::steady_clock::now();
    Image img = rt.render();
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    int black = 0, lost = 0;
    for (int y = 0; y < img.height; ++y)
        for (int x = 0; x < img.width; ++x) {
            const float* p = img.pixel(x, y);
            if (p[0] == 0 && p[1] == 0 && p[2] == 0) ++black;
            if (p[0] == 0.2f && p[1] == 0.2f && p[2] == 0.2f) ++lost;
        }
    const double fov = rt.camera.fovDeg * kPiD / 180.0;
    const double side = 2.0 * std::tan(0.5 * fov); // the image plane in units of the distance
    const double fraction = double(black) / double(img.width * img.height);
    const double radius = std::atan(std::sqrt(fraction * side * side / kPiD));
    std::printf("  shadow: %d x %d rays in %.0f ms, %d black pixels (%d lost): angular radius %.5f rad, Bardeen 3 sqrt(3) M / r0 = %.5f (%.2f%%)\n",
                img.width, img.height, ms, black, lost, radius, thetaShadow, 100.0 * (radius / thetaShadow - 1.0));
    if (std::getenv("RF_WRITE_IMAGES")) RayTracer::writePpm(img, "shadow.ppm");
    CHECK(lost == 0, "%d rays ended in neither the sky nor the hole", lost);
    CHECK(std::fabs(radius / thetaShadow - 1.0) < 0.03, "shadow radius %.5f vs %.5f", radius, thetaShadow);
}

void testHorizonPenetration() {
    // A radially infalling photon. In ingoing Eddington-Finkelstein coordinates the horizon is an
    // ordinary place: v stays constant and r runs through 2 M to the centre. In Boyer-Lindquist
    // (Schwarzschild) coordinates the same ray approaches r = 2 M with t growing without bound:
    // the coordinates, not the spacetime, break there - a GRMHD grid has to use the former.
    Metric ef;
    ef.coordinates = Metric::Coordinates::EddingtonFinkelstein;
    Geodesic geoEf(ef);
    GeodesicState s;
    s.x[1] = 10.0;
    s.x[2] = 0.5 * kPiD;
    s.p[0] = -1.0; // p_v = -E; p_r = 0 is the ingoing radial null direction (p^r = p_v = -E)
    Geodesic::Options opt;
    opt.tol = 1e-11;
    opt.rMin = 0.15;
    opt.rMax = 100.0;
    opt.hMax = 0.1;
    GeodesicResult r = geoEf.integrate(s, opt);
    std::printf("  Eddington-Finkelstein: r 10 M -> %.3f M in %d steps, v = %.2e (constant), lambda %.3f M, H %.1e\n", s.x[1], r.steps, s.x[0],
                s.lambda, r.H1);
    CHECK(r.end == GeodesicEnd::Captured && s.x[1] < 0.16, "the photon did not reach r = 0.15 M (r %.3f, end %d)", s.x[1], int(r.end));
    CHECK(std::isfinite(s.x[0]) && std::fabs(s.x[0]) < 1e-6, "v changed along the ingoing ray: %e", s.x[0]);

    Metric bl;
    Geodesic geoBl(bl);
    GeodesicState b;
    b.x[1] = 10.0;
    b.x[2] = 0.5 * kPiD;
    b.p[0] = -1.0;
    b.p[1] = -1.0 / (1.0 - 2.0 / 10.0); // ingoing: p^r = (1 - 2M/r) p_r < 0
    Geodesic::Options optBl = opt;
    optBl.rMin = 2.0 * (1.0 + 1e-6);
    optBl.lambdaMax = 1e3;
    // p_r = -E / (1 - 2M/r) diverges at the horizon: the step control stalls in ever smaller
    // steps, and once rounding wins the ray may even bounce - record the closest approach.
    double rClosest = 10.0, tClosest = 0.0;
    GeodesicResult rb = geoBl.integrate(b, optBl, [&](const GeodesicState&, const GeodesicState& after) {
        if (after.x[1] < rClosest) {
            rClosest = after.x[1];
            tClosest = after.x[0];
        }
        return false;
    });
    std::printf("  Boyer-Lindquist: r 10 M -> %.9f M closest in %d steps, t there %.1f M and growing like 2 M ln(1 / (r/2M - 1)) = %.1f; end %d\n",
                rClosest, rb.steps, tClosest, -2.0 * std::log(rClosest / 2.0 - 1.0), int(rb.end));
    CHECK(rClosest < 2.0 * (1.0 + 1e-4), "the Boyer-Lindquist ray stopped %e M above the horizon", rClosest - 2.0);
    CHECK(tClosest > 20.0, "t at the horizon should be large: %.1f", tClosest);
}
