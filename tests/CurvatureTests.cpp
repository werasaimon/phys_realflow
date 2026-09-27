// Curvature from the metric (relativity/Curvature.h) against the textbook answers: coordinates are
// not curvature (flat space in spherical coordinates), the 2-sphere, Schwarzschild's Christoffel
// symbols and Kretschmann scalar, Kerr (Henry 2000), Reissner-Nordstrom, de Sitter, Friedmann, the
// wormhole's negative energy, the tides of a static observer, and the geodesic equation with Gamma
// against Hamilton's equations of Geodesic.h.
//
// Tolerances come from the error budget in Curvature.h: with the default step the curvature is
// right to ~2e-10 relative (truncation and rounding meet there), Gamma to ~2e-11; checked at
// 1e-8 .. 1e-7 (a margin of 50 or more against another point, metric or machine).
//
// Registration (tests/Tests.h and tests/main.cpp):
//   run("curvature: flat space in spherical coordinates, the 2-sphere, loops vs einstein()", testCurvatureFlatAndSphere);
//   run("curvature: Schwarzschild, Kerr, Reissner-Nordstrom (Christoffels, Ricci, Kretschmann)", testCurvatureBlackHoles);
//   run("curvature: de Sitter, Friedmann, a wormhole's exotic matter", testCurvatureCosmology);
//   run("curvature: tides of a static observer, geodesics with Gamma vs Hamilton", testTidesAndGammaGeodesic);
#include "TestRunner.h"

#include "relativity/Curvature.h"
#include "relativity/Geodesic.h"
#include "relativity/Spacetime.h"

#include <cmath>

using namespace rf;

namespace {
const double kPiD = 3.14159265358979323846;

// Schwarzschild with a chosen finite-difference step (a fraction of the coordinate's size), to
// measure how the error of Gamma and of the curvature depends on it.
struct SchwarzschildWithStep : SchwarzschildSpacetime {
    SchwarzschildWithStep(double mass, double fraction) : SchwarzschildSpacetime(mass), fraction(fraction) {}
    double fraction;
    double differenceStep(int k, const double* x) const override { return fraction * std::max(1.0, std::fabs(x[k])); }
};

// The error of the finite differences against the step, at r = 7M: the largest relative error of
// four Christoffel symbols and of the Kretschmann scalar. Writes <RF_PLOT_DIR>/curvature_step.csv
// for tools/plot_curvature.py when RF_PLOT_DIR is set.
void measureStepConvergence(double& errAtDefault) {
    const double M = 1.0, r = 7.0, x[4] = {0.0, r, 1.1, 0.4}, f = r - 2 * M;
    const char* dir = std::getenv("RF_PLOT_DIR");
    FILE* csv = dir && *dir ? std::fopen((std::string(dir) + "/curvature_step.csv").c_str(), "w") : nullptr;
    if (csv) std::fprintf(csv, "step_fraction,gamma_rel_err,kretschmann_rel_err\n");
    std::printf("  step h / (coordinate size) -> relative error of Gamma, of K:");
    for (double h : {3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5}) {
        SchwarzschildWithStep st(M, h);
        const Tensor G = christoffel(st, x);
        const double eG = std::max({std::fabs(G(1, 0, 0) / (M * f / (r * r * r)) - 1), std::fabs(G(0, 0, 1) / (M / (r * f)) - 1),
                                    std::fabs(G(1, 1, 1) / (-M / (r * f)) - 1), std::fabs(G(2, 1, 2) * r - 1)});
        const double eK = std::fabs(curvatureAt(st, x).K / (48 * M * M / std::pow(r, 6)) - 1);
        if (h == 1e-3) errAtDefault = eK; // the default step of Spacetime
        std::printf(" %.0e: %.1e, %.1e;", h, eG, eK);
        if (csv) std::fprintf(csv, "%.9g,%.9g,%.9g\n", h, eG, eK);
    }
    std::printf("\n");
    if (csv) std::fclose(csv);
}
}

// Coordinates are not curvature: flat space in spherical coordinates has Christoffel symbols
// (Gamma^r_thth = -r) but no Riemann tensor. A sphere's surface is curved: R = 2 / radius^2. And
// the same Schwarzschild point by loops and by einstein(): the two ways must agree.
void testCurvatureFlatAndSphere() {
    MinkowskiSpherical flat;
    const double x[4] = {0.0, 3.0, 1.0, 0.5};
    const Tensor G = christoffel(flat, x);
    const Tensor R = riemann(flat, x);
    MinkowskiCartesian cart;
    const double xc[4] = {0.0, 1.0, 2.0, 3.0};
    const double cartR = riemann(cart, xc).maxAbs();
    TwoSphere sphere(2.0);
    const double xs[2] = {1.0, 0.3};
    const CurvatureAt cs = curvatureAt(sphere, xs);
    SchwarzschildSpacetime bh(1.0);
    const double xb[4] = {0.0, 7.0, 1.1, 0.4};
    const Tensor gL = christoffel(bh, xb), gE = christoffelEinstein(bh, xb);
    const Tensor rL = riemann(bh, xb), rE = riemannEinstein(bh, xb);
    const Tensor g = metricAt(bh, xb), gi = inverseMetricAt(bh, xb);
    const double kL = kretschmann(rL, g, gi), kE = kretschmannEinstein(rL, g, gi);
    const double ricDiff = (ricci(rL) - ricciEinstein(rL)).maxAbs();
    std::printf("  flat spherical: Gamma^r_thth %.10f (-r = -3), max |Riemann| %.1e; Cartesian %.1e; 2-sphere R = %.10f (2/R^2 = 0.5)\n",
                G(1, 2, 2), R.maxAbs(), cartR, cs.R);
    std::printf("  loops vs einstein(): Gamma %.1e, Riemann %.1e, Ricci %.1e, K %.1e (K = %.6e)\n", (gL - gE).maxAbs(),
                (rL - rE).maxAbs(), ricDiff, std::fabs(kL / kE - 1), kL);
    CHECK(std::fabs(G(1, 2, 2) + 3.0) < 1e-10, "Gamma^r_thth = %f, not -r", G(1, 2, 2));
    CHECK(R.maxAbs() < 1e-9 && cartR == 0.0, "flat space is curved: %e", R.maxAbs());
    CHECK(std::fabs(cs.R - 0.5) < 1e-8, "2-sphere scalar curvature %f, not 2/R^2", cs.R);
    CHECK((gL - gE).maxAbs() < 1e-13 && (rL - rE).maxAbs() < 1e-10 && ricDiff < 1e-15 && std::fabs(kL / kE - 1) < 1e-12,
          "the loops and einstein() disagree");
}

// Black holes: Schwarzschild's Christoffels against their closed forms, vacuum (Ricci = 0) and the
// Kretschmann scalar 48 M^2 / r^6; Kerr's Kretschmann (Henry 2000, ApJ 535, 350); Reissner-
// Nordstrom: R = 0 (the electromagnetic stress tensor is traceless) but Ricci != 0, G^t_t = -Q^2/r^4
// (8 pi times the field's energy density), K = 8 (6 M^2 r^2 - 12 M Q^2 r + 7 Q^4) / r^8.
void testCurvatureBlackHoles() {
    const double M = 1.0, r = 7.0, th = 1.1;
    const double x[4] = {0.0, r, th, 0.4};
    SchwarzschildSpacetime sw(M);
    const CurvatureAt c = curvatureAt(sw, x);
    const double f = r - 2 * M, s = std::sin(th), co = std::cos(th);
    struct Known { int l, m, n; double value; };
    const Known known[] = {{1, 0, 0, M * f / (r * r * r)}, {0, 0, 1, M / (r * f)}, {1, 1, 1, -M / (r * f)},
                           {2, 1, 2, 1 / r},               {1, 2, 2, -f},          {3, 2, 3, co / s},
                           {2, 3, 3, -s * co},             {1, 3, 3, -f * s * s},  {3, 1, 3, 1 / r}};
    double worstGamma = 0;
    for (const Known& k : known) worstGamma = std::max(worstGamma, std::fabs(c.gamma(k.l, k.m, k.n) - k.value) / std::fabs(k.value));
    const double K0 = 48 * M * M / std::pow(r, 6), scale = M / (r * r * r); // the curvature scale
    std::printf("  Schwarzschild r = 7M: 9 Christoffels vs closed form %.1e (relative); max |Ricci| %.1e; K %.10e vs 48M^2/r^6 %.10e\n",
                worstGamma, c.ricci.maxAbs(), c.K, K0);
    CHECK(worstGamma < 1e-8, "Christoffel symbols off by %e", worstGamma); // ~2e-11 at the default step (Curvature.h)
    double errDefault = 1;
    measureStepConvergence(errDefault);
    CHECK(errDefault < 1e-9, "the default step is not where the error is small: %e", errDefault);
    CHECK(c.ricci.maxAbs() < 1e-8 * scale * r * r && std::fabs(c.K / K0 - 1) < 1e-7, "Schwarzschild curvature wrong");
    Metric km;
    km.M = M;
    km.a = 0.9;
    KerrSpacetime kerr(km);
    const CurvatureAt ck = curvatureAt(kerr, x);
    const double a = km.a, ac = a * co, r2 = r * r;
    const double Kh = 48 * M * M * (std::pow(r, 6) - 15 * r2 * r2 * ac * ac + 15 * r2 * std::pow(ac, 4) - std::pow(ac, 6)) / std::pow(r2 + ac * ac, 6);
    std::printf("  Kerr a = 0.9, r = 7M, theta 1.1: K %.10e vs Henry 2000 %.10e (%.1e); max |Ricci| %.1e\n", ck.K, Kh,
                std::fabs(ck.K / Kh - 1), ck.ricci.maxAbs());
    CHECK(std::fabs(ck.K / Kh - 1) < 1e-7 && ck.ricci.maxAbs() < 1e-8 * scale * r * r, "Kerr curvature wrong");
    ReissnerNordstrom rn(M, 0.6);
    const CurvatureAt cr = curvatureAt(rn, x);
    const double Q = rn.Q, Krn = 8 * (6 * M * M * r2 - 12 * M * Q * Q * r + 7 * Q * Q * Q * Q) / std::pow(r, 8);
    const double Gtt = cr.ginv(0, 0) * cr.einstein(0, 0); // G^t_t (the metric is diagonal)
    std::printf("  Reissner-Nordstrom Q = 0.6: R %.1e, max |Ricci| %.3e (Q^2/r^4 = %.3e), G^t_t %.8e vs -Q^2/r^4 %.8e, K %.8e vs %.8e\n",
                cr.R, cr.ricci.maxAbs(), Q * Q / (r2 * r2), Gtt, -Q * Q / (r2 * r2), cr.K, Krn);
    CHECK(std::fabs(cr.R) < 1e-10 && cr.ricci.maxAbs() > 0.5 * Q * Q / (r2 * r2), "Reissner-Nordstrom: R must vanish, Ricci not");
    CHECK(std::fabs(Gtt * r2 * r2 / (Q * Q) + 1) < 1e-7 && std::fabs(cr.K / Krn - 1) < 1e-7, "Reissner-Nordstrom curvature wrong");
}

// Cosmology and exotica: de Sitter solves G_mn + Lambda g_mn = 0 with R = 4 Lambda; a flat universe
// a(t) = t^(2/3) expands at H = a'/a = 2/(3t) and obeys Friedmann's G_tt = 3 H^2 = 8 pi rho; the
// Ellis wormhole's throat needs negative energy density rho = G_tt / 8 pi = -1 / (8 pi b0^2)
// (Morris & Thorne 1988: no ordinary matter can hold a traversable wormhole open).
void testCurvatureCosmology() {
    const double Lambda = 0.03;
    DeSitterStatic ds(Lambda);
    const double x[4] = {0.0, 4.0, 1.0, 0.2};
    const CurvatureAt cd = curvatureAt(ds, x);
    const double vacuum = (cd.einstein + cd.g * Lambda).maxAbs();
    FlatFLRW flrw([](double t) { return std::pow(t, 2.0 / 3.0); });
    const double t = 2.0, xf[4] = {t, 0.3, 0.4, 0.5};
    const CurvatureAt cf = curvatureAt(flrw, xf);
    const double H = 2.0 / (3.0 * t);
    EllisWormhole wh(1.0);
    const double xw[4] = {0.0, 0.0, 1.2, 0.0}; // the throat l = 0
    const CurvatureAt cw = curvatureAt(wh, xw);
    const double rho = cw.einstein(0, 0) / (8 * kPiD), rhoExact = -1.0 / (8 * kPiD * wh.b0 * wh.b0);
    std::printf("  de Sitter Lambda = 0.03: R %.10f (4 Lambda = 0.12), max |G + Lambda g| %.1e\n", cd.R, vacuum);
    std::printf("  flat FLRW a = t^(2/3), t = 2: Gamma^x_tx %.10f vs H %.10f; G_tt %.10f vs 3H^2 %.10f\n", cf.gamma(1, 0, 1), H,
                cf.einstein(0, 0), 3 * H * H);
    std::printf("  Ellis wormhole b0 = 1, at the throat: rho = G_tt / 8pi %.10f vs -1/(8 pi b0^2) %.10f (negative: exotic matter)\n",
                rho, rhoExact);
    CHECK(std::fabs(cd.R / (4 * Lambda) - 1) < 1e-8 && vacuum < 1e-9, "de Sitter wrong");
    CHECK(std::fabs(cf.gamma(1, 0, 1) / H - 1) < 1e-10 && std::fabs(cf.einstein(0, 0) / (3 * H * H) - 1) < 1e-8, "Friedmann equation wrong");
    CHECK(rho < 0 && std::fabs(rho / rhoExact - 1) < 1e-8, "wormhole energy density %f vs %f", rho, rhoExact);
}

// Tides: two free fallers next to a static observer at r in Schwarzschild drift apart radially at
// +2M/r^3 per unit separation and together sideways at -M/r^3 (MTW 31.4). And the geodesic
// equation with Gamma (any spacetime) against Hamilton's equations (Geodesic.h) on a Kerr orbit.
void testTidesAndGammaGeodesic() {
    const double M = 1.0, r = 10.0;
    SchwarzschildSpacetime sw(M);
    const double x[4] = {0.0, r, 1.0, 0.0};
    const Tensor R = riemann(sw, x);
    const double f = 1 - 2 * M / r;
    const std::vector<double> u = {1 / std::sqrt(f), 0, 0, 0};
    const double radial = geodesicDeviation(R, u, {0, std::sqrt(f), 0, 0})[1] / std::sqrt(f); // orthonormal components
    const double side = geodesicDeviation(R, u, {0, 0, 1 / r, 0})[2] * r;
    const double tide = M / (r * r * r);
    std::printf("  tides at r = 10M: radial %.10e (2M/r^3 = %.10e), sideways %.10e (-M/r^3)\n", radial, 2 * tide, side);
    CHECK(std::fabs(radial / (2 * tide) - 1) < 1e-7 && std::fabs(side / tide + 1) < 1e-7, "tidal tensor wrong");
    // A bound Kerr orbit (that of testGeodesicInvariants) for lambda = 50 M with RK4, h = 0.05, in
    // both forms. u^m = g^mn p_n.
    Metric km;
    km.M = M;
    km.a = 0.9;
    Geodesic geo(km);
    GeodesicState s;
    s.x[1] = 10.0;
    s.x[2] = 0.5 * kPiD;
    const double dir[3] = {0.0, 0.3, std::sqrt(1.0 - 0.09)};
    Geodesic::momentumFromLocal(km, s.x, dir, 0.34, 1.0, s.p);
    double ginv[4][4];
    km.contravariant(s.x, ginv);
    std::vector<double> xs(s.x, s.x + 4), us(4, 0.0);
    for (int m = 0; m < 4; ++m)
        for (int n = 0; n < 4; ++n) us[size_t(m)] += ginv[m][n] * s.p[n];
    KerrSpacetime kerr(km);
    const double h = 0.05;
    for (int step = 0; step < 1000; ++step) {
        geo.stepRK4(s, h);
        geodesicStepRK4(kerr, xs, us, h);
    }
    double worst = 0;
    for (int m = 0; m < 4; ++m) worst = std::max(worst, std::fabs(xs[size_t(m)] - s.x[m]));
    std::printf("  Kerr orbit, lambda 50 M: Gamma form vs Hamilton form, max |dx| %.1e (r %.6f M, phi %.6f)\n", worst, s.x[1], s.x[3]);
    CHECK(worst < 1e-8, "the two forms of the geodesic equation disagree by %e", worst);
}
