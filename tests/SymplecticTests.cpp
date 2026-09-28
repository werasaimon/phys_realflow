// Symplectic geodesics in Kerr (src/relativity/SymplecticGeodesic.h) against the facts that make
// them worth having:
//   1. the exact gradient of KerrHamiltonian is the metric's own (the finite differences of
//      Geodesic.cpp agree to their error);
//   2. over thousands of orbits the energy of the symplectic methods only wobbles, while RK4's
//      creeps away linearly - measured as the least-squares trend of H over the run divided by the
//      largest excursion (1 = pure drift, 0 = pure wobble);
//   3. the orders: implicit midpoint and Tao-2 are of order 2, Tao-4 and RK4 of order 4;
//   4. reversibility: a symmetric method run forward and then back lands on its start.
// The orbit: a unit-mass particle around a Kerr hole of spin a = 0.9 M, eccentric (turning points
// r = 8 M and 15 M) and inclined (Carter's Q = 4 M^2); one radial period is ~ 250 M of proper time.
//
// Registration (tests/Tests.h and tests/main.cpp):
//   run("symplectic geodesics: the exact Kerr gradient matches the metric", testKerrHamiltonianGradient);
//   run("symplectic geodesics: energy bounded over thousands of orbits, RK4 drifts", testSymplecticLongOrbits);
//   run("symplectic geodesics: order 2 (midpoint, Tao-2) and 4 (Tao-4, RK4)", testSymplecticOrder);
//   run("symplectic geodesics: forward then back returns to the start", testSymplecticReversibility);
// The long run behind the plot of docs/08-relativity.md:
//   RF_LONG_ORBITS=1000000 RF_PLOT_DIR=plots RF_TEST="symplectic geodesics: energy" ./rf_tests
//   python tools/plot_symplectic.py plots docs/img
#include "TestRunner.h"

#include "core/Parallel.h"
#include "relativity/Geodesic.h"
#include "relativity/KerrHamiltonian.h"
#include "relativity/Metric.h"
#include "relativity/SymplecticGeodesic.h"

#include <cmath>

using namespace rf;

namespace {

constexpr double kPiExact = 3.14159265358979323846;
constexpr double kMass = 1.0, kSpin = 0.9;           // the hole
constexpr double kPeriapsis = 8.0, kApoapsis = 15.0; // the orbit's turning points [M]
constexpr double kCarterQ = 4.0;                     // its inclination [M^2]
constexpr double kOmega = 0.1;                       // Tao's binding: ~8x the orbital frequency, h^4 omega << 1

// The particle at its apoapsis on the equator: p_r = 0, p_theta = sqrt(Q) (cos theta = 0 there),
// p_t = -E, p_phi = L; H = -1/2 follows from R(r_a) = 0.
GeodesicState eccentricOrbit() {
    double E, L;
    KerrHamiltonian(kMass, kSpin).boundOrbit(kPeriapsis, kApoapsis, kCarterQ, E, L); // R(r_p) = R(r_a) = 0
    GeodesicState s;
    s.x[1] = kApoapsis;
    s.x[2] = 0.5 * kPiExact;
    s.p[0] = -E;
    s.p[2] = std::sqrt(kCarterQ);
    s.p[3] = L;
    return s;
}

// The largest difference of the r and theta motion between two states (what the orbit IS; t and
// phi only add up along it).
double orbitDifference(const GeodesicState& a, const GeodesicState& b) {
    double worst = 0;
    for (int i = 1; i <= 2; ++i) worst = std::max({worst, std::fabs(a.x[i] - b.x[i]), std::fabs(a.p[i] - b.p[i])});
    return worst;
}

// The slope of ln(error) against ln(h): the observed order of a method (least squares).
double observedOrder(const std::vector<double>& h, const std::vector<double>& error) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    const double n = double(h.size());
    for (size_t i = 0; i < h.size(); ++i) {
        const double x = std::log(h[i]), y = std::log(error[i]);
        sx += x, sy += y, sxx += x * x, sxy += x * y;
    }
    return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

// A CSV for the plots, written only when RF_PLOT_DIR names a folder (a normal run writes nothing).
class DriftCsv {
public:
    DriftCsv(const char* name, const char* header) {
        const char* dir = std::getenv("RF_PLOT_DIR");
        if (!dir || !*dir) return;
        file_ = std::fopen((std::string(dir) + "/" + name + ".csv").c_str(), "w");
        if (file_) std::fprintf(file_, "%s\n", header);
    }
    ~DriftCsv() {
        if (file_) std::fclose(file_);
    }
    DriftCsv(const DriftCsv&) = delete;
    DriftCsv& operator=(const DriftCsv&) = delete;
    void row(std::initializer_list<double> values) {
        if (!file_) return;
        const char* separator = "";
        for (double v : values) std::fprintf(file_, "%s%.9g", separator, v), separator = ",";
        std::fprintf(file_, "\n");
    }

private:
    FILE* file_ = nullptr;
};

} // namespace

void testKerrHamiltonianGradient() {
    Metric m;
    m.M = kMass;
    m.a = kSpin;
    const Geodesic geo(m);
    const KerrHamiltonian H(kMass, kSpin);
    // Three points: the orbit's start, a point off the plane with radial motion, close in (r = 3 M).
    const double pts[3][4] = {{15.0, 0.5 * kPiExact, 0.0, 2.0}, {11.3, 1.1, 0.12, 1.3}, {3.0, 0.7, -0.4, 0.9}};
    double worstValue = 0, worstGradient = 0;
    for (const auto& q : pts) {
        GeodesicState s = eccentricOrbit();
        s.x[1] = q[0], s.x[2] = q[1], s.p[1] = q[2], s.p[2] = q[3];
        worstValue = std::max(worstValue, std::fabs(H.value(s.x, s.p) - m.hamiltonian(s.x, s.p)));
        double dHdx[4], dHdp[4], dx[4], dp[4];
        H.gradient(s.x, s.p, dHdx, dHdp);
        geo.derivative(s, dx, dp); // 5-point finite differences of g^{mu nu}
        for (int i = 0; i < 4; ++i) {
            worstGradient = std::max(worstGradient, std::fabs(dHdp[i] - dx[i]) / (1.0 + std::fabs(dx[i])));
            worstGradient = std::max(worstGradient, std::fabs(-dHdx[i] - dp[i]) / (1.0 + std::fabs(dp[i])));
        }
    }
    std::printf("  Carter's form vs the metric: H off by %.1e, the gradient off by %.1e (the finite differences' own error)\n",
                worstValue, worstGradient);
    CHECK(worstValue < 1e-13, "H of Carter's form is not the metric's: %e", worstValue);
    CHECK(worstGradient < 1e-8, "the exact gradient differs from the metric's derivatives by %e", worstGradient);
}

namespace {

// One method of the long run: how it steps and its step h [M].
struct Method {
    const char* name;
    int kind; // 0 RK4, 1 implicit midpoint, 2 Tao-4
    double h;
};

// What a long run measured.
struct Drift {
    double largest = 0;           // max |H - H0| over every step of the run
    double trendRatio = 0;        // |least-squares slope of H - H0| * T / largest: 1 drift, 0 wobble
    double largestQ = 0;          // max |Q / Q0 - 1| (Carter's constant, sampled)
    bool keptEL = true;           // p_t = -E and p_phi = L still their starting bits
    double gradientsPerStep = 0;  // the cost of a step
    double seconds = 0;
    std::vector<double> envelope; // max |H - H0| so far, at each checkpoint (in orbits)
};

// The checkpoints of the plot: 24 per decade of orbits, from 1 to `orbits`.
std::vector<double> checkpoints(double orbits) {
    std::vector<double> c;
    for (double o = 1.0; o <= orbits * (1 + 1e-9); o *= std::pow(10.0, 1.0 / 24.0)) c.push_back(o);
    if (c.empty() || c.back() < orbits * (1 - 1e-9)) c.push_back(orbits);
    return c;
}

// The least-squares line through (lambda, H - H0), accumulated step by step.
struct Trend {
    double st = 0, sy = 0, stt = 0, sty = 0, n = 0;
    void add(double t, double y) { st += t, sy += y, stt += t * t, sty += t * y, n += 1; }
    double slope() const { return (n * sty - st * sy) / (n * stt - st * st); }
};

// Runs one method for `orbits` radial periods of length `period` and measures its drift.
Drift runMethod(const Method& method, double orbits, double period, const std::vector<double>& marks) {
    const KerrHamiltonian H(kMass, kSpin);
    const SymplecticGeodesic sg(H);
    Metric m;
    m.M = kMass, m.a = kSpin;
    const GeodesicState start = eccentricOrbit();
    GeodesicState s = start;
    TaoState z = SymplecticGeodesic::taoStart(start);
    double e0, l0, q0, e, l, q;
    m.constants(start.x, start.p, 1.0, e0, l0, q0);
    const double h0 = H.value(start.x, start.p);
    const long long steps = std::llround(orbits * period / method.h), everyQ = std::max(1LL, steps / 4000);
    Drift d;
    Trend trend;
    size_t mark = 0;
    const auto t0 = std::chrono::steady_clock::now();
    for (long long k = 1; k <= steps; ++k) {
        // 1. One step of the method.
        if (method.kind == 0) {
            sg.stepRK4(s, method.h);
        } else if (method.kind == 1) {
            sg.stepImplicitMidpoint(s, method.h);
        } else {
            sg.stepTao(z, method.h, 4, kOmega);
            s = z.s;
        }
        // 2. Its energy error, every step: the largest excursion and the trend line.
        const double lambda = double(k) * method.h, dH = H.value(s.x, s.p) - h0;
        d.largest = std::max(d.largest, std::fabs(dH));
        trend.add(lambda, dH);
        // 3. The envelope at the plot's checkpoints; Carter's Q now and then.
        for (; mark < marks.size() && lambda >= marks[mark] * period; ++mark) d.envelope.push_back(d.largest);
        if (k % everyQ == 0) {
            m.constants(s.x, s.p, 1.0, e, l, q);
            d.largestQ = std::max(d.largestQ, std::fabs(q / q0 - 1.0));
        }
    }
    for (; mark < marks.size(); ++mark) d.envelope.push_back(d.largest); // the last one, lost to rounding
    d.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    d.trendRatio = std::fabs(trend.slope() * double(steps) * method.h) / std::max(d.largest, 1e-300);
    d.keptEL = s.p[0] == start.p[0] && s.p[3] == start.p[3];
    // The H measured every step above is a value, not a gradient: not counted, not the method's cost.
    d.gradientsPerStep = double(sg.gradientEvaluations()) / double(steps);
    return d;
}

} // namespace

void testSymplecticLongOrbits() {
    const double period = SymplecticGeodesic(KerrHamiltonian(kMass, kSpin)).radialPeriod(eccentricOrbit());
    const char* env = std::getenv("RF_LONG_ORBITS");
    const double orbits = env && *env ? std::atof(env) : 2000.0;
    // RK4 twice: at the symplectic methods' step, and at a third of it - the same cost as Tao-4
    // (12 gradients a step against RK4's 4).
    const Method methods[4] = {{"RK4, h = 1 M", 0, 1.0}, {"RK4, h = 1/3 M (cost of Tao-4)", 0, 1.0 / 3.0},
                               {"implicit midpoint, h = 1 M", 1, 1.0}, {"Tao-4, h = 1 M, omega = 0.1", 2, 1.0}};
    const std::vector<double> marks = checkpoints(orbits);
    Drift d[4];
    parallelFor(4, [&](int i) { d[i] = runMethod(methods[i], orbits, period, marks); }, 1);
    std::printf("  Kerr a = 0.9, r 8..15 M, Q = 4, radial period %.2f M, %.0f orbits:\n", period, orbits);
    for (int i = 0; i < 4; ++i)
        std::printf("    %-32s max|dH| %.2e, trend/excursion %.2f, max|dQ/Q| %.2e, E and L %s, %4.1f gradients a step, %.1f s\n",
                    methods[i].name, d[i].largest, d[i].trendRatio, d[i].largestQ, d[i].keptEL ? "exact" : "CHANGED",
                    d[i].gradientsPerStep, d[i].seconds);
    // Where RK4 at the same cost overtakes Tao-4's bounded wobble: its error grows linearly.
    std::printf("    RK4 at Tao-4's cost passes Tao-4's excursion after ~%.0f orbits\n",
                orbits * d[3].largest / std::max(d[1].largest, 1e-300));
    DriftCsv csv("symplectic_drift", "orbits,rk4_h1,rk4_same_cost,implicit_midpoint,tao4");
    for (size_t k = 0; k < marks.size() && k < d[0].envelope.size() && k < d[1].envelope.size(); ++k)
        csv.row({marks[k], d[0].envelope[k], d[1].envelope[k], d[2].envelope[k], d[3].envelope[k]});
    CHECK(d[0].trendRatio > 0.9, "RK4's energy error should be a drift (trend/excursion %f)", d[0].trendRatio);
    CHECK(d[2].trendRatio < 0.2 && d[3].trendRatio < 0.2, "a symplectic energy error drifts: midpoint %f, Tao-4 %f",
          d[2].trendRatio, d[3].trendRatio);
    CHECK(d[2].largest < 5e-5 && d[3].largest < 3e-5, "the symplectic excursions are too large: %e, %e", d[2].largest, d[3].largest);
    CHECK(d[0].keptEL && d[1].keptEL && d[2].keptEL && d[3].keptEL, "E or L changed: t and phi are cyclic, they must not");
}

void testSymplecticOrder() {
    const KerrHamiltonian H(kMass, kSpin);
    const SymplecticGeodesic sg(H);
    const double T = 256.0; // about one radial period
    // 1. The reference: RK4 with a step 64 times finer than the finest tested (its error ~ 1e-15).
    GeodesicState reference = eccentricOrbit();
    for (int k = 0; k < 256 * 64; ++k) sg.stepRK4(reference, 1.0 / 64);
    // 2. Each method at h = 1, 1/2, 1/4 M; the error of the r and theta motion at T.
    const std::vector<double> steps = {1.0, 0.5, 0.25};
    std::vector<double> rk4, midpoint, tao2, tao4;
    for (double h : steps) {
        GeodesicState a = eccentricOrbit(), b = a;
        TaoState c = SymplecticGeodesic::taoStart(a), e = c;
        for (int k = 0; k < int(T / h); ++k) {
            sg.stepRK4(a, h);
            sg.stepImplicitMidpoint(b, h);
            sg.stepTao(c, h, 2, kOmega);
            sg.stepTao(e, h, 4, kOmega);
        }
        rk4.push_back(orbitDifference(a, reference));
        midpoint.push_back(orbitDifference(b, reference));
        tao2.push_back(orbitDifference(c.s, reference));
        tao4.push_back(orbitDifference(e.s, reference));
    }
    // 3. The orders: the slope of ln(error) against ln(h).
    const double oRK4 = observedOrder(steps, rk4), oMid = observedOrder(steps, midpoint);
    const double oTao2 = observedOrder(steps, tao2), oTao4 = observedOrder(steps, tao4);
    std::printf("  error at T = 256 M for h = 1, 1/2, 1/4 M:\n"
                "    implicit midpoint %.2e %.2e %.2e -> order %.2f (theory 2)\n"
                "    Tao-2             %.2e %.2e %.2e -> order %.2f (theory 2)\n"
                "    Tao-4             %.2e %.2e %.2e -> order %.2f (theory 4)\n"
                "    RK4               %.2e %.2e %.2e -> order %.2f (theory 4)\n",
                midpoint[0], midpoint[1], midpoint[2], oMid, tao2[0], tao2[1], tao2[2], oTao2, tao4[0], tao4[1], tao4[2], oTao4,
                rk4[0], rk4[1], rk4[2], oRK4);
    CHECK(std::fabs(oMid - 2) < 0.2 && std::fabs(oTao2 - 2) < 0.2, "second-order methods: midpoint %f, Tao-2 %f", oMid, oTao2);
    CHECK(std::fabs(oTao4 - 4) < 0.4 && std::fabs(oRK4 - 4) < 0.4, "fourth-order methods: Tao-4 %f, RK4 %f", oTao4, oRK4);
}

void testSymplecticReversibility() {
    const KerrHamiltonian H(kMass, kSpin);
    const SymplecticGeodesic sg(H);
    const GeodesicState start = eccentricOrbit();
    const int n = 2000; // about 8 orbits out, and 8 back
    GeodesicState rk4 = start, midpoint = start;
    TaoState tao2 = SymplecticGeodesic::taoStart(start), tao4 = tao2;
    for (int k = 0; k < n; ++k) {
        sg.stepRK4(rk4, 1.0);
        sg.stepImplicitMidpoint(midpoint, 1.0);
        sg.stepTao(tao2, 1.0, 2, kOmega);
        sg.stepTao(tao4, 1.0, 4, kOmega);
    }
    for (int k = 0; k < n; ++k) { // the same steps backwards: h = -1
        sg.stepRK4(rk4, -1.0);
        sg.stepImplicitMidpoint(midpoint, -1.0);
        sg.stepTao(tao2, -1.0, 2, kOmega);
        sg.stepTao(tao4, -1.0, 4, kOmega);
    }
    const double eRK4 = orbitDifference(rk4, start), eMid = orbitDifference(midpoint, start);
    const double eTao2 = orbitDifference(tao2.s, start), eTao4 = orbitDifference(tao4.s, start);
    std::printf("  2000 steps of 1 M out and back: midpoint %.1e, Tao-2 %.1e, Tao-4 %.1e from the start; RK4 %.1e (not symmetric)\n",
                eMid, eTao2, eTao4, eRK4);
    CHECK(eMid < 1e-11 && eTao2 < 1e-11 && eTao4 < 1e-11, "a symmetric method did not come back: %e %e %e", eMid, eTao2, eTao4);
}
