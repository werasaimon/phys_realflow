#pragma once
// Geodesics of a Metric: light rays and free particles as Hamilton's equations for the eight
// variables (x^mu, p_mu) with H = g^{mu nu} p_mu p_nu / 2 (Carter 1968; Misner, Thorne & Wheeler
// 25.2; the GYOTO and DNGR ray tracers integrate the same system):
//   dx^mu / dlambda = g^{mu nu} p_nu,      dp_mu / dlambda = -1/2 (d g^{ab} / d x^mu) p_a p_b.
// The affine parameter lambda is proper time for a unit-mass particle. Two integrators: classical
// RK4 with a fixed step (for the comparison) and Dormand-Prince RK45 with step control - the
// tests hold E, L and Carter's Q to 1e-8 over hundreds of orbits with tolerance 1e-10. The metric
// derivatives are 5-point central finite differences of g^{mu nu} in r and theta (the metric does
// not depend on t and phi): truncation ~ h^4, rounding ~ 1e-12 with h = 1e-4 r; simpler than the
// analytic Kerr derivatives and exact enough for the tolerances above. A symplectic integrator
// (bounded energy error over 1e6 orbits) is the next step.
#include "relativity/Metric.h"

#include <functional>

namespace rf {

struct GeodesicState {
    double x[4] = {0, 0, 0, 0}; // t (or v), r, theta, phi
    double p[4] = {0, 0, 0, 0}; // covariant momentum
    double lambda = 0;          // affine parameter
};

enum class GeodesicEnd { Running, Captured, Escaped, LambdaLimit, StepLimit, Stopped };

struct GeodesicResult {
    GeodesicEnd end = GeodesicEnd::Running;
    int steps = 0;
    double E0 = 0, L0 = 0, Q0 = 0, H0 = 0; // constants of motion at the start ...
    double E1 = 0, L1 = 0, Q1 = 0, H1 = 0; // ... and at the end: the drift is the integrator's error
};

class Geodesic {
public:
    struct Options {
        double rMin = -1;         // captured below this radius (< 0: the horizon times 1 + 1e-3)
        double rMax = 1000;       // escaped beyond this radius
        double lambdaMax = 1e5;   // and after this much affine parameter
        int maxSteps = 4000000;
        bool adaptive = true;     // RK45 with `tol`, else RK4 with `hFixed`
        double tol = 1e-10;       // RK45: absolute and relative tolerance per step
        double hInit = 0.1;       // RK45: first step; RK4: the step
        double hMax = 50.0;       // RK45: largest step
        double hFixed = 0.05;
        bool backward = false;    // integrate towards smaller lambda (ray tracing: where light came from)
        double mu = 0;            // 0 light, 1 unit-mass particle (for the constants of motion)
    };

    explicit Geodesic(const Metric& metric) : metric_(metric) {}
    const Metric& metric() const { return metric_; }

    // Hamilton's equations at s.
    void derivative(const GeodesicState& s, double dx[4], double dp[4]) const;
    // One RK4 step of length h (may be negative).
    void stepRK4(GeodesicState& s, double h) const;
    // One accepted Dormand-Prince step: |h| is shrunk until the error estimate is below tol, then
    // grown for the next call. Returns the number of attempts.
    int stepRK45(GeodesicState& s, double& h, double tol, double hMax) const;

    // Integrates until captured, escaped, the limits, or `stop(before, after)` returns true.
    GeodesicResult integrate(GeodesicState& s, const Options& opt,
                             const std::function<bool(const GeodesicState&, const GeodesicState&)>& stop = {}) const;

    // The momentum of a photon (or of a particle moving at `speed` < 1) that an observer at x
    // sees leaving along the unit direction dir = (n^r, n^theta, n^phi) of the observer's own
    // orthonormal frame, with local energy `energyLocal`. The frame is the ZAMO tetrad of the
    // Boyer-Lindquist coordinates (Bardeen, Press & Teukolsky 1972, eq. 3.3): the observer
    // co-rotating with the frame dragging, e_t = (d_t + omega d_phi) / alpha, e_i = d_i / sqrt(g_ii);
    // for a = 0 it is the static observer.
    static void momentumFromLocal(const Metric& m, const double x[4], const double dir[3], double speed,
                                  double energyLocal, double p[4]);

private:
    const Metric& metric_;
    void contravariantDerivatives(const double x[4], double ginv[4][4], double dg[4][4][4]) const;
};

} // namespace rf
