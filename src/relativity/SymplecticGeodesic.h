#pragma once
// Symplectic integrators for geodesics: the energy error stays bounded for millions of orbits.
//
// Why. The motion along a geodesic IS Hamiltonian - H = 1/2 g^{mu nu} p_mu p_nu - so its exact flow
// keeps the symplectic form dx ^ dp. Runge-Kutta methods (Geodesic.h) do not: each step makes a
// small error of the same sign, and the energy creeps away, linearly in time. A symplectic step is
// instead the EXACT flow of a slightly different "shadow" Hamiltonian H~ = H + O(h^order) (backward
// error analysis; Hairer, Lubich & Wanner, "Geometric Numerical Integration", 2006, ch. IX), so
// H~ is kept and H only wobbles around its start - over any number of orbits. Two conditions:
// the force must be an exact gradient (KerrHamiltonian gives one) and the step must stay FIXED
// (a step that changes with the state is a different map every time and breaks the argument).
//
// The difficulty of general relativity. The simplest symplectic method, leapfrog (Stoermer-Verlet),
// needs H = T(p) + U(x). The Kerr Hamiltonian mixes them - g^{mu nu}(x) p p - it is not separable.
// Two ways around that, both here:
//   1. Implicit midpoint: z1 = z0 + h J grad H((z0 + z1) / 2). Symplectic and symmetric for ANY H,
//      order 2; each step solves for z1 by fixed-point iteration (Hairer, Lubich & Wanner, VI.3).
//      The simple, readable baseline.
//   2. Tao 2016 (Phys. Rev. E 94, 043303, doi:10.1103/PhysRevE.94.043303): explicit. Two copies of
//      phase space, (q, p) and (x, y), bound by a spring: H^(q,p,x,y) = H(q,y) + H(x,p) + omega
//      (|q - x|^2 + |p - y|^2) / 2. Each of the three terms has an exact flow (below), so their
//      symmetric composition is explicit and symplectic; the triple jump (Yoshida 1990) makes it
//      order 4. Tao proves an error O(T h^order omega) for integrable H - our Kerr orbits are
//      integrable (E, L, Q) - with omega above a threshold and h^order omega << 1.
// Why Tao and not the explicit Kerr splitting of Wu, Wang, Sun & Liu (2021, ApJ 914, 63,
// doi:10.3847/1538-4357/abfc45): theirs is faster per step but specific - a time transformation
// dtau = dlambda / Sigma and five hand-solved sub-flows of Kerr in Boyer-Lindquist form. Tao needs
// only H and its gradient, so it works for any spacetime we can write down, in about fifty lines.
//
// RK4 with the same exact gradient is here too: the fair comparison (same force, same cost).
#include "relativity/Geodesic.h"
#include "relativity/KerrHamiltonian.h"

namespace rf {

// Tao's extended phase space: the physical copy (q, p) = (s.x, s.p) and the shadow copy (x, y).
// Both start equal; the binding spring keeps them within O(1 / sqrt(omega)).
struct TaoState {
    GeodesicState s;
    double x[4] = {0, 0, 0, 0};
    double y[4] = {0, 0, 0, 0};
};

class SymplecticGeodesic {
public:
    explicit SymplecticGeodesic(const KerrHamiltonian& h) : H_(h) {} // a copy: two numbers
    const KerrHamiltonian& hamiltonian() const { return H_; }

    // One implicit-midpoint step of length h (may be negative). Returns the fixed-point iterations
    // it took to settle to rounding (each is one evaluation of the gradient).
    int stepImplicitMidpoint(GeodesicState& s, double h) const;

    // Tao's extended state with both copies at s.
    static TaoState taoStart(const GeodesicState& s);
    // One step of Tao's method of even `order` (2, 4, 6, ...) with binding constant `omega`.
    void stepTao(TaoState& z, double h, int order, double omega) const;

    // One classical Runge-Kutta step with the same exact gradient (not symplectic: the comparison).
    void stepRK4(GeodesicState& s, double h) const;

    // One radial period of a bound orbit that starts at its apoapsis (p_r = 0): the proper time to
    // fall in through the periapsis and climb back out, by fine RK4 steps.
    double radialPeriod(const GeodesicState& apoapsis) const;

    // Gradient evaluations so far: the cost of a run, whatever the method.
    long long gradientEvaluations() const { return evaluations_; }

private:
    // dz/dlambda = J grad H: dx = dH/dp, dp = -dH/dx.
    void field(const double x[4], const double p[4], double dx[4], double dp[4]) const;
    // The exact flows of Tao's three terms for a time d.
    void flowA(TaoState& z, double d) const;                 // H(q, y): moves p and x
    void flowB(TaoState& z, double d) const;                 // H(x, p): moves q and y
    static void flowC(TaoState& z, double d, double omega);  // the spring: rotates the difference
    void stepTao2(TaoState& z, double h, double omega) const;

    const KerrHamiltonian H_;
    mutable long long evaluations_ = 0;
};

} // namespace rf
