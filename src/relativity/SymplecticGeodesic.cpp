// Implicit midpoint, Tao's explicit symplectic method and RK4 for Kerr geodesics, all with the
// exact gradient of KerrHamiltonian. See SymplecticGeodesic.h.
#include "relativity/SymplecticGeodesic.h"

#include <algorithm>
#include <cmath>

namespace rf {

void SymplecticGeodesic::field(const double x[4], const double p[4], double dx[4], double dp[4]) const {
    double dHdx[4], dHdp[4];
    H_.gradient(x, p, dHdx, dHdp);
    ++evaluations_;
    for (int i = 0; i < 4; ++i) {
        dx[i] = dHdp[i];  // dx^mu / dlambda =  dH / dp_mu
        dp[i] = -dHdx[i]; // dp_mu / dlambda = -dH / dx^mu
    }
}

// ---------------------------------------------------------------------------------------------
// Implicit midpoint: the end z1 is where the field at the middle (z0 + z1) / 2 carries z0 in h.
int SymplecticGeodesic::stepImplicitMidpoint(GeodesicState& s, double h) const {
    constexpr int kMaxIterations = 60;
    constexpr double kSettled = 1e-15; // relative change of the guess: rounding, nothing more
    // 1. The first guess of the end: an explicit Euler step.
    double dx[4], dp[4], x1[4], p1[4];
    field(s.x, s.p, dx, dp);
    for (int i = 0; i < 4; ++i) {
        x1[i] = s.x[i] + h * dx[i];
        p1[i] = s.p[i] + h * dp[i];
    }
    // 2. Improve it: z1 <- z0 + h f((z0 + z1) / 2), until it stops changing. The map contracts
    //    with a factor ~ h |df/dz| / 2, a few digits per sweep at orbital steps.
    int iterations = 0;
    for (double change = 1.0; change > kSettled && iterations < kMaxIterations; ++iterations) {
        double xm[4], pm[4];
        for (int i = 0; i < 4; ++i) {
            xm[i] = 0.5 * (s.x[i] + x1[i]);
            pm[i] = 0.5 * (s.p[i] + p1[i]);
        }
        field(xm, pm, dx, dp);
        change = 0.0;
        for (int i = 0; i < 4; ++i) {
            const double nx = s.x[i] + h * dx[i], np = s.p[i] + h * dp[i];
            change = std::max(change, std::fabs(nx - x1[i]) / (1.0 + std::fabs(nx)));
            change = std::max(change, std::fabs(np - p1[i]) / (1.0 + std::fabs(np)));
            x1[i] = nx;
            p1[i] = np;
        }
    }
    // 3. The settled guess is the step.
    for (int i = 0; i < 4; ++i) {
        s.x[i] = x1[i];
        s.p[i] = p1[i];
    }
    s.lambda += h;
    return iterations;
}

// ---------------------------------------------------------------------------------------------
// Tao 2016, eq. (1): the exact flows of the three terms of H^ = H(q, y) + H(x, p) + omega H_C.

TaoState SymplecticGeodesic::taoStart(const GeodesicState& s) {
    TaoState z;
    z.s = s;
    for (int i = 0; i < 4; ++i) {
        z.x[i] = s.x[i];
        z.y[i] = s.p[i];
    }
    return z;
}

// H_A = H(q, y) does not depend on p or x, so for a time d its flow keeps q and y and moves
//   p <- p - d dH/dq (q, y),     x <- x + d dH/dy (q, y).
void SymplecticGeodesic::flowA(TaoState& z, double d) const {
    double dx[4], dp[4];
    field(z.s.x, z.y, dx, dp);
    for (int i = 0; i < 4; ++i) {
        z.s.p[i] += d * dp[i];
        z.x[i] += d * dx[i];
    }
}

// H_B = H(x, p) keeps x and p and moves
//   q <- q + d dH/dp (x, p),     y <- y - d dH/dx (x, p).
void SymplecticGeodesic::flowB(TaoState& z, double d) const {
    double dx[4], dp[4];
    field(z.x, z.s.p, dx, dp);
    for (int i = 0; i < 4; ++i) {
        z.s.x[i] += d * dx[i];
        z.y[i] += d * dp[i];
    }
}

// omega H_C = omega (|q - x|^2 + |p - y|^2) / 2 is a harmonic spring between the copies: the sums
// q + x, p + y stay, the differences (q - x, p - y) turn by the angle 2 omega d:
//   (q - x, p - y) <- R (q - x, p - y),   R = [cos 2wd, sin 2wd; -sin 2wd, cos 2wd].
// The spring binds only r and theta. H does not depend on t and phi, so their momenta (-E and L)
// never change in either copy, and nothing flows back from t and phi into the motion: they need no
// binding. Bound as well, the spring would turn their growing difference into p_t and p_phi and
// spoil E and L - which, unbound, stay exact to the last bit (the extended Hamiltonian stays
// Hamiltonian and its exact flow still keeps the two copies equal).
void SymplecticGeodesic::flowC(TaoState& z, double d, double omega) {
    const double c = std::cos(2.0 * omega * d), s = std::sin(2.0 * omega * d);
    for (int i = 1; i <= 2; ++i) { // r and theta
        const double sumQ = z.s.x[i] + z.x[i], sumP = z.s.p[i] + z.y[i];
        const double difQ = z.s.x[i] - z.x[i], difP = z.s.p[i] - z.y[i];
        const double turnQ = c * difQ + s * difP, turnP = -s * difQ + c * difP;
        z.s.x[i] = 0.5 * (sumQ + turnQ);
        z.s.p[i] = 0.5 * (sumP + turnP);
        z.x[i] = 0.5 * (sumQ - turnQ);
        z.y[i] = 0.5 * (sumP - turnP);
    }
}

// Tao 2016, eq. (2): the symmetric (Strang) composition A(h/2) B(h/2) C(h) B(h/2) A(h/2), order 2.
void SymplecticGeodesic::stepTao2(TaoState& z, double h, double omega) const {
    flowA(z, 0.5 * h);
    flowB(z, 0.5 * h);
    flowC(z, h, omega);
    flowB(z, 0.5 * h);
    flowA(z, 0.5 * h);
}

// The triple jump (Yoshida 1990; Tao 2016, eq. 3): three symmetric steps of order p = l - 2 with
// lengths gamma h, (1 - 2 gamma) h, gamma h make one of order l when their leading errors cancel,
//   2 gamma^(p+1) + (1 - 2 gamma)^(p+1) = 0   =>   gamma = 1 / (2 - 2^(1 / (p + 1))),
// 1 / (2 - 2^(1/3)) = 1.3512 for order 4. The middle step goes backwards (1 - 2 gamma = -1.70).
// (Tao's eq. 3 prints the exponent as 1/(l+1); the condition above, and a measured order of 2
// instead of 4 with that exponent, say 1/(p+1) = 1/(l-1).)
void SymplecticGeodesic::stepTao(TaoState& z, double h, int order, double omega) const {
    if (order <= 2) {
        stepTao2(z, h, omega);
        z.s.lambda += h; // the innermost steps count the time; their lengths add up to h
        return;
    }
    const int base = order - 2; // the order of the three steps composed
    const double gamma = 1.0 / (2.0 - std::pow(2.0, 1.0 / (base + 1)));
    stepTao(z, gamma * h, base, omega);
    stepTao(z, (1.0 - 2.0 * gamma) * h, base, omega);
    stepTao(z, gamma * h, base, omega);
}

// ---------------------------------------------------------------------------------------------
// Classical Runge-Kutta: four evaluations of the field, order 4, not symplectic.
void SymplecticGeodesic::stepRK4(GeodesicState& s, double h) const {
    double kx[4][4], kp[4][4];
    const double weight[4] = {0.0, 0.5, 0.5, 1.0};
    for (int stage = 0; stage < 4; ++stage) {
        double x[4], p[4];
        for (int i = 0; i < 4; ++i) {
            x[i] = s.x[i] + (stage ? weight[stage] * h * kx[stage - 1][i] : 0.0);
            p[i] = s.p[i] + (stage ? weight[stage] * h * kp[stage - 1][i] : 0.0);
        }
        field(x, p, kx[stage], kp[stage]);
    }
    for (int i = 0; i < 4; ++i) {
        s.x[i] += h / 6.0 * (kx[0][i] + 2.0 * kx[1][i] + 2.0 * kx[2][i] + kx[3][i]);
        s.p[i] += h / 6.0 * (kp[0][i] + 2.0 * kp[1][i] + 2.0 * kp[2][i] + kp[3][i]);
    }
    s.lambda += h;
}

// p_r turns negative (falling in), crosses to positive at the periapsis and back to negative at
// the next apoapsis; that last crossing, interpolated inside its step, ends the period.
double SymplecticGeodesic::radialPeriod(const GeodesicState& apoapsis) const {
    constexpr double kStep = 0.01;
    constexpr int kMaxSteps = 100000000; // 10^6 M of proper time: far longer than any bound orbit here
    GeodesicState s = apoapsis;
    s.lambda = 0.0;
    bool climbing = false;
    for (int k = 0; k < kMaxSteps; ++k) {
        const GeodesicState before = s;
        stepRK4(s, kStep);
        if (before.p[1] < 0 && s.p[1] >= 0) climbing = true;
        if (climbing && before.p[1] > 0 && s.p[1] <= 0) return before.lambda + kStep * before.p[1] / (before.p[1] - s.p[1]);
    }
    return 0.0; // not a bound orbit
}

} // namespace rf
