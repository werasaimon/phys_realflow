// The Kerr Hamiltonian in Carter's separated form and its exact gradient. See KerrHamiltonian.h.
#include "relativity/KerrHamiltonian.h"

#include <algorithm>
#include <cmath>

namespace rf {

namespace {

// The building blocks of Carter's form at one point of phase space, computed once and shared by
// the value and the gradient.
struct Carter {
    double sinTheta, cosTheta; // of the polar angle
    double sigma;              // Sigma = r^2 + a^2 cos^2 theta
    double delta;              // Delta = r^2 - 2 M r + a^2 (zero on the horizons)
    double w;                  // W = (r^2 + a^2) E - a L: the radial "energy" of the orbit
    double v;                  // V = L - a E sin^2 theta: the polar "angular momentum"
    double n;                  // N = 2 Sigma H = Delta p_r^2 + p_theta^2 - W^2 / Delta + V^2 / sin^2 theta
};

Carter carter(double M, double a, const double x[4], const double p[4]) {
    Carter c;
    const double r = x[1], energy = -p[0], angular = p[3];
    c.sinTheta = std::sin(x[2]);
    c.cosTheta = std::cos(x[2]);
    const double sin2 = c.sinTheta * c.sinTheta;
    c.sigma = r * r + a * a * c.cosTheta * c.cosTheta;
    c.delta = r * r - 2.0 * M * r + a * a;
    c.w = (r * r + a * a) * energy - a * angular;
    c.v = angular - a * energy * sin2;
    c.n = c.delta * p[1] * p[1] + p[2] * p[2] - c.w * c.w / c.delta + c.v * c.v / sin2;
    return c;
}

} // namespace

double KerrHamiltonian::value(const double x[4], const double p[4]) const {
    const Carter c = carter(M, a, x, p);
    return c.n / (2.0 * c.sigma); // H = N / (2 Sigma)
}

void KerrHamiltonian::gradient(const double x[4], const double p[4], double dHdx[4], double dHdp[4]) const {
    const Carter c = carter(M, a, x, p);
    const double r = x[1], energy = -p[0], radial = p[1], polar = p[2];
    const double sin2 = c.sinTheta * c.sinTheta;
    const double twoSigma2 = 2.0 * c.sigma * c.sigma;

    // 1. The coordinates t and phi do not appear: their forces are exactly zero.
    dHdx[0] = dHdx[3] = 0.0;

    // 2. Along r. H = N / (2 Sigma), so dH/dr = (N_r Sigma - N Sigma_r) / (2 Sigma^2), with
    //    Delta_r = 2r - 2M, W_r = 2 r E, Sigma_r = 2r and
    //    N_r = Delta_r p_r^2 - 2 W W_r / Delta + W^2 Delta_r / Delta^2.
    const double deltaR = 2.0 * r - 2.0 * M, wR = 2.0 * r * energy;
    const double nR = deltaR * radial * radial - 2.0 * c.w * wR / c.delta + c.w * c.w * deltaR / (c.delta * c.delta);
    dHdx[1] = (nR * c.sigma - c.n * 2.0 * r) / twoSigma2;

    // 3. Along theta: V_theta = -2 a E sin cos, Sigma_theta = -2 a^2 sin cos and
    //    N_theta = d(V^2 / sin^2)/dtheta = 2 V V_theta / sin^2 - 2 V^2 cos / sin^3.
    const double sinCos = c.sinTheta * c.cosTheta;
    const double vTheta = -2.0 * a * energy * sinCos, sigmaTheta = -2.0 * a * a * sinCos;
    const double nTheta = 2.0 * c.v * vTheta / sin2 - 2.0 * c.v * c.v * c.cosTheta / (sin2 * c.sinTheta);
    dHdx[2] = (nTheta * c.sigma - c.n * sigmaTheta) / twoSigma2;

    // 4. The velocities dx/dlambda = dH/dp. For p_r and p_theta H is a plain quadratic; p_t and
    //    p_phi enter through W and V (E = -p_t, L = p_phi):
    //      dH/dp_t   =  (W (r^2 + a^2) / Delta + a V) / Sigma     (so dt/dlambda > 0 outside),
    //      dH/dp_phi =  (a W / Delta + V / sin^2 theta) / Sigma.
    dHdp[1] = c.delta * radial / c.sigma;
    dHdp[2] = polar / c.sigma;
    dHdp[0] = (c.w * (r * r + a * a) / c.delta + a * c.v) / c.sigma;
    dHdp[3] = (a * c.w / c.delta + c.v / sin2) / c.sigma;
}

double KerrHamiltonian::radialPotential(double E, double L, double Q, double r) const {
    const double delta = r * r - 2.0 * M * r + a * a, w = (r * r + a * a) * E - a * L;
    return w * w - delta * (r * r + (L - a * E) * (L - a * E) + Q);
}

// Two equations R(r_p) = R(r_a) = 0 in two unknowns (E, L). Newton's method with a Jacobian from
// forward differences; the guess is the Newtonian orbit (E^2 ~ 1 - M / a_semi, L^2 ~ M p).
bool KerrHamiltonian::boundOrbit(double rPeri, double rApo, double Q, double& E, double& L) const {
    const double semiMajor = 0.5 * (rPeri + rApo), semiLatus = 2.0 * rPeri * rApo / (rPeri + rApo);
    E = std::sqrt(1.0 - M / semiMajor);
    L = std::sqrt(M * semiLatus);
    for (int iteration = 0; iteration < 60; ++iteration) {
        const double f1 = radialPotential(E, L, Q, rPeri), f2 = radialPotential(E, L, Q, rApo);
        const double d = 1e-7 * std::max(1.0, std::fabs(L));
        const double a11 = (radialPotential(E + d, L, Q, rPeri) - f1) / d, a12 = (radialPotential(E, L + d, Q, rPeri) - f1) / d;
        const double a21 = (radialPotential(E + d, L, Q, rApo) - f2) / d, a22 = (radialPotential(E, L + d, Q, rApo) - f2) / d;
        const double det = a11 * a22 - a12 * a21;
        if (det == 0.0) return false;
        const double dE = (a22 * f1 - a12 * f2) / det, dL = (a11 * f2 - a21 * f1) / det;
        E -= dE;
        L -= dL;
        if (std::fabs(dE) + std::fabs(dL) < 1e-15) return true;
    }
    return std::fabs(radialPotential(E, L, Q, rPeri)) + std::fabs(radialPotential(E, L, Q, rApo)) < 1e-9;
}

} // namespace rf
