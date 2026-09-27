#pragma once
// Spacetime of a rotating black hole: the Kerr metric (Kerr 1963) in Boyer-Lindquist coordinates
// (t, r, theta, phi), and for the non-rotating case (Schwarzschild) also the ingoing
// Eddington-Finkelstein coordinates (v, r, theta, phi) that are regular on the horizon.
// Geometric units G = c = 1: lengths and times in units of the mass M (for the Sun M = 1.48 km =
// 4.9 us; for M87* 6.5e9 solar masses, M = 9.6e12 m). Everything is double precision: the tests
// hold the constants of motion to 1e-8 over hundreds of orbits.
//
// Boyer-Lindquist form (Bardeen, Press & Teukolsky 1972, ApJ 178, 347; Misner, Thorne & Wheeler
// "Gravitation" 33.2), with Sigma = r^2 + a^2 cos^2 theta and Delta = r^2 - 2 M r + a^2:
//   ds^2 = -(1 - 2 M r / Sigma) dt^2 - (4 M a r sin^2 theta / Sigma) dt dphi + (Sigma / Delta) dr^2
//          + Sigma dtheta^2 + (r^2 + a^2 + 2 M a^2 r sin^2 theta / Sigma) sin^2 theta dphi^2.
// The horizon is the larger root of Delta: r+ = M + sqrt(M^2 - a^2); it is a coordinate
// singularity of these coordinates (t runs to infinity there), which is why the infalling test
// uses Eddington-Finkelstein: v = t + r + 2 M ln(r / 2M - 1) turns the Schwarzschild metric into
//   ds^2 = -(1 - 2 M / r) dv^2 + 2 dv dr + r^2 dOmega^2,
// finite through r = 2 M (Finkelstein 1958) - the kind of coordinates a GRMHD code needs.
#include <cmath>

namespace rf {

class Metric {
public:
    enum class Coordinates { BoyerLindquist, EddingtonFinkelstein };

    double M = 1.0; // mass
    double a = 0.0; // spin parameter J / M, |a| <= M; 0 = Schwarzschild
    Coordinates coordinates = Coordinates::BoyerLindquist;

    // g_{mu nu} and g^{mu nu} at the point x = (t, r, theta, phi) (or (v, r, theta, phi)).
    void covariant(const double x[4], double g[4][4]) const;
    void contravariant(const double x[4], double ginv[4][4]) const;

    double horizonRadius() const { return M + std::sqrt(std::max(M * M - a * a, 0.0)); } // r+
    // Circular photon orbit in the equatorial plane (Bardeen, Press & Teukolsky 1972, eq. 2.18):
    // r_ph = 2 M (1 + cos(2/3 arccos(-+ a / M))), 3 M for a = 0, M for a = M prograde.
    double photonSphereRadius(bool prograde = true) const;
    // Innermost stable circular orbit (ibid. eq. 2.21): 6 M for a = 0, M (prograde) / 9 M
    // (retrograde) for a = M.
    double iscoRadius(bool prograde = true) const;
    // Radius of the shadow of a non-rotating hole seen from far away (Bardeen 1973): the critical
    // impact parameter 3 sqrt(3) M of the photon sphere.
    static double shadowRadiusSchwarzschild(double M) { return 3.0 * std::sqrt(3.0) * M; }
    // Energy and angular momentum of the circular equatorial orbit of a massive particle at r
    // (ibid. eq. 2.12-2.13); prograde. Exact inputs for the ISCO and precession tests.
    void circularOrbit(double r, double& E, double& L, bool prograde = true) const;
    // Keplerian angular velocity of that orbit, Omega = +-sqrt(M) / (r^{3/2} +- a sqrt(M)).
    double keplerianOmega(double r, bool prograde = true) const;

    // The constants of motion of a geodesic with covariant momentum p at x: the energy E = -p_t,
    // the axial angular momentum L = p_phi and Carter's constant (Carter 1968)
    //   Q = p_theta^2 + cos^2 theta (a^2 (mu^2 - E^2) + L^2 / sin^2 theta),
    // mu = 0 for light, 1 for a unit-mass particle. All three are conserved exactly along the
    // geodesic: E and L by the symmetries in t and phi (Noether), Q by the hidden symmetry of
    // Kerr - the drift of the three numbers measures the integrator, nothing else.
    void constants(const double x[4], const double p[4], double mu, double& E, double& L, double& Q) const;
    // The Hamiltonian H = g^{mu nu} p_mu p_nu / 2: 0 for light, -1/2 for a unit-mass particle.
    double hamiltonian(const double x[4], const double p[4]) const;

    // SI helpers: the Schwarzschild radius 2 G M / c^2 of a mass in kilograms, and the unit of
    // time G M / c^3.
    static double schwarzschildRadiusSI(double massKg);
    static double timeUnitSI(double massKg);
};

// A 4x4 symmetric matrix inverse (Gauss-Jordan); false if singular.
bool invert4(const double m[4][4], double out[4][4]);

} // namespace rf
