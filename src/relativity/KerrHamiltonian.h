#pragma once
// The Hamiltonian of a free particle or a light ray in the Kerr spacetime, in Boyer-Lindquist
// coordinates (t, r, theta, phi), with its EXACT gradient - what a symplectic integrator needs
// (SymplecticGeodesic.h).
//
// A geodesic is Hamilton's flow of H(x, p) = 1/2 g^{mu nu}(x) p_mu p_nu (Misner, Thorne & Wheeler,
// "Gravitation", 25.2). In Kerr the inverse metric has a closed form, and Carter (1968) wrote H so
// that it shows how the motion separates. With E = -p_t and L = p_phi (units G = c = 1):
//
//     Sigma = r^2 + a^2 cos^2 theta,          Delta = r^2 - 2 M r + a^2,
//     W     = (r^2 + a^2) E - a L,            V     = L - a E sin^2 theta,
//     2 Sigma H = Delta p_r^2 + p_theta^2 - W^2 / Delta + V^2 / sin^2 theta.
//
// Every partial derivative of this line is written out by hand in KerrHamiltonian.cpp - the chain
// rule and nothing else. Why not the finite differences Geodesic.cpp uses: a symplectic method keeps
// the energy bounded because each of its steps is the exact flow of a nearby Hamiltonian (backward
// error analysis; Hairer, Lubich & Wanner, "Geometric Numerical Integration", 2006, ch. IX). A force
// taken from differences is not exactly the gradient of any function, and that guarantee would be
// lost at the level of the differences' error. Exact derivatives keep it.
//
// t and phi do not appear in H: dp_t/dlambda = dp_phi/dlambda = 0 exactly, so E and L are kept to
// the last bit by every integrator (Noether: symmetry in time and about the axis). Valid outside the
// horizon (Delta > 0) and off the axis (sin theta != 0) - where Boyer-Lindquist coordinates are.

namespace rf {

class KerrHamiltonian {
public:
    KerrHamiltonian(double mass, double spin) : M(mass), a(spin) {}

    double M = 1.0; // mass
    double a = 0.0; // spin J / M, |a| <= M

    // H(x, p): 0 for light, -1/2 for a particle of unit mass (then lambda is its proper time).
    double value(const double x[4], const double p[4]) const;

    // The gradient: dHdx[mu] = dH/dx^mu and dHdp[mu] = dH/dp_mu. Hamilton's equations are then
    //   dx^mu/dlambda = dH/dp_mu,      dp_mu/dlambda = -dH/dx^mu.
    void gradient(const double x[4], const double p[4], double dHdx[4], double dHdp[4]) const;

    // Carter's radial potential of a unit mass, R(r) = W^2 - Delta (r^2 + (L - a E)^2 + Q): the
    // particle moves where R >= 0 and turns where R = 0 (Carter 1968; Bardeen, Press & Teukolsky
    // 1972). Q is Carter's constant.
    double radialPotential(double E, double L, double Q, double r) const;
    // E and L of the bound orbit of Carter's constant Q that turns at the radii rPeri and rApo
    // (R = 0 at both): Newton's method from a Newtonian guess. False if it did not converge.
    bool boundOrbit(double rPeri, double rApo, double Q, double& E, double& L) const;
};

} // namespace rf
