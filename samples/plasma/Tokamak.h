#pragma once
// A tokamak: a ring of plasma held by a magnetic field inside a toroidal vessel (Artsimovich 1972;
// Wesson, "Tokamaks", 4th ed. 2011, ch. 3 and 6; Freidberg, "Ideal MHD", 2014, ch. 11).
// Coordinates about the vertical axis y: the major radius R (distance from the axis), the
// toroidal angle phi (around the axis); in the poloidal plane the minor radius r (distance from
// the magnetic axis, the circle R = R0, y = 0) and the poloidal angle theta. Three fields make
// the ring:
//   * the toroidal field of the coils,  B_phi = B0 R0 / R  (current-free between the coils);
//   * the poloidal field B_theta(r) of the plasma current I_p flowing along phi, here with the
//     constant profile of the classic kink analysis, j = I_p / (pi a^2) in r < a:
//     B_theta = mu0 I_p r / (2 pi a^2) inside, mu0 I_p / (2 pi r) outside (the straight-column
//     values; the torus bends them by ~r / R0). Its vector potential is the sum of circular
//     current loops filling the channel, each the exact loop potential in elliptic integrals
//     (Jackson, "Classical Electrodynamics", eq. 5.37) - the straight column's A(r) phi_hat wrapped
//     into a torus is not current-free outside the channel and carries a wrong vertical field;
//   * the vertical field B_v that holds the ring against its own outward hoop force. The vessel
//     wall is a perfectly conducting shell, and its image currents hold the ring too: without
//     B_v the ring sits shifted outwards by (Mukhovatov & Shafranov 1971, Nucl. Fusion 11)
//       Delta = b^2 / (2 R0) [ ln(b / a) + (Lambda + 1/2) (1 - a^2 / b^2) ],  Lambda = beta_p + l_i / 2 - 1,
//     and the vertical field that centres it, from the shell's restoring force mu0 I_p^2 Delta / (2 pi b^2)
//     per unit length, is  B_v = mu0 I_p Delta / (2 pi b^2)  (l_i = 1/2 for the constant profile).
//     As in a real machine the vertical field is under feedback (radial position control): the
//     coils add gain x shift + damping x its rate to B_v, so the ring, never exactly in the
//     equilibrium of the formulas on a grid, does not swing in and out.
// Between the current channel (r < a) and the wall (r = b) lies a resistive "vacuum": the same
// fluid, but with a resistivity so large that the field slips through it within an Alfven time,
// as in the halo region of the tokamak codes (NIMROD, M3D-C1). Without it the perfectly conducting
// fluid there would hold the field lines and stabilise the kink.
// The field lines wind helically; the safety factor q = r B_phi / (R0 B_theta) - constant q_a
// inside the channel - is the number of toroidal turns a line makes per poloidal turn. The m = 1,
// n = 1 kink (Kruskal 1954, Shafranov 1956): for the constant profile with a conducting wall at b,
//   delta W = W0 (1 - q_a) [ (1 - q_a) (b^2 + a^2) / (b^2 - a^2) - 1 ],   W0 = 2 pi^2 R0 B_theta(a)^2 xi^2 / mu0,
// negative - unstable - for  2 a^2 / (a^2 + b^2) < q_a < 1: above q_a = 1 the line tension wins
// (the Kruskal-Shafranov limit), below the wall limit the image currents in the wall hold the
// column. The column shifts rigidly, the fluid outside flows around it (added mass with the wall
// (b^2 + a^2) / (b^2 - a^2) of the displaced mass), so
//   gamma^2 = (v_A,theta / a)^2 (1 - q_a) (q_a (b^2 + a^2) - 2 a^2) / b^2,  v_A,theta = B_theta(a) / sqrt(mu0 rho).
// This already counts the gas around the column (the scene's "vacuum" is gas of the same density):
// with a massless vacuum the kink would grow 1.63x faster (kinkGrowthRateWithOuterGas). The
// measured kink grows at 0.4-0.6 of this gamma; what explains the gap is in docs/06-mhd-plasma.md.
// Everything is scaled down (millitesla, metres, Alfven speeds of m/s): the equations are the
// same, only the time runs 10^6 times slower than in a real machine.

#include "plasma/MagneticField.h"
#include "math/Math.h"

#include <vector>

namespace rf {

class Tokamak {
public:
    Vector3 centre{0.0f};
    float majorRadius = 0.4f;      // R0 [m]
    float minorRadius = 0.12f;     // a: radius of the current channel [m]
    float vesselRadius = 0.24f;    // b: the wall of the vessel in the poloidal plane [m]
    float toroidalField = 2.5e-3f; // B0 = B_phi on the magnetic axis [T]: v_A = 2.2 m/s at density 1
    float safetyFactorEdge = 0.7f; // q_a: sets the plasma current; 0.4 < q_a < 1 kinks with b = 2a
    bool verticalField = true;     // the equilibrium B_v (false: the wall's image currents alone hold the ring)
    bool positionControl = true;   // feedback on B_v holding the ring on the magnetic axis
    // [m^2/s] of the fluid between the channel and the wall. A vacuum holds no field lines: the
    // halo's magnetic Reynolds number Rm = a v_A,theta / eta must be << 1, so the field slips
    // through it faster than the column moves. 2 m^2/s gives Rm = 0.06 at q_a = 0.7; the kink then
    // grows as with 20 m^2/s (Rm = 0.006) to 1 %, while 0.2 (Rm = 0.6, the old value) held the
    // lines and slowed it by a fifth (docs/06-mhd-plasma.md, the tokamak section).
    float vacuumResistivity = 2.0f;
    // The helical m = 1, n = 1 offset of the current channel at the start, as a fraction of a: the
    // seed of the kink (a perfectly symmetric start would stay symmetric to rounding for a while).
    float seedDisplacement = 0.02f;

    float plasmaCurrent() const;          // I_p [A] from q_a: 2 pi a^2 B0 / (mu0 R0 q_a)
    float poloidalField(float r) const;   // B_theta(r) [T]
    float safetyFactor(float r) const;    // q(r): q_a inside the channel, q_a r^2 / a^2 outside
    float equilibriumShift() const;       // Delta [m]: where the ring sits in the shell without B_v
    float verticalFieldStrength() const;  // B_v [T] that centres the ring, upwards for a current along +phi
    // B_v [T] that holds a ring with no shell at all (Shafranov 1966; Wesson, "Tokamaks", eq. 3.9.2):
    //   B_v = mu0 I_p / (4 pi R0) [ln(8 R0 / a) + beta_p + l_i / 2 - 3/2]
    // - the most the vertical field ever needs; the position control's integral stays below it.
    float freeRingVerticalField() const;
    float internalInductance() const { return 0.5f; } // l_i of the constant profile
    float wallLimit() const;              // 2 a^2 / (a^2 + b^2): below it the wall holds the kink
    bool kinkUnstable() const { return safetyFactorEdge > wallLimit() && safetyFactorEdge < 1.0f; }
    float kinkGrowthRate(float density) const; // gamma [1/s] of the ideal m = 1 kink, 0 when stable
    // The same kink with the gap between the column and the wall filled with gas of density
    // outerDensityRatio x `density` (the column's). The column shifts rigidly and pushes that gas
    // around it: 2D potential flow between coaxial cylinders adds the mass per unit length
    //   m_added = rho_out pi a^2 (b^2 + a^2) / (b^2 - a^2)   (Lamb, "Hydrodynamics", 1932;
    //   Brennen 1982, "A review of added mass and fluid inertial forces", NCEL CR 82.010),
    // while the energy principle's delta W (Tokamak.h, the head) does not change, so
    //   gamma^2 = gamma_vac^2 / (1 + (rho_out / rho_in) (b^2 + a^2) / (b^2 - a^2)),
    //   gamma_vac^2 = 2 (v_A,theta / a)^2 (1 - q_a) (q_a (b^2 + a^2) - 2 a^2) / (b^2 - a^2),
    // gamma_vac being the textbook kink with a massless vacuum (ratio 0). Ratio 1 - the same gas
    // everywhere, as in the scene - is kinkGrowthRate: the added mass is already in it (5/3 of
    // the column's at b = 2a, 2.67x the inertia, gamma 0.61 gamma_vac). Assumptions: the long
    // wavelength limit (n = 1 along the torus, k a = a / R0 << 1), a rigid m = 1 shift of the
    // constant-current column, incompressible flow. At k a = 0.3 (this device) the 3D flow
    // around the helix (Bessel I1, K1) adds only 1.44 instead of 1.67 column masses: gamma is
    // 4.5 % higher than this long-wavelength value.
    float kinkGrowthRateWithOuterGas(float density, float outerDensityRatio) const;
    bool inside(const Vector3& x) const;  // inside the vessel

    // Position control: the gains of the PID law
    //   B_v = B_v(equilibrium) + gain shift + damping rate + integral * (time integral of the shift),
    // the control's stiffness 4x the shell's, critically damped (gas of `density`), the integral
    // part four times slower than the loop. The integral removes the offset a PD law leaves when
    // B_v(equilibrium) is off: that formula is the large-aspect-ratio one, and here R0 / b = 1.67
    // (as the radial position control of real machines: Ariola & Pironti, "Magnetic Control of
    // Tokamak Plasmas", 2nd ed. 2016, ch. 5).
    void controlGains(float density, float& gain, float& damping, float& integral) const;

    // Vector potentials [T m] of the fields: the coils (toroidal + vertical, current-free: the
    // background B0 of MagneticField) and the plasma current (the evolving B1 at the start).
    Vector3 coilPotential(const Vector3& x) const { return coilPotential(x, verticalFieldStrength()); }
    Vector3 coilPotential(const Vector3& x, float verticalFieldT) const;
    Vector3 plasmaPotential(const Vector3& x) const;
    // A_phi [T m] of one circular current loop (radius Rl at height yl, current I) at (R, y).
    static double loopPotential(double R, double y, double Rl, double yl, double I);
    // The glowing plasma: 1 inside the current channel, 0 outside (a soft edge of a few cm).
    float tracer(const Vector3& x) const;
    // Resistivity [m^2/s] at a point: that of the plasma inside the channel, of the "vacuum" outside.
    float resistivityAt(const Vector3& x, float plasmaResistivity) const;

    // Diagnostics on the grid (dx = its spacing):
    // the plasma current [A] by Ampere's law, the circulation of B around the channel in the
    // poloidal plane phi = 0 (a Rogowski coil);
    float measuredCurrent(const MagneticField& m, float dx) const;
    // the outward shift of the current channel's centroid from the magnetic axis [m], averaged
    // around the torus (the n = 0 part: the equilibrium position in the shell);
    float measuredShift(const MagneticField& m, float dx) const;
    // the kink: amplitude of the n = 1 harmonic, around the torus, of the current channel's
    // centroid displacement from the magnetic axis [m] - the seed for a symmetric ring, growing
    // when the column winds into a helix.
    float kinkAmplitude(const MagneticField& m, float dx) const;
    // The same n = 1 amplitude for the plasma itself: the centroid of the tracer (the gas's smoke
    // field, carried by the flow; gridOrigin and dx of the gas grid). The current can lag the
    // moving plasma (it diffuses where the resistivity ramps up), so the two differ; the plasma's
    // is the displacement xi of the energy principle.
    float plasmaKinkAmplitude(const Field3& tracer, const Vector3& gridOrigin, float dx) const;

private:
    // Centroid (dR, dy) of the toroidal current in N poloidal planes around the torus.
    void currentCentroids(const MagneticField& m, float dx, std::vector<Vector2>& out) const;
    // Centroid (dR, dy) of the tracer in the same planes.
    void tracerCentroids(const Field3& tracer, const Vector3& gridOrigin, float dx, std::vector<Vector2>& out) const;
    // The n = 1 Fourier amplitude of the centroids' displacement dR + i dy around the torus.
    static float helicalAmplitude(const std::vector<Vector2>& centroids);
    // The current channel's centre in the poloidal plane at toroidal angle phi (with the seed).
    void channelCentre(float phi, float& R, float& y) const;
    float channelRadius(const Vector3& x) const; // r from the (seeded) channel centre
};

} // namespace rf
