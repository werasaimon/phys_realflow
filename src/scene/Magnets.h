#pragma once
// Permanent magnets on rigid bodies, as point dipoles: every magnet makes the field of a dipole,
// and every other magnet feels a force and a torque in it. This is the physics behind the
// editor's "magnet" role (SceneGraph.h): two bar magnets pull head to tail, push side by side, and
// a compass needle turns along the field.
//
// Formulas (SI units, mu0 = 4 pi 1e-7 T m / A):
//   field of dipole m at offset r:   B = mu0 / (4 pi |r|^3) (3 r_hat (m . r_hat) - m)
//                                    (Jackson, Classical Electrodynamics, 3rd ed., eq. 5.56)
//   force on dipole 2 from dipole 1: F = 3 mu0 / (4 pi |r|^4) [ (m1 . r_hat) m2 + (m2 . r_hat) m1
//                                        + (m1 . m2) r_hat - 5 (m1 . r_hat)(m2 . r_hat) r_hat ],
//                                    r = x2 - x1 (Yung, Landecker & Villani 1998, Magnetic and
//                                    Electrical Separation 9, "An analytic solution for the force
//                                    between two magnetic dipoles")
//   torque on dipole 2:              tau = m2 x B1(r)
// The force pair is equal and opposite (Newton's third law holds exactly, the test checks it).
#include "math/Math.h"

#include <vector>

namespace rf {

class RigidWorld;

// The field [T] of a dipole with moment m [A m^2] at the offset r from it.
Vector3 dipoleField(const Vector3& m, const Vector3& r);
// The force [N] on dipole m2 from dipole m1, r = position of m2 - position of m1.
Vector3 dipoleForce(const Vector3& m1, const Vector3& m2, const Vector3& r);

// One frame of magnet forces on rigid bodies: every pair pushes or pulls, every magnet is turned
// by the field of the others. `momentsBody` are the moments in each body's own frame (they turn
// with the body). Applied as impulses force * dt and torque * dt. Returns the largest force [N].
float applyMagnetForces(RigidWorld& world, const std::vector<int>& bodies, const std::vector<Vector3>& momentsBody, float dt);

} // namespace rf
