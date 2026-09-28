// The turn of a rigid body over one step: how its orientation and its spin move when nothing pushes
// it - Euler's equations of the torque-free body - by a symplectic splitting (McLachlan 1993,
// "Explicit Lie-Poisson integration and the Euler equations", Phys. Rev. Lett. 71, 3043; Dullweber,
// Leimkuhler & McLachlan 1997, "Symplectic splitting methods for rigid body molecular dynamics",
// J. Chem. Phys. 107, 5840). The pushes - torques, contact impulses - have already changed the spin
// before this drift (a kick), so a step is "kick, then turn freely", as in the symplectic Euler of
// the positions.
//
// Why a splitting. In the body's principal frame the energy of the turn is a sum of three parts,
//
//     T = P1^2 / 2 I1  +  P2^2 / 2 I2  +  P3^2 / 2 I3,        P = I w  (angular momentum, body frame)
//
// and each part alone is a motion we can do EXACTLY: a spin about one principal axis at the
// constant rate Pi / Ii. Over a time t the body turns by the angle theta = t Pi / Ii about that
// axis, and the momentum, seen from the turning body, turns back by -theta about the same axis
// (its own component Pi does not change, so theta is exact, not an estimate). The whole motion is
// these exact turns one after the other, in the symmetric order
//
//     R1(h/2)  R2(h/2)  R3(h)  R2(h/2)  R1(h/2)
//
// which makes the method second order and symplectic: the energy error stays bounded forever
// instead of drifting, and because every sub-step is a rotation, |P| and the world angular momentum
// L = R P are kept to rounding. The implicit midpoint rule it replaces kept |P| and T for the spin
// but advanced the orientation with the spin of the end of the step - first order in the time
// step: the world L drifted, and a tennis racket flipped 4 times in 20 s instead of 3.
//
// In double precision, inside: a float round trip (spin -> momentum -> five turns -> spin) every
// step added its rounding, ~1e-7, to |L| each time - a drift that grew as the step shrank. The
// body keeps float; only this arithmetic is double.
//
// A body with three equal moments (a cube, a sphere) spins about a fixed axis at a constant rate:
// that motion is one exact rotation, done in one step as before.
#include "rigid/RigidWorld.h"

#include <cmath>

namespace rf {

namespace {

// A unit quaternion in double: just what the free turn needs.
struct TurnD {
    double w = 1, x = 0, y = 0, z = 0;

    TurnD operator*(const TurnD& q) const { // this rotation after q
        return {w * q.w - x * q.x - y * q.y - z * q.z, w * q.x + x * q.w + y * q.z - z * q.y,
                w * q.y - x * q.z + y * q.w + z * q.x, w * q.z + x * q.y - y * q.x + z * q.w};
    }
    // v turned by this rotation: v + 2 u x (u x v + w v), u = (x, y, z).
    void rotate(const double v[3], double out[3]) const {
        const double t[3] = {y * v[2] - z * v[1] + w * v[0], z * v[0] - x * v[2] + w * v[1], x * v[1] - y * v[0] + w * v[2]};
        out[0] = v[0] + 2 * (y * t[2] - z * t[1]);
        out[1] = v[1] + 2 * (z * t[0] - x * t[2]);
        out[2] = v[2] + 2 * (x * t[1] - y * t[0]);
    }
    TurnD conjugate() const { return {w, -x, -y, -z}; }
};

// One exact sub-turn about principal axis `axis` for the time `tau`: the body turns by
// theta = tau Pi / Ii about its own axis (a rotation on the right of `rot`: in the body frame), and
// P turns by -theta about the same axis, so that R P - the world angular momentum - stays as it was.
void turnAboutAxis(TurnD& rot, double P[3], const double inverseMoments[3], int axis, double tau) {
    const double theta = tau * P[axis] * inverseMoments[axis];
    const double c = rf::cos(theta), s = rf::sin(theta);
    TurnD step{rf::cos(0.5 * theta), 0, 0, 0};
    (axis == 0 ? step.x : axis == 1 ? step.y : step.z) = rf::sin(0.5 * theta);
    rot = rot * step;
    // The other two components (j, k in cyclic order after `axis`) turn by -theta.
    const int j = (axis + 1) % 3, k = (axis + 2) % 3;
    const double Pj = P[j], Pk = P[k];
    P[j] = c * Pj + s * Pk;
    P[k] = -s * Pj + c * Pk;
}

} // namespace

// The drift of one body over the step h (see the file head):
//   1. a body with three equal moments: one exact turn about its spin axis, as before;
//   2. otherwise the spin as angular momentum in the principal frame, P = I (R^-1 w), in double;
//   3. the five exact sub-turns R1(h/2) R2(h/2) R3(h) R2(h/2) R1(h/2);
//   4. the spin back in the world, w = R I^-1 P, and the orientation back in float;
//   5. the pseudo spin of the split impulse (penetration recovery), a small turn of its own that
//      moves the pose but is no part of the motion.
void RigidWorld::turnFreely(RigidBody& b, float h) {
    const Vector3& inv = b.invInertiaLocal;
    if (inv.x == inv.y && inv.y == inv.z) { // 1.
        b.rot = b.rot.integrated(b.angVel + b.biasAngVel, h);
        return;
    }
    const double inverseMoments[3] = {inv.x, inv.y, inv.z};
    TurnD rot{b.rot.w, b.rot.x, b.rot.y, b.rot.z}; // 2.
    const double norm = std::sqrt(rot.w * rot.w + rot.x * rot.x + rot.y * rot.y + rot.z * rot.z);
    rot = {rot.w / norm, rot.x / norm, rot.y / norm, rot.z / norm};
    const double worldSpin[3] = {b.angVel.x, b.angVel.y, b.angVel.z};
    double bodySpin[3], P[3];
    rot.conjugate().rotate(worldSpin, bodySpin);
    for (int a = 0; a < 3; ++a) P[a] = bodySpin[a] / inverseMoments[a];
    turnAboutAxis(rot, P, inverseMoments, 0, 0.5 * h); // 3.
    turnAboutAxis(rot, P, inverseMoments, 1, 0.5 * h);
    turnAboutAxis(rot, P, inverseMoments, 2, h);
    turnAboutAxis(rot, P, inverseMoments, 1, 0.5 * h);
    turnAboutAxis(rot, P, inverseMoments, 0, 0.5 * h);
    for (int a = 0; a < 3; ++a) bodySpin[a] = P[a] * inverseMoments[a]; // 4.
    double spin[3];
    rot.rotate(bodySpin, spin);
    b.angVel = Vector3(float(spin[0]), float(spin[1]), float(spin[2]));
    b.rot = Quaternion{float(rot.w), float(rot.x), float(rot.y), float(rot.z)}.normalized();
    if (length2(b.biasAngVel) > 0) b.rot = b.rot.integrated(b.biasAngVel, h); // 5.
}

} // namespace rf
