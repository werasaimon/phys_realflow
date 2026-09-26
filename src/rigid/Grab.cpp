// The mouse joint of RigidWorld (Box2D's b2MouseJoint): a soft spring from a point on a body
// to a target, solved as a constraint with the contacts.
#include "rigid/RigidWorld.h"

#include <algorithm>

namespace rf {

void RigidWorld::grab(int i, const Vector3& p) {
    if (i < 0 || i >= int(bodies_.size()) || bodies_[i].invMass == 0) return;
    grab_ = GrabJoint();
    grab_.active = true;
    grab_.body = i;
    grab_.localAnchor = bodies_[i].rotation().transposed() * (p - bodies_[i].pos);
    grab_.target = p;
    wake(i);
}

// Soft constraint (Box2D mouse joint): spring-damper with the given frequency and damping ratio,
// expressed through gamma (softness) and beta (position feedback), solved as a 3D point constraint.
void RigidWorld::prepareGrab(float dt) {
    if (!grab_.active) return;
    if (grab_.body >= int(bodies_.size())) { releaseGrab(); return; }
    RigidBody& b = bodies_[grab_.body];
    b.sleepTimer = 0; // a held body never falls asleep
    if (b.invMass == 0) return;
    const float omega = 2.0f * kPi * grab_.frequency;
    const float d = 2.0f * b.mass * grab_.damping * omega;
    const float k = b.mass * omega * omega;
    grabGamma_ = dt * (d + dt * k);
    grabGamma_ = grabGamma_ > 0 ? 1.0f / grabGamma_ : 0.0f;
    const float beta = dt * k * grabGamma_;
    Vector3 r = b.rotation() * grab_.localAnchor;
    Matrix3x3 S = Matrix3x3::skew(r);
    Matrix3x3 K = S * b.invInertiaWorld * S.transposed();
    for (int i = 0; i < 3; ++i) K.m[i][i] += b.invMass + grabGamma_;
    grabMass_ = K.inverse();
    grabBias_ = (b.pos + r - grab_.target) * beta;
    // Warm start and extra angular damping so the held body does not spin forever.
    b.vel += grab_.impulse * b.invMass;
    b.angVel += b.applyInvInertiaWorld(cross(r, grab_.impulse));
    b.angVel *= std::max(0.0f, 1.0f - 2.0f * dt);
}

void RigidWorld::solveGrab(float dt) {
    if (!grab_.active) return;
    RigidBody& b = bodies_[grab_.body];
    if (b.invMass == 0) return;
    Vector3 r = b.rotation() * grab_.localAnchor;
    Vector3 cdot = b.vel + cross(b.angVel, r);
    Vector3 imp = grabMass_ * (-(cdot + grabBias_ + grab_.impulse * grabGamma_));
    Vector3 old = grab_.impulse;
    grab_.impulse += imp;
    float maxImpulse = grab_.maxForce * b.mass * 9.81f * dt; // weights at standard gravity, whatever params.gravity is
    float l = length(grab_.impulse);
    if (l > maxImpulse) grab_.impulse *= maxImpulse / l;
    imp = grab_.impulse - old;
    b.vel += imp * b.invMass;
    b.angVel += b.applyInvInertiaWorld(cross(r, imp));
}

} // namespace rf
