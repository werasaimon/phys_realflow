// XPBD rigid-body solver: M. Mueller, M. Macklin, N. Chentanez, S. Jeschke, T.-Y. Kim,
// "Detailed Rigid Body Simulation with Extended Position Based Dynamics", SCA 2020.
//
// Per substep h:  integrate -> solve positions (contacts, static friction) -> derive velocities
//                 -> solve velocities (dynamic friction, restitution, resting damping)
// Contacts keep their anchor points in the bodies' local frames, so the penetration depth is
// re-evaluated every substep while collision detection runs only every few substeps.

#include "rigid/RigidWorld.h"

#include "core/Probe.h"

#include <algorithm>

namespace rf {

Vector3 RigidWorld::anchorA(const XContact& c, bool prev) const {
    const RigidBody& A = bodies_[c.a];
    return prev ? A.prevPos + A.prevRot.rotate(c.rA) : A.pos + A.rot.rotate(c.rA);
}

Vector3 RigidWorld::anchorB(const XContact& c, bool prev) const {
    if (c.b < 0) return c.rB;
    const RigidBody& B = bodies_[c.b];
    return prev ? B.prevPos + B.prevRot.rotate(c.rB) : B.pos + B.rot.rotate(c.rB);
}

void RigidWorld::buildXContacts() {
    xcontacts_.clear();
    for (const Manifold& m : manifolds_) {
        const RigidBody& A = bodies_[m.a];
        const RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
        for (const SolverPoint& p : m.points) {
            XContact c;
            c.a = m.a;
            c.b = m.b;
            c.n = p.normal;
            // Surface points: A's point lies on the -n side of the midpoint, B's on the +n side.
            Vector3 pA = p.position - p.normal * (0.5f * p.depth);
            Vector3 pB = p.position + p.normal * (0.5f * p.depth);
            c.rA = A.rot.conjugate().rotate(pA - A.pos);
            c.rB = B ? B->rot.conjugate().rotate(pB - B->pos) : pB;
            c.friction = m.friction;
            c.staticFriction = m.staticFriction;
            c.restitution = m.restitution;
            xcontacts_.push_back(c);
        }
    }
    // Ground-up order (Gauss-Seidel): support reaches the top of a stack within the substep.
    Vector3 up = normalize(-params.gravity);
    if (length2(up) < 0.5f) up = Vector3(0, 1, 0);
    auto height = [&](const XContact& c) {
        float h = dot(bodies_[c.a].pos, up);
        return c.b >= 0 ? std::min(h, dot(bodies_[c.b].pos, up)) : h - 1e3f;
    };
    std::stable_sort(xcontacts_.begin(), xcontacts_.end(),
                     [&](const XContact& x, const XContact& y) { return height(x) < height(y); });
    contactCount_ = xcontacts_.size();
}

// Moves A by +corr and B by -corr at the given points, distributed by generalised inverse masses
// w = 1/m + (r x n)^T I^-1 (r x n). Returns the positional multiplier |lambda|.
float RigidWorld::applyPositionalCorrection(int a, int b, const Vector3& corr, const Vector3& pA, const Vector3& pB) {
    float mag = length(corr);
    if (mag < 1e-12f) return 0;
    Vector3 n = corr / mag;
    RigidBody& A = bodies_[a];
    RigidBody* B = b >= 0 ? &bodies_[b] : nullptr;
    Vector3 rA = pA - A.pos, rB = B ? pB - B->pos : Vector3(0.0f);
    Vector3 cA = cross(rA, n), cB = cross(rB, n);
    float w1 = A.invMass > 0 ? A.invMass + dot(cA, A.applyInvInertiaWorld(cA)) : 0.0f;
    float w2 = (B && B->invMass > 0) ? B->invMass + dot(cB, B->applyInvInertiaWorld(cB)) : 0.0f;
    if (w1 + w2 < 1e-12f) return 0;
    float lambda = mag / (w1 + w2);
    Vector3 P = n * lambda;
    if (A.invMass > 0) {
        A.pos += P * A.invMass;
        A.rot = A.rot.integrated(A.applyInvInertiaWorld(cross(rA, P)), 1.0f);
    }
    if (B && B->invMass > 0) {
        B->pos -= P * B->invMass;
        B->rot = B->rot.integrated(B->applyInvInertiaWorld(cross(rB, -P)), 1.0f);
    }
    return lambda;
}

void RigidWorld::applyVelocityChange(int a, int b, const Vector3& dv, const Vector3& pA, const Vector3& pB) {
    float mag = length(dv);
    if (mag < 1e-12f) return;
    Vector3 n = dv / mag;
    RigidBody& A = bodies_[a];
    RigidBody* B = b >= 0 ? &bodies_[b] : nullptr;
    Vector3 rA = pA - A.pos, rB = B ? pB - B->pos : Vector3(0.0f);
    Vector3 cA = cross(rA, n), cB = cross(rB, n);
    float w1 = A.invMass > 0 ? A.invMass + dot(cA, A.applyInvInertiaWorld(cA)) : 0.0f;
    float w2 = (B && B->invMass > 0) ? B->invMass + dot(cB, B->applyInvInertiaWorld(cB)) : 0.0f;
    if (w1 + w2 < 1e-12f) return;
    Vector3 P = n * (mag / (w1 + w2));
    if (A.invMass > 0) {
        A.vel += P * A.invMass;
        A.angVel += A.applyInvInertiaWorld(cross(rA, P));
    }
    if (B && B->invMass > 0) {
        B->vel -= P * B->invMass;
        B->angVel -= B->applyInvInertiaWorld(cross(rB, P));
    }
}

void RigidWorld::solveXContactPosition(XContact& c, float h) {
    Vector3 pA = anchorA(c), pB = anchorB(c);
    float d = dot(pB - pA, c.n); // penetration along the normal (> 0: overlapping)
    if (d <= 0) return;
    c.active = true;
    c.lambdaN += applyPositionalCorrection(c.a, c.b, c.n * d, pA, pB);

    // Static friction as a position constraint: undo the tangential slip of the contact points
    // during this substep if the required multiplier stays inside the static friction cone.
    pA = anchorA(c);
    pB = anchorB(c);
    Vector3 dp = (pA - anchorA(c, true)) - (pB - anchorB(c, true));
    Vector3 dpt = dp - c.n * dot(dp, c.n);
    float slip = length(dpt);
    if (slip < 1e-9f) return;
    // Multiplier the correction would need (evaluated without applying it).
    RigidBody& A = bodies_[c.a];
    RigidBody* B = c.b >= 0 ? &bodies_[c.b] : nullptr;
    Vector3 t = dpt / slip;
    Vector3 rA = pA - A.pos, rB = B ? pB - B->pos : Vector3(0.0f);
    Vector3 cA = cross(rA, t), cB = cross(rB, t);
    float w = (A.invMass > 0 ? A.invMass + dot(cA, A.applyInvInertiaWorld(cA)) : 0.0f) +
              ((B && B->invMass > 0) ? B->invMass + dot(cB, B->applyInvInertiaWorld(cB)) : 0.0f);
    if (w < 1e-12f) return;
    float lambdaT = slip / w;
    if (c.lambdaT + lambdaT < c.staticFriction * c.lambdaN) {
        c.lambdaT += applyPositionalCorrection(c.a, c.b, -dpt, pA, pB);
    }
    (void)h;
}

void RigidWorld::solveXContactVelocity(XContact& c, float h) {
    if (!c.active) return;
    RigidBody& A = bodies_[c.a];
    RigidBody* B = c.b >= 0 ? &bodies_[c.b] : nullptr;
    Vector3 pA = anchorA(c), pB = anchorB(c);
    Vector3 v = A.velocityAt(pA) - (B ? B->velocityAt(pB) : Vector3(0.0f));
    float vn = dot(v, c.n);
    Vector3 vt = v - c.n * vn;
    float vtl = length(vt);
    Vector3 dv(0.0f);
    // Dynamic (kinetic) friction: |dv| <= h * mu_k * |f_n|, with f_n = lambda_n / h^2.
    if (vtl > 1e-9f) {
        float fn = c.lambdaN / (h * h);
        dv -= vt * (std::min(h * c.friction * fn, vtl) / vtl);
    }
    // Restitution (switched off for slow contacts so resting bodies do not jitter).
    float e = std::fabs(c.vnPrev) <= 2.0f * length(params.gravity) * h ? 0.0f : c.restitution;
    float target = std::max(-e * c.vnPrev, 0.0f);
    if (vn < target) dv += c.n * (target - vn);
    applyVelocityChange(c.a, c.b, dv, pA, pB);
}

// One XPBD substep (Macklin, Müller & Chentanez 2016; Müller et al. 2020 "Detailed rigid body
// simulation with XPBD"): predict the poses, correct them by the contact and joint constraints,
// derive the velocities from the corrected poses, then friction, restitution and damping at
// velocity level.
void RigidWorld::stepXPBD(float h) {
    lastDt_ = h;
    {
        Probe::Timer t("rigid/collide ms");
        xpbdDetectContacts(h);
    }
    {
        Probe::Timer t("rigid/integrate ms");
        xpbdIntegrate(h);
    }
    Probe::Timer t("rigid/solve ms"); // the rest of the substep: the constraints and the velocities they leave
    // 2) Positions: contacts with static friction.
    for (int it = 0; it < params.positionIterations; ++it) {
        for (XContact& c : xcontacts_) solveXContactPosition(c, h);
        for (RigidBody& b : bodies_) b.updateInertia();
    }
    solveJointPositions();
    xpbdVelocitiesFromPoses(h);
    // 4) Velocity level: dynamic friction, restitution, mouse joint.
    for (XContact& c : xcontacts_) solveXContactVelocity(c, h);
    prepareGrab(h);
    for (auto& j : joints_) j->prepare(bodies_, h, params.warmStarting);
    for (int it = 0; it < 4; ++it) {
        solveGrab(h);
        for (auto& j : joints_) j->solveVelocity(bodies_);
    }
    xpbdDamp(h);
}

// Collision detection every few substeps; the margin covers the motion until the next one, so a
// contact found now is still valid when the bodies have moved on.
void RigidWorld::xpbdDetectContacts(float h) {
    if (substepCounter_++ % std::max(1, params.collisionInterval) != 0 && !xcontacts_.empty()) return;
    float vmax = 0;
    for (const RigidBody& b : bodies_)
        if (b.invMass > 0) vmax = std::max(vmax, length(b.vel) + length(b.angVel) * b.boundingRadius());
    float saved = params.contactMargin;
    params.contactMargin = std::max(saved, 1.5f * vmax * h * std::max(1, params.collisionInterval));
    for (RigidBody& b : bodies_) b.updateInertia();
    collide();
    params.contactMargin = saved;
    buildXContacts();
}

// 1) Integrate (explicit, gyroscopic term included): the predicted poses the constraints will
// correct; every contact remembers its normal velocity before the substep, for the restitution.
void RigidWorld::xpbdIntegrate(float h) {
    for (RigidBody& b : bodies_) {
        b.prevPos = b.pos;
        b.prevRot = b.rot;
        if (b.invMass == 0) continue;
        b.vel += (params.gravity + b.force * b.invMass) * h;
        Matrix3x3 R = b.rotation();
        Vector3 wl = R.transposed() * b.angVel;
        Vector3 I(1.0f / b.invInertiaLocal.x, 1.0f / b.invInertiaLocal.y, 1.0f / b.invInertiaLocal.z);
        Vector3 gyro = R * cross(wl, I * wl);
        b.angVel += b.applyInvInertiaWorld(b.torque - gyro) * h;
        b.pos += b.vel * h;
        b.rot = b.rot.integrated(b.angVel, h);
        b.updateInertia();
    }
    // Normal velocity before the substep, for restitution.
    for (XContact& c : xcontacts_) {
        c.lambdaN = c.lambdaT = 0;
        c.active = false;
        const RigidBody& A = bodies_[c.a];
        const RigidBody* B = c.b >= 0 ? &bodies_[c.b] : nullptr;
        Vector3 pA = anchorA(c, true), pB = anchorB(c, true);
        c.vnPrev = dot(A.velocityAt(pA) - (B ? B->velocityAt(pB) : Vector3(0.0f)), c.n);
    }
}

// 3) Velocities from the position change: this is what makes XPBD stable - the velocity is
// whatever the constraints allowed the body to move.
void RigidWorld::xpbdVelocitiesFromPoses(float h) {
    for (RigidBody& b : bodies_) {
        if (b.invMass == 0) continue;
        b.vel = (b.pos - b.prevPos) / h;
        Quaternion dq = b.rot * b.prevRot.conjugate();
        Vector3 w(dq.x, dq.y, dq.z);
        b.angVel = w * (2.0f / h) * (dq.w >= 0 ? 1.0f : -1.0f);
        b.updateInertia();
    }
}

// Damping: global, plus the capped anti-phase damping of slow bodies in contact; the forces of
// this substep are cleared.
void RigidWorld::xpbdDamp(float h) {
    const float ld = std::max(0.0f, 1.0f - params.linearDamping * h);
    const float ad = std::max(0.0f, 1.0f - params.angularDamping * h);
    std::vector<char> touching(bodies_.size(), 0);
    for (const XContact& c : xcontacts_)
        if (c.active) {
            touching[c.a] = 1;
            if (c.b >= 0) touching[c.b] = 1;
        }
    const float cap = params.restDamping * h;
    for (size_t i = 0; i < bodies_.size(); ++i) {
        RigidBody& b = bodies_[i];
        b.force = Vector3(0.0f);
        b.torque = Vector3(0.0f);
        if (b.invMass == 0) continue;
        b.vel *= ld;
        b.angVel *= ad;
        if (touching[i] && params.restDamping > 0) {
            float v = length(b.vel), w = length(b.angVel);
            if (v < params.restLinearThreshold) b.vel *= v > cap ? 1.0f - cap / v : 0.0f;
            float wcap = cap / std::max(b.boundingRadius(), 1e-3f);
            if (w < params.restAngularThreshold) b.angVel *= w > wcap ? 1.0f - wcap / w : 0.0f;
        }
    }
}

} // namespace rf
