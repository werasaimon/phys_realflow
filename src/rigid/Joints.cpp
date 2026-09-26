#include "rigid/Joints.h"

#include <algorithm>

namespace rf {

namespace {

// The static world as a body: infinite mass, identity pose. Never written to: applyRowImpulse and
// moveBody skip bodies with invMass == 0 (so a bad impulse cannot poison this shared object).
RigidBody& worldBody() {
    static RigidBody w = [] {
        RigidBody b;
        b.mass = 0;
        b.invMass = 0;
        b.invInertiaLocal = Vector3(0.0f);
        b.invInertiaWorld = Matrix3x3::zero();
        return b;
    }();
    return w;
}

RigidBody& bodyOf(std::vector<RigidBody>& bodies, int i) { return i >= 0 ? bodies[i] : worldBody(); }
const RigidBody& bodyOf(const std::vector<RigidBody>& bodies, int i) { return i >= 0 ? bodies[i] : worldBody(); }

// Velocity change of a body by the impulse lambda of one Jacobian row (its linear and angular part).
void applyRowImpulse(RigidBody& b, const Vector3& lin, const Vector3& ang, float lambda) {
    if (b.invMass == 0) return;
    b.vel += lin * (b.invMass * lambda);
    b.angVel += b.applyInvInertiaWorld(ang) * lambda;
}

// Position-stage step limits (keep the nonlinear Gauss-Seidel robust for large errors).
constexpr float kMaxLinearCorrection = 0.2f;  // m
constexpr float kMaxAngularCorrection = 0.5f; // rad

void moveBody(RigidBody& b, const Vector3& dp, const Vector3& dRot) {
    if (b.invMass == 0) return;
    b.pos += dp;
    b.rot = b.rot.integrated(dRot, 1.0f);
    b.updateInertia();
}

} // namespace

// ---------------------------------------------------------------------------
// Generic Jacobian-row machinery
// ---------------------------------------------------------------------------
Vector3 Joint::worldAnchorA(const std::vector<RigidBody>& bodies) const {
    const RigidBody& A = bodyOf(bodies, a);
    return A.pos + A.rot.rotate(localAnchorA);
}

Vector3 Joint::worldAnchorB(const std::vector<RigidBody>& bodies) const {
    const RigidBody& B = bodyOf(bodies, b);
    return B.pos + B.rot.rotate(localAnchorB);
}

Vector3 Joint::worldAxis(const std::vector<RigidBody>& bodies) const {
    return bodyOf(bodies, a).rot.rotate(localAxisA);
}

float Joint::appliedImpulse() const {
    float s = 0;
    for (const JacobianRow& r : rows_) s += r.lambda * r.lambda;
    return std::sqrt(s);
}

void Joint::addRow(const std::vector<RigidBody>& bodies, JacobianRow r) {
    const RigidBody& A = bodyOf(bodies, a);
    const RigidBody& B = bodyOf(bodies, b);
    float k = A.invMass * dot(r.linA, r.linA) + dot(r.angA, A.applyInvInertiaWorld(r.angA)) +
              B.invMass * dot(r.linB, r.linB) + dot(r.angB, B.applyInvInertiaWorld(r.angB)) + r.softness;
    r.effMass = k > 1e-12f ? 1.0f / k : 0.0f;
    rows_.push_back(r);
}

void Joint::prepare(std::vector<RigidBody>& bodies, float h, bool warmStart) {
    rows_.clear();
    buildRows(bodies, h);
    RigidBody& A = bodyOf(bodies, a);
    RigidBody& B = bodyOf(bodies, b);
    const bool warm = warmStart && warm_.size() == rows_.size();
    for (size_t i = 0; i < rows_.size(); ++i) {
        JacobianRow& r = rows_[i];
        r.lambda = warm ? clampv(warm_[i], r.lo, r.hi) : 0.0f;
        if (r.lambda == 0) continue;
        applyRowImpulse(A, r.linA, r.angA, r.lambda);
        applyRowImpulse(B, r.linB, r.angB, r.lambda);
    }
}

void Joint::solveVelocity(std::vector<RigidBody>& bodies) {
    RigidBody& A = bodyOf(bodies, a);
    RigidBody& B = bodyOf(bodies, b);
    for (JacobianRow& r : rows_) {
        float jv = dot(r.linA, A.vel) + dot(r.angA, A.angVel) + dot(r.linB, B.vel) + dot(r.angB, B.angVel);
        float d = -r.effMass * (jv + r.bias + r.softness * r.lambda);
        float old = r.lambda;
        r.lambda = clampv(old + d, r.lo, r.hi);
        d = r.lambda - old;
        if (d == 0) continue;
        applyRowImpulse(A, r.linA, r.angA, d);
        applyRowImpulse(B, r.linB, r.angB, d);
    }
    warm_.resize(rows_.size());
    for (size_t i = 0; i < rows_.size(); ++i) warm_[i] = rows_[i].lambda;
}

// Position stage helpers ------------------------------------------------------------------------
// Point constraint C = pA - pB: lambda = -K^-1 C with K = sum(1/m I - [r]x I^-1 [r]x).
float Joint::correctPoint(std::vector<RigidBody>& bodies, const Vector3& C, const Vector3& pA, const Vector3& pB, float beta) {
    RigidBody& A = bodyOf(bodies, a);
    RigidBody& B = bodyOf(bodies, b);
    float err = length(C);
    if (err < 1e-7f) return err;
    Vector3 rA = pA - A.pos, rB = pB - B.pos;
    Matrix3x3 SA = Matrix3x3::skew(rA), SB = Matrix3x3::skew(rB);
    Matrix3x3 K = Matrix3x3::diag(Vector3(A.invMass + B.invMass)) + SA * A.invInertiaWorld * SA.transposed() +
                  SB * B.invInertiaWorld * SB.transposed();
    Matrix3x3 Ki = K.inverse(); // zero for two static bodies: no correction
    Vector3 corr = C * beta;
    if (length(corr) > kMaxLinearCorrection) corr *= kMaxLinearCorrection / length(corr);
    Vector3 P = Ki * (-corr);
    moveBody(A, P * A.invMass, A.applyInvInertiaWorld(cross(rA, P)));
    moveBody(B, -P * B.invMass, B.applyInvInertiaWorld(cross(rB, -P)));
    return err;
}

// Angular constraint with rotation error e (A relative to B): L = -(IA^-1 + IB^-1)^-1 e.
float Joint::correctAngle(std::vector<RigidBody>& bodies, Vector3 e, float beta) {
    RigidBody& A = bodyOf(bodies, a);
    RigidBody& B = bodyOf(bodies, b);
    float err = length(e);
    if (err < 1e-7f) return err;
    Matrix3x3 Ki = (A.invInertiaWorld + B.invInertiaWorld).inverse(); // zero if neither can turn
    e *= beta;
    if (length(e) > kMaxAngularCorrection) e *= kMaxAngularCorrection / length(e);
    Vector3 L = Ki * (-e);
    moveBody(A, Vector3(0.0f), A.applyInvInertiaWorld(L));
    moveBody(B, Vector3(0.0f), B.applyInvInertiaWorld(-L));
    return err;
}

// Scalar constraint C along n (C = n.(pA - pB) - target): moves A by +lambda n, B by -lambda n.
float Joint::correctAlong(std::vector<RigidBody>& bodies, const Vector3& n, float C, const Vector3& pA, const Vector3& pB, float beta) {
    RigidBody& A = bodyOf(bodies, a);
    RigidBody& B = bodyOf(bodies, b);
    Vector3 rA = pA - A.pos, rB = pB - B.pos;
    Vector3 cA = cross(rA, n), cB = cross(rB, n);
    float k = A.invMass + dot(cA, A.applyInvInertiaWorld(cA)) + B.invMass + dot(cB, B.applyInvInertiaWorld(cB));
    if (k < 1e-12f) return std::fabs(C);
    float c = clampv(C * beta, -kMaxLinearCorrection, kMaxLinearCorrection);
    Vector3 P = n * (-c / k);
    moveBody(A, P * A.invMass, A.applyInvInertiaWorld(cross(rA, P)));
    moveBody(B, -P * B.invMass, B.applyInvInertiaWorld(cross(rB, -P)));
    return std::fabs(C);
}

// Orientation error on SO(3): E = (qB^-1 qA) refRel^-1, log(E) mapped to world through B.
Vector3 Joint::rotationError(const std::vector<RigidBody>& bodies) const {
    const RigidBody& A = bodyOf(bodies, a);
    const RigidBody& B = bodyOf(bodies, b);
    Quaternion qRel = B.rot.conjugate() * A.rot;
    return B.rot.rotate((qRel * refRel.conjugate()).log());
}

// ---------------------------------------------------------------------------
// Ball
// ---------------------------------------------------------------------------
void BallJoint::buildRows(const std::vector<RigidBody>& bodies, float) {
    const RigidBody& A = bodyOf(bodies, a);
    const RigidBody& B = bodyOf(bodies, b);
    Vector3 rA = A.rot.rotate(localAnchorA), rB = B.rot.rotate(localAnchorB);
    const Vector3 axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (const Vector3& e : axes) {
        JacobianRow r;
        r.linA = e;
        r.angA = cross(rA, e);
        r.linB = -e;
        r.angB = -cross(rB, e);
        addRow(bodies, r);
    }
}

float BallJoint::solvePosition(std::vector<RigidBody>& bodies) {
    Vector3 pA = worldAnchorA(bodies), pB = worldAnchorB(bodies);
    return correctPoint(bodies, pA - pB, pA, pB, 0.5f);
}

// ---------------------------------------------------------------------------
// Hinge
// ---------------------------------------------------------------------------
float HingeJoint::angle(const std::vector<RigidBody>& bodies) const {
    const RigidBody& A = bodyOf(bodies, a);
    const RigidBody& B = bodyOf(bodies, b);
    Vector3 axis = A.rot.rotate(localAxisA);
    Vector3 ra = A.rot.rotate(localRefA), rb = B.rot.rotate(localRefB);
    return std::atan2(dot(cross(rb, ra), axis), dot(rb, ra));
}

void HingeJoint::buildRows(const std::vector<RigidBody>& bodies, float h) {
    BallJoint::buildRows(bodies, h);
    const RigidBody& A = bodyOf(bodies, a);
    Vector3 axis = A.rot.rotate(localAxisA);
    Vector3 t1 = anyPerpendicular(axis), t2 = cross(axis, t1);
    for (const Vector3& t : {t1, t2}) { // relative rotation only about the axis
        JacobianRow r;
        r.angA = t;
        r.angB = -t;
        addRow(bodies, r);
    }
    if (motorEnabled) { // drive (wA - wB).axis to motorSpeed with a bounded torque
        JacobianRow r;
        r.angA = axis;
        r.angB = -axis;
        r.bias = -motorSpeed;
        r.lo = -maxMotorTorque * h;
        r.hi = maxMotorTorque * h;
        addRow(bodies, r);
    }
    if (limitEnabled) {
        float th = angle(bodies);
        if (th <= lower || th >= upper) {
            JacobianRow r;
            r.angA = axis;
            r.angB = -axis;
            if (th <= lower) r.lo = 0;
            else r.hi = 0;
            addRow(bodies, r);
        }
    }
}

float HingeJoint::solvePosition(std::vector<RigidBody>& bodies) {
    float err = BallJoint::solvePosition(bodies);
    const RigidBody& A = bodyOf(bodies, a);
    const RigidBody& B = bodyOf(bodies, b);
    Vector3 aA = A.rot.rotate(localAxisA), aB = B.rot.rotate(localAxisB);
    err += correctAngle(bodies, cross(aB, aA), 0.5f); // tilt of A's axis away from B's axis
    if (limitEnabled) {
        float th = angle(bodies);
        Vector3 axis = bodyOf(bodies, a).rot.rotate(localAxisA);
        if (th < lower) err += correctAngle(bodies, axis * (th - lower), 0.5f);
        else if (th > upper) err += correctAngle(bodies, axis * (th - upper), 0.5f);
    }
    return err;
}

// ---------------------------------------------------------------------------
// Slider
// ---------------------------------------------------------------------------
float SliderJoint::translation(const std::vector<RigidBody>& bodies) const {
    return dot(worldAnchorB(bodies) - worldAnchorA(bodies), worldAxis(bodies));
}

void SliderJoint::buildRows(const std::vector<RigidBody>& bodies, float) {
    const RigidBody& A = bodyOf(bodies, a);
    const RigidBody& B = bodyOf(bodies, b);
    const Vector3 axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (const Vector3& e : axes) { // no relative rotation
        JacobianRow r;
        r.angA = e;
        r.angB = -e;
        addRow(bodies, r);
    }
    Vector3 rA = A.rot.rotate(localAnchorA), rB = B.rot.rotate(localAnchorB);
    Vector3 d = (B.pos + rB) - (A.pos + rA);
    Vector3 axis = A.rot.rotate(localAxisA);
    Vector3 t1 = anyPerpendicular(axis), t2 = cross(axis, t1);
    auto linearRow = [&](const Vector3& t) { // C = t.(pB - pA), t rotating with A
        JacobianRow r;
        r.linA = -t;
        r.angA = -cross(rA + d, t);
        r.linB = t;
        r.angB = cross(rB, t);
        return r;
    };
    addRow(bodies, linearRow(t1));
    addRow(bodies, linearRow(t2));
    if (limitEnabled) {
        float s = dot(d, axis);
        if (s <= lower || s >= upper) {
            JacobianRow r = linearRow(axis);
            if (s <= lower) r.lo = 0;
            else r.hi = 0;
            addRow(bodies, r);
        }
    }
}

float SliderJoint::solvePosition(std::vector<RigidBody>& bodies) {
    float err = correctAngle(bodies, rotationError(bodies), 0.5f);
    Vector3 pA = worldAnchorA(bodies), pB = worldAnchorB(bodies);
    Vector3 axis = worldAxis(bodies);
    Vector3 d = pB - pA;
    Vector3 perp = d - axis * dot(axis, d);
    err += correctPoint(bodies, -perp, pA, pB, 0.5f);
    if (limitEnabled) {
        pA = worldAnchorA(bodies);
        pB = worldAnchorB(bodies);
        float s = dot(pB - pA, axis);
        // correctAlong changes axis.(pA - pB) = -s; its target is -lower / -upper.
        if (s < lower) err += correctAlong(bodies, axis, -s + lower, pA, pB, 0.5f);
        else if (s > upper) err += correctAlong(bodies, axis, -s + upper, pA, pB, 0.5f);
    }
    return err;
}

// ---------------------------------------------------------------------------
// Fixed (weld)
// ---------------------------------------------------------------------------
void FixedJoint::buildRows(const std::vector<RigidBody>& bodies, float h) {
    BallJoint::buildRows(bodies, h);
    const Vector3 axes[3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (const Vector3& e : axes) {
        JacobianRow r;
        r.angA = e;
        r.angB = -e;
        addRow(bodies, r);
    }
}

float FixedJoint::solvePosition(std::vector<RigidBody>& bodies) {
    float err = correctAngle(bodies, rotationError(bodies), 0.5f);
    return err + BallJoint::solvePosition(bodies);
}

// ---------------------------------------------------------------------------
// Distance: rod / rope / spring
// ---------------------------------------------------------------------------
void DistanceJoint::buildRows(const std::vector<RigidBody>& bodies, float h) {
    const RigidBody& A = bodyOf(bodies, a);
    const RigidBody& B = bodyOf(bodies, b);
    Vector3 rA = A.rot.rotate(localAnchorA), rB = B.rot.rotate(localAnchorB);
    Vector3 dvec = (A.pos + rA) - (B.pos + rB);
    float len = rf::length(dvec);
    if (len < 1e-6f) return;
    Vector3 n = dvec / len;
    float C = len - length;
    if (rope && C < 0) return; // slack rope carries nothing
    JacobianRow r;
    r.linA = n;
    r.angA = cross(rA, n);
    r.linB = -n;
    r.angB = -cross(rB, n);
    if (rope) r.hi = 0; // a rope only pulls
    if (frequency > 0) {
        // Soft constraint (spring-damper) in the Box2D formulation.
        float k0 = A.invMass + dot(r.angA, A.applyInvInertiaWorld(r.angA)) + B.invMass + dot(r.angB, B.applyInvInertiaWorld(r.angB));
        float mEff = k0 > 1e-12f ? 1.0f / k0 : 0.0f;
        float omega = 2.0f * kPi * frequency;
        float stiff = mEff * omega * omega, damp = 2.0f * mEff * dampingRatio * omega;
        float g = h * (damp + h * stiff);
        r.softness = g > 0 ? 1.0f / g : 0.0f;
        r.bias = C * h * stiff * r.softness;
    }
    addRow(bodies, r);
}

float DistanceJoint::solvePosition(std::vector<RigidBody>& bodies) {
    if (frequency > 0) return 0; // springs are meant to stretch
    Vector3 pA = worldAnchorA(bodies), pB = worldAnchorB(bodies);
    Vector3 dvec = pA - pB;
    float len = rf::length(dvec);
    if (len < 1e-6f) return 0;
    float C = len - length;
    if (rope && C < 0) return 0;
    return correctAlong(bodies, dvec / len, C, pA, pB, 0.5f);
}

} // namespace rf
