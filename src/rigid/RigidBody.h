#pragma once
// Rigid body state shared by the world, the contact solver and the joints.

#include "math/Math.h"
#include "rigid/GjkEpa.h"
#include "rigid/Shapes.h"

#include <memory>

namespace rf {

struct RigidBody {
    std::shared_ptr<const ConvexShape> shape;
    Vector3 pos, vel, angVel;
    Quaternion rot;
    float mass = 1, invMass = 1;
    Vector3 invInertiaLocal{1.0f}; // diagonal, principal body frame
    float restitution = 0.2f;
    float friction = 0.5f;       // kinetic (sliding) friction coefficient
    float staticFriction = 0.7f; // static (sticking) friction coefficient, >= kinetic
    Vector3 color{0.8f, 0.5f, 0.2f};
    // What the body looks like, apart from what it collides with (its `shape`): a mesh in the BODY
    // frame (vertex_world = pos + rot * v), e.g. a horse model on a capsule collider. Null: the
    // collider itself is drawn. The solver never reads it; it is for viewers only.
    std::shared_ptr<const TriMesh> visualMesh;
    Vector3 force, torque; // accumulated external loads for the next step
    // Split pseudo velocities correct poses without directly changing physical velocities.
    Vector3 biasVel, biasAngVel;
    // Island sleeping (Box2D-style): time spent below the sleep thresholds, and the state.
    float sleepTimer = 0;
    bool sleeping = false;
    int sleepIsland = -1; // island the body fell asleep with: it wakes as a whole
    // Pose at the start of the (sub)step: XPBD velocity update and CCD sweeps.
    Vector3 prevPos;
    Quaternion prevRot;
    // false: the slot of a destroyed body (RigidWorld::destroyBody) - static, parked far away.
    bool alive = true;

    ShapeType type() const { return shape->type(); }
    float radius() const { return shape->boundingRadius(); }
    Vector3 halfExtents() const {
        return type() == ShapeType::Box ? static_cast<const BoxShape*>(shape.get())->halfExtents() : Vector3(radius());
    }
    Matrix3x3 rotation() const { return rot.toMatrix3x3(); }
    PosedShape posed() const { return {shape.get(), rotation(), pos}; }
    // World-space inverse inertia, cached once per step (the solver applies it thousands of times).
    Matrix3x3 invInertiaWorld = Matrix3x3::zero();
    void updateInertia() {
        Matrix3x3 R = rotation();
        invInertiaWorld = R * Matrix3x3::diag(invInertiaLocal) * R.transposed();
    }
    Vector3 applyInvInertiaWorld(const Vector3& v) const { return invInertiaWorld * v; }
    Vector3 velocityAt(const Vector3& p) const { return vel + cross(angVel, p - pos); }
    AABB worldBounds() const { return shape->boundsAt(rotation(), pos); }
    // Signed distance from world point p to the body surface (negative inside) and outward normal.
    float signedDistance(const Vector3& p, Vector3& normal) const;
    float boundingRadius() const { return shape->boundingRadius(); }
};

} // namespace rf
