#pragma once
// Continuous collision detection: conservative advancement (B. Mirtich, 1996) on GJK distances.
//
// The motion of each shape over the step is parameterised by s in [0, 1]: position linear, the
// orientation by the exponential map of its rotation vector. GJK supplies a search direction;
// a fresh support plane bounds the gap from below. The gap cannot close faster than
//     bound = |dP_rel| + |dTheta_A| r_A + |dTheta_B| r_B      (per unit s),
// so a positive support gap permits conservative advancement in exact arithmetic (Mirtich 1996).
// ConservativeAdvancement owns this search. Float support/rotation arithmetic is not interval
// certified; see docs/book/09-ccd-contract.md for the implementation's explicit compromises.

#include "rigid/GjkEpa.h"

#include <cstddef>

namespace rf {

struct SweptPose {
    const ConvexShape* shape = nullptr;
    Vector3 p0, p1;   // positions of the body at s = 0 and s = 1
    Quaternion q0;       // orientation at s = 0
    Vector3 dTheta;   // rotation vector over the whole step (q1 = exp(dTheta) q0)
    // A convex part of a compound body: its pose in the body frame (identity: the whole body).
    Matrix3x3 partR;
    Vector3 partT;
    // Quadratic translation: p(s) = lerp(p0,p1,s) + translationCurve*s*(s-1).
    // Constant acceleration a over duration h gives translationCurve = a*h*h/2.
    Vector3 translationCurve;
    PosedShape at(float s) const {
        const Matrix3x3 R = q0.integrated(dTheta, s).toMatrix3x3();
        return {shape, R * partR, p0 + (p1 - p0) * s + translationCurve * (s * (s - 1)) + R * partT};
    }
    // How far a surface point can move by the rotation; a centred sphere looks the same at any orientation.
    float angularReach() const {
        if (shape->type() == ShapeType::Sphere && length2(partT) == 0) return 0.0f;
        float radius = shape->boundingRadius();
        if (shape->type() == ShapeType::Triangle) {
            const auto& triangle = static_cast<const TriangleShape&>(*shape);
            for (int i = 0; i < 3; ++i) radius = std::max(radius, length(triangle.vertex(i)));
        }
        return length(dTheta) * (length(partT) + radius);
    }
};

enum class ToiStatus { Separated, Impact, InitialContact, Unresolved };
enum class ToiReason { None, IterationLimit, NoProgress, NoSeparatingPlane, InvalidInput };

struct ToiResult {
    bool hit = false;   // compatibility: true only for Impact; inspect status for uncertainty
    float s = 0;       // stop fraction for Impact/Unresolved; 1 for Separated, 0 for InitialContact
    int iterations = 0;
    ToiStatus status = ToiStatus::Unresolved;
    ToiReason reason = ToiReason::IterationLimit;
    bool needsClamping() const { return status == ToiStatus::Impact || status == ToiStatus::Unresolved; }
};

// Counts convex/plane queries over all trials of the last requested impulse-solver step.
struct CcdDiagnostics {
    size_t queries = 0;
    size_t unresolved = 0;
    size_t initialContacts = 0;
    bool passLimitReached = false; // legacy clamping diagnostic; retries use RigidStepResult instead
    float observe(const ToiResult& result) {
        ++queries;
        unresolved += result.status == ToiStatus::Unresolved;
        initialContacts += result.status == ToiStatus::InitialContact;
        return result.needsClamping() ? result.s : 1.0f;
    }
};

// Contains the entire interpolated motion, including the arc between the endpoint orientations.
AABB sweptBounds(const SweptPose& sweep, float tolerance = 0);

// Time of impact between two swept convex shapes (either may be static: p0 == p1, dTheta == 0).
// Initially touching/overlapping pairs report InitialContact, not Separated. RigidWorld delegates
// those to its discrete solver. An exhausted/invalid query reports Unresolved, never Separated.
ToiResult timeOfImpact(const SweptPose& A, const SweptPose& B, float tolerance, int maxIterations = 64);

// Swept shape against a static plane {x : n.x = d} (domain walls), exact support distance.
ToiResult timeOfImpactPlane(const SweptPose& A, const Vector3& n, float d, float tolerance, int maxIterations = 64);

} // namespace rf
