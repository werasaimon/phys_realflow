#pragma once
// Continuous collision detection: conservative advancement (B. Mirtich, 1996) on GJK distances.
//
// The motion of each shape over the step is parameterised by s in [0, 1]: position linear, the
// orientation by the exponential map of its rotation vector. At the current s GJK gives the exact
// distance d between the posed shapes; the relative motion cannot close a gap faster than
//     bound = |dP_rel| + |dTheta_A| r_A + |dTheta_B| r_B      (per unit s),
// so s can safely advance by d / bound. The loop stops when d < tolerance (time of impact) or
// s passes 1 (no collision within the step). Works for every convex shape (support mapping).

#include "rigid/GjkEpa.h"

namespace rf {

struct SweptPose {
    const ConvexShape* shape = nullptr;
    Vector3 p0, p1;   // positions at s = 0 and s = 1
    Quaternion q0;       // orientation at s = 0
    Vector3 dTheta;   // rotation vector over the whole step (q1 = exp(dTheta) q0)
    PosedShape at(float s) const {
        return {shape, q0.integrated(dTheta, s).toMatrix3x3(), p0 + (p1 - p0) * s};
    }
    // How far a surface point can move by the rotation; a sphere looks the same at any orientation.
    float angularReach() const { return shape->type() == ShapeType::Sphere ? 0.0f : length(dTheta) * shape->boundingRadius(); }
};

struct ToiResult {
    bool hit = false;
    float s = 1;        // fraction of the step at the time of impact
    int iterations = 0;
};

// Time of impact between two swept convex shapes (either may be static: p0 == p1, dTheta == 0).
// Pairs already overlapping at s = 0 report no hit (resting contacts belong to the discrete solver).
ToiResult timeOfImpact(const SweptPose& A, const SweptPose& B, float tolerance, int maxIterations = 64);

// Swept shape against a static plane {x : n.x = d} (domain walls), exact support distance.
ToiResult timeOfImpactPlane(const SweptPose& A, const Vector3& n, float d, float tolerance, int maxIterations = 64);

} // namespace rf
