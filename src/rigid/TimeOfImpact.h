#pragma once
// Continuous collision detection: conservative advancement (B. Mirtich, 1996) on GJK distances.
//
// The motion of each shape over the step is parameterised by s in [0, 1]: position linear, the
// orientation by the exponential map of its rotation vector. At the current s GJK estimates the
// distance d between the posed shapes; the relative motion cannot close a gap faster than
//     bound = |dP_rel| + |dTheta_A| r_A + |dTheta_B| r_B      (per unit s),
// so an exact distance permits advancing by d / bound. A finite GJK iterate need not be a lower
// bound: this implementation is not a certified intersection-free trajectory algorithm. The
// remaining advancement/budget and multi-body clamping limitations are recorded in docs/22.

#include "rigid/GjkEpa.h"

namespace rf {

struct SweptPose {
    const ConvexShape* shape = nullptr;
    Vector3 p0, p1;   // positions of the body at s = 0 and s = 1
    Quaternion q0;       // orientation at s = 0
    Vector3 dTheta;   // rotation vector over the whole step (q1 = exp(dTheta) q0)
    // A convex part of a compound body: its pose in the body frame (identity: the whole body).
    Matrix3x3 partR;
    Vector3 partT;
    PosedShape at(float s) const {
        const Matrix3x3 R = q0.integrated(dTheta, s).toMatrix3x3();
        return {shape, R * partR, p0 + (p1 - p0) * s + R * partT};
    }
    // How far a surface point can move by the rotation; a centred sphere looks the same at any orientation.
    float angularReach() const {
        if (shape->type() == ShapeType::Sphere && length2(partT) == 0) return 0.0f;
        return length(dTheta) * (length(partT) + shape->boundingRadius());
    }
};

struct ToiResult {
    bool hit = false;
    float s = 1;        // fraction of the step at the time of impact
    int iterations = 0;
};

// Contains the entire interpolated motion, including the arc between the endpoint orientations.
AABB sweptBounds(const SweptPose& sweep, float tolerance = 0);

// Time of impact between two swept convex shapes (either may be static: p0 == p1, dTheta == 0).
// Pairs already overlapping at s = 0 report no hit (resting contacts belong to the discrete solver).
ToiResult timeOfImpact(const SweptPose& A, const SweptPose& B, float tolerance, int maxIterations = 64);

// Swept shape against a static plane {x : n.x = d} (domain walls), exact support distance.
ToiResult timeOfImpactPlane(const SweptPose& A, const Vector3& n, float d, float tolerance, int maxIterations = 64);

} // namespace rf
