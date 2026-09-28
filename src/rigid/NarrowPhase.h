#pragma once
// Narrow phase: contact manifolds between posed convex shapes.
//
//   sphere-sphere, sphere-box      analytic
//   box-box                        SAT (15 axes) + Sutherland-Hodgman clipping / edge-edge
//   sphere-convex                  GJK distance (point vs shape) + EPA when the centre is inside
//   convex-convex (hull, triangle) GJK + EPA for the normal/depth, then the supporting faces of both
//                                  shapes ("boundary simplices") are clipped against each other (as in
//                                  Jolt); vertex contacts fall back to the perturbation method (Bullet)
//   capsule-anything               the capsule is its core segment plus a radius: segment-point,
//                                  segment-segment analytic; segment-convex by GJK/EPA on the core,
//                                  and a capsule lying on a face touches it at its two ends
//
// Algorithms are chosen through a double-dispatch table indexed by the two shape types.

#include "rigid/GjkEpa.h"

#include <vector>

namespace rf {

struct ContactPoint {
    Vector3 position; // world, midway between the two surfaces
    Vector3 normal;   // unit, from B towards A (pushes A out of B)
    float depth = 0; // > 0 penetration, < 0 speculative gap (within the contact margin)
};

struct ContactManifold {
    std::vector<ContactPoint> points;
    void add(const Vector3& p, const Vector3& n, float depth) { points.push_back({p, n, depth}); }
};

// Keeps at most `maxPoints`: the deepest point plus the points spanning the largest area.
void reduceManifold(std::vector<ContactPoint>& pts, size_t maxPoints = 4);

class NarrowPhase {
public:
    using Algorithm = bool (*)(const PosedShape& A, const PosedShape& B, ContactManifold& m);

    NarrowPhase();

    // Shapes closer than this produce (speculative) contacts with negative depth, so contacts
    // persist through tiny separations and warm starting keeps working.
    static inline float margin = 0.01f;
    // Appends contacts (normals from B to A) and returns true when the shapes touch.
    bool collide(const PosedShape& A, const PosedShape& B, ContactManifold& m) const;

    static bool sphereSphere(const PosedShape& A, const PosedShape& B, ContactManifold& m);
    static bool sphereBox(const PosedShape& A, const PosedShape& B, ContactManifold& m);
    static bool sphereConvex(const PosedShape& A, const PosedShape& B, ContactManifold& m);
    static bool boxBox(const PosedShape& A, const PosedShape& B, ContactManifold& m);
    static bool convexConvex(const PosedShape& A, const PosedShape& B, ContactManifold& m);
    static bool sphereCapsule(const PosedShape& A, const PosedShape& B, ContactManifold& m);
    static bool capsuleCapsule(const PosedShape& A, const PosedShape& B, ContactManifold& m);
    static bool capsuleConvex(const PosedShape& A, const PosedShape& B, ContactManifold& m);
    // Manifold from the supporting faces along n (from B to A), the deepest point kept within
    // [depthLow, depthHigh] (EPA's bracket of the pair's depth; infinite: none); false for
    // vertex/edge-face cases.
    static bool faceManifold(const PosedShape& A, const PosedShape& B, const Vector3& n, float depthLow, float depthHigh,
                             std::vector<ContactPoint>& pts);

private:
    static constexpr int kTypes = 6; // Compound is dispatched part by part before the table
    Algorithm table_[kTypes][kTypes] = {};

    // The steps of boxBox(): the best separating axis of the 15 (SAT), and the contact it gives.
    struct BoxAxis {
        float sep = -kInf; // separation along the axis (negative: overlap)
        Vector3 axis;      // from A towards B
        int kind = 0;      // 0: a face of A, 1: a face of B, 2: an edge of A against an edge of B
        int i = 0, j = 0;  // the face / edge indices
    };
    static bool boxSeparatingAxes(const PosedShape& A, const PosedShape& B, BoxAxis& best, BoxAxis& bestEdge);
    static void boxEdgeContact(const PosedShape& A, const PosedShape& B, const BoxAxis& e, ContactManifold& m);
    static bool boxFaceContact(const PosedShape& A, const PosedShape& B, const BoxAxis& best, std::vector<ContactPoint>& pts);
    // The steps of faceManifold() and convexConvex().
    static void clipIncidentFace(const std::vector<Vector3>& ref, const Vector3& refNormal, std::vector<Vector3>& inc,
                                 std::vector<Vector3>& clipped);
    static void perturbationManifold(const PosedShape& A, const PosedShape& B, const Vector3& n, std::vector<ContactPoint>& pts);
    // The step of capsuleConvex(): a capsule lying on B's face touches it at both ends of its side.
    static bool capsuleOnFace(const PosedShape& A, const PosedShape& B, const Vector3& n, ContactManifold& m);
};

} // namespace rf
