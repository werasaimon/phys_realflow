#pragma once
// Narrow phase: contact manifolds between posed convex shapes.
//
//   sphere-sphere, sphere-box      analytic
//   box-box                        SAT (15 axes) + Sutherland-Hodgman clipping / edge-edge
//   sphere-convex                  GJK distance (point vs shape) + EPA when the centre is inside
//   convex-convex (hull, triangle) GJK + EPA for the normal/depth, then the supporting faces of both
//                                  shapes ("boundary simplices") are clipped against each other (as in
//                                  Jolt); vertex contacts fall back to the perturbation method (Bullet)
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
    // Manifold from the supporting faces along n (from B to A); false for vertex/edge-face cases.
    static bool faceManifold(const PosedShape& A, const PosedShape& B, const Vector3& n, std::vector<ContactPoint>& pts);

private:
    static constexpr int kTypes = 5; // Compound is dispatched part by part before the table
    Algorithm table_[kTypes][kTypes] = {};
};

} // namespace rf
