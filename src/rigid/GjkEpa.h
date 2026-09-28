#pragma once
// GJK (Gilbert-Johnson-Keerthi) distance / intersection and EPA (Expanding Polytope Algorithm)
// penetration depth for two posed convex shapes, working on the Minkowski difference A - B.

#include "rigid/Shapes.h"

#include <vector>

namespace rf {

// A convex shape placed in the world: x_world = R * x_local + p.
struct PosedShape {
    const ConvexShape* shape = nullptr;
    Matrix3x3 R;
    Vector3 p;
    Vector3 support(const Vector3& d) const { return p + R * shape->support(R.transposed() * d); }
    Vector3 center() const { return p; }
    // Supporting feature in world space (see ConvexShape::supportFeature).
    void feature(const Vector3& d, std::vector<Vector3>& out) const {
        shape->supportFeature(R.transposed() * d, out);
        for (Vector3& v : out) v = p + R * v;
    }
};

struct GjkResult {
    bool intersect = false;
    float distance = 0;   // valid when !intersect
    Vector3 pointA, pointB;  // closest points (world) when separated
    // Final simplex, handed to EPA when intersecting.
    int simplexSize = 0;
    Vector3 w[4], a[4], b[4];
};

// maxDistance: stop as soon as the shapes are proven farther apart than this (the support point
// along -v bounds the distance from below: dist >= dot(v, w) / |v|); `distance` then holds that
// lower bound (> maxDistance) and the closest points are not computed. As Jolt's inMaxDistSq.
GjkResult gjk(const PosedShape& A, const PosedShape& B, float maxDistance = kInf);

struct PenetrationResult {
    bool valid = false;
    Vector3 normal;      // unit, from B towards A (moving A along +normal by depth separates them)
    float depth = 0;
    // How far A must go along the normal lies between depth and depthMax: EPA's polytope lies
    // inside the Minkowski difference A - B (its face is no farther than the boundary), and the
    // support of A - B along the normal is as far as the boundary can be. Converged, the two agree
    // to EPA's tolerance; stopped early (a touch it cannot resolve), depth is a bound from below.
    float depthMax = 0;
    Vector3 pointA;      // deepest point of A (inside B)
    Vector3 pointB;      // deepest point of B (inside A)
};

// EPA starting from an intersecting GJK simplex.
PenetrationResult epa(const PosedShape& A, const PosedShape& B, const GjkResult& g);

// Convenience: GJK, then EPA when the shapes overlap. Returns false when they are separated.
bool penetration(const PosedShape& A, const PosedShape& B, PenetrationResult& out);

// A record of one GJK (+ EPA) run, for the research view of a watched pair: the simplex after
// every GJK iteration (each vertex a point w = a - b of the Minkowski difference with the support
// points a on A and b on B that made it) and the final EPA polytope. Recording happens only on the
// thread that set a trace, and only for the runs between setGjkTrace(&t) and setGjkTrace(nullptr):
// the hot path pays one pointer test.
struct GjkTrace {
    struct Vertex { Vector3 w, a, b; };
    std::vector<std::vector<Vertex>> simplices; // one per GJK iteration (the last run only)
    std::vector<Vertex> polytope;               // the EPA polytope's vertices (the last run)
    std::vector<int> faces;                     // its final faces, three indices each
    int epaIterations = 0;
    void clear() { simplices.clear(); polytope.clear(); faces.clear(); epaIterations = 0; }
};
void setGjkTrace(GjkTrace* trace);

} // namespace rf
