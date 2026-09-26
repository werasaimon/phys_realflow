#include "rigid/Ccd.h"

namespace rf {

ToiResult timeOfImpact(const SweptPose& A, const SweptPose& B, float tol, int maxIt) {
    ToiResult r;
    const float bound = length((A.p1 - A.p0) - (B.p1 - B.p0)) + A.angularReach() + B.angularReach();
    if (bound < 1e-9f) return r;
    float s = 0;
    for (int it = 0; it < maxIt; ++it) {
        r.iterations = it + 1;
        GjkResult g = gjk(A.at(s), B.at(s));
        if (g.intersect) {
            if (s == 0) return r; // already touching: the discrete solver owns this pair
            r.hit = true;         // overshot inside the tolerance band: report the last safe s
            r.s = s;
            return r;
        }
        if (g.distance < tol) {
            // Touching at the start of the step: a resting/sliding contact owned by the discrete
            // (speculative) solver - never freeze such pairs.
            if (s == 0) return r;
            r.hit = true;
            r.s = s;
            return r;
        }
        // Conservative step: the gap cannot close faster than `bound` per unit s. Aim slightly
        // short of contact so the next GJK still sees separated shapes.
        s += std::max((g.distance - 0.5f * tol) / bound, 1e-6f);
        if (s > 1.0f) return r;
    }
    return r;
}

ToiResult timeOfImpactPlane(const SweptPose& A, const Vector3& n, float d, float tol, int maxIt) {
    ToiResult r;
    const float bound = std::fabs(dot(A.p1 - A.p0, n)) + A.angularReach();
    if (bound < 1e-9f) return r;
    float s = 0;
    for (int it = 0; it < maxIt; ++it) {
        r.iterations = it + 1;
        PosedShape ps = A.at(s);
        float dist = dot(ps.support(-n), n) - d; // exact distance of the shape to the plane
        if (dist < 0) {
            if (s == 0) return r;
            r.hit = true;
            r.s = s;
            return r;
        }
        if (dist < tol) {
            if (s == 0) return r; // touching the plane already: discrete contact
            r.hit = true;
            r.s = s;
            return r;
        }
        s += std::max((dist - 0.5f * tol) / bound, 1e-6f);
        if (s > 1.0f) return r;
    }
    return r;
}

} // namespace rf
