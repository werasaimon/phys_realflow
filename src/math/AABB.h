#pragma once
// Axis-aligned bounding box (float). Empty (invalid) until something is added.

#include "math/Vector3.h"

namespace rf {

struct AABB {
    Vector3 lo{kInf, kInf, kInf};
    Vector3 hi{-kInf, -kInf, -kInf};

    AABB() = default;
    AABB(const Vector3& l, const Vector3& h) : lo(l), hi(h) {}

    bool valid() const { return lo.x <= hi.x && lo.y <= hi.y && lo.z <= hi.z; }
    void expand(const Vector3& p) { lo = vmin(lo, p); hi = vmax(hi, p); }
    void expand(const AABB& b) { lo = vmin(lo, b.lo); hi = vmax(hi, b.hi); }
    Vector3 center() const { return (lo + hi) * 0.5f; }
    Vector3 extent() const { return hi - lo; }
    float surfaceArea() const {
        Vector3 e = extent();
        return 2.0f * (e.x * e.y + e.y * e.z + e.z * e.x);
    }
    bool overlaps(const AABB& b) const {
        return lo.x <= b.hi.x && hi.x >= b.lo.x && lo.y <= b.hi.y && hi.y >= b.lo.y && lo.z <= b.hi.z && hi.z >= b.lo.z;
    }
    bool contains(const Vector3& p) const {
        return p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y && p.z >= lo.z && p.z <= hi.z;
    }
    float distance2(const Vector3& p) const {
        Vector3 d = vmax(vmax(lo - p, p - hi), Vector3(0.0f));
        return length2(d);
    }
    // Slab test; returns entry distance or +inf when missed.
    float rayHit(const Vector3& o, const Vector3& invDir, float tmax) const {
        float t0 = 0, t1 = tmax;
        for (int a = 0; a < 3; ++a) {
            float tn = (lo[a] - o[a]) * invDir[a];
            float tf = (hi[a] - o[a]) * invDir[a];
            if (tn > tf) std::swap(tn, tf);
            t0 = tn > t0 ? tn : t0;
            t1 = tf < t1 ? tf : t1;
            if (t0 > t1) return kInf;
        }
        return t0;
    }
};

} // namespace rf
