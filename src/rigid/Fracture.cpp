// Voronoi cells of a convex hull (Fracture.h): a cell starts as the whole hull, kept as a list
// of flat polygons, and is cut by one plane per other seed. Cutting a convex polytope by a plane
// is Sutherland-Hodgman on every face plus one new face - the cap - made of the cut edges. At
// the end the polygons are fanned into triangles and the vertices welded, so each cell is a
// closed convex triangle mesh that ConvexHullShape accepts.
#include "rigid/Fracture.h"

#include <algorithm>
#include <cmath>

namespace rf {

namespace {

// A convex polytope as its faces, each a loop of vertices wound counter-clockwise seen from
// outside (so the loop's normal points out).
struct Polygon {
    std::vector<Vector3> loop;
};
using Polytope = std::vector<Polygon>;

Polytope polytopeFromMesh(const TriMesh& mesh) {
    Polytope faces;
    faces.reserve(mesh.triangles.size());
    for (const auto& t : mesh.triangles)
        faces.push_back({{mesh.positions[t[0]], mesh.positions[t[1]], mesh.positions[t[2]]}});
    return faces;
}

// A cut edge: where the plane crossed one face (two points). The cap is built from them.
struct CutEdge {
    Vector3 a, b;
};

// Sutherland-Hodgman: the part of the loop with n . x <= d. Returns false when nothing is left.
// `eps` is the distance within which a vertex counts as on the plane (no sliver triangles).
bool clipPolygon(const std::vector<Vector3>& in, const Vector3& n, float d, float eps, std::vector<Vector3>& out, CutEdge& cut, bool& wasCut) {
    out.clear();
    wasCut = false;
    int crossings = 0;
    const size_t m = in.size();
    for (size_t i = 0; i < m; ++i) {
        const Vector3& p = in[i];
        const Vector3& q = in[(i + 1) % m];
        const float sp = dot(n, p) - d, sq = dot(n, q) - d;
        const bool pIn = sp <= eps, qIn = sq <= eps;
        if (pIn) out.push_back(p);
        // The edge crosses the plane strictly: add the crossing point and remember it for the cap.
        if (pIn != qIn && std::fabs(sp - sq) > 1e-20f) {
            const Vector3 x = p + (q - p) * (sp / (sp - sq));
            out.push_back(x);
            (crossings == 0 ? cut.a : cut.b) = x;
            ++crossings;
            wasCut = true;
        }
    }
    if (crossings != 2) wasCut = false; // touching a vertex only, or nothing crossed
    return out.size() >= 3;
}

// The cap: the cut points sorted by angle around their centre in the plane, wound so that the
// cap's normal is n (n points out of the kept half-space, so the cap faces outward).
Polygon buildCap(const std::vector<CutEdge>& cuts, const Vector3& n, float eps) {
    std::vector<Vector3> pts;
    for (const CutEdge& c : cuts) {
        for (const Vector3& x : {c.a, c.b}) {
            bool dup = false;
            for (const Vector3& y : pts) dup |= length2(x - y) < eps * eps;
            if (!dup) pts.push_back(x);
        }
    }
    Polygon cap;
    if (pts.size() < 3) return cap;
    Vector3 centre(0.0f);
    for (const Vector3& p : pts) centre += p;
    centre /= float(pts.size());
    const Vector3 u = normalize(std::fabs(n.x) < 0.9f ? cross(n, Vector3(1, 0, 0)) : cross(n, Vector3(0, 1, 0)));
    const Vector3 v = cross(n, u); // u x v = n: counter-clockwise in (u, v) is seen from +n
    std::sort(pts.begin(), pts.end(), [&](const Vector3& a, const Vector3& b) {
        return std::atan2(dot(a - centre, v), dot(a - centre, u)) < std::atan2(dot(b - centre, v), dot(b - centre, u));
    });
    cap.loop = std::move(pts);
    return cap;
}

// Cuts the polytope by the plane n . x <= d in place. False when the polytope vanished.
bool clipPolytope(Polytope& faces, const Vector3& n, float d, float eps) {
    Polytope kept;
    std::vector<CutEdge> cuts;
    std::vector<Vector3> out;
    for (const Polygon& f : faces) {
        CutEdge cut;
        bool wasCut = false;
        if (clipPolygon(f.loop, n, d, eps, out, cut, wasCut)) kept.push_back({out});
        if (wasCut) cuts.push_back(cut);
    }
    if (kept.empty()) return false;
    if (cuts.size() >= 3) {
        Polygon cap = buildCap(cuts, n, eps);
        if (cap.loop.size() >= 3) kept.push_back(std::move(cap));
    }
    faces.swap(kept);
    return true;
}

// Fans a polygon into triangles from the vertex whose fan has the largest smallest triangle.
// Several cuts through the same edge leave collinear vertices on a loop; they must stay (the
// neighbouring face has them too, and the mesh is only watertight if both sides share them),
// so the fan apex is chosen off their line - a fan from a collinear vertex would make a
// zero-area triangle with a meaningless normal.
void fanPolygon(const std::vector<Vector3>& loop, uint32_t base, std::vector<std::array<uint32_t, 3>>& triangles) {
    const uint32_t n = uint32_t(loop.size());
    uint32_t apex = 0;
    float best = -1.0f;
    for (uint32_t a = 0; a < n; ++a) {
        float smallest = 1e30f;
        for (uint32_t k = 1; k + 1 < n; ++k) {
            const uint32_t i = (a + k) % n, j = (a + k + 1) % n;
            smallest = std::min(smallest, length(cross(loop[i] - loop[a], loop[j] - loop[a])));
        }
        if (smallest > best) {
            best = smallest;
            apex = a;
        }
    }
    if (best <= 0.0f) return; // every vertex on one line: no area, nothing to draw
    for (uint32_t k = 1; k + 1 < n; ++k) triangles.push_back({base + apex, base + (apex + k) % n, base + (apex + k + 1) % n});
}

// Every polygon into triangles, then the shared vertices welded and exact slivers dropped.
TriMesh meshFromPolytope(const Polytope& faces, float weldEps) {
    TriMesh m;
    for (const Polygon& f : faces) {
        const uint32_t base = uint32_t(m.positions.size());
        for (const Vector3& p : f.loop) m.positions.push_back(p);
        fanPolygon(f.loop, base, m.triangles);
    }
    m.weld(weldEps);
    m.removeDegenerate();
    m.orientOutward();
    return m;
}

struct Random {
    uint32_t seed;
    float next() {
        seed = seed * 1664525u + 1013904223u;
        return (seed >> 8) * (1.0f / 16777216.0f);
    }
};

} // namespace

bool insideConvex(const TriMesh& hull, const Vector3& p, float tolerance) {
    for (size_t t = 0; t < hull.triangles.size(); ++t) {
        const Vector3 n = hull.faceNormal(t);
        if (dot(n, p - hull.positions[hull.triangles[t][0]]) > tolerance) return false;
    }
    return true;
}

std::vector<TriMesh> voronoiCells(const TriMesh& hull, const std::vector<Vector3>& seeds) {
    const float size = maxComp(hull.bounds().extent());
    const float eps = 1e-6f * size;
    std::vector<TriMesh> cells;
    cells.reserve(seeds.size());
    for (size_t i = 0; i < seeds.size(); ++i) {
        Polytope cell = polytopeFromMesh(hull);
        bool alive = true;
        for (size_t j = 0; j < seeds.size() && alive; ++j) {
            if (j == i) continue;
            // The bisector: keep the half-space nearer to seed i. n . x <= n . midpoint.
            const Vector3 n = normalize(seeds[j] - seeds[i]);
            const float d = dot(n, (seeds[i] + seeds[j]) * 0.5f);
            alive = clipPolytope(cell, n, d, eps);
        }
        if (!alive) continue;
        TriMesh m = meshFromPolytope(cell, 4.0f * eps);
        if (!m.empty()) cells.push_back(std::move(m));
    }
    return cells;
}

std::vector<Vector3> uniformSeeds(const TriMesh& hull, int count, uint32_t seed) {
    Random rnd{seed};
    const AABB box = hull.bounds();
    std::vector<Vector3> out;
    for (int tries = 0; int(out.size()) < count && tries < 100 * count; ++tries) {
        const Vector3 p = box.lo + box.extent() * Vector3(rnd.next(), rnd.next(), rnd.next());
        if (insideConvex(hull, p)) out.push_back(p);
    }
    return out;
}

std::vector<Vector3> impactSeeds(const TriMesh& hull, const Vector3& impact, int count, float radius, uint32_t seed) {
    Random rnd{seed};
    std::vector<Vector3> out;
    for (int tries = 0; int(out.size()) < count && tries < 100 * count; ++tries) {
        // A random direction (uniform on the sphere) and a distance crowded towards the impact.
        const float z = 2.0f * rnd.next() - 1.0f, phi = 2.0f * kPi * rnd.next();
        const float s = std::sqrt(std::max(0.0f, 1.0f - z * z));
        const Vector3 dir(s * std::cos(phi), s * std::sin(phi), z);
        const float u = rnd.next();
        const Vector3 p = impact + dir * (radius * u * u);
        if (insideConvex(hull, p)) out.push_back(p);
    }
    return out;
}

} // namespace rf
