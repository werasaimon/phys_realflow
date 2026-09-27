// Distance and penetration between two convex shapes: GJK (Gilbert, Johnson & Keerthi 1988) walks
// a simplex in the Minkowski difference towards the origin; when the shapes overlap, EPA (van
// den Bergen 2001) expands a polytope to the boundary and reads the penetration off its closest
// face. Declared in GjkEpa.h.
#include "rigid/GjkEpa.h"

#include "core/Parallel.h"

#include <vector>

namespace rf {

namespace {

struct SV {
    Vector3 w, a, b; // w = a - b
};

SV supportAB(const PosedShape& A, const PosedShape& B, const Vector3& d) {
    Vector3 a = A.support(d), b = B.support(-d);
    return {a - b, a, b};
}

struct Simplex {
    SV v[4];
    float lam[4] = {1, 0, 0, 0};
    int n = 0;
};

// Closest point to the origin on segment / triangle, with barycentric weights.
Vector3 closestSegment(const Vector3& a, const Vector3& b, float bary[2]) {
    Vector3 ab = b - a;
    float l2 = dot(ab, ab);
    float t = l2 > 1e-20f ? clampv(-dot(a, ab) / l2, 0.0f, 1.0f) : 0.0f;
    bary[0] = 1 - t;
    bary[1] = t;
    return a + ab * t;
}

Vector3 closestTriangle(const Vector3& a, const Vector3& b, const Vector3& c, float bary[3]) {
    Vector3 ab = b - a, ac = c - a, ap = -a;
    float d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) { bary[0] = 1; bary[1] = bary[2] = 0; return a; }
    Vector3 bp = -b;
    float d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) { bary[1] = 1; bary[0] = bary[2] = 0; return b; }
    float vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        float v = d1 / (d1 - d3);
        bary[0] = 1 - v; bary[1] = v; bary[2] = 0;
        return a + ab * v;
    }
    Vector3 cp = -c;
    float d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) { bary[2] = 1; bary[0] = bary[1] = 0; return c; }
    float vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        float w = d2 / (d2 - d6);
        bary[0] = 1 - w; bary[1] = 0; bary[2] = w;
        return a + ac * w;
    }
    float va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        bary[0] = 0; bary[1] = 1 - w; bary[2] = w;
        return b + (c - b) * w;
    }
    float denom = va + vb + vc;
    if (std::fabs(denom) < 1e-30f) { bary[0] = 1; bary[1] = bary[2] = 0; return a; }
    float v = vb / denom, w = vc / denom;
    bary[0] = 1 - v - w; bary[1] = v; bary[2] = w;
    return a + ab * v + ac * w;
}

// Keeps only the vertices with non-zero weight.
void compact(Simplex& s) {
    int k = 0;
    for (int i = 0; i < s.n; ++i)
        if (s.lam[i] > 0) {
            s.v[k] = s.v[i];
            s.lam[k] = s.lam[i];
            ++k;
        }
    s.n = std::max(k, 1);
}

// The tetrahedron case of solveSimplex: the origin is inside when it lies on the inner side of
// all four faces; otherwise the simplex shrinks to the face closest to the origin. Side tests in
// double precision with a relative tolerance: thin shapes give nearly flat Minkowski tetrahedra
// where float rounding flips the sign and fakes an "inside".
Vector3 closestOnTetrahedron(Simplex& s, bool& inside) {
    const int faces[4][4] = {{0, 1, 2, 3}, {0, 3, 1, 2}, {0, 2, 3, 1}, {1, 3, 2, 0}}; // i, j, k, opposite
    bool anyOutside = false;
    float best = kInf;
    Simplex bestS;
    Vector3 bestV;
    double scale = 0;
    for (int i = 0; i < 4; ++i) scale = std::max(scale, double(length(s.v[i].w)));
    for (auto& f : faces) {
        const Vector3 &a = s.v[f[0]].w, &b = s.v[f[1]].w, &c = s.v[f[2]].w, &d = s.v[f[3]].w;
        double ab[3] = {double(b.x) - a.x, double(b.y) - a.y, double(b.z) - a.z};
        double ac[3] = {double(c.x) - a.x, double(c.y) - a.y, double(c.z) - a.z};
        double n[3] = {ab[1] * ac[2] - ab[2] * ac[1], ab[2] * ac[0] - ab[0] * ac[2], ab[0] * ac[1] - ab[1] * ac[0]};
        double nl = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        double sO = -(n[0] * a.x + n[1] * a.y + n[2] * a.z);
        double sD = n[0] * (double(d.x) - a.x) + n[1] * (double(d.y) - a.y) + n[2] * (double(d.z) - a.z);
        const double eps = 1e-6 * nl * std::max(scale, 1e-9);
        bool flat = std::fabs(sD) <= eps;
        // "Inside" of this face only when the origin is clearly on the side of the opposite vertex.
        bool insideFace = !flat && sO * sD > 0 && std::fabs(sO) > eps;
        if (insideFace) continue;
        anyOutside = true;
        float bary[3];
        Vector3 v = closestTriangle(a, b, c, bary);
        float l2 = dot(v, v);
        if (l2 < best) {
            best = l2;
            bestV = v;
            bestS.n = 3;
            bestS.v[0] = s.v[f[0]]; bestS.v[1] = s.v[f[1]]; bestS.v[2] = s.v[f[2]];
            bestS.lam[0] = bary[0]; bestS.lam[1] = bary[1]; bestS.lam[2] = bary[2];
        }
    }
    if (!anyOutside) {
        inside = true;
        return Vector3(0.0f);
    }
    s = bestS;
    compact(s);
    return bestV;
}

// Replaces the simplex by the smallest sub-simplex containing its point closest to the origin.
// Returns that point; `inside` is set when the tetrahedron contains the origin.
Vector3 solveSimplex(Simplex& s, bool& inside) {
    inside = false;
    switch (s.n) {
    case 1:
        s.lam[0] = 1;
        return s.v[0].w;
    case 2: {
        float b[2];
        Vector3 v = closestSegment(s.v[0].w, s.v[1].w, b);
        s.lam[0] = b[0]; s.lam[1] = b[1];
        compact(s);
        return v;
    }
    case 3: {
        float b[3];
        Vector3 v = closestTriangle(s.v[0].w, s.v[1].w, s.v[2].w, b);
        s.lam[0] = b[0]; s.lam[1] = b[1]; s.lam[2] = b[2];
        compact(s);
        return v;
    }
    default:
        return closestOnTetrahedron(s, inside);
    }
}

} // namespace

// ---------------------------------------------------------------------------
// GJK
// ---------------------------------------------------------------------------
GjkResult gjk(const PosedShape& A, const PosedShape& B, float maxDistance) {
    GjkResult r;
    Simplex s;
    Vector3 d = A.center() - B.center();
    if (length2(d) < 1e-12f) d = Vector3(1, 0, 0);
    s.v[0] = supportAB(A, B, d);
    s.n = 1;
    s.lam[0] = 1;
    Vector3 v = s.v[0].w;
    bool inside = false;
    for (int iter = 0; iter < 64; ++iter) {
        float dist2 = dot(v, v);
        if (dist2 < 1e-12f) { inside = true; break; }
        SV w = supportAB(A, B, -v);
        // Early out: the whole Minkowski difference lies beyond the plane dot(v, x) = dot(v, w).
        const float vw = dot(v, w.w);
        if (vw > 0 && vw * vw > maxDistance * maxDistance * dist2) {
            r.intersect = false;
            r.distance = vw / std::sqrt(dist2);
            r.simplexSize = 0;
            return r;
        }
        // No further progress towards the origin: v is (numerically) the closest point.
        if (dist2 - dot(v, w.w) <= 1e-6f * dist2 + 1e-10f) break;
        bool duplicate = false;
        for (int i = 0; i < s.n; ++i) duplicate |= length2(s.v[i].w - w.w) < 1e-14f;
        if (duplicate) break;
        Simplex prev = s;
        s.v[s.n++] = w;
        Vector3 vNew = solveSimplex(s, inside);
        if (inside) break;
        if (dot(vNew, vNew) >= dist2) { s = prev; break; } // no progress: keep the best simplex found
        v = vNew;
    }
    r.intersect = inside;
    r.simplexSize = s.n;
    for (int i = 0; i < s.n; ++i) { r.w[i] = s.v[i].w; r.a[i] = s.v[i].a; r.b[i] = s.v[i].b; }
    if (!inside) {
        r.distance = length(v);
        r.pointA = r.pointB = Vector3(0.0f);
        for (int i = 0; i < s.n; ++i) {
            r.pointA += s.v[i].a * s.lam[i];
            r.pointB += s.v[i].b * s.lam[i];
        }
    }
    return r;
}

// ---------------------------------------------------------------------------
// EPA
// ---------------------------------------------------------------------------
// A face of the EPA polytope.
struct EpaFace {
    int a, b, c;
    Vector3 n;
    float d;
};

// Scratch of EPA, one per thread of the pool, cleared before use and never freed: a deep
// contact costs no allocation (indexed by the worker number, not thread_local objects - see
// ThreadPool::workerIndex).
namespace {
struct EpaScratch {
    std::vector<SV> V;                          // the polytope's vertices
    std::vector<EpaFace> faces, kept;           // its faces, and those that survive an expansion
    std::vector<std::pair<int, int>> horizon;   // the edges the new vertex is joined to
};
EpaScratch& epaScratch() {
    static std::vector<EpaScratch> all(size_t(ThreadPool::instance().threadCount()));
    return all[size_t(ThreadPool::workerIndex())];
}
} // namespace

namespace {

// A face of the polytope through three of its vertices, with its outward unit normal and the
// distance of its plane from the origin (a degenerate face gets an infinite distance: never chosen).
EpaFace makeEpaFace(const std::vector<SV>& V, int a, int b, int c) {
    EpaFace f{a, b, c, cross(V[b].w - V[a].w, V[c].w - V[a].w), 0};
    float l = length(f.n);
    if (l < 1e-12f) { f.n = Vector3(0.0f); f.d = kInf; return f; } // degenerate: never chosen
    f.n /= l;
    f.d = dot(f.n, V[a].w);
    return f;
}

// GJK may end with a point, a segment or a triangle containing the origin; EPA needs a tetrahedron
// with volume, so support points in other directions are added until it has one.
void inflateToTetrahedron(const PosedShape& A, const PosedShape& B, std::vector<SV>& V) {
    const float eps = 1e-6f;
    if (V.size() == 1) {
        const Vector3 dirs[6] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
        for (const Vector3& d : dirs) {
            SV w = supportAB(A, B, d);
            if (length(w.w - V[0].w) > eps) { V.push_back(w); break; }
        }
    }
    if (V.size() == 2) {
        Vector3 line = normalize(V[1].w - V[0].w);
        Vector3 axis = std::fabs(line.x) < 0.57f ? Vector3(1, 0, 0) : (std::fabs(line.y) < 0.57f ? Vector3(0, 1, 0) : Vector3(0, 0, 1));
        Vector3 p0 = normalize(cross(line, axis)), p1 = cross(line, p0);
        for (int k = 0; k < 6; ++k) {
            float t = kPi / 3.0f * k;
            SV w = supportAB(A, B, p0 * std::cos(t) + p1 * std::sin(t));
            Vector3 rel = w.w - V[0].w;
            if (length(rel - line * dot(rel, line)) > eps) { V.push_back(w); break; }
        }
    }
    if (V.size() == 3) {
        Vector3 n = normalize(cross(V[1].w - V[0].w, V[2].w - V[0].w));
        SV w = supportAB(A, B, n);
        if (std::fabs(dot(w.w - V[0].w, n)) > eps) V.push_back(w);
        else {
            w = supportAB(A, B, -n);
            if (std::fabs(dot(w.w - V[0].w, n)) > eps) V.push_back(w);
        }
    }
}

// The four faces of the starting tetrahedron, oriented outwards.
void buildInitialPolytope(EpaScratch& S) {
    const std::vector<SV>& V = S.V;
    std::vector<EpaFace>& faces = S.faces;
    faces.clear();
    const int tf[4][4] = {{0, 1, 2, 3}, {0, 3, 1, 2}, {0, 2, 3, 1}, {1, 3, 2, 0}};
    for (auto& t : tf) {
        int a = t[0], b = t[1], c = t[2];
        Vector3 n = cross(V[b].w - V[a].w, V[c].w - V[a].w);
        if (dot(n, V[t[3]].w - V[a].w) > 0) std::swap(b, c); // outward orientation
        faces.push_back(makeEpaFace(V, a, b, c));
    }
}

// The expansion (van den Bergen 2001): take the face closest to the origin, ask the Minkowski
// difference for its support point in that direction; if the point is not beyond the face the
// polytope has reached the boundary. Otherwise remove the faces seen from the point and join its
// horizon edges to it. Bounded by 96 rounds and 512 faces.
void expandPolytope(const PosedShape& A, const PosedShape& B, EpaScratch& S) {
    std::vector<SV>& V = S.V;
    std::vector<EpaFace>& faces = S.faces;
    for (int iter = 0; iter < 96; ++iter) {
        int closest = 0;
        for (int i = 1; i < int(faces.size()); ++i)
            if (faces[i].d < faces[closest].d) closest = i;
        const EpaFace f = faces[closest];
        if (f.d == kInf) return;
        SV w = supportAB(A, B, f.n);
        float dist = dot(w.w, f.n);
        if (dist - f.d < 1e-5f * std::max(1.0f, std::fabs(dist))) break; // converged
        int wi = int(V.size());
        V.push_back(w);
        // Remove the faces seen from w and collect the horizon (edges used by exactly one of them).
        std::vector<std::pair<int, int>>& horizon = S.horizon;
        std::vector<EpaFace>& kept = S.kept;
        horizon.clear();
        kept.clear();
        for (const EpaFace& fc : faces) {
            if (fc.d != kInf && dot(fc.n, w.w - V[fc.a].w) > 1e-7f) {
                const int e[3][2] = {{fc.a, fc.b}, {fc.b, fc.c}, {fc.c, fc.a}};
                for (auto& ed : e) {
                    auto it = std::find(horizon.begin(), horizon.end(), std::make_pair(ed[1], ed[0]));
                    if (it != horizon.end()) horizon.erase(it);
                    else horizon.emplace_back(ed[0], ed[1]);
                }
            } else {
                kept.push_back(fc);
            }
        }
        if (horizon.empty()) break;
        faces.swap(kept);
        for (auto& ed : horizon) faces.push_back(makeEpaFace(V, ed.first, ed.second, wi));
        if (faces.size() > 512) break;
    }
}

// The penetration from the closest face of the finished polytope: its distance is the depth, its
// normal the direction, and the witness points on A and B come from the same barycentric weights
// (the face's vertices remember which support points of A and B made them).
PenetrationResult penetrationFromPolytope(const EpaScratch& S) {
    PenetrationResult res;
    const std::vector<SV>& V = S.V;
    const std::vector<EpaFace>& faces = S.faces;
    // The loop may end right after rebuilding the polytope (face limit, iteration limit): the
    // index found at the top of the last iteration then points into the old face list.
    int closest = 0;
    for (int i = 1; i < int(faces.size()); ++i)
        if (faces[i].d < faces[closest].d) closest = i;
    if (faces[closest].d == kInf) return res;
    const EpaFace& f = faces[closest];
    // Barycentric coordinates of the origin's projection onto the closest face.
    Vector3 p = f.n * f.d;
    Vector3 a = V[f.a].w, b = V[f.b].w, c = V[f.c].w;
    Vector3 v0 = b - a, v1 = c - a, v2 = p - a;
    float d00 = dot(v0, v0), d01 = dot(v0, v1), d11 = dot(v1, v1), d20 = dot(v2, v0), d21 = dot(v2, v1);
    float den = d00 * d11 - d01 * d01;
    float bv = 0, bw = 0;
    if (std::fabs(den) > 1e-20f) {
        bv = (d11 * d20 - d01 * d21) / den;
        bw = (d00 * d21 - d01 * d20) / den;
    }
    float bu = 1 - bv - bw;
    res.pointA = V[f.a].a * bu + V[f.b].a * bv + V[f.c].a * bw;
    res.pointB = V[f.a].b * bu + V[f.b].b * bv + V[f.c].b * bw;
    // f.n points out of A-B; A must move along -f.n to separate.
    res.normal = -f.n;
    res.depth = std::max(f.d, 0.0f);
    res.valid = true;
    return res;
}

} // namespace

// EPA, the expanding polytope algorithm: from the simplex GJK left inside the Minkowski
// difference, a polytope grows outwards until it touches the difference's boundary; the face it
// touches with gives the penetration depth, the normal and the witness points.
PenetrationResult epa(const PosedShape& A, const PosedShape& B, const GjkResult& g) {
    EpaScratch& S = epaScratch();
    std::vector<SV>& V = S.V;
    V.clear();
    for (int i = 0; i < g.simplexSize; ++i) V.push_back({g.w[i], g.a[i], g.b[i]});
    inflateToTetrahedron(A, B, V);
    if (V.size() < 4) return PenetrationResult(); // touching with zero volume: no meaningful penetration
    buildInitialPolytope(S);
    expandPolytope(A, B, S);
    return penetrationFromPolytope(S);
}

bool penetration(const PosedShape& A, const PosedShape& B, PenetrationResult& out) {
    GjkResult g = gjk(A, B);
    if (!g.intersect) return false;
    out = epa(A, B, g);
    return out.valid;
}

} // namespace rf
