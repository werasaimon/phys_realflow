#include "rigid/ConvexDecomposition.h"

#include "spatial/BVH.h"
#include "core/Parallel.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <queue>
#include <unordered_map>

namespace rf {

// ---------------------------------------------------------------------------
// Convex hull: Quickhull (Barber, Dobkin & Huhdanpaa 1996) in double precision
// ---------------------------------------------------------------------------
// Every face keeps the points outside it (conflict list); the farthest outside point is added
// first, the visible region is grown as a connected patch from that face (so the horizon is one
// closed loop), and the orphaned points are handed to the new faces. Adding extreme points first
// keeps the new faces well shaped - this is what makes the hull robust for dense, nearly coplanar
// samples of smooth surfaces.
TriMesh buildConvexHull(const std::vector<Vector3>& input, int maxVertices) {
    struct D3 { double x, y, z; };
    auto sub = [](const D3& a, const D3& b) { return D3{a.x - b.x, a.y - b.y, a.z - b.z}; };
    auto crs = [](const D3& a, const D3& b) { return D3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; };
    auto dt = [](const D3& a, const D3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; };

    TriMesh out;
    if (input.size() < 4) return out;
    AABB bb;
    for (const Vector3& p : input)
        if (std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z)) bb.expand(p);
    if (!bb.valid()) return out;
    const double scale = std::max(double(length(bb.extent())), 1e-9);
    const double maxAbs = std::max({std::fabs(bb.lo.x), std::fabs(bb.lo.y), std::fabs(bb.lo.z), std::fabs(bb.hi.x),
                                    std::fabs(bb.hi.y), std::fabs(bb.hi.z)});
    // Distance tolerance: the coordinates are only float accurate.
    const double eps = std::max(3.0 * (maxAbs + scale) * 1.2e-7, 1e-12);

    // Deduplicate (points closer than ~1e-5 of the size add nothing but degeneracy).
    std::vector<D3> P;
    {
        std::unordered_map<uint64_t, int> seen;
        const double q = 1.0 / (1e-5 * scale);
        auto k = [&](double v) { return uint64_t(uint32_t(int32_t(std::floor(v * q)))) & 0x1FFFFF; };
        for (const Vector3& p : input) {
            if (!(std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z))) continue;
            uint64_t key = (k(p.x) << 42) | (k(p.y) << 21) | k(p.z);
            if (seen.emplace(key, 1).second) P.push_back({p.x, p.y, p.z});
        }
    }
    const int n = int(P.size());
    if (n < 4) return out;

    // Initial tetrahedron from extreme points.
    int i0 = 0;
    for (int i = 1; i < n; ++i)
        if (P[i].x < P[i0].x) i0 = i;
    int i1 = i0;
    double best = -1;
    for (int i = 0; i < n; ++i) {
        D3 v = sub(P[i], P[i0]);
        if (dt(v, v) > best) { best = dt(v, v); i1 = i; }
    }
    int i2 = i0;
    best = -1;
    const D3 e01 = sub(P[i1], P[i0]);
    for (int i = 0; i < n; ++i) {
        D3 c = crs(e01, sub(P[i], P[i0]));
        if (dt(c, c) > best) { best = dt(c, c); i2 = i; }
    }
    const D3 nrm = crs(e01, sub(P[i2], P[i0]));
    const double nl = std::sqrt(dt(nrm, nrm));
    if (nl <= 0) return out;
    int i3 = i0;
    best = -1;
    for (int i = 0; i < n; ++i) {
        double d = std::fabs(dt(sub(P[i], P[i0]), nrm)) / nl;
        if (d > best) { best = d; i3 = i; }
    }
    if (best < 10 * eps) return out; // flat point set

    struct Face {
        int a, b, c;
        D3 n;
        double d;
        bool alive = true;
        std::vector<int> outside;
        int farthest = -1;
        double farDist = 0;
    };
    std::vector<Face> faces;
    auto dist = [&](const Face& f, int pi) { return dt(f.n, P[pi]) - f.d; };
    auto ek = [](int a, int b) { return (uint64_t(uint32_t(a)) << 32) | uint32_t(b); };
    std::unordered_map<uint64_t, int> edgeFace; // directed edge -> alive face owning it
    auto addFace = [&](int a, int b, int c) {
        Face f;
        f.a = a, f.b = b, f.c = c;
        D3 cr = crs(sub(P[b], P[a]), sub(P[c], P[a]));
        double l = std::sqrt(dt(cr, cr));
        f.n = l > 0 ? D3{cr.x / l, cr.y / l, cr.z / l} : D3{0, 0, 0};
        f.d = dt(f.n, P[a]);
        faces.push_back(std::move(f));
        int fi = int(faces.size()) - 1;
        edgeFace[ek(a, b)] = fi;
        edgeFace[ek(b, c)] = fi;
        edgeFace[ek(c, a)] = fi;
        return fi;
    };
    auto assign = [&](int pi, const std::vector<int>& cand) {
        for (int fi : cand) {
            double d = dist(faces[fi], pi);
            if (d > eps) {
                Face& f = faces[fi];
                f.outside.push_back(pi);
                if (d > f.farDist) { f.farDist = d; f.farthest = pi; }
                return;
            }
        }
    };
    {
        const D3 centroid{(P[i0].x + P[i1].x + P[i2].x + P[i3].x) * 0.25, (P[i0].y + P[i1].y + P[i2].y + P[i3].y) * 0.25,
                          (P[i0].z + P[i1].z + P[i2].z + P[i3].z) * 0.25};
        const int tet[4][3] = {{i0, i1, i2}, {i0, i1, i3}, {i0, i2, i3}, {i1, i2, i3}};
        std::vector<int> init;
        for (auto& t : tet) {
            int a = t[0], b = t[1], c = t[2];
            if (dt(crs(sub(P[b], P[a]), sub(P[c], P[a])), sub(P[a], centroid)) < 0) std::swap(b, c);
            init.push_back(addFace(a, b, c));
        }
        for (int pi = 0; pi < n; ++pi)
            if (pi != i0 && pi != i1 && pi != i2 && pi != i3) assign(pi, init);
    }

    // Globally farthest point first: with a vertex limit the hull stops as the best approximation
    // of that size (greedy reduction of the Hausdorff distance, as PhysX's vertex-limited cooking).
    std::priority_queue<std::pair<double, int>> pending;
    for (int fi = 0; fi < 4; ++fi)
        if (faces[fi].farthest >= 0) pending.push({faces[fi].farDist, fi});
    std::vector<int> mark;
    int iter = 0, vertices = 4;
    while (!pending.empty()) {
        const int fi = pending.top().second;
        pending.pop();
        if (!faces[fi].alive || faces[fi].farthest < 0) continue;
        if (maxVertices >= 4 && vertices >= maxVertices) break;
        ++vertices;
        const int pi = faces[fi].farthest;
        ++iter;
        if (mark.size() < faces.size()) mark.resize(faces.size() * 2, -1);

        // Connected visible region.
        std::vector<int> region = {fi}, stack = {fi};
        mark[fi] = iter;
        while (!stack.empty()) {
            const int f = stack.back();
            stack.pop_back();
            const int e[3][2] = {{faces[f].a, faces[f].b}, {faces[f].b, faces[f].c}, {faces[f].c, faces[f].a}};
            for (auto& ed : e) {
                auto it = edgeFace.find(ek(ed[1], ed[0]));
                if (it == edgeFace.end()) continue;
                const int nb = it->second;
                if (mark[nb] == iter || dist(faces[nb], pi) <= eps) continue;
                mark[nb] = iter;
                region.push_back(nb);
                stack.push_back(nb);
            }
        }
        std::vector<std::pair<int, int>> horizon;
        for (int f : region) {
            const int e[3][2] = {{faces[f].a, faces[f].b}, {faces[f].b, faces[f].c}, {faces[f].c, faces[f].a}};
            for (auto& ed : e) {
                auto it = edgeFace.find(ek(ed[1], ed[0]));
                if (it == edgeFace.end() || mark[it->second] != iter) horizon.emplace_back(ed[0], ed[1]);
            }
        }
        std::vector<int> orphans;
        for (int f : region) {
            Face& F = faces[f];
            F.alive = false;
            edgeFace.erase(ek(F.a, F.b));
            edgeFace.erase(ek(F.b, F.c));
            edgeFace.erase(ek(F.c, F.a));
            for (int q : F.outside)
                if (q != pi) orphans.push_back(q);
            std::vector<int>().swap(F.outside);
        }
        std::vector<int> created;
        for (auto& ed : horizon) created.push_back(addFace(ed.first, ed.second, pi));
        for (int q : orphans) assign(q, created);
        for (int f : created)
            if (faces[f].farthest >= 0) pending.push({faces[f].farDist, f});
    }

    std::unordered_map<int, uint32_t> remap;
    for (const Face& f : faces) {
        if (!f.alive) continue;
        std::array<uint32_t, 3> t;
        const int idx[3] = {f.a, f.b, f.c};
        for (int k = 0; k < 3; ++k) {
            auto it = remap.find(idx[k]);
            if (it == remap.end()) {
                it = remap.emplace(idx[k], uint32_t(out.positions.size())).first;
                out.positions.push_back(Vector3(float(P[idx[k]].x), float(P[idx[k]].y), float(P[idx[k]].z)));
            }
            t[k] = it->second;
        }
        out.triangles.push_back(t);
    }
    return out;
}

// ---------------------------------------------------------------------------
// Voxel-based approximate convex decomposition
// ---------------------------------------------------------------------------
std::vector<TriMesh> convexDecomposition(const std::vector<TriMesh>& parts, const DecompositionParams& prm) {
    std::vector<TriMesh> result;
    AABB bb;
    for (const TriMesh& m : parts) bb.expand(m.bounds());
    if (!bb.valid()) return result;
    const float h = maxComp(bb.extent()) / float(std::max(prm.resolution, 4));
    const int nx = std::max(1, int(std::ceil(bb.extent().x / h))), ny = std::max(1, int(std::ceil(bb.extent().y / h))),
              nz = std::max(1, int(std::ceil(bb.extent().z / h)));

    std::vector<MeshBVH> bvhs(parts.size());
    for (size_t i = 0; i < parts.size(); ++i) bvhs[i].build(parts[i]);
    std::vector<uint8_t> occ(size_t(nx) * ny * nz, 0);
    parallelFor(nz, [&](int k) {
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                Vector3 c = bb.lo + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * h;
                for (size_t m = 0; m < bvhs.size(); ++m)
                    if (bvhs[m].bounds().contains(c) && bvhs[m].isInside(c)) { occ[size_t(i) + size_t(nx) * (j + size_t(ny) * k)] = 1; break; }
            }
    }, 1);

    struct V { int i, j, k; };
    std::vector<V> all;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                if (occ[size_t(i) + size_t(nx) * (j + size_t(ny) * k)]) all.push_back({i, j, k});
    if (all.empty()) return result;
    const float totalVolume = float(all.size()) * h * h * h;

    // Hull of a voxel set: only voxels with a face exposed within the set can contribute extreme
    // corners, and shared corners are emitted once.
    std::vector<int> stamp(occ.size(), 0);
    int stampId = 0;
    auto hullOf = [&](const std::vector<V>& vs) {
        ++stampId;
        auto id = [&](int i, int j, int k) { return size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * size_t(k)); };
        for (const V& v : vs) stamp[id(v.i, v.j, v.k)] = stampId;
        auto member = [&](int i, int j, int k) {
            return i >= 0 && j >= 0 && k >= 0 && i < nx && j < ny && k < nz && stamp[id(i, j, k)] == stampId;
        };
        std::vector<Vector3> pts;
        for (const V& v : vs) {
            if (member(v.i + 1, v.j, v.k) && member(v.i - 1, v.j, v.k) && member(v.i, v.j + 1, v.k) &&
                member(v.i, v.j - 1, v.k) && member(v.i, v.j, v.k + 1) && member(v.i, v.j, v.k - 1))
                continue;
            for (int c = 0; c < 8; ++c)
                pts.push_back(bb.lo + Vector3(float(v.i + (c & 1)), float(v.j + ((c >> 1) & 1)), float(v.k + ((c >> 2) & 1))) * h);
        }
        return buildConvexHull(pts);
    };

    // 1) Recursive splitting of the voxel set into nearly convex regions.
    std::vector<std::vector<V>> leaves;
    std::function<void(std::vector<V>&, int)> split = [&](std::vector<V>& vs, int depth) {
        TriMesh hull = hullOf(vs);
        float vox = float(vs.size()) * h * h * h;
        float excess = hull.signedVolume() - vox;
        bool leaf = depth >= prm.maxDepth || vs.size() < 8 || excess <= prm.concavity * totalVolume ||
                    int(leaves.size()) >= prm.maxParts * 2;
        if (!leaf) {
            // Best axis-aligned cut: minimise the empty space in the bounding boxes of the halves.
            int bestAxis = -1, bestCut = 0;
            float bestCost = kInf;
            for (int axis = 0; axis < 3; ++axis) {
                int lo = 1 << 30, hi = -(1 << 30);
                for (const V& v : vs) {
                    int c = axis == 0 ? v.i : axis == 1 ? v.j : v.k;
                    lo = std::min(lo, c);
                    hi = std::max(hi, c);
                }
                for (int cut = lo + 1; cut <= hi; ++cut) {
                    AABB a, b;
                    int na = 0, nb = 0;
                    for (const V& v : vs) {
                        int c = axis == 0 ? v.i : axis == 1 ? v.j : v.k;
                        Vector3 p(float(v.i), float(v.j), float(v.k));
                        if (c < cut) { a.expand(p); ++na; }
                        else { b.expand(p); ++nb; }
                    }
                    if (!na || !nb) continue;
                    auto vol = [](const AABB& x) { Vector3 e = x.extent() + Vector3(1.0f); return e.x * e.y * e.z; };
                    float cost = vol(a) - na + vol(b) - nb;
                    if (cost < bestCost) { bestCost = cost; bestAxis = axis; bestCut = cut; }
                }
            }
            if (bestAxis >= 0) {
                std::vector<V> a, b;
                for (const V& v : vs) {
                    int c = bestAxis == 0 ? v.i : bestAxis == 1 ? v.j : v.k;
                    (c < bestCut ? a : b).push_back(v);
                }
                split(a, depth + 1);
                split(b, depth + 1);
                return;
            }
        }
        leaves.push_back(vs);
    };
    split(all, 0);

    // 2) Fit to the real surface: each part's hull is built from the points of the original smooth
    //    surface inside its region, plus the corners of fully interior voxels (inner cut faces).
    auto cellIdx = [&](int i, int j, int k) { return size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * size_t(k)); };
    std::vector<int> owner(occ.size(), -1);
    for (int p = 0; p < int(leaves.size()); ++p)
        for (const V& v : leaves[p]) owner[cellIdx(v.i, v.j, v.k)] = p;
    std::vector<std::vector<Vector3>> pts(leaves.size());
    auto ownerOf = [&](const Vector3& x) {
        Vector3 g = (x - bb.lo) / h;
        int i = clampv(int(g.x), 0, nx - 1), j = clampv(int(g.y), 0, ny - 1), k = clampv(int(g.z), 0, nz - 1);
        if (owner[cellIdx(i, j, k)] >= 0) return owner[cellIdx(i, j, k)];
        int best = -1;
        float bd = kInf;
        for (int dz = -1; dz <= 1; ++dz)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx) {
                    int a = i + dx, b = j + dy, c = k + dz;
                    if (a < 0 || b < 0 || c < 0 || a >= nx || b >= ny || c >= nz) continue;
                    int o = owner[cellIdx(a, b, c)];
                    if (o < 0) continue;
                    float d = length2(bb.lo + Vector3(a + 0.5f, b + 0.5f, c + 0.5f) * h - x);
                    if (d < bd) { bd = d; best = o; }
                }
        return best;
    };
    // Dense surface samples (spacing ~ h/3) of every part mesh.
    for (const TriMesh& m : parts)
        for (const auto& t : m.triangles) {
            const Vector3 &a = m.positions[t[0]], &b = m.positions[t[1]], &c = m.positions[t[2]];
            float L = std::max({length(b - a), length(c - a), length(c - b)});
            int n = std::max(1, int(std::ceil(L / (h / 2.0f))));
            for (int u = 0; u <= n; ++u)
                for (int v = 0; u + v <= n; ++v) {
                    Vector3 x = a + (b - a) * (float(u) / n) + (c - a) * (float(v) / n);
                    // Keep only points on the outer surface of the union (not inside another part).
                    bool buried = false;
                    for (size_t q = 0; q < bvhs.size() && !buried; ++q)
                        if (&parts[q] != &m && bvhs[q].bounds().contains(x) && bvhs[q].signedDistance(x, h) < -0.25f * h) buried = true;
                    if (buried) continue;
                    int o = ownerOf(x);
                    if (o >= 0) pts[o].push_back(x);
                }
        }
    for (int p = 0; p < int(leaves.size()); ++p)
        for (const V& v : leaves[p]) {
            bool interior = true;
            const int nb[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
            for (auto& d : nb) {
                int a = v.i + d[0], b = v.j + d[1], c = v.k + d[2];
                if (a < 0 || b < 0 || c < 0 || a >= nx || b >= ny || c >= nz || !occ[cellIdx(a, b, c)]) { interior = false; break; }
            }
            if (!interior) continue;
            for (int c = 0; c < 8; ++c)
                pts[p].push_back(bb.lo + Vector3(float(v.i + (c & 1)), float(v.j + ((c >> 1) & 1)), float(v.k + ((c >> 2) & 1))) * h);
        }

    // 3) Minimal set: greedily merge the pair whose merged hull adds the least volume, while the
    //    excess of the merged part stays within the concavity budget.
    struct Part { std::vector<Vector3> pts; float solidVolume; TriMesh hull; };
    std::vector<Part> P;
    for (int p = 0; p < int(leaves.size()); ++p) {
        if (pts[p].size() < 4) continue;
        Part part{{}, float(leaves[p].size()) * h * h * h, buildConvexHull(pts[p])};
        if (part.hull.empty()) continue;
        part.pts = part.hull.positions; // hull(A u B) = hull(verts(hull A) u verts(hull B))
        P.push_back(std::move(part));
    }
    const float budget = prm.concavity * totalVolume;
    for (;;) {
        int bi = -1, bj = -1;
        float bestExcess = kInf;
        TriMesh bestHull;
        for (int i = 0; i < int(P.size()); ++i)
            for (int j = i + 1; j < int(P.size()); ++j) {
                AABB a = P[i].hull.bounds(), b = P[j].hull.bounds();
                a.lo -= Vector3(h); a.hi += Vector3(h);
                if (!a.overlaps(b)) continue; // only neighbours
                std::vector<Vector3> u = P[i].pts;
                u.insert(u.end(), P[j].pts.begin(), P[j].pts.end());
                TriMesh hu = buildConvexHull(u);
                float excess = hu.signedVolume() - (P[i].solidVolume + P[j].solidVolume);
                if (excess < bestExcess) { bestExcess = excess; bi = i; bj = j; bestHull = std::move(hu); }
            }
        bool mustMerge = int(P.size()) > prm.maxParts;
        if (bi < 0 || (!mustMerge && bestExcess > budget)) break;
        P[bi].solidVolume += P[bj].solidVolume;
        P[bi].hull = std::move(bestHull);
        P[bi].pts = P[bi].hull.positions;
        P.erase(P.begin() + bj);
    }
    for (Part& part : P) {
        if (prm.maxHullVertices >= 4 && int(part.hull.positions.size()) > prm.maxHullVertices)
            part.hull = buildConvexHull(part.hull.positions, prm.maxHullVertices);
        result.push_back(std::move(part.hull));
    }
    return result;
}

} // namespace rf
