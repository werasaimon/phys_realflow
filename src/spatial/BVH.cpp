#include "spatial/BVH.h"

#include <algorithm>
#include <numeric>
#include <unordered_map>

namespace rf {

// ---------------------------------------------------------------------------
// Generic BVH
// ---------------------------------------------------------------------------
void BVH::build(const std::vector<AABB>& boxes, int maxLeaf) {
    nodes_.clear();
    prims_.resize(boxes.size());
    std::iota(prims_.begin(), prims_.end(), 0u);
    if (boxes.empty()) return;
    std::vector<Vector3> centroids(boxes.size());
    for (size_t i = 0; i < boxes.size(); ++i) centroids[i] = boxes[i].center();
    nodes_.reserve(boxes.size() * 2);
    BVHNode root;
    root.first = 0;
    root.count = int32_t(boxes.size());
    nodes_.push_back(root);
    subdivide(0, boxes, centroids, std::max(1, maxLeaf), 0);
}

void BVH::subdivide(int nodeIdx, const std::vector<AABB>& boxes, const std::vector<Vector3>& centroids, int maxLeaf,
                    int depth) {
    const int first = nodes_[nodeIdx].first;
    const int count = nodes_[nodeIdx].count;
    AABB nodeBox, cBox;
    for (int i = first; i < first + count; ++i) {
        nodeBox.expand(boxes[prims_[i]]);
        cBox.expand(centroids[prims_[i]]);
    }
    nodes_[nodeIdx].box = nodeBox;
    if (count <= maxLeaf || depth >= 48) return;

    constexpr int kBins = 12;
    // Bin of a centroid coordinate. The float is range-checked before the int conversion, so a NaN
    // or infinite one (a body that blew up, a broken mesh) lands in an end bin instead of indexing
    // with int(NaN).
    auto binOf = [&](float c, int axis, float scale) {
        const float f = (c - cBox.lo[axis]) * scale;
        if (!(f > 0.0f)) return 0;
        return f < float(kBins - 1) ? int(f) : kBins - 1;
    };
    Vector3 cExt = cBox.extent();
    int bestAxis = -1, bestSplit = -1;
    float bestCost = kInf;
    for (int axis = 0; axis < 3; ++axis) {
        if (cExt[axis] <= 1e-12f) continue;
        AABB binBox[kBins];
        int binCount[kBins] = {};
        float scale = kBins / cExt[axis];
        for (int i = first; i < first + count; ++i) {
            int b = binOf(centroids[prims_[i]][axis], axis, scale);
            binCount[b]++;
            binBox[b].expand(boxes[prims_[i]]);
        }
        float leftArea[kBins - 1];
        int leftCount[kBins - 1];
        AABB acc;
        int n = 0;
        for (int b = 0; b < kBins - 1; ++b) {
            acc.expand(binBox[b]);
            n += binCount[b];
            leftArea[b] = acc.valid() ? acc.surfaceArea() : 0.0f;
            leftCount[b] = n;
        }
        acc = AABB();
        n = 0;
        for (int b = kBins - 1; b > 0; --b) {
            acc.expand(binBox[b]);
            n += binCount[b];
            float cost = leftCount[b - 1] * leftArea[b - 1] + n * (acc.valid() ? acc.surfaceArea() : 0.0f);
            if (leftCount[b - 1] > 0 && n > 0 && cost < bestCost) {
                bestCost = cost;
                bestAxis = axis;
                bestSplit = b;
            }
        }
    }

    int mid;
    if (bestAxis < 0) {
        if (count <= 16) return; // all centroids coincide
        mid = first + count / 2;
    } else {
        float leafCost = count * nodeBox.surfaceArea();
        if (bestCost >= leafCost && count <= 16) return;
        float scale = kBins / cExt[bestAxis];
        auto it = std::partition(prims_.begin() + first, prims_.begin() + first + count, [&](uint32_t p) {
            return binOf(centroids[p][bestAxis], bestAxis, scale) < bestSplit;
        });
        mid = int(it - prims_.begin());
        if (mid == first || mid == first + count) mid = first + count / 2;
    }

    int left = int(nodes_.size());
    BVHNode l, r;
    l.first = first;
    l.count = mid - first;
    r.first = mid;
    r.count = first + count - mid;
    nodes_.push_back(l);
    nodes_.push_back(r);
    nodes_[nodeIdx].first = left;
    nodes_[nodeIdx].count = 0;
    subdivide(left, boxes, centroids, maxLeaf, depth + 1);
    subdivide(left + 1, boxes, centroids, maxLeaf, depth + 1);
}

// ---------------------------------------------------------------------------
// Mesh BVH
// ---------------------------------------------------------------------------
void MeshBVH::clear() {
    pos_.clear();
    tris_.clear();
    faceN_.clear();
    vertexN_.clear();
    edgeN_.clear();
    bvh_.clear();
}

void MeshBVH::build(const TriMesh& mesh) {
    clear();
    pos_ = mesh.positions;
    for (const auto& t : mesh.triangles) {
        Vector3 c = cross(pos_[t[1]] - pos_[t[0]], pos_[t[2]] - pos_[t[0]]);
        if (length2(c) <= 1e-30f) continue;
        tris_.push_back(t);
        faceN_.push_back(normalize(c));
    }
    vertexN_.assign(pos_.size(), Vector3(0.0f));
    std::unordered_map<uint64_t, Vector3> edgeSum;
    auto ekey = [](uint32_t a, uint32_t b) {
        if (a > b) std::swap(a, b);
        return (uint64_t(a) << 32) | b;
    };
    std::vector<AABB> boxes(tris_.size());
    for (size_t t = 0; t < tris_.size(); ++t) {
        const auto& tri = tris_[t];
        for (int c = 0; c < 3; ++c) {
            Vector3 e1 = normalize(pos_[tri[(c + 1) % 3]] - pos_[tri[c]]);
            Vector3 e2 = normalize(pos_[tri[(c + 2) % 3]] - pos_[tri[c]]);
            vertexN_[tri[c]] += faceN_[t] * std::acos(clampv(dot(e1, e2), -1.0f, 1.0f));
            edgeSum[ekey(tri[c], tri[(c + 1) % 3])] += faceN_[t];
            boxes[t].expand(pos_[tri[c]]);
        }
    }
    for (Vector3& n : vertexN_) n = normalize(n);
    edgeN_.resize(tris_.size() * 3);
    for (size_t t = 0; t < tris_.size(); ++t)
        for (int c = 0; c < 3; ++c) edgeN_[t * 3 + c] = normalize(edgeSum[ekey(tris_[t][c], tris_[t][(c + 1) % 3])]);
    bvh_.build(boxes, 4);
}

bool MeshBVH::raycast(const Vector3& o, const Vector3& d, float tmax, RayHit& hit) const {
    if (tris_.empty()) return false;
    Vector3 inv(1.0f / (std::fabs(d.x) > 1e-12f ? d.x : 1e-12f), 1.0f / (std::fabs(d.y) > 1e-12f ? d.y : 1e-12f),
             1.0f / (std::fabs(d.z) > 1e-12f ? d.z : 1e-12f));
    hit.t = tmax;
    hit.triangle = -1;
    const auto& nodes = bvh_.nodes();
    const auto& prims = bvh_.prims();
    int stack[64];
    int sp = 0;
    stack[sp++] = 0;
    while (sp) {
        const BVHNode& n = nodes[stack[--sp]];
        if (n.box.rayHit(o, inv, hit.t) == kInf) continue;
        if (n.isLeaf()) {
            for (int i = 0; i < n.count; ++i) {
                uint32_t t = prims[n.first + i];
                const Vector3& a = pos_[tris_[t][0]];
                Vector3 e1 = pos_[tris_[t][1]] - a, e2 = pos_[tris_[t][2]] - a;
                Vector3 pv = cross(d, e2);
                float det = dot(e1, pv);
                if (std::fabs(det) < 1e-14f) continue;
                float id = 1.0f / det;
                Vector3 tv = o - a;
                float u = dot(tv, pv) * id;
                if (u < 0 || u > 1) continue;
                Vector3 qv = cross(tv, e1);
                float v = dot(d, qv) * id;
                if (v < 0 || u + v > 1) continue;
                float tt = dot(e2, qv) * id;
                if (tt > 1e-7f && tt < hit.t) {
                    hit.t = tt;
                    hit.triangle = int(t);
                    hit.u = u;
                    hit.v = v;
                }
            }
        } else {
            stack[sp++] = n.first;
            stack[sp++] = n.first + 1;
        }
    }
    return hit.triangle >= 0;
}

// Ericson, Real-Time Collision Detection 5.1.5. region: 0 face, 1..3 vertex a,b,c, 4..6 edge ab,bc,ca.
static Vector3 closestPtTri(const Vector3& p, const Vector3& a, const Vector3& b, const Vector3& c, int& region) {
    Vector3 ab = b - a, ac = c - a, ap = p - a;
    float d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) { region = 1; return a; }
    Vector3 bp = p - b;
    float d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) { region = 2; return b; }
    float vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) { region = 4; return a + ab * (d1 / (d1 - d3)); }
    Vector3 cp = p - c;
    float d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) { region = 3; return c; }
    float vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) { region = 6; return a + ac * (d2 / (d2 - d6)); }
    float va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        region = 5;
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    float denom = 1.0f / (va + vb + vc);
    region = 0;
    return a + ab * (vb * denom) + ac * (vc * denom);
}

bool MeshBVH::closestPoint(const Vector3& p, float maxDist, ClosestHit& out) const {
    if (tris_.empty()) return false;
    float best2 = maxDist == kInf ? kInf : maxDist * maxDist;
    int bestTri = -1, bestRegion = 0;
    Vector3 bestPt;
    const auto& nodes = bvh_.nodes();
    const auto& prims = bvh_.prims();
    int stack[64];
    int sp = 0;
    stack[sp++] = 0;
    while (sp) {
        const BVHNode& n = nodes[stack[--sp]];
        if (n.box.distance2(p) > best2) continue;
        if (n.isLeaf()) {
            for (int i = 0; i < n.count; ++i) {
                uint32_t t = prims[n.first + i];
                int region;
                Vector3 q = closestPtTri(p, pos_[tris_[t][0]], pos_[tris_[t][1]], pos_[tris_[t][2]], region);
                float d2 = length2(p - q);
                if (d2 < best2) {
                    best2 = d2;
                    bestTri = int(t);
                    bestRegion = region;
                    bestPt = q;
                }
            }
        } else {
            // Visit the nearer child first (pushed last).
            float dl = nodes[n.first].box.distance2(p);
            float dr = nodes[n.first + 1].box.distance2(p);
            if (dl < dr) {
                stack[sp++] = n.first + 1;
                stack[sp++] = n.first;
            } else {
                stack[sp++] = n.first;
                stack[sp++] = n.first + 1;
            }
        }
    }
    if (bestTri < 0) return false;
    Vector3 N;
    const auto& tri = tris_[bestTri];
    switch (bestRegion) {
    case 0: N = faceN_[bestTri]; break;
    case 1: case 2: case 3: N = vertexN_[tri[bestRegion - 1]]; break;
    default: N = edgeN_[size_t(bestTri) * 3 + (bestRegion - 4)]; break;
    }
    out.point = bestPt;
    out.triangle = bestTri;
    out.distance = std::sqrt(best2);
    Vector3 diff = p - bestPt;
    bool inside = dot(diff, N) < 0;
    out.signedDistance = inside ? -out.distance : out.distance;
    // Outward direction: along the offset where it is well defined, otherwise the pseudo-normal.
    if (out.distance > 1e-6f) out.normal = inside ? -diff / out.distance : diff / out.distance;
    else out.normal = N;
    return true;
}

float MeshBVH::signedDistance(const Vector3& p, float maxDist, Vector3* normal) const {
    ClosestHit h;
    if (!closestPoint(p, maxDist, h)) return maxDist;
    if (normal) *normal = h.normal;
    return h.signedDistance;
}

bool MeshBVH::isInside(const Vector3& p) const {
    ClosestHit h;
    if (!closestPoint(p, kInf, h)) return false;
    return h.signedDistance < 0;
}

} // namespace rf
