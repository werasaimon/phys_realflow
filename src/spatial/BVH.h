#pragma once
// Bounding volume hierarchy (binned SAH) used for
//  * triangle meshes: ray casts, closest point, signed distance (angle-weighted pseudo-normals)
//  * rigid-body broad phase: AABB overlap queries.

#include "math/Math.h"
#include "core/Mesh.h"

#include <vector>

namespace rf {

struct BVHNode {
    AABB box;
    int32_t first = 0;  // interior: index of left child (right = first + 1); leaf: first primitive slot
    int32_t count = 0;  // 0 for interior nodes
    bool isLeaf() const { return count > 0; }
};

class BVH {
public:
    void build(const std::vector<AABB>& primBoxes, int maxLeafSize = 4);
    void clear() { nodes_.clear(); prims_.clear(); }
    bool empty() const { return nodes_.empty(); }

    const std::vector<BVHNode>& nodes() const { return nodes_; }
    const std::vector<uint32_t>& prims() const { return prims_; }
    AABB bounds() const { return nodes_.empty() ? AABB() : nodes_[0].box; }

    // Calls fn(primIndex) for every primitive whose box overlaps q.
    template <class F> void queryAABB(const AABB& q, F&& fn) const {
        if (nodes_.empty()) return;
        int stack[64];
        int sp = 0;
        stack[sp++] = 0;
        while (sp) {
            const BVHNode& n = nodes_[stack[--sp]];
            if (!n.box.overlaps(q)) continue;
            if (n.isLeaf()) {
                for (int i = 0; i < n.count; ++i) fn(prims_[n.first + i]);
            } else {
                stack[sp++] = n.first;
                stack[sp++] = n.first + 1;
            }
        }
    }

private:
    void subdivide(int nodeIdx, const std::vector<AABB>& boxes, const std::vector<Vector3>& centroids, int maxLeaf, int depth);

    std::vector<BVHNode> nodes_;
    std::vector<uint32_t> prims_;
};

struct RayHit {
    float t = kInf;
    int triangle = -1;
    float u = 0, v = 0;
};

struct ClosestHit {
    Vector3 point;
    Vector3 normal;           // pseudo-normal at the closest feature (outward)
    float distance = kInf; // unsigned
    float signedDistance = kInf;
    int triangle = -1;
};

class MeshBVH {
public:
    void build(const TriMesh& mesh);
    void clear();
    bool empty() const { return tris_.empty(); }
    AABB bounds() const { return bvh_.bounds(); }
    const BVH& bvh() const { return bvh_; }
    size_t triangleCount() const { return tris_.size(); }

    void triangle(size_t t, Vector3& a, Vector3& b, Vector3& c) const {
        a = pos_[tris_[t][0]]; b = pos_[tris_[t][1]]; c = pos_[tris_[t][2]];
    }
    const Vector3& faceNormal(size_t t) const { return faceN_[t]; }

    bool raycast(const Vector3& origin, const Vector3& dir, float tmax, RayHit& hit) const;
    // Nearest surface point within maxDist; false if nothing is that close.
    bool closestPoint(const Vector3& p, float maxDist, ClosestHit& out) const;
    // Signed distance (negative inside), clamped to +/-maxDist when farther than maxDist
    // (in which case the sign is +, i.e. assumed outside; use isInside() for a definitive answer).
    float signedDistance(const Vector3& p, float maxDist, Vector3* normal = nullptr) const;
    bool isInside(const Vector3& p) const;

private:
    std::vector<Vector3> pos_;
    std::vector<std::array<uint32_t, 3>> tris_;
    std::vector<Vector3> faceN_;
    std::vector<Vector3> vertexN_;  // angle weighted
    std::vector<Vector3> edgeN_;    // 3 per triangle: ab, bc, ca
    BVH bvh_;
};

} // namespace rf
