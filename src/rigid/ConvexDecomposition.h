#pragma once
// Convex hulls and approximate convex decomposition of non-convex solids.
//
//  * buildConvexHull: Quickhull in double precision, optional vertex limit (farthest first).
//  * convexDecomposition: in the spirit of V-HACD - the solid is voxelised (inside test through
//    MeshBVH), then the voxel set is split recursively by axis-aligned planes until every part is
//    nearly convex: (volume of its hull - its voxel volume) <= concavity * total volume. Each part
//    becomes the convex hull of its voxels' corners.

#include "core/Mesh.h"

#include <vector>

namespace rf {

// maxVertices > 0 stops after that many vertices (farthest points first).
TriMesh buildConvexHull(const std::vector<Vector3>& points, int maxVertices = 0);

struct DecompositionParams {
    int resolution = 32;      // voxels along the largest extent
    int maxParts = 24;
    int maxHullVertices = 64; // per convex part (collision cost of support mapping)
    float concavity = 0.015f; // allowed (hull - part) volume, as a fraction of the total volume
    int maxDepth = 10;
};

// `solidParts`: closed meshes whose union is the solid (overlaps allowed).
std::vector<TriMesh> convexDecomposition(const std::vector<TriMesh>& solidParts, const DecompositionParams& params = {});

} // namespace rf
