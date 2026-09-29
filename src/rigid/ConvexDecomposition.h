#pragma once
// Convex hulls and approximate convex decomposition of non-convex solids.
//
//  * buildConvexHull: Quickhull in double precision, optional vertex limit (farthest first).
//  * convexDecomposition: in the spirit of V-HACD and CoACD - the solid is voxelised (inside test
//    through MeshBVH), then the least convex region is split by an axis-aligned plane, again and
//    again, until every part is nearly convex: its hull reaches no deeper than `gap` into empty
//    space (CoACD's collision-aware concavity) and adds at most `concavity` of the volume
//    (V-HACD's). The hulls are fitted to the smooth surface and neighbours merged while the
//    merged part stays nearly convex.

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
    float gap = 0.01f;        // how far a hull may reach into empty space, as a fraction of the diagonal
    int maxDepth = 16;
};

// `solidParts`: closed meshes whose union is the solid (overlaps allowed).
std::vector<TriMesh> convexDecomposition(const std::vector<TriMesh>& solidParts, const DecompositionParams& params = {});

} // namespace rf
