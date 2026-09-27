#pragma once
// Breaking a body into pieces: the Voronoi cells of a convex hull. Seed points are scattered
// inside the body; the cell of a seed is everything closer to it than to any other seed - the
// hull cut by the bisector planes to all the other seeds, so every cell is convex and closed,
// and the cells fill the body exactly. This is how Blender's Cell Fracture, PhysX Blast and
// Unreal's Chaos pre-fracture a mesh (Muller, Chentanez & Kim 2013, "Real time dynamic fracture
// with volumetric approximate convex decompositions", use the same cells at run time). Each cell
// becomes a convex body (RigidWorld::addConvex); gluing the cells with breakable joints so that
// the body shatters on impact is the run-time part, built on top of this file.
#include "core/Mesh.h"
#include "math/Math.h"

#include <cstdint>
#include <vector>

namespace rf {

// The Voronoi cells of the seeds inside a closed convex mesh (outward oriented). Every cell is
// a closed, convex, outward-oriented triangle mesh with welded vertices; a seed outside the
// hull gets no cell (it is skipped), so the result can be shorter than `seeds`.
std::vector<TriMesh> voronoiCells(const TriMesh& convexHull, const std::vector<Vector3>& seeds);

// `count` seeds spread evenly through the hull (uniform in the bounding box, kept when inside).
std::vector<Vector3> uniformSeeds(const TriMesh& convexHull, int count, uint32_t seed);

// `count` seeds crowded around an impact point: the distance from it is radius * u^2 for a
// uniform u, so the pieces are small where the body was hit and large far away - the pattern
// of brittle fracture (glass, stone). Seeds that would fall outside the hull are redrawn.
std::vector<Vector3> impactSeeds(const TriMesh& convexHull, const Vector3& impactPoint, int count, float radius, uint32_t seed);

// True when the point is on the inner side of every face of a convex, outward-oriented mesh.
bool insideConvex(const TriMesh& convexHull, const Vector3& point, float tolerance = 0.0f);

} // namespace rf
