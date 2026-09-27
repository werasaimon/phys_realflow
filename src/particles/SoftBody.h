#pragma once
// Soft body made of particles, held together by shape matching (Mueller, Heidelberger, Teschner,
// Gross 2005, "Meshless deformations based on shape matching") on overlapping clusters, as in
// NVIDIA FleX: every cluster finds the rigid motion that best fits its particles now - the centre
// of mass and the rotation of the deformation (Mueller et al. 2016) - and pulls its particles
// towards where that rigid motion would put them. One big cluster = rigid; many small ones = the
// body can bend and squash, and springs back.
//
// For drawing, the body's surface mesh is skinned to the clusters (as FleX does): every vertex
// follows the rigid motions of the clusters of its nearest particle.

#include "core/Mesh.h"
#include "math/Math.h"

#include <vector>

namespace rf {

struct SoftCluster {
    std::vector<int> particles;        // global particle indices
    std::vector<Vector3> restOffsets;  // rest position - rest centre of mass (same order)
    Vector3 restCentre;                // centre of mass at rest
    Vector3 centre;                    // centre of mass now (last shape matching)
    Quaternion rotation;               // rotation now (also the warm start of the next extraction)
};

struct SoftBody {
    int object = -1;                   // particle object id (for self-collision filtering)
    int group = -1;                    // particle group (ParticleSystem::removeGroup removes it whole)
    std::vector<int> particles;
    std::vector<SoftCluster> clusters;
    // Fraction of the way back to the rest shape per substep (0.05 jelly, 0.3 rubber, 1 rigid).
    // Spread over the solver's passes as k' = 1 - (1 - k)^(1/n) (Mueller et al. 2007, PBD), so the
    // softness does not depend on how many passes a substep makes.
    float stiffness = 0.3f;
    Vector3 color{0.9f, 0.4f, 0.4f};
    float clusterRadius = 0;           // radius the clusters were built with [m] (the skinning blends within it)
    // Surface for drawing: rest mesh and, per vertex, the clusters it follows with their weights
    // (smooth linear blend skinning: every cluster within clusterRadius, weighted by closeness).
    TriMesh surface;
    std::vector<std::vector<int>> vertexClusters;
    std::vector<std::vector<float>> vertexWeights;
};

// Clusters for particles at rest positions `rest` (global indices `ids`): cluster centres on a
// lattice with the given spacing, each taking the particles within `radius` of it.
std::vector<SoftCluster> buildClusters(const std::vector<int>& ids, const std::vector<Vector3>& rest, float spacing,
                                       float radius);

// Binds the rest surface mesh to the clusters (vertex -> clusters of the nearest particle).
void bindSurface(SoftBody& body, const TriMesh& restSurface, const std::vector<Vector3>& particleRest);

// Current positions of the surface vertices: the average of where each of the vertex's clusters
// carries it (cluster centre + rotation of the rest offset).
void skinSurface(const SoftBody& body, std::vector<Vector3>& out);

// One shape-matching pass over all soft bodies: moves the predicted positions p towards the goals
// (average over the clusters of each particle). Particles with invMass 0 stay where they are.
// passesPerStep: how many such passes the substep makes (the stiffness is spread over them).
void solveShapeMatching(std::vector<SoftBody>& bodies, std::vector<Vector3>& p, const std::vector<float>& invMass, int passesPerStep);

} // namespace rf
