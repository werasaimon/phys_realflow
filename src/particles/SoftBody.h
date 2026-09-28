#pragma once
// Soft body made of particles on a lattice inside its mesh. Two models hold the particles together:
//
//  * NeoHookean (the default, SoftTets.cpp) - an elastic solid, as Houdini Vellum and the FEM codes
//    do it: the lattice's cubes are cut into tetrahedra, and every tetrahedron stores the energy of
//    a stable Neo-Hookean material (Smith, de Goes, Kim 2018), solved as two XPBD constraints per
//    tetrahedron (Macklin & Mueller 2021) in small steps (Macklin et al. 2019). The material is
//    given by numbers from a handbook: Young's modulus E, Poisson's ratio nu, density, friction.
//  * ShapeMatching (legacy, kept for comparison) - overlapping clusters pulled towards their
//    best-fit rigid pose (Mueller, Heidelberger, Teschner, Gross 2005, "Meshless deformations based
//    on shape matching"), as in NVIDIA FleX. Its "stiffness" is a solver fraction, not a material.
//
// For drawing, the body's surface mesh follows the particles: with tetrahedra every vertex keeps
// its barycentric coordinates in one tetrahedron (an embedding); with clusters it is skinned to
// the clusters' best-fit linear deformations (Mueller et al. 2005, sec. 4.3).

#include "core/Mesh.h"
#include "math/Math.h"

#include <vector>

namespace rf {

// How a soft body keeps its shape (see the head of this file).
enum class SoftModel : uint8_t { NeoHookean, ShapeMatching };

// A soft body's material, in the units a materials handbook uses.
struct SoftMaterial {
    SoftModel model = SoftModel::NeoHookean;
    float density = 1050.0f;     // [kg/m^3]
    float youngModulus = 1.5e4f; // E [Pa]: the stress that stretches it by 1 % is E / 100
    float poissonRatio = 0.45f;  // nu: how much it narrows when stretched; 0.5 keeps its volume exactly
    float friction = 0.5f;       // Coulomb coefficient: floors, bodies, other soft bodies
    // Rayleigh damping, the stiffness-proportional part (a damping matrix beta K): the time beta
    // [s]. A vibration of angular frequency w loses the share beta w / 2 of its amplitude per
    // radian - the body's fast jitter dies, its slow swing hardly notices. 0 = none.
    float damping = 1e-3f;
    float stiffness = 0.3f; // ShapeMatching only: the fraction of the way back to the rest shape per substep
};

// The Lame parameters of the material (Landau & Lifshitz, "Theory of Elasticity", sec. 5):
// shear modulus mu = E / (2 (1 + nu)), lambda_LE = E nu / ((1 + nu)(1 - 2 nu)). nu is kept in
// 0.05 .. 0.495: at 0.5 lambda is infinite (a truly incompressible solid). `lambda` is the one the
// stable Neo-Hookean energy takes, lambda_LE + mu: that energy, expanded for small strains, is
// linear elasticity with mu and lambda - mu (Smith, de Goes, Kim 2018, sec. 3.4, where the same
// shift of the parameters is made) - with lambda_LE alone a beam of nu = 0.3 came out 10 % soft.
struct Lame {
    float mu = 0, lambda = 0;
};
Lame lameParameters(const SoftMaterial& m);

// Young's modulus for an old scene's shape-matching "stiffness" 0..1 (files from before the
// tetrahedra): log-linear through 0.05 -> 5.6 kPa (a loose jelly), 0.3 -> 100 kPa (foam rubber),
// 0.5 -> 1 MPa (rubber), 0.7 and above -> 10 MPa (stiffer is beyond the small steps), so an old
// scene keeps its character: a stack of old "0.5" barrels stands as it did.
float youngFromStiffness(float stiffness);

// One tetrahedron of the elastic model: four particles, the inverse of its rest edge matrix
// D_m = [x1 - x0, x2 - x0, x3 - x0], its rest volume V, and the XPBD multipliers of its two
// constraints (deviatoric and hydrostatic) - their sum over a step is the force times dt^2.
struct SoftTet {
    int v[4] = {0, 0, 0, 0};
    Matrix3x3 restInverse;
    float restVolume = 0;
    float lambdaD = 0, lambdaH = 0;
};

// A skin vertex embedded in a tetrahedron: x = sum_k weight[k] x_{tet.v[k]} (barycentric; a vertex
// outside the tetrahedra - the half spacing between the particles and the mesh - extrapolates).
struct SkinBinding {
    int tet = 0;
    float weight[4] = {1, 0, 0, 0};
};

struct SoftCluster {
    std::vector<int> particles;        // global particle indices
    std::vector<Vector3> restOffsets;  // rest position - rest centre of mass (same order)
    Vector3 restCentre;                // centre of mass at rest
    Vector3 centre;                    // centre of mass now (last shape matching)
    Quaternion rotation;               // rotation now (also the warm start of the next extraction)
    // The best-fit linear map from rest to now, F = A_pq A_qq^-1 (Mueller et al. 2005, eq. 7):
    // A_pq = sum (p - c) q^T as in the rotation fit, A_qq = sum q q^T of the rest offsets (fixed,
    // so its inverse is kept). Identity at rest; used only to draw the skin.
    Matrix3x3 restInverseQQ = Matrix3x3::zero(); // A_qq^-1 (zero: a flat cluster has no linear fit)
    Matrix3x3 deformation;                        // F now
};

struct SoftBody {
    int object = -1;                   // particle object id (for self-collision filtering)
    int group = -1;                    // particle group (ParticleSystem::removeGroup removes it whole)
    std::vector<int> particles;
    // The material. With ShapeMatching, material.stiffness is the fraction of the way back to the
    // rest shape per substep (0.05 jelly, 0.3 rubber, 1 rigid), spread over the solver's passes as
    // k' = 1 - (1 - k)^(1/n) (Mueller et al. 2007, PBD).
    SoftMaterial material;
    Vector3 color{0.9f, 0.4f, 0.4f};
    bool hasTets() const { return !tets.empty(); } // the NeoHookean model (else: clusters)
    // NeoHookean: the tetrahedra, grouped by the cube of the lattice they cut (cubeStart[c] ..
    // cubeStart[c + 1]) and the cubes by colour (colourStart[k] .. colourStart[k + 1]): cubes of
    // one colour share no particle and are solved in parallel (SoftTets.cpp).
    std::vector<SoftTet> tets;
    std::vector<int> cubeStart;
    int colourStart[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<SkinBinding> skin; // per surface vertex (NeoHookean)
    // ShapeMatching: the clusters and the radius they were built with [m] (the skinning blends within it).
    std::vector<SoftCluster> clusters;
    float clusterRadius = 0;
    // Surface for drawing: rest mesh and, per vertex, the clusters it follows with their weights
    // (smooth linear blend skinning: every cluster within clusterRadius, weighted by closeness).
    TriMesh surface;
    std::vector<std::vector<int>> vertexClusters;
    std::vector<std::vector<float>> vertexWeights;
    std::vector<int> vertexAnchor;             // per vertex: its nearest particle at rest (slot in `particles`)
    std::vector<Vector3> vertexAnchorOffset;   // per vertex: rest position - that particle's rest position
};

// Clusters for particles at rest positions `rest` (global indices `ids`): cluster centres on a
// lattice with the given spacing, each taking the particles within `radius` of it.
std::vector<SoftCluster> buildClusters(const std::vector<int>& ids, const std::vector<Vector3>& rest, float spacing,
                                       float radius);

// Binds the rest surface mesh to the body: every vertex to its nearest particle (its anchor) and to
// the clusters around it, with smooth weights.
void bindSurface(SoftBody& body, const TriMesh& restSurface, const std::vector<Vector3>& particleRest);

// Current positions of the surface vertices, from the particles' current `positions`. Tetrahedra:
// each vertex at its barycentric coordinates in its tetrahedron. Clusters: each vertex rides on
// its anchor particle, its rest offset turned and squashed by the blend of its clusters' linear
// deformations.
void skinSurface(const SoftBody& body, const std::vector<Vector3>& positions, std::vector<Vector3>& out);

// The outward surface normals of a body's particles now (the contacts' signed distance field,
// Macklin et al. 2014, sec. 5.1): every particle's rest normal turned by the rotations of the
// clusters it belongs to, averaged. Both arrays are indexed by the global particle index.
void turnSurfaceNormals(const SoftBody& body, const std::vector<Vector3>& restNormal, std::vector<Vector3>& normal);

// One shape-matching pass over all soft bodies: moves the predicted positions p towards the goals
// (average over the clusters of each particle). Particles with invMass 0 stay where they are.
// passesPerStep: how many such passes the substep makes (the stiffness is spread over them).
void solveShapeMatching(std::vector<SoftBody>& bodies, std::vector<Vector3>& p, const std::vector<float>& invMass, int passesPerStep);

// ---------------------------------------------------------------------------------------------
// The Neo-Hookean tetrahedra (SoftTets.cpp)
// ---------------------------------------------------------------------------------------------

// The particle lattice of a body being built: point (i, j, k) sits at origin + (i, j, k) spacing;
// node numbers the points inside the body (-1: outside). buildLatticeTets fills cubeIndex: per
// cube of the lattice (the one whose lowest corner is point (i, j, k)) its number in the body,
// -1 if it holds no tetrahedron.
struct SoftLattice {
    int nx = 0, ny = 0, nz = 0;
    Vector3 origin;
    float spacing = 0;
    std::vector<int> node;
    std::vector<int> cubeIndex;
    int point(int i, int j, int k) const { return i + nx * (j + ny * k); }
    int nodeAt(int i, int j, int k) const {
        const bool outside = i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz;
        return outside ? -1 : node[size_t(point(i, j, k))];
    }
};

// Cuts every cube of the lattice into five tetrahedra - the ones whose four corners are all inside
// the body are kept - with their corners as lattice node numbers; fills body.tets (no rest shape
// yet), body.cubeStart, body.colourStart and lattice.cubeIndex.
void buildLatticeTets(SoftLattice& lattice, SoftBody& body);

// The rest shape of every tetrahedron from the particles' rest positions (indexed like the corners).
void setRestShape(SoftBody& body, const std::vector<Vector3>& rest);

// Embeds the rest surface in the tetrahedra (body.skin): every vertex in the tetrahedron that
// contains it, or - outside them all - the nearest one, found among the cubes around it.
void bindSurfaceToTets(SoftBody& body, const SoftLattice& lattice, const TriMesh& restSurface,
                       const std::vector<Vector3>& rest);

// One XPBD solve of the tetrahedra of one cube (SoftTets.cpp): step length h; with `stepStart`
// (the positions at the start of the step) the Rayleigh damping acts too.
void solveCubeTets(SoftBody& body, int cube, float h, std::vector<Vector3>& p, const std::vector<float>& invMass,
                   const std::vector<Vector3>* stepStart);

// The body's elastic energy now [J]: the Neo-Hookean energy density integrated over the tetrahedra.
double elasticEnergy(const SoftBody& body, const std::vector<Vector3>& x);

// The outward surface normals of a tetrahedral body's particles now: a normal turns with the
// deformation as n' ~ cof(F) n = det(F) F^-T n (Nanson's formula), F of the tetrahedra around it.
void turnNormalsWithTets(const SoftBody& body, const std::vector<Vector3>& x, const std::vector<Vector3>& restNormal,
                         std::vector<Vector3>& normal);

// The highest vibration frequency of the body's lattice [rad/s]: omega = 2 c / s, with c =
// sqrt((lambda + 2 mu) / rho) the speed of sound (the pressure wave) and s the spacing. The
// small steps must be short against 1 / omega (ParticleSystem::softSmallSteps).
float highestFrequency(const SoftMaterial& material, float spacing);

} // namespace rf
