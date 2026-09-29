#pragma once
// Soft body made of particles and tetrahedra. The particles fill the body's mesh on the particle
// solver's lattice (spacing 2r), as before, so they collide with liquid, cloth, bodies and each
// other like every other particle. Every cube of eight lattice particles is cut into six
// tetrahedra around a diagonal (Kuhn 1960, Freudenthal 1942), the diagonal mirrored from cube to
// cube so the material has no preferred direction; neighbouring cubes share whole faces, and a
// cube cut at the body's edge keeps the tetrahedra whose four corners are there.
// The tetrahedra carry the shape and the particles the volume of the material - a Young's modulus
// and a Poisson ratio (SoftBodySolver.cpp).
//
// For drawing, every vertex of the body's surface mesh rides in one tetrahedron: its barycentric
// coordinates there, fixed at rest, place it among the four particles now (a vertex just outside
// the lattice extrapolates the same linear map - the skin squashes and stretches with the body).

#include "core/Mesh.h"
#include "math/Math.h"

#include <array>
#include <vector>

namespace rf {

// What a soft body is made of. Typical Young's moduli: jelly 1-50 kPa, soft foam 50-500 kPa,
// silicone and soft rubber 0.5-5 MPa, tyre rubber ~10 MPa. Poisson's ratio: 0 = a cork that does
// not bulge when squeezed, 0.3 = most foams, 0.45-0.49 = rubber and jelly (nearly incompressible).
struct SoftMaterial {
    float density = 150.0f;      // [kg/m^3]
    float youngModulus = 5e4f;   // E [Pa]
    float poissonRatio = 0.3f;   // nu, 0 .. 0.49
    // Internal friction [1/s]: the wobble - every particle's velocity less the motion of the body
    // as a whole - decays as exp(-damping t). 3: a jelly's wobble is down to 5 % in a second;
    // 0: it rings until the contacts stop it.
    float damping = 3.0f;
};

// Lamé's first parameter mu (the shear modulus) and the second, lambda, from E and nu.
inline float shearModulus(const SoftMaterial& m) { return m.youngModulus / (2.0f * (1.0f + m.poissonRatio)); }
inline float lameLambda(const SoftMaterial& m) {
    const float nu = std::min(m.poissonRatio, 0.499f);
    return m.youngModulus * nu / ((1.0f + nu) * (1.0f - 2.0f * nu));
}

struct SoftTet {
    std::array<int, 4> v{};            // global particle indices; positive volume at rest
    Matrix3x3 restInverse;             // D_m^-1, D_m = [x1 - x0, x2 - x0, x3 - x0] at rest
    float restVolume = 0;              // [m^3]
    // XPBD multiplier of the shape constraint F - R, nine components in world axes (SoftBodySolver.cpp).
    Matrix3x3 lambdaShape = Matrix3x3::zero();
    Quaternion rotation;               // the rotation of F last found (the next search starts there)
};

// The volume is held per particle, not per tetrahedron (Bonet & Burton 1998, "A simple average
// nodal pressure tetrahedral element"): particle i's volume ratio is the mean of det F over the
// tetrahedra around it, J_i = sum_t (V_t / 4) det F_t / V_i, V_i = sum_t V_t / 4, and the constraint
// J_i - 1 has compliance 1 / (lambda V_i) (why: SoftBodySolver.cpp).
struct SoftNode {
    float restVolume = 0;       // V_i [m^3]
    float lambda = 0;           // XPBD multiplier of J_i - 1
    int starBegin = 0, starEnd = 0; // its tetrahedra: SoftBody::nodeStar[starBegin .. starEnd)
    int nearBegin = 0, nearEnd = 0; // their particles, each once: SoftBody::nodeNear[nearBegin .. nearEnd), at most 27
};

// What one pass needs of a body's material and step (SoftBodySolver.cpp).
struct SoftPassConstants {
    float shapeCompliance = 0;  // 1 / (2 mu) [1/Pa]; divided by the tetrahedron's volume
    float volumeCompliance = 0; // 1 / lambda [1/Pa]; divided by the node's volume; 0: nu = 0, no volume term
    float invDt2 = 0;           // 1 / dt^2 of the step the multipliers belong to
};

// A run of one colour of one body in a pass over several bodies (solveSoftBodies' work space).
struct SoftColourRun {
    int body = 0, begin = 0, offset = 0; // the body, its first index of the colour, the run's first place in the loop
};

struct SoftBody {
    int object = -1;                   // particle object id (for self-collision filtering)
    int group = -1;                    // particle group (ParticleSystem::removeGroup removes it whole)
    std::vector<int> particles;        // contiguous: particles[k] = particles.front() + k (slot k)
    // Tetrahedra sorted by colour: no two of one colour share a particle, so a colour is solved in
    // parallel and the result is the same on any number of threads. Colour c is
    // tets[colourStart[c] .. colourStart[c + 1]).
    std::vector<SoftTet> tets;
    std::vector<int> colourStart;
    // One node per particle (slot k = particles[k]); nodeOrder lists them by colour - the
    // tetrahedra around two nodes of one colour share no particle - colour c being
    // nodeOrder[nodeColourStart[c] .. nodeColourStart[c + 1]).
    std::vector<SoftNode> nodes;
    std::vector<int> nodeStar, nodeOrder, nodeColourStart;
    // Per node: the place among its particles (nodeNear) of every corner of every tetrahedron of its
    // star - nodeSlot[4 s + c], s the place in nodeStar - found once, so an update sums its gradient
    // without searching; and that gradient as its last update found it, for the next warm start.
    std::vector<int> nodeNear;
    std::vector<uint8_t> nodeSlot;
    std::vector<Vector3> nodeGradient;
    SoftMaterial material;
    Vector3 color{0.9f, 0.4f, 0.4f};
    float multiplierStep = 0;       // the small step the multipliers belong to (0: none yet)
    bool touchesOthers = false;     // next to liquid or cloth in this step: its material is solved in the passes too
    // Work space of a pass (SoftBodySolver.cpp): its constants, the moves it found the particles at,
    // the warm start's moves.
    SoftPassConstants pass;
    std::vector<Vector3> passStart, warmMove;
    // Surface for drawing: rest mesh; per vertex the tetrahedron it rides in and its barycentric
    // weights for the corners v[1], v[2], v[3] (v[0] takes the rest).
    TriMesh surface;
    std::vector<int> vertexTet;
    std::vector<Vector3> vertexWeights;
};

// The Kuhn tetrahedra of a set of lattice points: `cells` are the integer lattice coordinates of
// the points; returned are the four point indices of every tetrahedron whose corners all exist.
std::vector<std::array<int, 4>> latticeTetrahedra(const std::vector<std::array<int, 3>>& cells);

// The body's tetrahedra from point quadruples (global particle indices) at the rest positions:
// oriented to positive volume, D_m^-1 and the volume kept, then coloured and sorted by colour.
void buildTetrahedra(SoftBody& body, const std::vector<std::array<int, 4>>& quads, const std::vector<Vector3>& rest);
void buildNodes(SoftBody& body); // the nodes of the volume constraints (called by buildTetrahedra)

// Binds the rest surface mesh to the body: every vertex to the tetrahedron around it (or, outside
// the lattice, the nearest one), by barycentric weights. `rest` is indexed by the global particle index.
void bindSurface(SoftBody& body, const TriMesh& restSurface, const std::vector<Vector3>& rest);

// Current positions of the surface vertices from the particles' current positions.
void skinSurface(const SoftBody& body, const std::vector<Vector3>& positions, std::vector<Vector3>& out);

// The deformation gradient of a tetrahedron now, F = D_s D_m^-1.
Matrix3x3 deformationGradient(const SoftTet& t, const std::vector<Vector3>& positions);
// The same at the positions base + u, every edge summed from the edge of the base and the edge of u:
// a move u far below the float step of the base still counts (SoftBodySolver.cpp).
Matrix3x3 deformationGradient(const SoftTet& t, const std::vector<Vector3>& base, const std::vector<Vector3>& u);

// One XPBD pass of the material of several bodies (SoftBodySolver.cpp): every tetrahedron's shape
// and every particle's volume, colour by colour - one colour of all the bodies in one parallel loop
// (they share no particle) - for a step of length dt, at the positions base + u; it moves u only.
// withWarmStart: the multipliers carried from the last small step are applied first (the bodies' own
// small steps); without, the pass goes on from the multipliers as they stand (the unified passes of
// the particles). `runs` is work space.
void solveSoftBodies(const std::vector<SoftBody*>& bodies, const std::vector<Vector3>& base, std::vector<Vector3>& u,
                     const std::vector<float>& invMass, float dt, bool withWarmStart, std::vector<SoftColourRun>& runs);

// Every multiplier of the body times f (the same forces over a step of another length: lambda ~ h^2).
void scaleSoftMultipliers(SoftBody& body, float f);

// The body's internal friction on the velocities v at the positions x, at the end of a step dt
// (SoftBodySolver.cpp): its wobble decays, its motion as a whole stays.
void dampSoftBody(const SoftBody& b, std::vector<Vector3>& v, const std::vector<Vector3>& x, const std::vector<float>& invMass, float dt);

// The outward surface normals of a body's particles now (the contacts' signed distance field,
// Macklin et al. 2014, sec. 5.1): a normal is carried by the cofactor of the deformation, n ~ cof(F) n0
// (Nanson's formula), with F summed over the tetrahedra around the particle, weighted by volume.
// Both normal arrays are indexed by the global particle index.
void turnSurfaceNormals(const SoftBody& body, const std::vector<Vector3>& positions, const std::vector<Vector3>& restNormal,
                        std::vector<Vector3>& normal);

} // namespace rf
