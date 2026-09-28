// Soft bodies in ParticleSystem: building one from a closed mesh (the particle lattice inside it,
// cut into tetrahedra or grouped into clusters), and the tetrahedra's small steps. The material
// model itself - the stable Neo-Hookean constraints - is in SoftTets.cpp, the legacy clusters in
// SoftBody.cpp, the contacts with everything else in ParticleContacts.cpp.
#include "particles/ParticleSystem.h"

#include "core/Parallel.h"

#include <algorithm>

namespace rf {

// ---------------------------------------------------------------------------------------------
// Building
// ---------------------------------------------------------------------------------------------

int ParticleSystem::addSoftBody(const TriMesh& shape, float density, float stiffness, const Vector3& color,
                                const Vector3& velocity) {
    SoftMaterial material;
    material.density = density;
    material.youngModulus = youngFromStiffness(stiffness);
    material.stiffness = stiffness;
    return addSoftBody(shape, material, color, velocity);
}

int ParticleSystem::addSoftBody(const TriMesh& shape, const SoftMaterial& material, const Vector3& color,
                                const Vector3& velocity) {
    MeshBVH bvh;
    bvh.build(shape);
    SoftLattice lattice = fillLattice(bvh, shape.bounds());
    SoftBody body;
    body.object = nextObject_++;
    body.group = nextGroup_++;
    body.material = material;
    body.color = color;
    // The skin: the mesh cut finer (about a spacing and a half per edge), so it can bend with the body.
    const TriMesh skin = primitives::subdivided(shape, 1.5f * spacing());
    const bool tetrahedra = material.model == SoftModel::NeoHookean && addTetBody(body, lattice, skin, velocity);
    if (!tetrahedra) {
        body.material.model = SoftModel::ShapeMatching;
        addClusterBody(body, lattice, skin, velocity);
    }
    if (body.particles.empty()) return -1;
    measureSoftSurface(body, bvh);
    softBodies_.push_back(std::move(body));
    tetBatchesStale_ = true;
    return int(softBodies_.size()) - 1;
}

// The lattice with spacing s = 2r over the mesh's box, its first point half a spacing in from the
// corner: the particles' spheres then reach exactly to a face of a box-shaped mesh. The points
// inside the mesh (and the walls) are numbered in order, x fastest.
SoftLattice ParticleSystem::fillLattice(const MeshBVH& bvh, const AABB& bounds) const {
    const float s = spacing();
    SoftLattice lattice;
    lattice.spacing = s;
    lattice.origin = bounds.lo + Vector3(0.5f * s);
    const Vector3 reach = bounds.hi - lattice.origin;
    lattice.nx = std::max(0, int(std::ceil(reach.x / s)));
    lattice.ny = std::max(0, int(std::ceil(reach.y / s)));
    lattice.nz = std::max(0, int(std::ceil(reach.z / s)));
    lattice.node.assign(size_t(lattice.nx) * size_t(lattice.ny) * size_t(lattice.nz), -1);
    int count = 0;
    for (int k = 0; k < lattice.nz; ++k)
        for (int j = 0; j < lattice.ny; ++j)
            for (int i = 0; i < lattice.nx; ++i) {
                const Vector3 p = lattice.origin + Vector3(float(i), float(j), float(k)) * s;
                if (walls_.contains(p) && bvh.isInside(p)) lattice.node[size_t(lattice.point(i, j, k))] = count++;
            }
    return lattice;
}

// Particles for the lattice nodes marked in `keep`, in the lattice's order: mass = density s^3
// each (every particle stands for its cube of the lattice, so the body weighs what its mesh
// does), the material's friction. globalOfNode: the particle index of every node (-1: none).
void ParticleSystem::addSoftParticles(SoftBody& body, const SoftLattice& lattice, const std::vector<bool>& keep,
                                      const Vector3& velocity, std::vector<int>& globalOfNode) {
    const float s = spacing();
    const float invMass = 1.0f / (std::max(body.material.density, 1e-3f) * s * s * s);
    globalOfNode.assign(keep.size(), -1);
    for (int k = 0; k < lattice.nz; ++k)
        for (int j = 0; j < lattice.ny; ++j)
            for (int i = 0; i < lattice.nx; ++i) {
                const int node = lattice.nodeAt(i, j, k);
                if (node < 0 || !keep[size_t(node)] || int(x_.size()) >= params.maxParticles) continue;
                globalOfNode[size_t(node)] = int(x_.size());
                body.particles.push_back(int(x_.size()));
                addParticle(lattice.origin + Vector3(float(i), float(j), float(k)) * s, velocity, ParticlePhase::Soft, body.object,
                            body.group, invMass);
                friction_.back() = std::max(0.0f, body.material.friction);
            }
}

// The Neo-Hookean body: the lattice's cubes cut into tetrahedra (SoftTets.cpp); only the nodes of
// some tetrahedron become particles - a lone node (the tip of a thin spike of the mesh) would have
// nothing holding it. False (nothing added) when no tetrahedron fits or the particles would not.
bool ParticleSystem::addTetBody(SoftBody& body, SoftLattice& lattice, const TriMesh& skin, const Vector3& velocity) {
    buildLatticeTets(lattice, body);
    const int nodes = int(std::count_if(lattice.node.begin(), lattice.node.end(), [](int n) { return n >= 0; }));
    std::vector<bool> used(size_t(nodes), false);
    for (const SoftTet& t : body.tets)
        for (int v : t.v) used[size_t(v)] = true;
    const size_t needed = size_t(std::count(used.begin(), used.end(), true));
    if (body.tets.empty() || x_.size() + needed > size_t(params.maxParticles)) {
        body.tets.clear();
        body.cubeStart.clear();
        return false;
    }
    std::vector<int> global;
    addSoftParticles(body, lattice, used, velocity, global);
    for (SoftTet& t : body.tets)
        for (int& v : t.v) v = global[size_t(v)];
    setRestShape(body, rest_);
    bindSurfaceToTets(body, lattice, skin, rest_);
    return true;
}

// The legacy body of clusters on every lattice node. Clusters every 1.5 particle spacings, each 2
// spacings in radius (as FleX): a cluster spans ~4 particles, so a body a few particles thick bends
// and squashes between its clusters. Bigger clusters (3 / 4 spacings) covered a small body whole and
// made it rigid. A body must be at least 3 particles across: a cluster of a flat sheet of
// particles has no definite rotation, and its skin flies apart.
void ParticleSystem::addClusterBody(SoftBody& body, const SoftLattice& lattice, const TriMesh& skin, const Vector3& velocity) {
    const float s = spacing();
    const int nodes = int(std::count_if(lattice.node.begin(), lattice.node.end(), [](int n) { return n >= 0; }));
    std::vector<int> global;
    addSoftParticles(body, lattice, std::vector<bool>(size_t(nodes), true), velocity, global);
    if (body.particles.empty()) return;
    std::vector<Vector3> rest;
    for (int i : body.particles) rest.push_back(rest_[size_t(i)]);
    body.clusterRadius = 2.0f * s;
    body.clusters = buildClusters(body.particles, rest, 1.5f * s, body.clusterRadius);
    bindSurface(body, skin, rest);
}

// A soft body's signed distance field at a point p inside its mesh (Macklin et al. 2014, sec.
// 5.1): the depth under the surface and the outward gradient. The gradient is the central
// difference of the signed distance over a particle radius each way, not the normal of the nearest
// triangle. The two differ along an edge and at a corner, where the nearest face is any one of two
// or three: the difference quotient leans out between them - the diagonal arrows at the corners
// of the paper's Fig. 7. A particle on the rim of a face must count as "up and out", or a body
// resting on that face is pushed sideways by the one-sided contacts at the rim (eq. 20). Where
// the differences cancel - the middle of the body, its medial axis - the nearest face decides.
void ParticleSystem::measureSurface(const MeshBVH& bvh, const Vector3& p, float& depth, Vector3& normal) const {
    ClosestHit nearest;
    bvh.closestPoint(p, kInf, nearest);
    depth = std::max(0.0f, -nearest.signedDistance);
    const float h = params.particleRadius;
    auto slope = [&](const Vector3& e) { return bvh.signedDistance(p + e * h, kInf) - bvh.signedDistance(p - e * h, kInf); };
    const Vector3 gradient(slope(Vector3(1, 0, 0)), slope(Vector3(0, 1, 0)), slope(Vector3(0, 0, 1)));
    normal = length2(gradient) > 1e-6f * h * h ? normalize(gradient) : normalize(nearest.normal);
}

// The body's signed distance field, sampled at its particles (Macklin et al. 2014, sec. 5.1,
// Fig. 7): how deep under the mesh's surface each one sits and which way is out. The contacts
// with other bodies take their normal from it, so a deep overlap comes apart the way out of the
// body, not along whichever neighbour happens to be nearest.
void ParticleSystem::measureSoftSurface(const SoftBody& body, const MeshBVH& bvh) {
    for (int i : body.particles) {
        measureSurface(bvh, rest_[size_t(i)], surfaceDepth_[size_t(i)], restSurfaceNormal_[size_t(i)]);
        surfaceNormal_[size_t(i)] = restSurfaceNormal_[size_t(i)];
    }
}

double ParticleSystem::softElasticEnergy() const {
    double energy = 0;
    for (const SoftBody& b : softBodies_)
        if (b.hasTets()) energy += elasticEnergy(b, x_);
    return energy;
}

// ---------------------------------------------------------------------------------------------
// The tetrahedra's small steps
// ---------------------------------------------------------------------------------------------

void ParticleSystem::refreshTetBatches() {
    if (!tetBatchesStale_ && tetParticle_.size() == x_.size()) return; // (liquid added: the flags grow)
    tetBatchesStale_ = false;
    tetParticles_.clear();
    tetParticle_.assign(x_.size(), 0);
    for (std::vector<TetBatch>& batches : tetBatches_) batches.clear();
    for (int b = 0; b < int(softBodies_.size()); ++b) {
        const SoftBody& body = softBodies_[size_t(b)];
        if (!body.hasTets()) continue;
        tetParticles_.insert(tetParticles_.end(), body.particles.begin(), body.particles.end());
        for (int i : body.particles) tetParticle_[size_t(i)] = 1;
        for (int colour = 0; colour < 8; ++colour)
            for (int c = body.colourStart[colour]; c < body.colourStart[colour + 1]; ++c) tetBatches_[colour].push_back({b, c});
    }
}

// How many small steps a substep dt is cut into: enough that h omega <= softStepFraction for the
// stiffest body (omega its highest frequency). XPBD with one pass per step converges to the
// material's true stiffness only while h omega is small: a constraint then acts with the
// compliance alpha + h^2 w |grad C|^2 in place of alpha (Macklin et al. 2016, sec. 4; the ratio
// of the two is (h omega)^2) - a step too long makes a stiff body softer than its E.
int ParticleSystem::softSmallSteps(float dt) const {
    float omega = 0;
    for (const SoftBody& b : softBodies_)
        if (b.hasTets()) omega = std::max(omega, highestFrequency(b.material, spacing()));
    const int needed = int(std::ceil(dt * omega / std::max(params.softStepFraction, 1e-3f)));
    const int least = std::max(1, params.softSubsteps);
    return std::clamp(needed, least, std::max(least, params.softMaxSubsteps));
}

// One pass over every tetrahedron of every body, colour by colour; the cubes of a colour in
// parallel (they share no particle). The Rayleigh damping acts against the small step's start.
void ParticleSystem::solveTetPass(float h) {
    refreshTetBatches();
    for (const std::vector<TetBatch>& batches : tetBatches_)
        parallelFor(int(batches.size()), [&](int k) {
            const TetBatch& b = batches[size_t(k)];
            solveCubeTets(softBodies_[size_t(b.body)], b.cube, h, p_, invMass_, &softStart_);
        }, 8);
}

// The tetrahedral bodies' motion in m small steps of h = dt / m, each a recipe of four (as the
// cloth's, stepClothsInSmallSteps):
//   1. move: a free particle flies on under gravity, a held one slides in a straight line towards
//      where this substep puts it;
//   2. one XPBD pass over the tetrahedra, the multipliers from zero (a fresh step), then the
//      substep's contacts between elastic bodies and with cloth (ParticleContact::smallSteps);
//   3. the walls, the obstacle and the fixed bodies, with Coulomb friction (a moving body is met
//      in the main passes, where it is pushed back);
//   4. the velocity of this small step from its motion.
// What the main passes do to these particles afterwards (the liquid's push, the moving bodies,
// contacts with anything else) becomes a change of their velocity at the end of the substep
// (finishStep), and the tetrahedra answer it in the next substep's small steps. (Solving the
// tetrahedra in the main passes too, with the whole substep, only cost time: once the contacts
// between elastic bodies were in the small steps, every test gave the same numbers without it.)
void ParticleSystem::stepSoftTetsInSmallSteps(float dt) {
    refreshTetBatches();
    softSmallSteps_ = 0;
    if (tetParticles_.empty()) return;
    const int m = softSmallSteps_ = softSmallSteps(dt);
    const float h = dt / float(m);
    const Vector3 g = params.gravity;
    const int count = int(tetParticles_.size());
    softStart_.resize(x_.size());
    softVelocity_.resize(x_.size());
    softEnd_.resize(x_.size());
    parallelFor(count, [&](int k) { // the substep's start, the velocity before its gravity, a held particle's end
        const int i = tetParticles_[size_t(k)];
        softStart_[i] = x_[i];
        softVelocity_[i] = invMass_[i] > 0 ? v_[i] - g * dt : Vector3(0.0f);
        softEnd_[i] = p_[i];
    });
    prepareBodyQuery(false);
    for (int s = 0; s < m; ++s) {
        parallelFor(count, [&](int k) { // 1. move
            const int i = tetParticles_[size_t(k)];
            if (invMass_[i] == 0) {
                p_[i] = softStart_[i] + (softEnd_[i] - softStart_[i]) * (1.0f / float(m - s));
                return;
            }
            softVelocity_[i] += g * h;
            p_[i] = softStart_[i] + softVelocity_[i] * h;
        });
        for (SoftBody& b : softBodies_)
            for (SoftTet& t : b.tets) t.lambdaD = t.lambdaH = 0;
        solveTetPass(h);             // 2. the tetrahedra,
        solveParticleContacts(true); //    the contacts between elastic bodies (and cloth)
        parallelFor(count, [&](int k) { // 3. collide, 4. the velocity
            const int i = tetParticles_[size_t(k)];
            if (invMass_[i] > 0) {
                collide(i, p_[i], softStart_[i], false, h, true); // (the moving bodies: the main passes')
                softVelocity_[i] = (p_[i] - softStart_[i]) / h;
            }
            softStart_[i] = p_[i];
        });
    }
}

} // namespace rf
