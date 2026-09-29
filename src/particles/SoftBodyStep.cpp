// The soft bodies' own step inside every substep of the particles (ParticleSystem's): small steps,
// each the whole recipe - flight, material, walls, rigid bodies, other soft bodies - as Jolt's soft
// bodies and PhysX's TGS substeps do it ("Small Steps in Physics Simulation", Macklin et al. 2019).
// The contacts are found once per step, a plane per particle against a rigid body and a pair per
// two particles of different soft bodies, and solved in every small step: positions pushed out
// without a speed, then the velocities by sequential impulses. The material is in SoftBodySolver.cpp.
#include "particles/ParticleSystem.h"

#include "core/Parallel.h"

#include <algorithm>
#include <array>

namespace rf {

// How many small steps a soft body needs in a step dt: enough that a shear wave, c = sqrt(mu / rho),
// crosses at most params.softCourant lattice cells in one (a Courant number). The static answer
// does not depend on it (the multipliers are carried, SoftBodySolver.cpp); the dynamics do - a
// jelly cube (40 kPa, 150 kg/m^3, 3 cm cells) at 180 Hz takes 2, soft rubber (5 MPa, 1000 kg/m^3,
// 1.7 cm) 15.
int ParticleSystem::softSmallSteps(const SoftBody& b, float dt) const {
    const float wave = std::sqrt(shearModulus(b.material) / std::max(b.material.density, 1e-3f));
    const int needed = int(std::ceil(wave * dt / (std::max(params.softCourant, 1e-3f) * spacing())));
    return std::clamp(needed, std::max(1, params.softSubsteps), std::max(1, params.maxSoftSubsteps));
}

// The contacts a soft particle may meet in this step, found once from where the step starts (as
// Jolt's soft bodies and "Small Steps", sec. 4.2, do), within the reach of how far it can move in
// the step.
void ParticleSystem::findSoftContacts(float dt) {
    softBodyContacts_.clear();
    softPairs_.clear();
    findRigidPlanes(dt);
    findSoftPairs(dt);
}

// How far a particle can go in a step: its speed and gravity's half g dt^2.
float ParticleSystem::softTravel(int i, float dt) const {
    return length(v_[size_t(i)]) * dt + 0.5f * length(params.gravity) * dt * dt;
}

// Against a rigid body: the plane of the body's surface nearest to the particle, the body where it
// ends this step (the rigid solver steps first). The plane stays for the whole step, so a particle
// pressed into a box goes back out the way it came in - never out through the far face, the way
// the edge of a box used to cut into a jelly.
void ParticleSystem::findRigidPlanes(float dt) {
    if (!rigid_) return;
    prepareBodyQuery(false);
    std::vector<int>& candidates = bodyCandidates_[size_t(ThreadPool::workerIndex())];
    const auto& bodies = rigid_->bodies();
    const float r = params.particleRadius;
    for (const SoftBody& b : softBodies_)
        for (int i : b.particles) {
            if (invMass_[size_t(i)] == 0) continue;
            const Vector3& x = x_[size_t(i)];
            const float reach = 1.5f * r + softTravel(i, dt);
            rigid_->queryBodies(AABB(x - Vector3(reach), x + Vector3(reach)), candidates);
            for (int k : candidates) {
                Vector3 n;
                const float d = bodies[size_t(k)].signedDistance(x, n);
                if (d < reach) softBodyContacts_.push_back({i, k, n, x - n * d});
            }
        }
}

// Between two soft bodies: every pair of their particles that can come within d0 in the step (a
// sorted grid of the soft particles, cells as big as the farthest two can close in a step), pushed
// apart along the line of their centres at the step's start. Two bodies found inside each other -
// an overlap deeper than a quarter spacing, or one the pre-stabilization is still pulling apart -
// are only kept from going deeper: they come apart as wholes, without a speed
// (ParticleContacts.cpp); pushed apart particle by particle they would dent, and the dents spring
// back into a jump.
void ParticleSystem::findSoftPairs(float dt) {
    if (softBodies_.size() < 2) return;
    const float r = params.particleRadius, d0 = spacing();
    float fastest = 0;
    for (const SoftBody& b : softBodies_)
        for (int i : b.particles) fastest = std::max(fastest, softTravel(i, dt));
    const float cell = d0 + 2.0f * fastest + 0.25f * r;
    auto cellOf = [&](const Vector3& p) { return std::array<int, 3>{int(std::floor(p.x / cell)), int(std::floor(p.y / cell)), int(std::floor(p.z / cell))}; };
    auto key = [](int x, int y, int z) { return ((long long)(x & 0x1fffff) << 42) | ((long long)(y & 0x1fffff) << 21) | (long long)(z & 0x1fffff); };
    softCells_.clear();
    for (const SoftBody& b : softBodies_)
        for (int i : b.particles) {
            const auto c = cellOf(x_[size_t(i)]);
            softCells_.push_back({key(c[0], c[1], c[2]), i});
        }
    std::sort(softCells_.begin(), softCells_.end());
    auto pairWith = [&](int i, int j) {
        if (j <= i || object_[size_t(j)] == object_[size_t(i)] || invMass_[size_t(i)] + invMass_[size_t(j)] == 0) return;
        const Vector3 centres = x_[size_t(i)] - x_[size_t(j)];
        const float apart = length(centres), reach = d0 + softTravel(i, dt) + softTravel(j, dt) + 0.25f * r;
        if (apart >= reach || apart < 1e-9f) return;
        // softBodyOf_ is the last pre-stabilization's: particles appended since (the emitter's) are not in it.
        bool inside = apart < 0.75f * d0;
        if (!inside && !unresolved_.empty() && size_t(j) < softBodyOf_.size())
            inside = std::binary_search(unresolved_.begin(), unresolved_.end(), moverPair({i, j, Vector3(0.0f), 1.0f, d0, false}));
        softPairs_.push_back({i, j, centres / apart, inside ? std::min(apart, d0) : d0});
    };
    for (const SoftBody& b : softBodies_)
        for (int i : b.particles) {
            const auto c = cellOf(x_[size_t(i)]);
            for (int dz = -1; dz <= 1; ++dz)
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dx = -1; dx <= 1; ++dx) {
                        const long long k = key(c[0] + dx, c[1] + dy, c[2] + dz);
                        for (auto it = std::lower_bound(softCells_.begin(), softCells_.end(), std::make_pair(k, -1));
                             it != softCells_.end() && it->first == k; ++it)
                            pairWith(i, it->second);
                    }
        }
}

namespace {

// A contact's answer in the velocities, after its positions were put right (Jolt's order: the push
// out of the other side does not change the velocity, this does), one visit per small step as
// Jolt's soft bodies do: the approach along n stops - inelastic, a soft body's bounce is its
// material's - and Coulomb's friction takes back the slide, at most mu times the normal impulse (in
// a small step a body at rest is pressed in by its weight, g h, and mu >= tan(slope) holds it on a
// slope). `relative` is the particle's velocity against the other side, `inverseMass(d)` the pair's
// inverse mass along the direction d, `apply(J)` gives the particle the impulse J and the other side
// -J. (Four sweeps of accumulated impulses, as Box2D's, were tried: they did not settle a stack of
// soft barrels better, only differently. What had thrown a box off a jelly it landed on was the
// jelly's lattice lying off-centre in its box, ParticleSystem::addSoftBody.)
template <class InverseMass, class Apply>
void contactImpulse(const Vector3& relative, const Vector3& n, float mu, InverseMass inverseMass, Apply apply) {
    const float vn = dot(relative, n);
    if (vn >= 0) return;
    const float pushed = -vn / inverseMass(n);
    const Vector3 slide = relative - n * vn;
    const float along = length(slide);
    Vector3 friction(0.0f);
    if (along > 1e-9f) {
        friction = slide * (-1.0f / inverseMass(slide / along));
        const float limit = mu * pushed, size = length(friction);
        if (size > limit) friction *= limit / size;
    }
    apply(n * pushed + friction);
}

// Sweeps of the pairs' positions per small step, alternating their direction. A blow between two
// bodies is dozens of pairs pushing on each other; one sweep left them 0.28 d0 deep, past the
// quarter spacing that marks two bodies inside each other, the pairs were then only held, and
// falling soft barrels sank 0.47 d0 into each other (four sweeps: 0.02 d0).
constexpr int kPairSweeps = 4;

} // namespace

// The contacts of the soft particles in one small step: every particle pushed out of the planes of
// the rigid bodies and apart from the other soft bodies' particles (positions, no speed; the pairs
// in kPairSweeps sweeps), then the velocities (contactImpulse). A rigid body keeps its pose; its velocity changes in the step's copy
// (bodyV_, bodyW_), and what the step gave it is added to the body at the end (bodyDv_, bodyDw_) -
// Jolt's way. The lever arm runs to the contact point on the plane, r below the particle's centre.
void ParticleSystem::solveSoftContacts() {
    const float r = params.particleRadius, mu = params.solidFriction;
    std::vector<Vector3>& u = softMove_; // the particles are at x_ + u (stepSoftBodies)
    for (SoftBodyContact& c : softBodyContacts_) {
        const size_t i = size_t(c.particle);
        const float depth = r - dot((x_[i] - c.point) + u[i], c.normal);
        c.touching = depth > 0;
        if (c.touching) u[i] += c.normal * depth;
    }
    const int np = int(softPairs_.size());
    for (SoftPair& c : softPairs_) c.touching = false;
    for (int sweep = 0; sweep < kPairSweeps; ++sweep)
        for (int q = 0; q < np; ++q) {
            SoftPair& c = softPairs_[size_t(sweep % 2 ? np - 1 - q : q)];
            const size_t i = size_t(c.i), j = size_t(c.j);
            const float depth = pushAlong((x_[i] - x_[j]) + (u[i] - u[j]), c.normal, c.target), wi = invMass_[i], wj = invMass_[j];
            if (depth <= 0) continue;
            c.touching = true;
            u[i] += c.normal * (depth * wi / (wi + wj));
            u[j] -= c.normal * (depth * wj / (wi + wj));
        }
    if (rigid_)
        for (const SoftBodyContact& c : softBodyContacts_) {
            if (!c.touching) continue;
            const size_t i = size_t(c.particle), k = size_t(c.body);
            const RigidBody& body = rigid_->bodies()[k];
            const Vector3 arm = (x_[i] - body.pos) + u[i] - c.normal * r;
            const float wp = invMass_[i];
            contactImpulse(v_[i] - (bodyV_[k] + cross(bodyW_[k], arm)), c.normal, mu,
                           [&](const Vector3& d) { const Vector3 rd = cross(arm, d); return wp + body.invMass + dot(rd, body.applyInvInertiaWorld(rd)); },
                           [&](const Vector3& J) {
                               v_[i] += J * wp;
                               const Vector3 dV = J * -body.invMass, dW = body.applyInvInertiaWorld(cross(arm, J)) * -1.0f;
                               bodyV_[k] += dV, bodyW_[k] += dW, bodyDv_[k] += dV, bodyDw_[k] += dW;
                           });
        }
    for (const SoftPair& c : softPairs_) {
        if (!c.touching) continue;
        const size_t i = size_t(c.i), j = size_t(c.j);
        const float wi = invMass_[i], wj = invMass_[j];
        contactImpulse(v_[i] - v_[j], c.normal, mu, [&](const Vector3&) { return wi + wj; },
                       [&](const Vector3& J) { v_[i] += J * wi, v_[j] -= J * wj; });
    }
}

// The small steps of all soft bodies, m of h = dt / m, each:
//   1. fly: v += g h, x += v h (a held particle - pinned, grabbed - slides to where the step puts it);
//   2. the material: one XPBD pass, warm-started (SoftBodySolver.cpp); v += its move / h;
//   3. the walls and the obstacle, then the rigid bodies and the other soft bodies (solveSoftContacts),
//      each contact pushed out and its approach stopped with friction (contactImpulse).
// The body leaves the step with its last small step's velocity: it falls at exactly g, rests with
// none, and sags as its material says. (Contacts solved after the small steps instead, in the
// passes, got no velocity right: the mean over the step lagged gravity - a body fell at 0.63 g -,
// the last small step's left a supported body a phantom 20 mm/s.) The liquid and the cloth meet the
// soft particles in the passes that follow; what they push a soft particle by is added to its
// position and, over dt, to its velocity (finishStep).
void ParticleSystem::stepSoftBodies(float dt) {
    if (softBodies_.empty()) return;
    int m = 1;
    for (const SoftBody& b : softBodies_) m = std::max(m, softSmallSteps(b, dt));
    lastSoftSmallSteps_ = m;
    const float h = dt / float(m), mu = params.solidFriction;
    const Vector3 g = params.gravity;
    findSoftContacts(dt);
    if (rigid_) {
        const auto& bodies = rigid_->bodies();
        bodyV_.resize(bodies.size());
        bodyW_.resize(bodies.size());
        for (size_t k = 0; k < bodies.size(); ++k) bodyV_[k] = bodies[k].vel, bodyW_[k] = bodies[k].angVel;
    }
    std::vector<Vector3>& u = softMove_; // every soft particle is at x_ + u (SoftBodySolver.cpp: why)
    for (SoftBody& b : softBodies_) {
        for (int i : b.particles) {
            u[size_t(i)] = Vector3(0.0f);
            if (invMass_[size_t(i)] == 0) softFlight_[size_t(i)] = p_[size_t(i)] - x_[size_t(i)]; // where the held ones end
        }
        // The multipliers are force x h^2: carried from the last small step, rescaled if h changed.
        const float scale = b.multiplierStep > 0 ? sqr(h / b.multiplierStep) : 0.0f;
        if (scale != 1.0f) scaleSoftMultipliers(b, scale);
        b.multiplierStep = h;
    }
    for (int s = 0; s < m; ++s) {
        for (SoftBody& b : softBodies_) {
            for (int i : b.particles) { // 1. fly
                const size_t k = size_t(i);
                if (invMass_[k] == 0) {
                    u[k] = softFlight_[k] * (float(s + 1) / float(m));
                    continue;
                }
                v_[k] += g * h;
                u[k] += v_[k] * h;
                softFlight_[k] = u[k];
            }
            solveSoftBody(b, x_, u, invMass_, h, true); // 2. the material
            for (int i : b.particles) { // 3. the walls and the obstacle
                const size_t k = size_t(i);
                if (invMass_[k] == 0) continue;
                v_[k] += (u[k] - softFlight_[k]) / h;
                Vector3 at = x_[k] + u[k];
                const Vector3 before = at;
                collideWallsAndMesh(at, at); // pushed out without friction; the friction is in the velocities
                const Vector3 push = at - before;
                const float depth = length(push), w = invMass_[k];
                if (depth == 0) continue;
                u[k] += push;
                contactImpulse(v_[k], push / depth, mu, [&](const Vector3&) { return w; }, [&](const Vector3& J) { v_[k] += J * w; });
            }
        }
        solveSoftContacts(); // 3. the rigid bodies and the other soft bodies
    }
    for (const SoftBody& b : softBodies_)
        for (int i : b.particles) p_[size_t(i)] = x_[size_t(i)] + u[size_t(i)];
}

} // namespace rf
