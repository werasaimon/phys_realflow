// The soft bodies' own step inside every substep of the particles (ParticleSystem's): small steps,
// each the whole recipe - flight, material, walls, rigid bodies, other soft bodies - as Jolt's soft
// bodies and PhysX's TGS substeps do it ("Small Steps in Physics Simulation", Macklin et al. 2019).
// The contacts are found once per step, a plane per particle against a rigid body and a pair per
// two particles of different soft bodies, and solved in every small step: positions pushed out
// without a speed, then the velocities by sequential impulses. The material is in SoftBodySolver.cpp.
#include "particles/ParticleSystem.h"

#include "core/Parallel.h"
#include "core/Probe.h"

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
// the step. Two levels, as the broad phases of Bullet and PhysX: the soft bodies' boxes first -
// only the particles of a body whose box meets a rigid body or another soft body go on - then the
// particles, in parallel. Every search counts first and fills second, each particle at its place in
// the list: the list comes out the same on any number of threads.
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

// The box of every soft body where its particles start the step, grown by how far the farthest of
// them can go and by the reach of a contact (1.5 r).
void ParticleSystem::softBodyBoxes(float dt) {
    softBoxes_.resize(softBodies_.size());
    const float r = params.particleRadius;
    parallelFor(int(softBodies_.size()), [&](int b) {
        AABB box;
        float travel = 0;
        for (int i : softBodies_[size_t(b)].particles) box.expand(x_[size_t(i)]), travel = std::max(travel, softTravel(i, dt));
        const Vector3 grow(1.5f * r + travel);
        softBoxes_[size_t(b)] = AABB(box.lo - grow, box.hi + grow);
    }, 1);
}

// The first pass of a search over n items: count(q) says how many results item q has; offsets[q]
// becomes the place its results start at (the second pass writes them there), results in item order.
template <class Count>
static int countResults(int n, std::vector<int>& offsets, Count count) {
    offsets.assign(size_t(n) + 1, 0);
    parallelFor(n, [&](int q) { offsets[size_t(q) + 1] = count(q); }, 64);
    for (int q = 0; q < n; ++q) offsets[size_t(q) + 1] += offsets[size_t(q)];
    return offsets[size_t(n)];
}

// Against a rigid body: the plane of the body's surface nearest to the particle, the body where it
// ends this step (the rigid solver steps first). The plane stays for the whole step, so a particle
// pressed into a box goes back out the way it came in - never out through the far face, the way
// the edge of a box used to cut into a jelly. The world tree is asked once per soft body, with its box.
void ParticleSystem::findRigidPlanes(float dt) {
    if (!rigid_) return;
    prepareBodyQuery(false);
    const auto& bodies = rigid_->bodies();
    const float r = params.particleRadius;
    softNearRigid_.resize(softBodies_.size());
    softSearch_.clear(), softSearchBody_.clear();
    for (size_t b = 0; b < softBodies_.size(); ++b) {
        softNearRigid_[b].clear();
        rigid_->queryBodies(softBoxes_[b], softNearRigid_[b]); // a sleeping body too: it still holds what lies on it
        if (softNearRigid_[b].empty()) continue;
        for (int i : softBodies_[b].particles)
            if (invMass_[size_t(i)] != 0) softSearch_.push_back(i), softSearchBody_.push_back(int(b));
    }
    auto planes = [&](int q, auto emit) {
        const int i = softSearch_[size_t(q)];
        const Vector3& x = x_[size_t(i)];
        const float reach = 1.5f * r + softTravel(i, dt);
        int found = 0;
        for (int k : softNearRigid_[size_t(softSearchBody_[size_t(q)])]) {
            const RigidBody& body = bodies[size_t(k)];
            if (length2(x - body.pos) > sqr(body.boundingRadius() + reach)) continue;
            Vector3 n;
            const float d = body.signedDistance(x, n);
            if (d < reach) emit(found++, i, k, n, x - n * d);
        }
        return found;
    };
    const int n = int(softSearch_.size());
    softBodyContacts_.resize(size_t(countResults(n, softCount_, [&](int q) { return planes(q, [](int, int, int, const Vector3&, const Vector3&) {}); })));
    parallelFor(n, [&](int q) {
        planes(q, [&](int k, int i, int body, const Vector3& nrm, const Vector3& point) {
            softBodyContacts_[size_t(softCount_[size_t(q)] + k)] = {i, body, nrm, point};
        });
    }, 64);
}

// Which soft bodies' boxes meet another's: sweep and prune along x (Bullet's AxisSweep3 on one
// axis) - sorted by the low end of their boxes, a body meets only those that start before its box ends.
void ParticleSystem::softBodiesThatMeet() {
    const int B = int(softBodies_.size());
    softOrder_.resize(size_t(B));
    for (int b = 0; b < B; ++b) softOrder_[size_t(b)] = b;
    auto box = [&](int k) -> const AABB& { return softBoxes_[size_t(softOrder_[size_t(k)])]; };
    std::sort(softOrder_.begin(), softOrder_.end(), [&](int a, int b) { return softBoxes_[size_t(a)].lo.x < softBoxes_[size_t(b)].lo.x; });
    softMeets_.assign(size_t(B), 0);
    for (int a = 0; a < B; ++a)
        for (int c = a + 1; c < B && box(c).lo.x <= box(a).hi.x; ++c)
            if (box(a).overlaps(box(c))) softMeets_[size_t(softOrder_[size_t(a)])] = softMeets_[size_t(softOrder_[size_t(c)])] = 1;
}

// Between two soft bodies: every pair of their particles that can come within d0 in the step,
// pushed apart along the line of their centres at the step's start. The bodies first
// (softBodiesThatMeet), then the particles of the bodies that meet another, hashed into a grid of cells as big as the farthest two can close in a step (Teschner et
// al. 2003, "Optimized Spatial Hashing for Collision Detection of Deformable Objects"), sorted into
// the cells by counting. Two bodies found inside each other - an overlap deeper than a quarter
// spacing, or one the pre-stabilization is still pulling apart - are only kept from going deeper:
// they come apart as wholes, without a speed (ParticleContacts.cpp); pushed apart particle by
// particle they would dent, and the dents spring back into a jump.
void ParticleSystem::findSoftPairs(float dt) {
    const int B = int(softBodies_.size());
    if (B < 2) return;
    softBodiesThatMeet();
    softSearch_.clear();
    float fastest = 0;
    for (int b = 0; b < B; ++b)
        if (softMeets_[size_t(b)])
            for (int i : softBodies_[size_t(b)].particles) softSearch_.push_back(i), fastest = std::max(fastest, softTravel(i, dt));
    const int n = int(softSearch_.size());
    if (n == 0) return;
    const float r = params.particleRadius, d0 = spacing(), cell = d0 + 2.0f * fastest + 0.25f * r;
    size_t buckets = 1;
    while (buckets < size_t(2 * n)) buckets *= 2;
    auto cellOf = [&](const Vector3& p) { return std::array<int, 3>{int(std::floor(p.x / cell)), int(std::floor(p.y / cell)), int(std::floor(p.z / cell))}; };
    auto bucket = [&](int x, int y, int z) { return int((unsigned(x) * 73856093u ^ unsigned(y) * 19349663u ^ unsigned(z) * 83492791u) & unsigned(buckets - 1)); };
    softCellStart_.assign(buckets + 1, 0);
    for (int i : softSearch_) {
        const auto c = cellOf(x_[size_t(i)]);
        ++softCellStart_[size_t(bucket(c[0], c[1], c[2])) + 1];
    }
    for (size_t k = 0; k < buckets; ++k) softCellStart_[k + 1] += softCellStart_[k];
    softCellFill_.assign(softCellStart_.begin(), softCellStart_.end() - 1);
    softCellItems_.resize(size_t(n));
    for (int i : softSearch_) {
        const auto c = cellOf(x_[size_t(i)]);
        softCellItems_[size_t(softCellFill_[size_t(bucket(c[0], c[1], c[2]))]++)] = i;
    }
    auto pairsOf = [&](int q, auto emit) {
        const int i = softSearch_[size_t(q)];
        const auto c = cellOf(x_[size_t(i)]);
        int seen[27], count = 0, found = 0;
        for (int dz = -1; dz <= 1; ++dz)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx) seen[count++] = bucket(c[0] + dx, c[1] + dy, c[2] + dz);
        std::sort(seen, seen + count);
        count = int(std::unique(seen, seen + count) - seen); // two neighbour cells may share a bucket
        for (int k = 0; k < count; ++k)
            for (int e = softCellStart_[size_t(seen[k])]; e < softCellStart_[size_t(seen[k]) + 1]; ++e) {
                const int j = softCellItems_[size_t(e)];
                if (j <= i || object_[size_t(j)] == object_[size_t(i)] || invMass_[size_t(i)] + invMass_[size_t(j)] == 0) continue;
                if (softAsleep_[size_t(i)] && softAsleep_[size_t(j)]) continue; // two fixed ones
                const Vector3 centres = x_[size_t(i)] - x_[size_t(j)];
                const float apart = length(centres), reach = d0 + softTravel(i, dt) + softTravel(j, dt) + 0.25f * r;
                if (apart < reach && apart >= 1e-9f) emit(found++, i, j, centres, apart);
            }
        return found;
    };
    softPairs_.resize(size_t(countResults(n, softCount_, [&](int q) { return pairsOf(q, [](int, int, int, const Vector3&, float) {}); })));
    parallelFor(n, [&](int q) {
        pairsOf(q, [&](int k, int i, int j, const Vector3& centres, float apart) {
            // softBodyOf_ is the last pre-stabilization's: particles appended since (the emitter's) are not in it.
            bool inside = apart < 0.75f * d0;
            if (!inside && !unresolved_.empty() && size_t(j) < softBodyOf_.size())
                inside = std::binary_search(unresolved_.begin(), unresolved_.end(), moverPair({i, j, Vector3(0.0f), 1.0f, d0, false}));
            softPairs_[size_t(softCount_[size_t(q)] + k)] = {i, j, centres / apart, inside ? std::min(apart, d0) : d0};
        });
    }, 64);
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
void ParticleSystem::solveSoftContacts(float dt) {
    const float r = params.particleRadius, mu = params.solidFriction;
    std::vector<Vector3>& u = softMove_; // the particles are at x_ + u (stepSoftBodies)
    for (SoftBodyContact& c : softBodyContacts_) {
        const size_t i = size_t(c.particle);
        c.depth = r - dot((x_[i] - c.point) + u[i], c.normal);
        c.touching = c.depth > 0;
        if (c.touching && !softAsleep_[i]) u[i] += c.normal * c.depth, c.depth = 0;
    }
    const int np = int(softPairs_.size());
    for (SoftPair& c : softPairs_) c.touching = false;
    for (int sweep = 0; sweep < kPairSweeps; ++sweep)
        for (int q = 0; q < np; ++q) {
            SoftPair& c = softPairs_[size_t(sweep % 2 ? np - 1 - q : q)];
            const size_t i = size_t(c.i), j = size_t(c.j);
            const float depth = pushAlong((x_[i] - x_[j]) + (u[i] - u[j]), c.normal, c.target);
        const float wi = softAsleep_[i] ? 0.0f : invMass_[i], wj = softAsleep_[j] ? 0.0f : invMass_[j];
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
            // A sleeping particle is a fixed one: not pushed out, it pushes the body out instead, to
            // leave by the next step (a box lying on a sleeping jelly stays on top: without it, it
            // would sink by g dt^2 every step, nothing in a sleeping body pushing back).
            const float wp = softAsleep_[i] ? 0.0f : invMass_[i];
            contactImpulse(v_[i] - (bodyV_[k] + cross(bodyW_[k], arm)) - c.normal * (c.depth / dt), c.normal, mu,
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
        const float wi = softAsleep_[i] ? 0.0f : invMass_[i], wj = softAsleep_[j] ? 0.0f : invMass_[j];
        contactImpulse(v_[i] - v_[j], c.normal, mu, [&](const Vector3&) { return wi + wj; },
                       [&](const Vector3& J) { v_[i] += J * wi, v_[j] -= J * wj; });
    }
}

// The small steps of all soft bodies, m of h = dt / m, each:
//   1. fly (softFlight);
//   2. the material: one XPBD pass, warm-started (SoftBodySolver.cpp);
//   3. the walls and the obstacle (softWalls), then the rigid bodies and the other soft bodies
//      (solveSoftContacts), each contact pushed out and its approach stopped with friction.
// The body leaves the step with its last small step's velocity: it falls at exactly g, rests with
// none, and sags as its material says. (Contacts solved after the small steps instead, in the
// passes, got no velocity right: the mean over the step lagged gravity - a body fell at 0.63 g -,
// the last small step's left a supported body a phantom 20 mm/s.) The liquid and the cloth meet the
// soft particles in the passes that follow; what they push a soft particle by is added to its
// position and, over dt, to its velocity (finishStep).
// A body's flight in small step s of m (of length h): v += g h, u += v h; a held particle slides to
// where the step puts it.
void ParticleSystem::softFlight(const SoftBody& b, int s, int m, float h) {
    std::vector<Vector3>& u = softMove_;
    const Vector3 g = params.gravity;
    for (int i : b.particles) {
        const size_t k = size_t(i);
        if (invMass_[k] == 0) {
            u[k] = softFlight_[k] * (float(s + 1) / float(m));
            continue;
        }
        v_[k] += g * h;
        u[k] += v_[k] * h;
        softFlight_[k] = u[k];
    }
}

// After the material: v += its move / h; then the walls and the obstacle push the particles out and
// stop their approach with friction (contactImpulse).
void ParticleSystem::softWalls(const SoftBody& b, float h) {
    std::vector<Vector3>& u = softMove_;
    const float mu = params.solidFriction;
    for (int i : b.particles) {
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

void ParticleSystem::stepSoftBodies(float dt) {
    if (softBodies_.empty()) return;
    softBodyBoxes(dt);
    wakeSleepingSoftBodies(dt);
    int m = 1;
    softAll_.clear();
    for (SoftBody& b : softBodies_)
        if (!b.asleep) softAll_.push_back(&b), m = std::max(m, softSmallSteps(b, dt));
    lastSoftSmallSteps_ = m;
    Probe::set("particles/soft bodies asleep", double(softBodies_.size() - softAll_.size()));
    const float h = dt / float(m);
    std::vector<Vector3>& u = softMove_; // every soft particle is at x_ + u (SoftBodySolver.cpp: why)
    for (SoftBody& b : softBodies_)
        for (int i : b.particles) {
            u[size_t(i)] = Vector3(0.0f);
            softAsleep_[size_t(i)] = uint8_t(b.asleep);
            if (b.asleep) v_[size_t(i)] = Vector3(0.0f); // what pushed it too softly to wake it is dropped
            if (invMass_[size_t(i)] == 0) softFlight_[size_t(i)] = p_[size_t(i)] - x_[size_t(i)]; // where the held ones end
        }
    findSoftContacts(dt);
    if (rigid_) {
        const auto& bodies = rigid_->bodies();
        bodyV_.resize(bodies.size());
        bodyW_.resize(bodies.size());
        for (size_t k = 0; k < bodies.size(); ++k) bodyV_[k] = bodies[k].vel, bodyW_[k] = bodies[k].angVel;
    }
    if (softAll_.empty()) { // everything sleeps (p_ = x_ already), but still holds up what lies on it
        solveSoftContacts(dt);
        return;
    }
    for (SoftBody* b : softAll_) {
        // The multipliers are force x h^2: carried from the last small step, rescaled if h changed.
        const float scale = b->multiplierStep > 0 ? sqr(h / b->multiplierStep) : 0.0f;
        if (scale != 1.0f) scaleSoftMultipliers(*b, scale);
        b->multiplierStep = h;
    }
    const int awake = int(softAll_.size());
    for (int s = 0; s < m; ++s) {
        // The bodies share no particle: their flights and walls run a body per thread, their material
        // one colour of all of them per parallel loop.
        parallelFor(awake, [&](int b) { softFlight(*softAll_[size_t(b)], s, m, h); }, 1);
        solveSoftBodies(softAll_, x_, u, invMass_, h, true, softRuns_);
        parallelFor(awake, [&](int b) { softWalls(*softAll_[size_t(b)], h); }, 1);
        solveSoftContacts(dt); // 3. the rigid bodies and the other soft bodies
    }
    for (const SoftBody* b : softAll_)
        for (int i : b->particles) p_[size_t(i)] = x_[size_t(i)] + u[size_t(i)];
}

// Sleeping, as the rigid world's islands do it: a soft body whose particles have all moved slower
// than params.softSleepSpeed (by their move in the step) for params.softSleepTime - not held by the mouse, not in liquid or cloth -
// sleeps: it is not stepped, its particles hold still and meet the other bodies as fixed ones. Its
// multipliers keep the forces that held it, so it wakes exactly in the balance it fell asleep in,
// without a jolt. A lying jelly costs nothing; a scene that has come to rest, next to nothing.
void ParticleSystem::sleepStillSoftBodies(float dt) {
    for (size_t k = 0; k < softBodies_.size(); ++k) {
        SoftBody& b = softBodies_[k];
        if (!params.softSleeping) b.asleep = false, b.stillTime = 0;
        if (b.asleep || !params.softSleeping) continue;
        // How far the particles went in the step, not their velocities: under a rigid body lying on
        // it a still jelly's top particles carry the body's g dt of every step (it gets gravity in
        // its own step, the support in this one) and would never count as still.
        float farthest = 0;
        for (int i : b.particles) farthest = std::max(farthest, length2(p_[size_t(i)] - x_[size_t(i)]));
        const bool held = grabbing() && object_[size_t(grab_.particles[0])] == b.object;
        const bool still = farthest < sqr(params.softSleepSpeed * dt) && !held && !b.touchesOthers;
        b.stillTime = still ? b.stillTime + dt : 0.0f;
        if (b.stillTime < params.softSleepTime) continue;
        b.asleep = true;
        for (int i : b.particles) v_[size_t(i)] = Vector3(0.0f);
        b.sleepRigid.clear();
        if (rigid_ && k < softBoxes_.size()) rigid_->queryBodies(softBoxes_[k], b.sleepRigid);
    }
}

// What wakes a sleeping body: liquid or cloth reaching it (last step), an awake soft body that moved
// in its last step and whose box meets its box, a rigid body in its box moving faster than the
// sleep speed plus twice gravity's kick of a step (what lies on the body at rest moves by g dt
// between its own step and the support), or a change among the rigid bodies in its box (one that
// held it up is gone). Woken, it counts as moving for a step, so it wakes a sleeping stack one body
// per step from where it was touched.
void ParticleSystem::wakeSleepingSoftBodies(float dt) {
    const size_t count = softBodies_.size();
    for (size_t k = 0; k < count; ++k) {
        SoftBody& b = softBodies_[k];
        if (!b.asleep) continue;
        bool wake = b.touchesOthers || !params.softSleeping ||
                    std::find(softWakeObjects_.begin(), softWakeObjects_.end(), b.object) != softWakeObjects_.end();
        for (size_t a = 0; a < count && !wake; ++a)
            wake = !softBodies_[a].asleep && softBodies_[a].stillTime == 0 && softBoxes_[a].overlaps(softBoxes_[k]);
        if (!wake && rigid_) {
            rigid_->queryBodies(softBoxes_[k], softRigidNow_);
            wake = softRigidNow_ != b.sleepRigid;
            const float quick = params.softSleepSpeed + 2.0f * length(params.gravity) * dt;
            for (int r : softRigidNow_) {
                const RigidBody& body = rigid_->bodies()[size_t(r)];
                const float speed = length(body.vel) + length(body.angVel) * body.boundingRadius();
                wake = wake || (body.invMass > 0 && !body.sleeping && speed > quick);
            }
        }
        if (wake) b.asleep = false, b.stillTime = 0;
    }
    softWakeObjects_.clear();
}

size_t ParticleSystem::sleepingSoftBodies() const {
    size_t n = 0;
    for (const SoftBody& b : softBodies_) n += b.asleep ? 1 : 0;
    return n;
}

void ParticleSystem::wakeSoftBodies() {
    for (SoftBody& b : softBodies_) b.asleep = false, b.stillTime = 0;
}

} // namespace rf
