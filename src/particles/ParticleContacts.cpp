// Contacts of ParticleSystem: solid particles against each other, particles against the walls
// and the static mesh, and the two-way XPBD contacts with rigid bodies (Mueller et al. 2020).
//
// Particles against particles follow FleX (Macklin, Mueller, Chentanez, Kim 2014, "Unified
// Particle Physics for Real-Time Applications", ACM TOG 33(4)), in three steps a substep:
//   1. findParticleContacts - the pairs that may touch within the step;
//   2. preStabilizeContacts - bodies found INSIDE each other at the start of the step (made so,
//      pushed in by the mouse, come through at speed) are pulled apart before anything else, the
//      start and the predicted positions alike, so the overlap never becomes a speed (sec. 4.4);
//      they come apart the way out of the bodies, told by the bodies' signed distance fields
//      (sec. 5.1), each soft body moved as a whole;
//   3. solveParticleContacts - the main solve: every pair two hard spheres kept d0 = 2r apart,
//      the upper particle of a pile counting as the lighter (sec. 5.2).
#include "particles/ParticleSystem.h"

#include "core/Parallel.h"
#include "core/Probe.h"

#include <algorithm>

namespace rf {

// An overlap deeper than this share of the spacing d0 is no longer a touch but an intersection -
// the two are inside each other. The main solve leaves an intersection to the pre-stabilization
// (see findIntersections for why a quarter).
constexpr float kIntersectionDepth = 0.25f;
// Overlaps shallower than this share of d0 count as touching (the solvers' round-off).
constexpr float kOverlapTolerance = 0.01f;
// Jacobi passes of the pre-stabilization per substep (it stops sooner when nothing is left).
constexpr int kPreStabilizationPasses = 8;

// How far two hard spheres whose centres are r = x_i - x_j apart must move apart along the unit
// direction n until they are `target` apart; 0 if they already are. It is |r + t n| = target
// solved for t: with a = r . n the part of r along n and s^2 = |r|^2 - a^2 the part across,
// t = sqrt(target^2 - s^2) - a. Along the line of centres (n = r / |r|) it is the familiar
// target - |r|; a pair side by side needs less, a pair that has swapped sides along n needs more.
float pushAlong(const Vector3& r, const Vector3& n, float target) {
    const float r2 = length2(r);
    if (r2 >= target * target) return 0.0f;
    const float along = dot(r, n);
    const float across2 = std::max(0.0f, r2 - along * along);
    return std::sqrt(target * target - across2) - along;
}

void ParticleSystem::findParticleContacts() {
    const int n = int(p_.size());
    const float d0 = spacing(), reach2 = sqr(1.5f * d0);
    // Per particle in parallel, then concatenated in particle order (deterministic).
    std::vector<std::vector<ParticleContact>> local(n);
    parallelFor(n, [&](int i) {
        const int* nb = &nbr_[size_t(i) * kMaxNeighbors];
        for (int m = 0; m < nbrCount_[i]; ++m) {
            const int j = nb[m];
            if (j <= i) continue;                   // each pair once
            if (isFluid(i) && isFluid(j)) continue; // fluid-fluid: the density constraint
            if (invMass_[i] + invMass_[j] == 0) continue;
            if (object_[i] >= 0 && object_[i] == object_[j]) {
                if (phase_[i] == uint8_t(ParticlePhase::Soft)) continue;         // held by shape matching
                if (length2(rest_[i] - rest_[j]) < 2.25f * d0 * d0) continue;     // cloth neighbours at rest
            }
            if (length2(p_[i] - p_[j]) >= reach2) continue; // cannot touch within this substep
            // The normal and the target come later: the pre-stabilization may still move the pair.
            local[i].push_back({i, j, Vector3(0.0f), stackLift(i, j), d0, false});
        }
    });
    contacts_.clear();
    for (auto& l : local) contacts_.insert(contacts_.end(), l.begin(), l.end());
}

// FleX's stiff-stack mass scaling for one pair (Macklin et al. 2014, sec. 5.2, eq. 21): in the
// contacts a particle counts as lighter the higher it sits, m* = m exp(-k h). A correction travels
// one particle per pass, and a pile of soft bodies is dozens of particles tall - with true masses
// it sinks into itself before the support from below arrives. With the scaling a contact lifts
// the upper particle and leaves the lower one in place, as the ground does. Only the difference
// of heights matters, so it is split in half between the two: i's inverse mass is multiplied by
// exp(+k dh / 2), j's by exp(-k dh / 2), dh = how much higher i sits - no huge numbers however
// tall the pile. The heights are the start positions, taken once per step as the paper does
// (Algorithm 1, line 4), and the scaling lives only in the contacts (the masses themselves stay).
// The price is FleX's too: a contact between particles at different heights conserves momentum
// only as a support does.
float ParticleSystem::stackLift(int i, int j) const {
    const float g = length(params.gravity);
    if (g == 0 || params.stackMassScaling == 0 || invMass_[i] == 0 || invMass_[j] == 0) return 1.0f;
    const Vector3 up = params.gravity * (-1.0f / g);
    const float k = params.stackMassScaling / spacing(); // per metre of height
    return std::exp(0.5f * k * dot(x_[i] - x_[j], up));
}

// ---------------------------------------------------------------------------
// Pre-stabilization (Macklin et al. 2014, sec. 4.4, Algorithm 1 lines 10-15)
// ---------------------------------------------------------------------------
// A position solver that pushes apart two particles found inside each other turns every metre of
// the push into a speed of 1/dt: the paper's Fig. 4, a particle started in the ground "pops" out
// of it. Two jellies made a third inside each other flew apart at 0.36 m/s; a jelly started
// 10 cm inside another on the floor was thrown 2.6 m up (+272 J), the stack scaling lifting it a
// little more every step. So an intersection is undone before the main solve, from the start
// positions x, and the predicted positions p are shifted alike: the velocity (p - x) / dt never
// sees it.
// Who moves: FleX pushes single particles; its rigid bodies take their shape back in the main
// loop and the dent is gone. A soft body keeps the dent and springs back from it into the other
// body, and the overlap has become a speed after all. So here a soft body moves as a whole - one
// shift for all its particles, its shape and its elastic energy untouched - and a cloth or
// liquid particle moves alone, as in FleX.
// Which way: the bodies' signed distance fields (sec. 5.1, intersectionNormal) - a deep overlap
// comes apart the way out of the bodies, where the lines between particle centres point anywhere
// and the lattices of two bodies interlock (the paper's Fig. 6).
void ParticleSystem::preStabilizeContacts() {
    findMovers();
    const int intersecting = findIntersections();
    Probe::set("particles/intersecting pairs", intersecting);
    if (intersecting == 0) return;
    for (const SoftBody& body : softBodies_) turnSurfaceNormals(body, x_, restSurfaceNormal_, surfaceNormal_);
    for (ParticleContact& c : contacts_)
        if (c.intersecting) c.normal = intersectionNormal(c.i, c.j);
    // Jacobi passes, each: the contacts push the movers apart, then every shifted mover is pushed
    // back out of the supports (a body on the floor is not pushed into it).
    for (int pass = 0; pass < kPreStabilizationPasses; ++pass) {
        if (!pushMoversApart()) break;
        pushMoversOutOfSupports();
    }
    moveByMovers();
    rememberUnresolved();
}

// Who moves in the pre-stabilization: every particle alone, with its inverse mass, except the
// particles of a soft body - they share one mover, the body (numbered n + its index), with the
// inverse of the body's mass; 0 if any of its particles is pinned or held by the mouse: a held
// body is a support. The shifts start at zero.
void ParticleSystem::findMovers() {
    const int n = int(x_.size()), bodies = int(softBodies_.size());
    moverOf_.resize(size_t(n));
    softBodyOf_.assign(size_t(n), -1);
    moverInvMass_.resize(size_t(n + bodies));
    moverShift_.assign(size_t(n + bodies), Vector3(0.0f));
    moverMove_.resize(size_t(n + bodies));
    moverAsks_.resize(size_t(n + bodies));
    supportPush_.resize(size_t(n));
    for (int i = 0; i < n; ++i) {
        moverOf_[i] = i;
        moverInvMass_[i] = invMass_[i];
    }
    for (int b = 0; b < bodies; ++b) {
        float mass = 0;
        bool held = false;
        for (int i : softBodies_[b].particles) {
            moverOf_[i] = n + b;
            softBodyOf_[i] = b;
            if (invMass_[i] == 0) held = true;
            else mass += 1.0f / invMass_[i];
        }
        moverInvMass_[n + b] = held || mass == 0 ? 0.0f : 1.0f / mass;
    }
}

// A pair of movers as one number that stays the same from step to step: a soft body counts as
// -1 - its index, any other particle as its own index, the smaller of the two first.
long long ParticleSystem::moverPair(const ParticleContact& c) const {
    auto id = [&](int i) { return softBodyOf_[size_t(i)] >= 0 ? -1LL - softBodyOf_[size_t(i)] : (long long)i; };
    const long long offset = 1LL << 30, range = 1LL << 31;
    return (std::min(id(c.i), id(c.j)) + offset) * range + (std::max(id(c.i), id(c.j)) + offset);
}

// Which pairs of movers are inside each other: some pair of their particles overlaps by more than
// kIntersectionDepth at the start of the step, or the last step's pre-stabilization had not yet
// pulled them apart (unresolved_). Every contact of such a pair is the pre-stabilization's -
// pushed apart with its movers here, only kept from going deeper in the main solve. Every other
// contact is an ordinary touch, the main solve's alone.
// Why a quarter spacing: in a pile under load the main solve leaves its pairs up to ~0.2 d0 deep
// (six soft barrels: 0.21 d0 at worst). Such a pair is a touch; lifting a whole barrel out of it
// would move the barrel on every step - sideways as well, wherever the pairs lean - and toppled
// the column within a quarter of a second. Why an intersection is held until it is apart: let go
// at a quarter spacing, the rest would be pushed apart by the main solve, with speed (0.08 m/s for
// two rigid jellies). Returns how many contacts belong to an intersection.
int ParticleSystem::findIntersections() {
    const float d0 = spacing(), deep = (1.0f - kIntersectionDepth) * d0, touching = (1.0f - kOverlapTolerance) * d0;
    intersections_.clear();
    for (const ParticleContact& c : contacts_) {
        const float apart2 = length2(x_[c.i] - x_[c.j]);
        if (apart2 >= touching * touching) continue;
        const long long pair = moverPair(c);
        if (apart2 < deep * deep || std::binary_search(unresolved_.begin(), unresolved_.end(), pair))
            intersections_.push_back(pair);
    }
    unresolved_.clear(); // rebuilt after this step's passes (rememberUnresolved)
    std::sort(intersections_.begin(), intersections_.end());
    intersections_.erase(std::unique(intersections_.begin(), intersections_.end()), intersections_.end());
    if (intersections_.empty()) return 0;
    int count = 0;
    for (ParticleContact& c : contacts_) {
        c.intersecting = std::binary_search(intersections_.begin(), intersections_.end(), moverPair(c));
        count += c.intersecting;
    }
    return count;
}

// The intersections this step's passes have not pulled apart yet (a deep one takes a few steps,
// its pairs being found as the bodies move): they stay intersections in the next step whatever
// their depth. One that has come apart is forgotten - if its bodies go on pressing on each other
// (one resting on the other), that is an ordinary touch again.
void ParticleSystem::rememberUnresolved() {
    const float touching = (1.0f - kOverlapTolerance) * spacing();
    unresolved_.clear();
    for (const ParticleContact& c : contacts_)
        if (c.intersecting && length2(x_[c.i] - x_[c.j]) < touching * touching) unresolved_.push_back(moverPair(c));
    std::sort(unresolved_.begin(), unresolved_.end());
    unresolved_.erase(std::unique(unresolved_.begin(), unresolved_.end()), unresolved_.end());
}

// The way the pre-stabilization pushes particle i of an intersection (j goes the other way).
//  * A cloth or liquid particle has no surface of its own: the pair is two hard spheres, pushed
//    apart along the line between their centres.
//  * Two particles of two soft bodies (Macklin et al. 2014, sec. 5.1): the one nearer to its own
//    body's surface decides (eq. 17). The other particle is pushed out of that body along the
//    surface normal there - the shortest way out - and the deciding one into its body. Deep
//    inside (a particle more than d0 under its surface), that normal is the whole answer. On the
//    surface layer the line between the centres is kept, but mirrored where it would drive the
//    neighbour into the body (eq. 20, after Mueller & Chentanez 2011, "Solid simulation with
//    oriented particles"): a surface particle is a one-sided sphere.
Vector3 ParticleSystem::intersectionNormal(int i, int j) const {
    const Vector3 centres = x_[i] - x_[j];
    const Vector3 line = length2(centres) > 1e-12f ? normalize(centres) : Vector3(0.0f);
    if (!hasSurface(i) || !hasSurface(j)) return line;
    const bool iDecides = surfaceDepth_[i] < surfaceDepth_[j];
    const int k = iDecides ? i : j;
    const Vector3 way = iDecides ? -surfaceNormal_[k] : surfaceNormal_[k]; // as i moves: k into its body, or out of it
    if (surfaceDepth_[k] >= spacing() || length2(line) == 0) return way;
    const float against = dot(line, way);
    return against >= 0 ? line : line - way * (2.0f * against); // mirrored into the free side
}

// One Jacobi pass over the intersections (the paper's constraint averaging, eq. 12): every
// overlapping pair, where the shifts so far have put it, asks to be pushed apart along its normal
// until the two spheres are d0 apart, the push shared between its two movers by their (scaled)
// inverse masses; then every mover moves by the mean of what its pairs asked. A body in a hundred
// pairs moves by their mean, not by their sum, and pushes that disagree (the side faces of two
// boxes inside each other) cancel instead of throwing it about. False when nothing moved.
bool ParticleSystem::pushMoversApart() {
    const float d0 = spacing(), tolerance = kOverlapTolerance * d0;
    std::fill(moverMove_.begin(), moverMove_.end(), Vector3(0.0f));
    std::fill(moverAsks_.begin(), moverAsks_.end(), 0);
    for (const ParticleContact& c : contacts_) {
        if (!c.intersecting) continue;
        const int a = moverOf_[c.i], b = moverOf_[c.j];
        const Vector3 r = (x_[c.i] + moverShift_[a]) - (x_[c.j] + moverShift_[b]);
        const float depth = pushAlong(r, c.normal, d0);
        const float wa = moverInvMass_[a] * c.lift, wb = moverInvMass_[b] / c.lift;
        if (depth <= tolerance || wa + wb == 0) continue; // apart, or two held bodies
        const Vector3 corr = c.normal * (depth / (wa + wb));
        moverMove_[a] += corr * wa;
        moverMove_[b] -= corr * wb;
        ++moverAsks_[a];
        ++moverAsks_[b];
    }
    bool moved = false;
    for (size_t m = 0; m < moverShift_.size(); ++m) {
        if (moverAsks_[m] == 0) continue;
        moverShift_[m] += moverMove_[m] / float(moverAsks_[m]);
        moved = true;
    }
    return moved;
}

// Every shifted mover back out of the supports: each of its particles, where the shift has put
// it, is tested against the walls, the obstacle and the fixed bodies (in parallel), and the mover
// follows the deepest push among them - a body pushed into the floor is lifted by as much as its
// lowest particle went in. Of a soft body only the surface layer is tested: it meets a wall first.
void ParticleSystem::pushMoversOutOfSupports() {
    const int n = int(x_.size());
    const float d0 = spacing();
    prepareBodyQuery(false);
    parallelFor(n, [&](int i) {
        const int m = moverOf_[i];
        supportPush_[i] = Vector3(0.0f);
        if (moverInvMass_[m] == 0 || length2(moverShift_[m]) == 0) return;
        if (m >= n && surfaceDepth_[i] >= d0) return; // deep inside a body moved as a whole
        Vector3 q = x_[i] + moverShift_[m];
        const Vector3 before = q;
        pushOutOfSupports(q);
        supportPush_[i] = q - before;
    });
    std::fill(moverMove_.begin(), moverMove_.end(), Vector3(0.0f));
    for (int i = 0; i < n; ++i) {
        Vector3& deepest = moverMove_[moverOf_[i]];
        if (length2(supportPush_[i]) > length2(deepest)) deepest = supportPush_[i];
    }
    for (size_t m = 0; m < moverShift_.size(); ++m) moverShift_[m] += moverMove_[m];
}

// The supports that never move push a point out: the domain walls, the obstacle mesh and the
// fixed rigid bodies. A moving body is left to the main solve, where it is pushed back.
void ParticleSystem::pushOutOfSupports(Vector3& p) {
    collideWallsAndMesh(p, p);
    if (!rigid_) return;
    const float r = params.particleRadius;
    const auto& bodies = rigid_->bodies();
    std::vector<int>& candidates = bodyCandidates_[size_t(ThreadPool::workerIndex())];
    rigid_->queryBodies(AABB(p - Vector3(bodyReach_), p + Vector3(bodyReach_)), candidates);
    for (int b : candidates) {
        const RigidBody& body = bodies[b];
        if (body.invMass > 0 || length2(p - body.pos) > sqr(body.boundingRadius() + r)) continue;
        Vector3 n;
        const float sd = body.signedDistance(p, n);
        if (sd < r) p += n * (r - sd);
    }
}

// The shifts of the pre-stabilization into the start and the predicted positions alike.
void ParticleSystem::moveByMovers() {
    parallelFor(int(x_.size()), [&](int i) {
        const Vector3& shift = moverShift_[moverOf_[i]];
        x_[i] += shift;
        p_[i] += shift;
    });
}

// ---------------------------------------------------------------------------
// The main solve
// ---------------------------------------------------------------------------
// Every pair, from where the step now begins: two hard spheres, pushed apart along the line
// between their centres - fixed for the step, the pair may not swap sides within it, even if
// other constraints (a cloth snapping back, a heavy body pushing) move one particle past the
// other - until they are d0 apart. A pair of an intersection that still overlaps is only kept
// from going deeper (its target is the distance it starts at): the pre-stabilization undoes it
// in the next steps, without a speed.
// Why not the signed distance normals here as well: at the rim of a soft barrel squashed under a
// load the surface normal points sideways, and the one-sided contacts of eq. 20 turned the
// support of the barrel above into a sideways push - a column of six toppled within 0.25 s,
// where two hard spheres keep it standing (docs/03-particles.md, "Contacts").
void ParticleSystem::setMainSolveTargets() {
    const float d0 = spacing(), touching = (1.0f - kOverlapTolerance) * d0;
    for (ParticleContact& c : contacts_) {
        Vector3 centres = x_[c.i] - x_[c.j];
        if (length2(centres) < 1e-12f) centres = p_[c.i] - p_[c.j];
        const float apart = length(centres);
        c.normal = apart > 1e-6f ? centres / apart : Vector3(0.0f);
        c.target = c.intersecting && apart < touching ? apart : d0;
    }
}

// One pass of the main solve (in every solid pass): Gauss-Seidel, each contact resolved at once,
// so a stack of contacts (a body resting on cloth resting on ...) passes its correction through
// within one sweep. The pair's (scaled) inverse masses share the correction.
// Limitation of all position-based solvers (FleX included): with a large mass ratio between
// touching particles (a heavy body on a very light cloth, beyond ~1:10) the light side takes
// almost the whole correction and the support converges too slowly - use realistic materials
// (canvas ~1-2 kg/m^2 under foam-like bodies) or more solid iterations.
// Two solid particles (soft body, cloth) in contact also rub (Macklin et al. 2014, sec. 6.1, after
// Coulomb): their relative slide along the contact over the step is taken back whole while it is
// under mu times the push the contact just gave - they stick - and by mu times the push otherwise.
// Without it the lattices of two bodies nest like eggs in a tray and slide off each other on the
// slopes of the hollows: the top barrel of a stack of six rubber barrels (10 MPa) slid 0.3 m sideways
// in a second and fell. Liquid particles do not rub (their viscosity is XSPH's).
void ParticleSystem::solveParticleContacts() {
    const float mu = params.solidFriction;
    for (const ParticleContact& c : contacts_) {
        if (isSoft(c.i) && isSoft(c.j)) continue; // two soft bodies: their own step's (SoftBodySolver.cpp)
        const float depth = pushAlong(p_[c.i] - p_[c.j], c.normal, c.target);
        if (depth <= 0) continue;
        const float wi = invMass_[c.i] * c.lift, wj = invMass_[c.j] / c.lift;
        const Vector3 corr = c.normal * (depth / (wi + wj));
        p_[c.i] += corr * wi;
        p_[c.j] -= corr * wj;
        if (mu <= 0 || isFluid(c.i) || isFluid(c.j)) continue;
        const Vector3 slide = (p_[c.i] - x_[c.i]) - (p_[c.j] - x_[c.j]);
        const Vector3 tangent = slide - c.normal * dot(slide, c.normal);
        const float length = std::sqrt(length2(tangent)), limit = mu * depth;
        if (length < 1e-12f) continue;
        // By the true masses: the stack's mass scaling is for the support along the normal. Rubbing
        // by the scaled ones moved the upper particle more than the lower one's share, sideways
        // momentum out of nothing - the stack of barrels gained 23 J and threw the top one up 0.29 m.
        const float ui = invMass_[c.i], uj = invMass_[c.j];
        const Vector3 back = tangent * (length <= limit ? 1.0f : limit / length) * (1.0f / (ui + uj));
        p_[c.i] -= back * ui;
        p_[c.j] += back * uj;
    }
}

void ParticleSystem::solveBodyContacts(float dt) {
    if (!rigid_ || bodyShift_.empty()) return;
    const auto& bodies = rigid_->bodies();
    // Contacts of each body, solved one after another (Gauss-Seidel): every contact sees the body
    // already moved by the ones before it.
    std::vector<std::vector<int>> contacts(bodies.size());
    for (int i = 0; i < int(contactBody_.size()); ++i)
        if (contactBody_[i] >= 0) {
            contacts[contactBody_[i]].push_back(i);
            contactBody_[i] = -1;
        }
    for (size_t b = 0; b < bodies.size(); ++b) {
        if (contacts[b].empty()) continue;
        const RigidBody& body = bodies[b];
        const Vector3 shift0 = bodyShift_[b], turn0 = bodyTurn_[b]; // the pose the contacts were found at
        for (int i : contacts[b]) {
            const Vector3 n = contactNormal_[i], c = contactPoint_[i];
            const Vector3 rc = c - (body.pos + bodyShift_[b]);
            // Penetration now: what it was, minus how far the body has moved away from it since.
            const Vector3 moved = (bodyShift_[b] - shift0) + cross(bodyTurn_[b] - turn0, rc);
            const float depth = contactDepth_[i] + dot(moved, n);
            if (depth <= 0) continue;
            // Generalized inverse masses (Mueller et al. 2020, eq. 2-3): the particle, the body at c.
            const float wp = invMass_[i];
            const Vector3 rn = cross(rc, n);
            const float wb = body.invMass + dot(rn, body.applyInvInertiaWorld(rn));
            const float lambda = depth / (wp + wb); // [kg m]
            p_[i] += n * (lambda * wp);
            bodyShift_[b] -= n * (lambda * body.invMass);
            bodyTurn_[b] -= body.applyInvInertiaWorld(rn) * lambda;
            // Friction: the tangential slip of the particle against the body, shared the same way.
            if (params.wallFriction > 0 && wp > 0) {
                const Vector3 bodyMove = body.velocityAt(c) * dt + (bodyShift_[b] + cross(bodyTurn_[b], rc));
                const Vector3 slip = (p_[i] - x_[i]) - bodyMove;
                const Vector3 t = (slip - n * dot(slip, n)) * params.wallFriction;
                const Vector3 lt = t / (wp + wb);
                p_[i] -= lt * wp;
                bodyShift_[b] += lt * body.invMass;
                bodyTurn_[b] += body.applyInvInertiaWorld(cross(rc, lt));
            }
        }
    }
}

void ParticleSystem::prepareBodyQuery(bool coupled) {
    if (!rigid_) return;
    rigid_->updateWorldTree();
    const auto& bodies = rigid_->bodies();
    // The query box must reach every body the exact test could accept: a particle within r of the
    // surface, the body shifted by bodyShift_ and turned by bodyTurn_ (a point of the surface moves
    // at most |turn| times the bounding radius).
    float reach = params.particleRadius;
    if (coupled && bodyShift_.size() == bodies.size())
        for (size_t b = 0; b < bodies.size(); ++b)
            reach = std::max(reach, params.particleRadius + length(bodyShift_[b]) + length(bodyTurn_[b]) * bodies[b].boundingRadius());
    bodyReach_ = reach;
    if (bodyCandidates_.size() != size_t(ThreadPool::instance().threadCount()))
        bodyCandidates_.resize(size_t(ThreadPool::instance().threadCount()));
}

// Domain walls and the static obstacle mesh: the point is pushed out, and the mesh's friction
// takes back part of the slip along it since `start`.
void ParticleSystem::collideWallsAndMesh(Vector3& p, const Vector3& start) const {
    const float r = params.particleRadius;
    // Domain walls.
    p = vmax(domain_.lo + Vector3(r), vmin(p, domain_.hi - Vector3(r)));

    // Static obstacle mesh.
    if (mesh_ && !mesh_->empty()) {
        AABB mb = mesh_->bounds();
        mb.lo -= Vector3(r);
        mb.hi += Vector3(r);
        if (mb.contains(p)) {
            ClosestHit hit;
            if (mesh_->closestPoint(p, 3 * r, hit) && hit.signedDistance < r) {
                Vector3 n = hit.normal;
                p += n * (r - hit.signedDistance);
                Vector3 dx = p - start;
                p -= (dx - n * dot(dx, n)) * params.wallFriction;
            }
        }
    }
}

void ParticleSystem::collide(int i, Vector3& p, const Vector3& start, bool record, float dt) {
    const float r = params.particleRadius;
    collideWallsAndMesh(p, start);

    // Rigid bodies. Fixed ones (and every body in the cloth's small steps, record = false): the
    // particle is pushed out. Movable ones in the main passes: the contact is only recorded here
    // and solved together with the body in solveBodyContacts. A body is met where this substep
    // has moved it: its start pose shifted by bodyShift_ and turned by bodyTurn_ (small angle).
    if (rigid_) {
        auto& bodies = rigid_->bodies();
        const bool coupled = record && bodyShift_.size() == bodies.size();
        // The bodies whose boxes reach this particle, from the world tree (prepareBodyQuery set
        // the reach and refreshed the tree before the pass); in index order, as the full loop was.
        std::vector<int>& candidates = bodyCandidates_[size_t(ThreadPool::workerIndex())];
        rigid_->queryBodies(AABB(p - Vector3(bodyReach_), p + Vector3(bodyReach_)), candidates);
        for (int b : candidates) {
            const RigidBody& body = bodies[b];
            const Vector3 shift = coupled ? bodyShift_[b] : Vector3(0.0f), turn = coupled ? bodyTurn_[b] : Vector3(0.0f);
            const Vector3 centre = body.pos + shift;
            if (length2(p - centre) > sqr(body.boundingRadius() + r)) continue;
            Vector3 n;
            const float sd = body.signedDistance(body.pos + (p - centre) - cross(turn, p - centre), n);
            if (sd >= r) continue;
            n = normalize(n + cross(turn, n));
            if (coupled && body.invMass > 0) {
                contactBody_[i] = b;
                contactNormal_[i] = n;
                contactDepth_[i] = r - sd;
                contactPoint_[i] = p - n * r;
                continue;
            }
            p += n * (r - sd);
            const Vector3 rel = (p - start) - body.velocityAt(p) * dt;
            p -= (rel - n * dot(rel, n)) * params.wallFriction;
        }
    }
}

} // namespace rf
