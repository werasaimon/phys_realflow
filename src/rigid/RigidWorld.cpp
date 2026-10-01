// RigidWorld: the bodies, the forces on them and the step of the world. The other parts live in
// ContactSolver.cpp, Islands.cpp, ShockPropagation.cpp, Grab.cpp, Joints.cpp, TimeOfImpact.cpp
// and XpbdSolver.cpp.
#include "rigid/RigidWorld.h"

#include "core/Parallel.h"
#include "core/Probe.h"

#include <algorithm>
#include <chrono>

namespace rf {

// ---------------------------------------------------------------------------
// Body
// ---------------------------------------------------------------------------

float RigidBody::signedDistance(const Vector3& p, Vector3& normal) const {
    Matrix3x3 R = rotation();
    Vector3 nl;
    float d = shape->signedDistance(R.transposed() * (p - pos), nl);
    normal = R * nl;
    return d;
}

// ---------------------------------------------------------------------------
// World setup
// ---------------------------------------------------------------------------
int RigidWorld::addBody(std::shared_ptr<const ConvexShape> shape, const Vector3& pos, const Quaternion& rot, float density,
                        const Vector3& color) {
    RigidBody b;
    b.shape = std::move(shape);
    b.pos = pos;
    b.rot = rot;
    b.color = color;
    if (density > 0) {
        b.mass = density * b.shape->volume();
        b.invMass = 1.0f / b.mass;
        Vector3 I = b.shape->unitInertia() * b.mass;
        b.invInertiaLocal = Vector3(1.0f / I.x, 1.0f / I.y, 1.0f / I.z);
    } else {
        b.mass = 0;
        b.invMass = 0;
        b.invInertiaLocal = Vector3(0.0f);
    }
    b.updateInertia();
    return placeBody(std::move(b));
}

// A new body goes into the slot a destroyed one left (the last freed first), else at the end.
int RigidWorld::placeBody(RigidBody&& b) {
    if (freeBodies_.empty()) {
        bodies_.push_back(std::move(b));
        return int(bodies_.size()) - 1;
    }
    const int slot = freeBodies_.back();
    freeBodies_.pop_back();
    bodies_[size_t(slot)] = std::move(b);
    return slot; // its world-tree leaf is inserted again by the next updateWorldTree()
}

// Removes one body and leaves its slot as an inert tombstone: no mass, no velocity, parked far
// outside the world at a spot of its own (two tombstones never touch each other), out of the
// world tree. Everything that named the body - cached contacts, manifolds, joints, the mouse
// joint, sleeping records - is dropped, so the next step does not see it at all.
void RigidWorld::destroyBody(int i) {
    if (!isAlive(i)) return;
    forgetBody(i);
    RigidBody& b = bodies_[size_t(i)];
    b.alive = false;
    b.mass = b.invMass = 0;
    b.invInertiaLocal = Vector3(0.0f);
    b.vel = b.angVel = b.force = b.torque = b.biasVel = b.biasAngVel = Vector3(0.0f);
    b.sleeping = false;
    b.sleepTimer = 0;
    b.sleepIsland = -1;
    b.pos = b.prevPos = Vector3(1e6f + 1000.0f * float(i), 1e6f, 1e6f);
    b.rot = b.prevRot = Quaternion();
    b.updateInertia();
    if (size_t(i) < treeProxies_.size() && treeProxies_[size_t(i)] >= 0) {
        worldTree_.remove(treeProxies_[size_t(i)]);
        treeProxies_[size_t(i)] = -1;
    }
    freeBodies_.push_back(i);
}

void RigidWorld::forgetBody(int i) {
    // The bodies it touched may rest on it: they are woken so they notice it is gone (a sleeping
    // stack keeps the cached contacts of its pairs, so the cache names every neighbour).
    std::vector<int> neighbours;
    for (auto it = cache_.begin(); it != cache_.end();) {
        const int a = keyA(it->first), b = keyB(it->first);
        if (a == i || b == i) {
            const int other = a == i ? b : a;
            if (other >= 0) neighbours.push_back(other);
            it = cache_.erase(it);
        } else {
            ++it;
        }
    }
    auto names = [i](const Manifold& m) { return m.a == i || m.b == i; };
    manifolds_.erase(std::remove_if(manifolds_.begin(), manifolds_.end(), names), manifolds_.end());
    auto jointNames = [i](const std::unique_ptr<Joint>& j) { return j->a == i || j->b == i; };
    joints_.erase(std::remove_if(joints_.begin(), joints_.end(), jointNames), joints_.end());
    auto frozenNames = [i](const Frozen& f) { return f.body == i; };
    frozen_.erase(std::remove_if(frozen_.begin(), frozen_.end(), frozenNames), frozen_.end());
    held_.erase(std::remove_if(held_.begin(), held_.end(), frozenNames), held_.end());
    auto xNames = [i](const XContact& c) { return c.a == i || c.b == i; };
    xcontacts_.erase(std::remove_if(xcontacts_.begin(), xcontacts_.end(), xNames), xcontacts_.end());
    if (grab_.body == i) grab_ = GrabJoint();
    if (size_t(i) < ccdClamped_.size()) ccdClamped_[size_t(i)] = 0;
    for (int n : neighbours)
        if (n < int(bodies_.size()) && bodies_[size_t(n)].sleeping) wake(n);
}

int RigidWorld::addSphere(const Vector3& pos, float r, float density, const Vector3& color) {
    return addBody(std::make_shared<SphereShape>(r), pos, Quaternion(), density, color);
}

int RigidWorld::addBox(const Vector3& pos, const Vector3& h, const Quaternion& rot, float density, const Vector3& color) {
    return addBody(std::make_shared<BoxShape>(h), pos, rot, density, color);
}

int RigidWorld::addConvex(const TriMesh& mesh, const Vector3& pos, const Quaternion& rot, float density, const Vector3& color) {
    auto hull = std::make_shared<ConvexHullShape>(mesh);
    // The hull's local frame is its principal frame: world orientation = rot * principal.
    Quaternion q = (rot * Quaternion::fromMatrix3x3(hull->principalRotation())).normalized();
    return addBody(hull, pos, q, density, color);
}

int RigidWorld::addCompound(std::shared_ptr<const CompoundShape> shape, const Vector3& pos, const Quaternion& rot, float density,
                            const Vector3& color) {
    Quaternion q = (rot * Quaternion::fromMatrix3x3(shape->principalRotation())).normalized();
    return addBody(std::move(shape), pos, q, density, color);
}

std::vector<RigidWorld::DebugContact> RigidWorld::debugContacts() const {
    std::vector<DebugContact> out;
    for (const Manifold& m : manifolds_) {
        const float friction = length(m.t1 * m.jt1 + m.t2 * m.jt2); // one friction for the whole patch
        for (const SolverPoint& p : m.points) out.push_back({m.a, m.b, p.position, p.normal, p.depth, p.jn, friction});
    }
    return out;
}

// The solver gives body A of a manifold +n jn per point and the friction t1 jt1 + t2 jt2, body B the
// opposite (ContactSolver.cpp): so the impulse on a body is that sum with the sign of its side.
Vector3 RigidWorld::contactImpulseOn(int body, int& points) const {
    Vector3 J(0.0f);
    points = 0;
    for (const Manifold& m : manifolds_) {
        const float side = m.a == body ? 1.0f : (m.b == body ? -1.0f : 0.0f);
        if (side == 0) continue;
        for (const SolverPoint& p : m.points) J += p.normal * (p.jn * side);
        J += (m.t1 * m.jt1 + m.t2 * m.jt2) * side;
        points += m.points.size();
    }
    return J;
}

float RigidWorld::deepestPenetration() const {
    float deepest = 0;
    for (const Manifold& m : manifolds_)
        for (const SolverPoint& p : m.points) deepest = std::max(deepest, p.depth);
    return deepest;
}

// Outside impulses wake a sleeping body, but do not reset the sleep timer of an awake one: a body
// that stays slow despite a steady trickle of small impulses (gas pressure, fluid) can still fall
// asleep; a real push makes it fast and resets the timer by itself.
void RigidWorld::applyExternalImpulse(int i, const Vector3& J, const Vector3& p) {
    if (length2(J) > 0 && bodies_[i].sleeping) wake(i);
    applyImpulse(i, J, p);
}

void RigidWorld::applyVelocityChange(int i, const Vector3& dv, const Vector3& dw) {
    RigidBody& b = bodies_[i];
    if (b.invMass == 0) return;
    if (b.sleeping) wake(i);
    b.vel += dv;
    b.angVel += dw;
}

void RigidWorld::applyExternalWrench(int i, const Vector3& J, const Vector3& L) {
    RigidBody& b = bodies_[i];
    if (b.invMass == 0) return;
    if (b.sleeping) wake(i);
    b.vel += J * b.invMass;
    b.angVel += b.applyInvInertiaWorld(L);
}

void RigidWorld::clear() {
    bodies_.clear();
    freeBodies_.clear();
    worldTree_.clear();
    treeProxies_.clear();
    joints_.clear();
    grab_ = GrabJoint();
    cache_.clear();
    manifolds_.clear();
    pairs_.clear();
    colors_.clear();
    levels_.clear();
    frozen_.clear();
    held_.clear();
    islandParent_.clear();
    ccdClamped_.clear();
    ccdHits_ = 0;
    ccdDiagnostics_ = CcdDiagnostics();
    contactCount_ = 0;
    lastStepResult_ = RigidStepResult();
    trialMotorWork_ = motorWorkTotal_ = 0;
    stepInitial_.reset();
    stepTrial_.reset();
    // XPBD keeps its contacts between collision passes: they belong to the old bodies.
    xcontacts_.clear();
    substepCounter_ = 0;
}

void RigidWorld::updateWorldTree() const {
    // New bodies get a leaf; every body's leaf follows its box (the tree changes only when a body
    // leaves the fat box of its leaf). A destroyed body has no leaf (-1); a body placed into its
    // slot gets one again.
    // A sleeping body has not moved since it fell asleep (anything that moves it wakes it), so its
    // leaf is left alone - a pile at rest costs the tree nothing, as in Box2D.
    treeProxies_.resize(bodies_.size(), -1);
    for (size_t i = 0; i < bodies_.size(); ++i) {
        if (!bodies_[i].alive || (bodies_[i].sleeping && treeProxies_[i] >= 0)) continue;
        if (treeProxies_[i] < 0) treeProxies_[i] = worldTree_.insert(bodies_[i].worldBounds(), int(i));
        else worldTree_.update(treeProxies_[i], bodies_[i].worldBounds());
    }
}

void RigidWorld::queryBodies(const AABB& box, std::vector<int>& out) const {
    worldTree_.query(box, out);
    std::sort(out.begin(), out.end()); // index order: the callers' loops stay deterministic
}

bool RigidWorld::raycast(const Vector3& o, const Vector3& d, float maxT, int& body, float& t, Vector3& normal) const {
    body = -1;
    t = maxT;
    Vector3 inv(1.0f / (std::fabs(d.x) > 1e-12f ? d.x : 1e-12f), 1.0f / (std::fabs(d.y) > 1e-12f ? d.y : 1e-12f),
             1.0f / (std::fabs(d.z) > 1e-12f ? d.z : 1e-12f));
    updateWorldTree();
    worldTree_.raycast(o, d, maxT, rayCandidates_);
    std::sort(rayCandidates_.begin(), rayCandidates_.end()); // the nearest wins; ties go to the lower index, as before
    for (int i : rayCandidates_) {
        const RigidBody& b = bodies_[i];
        if (b.worldBounds().rayHit(o, inv, t) == kInf) continue;
        Matrix3x3 Rt = b.rotation().transposed();
        float ti;
        Vector3 nl;
        if (b.shape->raycast(Rt * (o - b.pos), Rt * d, t, ti, nl) && ti < t) {
            t = ti;
            body = i;
            normal = b.rotation() * nl;
        }
    }
    return body >= 0;
}

// ---------------------------------------------------------------------------
// Joints
// ---------------------------------------------------------------------------

float RigidWorld::kineticEnergy() const {
    double e = 0;
    for (const RigidBody& b : bodies_) {
        if (b.invMass == 0) continue;
        e += 0.5 * b.mass * length2(b.vel);
        Vector3 wl = b.rotation().transposed() * b.angVel;
        Vector3 I(1 / b.invInertiaLocal.x, 1 / b.invInertiaLocal.y, 1 / b.invInertiaLocal.z);
        e += 0.5 * dot(wl, I * wl);
    }
    return float(e);
}

// Before the iterations: the manifolds in bottom-up order, their colours, and each one prepared.
void RigidWorld::prepare(float dt) {
    contactCount_ = 0;
    sortManifoldsBottomUp();
    for (const Manifold& m : manifolds_) contactCount_ += size_t(m.points.size());
    buildColors();
    forEachManifold([&](Manifold& m) { prepareManifold(m, dt); });
    warmStats_ = {int(contactCount_), 0}; // summed after the parallel preparation, in manifold order
    for (const Manifold& m : manifolds_) warmStats_.matched += m.warmMatched;
}

// Gauss-Seidel order: bottom-up (along gravity), so support propagates through stacks within the
// first iterations instead of one level per iteration. A manifold is a big record (~0.6 KB: four
// points, the friction patch, the block matrix), and sorting the records themselves moved them
// around log n times - 12 ms a frame for a thousand cubes. So the records stay where the narrow
// phase put them, and the order is a list of their indices (solveOrder_) plus each one's place in
// it (solveRank_); whoever needs the order - the colouring, the bounces, the shock pass - reads it
// through them:
//   1. each manifold's height is computed once, as the sort key;
//   2. the indices are sorted by it - ties by index, which gives exactly the order of a stable sort
//      without its temporary buffer;
//   3. each manifold's rank in that order.
void RigidWorld::sortManifoldsBottomUp() {
    Vector3 up = normalize(-params.gravity);
    if (length2(up) < 0.5f) up = Vector3(0, 1, 0);
    const int n = int(manifolds_.size());
    sortKeys_.resize(size_t(n));
    solveOrder_.resize(size_t(n));
    for (int i = 0; i < n; ++i) { // 1.
        const Manifold& m = manifolds_[size_t(i)];
        const float h = dot(bodies_[m.a].pos, up);
        sortKeys_[size_t(i)] = m.b >= 0 ? std::min(h, dot(bodies_[m.b].pos, up)) : h - 1e3f; // static contacts first
        solveOrder_[size_t(i)] = i;
    }
    std::sort(solveOrder_.begin(), solveOrder_.end(), [&](int x, int y) { // 2.
        const float kx = sortKeys_[size_t(x)], ky = sortKeys_[size_t(y)];
        return kx < ky || (kx == ky && x < y);
    });
    solveRank_.resize(size_t(n)); // 3.
    for (int r = 0; r < n; ++r) solveRank_[size_t(solveOrder_[size_t(r)])] = r;
}

void RigidWorld::solve() {
    forEachManifold([&](Manifold& m) { solveManifold(m); });
}

void RigidWorld::stepUnchecked(float dt, float correctionDt) {
    correctionDt_ = correctionDt;
    ccdHits_ = 0;
    ccdDiagnostics_ = CcdDiagnostics();
    if (bodies_.empty()) return;
    if (params.solver == RigidSolver::XPBD) {
        // XPBD has no islands: clear old sleep flags so switching solvers cannot freeze moving bodies.
        for (RigidBody& b : bodies_) {
            b.sleeping = false;
            b.sleepTimer = 0;
            b.sleepIsland = -1;
        }
        stepXPBD(dt);
        return;
    }
    // One step of the impulse solver, as a list of what happens to the bodies in order. Every
    // named step is a member below; the probe measures the time and the memory of each stage.
    // The timers follow one another and never nest, because the Laboratory's profiler adds them
    // up: integrate (the kick, later the poses), collide, solve, ccd, islands.
    lastDt_ = dt;
    const auto tStart = std::chrono::steady_clock::now();
    long long allocations = Probe::allocations.load(); // the census: which stage churns memory
    auto countAllocations = [&](const char* channel) {
        const long long now = Probe::allocations.load();
        Probe::add(channel, double(now - allocations));
        allocations = now;
    };
    {
        Probe::Timer t("rigid/integrate ms");
        beginStep();
        integrateVelocities(dt);
    }
    {
        Probe::Timer t("rigid/collide ms");
        collide();
        if (params.sleeping && updateIslands(false, dt)) collide(); // woken island: contacts among its bodies
    }
    countAllocations("memory/rigid collide");
    {
        Probe::Timer t("rigid/solve ms"); // prepare, the iterations, the bounces, the shock pass
        solveContacts(dt);
    }
    countAllocations("memory/rigid solve");
    if (Probe::drawEnabled()) {
        Probe::Timer t("rigid/debug draw ms");
        drawDebug();
    }
    {
        Probe::Timer t("rigid/solve ms"); // the contact cache for the next warm start, the resting damping
        rememberContactImpulses();
        countAllocations("memory/rigid cache");
        dampRestingBodies(dt);
    }
    {
        Probe::Timer t("rigid/integrate ms");
        integratePoses(dt);
    }
    finishStep(dt, tStart); // the ccd, the joints' positions, the islands: each has its own timer
    countAllocations("memory/rigid integrate"); // poses, CCD, joints, islands
}

// Before anything moves: the inertia tensors follow the orientations (which may have been edited
// from outside), the sleep flags obey the parameter, a grabbed body is awake, sleeping bodies are
// frozen into static ones for this step, and every body remembers where the step started (for the
// continuous collision and the joints' position stage).
void RigidWorld::beginStep() {
    for (RigidBody& b : bodies_) // orientation may have been edited externally (a sleeper has not turned)
        if (!b.sleeping) b.updateInertia();
    if (!params.sleeping)
        for (RigidBody& b : bodies_) b.sleeping = false;
    if (grab_.active && grab_.body < int(bodies_.size())) wake(grab_.body);
    freezeSleepers();
    for (RigidBody& b : bodies_) {
        b.prevPos = b.pos;
        b.prevRot = b.rot;
    }
}

// Newton's second law for one step, as a kick: gravity and the applied forces change the linear
// velocity, the applied torques the angular one. The tumble of a body with unequal moments of
// inertia (Euler's equations) is part of its free turn, later in the step (turnFreely,
// FreeRotation.cpp). Symplectic Euler: velocities first, positions later.
void RigidWorld::integrateVelocities(float dt) {
    for (RigidBody& b : bodies_) {
        if (b.invMass == 0) continue;
        b.vel += (params.gravity + b.force * b.invMass) * dt;
        b.angVel += b.applyInvInertiaWorld(b.torque) * dt;
    }
}

// The contact and joint solve at velocity level: prepare the manifolds (effective masses, warm
// start), then Gauss-Seidel iterations over contacts, joints and the mouse joint, then the
// bounces, then the shock propagation pass that lets tall stacks stand. Its time is timings_.solve.
void RigidWorld::solveContacts(float dt) {
    auto ts = std::chrono::steady_clock::now();
    if (params.measureMotorWork) for (auto& j : joints_) j->captureMotorMotion(bodies_);
    prepare(dt);
    prepareGrab(dt);
    for (auto& j : joints_) j->prepare(bodies_, dt, params.warmStarting, params.measureMotorWork);
    for (int it = 0; it < params.iterations; ++it) {
        solve();
        for (auto& j : joints_) j->solveVelocity(bodies_);
        solveGrab(dt);
    }
    // Bounces before the shock pass: it is one-sided with its own accumulators, so it cannot take a
    // separation back, but it does absorb the downward half of a bounce inside a stack.
    applyRestitution();
    if (params.shockPropagation && !manifolds_.empty()) propagateShock();
    trialMotorWork_ = 0;
    if (params.measureMotorWork) for (auto& j : joints_) {
        j->finishMotorWork(bodies_); trialMotorWork_ += j->motorWork();
    }
    timings_.solve = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - ts).count();
}

// The contact cache for the next step's warm start. Pairs of sleeping bodies produce no
// manifolds; keep their last impulses so a woken stack carries its weight immediately instead of
// sagging and being thrown apart. The cache is pruned in place (erase allocates nothing) and a
// persistent pair keeps its node: only a new pair costs an allocation - rebuilding the map every
// step was a thousand allocations per frame for a sleeping tower.
//   1. this step's manifolds are written into the cache, stamped with this step's number;
//   2. one pass over the cache drops every pair this step did not refresh - unless both its bodies
//      sleep (sleeping pairs make no manifolds, and a woken stack needs their impulses).
// (It took two passes over the whole cache before - in a sleeping pile of a thousand cubes that
// is thousands of pairs twice a substep, for the few that changed.)
void RigidWorld::rememberContactImpulses() {
    ++cacheStamp_;
    for (const Manifold& m : manifolds_) { // 1.
        CachedPair& c = cache_[key(m.a, m.b, m.sub)];
        c.stamp = cacheStamp_;
        c.points = m.points;
        c.friction = m.t1 * m.jt1 + m.t2 * m.jt2;
        c.twist = m.jtwist;
        c.roll = m.jroll;
        c.locked = m.locked;
        c.lockRef = m.lockRef;
        c.lock = m.jlock;
    }
    auto asleep = [&](int body) { // static (the walls: body < 0) or sleeping
        return body < 0 || (body < int(bodies_.size()) && (bodies_[body].sleeping || bodies_[body].mass <= 0));
    };
    for (auto it = cache_.begin(); it != cache_.end();) { // 2.
        const uint64_t k = it->first;
        const int a = keyA(k), b = keyB(k);
        const bool keep = it->second.stamp == cacheStamp_ || (a < int(bodies_.size()) && asleep(a) && asleep(b));
        it = keep ? std::next(it) : cache_.erase(it);
    }
}

// Resting-contact damping (cf. Bullet's additional damping): bodies that touch something and
// move slower than the thresholds get an opposing velocity change, capped so that it only eats
// the residual jitter and never real motion. Off by default (params.restDamping = 0).
void RigidWorld::dampRestingBodies(float dt) {
    if (params.restDamping <= 0) return;
    std::vector<char>& touching = touching_; // kept between steps
    touching.assign(bodies_.size(), 0);
    for (const Manifold& m : manifolds_) {
        bool active = false;
        for (const SolverPoint& p : m.points) active |= p.jn > 0;
        if (!active) continue;
        touching[m.a] = 1;
        if (m.b >= 0) touching[m.b] = 1;
    }
    const float cap = params.restDamping * dt;
    for (size_t i = 0; i < bodies_.size(); ++i) {
        RigidBody& b = bodies_[i];
        if (!touching[i] || b.invMass == 0) continue;
        float v = length(b.vel), w = length(b.angVel);
        if (v < params.restLinearThreshold) b.vel *= v > cap ? 1.0f - cap / v : 0.0f;
        float wcap = cap / std::max(b.boundingRadius(), 1e-3f);
        if (w < params.restAngularThreshold) b.angVel *= w > wcap ? 1.0f - wcap / w : 0.0f;
    }
}

// Positions and orientations follow the solved velocities (plus the pseudo velocities of the
// split impulse, which recover penetration without adding momentum); the global damping takes
// its share, the forces are cleared, and the sleep timer counts the time a body stays slow.
void RigidWorld::integratePoses(float dt) {
    const float ld = std::max(0.0f, 1.0f - params.linearDamping * dt);
    const float ad = std::max(0.0f, 1.0f - params.angularDamping * dt);
    for (RigidBody& b : bodies_) {
        b.force = Vector3(0.0f);
        b.torque = Vector3(0.0f);
        if (b.invMass == 0) continue;
        b.vel *= ld;
        b.angVel *= ad;
        b.pos += (b.vel + b.biasVel) * dt;
        turnFreely(b, dt); // Euler's equations of the free turn, by an exact splitting (FreeRotation.cpp)
        b.biasVel = Vector3(0.0f);
        b.biasAngVel = Vector3(0.0f);
        b.updateInertia();
        // Sleep timer: time continuously spent below the thresholds.
        bool slow = length(b.vel) < params.sleepLinear && length(b.angVel) < params.sleepAngular;
        b.sleepTimer = slow ? b.sleepTimer + dt : 0.0f;
    }
}

// After the poses moved: continuous collision checks selected sweeps, then joints correct the
// new poses (these later corrections are not covered by CCD). The frozen sleepers are
// released, and the islands decide who falls asleep. The timings of these stages are recorded.
void RigidWorld::finishStep(float dt, std::chrono::steady_clock::time_point tStart) {
    auto tc = std::chrono::steady_clock::now();
    continuousCollision(); // its probe timer is "rigid/ccd ms"
    auto ti = std::chrono::steady_clock::now();
    timings_.ccd = std::chrono::duration<float, std::milli>(ti - tc).count();
    {
        // Joint position stage (nonlinear Gauss-Seidel on the new poses), frozen sleepers still static.
        Probe::Timer t("rigid/solve ms");
        solveJointPositions();
    }
    {
        Probe::Timer t("rigid/islands ms"); // the sleepers let go, then who falls asleep
        unfreezeAll();
        if (params.sleeping) updateIslands(true, dt);
    }
    auto tEnd = std::chrono::steady_clock::now();
    timings_.islands = std::chrono::duration<float, std::milli>(tEnd - ti).count();
    timings_.total = std::chrono::duration<float, std::milli>(tEnd - tStart).count();
}

// What this step was made of, for the probe (the counts of the last substep, the hits summed).
void RigidWorld::reportStep() const {
    Probe::set("rigid/bodies", double(bodies_.size()));
    Probe::set("rigid/bodies awake", double(bodies_.size() - sleepingCount()));
    Probe::set("rigid/contacts", double(contactCount_));
    Probe::set("rigid/manifolds", double(manifolds_.size()));
    Probe::add("rigid/ccd hits", double(ccdHits_));
    Probe::add("rigid/ccd queries", double(ccdDiagnostics_.queries));
    Probe::add("rigid/ccd unresolved", double(ccdDiagnostics_.unresolved));
    Probe::add("rigid/ccd initial contacts", double(ccdDiagnostics_.initialContacts));
    Probe::add("rigid/ccd pass limit", ccdDiagnostics_.passLimitReached ? 1.0 : 0.0);
}


} // namespace rf
