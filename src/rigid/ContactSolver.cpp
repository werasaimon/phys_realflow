// Contacts of RigidWorld: the manifolds of the collision passes and their cache, the graph
// colouring for parallel solving, and the solver of one manifold - the per-point sequential
// impulses with the block LCP of its normal impulses and its friction (Catto 2005-2011).
#include "rigid/RigidWorld.h"

#include "core/Parallel.h"
#include "core/Probe.h"

#include <algorithm>
#include <bitset>
#include <chrono>

namespace rf {

// ---------------------------------------------------------------------------
// Collision detection
// ---------------------------------------------------------------------------
void RigidWorld::addManifold(std::vector<Manifold>& out, int a, int b, ContactManifold& cm) const {
    if (cm.points.empty()) return;
    reduceManifold(cm.points, 4);
    const RigidBody& A = bodies_[a];
    Manifold m;
    m.a = a;
    m.b = b;
    float fb = b >= 0 ? bodies_[b].friction : 0.6f;
    m.friction = std::sqrt(A.friction * fb); // kinetic (sliding) coefficient
    float fsb = b >= 0 ? bodies_[b].staticFriction : 0.8f;
    m.staticFriction = std::max(std::sqrt(A.staticFriction * fsb), m.friction);
    m.restitution = std::max(A.restitution, b >= 0 ? bodies_[b].restitution : 0.1f);
    Matrix3x3 Rt = A.rotation().transposed();
    for (const ContactPoint& c : cm.points) {
        SolverPoint p;
        p.position = c.position;
        p.normal = c.normal;
        p.depth = c.depth;
        p.localA = Rt * (c.position - A.pos);
        p.id = positionHash(p.localA, contactCell(A));
        m.points.push_back(p);
    }
    out.push_back(std::move(m));
}

// Contacts of one body with the static environment: the six domain walls and the static triangle
// mesh (a terrain, a vessel). Static bodies never collide with it.
void RigidWorld::collideStatic(int i, std::vector<Manifold>& out) const {
    if (bodies_[i].invMass == 0) return;
    if (params.collideWithDomain) collideWalls(i, out);
    if (mesh_ && !mesh_->empty()) collideStaticMesh(i, out);
}

// Domain walls: one manifold per plane (static ids -1 .. -6). A sphere touches a wall at one
// point; a box, hull or compound at every vertex within the contact margin of the plane.
void RigidWorld::collideWalls(int i, std::vector<Manifold>& out) const {
    const RigidBody& body = bodies_[i];
    const PosedShape ps = body.posed();
    const Vector3 normals[6] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    const Vector3 points[6] = {domain_.lo, domain_.hi, domain_.lo, domain_.hi, domain_.lo, domain_.hi};
    // The worker's scratch, kept between calls: no allocation per body and wall.
    CollideScratch& S = collideScratch_[size_t(ThreadPool::workerIndex())];
    std::vector<Vector3>& verts = S.verts;
    ContactManifold& wallContacts = S.wall;
    verts.clear();
    if (body.type() == ShapeType::Box) {
        Vector3 h = body.halfExtents();
        for (int k = 0; k < 8; ++k)
            verts.push_back(ps.p + ps.R * Vector3((k & 1) ? h.x : -h.x, (k & 2) ? h.y : -h.y, (k & 4) ? h.z : -h.z));
    } else if (body.type() == ShapeType::ConvexHull) {
        for (const Vector3& v : static_cast<const ConvexHullShape*>(body.shape.get())->vertices()) verts.push_back(ps.p + ps.R * v);
    } else if (body.type() == ShapeType::Compound) {
        for (const auto& c : static_cast<const CompoundShape*>(body.shape.get())->children())
            for (const Vector3& v : c.shape->vertices()) verts.push_back(ps.p + ps.R * (c.R * v + c.t));
    }
    for (int w = 0; w < 6; ++w) {
        const Vector3& n = normals[w];
        // Early out with the support point along -n.
        const float margin = params.contactMargin;
        if (dot(ps.support(-n) - points[w], n) >= margin) continue;
        ContactManifold& cm = wallContacts;
        cm.points.clear();
        if (body.type() == ShapeType::Sphere) {
            float d = dot(body.pos - points[w], n) - body.radius();
            cm.add(body.pos - n * (body.radius() + 0.5f * d), n, -d);
        } else {
            for (const Vector3& v : verts) {
                float d = dot(v - points[w], n);
                if (d < margin) cm.add(v - n * (0.5f * d), n, -d);
            }
        }
        addManifold(out, i, -1 - w, cm);
    }
}

// Static triangle mesh: BVH -> candidate triangles -> narrow phase (static id -7). The query box
// is widened by the contact margin, as for walls and body pairs: speculative contacts. The mesh
// is one-sided: a contact normal pointing against the face normal would pull the body through
// the surface and is dropped.
void RigidWorld::collideStaticMesh(int i, std::vector<Manifold>& out) const {
    const RigidBody& body = bodies_[i];
    const PosedShape ps = body.posed();
    AABB bb = body.worldBounds();
    bb.lo -= Vector3(params.contactMargin);
    bb.hi += Vector3(params.contactMargin);
    if (!mesh_->bounds().overlaps(bb)) return;
    // A terrain sends thousands of candidate triangles per step: the two manifolds are the
    // worker's scratch, cleared, never reallocated.
    CollideScratch& S = collideScratch_[size_t(ThreadPool::workerIndex())];
    ContactManifold &meshContacts = S.mesh, &triangleContacts = S.triangle;
    ContactManifold& cm = meshContacts;
    cm.points.clear();
    mesh_->bvh().queryAABB(bb, [&](uint32_t t) {
        Vector3 a, b, c;
        mesh_->triangle(t, a, b, c);
        TriangleShape tri(a, b, c);
        PosedShape pt{&tri, Matrix3x3(), Vector3(0.0f)};
        ContactManifold& local = triangleContacts;
        local.points.clear();
        if (!narrow_.collide(ps, pt, local)) return;
        const Vector3& fn = mesh_->faceNormal(t);
        for (const ContactPoint& p : local.points)
            if (dot(p.normal, fn) > 0.2f) cm.points.push_back(p); // one-sided: never pull through the surface
    });
    addManifold(out, i, -7, cm);
}

void RigidWorld::collide() {
    manifolds_.clear();
    NarrowPhase::margin = params.contactMargin;
    // 1) Broad phase: candidate pairs from fattened AABBs.
    boxes_.resize(bodies_.size());
    for (size_t i = 0; i < bodies_.size(); ++i) {
        AABB bb = bodies_[i].worldBounds();
        bb.lo -= Vector3(params.contactMargin);
        bb.hi += Vector3(params.contactMargin);
        boxes_[i] = bb;
    }
    auto t0 = std::chrono::steady_clock::now();
    broadphase_->update(boxes_);
    broadphase_->findPairs(pairs_);
    updateWorldTree(); // the bodies' leaves follow them: particles, rays and the viewer ask this tree
    Probe::set("rigid/tree height", worldTree_.height());
    auto t1 = std::chrono::steady_clock::now();
    timings_.broad = std::chrono::duration<float, std::milli>(t1 - t0).count();
    // 2) Narrow phase in parallel (static environment per body, then body pairs). Each task writes
    //    its own slot and the slots are concatenated in a fixed order -> deterministic results.
    //    The slots are kept between steps (cleared, not freed): the same vectors every step.
    const int nb = int(bodies_.size()), np = int(pairs_.size());
    if (collideScratch_.size() != size_t(ThreadPool::instance().threadCount()))
        collideScratch_.resize(size_t(ThreadPool::instance().threadCount()));
    if (slots_.size() < size_t(nb) + np) slots_.resize(size_t(nb) + np);
    for (auto& v : slots_) v.clear();
    // Small grains: pair costs differ by orders of magnitude (sphere-sphere vs compound-compound),
    // and the pool hands out chunks dynamically, so many small chunks balance the load.
    parallelFor(nb, [&](int i) { collideStatic(i, slots_[i]); }, 4);
    parallelFor(np, [&](int k) {
        auto [i, j] = pairs_[k];
        const RigidBody &A = bodies_[i], &B = bodies_[j];
        if (A.invMass == 0 && B.invMass == 0) return;
        ContactManifold& cm = collideScratch_[size_t(ThreadPool::workerIndex())].pair; // the worker's scratch
        cm.points.clear();
        if (narrow_.collide(A.posed(), B.posed(), cm)) addManifold(slots_[size_t(nb) + k], i, j, cm);
    }, 2);
    for (size_t s = 0; s < size_t(nb) + np; ++s)
        for (const Manifold& m : slots_[s]) manifolds_.push_back(m);
    timings_.narrow = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t1).count();
}

// ---------------------------------------------------------------------------
// Sequential impulse solver with warm starting
// ---------------------------------------------------------------------------
// Contact cache: points are identified by a hash of their quantised position in A's frame
// ---------------------------------------------------------------------------
float RigidWorld::contactCell(const RigidBody& A) { return 0.05f * A.boundingRadius() + 0.005f; }

uint64_t RigidWorld::positionHash(const Vector3& localA, float cell) {
    Vector3 q = localA / cell;
    auto k = [](float v) { return uint64_t(uint32_t(int32_t(std::floor(v))) & 0x1FFFFF); };
    return (k(q.x) << 42) | (k(q.y) << 21) | k(q.z);
}

const RigidWorld::SolverPoint* RigidWorld::findCached(const SolverPoints& old, const SolverPoint& p, float cell) {
    // Exact hash hit first (O(1) in practice), then the nearest point when the contact crossed a cell border.
    for (const SolverPoint& o : old)
        if (o.id == p.id && dot(o.normal, p.normal) > 0.95f) return &o;
    const SolverPoint* best = nullptr;
    float bestD = 1.5f * cell;
    for (const SolverPoint& o : old) {
        float d = length(o.localA - p.localA);
        if (d < bestD && dot(o.normal, p.normal) > 0.95f) { bestD = d; best = &o; }
    }
    return best;
}

// ---------------------------------------------------------------------------
// Sequential impulse solver: warm starting, split impulse, Coulomb friction
// ---------------------------------------------------------------------------
static float effMass(const RigidBody& A, const RigidBody* B, const Vector3& ra, const Vector3& rb, const Vector3& dir) {
    float k = A.invMass + dot(cross(A.applyInvInertiaWorld(cross(ra, dir)), ra), dir);
    if (B) k += B->invMass + dot(cross(B->applyInvInertiaWorld(cross(rb, dir)), rb), dir);
    return k > 1e-12f ? 1.0f / k : 0.0f;
}

void RigidWorld::buildColors() {
    // Graph colouring of the contact graph (bodies = vertices, manifolds = edges): manifolds of one
    // colour share no dynamic body and can be solved in parallel. Small scenes keep one sequential
    // batch in bottom-up order, which is the best Gauss-Seidel order for stacks.
    // The batches are kept between steps and only cleared: no allocation once they have grown.
    for (auto& batch : colors_) batch.clear();
    const int nm = int(manifolds_.size());
    if (nm < 256) {
        colors_.resize(1);
        colors_[0].resize(nm);
        for (int i = 0; i < nm; ++i) colors_[0][i] = i;
        parallelColors_ = false;
        return;
    }
    parallelColors_ = true;
    std::vector<uint64_t>& used = colorUsed_;
    used.assign(bodies_.size(), 0);
    colors_.resize(64);
    for (int i = 0; i < nm; ++i) {
        const Manifold& m = manifolds_[i];
        // Static bodies (infinite mass) never conflict.
        uint64_t mask = (bodies_[m.a].invMass > 0 ? used[m.a] : 0) | (m.b >= 0 && bodies_[m.b].invMass > 0 ? used[m.b] : 0);
        int c = 63;
        for (int k = 0; k < 63; ++k)
            if (!(mask & (uint64_t(1) << k))) { c = k; break; }
        colors_[c].push_back(i); // colour 63 is the overflow batch, solved sequentially
        used[m.a] |= uint64_t(1) << c;
        if (m.b >= 0) used[m.b] |= uint64_t(1) << c;
    }
    // Empty trailing batches stay (an empty batch costs nothing to visit, freeing it would cost
    // an allocation next step).
}

// One manifold before the iterations: what every contact point needs (effective mass, bias,
// bounce), the friction patch, the coupling matrix of the points, the rotational lock, and the
// warm start from last step's impulses. Each part is a step below.
void RigidWorld::prepareManifold(Manifold& m, float dt) {
    const RigidBody& A = bodies_[m.a];
    auto oldIt = params.warmStarting ? cache_.find(key(m.a, m.b)) : cache_.end();
    const CachedPair* old = oldIt != cache_.end() ? &oldIt->second : nullptr;
    const float cell = contactCell(A);
    prepareContactPoints(m, dt);
    prepareFrictionPatch(m);
    prepareNormalMassMatrix(m);
    prepareRotationalLock(m, old);
    if (!old) return;
    warmStartManifold(m, *old, cell);
}

// Every point of the manifold as one normal constraint: its effective mass along the normal, the
// target velocity of a speculative gap or the Baumgarte push of a penetration, and whether the
// touch is an impact that should bounce (the bounce itself is applied later by applyRestitution).
// The centre and mean normal of the patch are accumulated on the way.
void RigidWorld::prepareContactPoints(Manifold& m, float dt) {
    RigidBody& A = bodies_[m.a];
    RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
    m.center = Vector3(0.0f);
    m.normal = Vector3(0.0f);
    for (SolverPoint& p : m.points) {
        m.center += p.position;
        m.normal += p.normal;
        Vector3 ra = p.position - A.pos, rb = B ? p.position - B->pos : Vector3(0.0f);
        const Vector3& n = p.normal;
        p.massN = effMass(A, B, ra, rb, n);
        Vector3 dv = A.velocityAt(p.position) - (B ? B->velocityAt(p.position) : Vector3(0.0f));
        float vn = dot(dv, n);
        // The approach speed at the impact: the step integrated gravity into the velocities before
        // the solve, so vn carries one g dt more than the body had when it touched. Bouncing with
        // -e vn would gain 2 g dt / |vn| of energy per bounce (the Noether energy test: +0.74 % at
        // 600 Hz). Only a dynamic body falls; a static support does not.
        const Vector3 gRel = params.gravity * ((A.invMass > 0 ? 1.0f : 0.0f) - (B && B->invMass > 0 ? 1.0f : 0.0f));
        const float vnImpact = vn - dot(gRel, n) * dt;
        // An impact, not a resting touch: faster than the threshold, and the gap (if any) closes
        // within this step. The bounce is applied by applyRestitution() after the solve, never as
        // a target inside the iterations (see there for why).
        const bool impact = vn < -params.restitutionThreshold && (p.depth >= 0 || vn * dt < p.depth);
        p.bounce = impact ? -m.restitution * vnImpact : 0.0f;
        if (p.depth < 0) {
            // Speculative contact: allow closing the gap in this step, but not more. Gaps smaller
            // than the slop count as touching (dead zone): otherwise sub-millimetre differences
            // between corners become cm/s differences of target velocity and tilt landing boxes.
            // Inside the zone the gap closes softly (Baumgarte fraction per step); beyond it the
            // speculative bias allows closing everything but the zone. Continuous at -slop.
            const float slop = params.slop, beta = params.baumgarte;
            p.velocityBias = p.depth > -slop ? beta * p.depth / dt : (p.depth + slop * (1.0f - beta)) / dt;
            p.positionBias = 0;
        } else {
            p.velocityBias = 0;
            p.positionBias = params.baumgarte / dt * std::max(p.depth - params.slop, 0.0f);
            if (!params.splitImpulse) {
                p.velocityBias = std::max(p.velocityBias, p.positionBias);
                p.positionBias = 0;
            }
        }
        p.jn = p.jp = 0;
    }
}

// Friction is solved once per manifold at the centre of the contact patch (ReactPhysics3D): the
// mean normal, two tangents, the patch radius (lever of the twist), the effective masses along
// the tangents and about the normal, and the rolling mass; the friction impulses start at zero
// (the warm start refills them).
void RigidWorld::prepareFrictionPatch(Manifold& m) {
    const RigidBody& A = bodies_[m.a];
    const RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
    const float np = float(m.points.size());
    m.center /= np;
    m.normal = normalize(m.normal);
    const Vector3& n = m.normal;
    m.t1 = anyPerpendicular(n);
    m.t2 = cross(n, m.t1);
    m.patchRadius = 0;
    for (const SolverPoint& p : m.points) m.patchRadius += length(p.position - m.center) / np;
    // Lever of the twist and rolling limits: the size of the smaller body that can move (A is only
    // the lower index of the pair - often a large static platform created first).
    const float sizeA = A.invMass > 0 ? A.boundingRadius() : kInf;
    const float sizeB = B && B->invMass > 0 ? B->boundingRadius() : kInf;
    m.lever = std::min(sizeA, sizeB) < kInf ? std::min(sizeA, sizeB) : A.boundingRadius();

    // --- Friction at the patch centre (tangents, twist) and rolling resistance -------------------
    Vector3 ra = m.center - A.pos, rb = B ? m.center - B->pos : Vector3(0.0f);
    m.massT1 = effMass(A, B, ra, rb, m.t1);
    m.massT2 = effMass(A, B, ra, rb, m.t2);
    float kTwist = dot(n, A.applyInvInertiaWorld(n)) + (B ? dot(n, B->applyInvInertiaWorld(n)) : 0.0f);
    m.massTwist = kTwist > 1e-12f ? 1.0f / kTwist : 0.0f;
    Matrix3x3 Isum = A.invInertiaWorld;
    if (B) for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) Isum.m[i][j] += B->invInertiaWorld.m[i][j];
    m.rollMass = Isum.inverse(); // rolling impulse -> angular velocity; zero for two static bodies
    m.jt1 = m.jt2 = m.jtwist = 0;
    m.jroll = Vector3(0.0f);
    m.jlock = Vector3(0.0f);
}

// Effective-mass matrix of the normal constraints, K_ij = J_i M^-1 J_j^T: how an impulse at point j
// changes the normal velocity at point i through the two bodies. The block solver needs the
// whole matrix to load the (<= 4) points of a face contact together.
void RigidWorld::prepareNormalMassMatrix(Manifold& m) {
    const RigidBody& A = bodies_[m.a];
    const RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
    const int npt = std::min<int>(int(m.points.size()), 4);
    for (int i = 0; i < npt; ++i)
        for (int j = 0; j < npt; ++j) {
            const SolverPoint &pi = m.points[i], &pj = m.points[j];
            float k = 0;
            if (A.invMass > 0) {
                Vector3 ci = cross(pi.position - A.pos, pi.normal), cj = cross(pj.position - A.pos, pj.normal);
                k += A.invMass * dot(pi.normal, pj.normal) + dot(ci, A.applyInvInertiaWorld(cj));
            }
            if (B && B->invMass > 0) {
                Vector3 ci = cross(pi.position - B->pos, pi.normal), cj = cross(pj.position - B->pos, pj.normal);
                k += B->invMass * dot(pi.normal, pj.normal) + dot(ci, B->applyInvInertiaWorld(cj));
            }
            m.K[i][j] = k;
        }
}

// Rotational lock: a face contact (>= 3 points) that is (nearly) at rest relative to the other
// body keeps its relative orientation qB^-1 qA. The error E = q_rel q_ref^-1 lives on SO(3); its
// logarithm (a rotation vector in B's frame) drives an angular constraint like a fixed joint. A
// lock from the last step is kept (old) until a real relative rotation shows up.
void RigidWorld::prepareRotationalLock(Manifold& m, const CachedPair* old) {
    const RigidBody& A = bodies_[m.a];
    const RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
    const Quaternion qB = B ? B->rot : Quaternion();
    const Quaternion qRel = qB.conjugate() * A.rot;
    m.locked = false;
    if (params.rotationalLock && m.points.size() >= 3) {
        Vector3 wRel = A.angVel - (B ? B->angVel : Vector3(0.0f));
        if (old && old->locked) {
            m.locked = true;
            m.lockRef = old->lockRef;
        } else if (length(wRel) < 0.5f) {
            m.locked = true;
            m.lockRef = qRel;
        }
        if (m.locked) {
            Vector3 e = (qRel * m.lockRef.conjugate()).log(); // B frame
            m.lockError = (B ? B->rotation() : Matrix3x3()) * e;
            if (length(m.lockError) > 0.2f) m.locked = false; // really rotating: release the lock
        }
    }
}

// Warm start (Catto): last step's impulses are applied again before the iterations, so a resting
// stack starts every step already carrying its weight. Normal: matched points (same position hash
// / nearest) inherit their impulse. When the contact configuration changed (e.g. 4 corners -> 2
// edge points) the unmatched new points share what is left of last step's total, and the total is
// preserved, so support is neither lost nor doubled when points appear, vanish or jump. Friction
// lives on the manifold: independent of how the individual points moved.
void RigidWorld::warmStartManifold(Manifold& m, const CachedPair& old, float cell) {
    RigidBody& A = bodies_[m.a];
    RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
    const Vector3& n = m.normal;
    float oldTotal = 0;
    for (const SolverPoint& q : old.points) oldTotal += q.jn;
    float matched = 0;
    int unmatched = 0;
    char hit[8] = {0}; // addManifold keeps <= 4 points
    for (int k = 0; k < m.points.size() && k < 8; ++k) {
        SolverPoint& p = m.points[k];
        if (const SolverPoint* q = findCached(old.points, p, cell)) {
            p.jn = q->jn;
            matched += q->jn;
            hit[k] = 1;
        } else {
            ++unmatched;
        }
    }
    if (unmatched > 0 && oldTotal > matched) {
        float share = (oldTotal - matched) / float(unmatched);
        for (int k = 0; k < m.points.size() && k < 8; ++k)
            if (!hit[k]) m.points[k].jn = share;
    }
    float newTotal = 0;
    for (const SolverPoint& p : m.points) newTotal += p.jn;
    if (newTotal > oldTotal && newTotal > 0) {
        float sc = oldTotal / newTotal;
        for (SolverPoint& p : m.points) p.jn *= sc;
    }
    for (SolverPoint& p : m.points) {
        if (p.jn <= 0) continue;
        applyImpulse(m.a, p.normal * p.jn, p.position);
        if (B) applyImpulse(m.b, -p.normal * p.jn, p.position);
    }
    // Friction lives on the manifold: independent of how the individual points moved.
    m.jt1 = dot(old.friction, m.t1);
    m.jt2 = dot(old.friction, m.t2);
    m.jtwist = old.twist;
    m.jroll = old.roll;
    if (m.locked && old.locked) m.jlock = old.lock;
    Vector3 J = m.t1 * m.jt1 + m.t2 * m.jt2;
    applyImpulse(m.a, J, m.center);
    if (B) applyImpulse(m.b, -J, m.center);
    Vector3 L = n * m.jtwist + m.jroll + m.jlock; // pure angular impulses
    // (A static: not written at all - manifolds sharing a static body run in parallel batches.)
    if (A.invMass > 0) A.angVel += A.applyInvInertiaWorld(L);
    if (B && B->invMass > 0) B->angVel -= B->applyInvInertiaWorld(L);
}

// Block solver (Box2D's 2-point block solver generalised to <= 4 points): the normal impulses of
// the manifold are found together as the exact solution of the small LCP by enumerating active
// sets, largest first. A small CFM on the diagonal makes 4 coplanar points (rank-3 K) well posed
// and selects the minimum-energy, i.e. symmetric, load distribution.
// The small LCP of one manifold: K x + b >= 0, x >= 0, x^T (K x + b) = 0 for the (<= 4) normal
// impulses x. An active set is a guess of which points push (x > 0); the guess is right when the
// pushing points end at zero normal velocity and the others are not being pulled into contact.
namespace {
struct ContactLcp {
    int n = 0;
    float bb[4];     // the velocities the points would have with no impulse at all
    float Kc[4][4];  // the effective-mass matrix with the CFM on its diagonal
    float x[4] = {0, 0, 0, 0}; // the solution found by the last successful tryActiveSet

    // Solves the points of `mask` as equalities and checks the complementarity of the rest.
    bool tryActiveSet(int mask) {
        int idx[4], k = 0;
        for (int i = 0; i < n; ++i)
            if (mask & (1 << i)) idx[k++] = i;
        float xs4[4] = {0, 0, 0, 0};
        if (k > 0) {
            float Ms[4][4], rs[4], xs[4];
            for (int i = 0; i < k; ++i) {
                rs[i] = -bb[idx[i]];
                for (int j = 0; j < k; ++j) Ms[i][j] = Kc[idx[i]][idx[j]];
            }
            if (!solveSmall(k, Ms, rs, xs)) return false;
            for (int i = 0; i < k; ++i)
                if (xs[i] < 0.0f) return false;
            for (int i = 0; i < k; ++i) xs4[idx[i]] = xs[i];
        }
        for (int i = 0; i < n; ++i) {
            if (mask & (1 << i)) continue;
            float w = bb[i];
            for (int j = 0; j < n; ++j) w += Kc[i][j] * xs4[j];
            if (w < -1e-5f) return false;
        }
        for (int i = 0; i < n; ++i) x[i] = xs4[i];
        return true;
    }
};
} // namespace

void RigidWorld::blockNormalSolve(Manifold& m) {
    RigidBody& A = bodies_[m.a];
    RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
    ContactLcp lcp;
    const int n = lcp.n = int(m.points.size());
    float a[4], b[4];
    float trace = 0;
    for (int i = 0; i < n; ++i) trace += m.K[i][i];
    const float cfm = params.blockCfm * trace / float(n);
    for (int i = 0; i < n; ++i) {
        const SolverPoint& p = m.points[i];
        a[i] = p.jn;
        Vector3 v = A.velocityAt(p.position) - (B ? B->velocityAt(p.position) : Vector3(0.0f));
        b[i] = dot(v, p.normal) - p.velocityBias;
        for (int j = 0; j < n; ++j) lcp.Kc[i][j] = m.K[i][j] + (i == j ? cfm : 0.0f);
    }
    // Incremental form: v(x) = K (x - a) + v0  ->  K x + (v0 - K a).
    for (int i = 0; i < n; ++i) {
        lcp.bb[i] = b[i];
        for (int j = 0; j < n; ++j) lcp.bb[i] -= m.K[i][j] * a[j];
    }
    // Total enumeration of active sets (exact for <= 4 points). K + cfm I is positive definite, so
    // the LCP solution is unique and the search order only affects speed: first x = 0 (the whole
    // manifold separates - common for speculative points), then the active set of the previous
    // iteration (it rarely changes between iterations), then all sets, largest first.
    const int full = (1 << n) - 1;
    bool found = lcp.tryActiveSet(0);
    int foundMask = 0;
    if (!found && m.activeSet > 0 && m.activeSet <= full && lcp.tryActiveSet(m.activeSet)) found = true, foundMask = m.activeSet;
    for (int size = n; size >= 1 && !found; --size)
        for (int mask = full; mask >= 1 && !found; --mask) {
            if (int(std::bitset<32>(unsigned(mask)).count()) != size || mask == m.activeSet) continue;
            if (lcp.tryActiveSet(mask)) found = true, foundMask = mask;
        }
    if (found) m.activeSet = foundMask;
    // No valid active set (numerical corner case): x = 0. Apply the change of every impulse.
    for (int i = 0; i < n; ++i) {
        SolverPoint& p = m.points[i];
        float d = lcp.x[i] - p.jn;
        p.jn = lcp.x[i];
        if (d == 0) continue;
        applyImpulse(m.a, p.normal * d, p.position);
        if (B) applyImpulse(m.b, -p.normal * d, p.position);
    }
}

// Restitution as one pass after the solve (as Box2D v3 does), from the approach speed remembered
// before it. Inside the iterations a bounce would be an inequality held for the whole step: in a
// stack the lower contact keeps re-pushing the upper box to its bounce speed while the box above
// lands on it, so every level bounces faster than the one below and the stack explodes (100 boxes
// dropped from 1 cm reached 40 m/s). The bounce impulse is not added to jn either: jn is warm
// started, and an impact is a one-off - re-applying it next step as if it were the contact's
// steady load pumped the stack's energy by 15 % per substep.
void RigidWorld::applyRestitution() {
    for (Manifold& m : manifolds_) {
        bool impact = false, loaded = false;
        for (const SolverPoint& p : m.points) {
            impact |= p.bounce > 0;
            loaded |= p.jn > 0;
        }
        if (!impact || !loaded) continue; // a resting touch, or a contact that carried no load
        // The bounce is the normal solve of this manifold once more with the bounce as its target:
        // the block form couples the points, where one sweep point by point gives the first corner
        // the whole impulse and spins a plate hit flat (head-on plates re-hit each other at 30 m/s).
        // The accumulated impulses are put back afterwards: what warm starts the next step is the
        // contact's load, not the bounce.
        float saved[16], bias[16];
        const int np = std::min<int>(int(m.points.size()), 16);
        for (int k = 0; k < np; ++k) {
            SolverPoint& p = m.points[k];
            saved[k] = p.jn;
            bias[k] = p.velocityBias;
            p.velocityBias = std::max(p.velocityBias, p.bounce);
        }
        if (params.blockSolver && np >= 2 && np <= 4) {
            blockNormalSolve(m);
        } else {
            RigidBody& A = bodies_[m.a];
            RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
            for (int k = 0; k < np; ++k) {
                SolverPoint& p = m.points[k];
                const float vn = dot(A.velocityAt(p.position) - (B ? B->velocityAt(p.position) : Vector3(0.0f)), p.normal);
                const float old = p.jn;
                p.jn = std::max(old + p.massN * (p.velocityBias - vn), 0.0f);
                applyImpulse(m.a, p.normal * (p.jn - old), p.position);
                if (B) applyImpulse(m.b, p.normal * (old - p.jn), p.position);
            }
        }
        for (int k = 0; k < np; ++k) {
            m.points[k].jn = saved[k];
            m.points[k].velocityBias = bias[k];
        }
    }
}

// One Gauss-Seidel visit of a manifold: the normal impulses (no penetration), the split impulse
// (penetration recovery on the pseudo velocities), then the friction of the patch limited by the
// normal load just found.
void RigidWorld::solveManifold(Manifold& m) {
    const float total = solveNormalImpulses(m);
    solveSplitImpulse(m);
    solveFriction(m, total);
}

// Velocity of A relative to B at a point of the pair (B static: A's velocity).
Vector3 RigidWorld::relativeVelocity(const Manifold& m, const Vector3& point) const {
    const RigidBody& A = bodies_[m.a];
    return A.velocityAt(point) - (m.b >= 0 ? bodies_[m.b].velocityAt(point) : Vector3(0.0f));
}

// An impulse J on A at a point, and -J on B (Newton's third law).
void RigidWorld::applyPairImpulse(const Manifold& m, const Vector3& J, const Vector3& point) {
    applyImpulse(m.a, J, point);
    if (m.b >= 0) applyImpulse(m.b, -J, point);
}

// A pure angular impulse L on A and -L on B (twist, rolling, the rotational lock).
void RigidWorld::applyPairAngularImpulse(const Manifold& m, const Vector3& L) {
    RigidBody& A = bodies_[m.a];
    RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
    if (A.invMass > 0) A.angVel += A.applyInvInertiaWorld(L);
    if (B && B->invMass > 0) B->angVel -= B->applyInvInertiaWorld(L);
}

// Normal impulses of the manifold as one block: the (<= 4) points are relaxed together until they
// agree, instead of letting whichever corner is solved first take the whole load (which makes
// boxes landing flat start to rock). Points that cannot form a block (one point, or the block
// solver off) are relaxed one after another. Returns the total normal impulse - the load that
// limits the friction.
float RigidWorld::solveNormalImpulses(Manifold& m) {
    if (params.blockSolver && m.points.size() >= 2 && m.points.size() <= 4) {
        blockNormalSolve(m);
    } else {
        const int localIters = m.points.size() > 1 ? params.manifoldIterations : 1; // one point needs one pass
        for (int local = 0; local < localIters; ++local)
            for (SolverPoint& p : m.points) {
                float vn = dot(relativeVelocity(m, p.position), p.normal);
                float old = p.jn;
                p.jn = std::max(old + p.massN * (p.velocityBias - vn), 0.0f);
                applyPairImpulse(m, p.normal * (p.jn - old), p.position);
            }
    }
    float total = 0;
    for (const SolverPoint& p : m.points) total += p.jn;
    return total;
}

// Split impulse (Catto 2006): penetration is pushed out through separate pseudo velocities that
// move the bodies but never enter their real velocities - no energy from the position correction.
// Never warm started.
void RigidWorld::solveSplitImpulse(Manifold& m) {
    const RigidBody& A = bodies_[m.a];
    const RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
    for (SolverPoint& p : m.points) {
        if (p.positionBias <= 0) continue;
        Vector3 bvA = A.biasVel + cross(A.biasAngVel, p.position - A.pos);
        Vector3 bvB = B ? B->biasVel + cross(B->biasAngVel, p.position - B->pos) : Vector3(0.0f);
        float vb = dot(bvA - bvB, p.normal);
        float oldp = p.jp;
        p.jp = std::max(oldp + p.massN * (p.positionBias - vb), 0.0f);
        Vector3 Jp = p.normal * (p.jp - oldp);
        applyBiasImpulse(m.a, Jp, p.position);
        if (B) applyBiasImpulse(m.b, -Jp, p.position);
    }
}

// Friction at the patch centre: Coulomb limit from the total normal load (ReactPhysics3D). Static
// vs kinetic: while the patch sticks (sliding speed below the threshold) the static coefficient
// holds it; once it breaks loose the smaller kinetic coefficient applies. Then the twist about the
// normal, the rotational lock of a resting face contact, or the rolling resistance.
void RigidWorld::solveFriction(Manifold& m, float total) {
    const RigidBody& A = bodies_[m.a];
    const RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
    Vector3 vc = relativeVelocity(m, m.center);
    float slide = length(vc - m.normal * dot(vc, m.normal));
    const float mu = slide < params.stickVelocity ? m.staticFriction : m.friction;
    const float maxF = mu * total;
    float old = m.jt1;
    m.jt1 = clampv(old - m.massT1 * dot(relativeVelocity(m, m.center), m.t1), -maxF, maxF);
    applyPairImpulse(m, m.t1 * (m.jt1 - old), m.center);
    old = m.jt2;
    m.jt2 = clampv(old - m.massT2 * dot(relativeVelocity(m, m.center), m.t2), -maxF, maxF);
    applyPairImpulse(m, m.t2 * (m.jt2 - old), m.center);
    // Twist: relative spin about the normal, limited by the friction moment of the patch.
    Vector3 wRel = A.angVel - (B ? B->angVel : Vector3(0.0f));
    const float maxTwist = maxF * std::max(m.patchRadius, 0.25f * m.lever);
    old = m.jtwist;
    m.jtwist = clampv(old - m.massTwist * dot(wRel, m.normal), -maxTwist, maxTwist);
    applyPairAngularImpulse(m, m.normal * (m.jtwist - old));
    if (m.locked && total > 0) {
        // Angular constraint on SO(3) at velocity level: the relative angular velocity of a resting
        // face contact is driven to zero on all three axes (angular part of a fixed joint). It is
        // breakable - limited by the friction moment the patch can carry - and released in
        // prepareManifold() once log(E) shows a real relative rotation (toppling). No position-level
        // term: the flush orientation is defined by the contact geometry itself.
        const float limit = maxF * std::max(m.patchRadius, 0.25f * m.lever);
        wRel = A.angVel - (B ? B->angVel : Vector3(0.0f));
        Vector3 oldL = m.jlock;
        Vector3 jl = oldL - m.rollMass * wRel;
        float l = length(jl);
        if (l > limit) jl *= limit / l;
        m.jlock = jl;
        applyPairAngularImpulse(m, jl - oldL);
        return;
    }
    // Rolling resistance: damps the remaining relative rotation (tilting / rolling).
    if (params.rollingResistance > 0 && total > 0) {
        wRel = A.angVel - (B ? B->angVel : Vector3(0.0f));
        Vector3 wRoll = wRel - m.normal * dot(wRel, m.normal);
        Vector3 oldR = m.jroll;
        Vector3 jr = oldR - m.rollMass * wRoll;
        float limit = params.rollingResistance * total * m.lever;
        float l = length(jr);
        if (l > limit) jr *= limit / l;
        m.jroll = jr;
        applyPairAngularImpulse(m, jr - oldR);
    }
}

} // namespace rf
