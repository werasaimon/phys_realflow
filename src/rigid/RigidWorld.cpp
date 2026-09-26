#include "rigid/RigidWorld.h"

#include "core/Parallel.h"

#include <algorithm>
#include <chrono>

namespace rf {

// ---------------------------------------------------------------------------
// Body
// ---------------------------------------------------------------------------
AABB RigidBody::worldBounds() const { return shape->boundsAt(rotation(), pos); }

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
    bodies_.push_back(std::move(b));
    return int(bodies_.size()) - 1;
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

void RigidWorld::applyImpulse(int i, const Vector3& J, const Vector3& p) {
    RigidBody& b = bodies_[i];
    if (b.invMass == 0) return;
    b.vel += J * b.invMass;
    b.angVel += b.applyInvInertiaWorld(cross(p - b.pos, J));
}

void RigidWorld::applyBiasImpulse(int i, const Vector3& J, const Vector3& p) {
    RigidBody& b = bodies_[i];
    if (b.invMass == 0) return;
    b.biasVel += J * b.invMass;
    b.biasAngVel += b.applyInvInertiaWorld(cross(p - b.pos, J));
}

std::vector<RigidWorld::DebugContact> RigidWorld::debugContacts() const {
    std::vector<DebugContact> out;
    for (const Manifold& m : manifolds_)
        for (const SolverPoint& p : m.points) out.push_back({m.a, m.b, p.position, p.normal, p.depth, p.jn});
    return out;
}

bool RigidWorld::wakeIsland(int i) {
    const int island = bodies_[i].sleepIsland;
    bool any = false;
    for (RigidBody& b : bodies_) {
        bool member = &b == &bodies_[i] || (island >= 0 && b.sleepIsland == island);
        if (!member) continue;
        any |= b.sleeping;
        b.sleeping = false;
        b.sleepIsland = -1;
    }
    return any;
}

void RigidWorld::wake(int i) {
    wakeIsland(i);
    bodies_[i].sleepTimer = 0;
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

size_t RigidWorld::sleepingCount() const {
    size_t n = 0;
    for (const RigidBody& b : bodies_) n += b.sleeping;
    return n;
}

void RigidWorld::clear() {
    bodies_.clear();
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
    contactCount_ = 0;
    // XPBD keeps its contacts between collision passes: they belong to the old bodies.
    xcontacts_.clear();
    substepCounter_ = 0;
}

RigidWorld::Frozen RigidWorld::makeStatic(int i) {
    RigidBody& b = bodies_[i];
    const Frozen f{i, b.invMass, b.invInertiaLocal};
    b.invMass = 0;
    b.invInertiaLocal = Vector3(0.0f);
    b.updateInertia();
    return f;
}

void RigidWorld::restore(const Frozen& f) {
    RigidBody& b = bodies_[f.body];
    b.invMass = f.invMass;
    b.invInertiaLocal = f.invInertiaLocal;
    b.updateInertia();
}

void RigidWorld::hold(int i) {
    if (i < 0 || i >= int(bodies_.size()) || bodies_[i].invMass == 0) return; // static already
    held_.push_back(makeStatic(i));
}

void RigidWorld::releaseHeld() {
    for (const Frozen& f : held_) {
        if (f.body >= int(bodies_.size())) continue;
        restore(f);
        wake(f.body);
    }
    held_.clear();
}

// Sleeping bodies take part in the step as static bodies: their inverse mass is set to zero for
// the duration of the step (restored in unfreezeAll), so every solver path treats them as fixed.
void RigidWorld::freezeSleepers() {
    frozen_.clear();
    for (int i = 0; i < int(bodies_.size()); ++i) {
        RigidBody& b = bodies_[i];
        if (!b.sleeping || b.invMass == 0) continue;
        b.vel = b.angVel = Vector3(0.0f);
        frozen_.push_back(makeStatic(i));
    }
}

void RigidWorld::unfreezeAll(bool onlyAwake) {
    std::vector<Frozen> keep;
    for (const Frozen& f : frozen_) {
        if (onlyAwake && bodies_[f.body].sleeping) { keep.push_back(f); continue; }
        restore(f);
    }
    frozen_.swap(keep);
}

// Islands = connected components of the contact graph among dynamic bodies (union-find).
// Before the solve: an island containing any moving body wakes up entirely (and its frozen
// members get their mass back). After the solve: an island whose bodies have all been slow for
// sleepTime falls asleep.
bool RigidWorld::updateIslands(bool decideSleep, float dt) {
    const int n = int(bodies_.size());
    islandParent_.resize(n);
    for (int i = 0; i < n; ++i) islandParent_[i] = i;
    auto find = [&](int x) {
        while (islandParent_[x] != x) x = islandParent_[x] = islandParent_[islandParent_[x]];
        return x;
    };
    for (const Manifold& m : manifolds_) {
        if (m.b < 0 || bodies_[m.a].mass <= 0 || bodies_[m.b].mass <= 0) continue;
        bool touching = false;
        for (const SolverPoint& p : m.points) touching |= p.depth > -params.slop;
        if (touching) islandParent_[find(m.a)] = find(m.b);
    }
    for (const auto& j : joints_)
        if (j->b >= 0 && bodies_[j->a].mass > 0 && bodies_[j->b].mass > 0) islandParent_[find(j->a)] = find(j->b);
    std::vector<float> minTimer(n, kInf);
    for (int i = 0; i < n; ++i)
        if (bodies_[i].mass > 0) minTimer[find(i)] = std::min(minTimer[find(i)], bodies_[i].sleepTimer);
    bool woke = false;
    std::vector<int> newIsland(n, -1);
    for (int i = 0; i < n; ++i) {
        RigidBody& b = bodies_[i];
        if (b.mass <= 0) continue;
        bool islandRests = minTimer[find(i)] >= params.sleepTime;
        if (!decideSleep) {
            // Touched by a moving body: the whole sleeping island wakes (Box2D).
            if (b.sleeping && !islandRests) woke |= wakeIsland(i);
        } else if (islandRests && !b.sleeping) {
            int r = find(i);
            if (newIsland[r] < 0) newIsland[r] = nextIsland_++;
            b.sleeping = true;
            b.sleepIsland = newIsland[r];
            b.vel = b.angVel = Vector3(0.0f);
        }
    }
    if (woke) unfreezeAll(true); // woken bodies get their mass back for this very step
    return woke;
    (void)dt;
}

bool RigidWorld::raycast(const Vector3& o, const Vector3& d, float maxT, int& body, float& t, Vector3& normal) const {
    body = -1;
    t = maxT;
    Vector3 inv(1.0f / (std::fabs(d.x) > 1e-12f ? d.x : 1e-12f), 1.0f / (std::fabs(d.y) > 1e-12f ? d.y : 1e-12f),
             1.0f / (std::fabs(d.z) > 1e-12f ? d.z : 1e-12f));
    for (int i = 0; i < int(bodies_.size()); ++i) {
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

void RigidWorld::grab(int i, const Vector3& p) {
    if (i < 0 || i >= int(bodies_.size()) || bodies_[i].invMass == 0) return;
    grab_ = GrabJoint();
    grab_.active = true;
    grab_.body = i;
    grab_.localAnchor = bodies_[i].rotation().transposed() * (p - bodies_[i].pos);
    grab_.target = p;
    wake(i);
}

// Soft constraint (Box2D mouse joint): spring-damper with the given frequency and damping ratio,
// expressed through gamma (softness) and beta (position feedback), solved as a 3D point constraint.
void RigidWorld::prepareGrab(float dt) {
    if (!grab_.active) return;
    if (grab_.body >= int(bodies_.size())) { releaseGrab(); return; }
    RigidBody& b = bodies_[grab_.body];
    b.sleepTimer = 0; // a held body never falls asleep
    if (b.invMass == 0) return;
    const float omega = 2.0f * kPi * grab_.frequency;
    const float d = 2.0f * b.mass * grab_.damping * omega;
    const float k = b.mass * omega * omega;
    grabGamma_ = dt * (d + dt * k);
    grabGamma_ = grabGamma_ > 0 ? 1.0f / grabGamma_ : 0.0f;
    const float beta = dt * k * grabGamma_;
    Vector3 r = b.rotation() * grab_.localAnchor;
    Matrix3x3 S = Matrix3x3::skew(r);
    Matrix3x3 K = S * b.invInertiaWorld * S.transposed();
    for (int i = 0; i < 3; ++i) K.m[i][i] += b.invMass + grabGamma_;
    grabMass_ = K.inverse();
    grabBias_ = (b.pos + r - grab_.target) * beta;
    // Warm start and extra angular damping so the held body does not spin forever.
    b.vel += grab_.impulse * b.invMass;
    b.angVel += b.applyInvInertiaWorld(cross(r, grab_.impulse));
    b.angVel *= std::max(0.0f, 1.0f - 2.0f * dt);
}

void RigidWorld::solveGrab(float dt) {
    if (!grab_.active) return;
    RigidBody& b = bodies_[grab_.body];
    if (b.invMass == 0) return;
    Vector3 r = b.rotation() * grab_.localAnchor;
    Vector3 cdot = b.vel + cross(b.angVel, r);
    Vector3 imp = grabMass_ * (-(cdot + grabBias_ + grab_.impulse * grabGamma_));
    Vector3 old = grab_.impulse;
    grab_.impulse += imp;
    float maxImpulse = grab_.maxForce * b.mass * 9.81f * dt; // weights at standard gravity, whatever params.gravity is
    float l = length(grab_.impulse);
    if (l > maxImpulse) grab_.impulse *= maxImpulse / l;
    imp = grab_.impulse - old;
    b.vel += imp * b.invMass;
    b.angVel += b.applyInvInertiaWorld(cross(r, imp));
}

// ---------------------------------------------------------------------------
// Joints
// ---------------------------------------------------------------------------
template <class J>
J& RigidWorld::attach(std::unique_ptr<J> j, const Vector3& anchorA, const Vector3& anchorB, const Vector3& axisW) {
    static const RigidBody world = [] { RigidBody w; w.mass = 0; w.invMass = 0; return w; }();
    const RigidBody& A = bodies_[j->a];
    const RigidBody& B = j->b >= 0 ? bodies_[j->b] : world;
    Vector3 axis = normalize(axisW);
    j->localAnchorA = A.rot.conjugate().rotate(anchorA - A.pos);
    j->localAnchorB = B.rot.conjugate().rotate(anchorB - B.pos);
    j->localAxisA = A.rot.conjugate().rotate(axis);
    j->localAxisB = B.rot.conjugate().rotate(axis);
    j->refRel = B.rot.conjugate() * A.rot;
    if (auto* h = dynamic_cast<HingeJoint*>(j.get())) {
        Vector3 ref = anyPerpendicular(axis);
        h->localRefA = A.rot.conjugate().rotate(ref);
        h->localRefB = B.rot.conjugate().rotate(ref);
    }
    wake(j->a);
    if (j->b >= 0) wake(j->b);
    J& ref = *j;
    joints_.push_back(std::move(j));
    return ref;
}

BallJoint& RigidWorld::addBallJoint(int a, int b, const Vector3& p) {
    return attach(std::make_unique<BallJoint>(a, b), p, p, Vector3(1, 0, 0));
}

HingeJoint& RigidWorld::addHingeJoint(int a, int b, const Vector3& p, const Vector3& axis) {
    return attach(std::make_unique<HingeJoint>(a, b), p, p, axis);
}

SliderJoint& RigidWorld::addSliderJoint(int a, int b, const Vector3& axis) {
    Vector3 p = bodies_[a].pos; // anchor at A's centre
    return attach(std::make_unique<SliderJoint>(a, b), p, p, axis);
}

FixedJoint& RigidWorld::addFixedJoint(int a, int b) {
    Vector3 p = b >= 0 ? (bodies_[a].pos + bodies_[b].pos) * 0.5f : bodies_[a].pos;
    return attach(std::make_unique<FixedJoint>(a, b), p, p, Vector3(1, 0, 0));
}

DistanceJoint& RigidWorld::addDistanceJoint(int a, int b, const Vector3& pa, const Vector3& pb) {
    DistanceJoint& j = attach(std::make_unique<DistanceJoint>(a, b), pa, pb, Vector3(1, 0, 0));
    j.length = length(pa - pb);
    return j;
}

void RigidWorld::solveJointPositions() {
    for (int it = 0; it < params.jointPositionIterations; ++it)
        for (auto& j : joints_) j->solvePosition(bodies_);
}

// ---------------------------------------------------------------------------
// Continuous collision detection
// ---------------------------------------------------------------------------
void RigidWorld::continuousCollision() {
    ccdHits_ = 0;
    ccdClamped_.assign(bodies_.size(), 0);
    if (!params.ccd) return;
    const float tol = params.ccdTolerance;
    auto sweptOf = [&](const RigidBody& b) {
        SweptPose sp;
        sp.shape = b.shape.get();
        sp.p0 = b.prevPos;
        sp.p1 = b.pos;
        sp.q0 = b.prevRot;
        sp.dTheta = (b.rot * b.prevRot.conjugate()).log();
        return sp;
    };
    auto boundsAt = [](const SweptPose& sp, float s) {
        PosedShape ps = sp.at(s);
        return ps.shape->boundsAt(ps.R, ps.p);
    };
    auto isFast = [&](const RigidBody& b, const SweptPose& sp) {
        Vector3 ext = b.shape->localBounds().extent();
        float minExtent = std::max(0.5f * minComp(ext), 1e-3f);
        return length(sp.p1 - sp.p0) + sp.angularReach() >= params.ccdThreshold * minExtent;
    };
    // A compound body is swept part by part: the hull of its union is not its surface, a fast body
    // can fly into a concave region of that hull (a teapot's handle) without touching any part.
    auto partsOf = [](const SweptPose& sp, std::vector<SweptPose>& out) {
        out.clear();
        if (sp.shape->type() != ShapeType::Compound) {
            out.push_back(sp);
            return;
        }
        for (const CompoundShape::Child& c : static_cast<const CompoundShape*>(sp.shape)->children()) {
            SweptPose part = sp;
            part.shape = c.shape.get();
            part.partR = c.R;
            part.partT = c.t;
            out.push_back(part);
        }
    };
    // Earliest time of impact of any part of A with any part of B.
    auto partsTimeOfImpact = [&](const SweptPose& A, const SweptPose& B) {
        std::vector<SweptPose> pa, pb;
        partsOf(A, pa);
        partsOf(B, pb);
        ToiResult first;
        for (const SweptPose& a : pa)
            for (const SweptPose& b : pb) {
                ToiResult r = timeOfImpact(a, b, tol);
                if (r.hit && r.s < first.s) first = r;
            }
        return first;
    };
    const int n = int(bodies_.size());
    {
        // Nothing fast this step (the usual case): no sweeps to test.
        bool anyFast = false;
        for (int i = 0; i < n && !anyFast; ++i)
            anyFast = bodies_[i].invMass > 0 && !bodies_[i].sleeping && isFast(bodies_[i], sweptOf(bodies_[i]));
        if (!anyFast) return;
    }

    // Passes: find the times of impact of all fast bodies with the motions of this pass; both
    // bodies of a colliding pair stop at their common time of impact (motion clamping); a stopped
    // body stays put for the rest of the step (its sweep becomes static), and the pass repeats so
    // bodies that would now run into it are caught too (A hits B, B hits C, ...).
    for (int pass = 0; pass < 8; ++pass) {
        std::vector<SweptPose> sw(n);
        std::vector<AABB> swBox(n);
        for (int i = 0; i < n; ++i) {
            sw[i] = sweptOf(bodies_[i]);
            swBox[i] = boundsAt(sw[i], 0);
            swBox[i].expand(boundsAt(sw[i], 1));
            swBox[i].lo -= Vector3(tol);
            swBox[i].hi += Vector3(tol);
        }
        std::vector<float> sMin(n, 1.0f);
        // Broad phase for CCD: BVH over the swept bounds (O(n log n) instead of all pairs).
        BVH sweptTree;
        sweptTree.build(swBox, 2);
        for (int i = 0; i < n; ++i) {
            const RigidBody& body = bodies_[i];
            if (body.invMass == 0 || !isFast(body, sw[i])) continue;
            sweptTree.queryAABB(swBox[i], [&](uint32_t ju) {
                int j = int(ju);
                if (j == i || !swBox[j].overlaps(swBox[i])) return;
                ToiResult r = partsTimeOfImpact(sw[i], sw[j]);
                if (!r.hit) return;
                sMin[i] = std::min(sMin[i], r.s);
                if (bodies_[j].invMass > 0) sMin[j] = std::min(sMin[j], r.s); // the pair stops together
            });
            if (params.collideWithDomain) {
                const Vector3 normals[6] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
                const Vector3 points[6] = {domain_.lo, domain_.hi, domain_.lo, domain_.hi, domain_.lo, domain_.hi};
                for (int w = 0; w < 6; ++w) {
                    ToiResult r = timeOfImpactPlane(sw[i], normals[w], dot(normals[w], points[w]), tol);
                    if (r.hit) sMin[i] = std::min(sMin[i], r.s);
                }
            }
            if (mesh_ && !mesh_->empty() && mesh_->bounds().overlaps(swBox[i])) {
                mesh_->bvh().queryAABB(swBox[i], [&](uint32_t t) {
                    Vector3 a, b, c;
                    mesh_->triangle(t, a, b, c);
                    TriangleShape tri(a, b, c);
                    SweptPose T;
                    T.shape = &tri;
                    ToiResult r = partsTimeOfImpact(sw[i], T);
                    if (r.hit) sMin[i] = std::min(sMin[i], r.s);
                });
            }
        }
        bool any = false;
        for (int i = 0; i < n; ++i) {
            if (sMin[i] >= 1.0f) continue;
            RigidBody& b = bodies_[i];
            // Stop at the time of impact, keep the velocity: the speculative contact of the next
            // step resolves the impact (with restitution, see prepareManifold).
            b.pos = sw[i].p0 + (sw[i].p1 - sw[i].p0) * sMin[i];
            b.rot = sw[i].q0.integrated(sw[i].dTheta, sMin[i]);
            b.prevPos = b.pos; // static for the remaining passes of this step
            b.prevRot = b.rot;
            b.updateInertia();
            if (!ccdClamped_[i]) ++ccdHits_;
            ccdClamped_[i] = 1;
            any = true;
        }
        if (!any) break;
    }
}

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

void RigidWorld::collideStatic(int i, std::vector<Manifold>& out) const {
    const RigidBody& body = bodies_[i];
    if (body.invMass == 0) return;
    const PosedShape ps = body.posed();

    // Domain walls: one manifold per plane (static ids -1 .. -6).
    if (params.collideWithDomain) {
        const Vector3 normals[6] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
        const Vector3 points[6] = {domain_.lo, domain_.hi, domain_.lo, domain_.hi, domain_.lo, domain_.hi};
        std::vector<Vector3> verts;
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
            ContactManifold cm;
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
    // is widened by the contact margin, as for walls and body pairs: speculative contacts.
    if (mesh_ && !mesh_->empty()) {
        AABB bb = body.worldBounds();
        bb.lo -= Vector3(params.contactMargin);
        bb.hi += Vector3(params.contactMargin);
        if (!mesh_->bounds().overlaps(bb)) return;
        ContactManifold cm;
        mesh_->bvh().queryAABB(bb, [&](uint32_t t) {
            Vector3 a, b, c;
            mesh_->triangle(t, a, b, c);
            TriangleShape tri(a, b, c);
            PosedShape pt{&tri, Matrix3x3(), Vector3(0.0f)};
            ContactManifold local;
            if (!narrow_.collide(ps, pt, local)) return;
            const Vector3& fn = mesh_->faceNormal(t);
            for (const ContactPoint& p : local.points)
                if (dot(p.normal, fn) > 0.2f) cm.points.push_back(p); // one-sided: never pull through the surface
        });
        addManifold(out, i, -7, cm);
    }
}

void RigidWorld::collide() {
    manifolds_.clear();
    NarrowPhase::margin = params.contactMargin;
    // 1) Broad phase: candidate pairs from fattened AABBs.
    std::vector<AABB> boxes(bodies_.size());
    for (size_t i = 0; i < bodies_.size(); ++i) {
        AABB bb = bodies_[i].worldBounds();
        bb.lo -= Vector3(params.contactMargin);
        bb.hi += Vector3(params.contactMargin);
        boxes[i] = bb;
    }
    auto t0 = std::chrono::steady_clock::now();
    broadphase_->update(boxes);
    broadphase_->findPairs(pairs_);
    auto t1 = std::chrono::steady_clock::now();
    timings_.broad = std::chrono::duration<float, std::milli>(t1 - t0).count();
    // 2) Narrow phase in parallel (static environment per body, then body pairs). Each task writes
    //    its own slot and the slots are concatenated in a fixed order -> deterministic results.
    const int nb = int(bodies_.size()), np = int(pairs_.size());
    std::vector<std::vector<Manifold>> slots(size_t(nb) + np);
    // Small grains: pair costs differ by orders of magnitude (sphere-sphere vs compound-compound),
    // and the pool hands out chunks dynamically, so many small chunks balance the load.
    parallelFor(nb, [&](int i) { collideStatic(i, slots[i]); }, 4);
    parallelFor(np, [&](int k) {
        auto [i, j] = pairs_[k];
        const RigidBody &A = bodies_[i], &B = bodies_[j];
        if (A.invMass == 0 && B.invMass == 0) return;
        ContactManifold cm;
        if (narrow_.collide(A.posed(), B.posed(), cm)) addManifold(slots[size_t(nb) + k], i, j, cm);
    }, 2);
    for (auto& v : slots)
        for (auto& m : v) manifolds_.push_back(std::move(m));
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

const RigidWorld::SolverPoint* RigidWorld::findCached(const std::vector<SolverPoint>& old, const SolverPoint& p, float cell) {
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
    colors_.clear();
    const int nm = int(manifolds_.size());
    if (nm < 256) {
        colors_.emplace_back(nm);
        for (int i = 0; i < nm; ++i) colors_[0][i] = i;
        parallelColors_ = false;
        return;
    }
    parallelColors_ = true;
    std::vector<uint64_t> used(bodies_.size(), 0);
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
    while (!colors_.empty() && colors_.back().empty() && colors_.size() > 1) colors_.pop_back();
}

template <class F> void RigidWorld::forEachManifold(F&& f) {
    for (size_t c = 0; c < colors_.size(); ++c) {
        const auto& batch = colors_[c];
        if (parallelColors_ && c < 63 && batch.size() >= 64)
            parallelFor(int(batch.size()), [&](int k) { f(manifolds_[batch[k]]); }, 16);
        else
            for (int i : batch) f(manifolds_[i]);
    }
}

void RigidWorld::prepareManifold(Manifold& m, float dt) {
    RigidBody& A = bodies_[m.a];
    RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
    auto oldIt = params.warmStarting ? cache_.find(key(m.a, m.b)) : cache_.end();
    const CachedPair* old = oldIt != cache_.end() ? &oldIt->second : nullptr;
    const float cell = contactCell(A);

    // --- Normal constraints per point ------------------------------------------------------
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
        if (p.depth < 0) {
            // Speculative contact: allow closing the gap in this step, but not more. Gaps smaller
            // than the slop count as touching (dead zone): otherwise sub-millimetre differences
            // between corners become cm/s differences of target velocity and tilt landing boxes.
            // Inside the zone the gap closes softly (Baumgarte fraction per step); beyond it the
            // speculative bias allows closing everything but the zone. Continuous at -slop.
            const float slop = params.slop, beta = params.baumgarte;
            p.velocityBias = p.depth > -slop ? beta * p.depth / dt : (p.depth + slop * (1.0f - beta)) / dt;
            p.positionBias = 0;
            // The gap closes within this step at this approach speed: it is an impact, not a resting
            // touch - apply restitution to the approach velocity (after CCD clamping every fast hit
            // arrives here as a speculative contact).
            bool clamped = (m.a < int(ccdClamped_.size()) && ccdClamped_[m.a]) ||
                           (m.b >= 0 && m.b < int(ccdClamped_.size()) && ccdClamped_[m.b]);
            if (clamped && vn < -1.0f && vn * dt < p.depth) p.velocityBias = std::max(p.velocityBias, -m.restitution * vn);
        } else {
            p.velocityBias = vn < -1.0f ? -m.restitution * vn : 0.0f;
            p.positionBias = params.baumgarte / dt * std::max(p.depth - params.slop, 0.0f);
            if (!params.splitImpulse) {
                p.velocityBias = std::max(p.velocityBias, p.positionBias);
                p.positionBias = 0;
            }
        }
        p.jn = p.jp = 0;
    }
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

    // Effective-mass matrix of the normal constraints: K_ij = J_i M^-1 J_j^T.
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

    // Rotational lock: a face contact (>= 3 points) that is (nearly) at rest relative to the other
    // body keeps its relative orientation qB^-1 qA. The error E = q_rel q_ref^-1 lives on SO(3); its
    // logarithm (a rotation vector in B's frame) drives an angular constraint like a fixed joint.
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
    if (!old) return;

    // --- Warm start -------------------------------------------------------------------------------
    // Normal: matched points (same position hash / nearest) inherit their impulse. When the
    // contact configuration changed (e.g. 4 corners -> 2 edge points) the unmatched new points
    // share what is left of last step's total, and the total is preserved, so support is neither
    // lost nor doubled when points appear, vanish or jump.
    float oldTotal = 0;
    for (const SolverPoint& q : old->points) oldTotal += q.jn;
    float matched = 0;
    int unmatched = 0;
    char hit[8] = {0}; // addManifold keeps <= 4 points
    for (size_t k = 0; k < m.points.size() && k < 8; ++k) {
        SolverPoint& p = m.points[k];
        if (const SolverPoint* q = findCached(old->points, p, cell)) {
            p.jn = q->jn;
            matched += q->jn;
            hit[k] = 1;
        } else {
            ++unmatched;
        }
    }
    if (unmatched > 0 && oldTotal > matched) {
        float share = (oldTotal - matched) / float(unmatched);
        for (size_t k = 0; k < m.points.size() && k < 8; ++k)
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
    m.jt1 = dot(old->friction, m.t1);
    m.jt2 = dot(old->friction, m.t2);
    m.jtwist = old->twist;
    m.jroll = old->roll;
    if (m.locked && old->locked) m.jlock = old->lock;
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
void RigidWorld::blockNormalSolve(Manifold& m) {
    RigidBody& A = bodies_[m.a];
    RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
    const int n = int(m.points.size());
    float a[4], b[4], Kc[4][4];
    float trace = 0;
    for (int i = 0; i < n; ++i) trace += m.K[i][i];
    const float cfm = params.blockCfm * trace / float(n);
    for (int i = 0; i < n; ++i) {
        const SolverPoint& p = m.points[i];
        a[i] = p.jn;
        Vector3 v = A.velocityAt(p.position) - (B ? B->velocityAt(p.position) : Vector3(0.0f));
        b[i] = dot(v, p.normal) - p.velocityBias;
        for (int j = 0; j < n; ++j) Kc[i][j] = m.K[i][j] + (i == j ? cfm : 0.0f);
    }
    // Incremental form: v(x) = K (x - a) + v0  ->  K x + (v0 - K a).
    float bb[4];
    for (int i = 0; i < n; ++i) {
        bb[i] = b[i];
        for (int j = 0; j < n; ++j) bb[i] -= m.K[i][j] * a[j];
    }
    // Total enumeration of active sets (exact for <= 4 points). K + cfm I is positive definite, so
    // the LCP solution is unique and the search order only affects speed: first x = 0 (the whole
    // manifold separates - common for speculative points), then the active set of the previous
    // iteration (it rarely changes between iterations), then all sets, largest first.
    float xBest[4] = {0, 0, 0, 0};
    auto tryMask = [&](int mask) {
        int idx[4], k = 0;
        for (int i = 0; i < n; ++i)
            if (mask & (1 << i)) idx[k++] = i;
        float x[4] = {0, 0, 0, 0};
        if (k > 0) {
            float Ms[4][4], rs[4], xs[4];
            for (int i = 0; i < k; ++i) {
                rs[i] = -bb[idx[i]];
                for (int j = 0; j < k; ++j) Ms[i][j] = Kc[idx[i]][idx[j]];
            }
            if (!solveSmall(k, Ms, rs, xs)) return false;
            for (int i = 0; i < k; ++i)
                if (xs[i] < 0.0f) return false;
            for (int i = 0; i < k; ++i) x[idx[i]] = xs[i];
        }
        for (int i = 0; i < n; ++i) {
            if (mask & (1 << i)) continue;
            float w = bb[i];
            for (int j = 0; j < n; ++j) w += Kc[i][j] * x[j];
            if (w < -1e-5f) return false;
        }
        for (int i = 0; i < n; ++i) xBest[i] = x[i];
        return true;
    };
    const int full = (1 << n) - 1;
    bool found = tryMask(0);
    int foundMask = 0;
    if (!found && m.activeSet > 0 && m.activeSet <= full && tryMask(m.activeSet)) found = true, foundMask = m.activeSet;
    for (int size = n; size >= 1 && !found; --size)
        for (int mask = full; mask >= 1 && !found; --mask) {
            if (__builtin_popcount(unsigned(mask)) != size || mask == m.activeSet) continue;
            if (tryMask(mask)) found = true, foundMask = mask;
        }
    if (found) m.activeSet = foundMask;
    // No valid active set (numerical corner case): x = 0.
    for (int i = 0; i < n; ++i) {
        SolverPoint& p = m.points[i];
        float d = xBest[i] - p.jn;
        p.jn = xBest[i];
        if (d == 0) continue;
        applyImpulse(m.a, p.normal * d, p.position);
        if (B) applyImpulse(m.b, -p.normal * d, p.position);
    }
}

void RigidWorld::solveManifold(Manifold& m) {
    RigidBody& A = bodies_[m.a];
    RigidBody* B = m.b >= 0 ? &bodies_[m.b] : nullptr;
    auto relVel = [&](const Vector3& p) { return A.velocityAt(p) - (B ? B->velocityAt(p) : Vector3(0.0f)); };
    auto apply = [&](const Vector3& J, const Vector3& p) {
        applyImpulse(m.a, J, p);
        if (B) applyImpulse(m.b, -J, p);
    };
    auto applyAngular = [&](const Vector3& L) {
        if (A.invMass > 0) A.angVel += A.applyInvInertiaWorld(L);
        if (B && B->invMass > 0) B->angVel -= B->applyInvInertiaWorld(L);
    };
    // Normal impulses of the manifold as one block: the (<= 4) points are relaxed together until
    // they agree, instead of letting whichever corner is solved first take the whole load (which
    // makes boxes landing flat start to rock).
    float total = 0;
    if (params.blockSolver && m.points.size() >= 2 && m.points.size() <= 4) {
        blockNormalSolve(m);
    } else {
        const int localIters = m.points.size() > 1 ? params.manifoldIterations : 1; // one point needs one pass
        for (int local = 0; local < localIters; ++local)
            for (SolverPoint& p : m.points) {
                float vn = dot(relVel(p.position), p.normal);
                float old = p.jn;
                p.jn = std::max(old + p.massN * (p.velocityBias - vn), 0.0f);
                apply(p.normal * (p.jn - old), p.position);
            }
    }
    for (const SolverPoint& p : m.points) total += p.jn;
    // Split impulse: penetration recovery on the pseudo velocities only (never warm started).
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
    // Friction at the patch centre: Coulomb limit from the total normal load (ReactPhysics3D).
    // Static vs kinetic: while the patch sticks (sliding speed below the threshold) the static
    // coefficient holds it; once it breaks loose the smaller kinetic coefficient applies.
    Vector3 vc = relVel(m.center);
    float slide = length(vc - m.normal * dot(vc, m.normal));
    const float mu = slide < params.stickVelocity ? m.staticFriction : m.friction;
    const float maxF = mu * total;
    float old = m.jt1;
    m.jt1 = clampv(old - m.massT1 * dot(relVel(m.center), m.t1), -maxF, maxF);
    apply(m.t1 * (m.jt1 - old), m.center);
    old = m.jt2;
    m.jt2 = clampv(old - m.massT2 * dot(relVel(m.center), m.t2), -maxF, maxF);
    apply(m.t2 * (m.jt2 - old), m.center);
    // Twist: relative spin about the normal, limited by the friction moment of the patch.
    Vector3 wRel = A.angVel - (B ? B->angVel : Vector3(0.0f));
    const float maxTwist = maxF * std::max(m.patchRadius, 0.25f * m.lever);
    old = m.jtwist;
    m.jtwist = clampv(old - m.massTwist * dot(wRel, m.normal), -maxTwist, maxTwist);
    applyAngular(m.normal * (m.jtwist - old));
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
        applyAngular(jl - oldL);
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
        applyAngular(jr - oldR);
    }
}

void RigidWorld::computeLevels() {
    // BFS over the contact graph starting from bodies touching the static environment (level 0).
    const int n = int(bodies_.size());
    const int kFar = 1 << 28;
    levels_.assign(n, kFar);
    // Adjacency in compressed rows (one allocation instead of one vector per body).
    std::vector<int> start(n + 1, 0), adj, queue;
    auto dynamicPair = [&](const Manifold& m) {
        return m.b >= 0 && bodies_[m.a].invMass > 0 && bodies_[m.b].invMass > 0;
    };
    for (const Manifold& m : manifolds_) {
        bool bStatic = m.b < 0 || bodies_[m.b].invMass == 0;
        bool aStatic = bodies_[m.a].invMass == 0;
        if (bStatic && !aStatic && levels_[m.a] != 0) { levels_[m.a] = 0; queue.push_back(m.a); }
        if (aStatic && m.b >= 0 && !bStatic && levels_[m.b] != 0) { levels_[m.b] = 0; queue.push_back(m.b); }
        if (dynamicPair(m)) ++start[m.a + 1], ++start[m.b + 1];
    }
    for (int i = 0; i < n; ++i) start[i + 1] += start[i];
    adj.resize(start[n]);
    std::vector<int> fill(start.begin(), start.end() - 1);
    for (const Manifold& m : manifolds_)
        if (dynamicPair(m)) adj[fill[m.a]++] = m.b, adj[fill[m.b]++] = m.a;
    for (size_t h = 0; h < queue.size(); ++h) {
        int u = queue[h];
        for (int e = start[u]; e < start[u + 1]; ++e) {
            int v = adj[e];
            if (levels_[v] > levels_[u] + 1) { levels_[v] = levels_[u] + 1; queue.push_back(v); }
        }
    }
}

void RigidWorld::solveManifoldShock(Manifold& m) {
    // Level of each side: static environment counts as -1 (always "below").
    const int la = bodies_[m.a].invMass == 0 ? -1 : levels_[m.a];
    const int lb = (m.b < 0 || bodies_[m.b].invMass == 0) ? -1 : levels_[m.b];
    if (la == lb) { solveManifold(m); return; } // same level: ordinary two-sided solve
    const bool upperIsA = la > lb;
    const int ui = upperIsA ? m.a : m.b;
    RigidBody& U = bodies_[ui];
    const float s = upperIsA ? 1.0f : -1.0f; // impulses are expressed for A (normal from B to A)
    // Velocity correction of the upper body only. Uses its own accumulators (starting at zero) so the
    // warm-start impulses stay those of the two-sided solve; friction is left to the regular
    // iterations. The upper body against a frozen support is a tiny problem (<= 4 points, one
    // body): iterate it locally until converged.
    float acc[16] = {0};
    const int np = std::min<int>(int(m.points.size()), 16);
    for (int local = 0; local < 8; ++local)
        for (int k = 0; k < np; ++k) {
            SolverPoint& p = m.points[k];
            Vector3 r = p.position - U.pos;
            Vector3 vu = U.velocityAt(p.position);
            Vector3 vl(0.0f);
            if (upperIsA) { if (m.b >= 0) vl = bodies_[m.b].velocityAt(p.position); }
            else vl = bodies_[m.a].velocityAt(p.position);
            float vn = dot(upperIsA ? vu - vl : vl - vu, p.normal);
            float kk = U.invMass + dot(cross(U.applyInvInertiaWorld(cross(r, p.normal)), r), p.normal);
            if (kk < 1e-12f) continue;
            float target = std::min(p.velocityBias, 0.0f); // speculative gap may still close; never push apart
            float old = acc[k];
            acc[k] = std::max(old + (target - vn) / kk, 0.0f);
            applyImpulse(ui, p.normal * ((acc[k] - old) * s), p.position);
        }

    // Friction in the same one-directional pass: the upper body is dragged along by (or held on)
    // its frozen support up to the Coulomb limit, so a sideways pull reaches the top of a column in
    // one sweep instead of one level per iteration (which shears a stack into a "staircase").
    if (!shockFrictionPass_) return;
    float normalTotal = 0;
    Vector3 c(0.0f), nrm(0.0f);
    for (int k = 0; k < np; ++k) {
        normalTotal += m.points[k].jn + acc[k];
        c += m.points[k].position;
        nrm += m.points[k].normal;
    }
    if (normalTotal <= 0) return;
    c /= float(np);
    nrm = normalize(nrm);
    Vector3 vu = U.velocityAt(c);
    Vector3 vl(0.0f);
    if (upperIsA) { if (m.b >= 0) vl = bodies_[m.b].velocityAt(c); }
    else vl = bodies_[m.a].velocityAt(c);
    Vector3 vrel = vu - vl;                          // velocity of the upper body relative to its support
    Vector3 vt = vrel - nrm * dot(vrel, nrm);
    float vtl = length(vt);
    if (vtl < 1e-7f) return;
    // Pure translation (no torque): a correction pass must not spin bodies up; tipping torques are
    // left to the regular two-sided friction.
    Vector3 t = vt / vtl;
    // Two limits:
    //  * the support's friction has to move everything resting on it: in equilibrium the normal
    //    impulse equals the weight of the whole column above times h, so the Coulomb-limited
    //    velocity change of that column is mu * |g| * h per step, independent of its height (this
    //    also keeps the pure translation - it ignores the rotation - a small correction);
    //  * the Coulomb cone of this contact: mu times its normal impulse, minus the friction impulse
    //    the regular iterations have already applied. The pass only finishes what they could not
    //    carry up the column; it never adds friction beyond Coulomb (a sliding body keeps mu_k N).
    const float mu = vtl < params.stickVelocity ? m.staticFriction : m.friction;
    const float used = length(m.t1 * m.jt1 + m.t2 * m.jt2);
    const float budget = std::max(0.0f, mu * normalTotal - used);
    const float dv = std::min({vtl, mu * length(params.gravity) * lastDt_, budget * U.invMass});
    U.vel -= t * dv;
}

void RigidWorld::prepare(float dt) {
    contactCount_ = 0;
    // Gauss-Seidel order: bottom-up (along gravity), so support propagates through stacks within
    // the first iterations instead of one level per iteration.
    Vector3 up = normalize(-params.gravity);
    if (length2(up) < 0.5f) up = Vector3(0, 1, 0);
    auto height = [&](const Manifold& m) {
        float h = dot(bodies_[m.a].pos, up);
        return m.b >= 0 ? std::min(h, dot(bodies_[m.b].pos, up)) : h - 1e3f; // static contacts first
    };
    std::stable_sort(manifolds_.begin(), manifolds_.end(),
                     [&](const Manifold& x, const Manifold& y) { return height(x) < height(y); });
    for (const Manifold& m : manifolds_) contactCount_ += m.points.size();
    buildColors();
    forEachManifold([&](Manifold& m) { prepareManifold(m, dt); });
}

void RigidWorld::solve() {
    forEachManifold([&](Manifold& m) { solveManifold(m); });
}

void RigidWorld::step(float dt) {
    if (bodies_.empty()) return;
    if (params.solver == RigidSolver::XPBD) {
        // The XPBD path has no islands: every body is awake (a body that fell asleep under the
        // impulse solver would otherwise keep the flag while it moves, and be frozen mid-air by the
        // next impulse step; the gas also reuses the cells of "resting" bodies).
        for (RigidBody& b : bodies_) {
            b.sleeping = false;
            b.sleepTimer = 0;
            b.sleepIsland = -1;
        }
        stepXPBD(dt);
        return;
    }
    lastDt_ = dt;
    const auto tStart = std::chrono::steady_clock::now();
    for (RigidBody& b : bodies_) b.updateInertia(); // orientation may have been edited externally
    if (!params.sleeping)
        for (RigidBody& b : bodies_) b.sleeping = false;
    if (grab_.active && grab_.body < int(bodies_.size())) wake(grab_.body);
    freezeSleepers();
    for (RigidBody& b : bodies_) {
        b.prevPos = b.pos;
        b.prevRot = b.rot;
    }
    for (RigidBody& b : bodies_) {
        if (b.invMass == 0) continue;
        b.vel += (params.gravity + b.force * b.invMass) * dt;
        b.angVel += b.applyInvInertiaWorld(b.torque) * dt;
    }
    collide();
    if (params.sleeping && updateIslands(false, dt)) collide(); // woken island: contacts among its bodies
    auto ts = std::chrono::steady_clock::now();
    prepare(dt);
    prepareGrab(dt);
    for (auto& j : joints_) j->prepare(bodies_, dt, params.warmStarting);
    for (int it = 0; it < params.iterations; ++it) {
        solve();
        for (auto& j : joints_) j->solveVelocity(bodies_);
        solveGrab(dt);
    }
    if (params.shockPropagation && !manifolds_.empty()) {
        computeLevels();
        // Ground-up order: sort by the lower level of each manifold.
        std::vector<int> order(manifolds_.size());
        for (size_t i = 0; i < order.size(); ++i) order[i] = int(i);
        auto lvl = [&](const Manifold& m) {
            int la = bodies_[m.a].invMass == 0 ? -1 : levels_[m.a];
            int lb = (m.b < 0 || bodies_[m.b].invMass == 0) ? -1 : levels_[m.b];
            return std::min(la, lb);
        };
        std::stable_sort(order.begin(), order.end(), [&](int x, int y) { return lvl(manifolds_[x]) < lvl(manifolds_[y]); });
        for (int pass = 0; pass < params.shockIterations; ++pass) {
            shockFrictionPass_ = params.shockFriction && pass == params.shockIterations - 1; // once per step
            for (int i : order) solveManifoldShock(manifolds_[i]);
        }
    }
    timings_.solve = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - ts).count();

    // Pairs of sleeping bodies produce no manifolds; keep their last impulses so a woken stack
    // carries its weight immediately instead of sagging and being thrown apart.
    std::unordered_map<uint64_t, CachedPair> keep;
    for (auto& [k, c] : cache_) {
        int a = int(k >> 32), b = int(k & 0xffffffffu) - 64;
        bool aSleep = a < int(bodies_.size()) && (bodies_[a].sleeping || bodies_[a].mass <= 0);
        bool bSleep = b < 0 || (b < int(bodies_.size()) && (bodies_[b].sleeping || bodies_[b].mass <= 0));
        if (aSleep && bSleep) keep.emplace(k, std::move(c));
    }
    cache_.swap(keep);
    cache_.reserve(cache_.size() + manifolds_.size());
    for (const Manifold& m : manifolds_) {
        CachedPair& c = cache_[key(m.a, m.b)];
        c.points = m.points;
        c.friction = m.t1 * m.jt1 + m.t2 * m.jt2;
        c.twist = m.jtwist;
        c.roll = m.jroll;
        c.locked = m.locked;
        c.lockRef = m.lockRef;
        c.lock = m.jlock;
    }

    // Resting-contact damping (cf. Bullet's additional damping): bodies that touch something and
    // move slower than the thresholds get an opposing velocity change, capped so that it only
    // eats the residual jitter and never real motion.
    if (params.restDamping > 0) {
        std::vector<char> touching(bodies_.size(), 0);
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
    const float ld = std::max(0.0f, 1.0f - params.linearDamping * dt);
    const float ad = std::max(0.0f, 1.0f - params.angularDamping * dt);
    for (RigidBody& b : bodies_) {
        b.force = Vector3(0.0f);
        b.torque = Vector3(0.0f);
        if (b.invMass == 0) continue;
        b.vel *= ld;
        b.angVel *= ad;
        b.pos += (b.vel + b.biasVel) * dt;
        b.rot = b.rot.integrated(b.angVel + b.biasAngVel, dt);
        b.biasVel = Vector3(0.0f);
        b.biasAngVel = Vector3(0.0f);
        b.updateInertia();
        // Sleep timer: time continuously spent below the thresholds.
        bool slow = length(b.vel) < params.sleepLinear && length(b.angVel) < params.sleepAngular;
        b.sleepTimer = slow ? b.sleepTimer + dt : 0.0f;
    }
    // Continuous collision: fast bodies are clamped to their time of impact (no tunnelling).
    auto tc = std::chrono::steady_clock::now();
    continuousCollision();
    auto ti = std::chrono::steady_clock::now();
    timings_.ccd = std::chrono::duration<float, std::milli>(ti - tc).count();
    // Joint position stage (nonlinear Gauss-Seidel on the new poses), frozen sleepers still static.
    solveJointPositions();
    unfreezeAll();
    if (params.sleeping) updateIslands(true, dt);
    auto tEnd = std::chrono::steady_clock::now();
    timings_.islands = std::chrono::duration<float, std::milli>(tEnd - ti).count();
    timings_.total = std::chrono::duration<float, std::milli>(tEnd - tStart).count();
}

} // namespace rf
