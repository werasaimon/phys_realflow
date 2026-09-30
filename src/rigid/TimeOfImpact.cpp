// Continuous collision detection: the time of impact of two moving shapes by conservative
// advancement (the distance shrinks no faster than the bound on the relative speed), and the
// motion-clamping stage of RigidWorld for fast bodies. See TimeOfImpact.h for its limitations.
#include "rigid/TimeOfImpact.h"
#include "rigid/RigidWorld.h"

#include "core/Parallel.h"
#include "core/Probe.h"

#include <algorithm>

namespace rf {

AABB sweptBounds(const SweptPose& sweep, float tolerance) {
    const PosedShape start = sweep.at(0), end = sweep.at(1);
    AABB box = start.shape->boundsAt(start.R, start.p);
    box.expand(end.shape->boundsAt(end.R, end.p));
    // For each material point x(s), |x''(s)| <= |theta|^2 r. Linear interpolation
    // between its endpoints has error <= max|x''| s(1-s)/2 <= |theta|^2 r/8.
    // Translation is linear; r includes the offset of a compound child. See docs/22.
    const float arc = 0.125f * length(sweep.dTheta) * sweep.angularReach();
    box.lo -= Vector3(tolerance + arc);
    box.hi += Vector3(tolerance + arc);
    return box;
}

ToiResult timeOfImpact(const SweptPose& A, const SweptPose& B, float tol, int maxIt) {
    ToiResult r;
    const float bound = length((A.p1 - A.p0) - (B.p1 - B.p0)) + A.angularReach() + B.angularReach();
    if (bound < 1e-9f) return r;
    float s = 0;
    for (int it = 0; it < maxIt; ++it) {
        r.iterations = it + 1;
        GjkResult g = gjk(A.at(s), B.at(s));
        if (g.intersect) {
            if (s == 0) return r; // already touching: the discrete solver owns this pair
            r.hit = true;         // intersecting iterate; no earlier safe bracket is retained here
            r.s = s;
            return r;
        }
        if (g.distance < tol) {
            // Touching at the start of the step: a resting/sliding contact owned by the discrete
            // (speculative) solver - never freeze such pairs.
            if (s == 0) return r;
            r.hit = true;
            r.s = s;
            return r;
        }
        // Conservative step: the gap cannot close faster than `bound` per unit s. Aim slightly
        // short of contact so the next GJK still sees separated shapes.
        s += std::max((g.distance - 0.5f * tol) / bound, 1e-6f);
        if (s > 1.0f) return r;
    }
    return r;
}

ToiResult timeOfImpactPlane(const SweptPose& A, const Vector3& n, float d, float tol, int maxIt) {
    ToiResult r;
    const float bound = std::fabs(dot(A.p1 - A.p0, n)) + A.angularReach();
    if (bound < 1e-9f) return r;
    float s = 0;
    for (int it = 0; it < maxIt; ++it) {
        r.iterations = it + 1;
        PosedShape ps = A.at(s);
        float dist = dot(ps.support(-n), n) - d; // exact distance of the shape to the plane
        if (dist < 0) {
            if (s == 0) return r;
            r.hit = true;
            r.s = s;
            return r;
        }
        if (dist < tol) {
            if (s == 0) return r; // touching the plane already: discrete contact
            r.hit = true;
            r.s = s;
            return r;
        }
        s += std::max((dist - 0.5f * tol) / bound, 1e-6f);
        if (s > 1.0f) return r;
    }
    return r;
}

// ---------------------------------------------------------------------------
// Continuous collision in RigidWorld: conservative advancement of the fast bodies
// ---------------------------------------------------------------------------
namespace {

// The motion of a body over this step, from where it was to where the solver put it.
SweptPose sweptOf(const RigidBody& b) {
    SweptPose sp;
    sp.shape = b.shape.get();
    sp.p0 = b.prevPos;
    sp.p1 = b.pos;
    sp.q0 = b.prevRot;
    sp.dTheta = (b.rot * b.prevRot.conjugate()).log();
    return sp;
}

// Fast: the body moves (or its rim turns) further than a fraction of its own thickness in one
// step - it could pass through a thin obstacle between two discrete checks.
bool isFast(const RigidBody& b, const SweptPose& sp, float threshold) {
    Vector3 ext = b.shape->localBounds().extent();
    float minExtent = std::max(0.5f * minComp(ext), 1e-3f);
    return length(sp.p1 - sp.p0) + sp.angularReach() >= threshold * minExtent;
}

// A compound body is swept part by part: the hull of its union is not its surface, a fast body
// can fly into a concave region of that hull (a teapot's handle) without touching any part.
void partsOf(const SweptPose& sp, std::vector<SweptPose>& out) {
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
}

// Earliest time of impact of any part of A with any part of B.
ToiResult partsTimeOfImpact(const SweptPose& A, const SweptPose& B, float tol) {
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
}

} // namespace

// Fast bodies are stopped at the moment they would touch something, so nothing tunnels through
// a wall between two steps. Passes: find the times of impact of all fast bodies with the motions
// of this pass; both bodies of a colliding pair stop at their common time of impact (motion
// clamping); a stopped body stays put for the rest of the step (its sweep becomes static), and
// the pass repeats so bodies that would now run into it are caught too (A hits B, B hits C, ...).
//
// What it costs: only the fast bodies are swept against the others, as Box2D's "bullet" bodies and
// Bullet's ccdMotionThreshold do - the slow ones are left to the speculative contacts. The fast
// list is made once per step, not once per pass (a pass only refreshes the bodies it clamped). The
// candidates come from a plain scan while the fast bodies are few - a tree over a thousand swept
// boxes costs a millisecond to build, a scan of them for one fast body a microsecond - and a slow
// body gets its exact sweep only when its coarse box (coarseSweepBox) reaches a fast body's sweep.
// A pile of a thousand falling cubes spent 14 ms a frame here before.
void RigidWorld::continuousCollision() {
    Probe::Timer timer("rigid/ccd ms");
    ccdHits_ = 0;
    ccdClamped_.assign(bodies_.size(), 0);
    if (!params.ccd) return;
    // 1. Who is fast this step (the usual answer: nobody, and nothing else is done).
    if (!findFastBodies()) return;
    // 2. The sweeps: the fast bodies' now, the others' when a scan finds them near - or, with many
    //    fast bodies, every body's at once and a tree over them.
    const size_t n = bodies_.size();
    sweeps_.resize(n);
    sweepBoxes_.resize(n);
    const bool useTree = fastBodies_.size() * n > kCcdScanLimit;
    sweepReady_.assign(n, useTree ? 1 : 0);
    if (useTree) {
        sweepBodies();
        sweptTree_.build(sweepBoxes_, 2);
    } else {
        for (int i : fastBodies_) sweepBody(i), sweepReady_[size_t(i)] = 1;
    }
    // 3. The passes: times of impact, then the clamping, until nothing more is stopped.
    for (int pass = 0; pass < 8; ++pass) {
        impactTimes_.assign(bodies_.size(), 1.0f);
        findTimesOfImpact(impactTimes_, useTree);
        if (!clampToTimesOfImpact(impactTimes_)) break;
    }
}

// The bodies that move (or whose rim turns) further in this step than a fraction of their own
// thickness (isFast); sleeping and static ones never are. Returns whether there is any.
// A quick bound screens them first, without the logarithm of the turn: for the step's rotation
// dq = q1 q0^-1 with vector part v, the turn angle is 2 asin|v| <= pi |v|, so a body slow even by
// that bound is slow - only the others get the exact test (the same bodies come out).
bool RigidWorld::findFastBodies() {
    fastBodies_.clear();
    motionBound_.assign(bodies_.size(), 0.0f);
    for (int i = 0; i < int(bodies_.size()); ++i) {
        const RigidBody& b = bodies_[i];
        if (b.invMass == 0 || b.sleeping) continue;
        const Quaternion dq = b.rot * b.prevRot.conjugate();
        const float turnBound = kPi * std::sqrt(dq.x * dq.x + dq.y * dq.y + dq.z * dq.z);
        const float reachBound = length(b.pos - b.prevPos) + turnBound * b.shape->boundingRadius();
        motionBound_[size_t(i)] = reachBound; // no point of the body moved further this step
        const float minExtent = std::max(0.5f * minComp(b.shape->localBounds().extent()), 1e-3f);
        if (reachBound < params.ccdThreshold * minExtent) continue; // slow even by the bound
        if (isFast(b, sweptOf(b), params.ccdThreshold)) fastBodies_.push_back(i);
    }
    return !fastBodies_.empty();
}

// A box that surely holds body j's whole sweep over this step, made without the sweep: its box at
// the start of the step (made for the broad phase, already widened by the contact margin), grown
// by how far any of its points can have moved (motionBound_) and by the tolerance.
AABB RigidWorld::coarseSweepBox(int j) const {
    AABB box = boxes_[size_t(j)];
    const float grow = motionBound_[size_t(j)] + params.ccdTolerance;
    box.lo -= Vector3(grow);
    box.hi += Vector3(grow);
    return box;
}

// The sweep of one body over this step and the box that contains it from start to end, widened by
// the tolerance.
void RigidWorld::sweepBody(int i) {
    sweeps_[size_t(i)] = sweptOf(bodies_[i]);
    sweepBoxes_[size_t(i)] = sweptBounds(sweeps_[size_t(i)], params.ccdTolerance);
}

// The sweeps of all bodies (kept between steps: no allocation once grown).
void RigidWorld::sweepBodies() {
    sweeps_.resize(bodies_.size());
    sweepBoxes_.resize(bodies_.size());
    for (int i = 0; i < int(bodies_.size()); ++i) sweepBody(i);
}

// For every fast body not yet stopped, the earliest moment (0..1 of the step) it touches another
// body, a domain wall or the static mesh: conservative advancement (timeOfImpact) on the pairs
// whose swept boxes meet. The candidates come from the tree over the swept boxes (built once per
// step: a stopped body's new box lies inside its old one, so the tree stays conservative) or from a
// plain scan of all boxes. Both bodies of a pair stop together.
void RigidWorld::findTimesOfImpact(std::vector<float>& sMin, bool useTree) {
    const float tol = params.ccdTolerance;
    const std::vector<SweptPose>& sw = sweeps_;
    const std::vector<AABB>& swBox = sweepBoxes_;
    for (int i : fastBodies_) {
        if (ccdClamped_[size_t(i)]) continue; // stopped by an earlier pass: static now
        auto tryPair = [&](int j) {
            if (j == i || !swBox[size_t(j)].overlaps(swBox[size_t(i)])) return;
            ToiResult r = partsTimeOfImpact(sw[size_t(i)], sw[size_t(j)], tol);
            if (!r.hit) return;
            sMin[size_t(i)] = std::min(sMin[size_t(i)], r.s);
            if (bodies_[j].invMass > 0) sMin[size_t(j)] = std::min(sMin[size_t(j)], r.s); // the pair stops together
        };
        if (useTree) sweptTree_.queryAABB(swBox[size_t(i)], [&](uint32_t j) { tryPair(int(j)); });
        else
            for (int j = 0; j < int(bodies_.size()); ++j) {
                if (!sweepReady_[size_t(j)]) { // a slow body: its exact sweep only if it can be near
                    if (j == i || !coarseSweepBox(j).overlaps(swBox[size_t(i)])) continue;
                    sweepBody(j);
                    sweepReady_[size_t(j)] = 1;
                }
                tryPair(j);
            }
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
                ToiResult r = partsTimeOfImpact(sw[i], T, tol);
                if (r.hit) sMin[i] = std::min(sMin[i], r.s);
            });
        }
    }
}

// Every body with an impact before the end of the step is put back to that moment, keeping its
// velocity: the speculative contact of the next step resolves the impact (with restitution, see
// prepareManifold). Returns whether anything was clamped (then another pass is due).
bool RigidWorld::clampToTimesOfImpact(const std::vector<float>& sMin) {
    const std::vector<SweptPose>& sw = sweeps_;
    bool any = false;
    for (int i = 0; i < int(bodies_.size()); ++i) {
        if (sMin[i] >= 1.0f) continue;
        RigidBody& b = bodies_[i];
        b.pos = sw[i].p0 + (sw[i].p1 - sw[i].p0) * sMin[i];
        b.rot = sw[i].q0.integrated(sw[i].dTheta, sMin[i]);
        b.prevPos = b.pos; // static for the remaining passes of this step
        b.prevRot = b.rot;
        b.updateInertia();
        sweepBody(i); // its sweep is now a point: the next pass sees it standing there
        sweepReady_[size_t(i)] = 1;
        if (!ccdClamped_[i]) ++ccdHits_;
        ccdClamped_[i] = 1;
        any = true;
    }
    return any;
}

} // namespace rf
