// Continuous collision detection: the time of impact of two moving shapes by conservative
// advancement (the distance shrinks no faster than the bound on the relative speed), and the
// trial-checking stage of RigidWorld for fast bodies. See TimeOfImpact.h for its limitations.
#include "rigid/TimeOfImpact.h"
#include "rigid/ConservativeAdvancement.h"
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
    const Vector3 padding = Vector3(tolerance + arc) + vabs(sweep.translationCurve) * 0.25f;
    box.lo -= padding;
    box.hi += padding;
    return box;
}

ToiResult timeOfImpact(const SweptPose& A, const SweptPose& B, float tol, int maxIt) {
    return ConservativeAdvancement(tol, maxIt).between(A, B);
}

ToiResult timeOfImpactPlane(const SweptPose& A, const Vector3& n, float d, float tol, int maxIt) {
    return ConservativeAdvancement(tol, maxIt).againstPlane(A, n, d);
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

// Inspect every child pair, including uncertainty in a pair later than the earliest impact.
float partsTimeOfImpact(const SweptPose& A, const SweptPose& B, float tol, int maxIt, CcdDiagnostics& diagnostics) {
    std::vector<SweptPose> pa, pb;
    partsOf(A, pa);
    partsOf(B, pb);
    float first = 1;
    ConservativeAdvancement query(tol, maxIt);
    for (const SweptPose& a : pa)
        for (const SweptPose& b : pb)
            first = std::min(first, diagnostics.observe(query.between(a, b)));
    return first;
}

} // namespace

// Check the completed trial trajectory without changing any pose or velocity. A hit or an
// unresolved query rejects the entire trial in RigidStep, which restores and reintegrates it.
// Candidate selection still uses the motion threshold (or ccdAllMoving); initial contacts are
// handled by the discrete solver. This is not an interval-certified collision-free integrator.
void RigidWorld::continuousCollision() {
    Probe::Timer timer("rigid/ccd ms");
    ccdHits_ = 0;
    ccdDiagnostics_ = CcdDiagnostics();
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
    impactTimes_.assign(n, 1.0f);
    findTimesOfImpact(impactTimes_, useTree);
    for (size_t i = 0; i < n; ++i) {
        if (impactTimes_[i] >= 1) continue;
        ccdClamped_[i] = 1; // compatibility diagnostic: this body's trial needs rejection
        ++ccdHits_;
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
        if (!b.alive || b.invMass == 0 || b.sleeping) continue;
        const Quaternion dq = b.rot * b.prevRot.conjugate();
        const float turnBound = kPi * std::sqrt(dq.x * dq.x + dq.y * dq.y + dq.z * dq.z);
        const float reachBound = length(b.pos - b.prevPos) + turnBound * b.shape->boundingRadius();
        motionBound_[size_t(i)] = reachBound; // no point of the body moved further this step
        if (params.ccdAllMoving && reachBound > 0) { fastBodies_.push_back(i); continue; }
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

// Earliest impact/uncertainty for each selected body against every candidate body, domain
// wall or mesh triangle. Both bodies of a dynamic pair contribute to the rejection diagnostic.
void RigidWorld::findTimesOfImpact(std::vector<float>& sMin, bool useTree) {
    const float tol = params.ccdTolerance;
    const std::vector<SweptPose>& sw = sweeps_;
    const std::vector<AABB>& swBox = sweepBoxes_;
    for (int i : fastBodies_) {
        auto tryPair = [&](int j) {
            if (j == i || !swBox[size_t(j)].overlaps(swBox[size_t(i)])) return;
            if (!bodies_[j].alive) return;
            const float stop = partsTimeOfImpact(sw[size_t(i)], sw[size_t(j)], tol, params.ccdMaxIterations, ccdDiagnostics_);
            sMin[size_t(i)] = std::min(sMin[size_t(i)], stop);
            if (bodies_[j].invMass > 0) sMin[size_t(j)] = std::min(sMin[size_t(j)], stop);
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
                const ToiResult r = timeOfImpactPlane(sw[i], normals[w], dot(normals[w], points[w]), tol, params.ccdMaxIterations);
                sMin[i] = std::min(sMin[i], ccdDiagnostics_.observe(r));
            }
        }
        if (mesh_ && !mesh_->empty() && mesh_->bounds().overlaps(swBox[i])) {
            mesh_->bvh().queryAABB(swBox[i], [&](uint32_t t) {
                Vector3 a, b, c;
                mesh_->triangle(t, a, b, c);
                TriangleShape tri(a, b, c);
                SweptPose T;
                T.shape = &tri;
                const float stop = partsTimeOfImpact(sw[i], T, tol, params.ccdMaxIterations, ccdDiagnostics_);
                sMin[i] = std::min(sMin[i], stop);
            });
        }
    }
}

} // namespace rf
