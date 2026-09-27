// Continuous collision detection: the time of impact of two moving shapes by conservative
// advancement (the distance shrinks no faster than the bound on the relative speed), and the
// step of RigidWorld that stops fast bodies at their impact so nothing tunnels. See TimeOfImpact.h.
#include "rigid/TimeOfImpact.h"
#include "rigid/RigidWorld.h"

#include "core/Parallel.h"
#include "core/Probe.h"

#include <algorithm>

namespace rf {

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
            r.hit = true;         // overshot inside the tolerance band: report the last safe s
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

AABB boundsAt(const SweptPose& sp, float s) {
    PosedShape ps = sp.at(s);
    return ps.shape->boundsAt(ps.R, ps.p);
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
void RigidWorld::continuousCollision() {
    Probe::Timer timer("rigid/ccd ms");
    ccdHits_ = 0;
    ccdClamped_.assign(bodies_.size(), 0);
    if (!params.ccd) return;
    if (!anyFastBody()) return; // nothing fast this step (the usual case): no sweeps to test
    for (int pass = 0; pass < 8; ++pass) {
        sweepBodies();
        std::vector<float> sMin(bodies_.size(), 1.0f);
        findTimesOfImpact(sMin);
        if (!clampToTimesOfImpact(sMin)) break;
    }
}

bool RigidWorld::anyFastBody() const {
    for (const RigidBody& b : bodies_)
        if (b.invMass > 0 && !b.sleeping && isFast(b, sweptOf(b), params.ccdThreshold)) return true;
    return false;
}

// The sweep of every body over this step and the box that contains it from start to end,
// widened by the tolerance (kept between steps: no allocation).
void RigidWorld::sweepBodies() {
    const float tol = params.ccdTolerance;
    const int n = int(bodies_.size());
    std::vector<SweptPose>& sw = sweeps_;
    std::vector<AABB>& swBox = sweepBoxes_;
    sw.resize(n);
    swBox.resize(n);
    for (int i = 0; i < n; ++i) {
        sw[i] = sweptOf(bodies_[i]);
        swBox[i] = boundsAt(sw[i], 0);
        swBox[i].expand(boundsAt(sw[i], 1));
        swBox[i].lo -= Vector3(tol);
        swBox[i].hi += Vector3(tol);
    }
}

// For every fast body, the earliest moment (0..1 of the step) it touches another body, a domain
// wall or the static mesh: conservative advancement (timeOfImpact) on the pairs a BVH over the
// swept boxes proposes (O(n log n) instead of all pairs). Both bodies of a pair stop together.
void RigidWorld::findTimesOfImpact(std::vector<float>& sMin) {
    const float tol = params.ccdTolerance;
    const int n = int(bodies_.size());
    const std::vector<SweptPose>& sw = sweeps_;
    const std::vector<AABB>& swBox = sweepBoxes_;
    BVH sweptTree;
    sweptTree.build(swBox, 2);
    for (int i = 0; i < n; ++i) {
        const RigidBody& body = bodies_[i];
        if (body.invMass == 0 || !isFast(body, sw[i], params.ccdThreshold)) continue;
        sweptTree.queryAABB(swBox[i], [&](uint32_t ju) {
            int j = int(ju);
            if (j == i || !swBox[j].overlaps(swBox[i])) return;
            ToiResult r = partsTimeOfImpact(sw[i], sw[j], tol);
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
        if (!ccdClamped_[i]) ++ccdHits_;
        ccdClamped_[i] = 1;
        any = true;
    }
    return any;
}

} // namespace rf
