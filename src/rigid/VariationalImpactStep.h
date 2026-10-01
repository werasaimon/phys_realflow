#pragma once
// Experimental variational impact flow: exact constant-force action between events, coupled
// elastic normal impulses at one configuration, and the full remaining time after each event.
// Dynamic bodies must be frictionless elastic spheres with isotropic inertia, zero torque,
// damping and sleeping; joints/grab are unsupported. Static convex/compound/mesh geometry is
// queried through the SDK. Persistent contact is an explicit atomic rejection, not a freeze.
// Float CCD and a nonzero contact band approximate event geometry; no symplecticity certificate
// for this finite-tolerance implementation is claimed. See docs/book/14-variational-impact.md.
#include "rigid/RigidWorld.h"

namespace rf {

class VariationalImpactStep {
public:
    explicit VariationalImpactStep(RigidWorld& world) : world_(world) {}
    RigidStepResult advance(float dt);

private:
    RigidWorld& world_;
    RigidStepResult result_;
    RigidWorld::Timings profile_;
    std::vector<SweptPose> sweeps_;
    std::vector<AABB> bounds_;
    RigidStepStatus validate() const;
    RigidStepStatus run(float dt);
    bool solveImpact();
    bool validateImpulses(double before) const;
    double kinetic() const;
    Vector3 acceleration(int body) const;
    void fly(float h);
    void makeSweeps(float h);
    ToiResult firstEvent(float h);
    ToiResult pairEvent(const SweptPose& a, const SweptPose& b) const;
    bool consumeFlight(float& remaining);
};

} // namespace rf
