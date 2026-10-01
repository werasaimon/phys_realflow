#pragma once
// PhysRealFlow adapter: fixed declared steps, optionally atomic all-moving CCD retries.
#include "BouncingCase.h"
#include "rigid/RigidWorld.h"
#include <array>

namespace rf::benchmark {
class BouncingSdk {
public:
    BouncingSdk(int iterations, bool ccd, bool variational = false) : variational_(variational) {
        auto& p = world_.params;
        p.iterations = iterations;
        p.substeps = 1;
        p.collideWithDomain = false;
        p.linearDamping = p.angularDamping = p.rollingResistance = 0;
        p.baumgarte = p.blockCfm = p.slop = p.restitutionThreshold = 0;
        p.sleeping = p.shockPropagation = p.rotationalLock = false;
        p.ccd = ccd;
        p.ccdAllMoving = ccd;
        world_.addBox({0,-5,0}, {25,5,25}, rf::Quaternion(), 0, {1,1,1});
        for (int i = 0; i < count; ++i) {
            const int id = world_.addSphere(startPosition(i), radius, 1, {1,1,1});
            auto& body = world_.bodies()[size_t(id)];
            body.mass = mass;
            body.invMass = 0.1f;
            body.invInertiaLocal = {25,25,25};
            body.updateInertia();
        }
        for (auto& body : world_.bodies()) {
            body.friction = body.staticFriction = 0;
            body.restitution = 1;
        }
    }

    bool step(float dt, Work& work) {
        const bool completed = variational_ ? world_.tryVariationalStep(dt) : world_.tryStep(dt);
        const auto& result = world_.lastStepResult();
        work.accepted += result.acceptedSubsteps;
        work.rejected += result.rejectedTrials;
        work.time += result.acceptedDt;
        work.sdkImpactEvents += result.impactEvents;
        work.sdkCcdQueries += world_.ccdDiagnostics().queries;
        const auto& times = world_.timings();
        work.sdkBroadMs += times.broad; work.sdkNarrowMs += times.narrow;
        work.sdkSolveMs += times.solve; work.sdkCcdMs += times.ccd;
        return completed;
    }

    std::array<State, count> states() const {
        std::array<State, count> out;
        for (int i = 0; i < count; ++i) {
            const auto& body = world_.bodies()[size_t(i + 1)];
            out[size_t(i)] = {body.pos, body.vel, body.angVel};
        }
        return out;
    }

private:
    bool variational_ = false;
    rf::RigidWorld world_;
};

} // namespace rf::benchmark
