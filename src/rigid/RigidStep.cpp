// Reintegrate the complete dynamical state, never interpolate the pose while retaining an
// end-of-step spin. Mirtich, Timewarp Rigid Body Simulation (2000), section 2.1 motivates state
// rollback; this is synchronous binary subdivision, not his asynchronous Timewarp algorithm.
#include "rigid/RigidStep.h"
#include "rigid/RigidStepCheckpoint.h"
#include "core/Probe.h"

#include <cmath>
#include <array>
#include <stdexcept>

namespace rf {
namespace {

bool finiteState(const RigidWorld& w) {
    for (const RigidBody& b : w.bodies()) {
        if (!b.alive) continue;
        for (const Vector3& v : {b.pos, b.vel, b.angVel})
            if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z)) return false;
        if (!std::isfinite(b.rot.w) || !std::isfinite(b.rot.x) || !std::isfinite(b.rot.y) || !std::isfinite(b.rot.z)) return false;
    }
    return std::isfinite(w.kineticEnergy());
}

void addDiagnostics(CcdDiagnostics& total, const CcdDiagnostics& trial) {
    total.queries += trial.queries;
    total.unresolved += trial.unresolved;
    total.initialContacts += trial.initialContacts;
    total.passLimitReached |= trial.passLimitReached;
}

void addTimings(RigidWorld::Timings& total, const RigidWorld::Timings& trial) {
    total.broad += trial.broad; total.narrow += trial.narrow; total.solve += trial.solve;
    total.ccd += trial.ccd; total.islands += trial.islands;
}

} // namespace

void RigidWorld::step(float dt) {
    if (!tryStep(dt)) throw std::runtime_error("RigidWorld step rejected; state restored (inspect lastStepResult)");
}

bool RigidWorld::tryStep(float dt) {
    lastStepResult_ = RigidStep(*this).advance(dt);
    motorWorkTotal_ += lastStepResult_.motorWork;
    if (lastStepResult_.status == RigidStepStatus::InvalidInput) {
        ccdHits_ = 0; ccdDiagnostics_ = CcdDiagnostics(); timings_ = Timings();
    }
    reportStep();
    Probe::add("rigid/step retries", double(lastStepResult_.rejectedTrials));
    Probe::set("rigid/step rejected", lastStepResult_.completed() ? 0 : 1);
    return lastStepResult_.completed();
}

RigidStepResult RigidStep::advance(float dt) {
    RigidStepResult result; result.requestedDt = dt;
    if (!std::isfinite(dt) || dt <= 0) { result.status = RigidStepStatus::InvalidInput; return result; }
    RigidWorld& w = world_;
    if (!w.params.ccd || w.params.solver != RigidSolver::SequentialImpulse || w.bodies_.empty()) {
        w.stepUnchecked(dt, dt); result.acceptedDt = dt; result.acceptedSubsteps = 1;
        if (w.params.measureMotorWork && w.params.solver == RigidSolver::SequentialImpulse && !w.bodies_.empty())
            result.motorWork = w.trialMotorWork_;
        return result;
    }
    const auto start = std::chrono::steady_clock::now();
    if (!w.stepInitial_) w.stepInitial_ = std::make_shared<RigidStepCheckpoint>(w);
    else w.stepInitial_->save(w);
    RigidStepCheckpoint& initial = *w.stepInitial_;
    std::array<std::pair<float, int>, 512> pending;
    pending[0] = {dt, 0};
    size_t pendingCount = 1;
    CcdDiagnostics diagnostics;
    RigidWorld::Timings timings;
    size_t hits = 0;
    while (pendingCount) {
        if (result.acceptedSubsteps + result.rejectedTrials >= 1024) { result.status = RigidStepStatus::SubdivisionLimit; break; }
        const auto [h, depth] = pending[--pendingCount];
        initial.restoreLoads(w);
        const RigidStepCheckpoint* before = &initial;
        if (result.acceptedSubsteps > 0) {
            if (!w.stepTrial_) w.stepTrial_ = std::make_shared<RigidStepCheckpoint>(w);
            else w.stepTrial_->save(w);
            before = w.stepTrial_.get();
        }
        w.stepUnchecked(h, dt);
        addTimings(timings, w.timings_);
        addDiagnostics(diagnostics, w.ccdDiagnostics_); hits += w.ccdHits_;
        const bool finite = finiteState(w);
        if (w.ccdHits_ == 0 && finite) {
            ++result.acceptedSubsteps;
            if (w.params.measureMotorWork) result.motorWork += w.trialMotorWork_;
            result.maxContactDepth = std::max(result.maxContactDepth, w.deepestPenetration());
            continue;
        }
        ++result.rejectedTrials;
        before->restore(w);
        if (diagnostics.unresolved && w.params.ccdMaxIterations <= 0) { result.status = RigidStepStatus::Unresolved; break; }
        const float half = h * 0.5f;
        if (depth >= w.params.ccdMaxSubdivisions || pendingCount + 2 > pending.size() || half <= 0 || half + half != h) {
            result.status = finite ? RigidStepStatus::SubdivisionLimit : RigidStepStatus::NonFinite; break;
        }
        pending[pendingCount++] = {half, depth + 1}; pending[pendingCount++] = {half, depth + 1};
    }
    if (!result.completed()) {
        initial.restore(w); result.acceptedSubsteps = 0; result.maxContactDepth = 0; result.motorWork = 0;
    }
    else result.acceptedDt = dt;
    w.ccdDiagnostics_ = diagnostics; w.ccdHits_ = hits;
    w.timings_ = timings;
    w.timings_.total = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - start).count();
    return result;
}

} // namespace rf
