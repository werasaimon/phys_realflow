// Fetecau, Marsden, Ortiz & West (2003), sections 3.1-3.2: resolve an impact configuration/time
// and enforce tangential momentum and energy jump conditions. Here the exact discrete action
// is available for constant-force sphere translation (VariationalFlight), while CCD supplies
// an approximate event location. This is a restricted experimental path, not the paper's full
// rotating-body algorithm. Existing SI stepping remains the default.
#include "rigid/VariationalImpactStep.h"
#include "rigid/VariationalFlight.h"
#include "rigid/ConservativeAdvancement.h"
#include "rigid/RigidStepCheckpoint.h"
#include "core/Probe.h"

namespace rf {
namespace {
bool finite(const Vector3& v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

void includeEvent(ToiResult& first, const ToiResult& candidate) {
    if (candidate.status == ToiStatus::Unresolved) { first = candidate; return; }
    if (first.status == ToiStatus::Unresolved) return;
    if (candidate.status == ToiStatus::InitialContact || candidate.s < first.s) first = candidate;
}
} // namespace

bool RigidWorld::tryVariationalStep(float dt) {
    lastStepResult_ = VariationalImpactStep(*this).advance(dt);
    reportStep();
    Probe::set("rigid/variational events", double(lastStepResult_.impactEvents));
    return lastStepResult_.completed();
}

RigidStepStatus VariationalImpactStep::validate() const {
    const auto& w = world_;
    const auto& p = w.params;
    if (!finite(p.gravity) || !std::isfinite(p.ccdTolerance) || p.ccdTolerance <= 0
        || !std::isfinite(p.contactMargin) || p.contactMargin < 2 * p.ccdTolerance
        || p.ccdMaxIterations <= 0 || p.iterations <= 0) return RigidStepStatus::InvalidInput;
    if (p.solver != RigidSolver::SequentialImpulse || p.sleeping || p.linearDamping != 0
        || p.angularDamping != 0 || p.restDamping != 0 || p.rollingResistance != 0
        || !w.joints_.empty() || w.grab_.active) return RigidStepStatus::Unsupported;
    for (const auto& b : w.bodies_) {
        if (!b.alive) continue;
        if (!finite(b.pos) || !finite(b.vel) || !finite(b.angVel) || !finite(b.force)
            || !finite(b.torque) || !std::isfinite(b.mass) || !std::isfinite(b.invMass)
            || !finite(b.invInertiaLocal) || !std::isfinite(length2(b.angVel))) return RigidStepStatus::NonFinite;
        const double norm = double(b.rot.w) * b.rot.w + double(b.rot.x) * b.rot.x
                            + double(b.rot.y) * b.rot.y + double(b.rot.z) * b.rot.z;
        if (!b.shape || !std::isfinite(norm) || std::fabs(norm - 1) > 1e-3) return RigidStepStatus::InvalidInput;
        if (b.invMass == 0) {
            if (length2(b.vel) + length2(b.angVel) != 0) return RigidStepStatus::Unsupported;
            continue;
        }
        if (b.type() != ShapeType::Sphere || b.sleeping || b.restitution != 1
            || b.friction != 0 || b.staticFriction != 0 || length2(b.torque) != 0
            || length2(b.biasVel) + length2(b.biasAngVel) != 0
            || b.invInertiaLocal.x != b.invInertiaLocal.y || b.invInertiaLocal.y != b.invInertiaLocal.z)
            return RigidStepStatus::Unsupported;
        if (b.mass <= 0 || b.invMass <= 0 || b.invInertiaLocal.x <= 0
            || std::fabs(double(b.mass) * b.invMass - 1) > 1e-5) return RigidStepStatus::InvalidInput;
    }
    return RigidStepStatus::Completed;
}

RigidStepResult VariationalImpactStep::advance(float dt) {
    result_ = RigidStepResult(); result_.requestedDt = dt;
    if (!std::isfinite(dt) || dt <= 0) { result_.status = RigidStepStatus::InvalidInput; return result_; }
    result_.status = validate();
    if (!result_.completed()) return result_;
    const auto start = std::chrono::steady_clock::now();
    RigidStepCheckpoint initial(world_);
    profile_ = RigidWorld::Timings();
    world_.ccdDiagnostics_ = CcdDiagnostics(); world_.ccdHits_ = 0;
    result_.status = run(dt);
    const CcdDiagnostics diagnostics = world_.ccdDiagnostics_;
    if (!result_.completed()) {
        initial.restore(world_); result_.acceptedSubsteps = result_.impactEvents = 0; result_.maxContactDepth = 0;
    } else {
        result_.acceptedDt = dt;
        for (auto& b : world_.bodies_) b.force = b.torque = Vector3(0);
        world_.cache_.clear(); // event impulses are not persistent support forces
        world_.lastDt_ = dt;
        world_.updateWorldTree();
    }
    world_.ccdDiagnostics_ = diagnostics;
    world_.timings_ = profile_;
    world_.timings_.total = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - start).count();
    return result_;
}

Vector3 VariationalImpactStep::acceleration(int body) const {
    if (body < 0 || world_.bodies_[size_t(body)].invMass == 0) return Vector3(0);
    const auto& b = world_.bodies_[size_t(body)];
    return world_.params.gravity + b.force * b.invMass;
}

void VariationalImpactStep::fly(float h) {
    for (size_t i = 0; i < world_.bodies_.size(); ++i) {
        auto& b = world_.bodies_[i];
        if (b.alive && b.invMass > 0) VariationalFlight::advance(b, acceleration(int(i)), h);
    }
    ++result_.acceptedSubsteps;
}

RigidStepStatus VariationalImpactStep::run(float dt) {
    float remaining = dt;
    while (remaining > 0) {
        if (result_.acceptedSubsteps + result_.impactEvents >= 256) return RigidStepStatus::EventLimit;
        world_.collide();
        profile_.broad += world_.timings_.broad; profile_.narrow += world_.timings_.narrow;
        const auto start = std::chrono::steady_clock::now();
        const bool resolved = solveImpact();
        profile_.solve += std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - start).count();
        if (!resolved) return RigidStepStatus::Unresolved;
        if (!consumeFlight(remaining)) return result_.status;
        for (const auto& b : world_.bodies_) {
            if (!finite(b.pos) || !finite(b.vel) || !finite(b.angVel) || !std::isfinite(b.rot.w)
                || !std::isfinite(b.rot.x) || !std::isfinite(b.rot.y) || !std::isfinite(b.rot.z)) return RigidStepStatus::NonFinite;
        }
    }
    return RigidStepStatus::Completed;
}

double VariationalImpactStep::kinetic() const {
    double total = 0;
    for (const auto& b : world_.bodies_) if (b.alive && b.invMass > 0)
        for (int k = 0; k < 3; ++k) total += 0.5 * b.mass * double(b.vel[k]) * b.vel[k];
    return total;
}

bool VariationalImpactStep::solveImpact() {
    const float band = world_.params.ccdTolerance * 1.01f;
    bool approaching = false;
    world_.contactCount_ = 0;
    for (auto& m : world_.manifolds_) {
        int count = 0;
        for (auto p : m.points) {
            if (p.depth < -band) continue;
            if (p.depth > band) return false; // an invalid initial overlap is not repaired by this model
            const float vn = dot(world_.relativeVelocity(m, p.position), p.normal);
            p.velocityBias = std::max(-vn, 0.0f); // elastic Newton jump, one configuration, no force kick
            p.massN = 1 / (world_.bodies_[m.a].invMass + (m.b >= 0 ? world_.bodies_[m.b].invMass : 0));
            p.jn = p.jp = p.positionBias = 0;
            m.points[count++] = p;
            approaching |= vn < -1e-6f;
            result_.maxContactDepth = std::max(result_.maxContactDepth, p.depth);
        }
        m.points.count = count;
        world_.contactCount_ += size_t(count);
    }
    if (!approaching) return true;
    const double before = kinetic();
    for (int it = 0; it < world_.params.iterations; ++it) for (auto& m : world_.manifolds_) for (auto& p : m.points) {
        const Vector3 relative = world_.bodies_[m.a].vel - (m.b >= 0 ? world_.bodies_[m.b].vel : Vector3(0));
        const float old = p.jn;
        p.jn = std::max(0.0f, old + p.massN * (p.velocityBias - dot(relative, p.normal)));
        const Vector3 impulse = p.normal * (p.jn - old);
        world_.bodies_[m.a].vel += impulse * world_.bodies_[m.a].invMass;
        if (m.b >= 0) world_.bodies_[m.b].vel -= impulse * world_.bodies_[m.b].invMass;
    }
    ++result_.impactEvents;
    return validateImpulses(before);
}

bool VariationalImpactStep::validateImpulses(double before) const {
    if (!std::isfinite(kinetic()) || std::fabs(kinetic() - before) > 2e-5 * std::max(1.0, before)) return false;
    for (const auto& m : world_.manifolds_) for (const auto& p : m.points) {
        const Vector3 relative = world_.bodies_[m.a].vel - (m.b >= 0 ? world_.bodies_[m.b].vel : Vector3(0));
        const float residual = dot(relative, p.normal) - p.velocityBias;
        const float tolerance = 2e-5f * std::max(1.0f, length(relative));
        if (residual < -tolerance || (p.jn > 0 && std::fabs(residual) > tolerance)) return false;
    }
    return true;
}

void VariationalImpactStep::makeSweeps(float h) {
    sweeps_.resize(world_.bodies_.size()); bounds_.resize(sweeps_.size());
    for (size_t i = 0; i < sweeps_.size(); ++i) {
        const auto& b = world_.bodies_[i];
        SweptPose sweep; sweep.shape = b.shape.get(); sweep.p0 = sweep.p1 = b.pos; sweep.q0 = b.rot;
        if (b.alive && b.invMass > 0) {
            RigidBody end = b; VariationalFlight::advance(end, acceleration(int(i)), h);
            sweep.p1 = end.pos; sweep.translationCurve = acceleration(int(i)) * (0.5f * h * h);
        }
        sweeps_[i] = sweep;
        bounds_[i] = sweptBounds(sweep, world_.params.ccdTolerance);
    }
}

ToiResult VariationalImpactStep::pairEvent(const SweptPose& a, const SweptPose& b) const {
    ConservativeAdvancement query(world_.params.ccdTolerance, world_.params.ccdMaxIterations, true);
    if (b.shape->type() != ShapeType::Compound) return query.between(a, b);
    ToiResult first; first.status = ToiStatus::Separated; first.s = 1;
    for (const auto& child : static_cast<const CompoundShape*>(b.shape)->children()) {
        SweptPose part = b; part.shape = child.shape.get(); part.partR = child.R; part.partT = child.t;
        includeEvent(first, query.between(a, part));
    }
    return first;
}

ToiResult VariationalImpactStep::firstEvent(float h) {
    makeSweeps(h);
    ToiResult first; first.status = ToiStatus::Separated; first.s = 1;
    auto observe = [&](const ToiResult& r) { world_.ccdDiagnostics_.observe(r); includeEvent(first, r); };
    const auto& w = world_;
    for (size_t i = 0; i < sweeps_.size(); ++i) {
        if (!w.bodies_[i].alive || w.bodies_[i].invMass == 0) continue;
        for (size_t j = 0; j < sweeps_.size(); ++j) {
            if (j == i || !w.bodies_[j].alive || (w.bodies_[j].invMass > 0 && j < i)) continue;
            if (bounds_[i].overlaps(bounds_[j])) observe(pairEvent(sweeps_[i], sweeps_[j]));
        }
        if (w.params.collideWithDomain) {
            const Vector3 normals[] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
            for (int k = 0; k < 6; ++k) {
                const float offset = dot(normals[k], k % 2 ? w.domain_.hi : w.domain_.lo);
                ConservativeAdvancement query(w.params.ccdTolerance, w.params.ccdMaxIterations, true);
                observe(query.againstPlane(sweeps_[i], normals[k], offset));
            }
        }
        if (w.mesh_ && !w.mesh_->empty()) w.mesh_->bvh().queryAABB(bounds_[i], [&](uint32_t t) {
            Vector3 a, b, c; w.mesh_->triangle(t, a, b, c); TriangleShape triangle(a, b, c);
            SweptPose mesh; mesh.shape = &triangle;
            observe(pairEvent(sweeps_[i], mesh));
        });
    }
    return first;
}

bool VariationalImpactStep::consumeFlight(float& remaining) {
    float h = remaining;
    ToiStatus last = ToiStatus::Unresolved;
    for (int retry = 0; retry <= world_.params.ccdMaxSubdivisions; ++retry) {
        const auto start = std::chrono::steady_clock::now();
        const ToiResult event = firstEvent(h);
        profile_.ccd += std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - start).count();
        last = event.status;
        if (event.status == ToiStatus::Separated || (event.status == ToiStatus::Impact && event.s > 0)) {
            const float consumed = h * event.s;
            if (consumed > 0 && remaining - consumed < remaining) {
                fly(consumed); remaining -= consumed; return true;
            }
        }
        ++result_.rejectedTrials;
        h *= 0.5f;
    }
    result_.status = last == ToiStatus::InitialContact ? RigidStepStatus::PersistentContact : RigidStepStatus::Unresolved;
    return false;
}

} // namespace rf
