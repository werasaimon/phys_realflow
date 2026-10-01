// Atomic CCD retries: an analytic load impulse, full rollback after a successful first half,
// and the mixed-scene energy overflow reproduced at a finer caller step than the demo uses.
#include "TestRunner.h"
#include "Tests.h"

#include <limits>
#include <stdexcept>

namespace {

void collisionWorld(RigidWorld& w) {
    w.params.gravity = Vector3(0);
    w.params.sleeping = w.params.collideWithDomain = false;
    w.params.linearDamping = w.params.angularDamping = 0;
    w.addBox(Vector3(0), {0.01f, 1, 1}, Quaternion(), 0, Vector3(1));
    const int bullet = w.addSphere({-1.5f, 0, 0}, 0.02f, 1000, Vector3(1));
    w.bodies()[bullet].vel = {100, 0, 0};
}

Vector3 angularMomentum(const RigidBody& b) {
    const Matrix3x3 r = b.rotation();
    return r * ((r.transposed() * b.angVel) / b.invInertiaLocal);
}

uint64_t mechanicalHash(const RigidWorld& w) {
    StateHash h;
    for (const RigidBody& b : w.bodies()) {
        h.add(b.pos); h.add(b.rot); h.add(b.vel); h.add(b.angVel);
        h.add(b.force); h.add(b.torque); h.add(b.biasVel); h.add(b.biasAngVel);
        h.add(b.prevPos); h.add(b.prevRot); h.add(b.invMass); h.add(b.invInertiaLocal);
        h.add(b.sleepTimer); h.add(float(b.sleeping)); h.add(float(b.sleepIsland));
    }
    h.add(w.grabJoint().impulse);
    return h.h;
}

double mechanicalEnergy(const RigidWorld& w) {
    double energy = w.kineticEnergy();
    for (const RigidBody& b : w.bodies()) if (b.mass > 0) energy -= b.mass * dot(w.params.gravity, b.pos);
    return energy;
}

// Independent face check for the sample's axis-aligned wall. If a body's full bounding box
// lies within the wall's y/z span, any positive intrusion past its left plane is penetration.
// Bodies leaving that span are outside this oracle's scope (the regular manifold metric remains).
float wallIntrusion(const RigidWorld& w) {
    const AABB wall = w.bodies()[0].worldBounds();
    float worst = 0;
    for (size_t i = 1; i < w.bodies().size(); ++i) {
        const RigidBody& body = w.bodies()[i];
        if (body.invMass == 0) continue;
        const AABB b = body.worldBounds();
        if (b.lo.y < wall.lo.y || b.hi.y > wall.hi.y || b.lo.z < wall.lo.z || b.hi.z > wall.hi.z) continue;
        // Measure actual intersection of the finite slab. A body may legitimately travel
        // around the wall and return on its far side; side alone is not a tunnelling oracle.
        worst = std::max(worst, std::min(b.hi.x - wall.lo.x, wall.hi.x - b.lo.x));
    }
    return worst;
}

} // namespace

void testRigidStepLoads() {
    RigidWorld w;
    collisionWorld(w);
    const int free = w.addSphere({5, 5, 0}, 0.1f, 1000, Vector3(1));
    const int spin = w.addBox({8, 5, 0}, {0.1f, 0.2f, 0.3f}, Quaternion(), 1000, Vector3(1));
    w.bodies()[free].force = {2, -3, 4}; w.bodies()[free].torque = {0.03f, 0.02f, -0.01f};
    w.bodies()[spin].angVel = {1, 2, 3};
    const RigidBody before = w.bodies()[free];
    const Vector3 spinL = angularMomentum(w.bodies()[spin]);
    const float dt = 0.02f;
    CHECK(w.tryStep(dt), "collision step failed");
    const RigidStepResult& result = w.lastStepResult();
    CHECK(result.rejectedTrials > 0 && result.acceptedSubsteps > 1 && result.acceptedDt == dt,
          "fixture did not reintegrate the caller's interval");
    const RigidBody& after = w.bodies()[free];
    CHECK(length(after.vel * after.mass - before.force * dt) < 1e-6f, "external force lost or counted twice");
    CHECK(length(angularMomentum(after) - before.torque * dt) < 1e-7f, "external torque lost or counted twice");
    CHECK(length(angularMomentum(w.bodies()[spin]) - spinL) < 1e-5f, "free anisotropic spin lost world angular momentum");
    CHECK(length2(after.force) + length2(after.torque) == 0, "completed step did not consume loads");
}

void testRigidStepRollback() {
    RigidWorld w, reference;
    for (RigidWorld* p : {&w, &reference}) {
        collisionWorld(*p);
        const int jointBody = p->addSphere({5, 1, 0}, 0.1f, 1000, Vector3(1));
        p->addBallJoint(jointBody, -1, {5, 2, 0});
        p->addBox({8, -0.1f, 0}, {1, 0.1f, 1}, Quaternion(), 0, Vector3(1));
        p->addBox({8, 0.1f, 0}, Vector3(0.1f), Quaternion(), 1000, Vector3(1));
        p->params.gravity = {0, -9.81f, 0};
        p->bodies()[jointBody].vel = {1, 0, 0};
        p->bodies()[1].vel = Vector3(0);
        for (int i = 0; i < 10; ++i) p->step(0.001f); // nonzero joint and contact warm starts
        p->bodies()[1].vel = {100, 0, 0};
        p->bodies()[jointBody].force = {2, 3, 4};
        p->params.ccdMaxSubdivisions = 1;
    }
    const Joint* jointIdentity = w.joints()[0].get();
    const uint64_t before = mechanicalHash(w);
    CHECK(!w.tryStep(0.02f), "insufficient subdivision budget accepted the step");
    CHECK(w.lastStepResult().status == RigidStepStatus::SubdivisionLimit && w.lastStepResult().rejectedTrials == 2,
          "fixture must accept its first half and reject its second half");
    CHECK(w.lastStepResult().acceptedDt == 0 && mechanicalHash(w) == before, "root rollback did not restore all bodies/loads");
    CHECK(w.joints()[0].get() == jointIdentity, "rollback invalidates caller-held joint reference");
    for (float dt : {0.0f, -1.0f, std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN()}) {
        CHECK(!w.tryStep(dt) && w.lastStepResult().status == RigidStepStatus::InvalidInput, "invalid dt accepted");
        CHECK(mechanicalHash(w) == before, "invalid dt changed the state");
    }
    bool threw = false;
    try { w.step(0.02f); } catch (const std::runtime_error&) { threw = true; }
    CHECK(threw && mechanicalHash(w) == before, "void step must report failure without advancing state");
    w.params.ccdMaxSubdivisions = reference.params.ccdMaxSubdivisions = 16;
    for (int i = 0; i < 30; ++i) {
        w.step(0.001f); reference.step(0.001f);
        CHECK(mechanicalHash(w) == mechanicalHash(reference), "rejected trial polluted subsequent dynamics at %d", i);
    }
}

void testRigidStepMixedEnergy() {
    for (bool allMoving : {false, true}) for (int factor : {1, 2, 4}) {
        Simulation sim;
        loadSample(sim, Preset::RigidCcd);
        RigidWorld& w = sim.rigid;
        w.params.ccdAllMoving = allMoving;
        const float dt = 1.0f / (600 * factor);
        const double initial = mechanicalEnergy(w);
        double peak = initial;
        size_t retries = 0, accepted = 0;
        float depth = 0, wallDepth = 0;
        for (int i = 0; i < 1200 * factor; ++i) {
            const bool complete = w.tryStep(dt);
            CHECK(complete, "mixed scene stalled: factor %d step %d status %d", factor, i, int(w.lastStepResult().status));
            if (!complete) break;
            const double energy = mechanicalEnergy(w);
            CHECK(std::isfinite(energy) && energy <= initial * 1.05, "energy overflow/gain: factor %d step %d E=%g", factor, i, energy);
            if (!std::isfinite(energy)) break;
            peak = std::max(peak, energy); depth = std::max(depth, w.lastStepResult().maxContactDepth);
            wallDepth = std::max(wallDepth, wallIntrusion(w));
            retries += w.lastStepResult().rejectedTrials; accepted += w.lastStepResult().acceptedSubsteps;
        }
        CHECK(retries > 0, "the energy regression did not exercise retries");
        if (allMoving) CHECK(depth < 0.002f && wallDepth < 0.002f, "all-moving sample exceeds 2 mm: manifold %g wall %g", depth, wallDepth);
        std::printf("  allMoving=%d dt=1/%d: E0=%.9g peak=%.9g final=%.9g J; manifold=%.6g wall=%.6g mm; accepted=%zu rejected=%zu\n",
                    int(allMoving), 600 * factor, initial, peak, mechanicalEnergy(w), depth * 1000, wallDepth * 1000, accepted, retries);
    }
}

// Baumgarte's per-step correction must decrease with a rejected trial's duration; otherwise
// an initially overlapping sphere is pushed into its neighbour by the same amount at every retry.
void testRigidStepRecovery() {
    RigidWorld w;
    w.params.gravity = Vector3(0);
    w.params.sleeping = w.params.collideWithDomain = false;
    w.params.linearDamping = w.params.angularDamping = 0;
    w.params.ccdAllMoving = true;
    w.addSphere({0, 0, 0}, 0.1f, 0, Vector3(1));
    const int moving = w.addSphere({0.1f, 0, 0}, 0.1f, 1000, Vector3(1));
    w.addSphere({0.32f, 0, 0}, 0.1f, 0, Vector3(1));
    CHECK(w.tryStep(0.01f), "initial-overlap recovery cannot shrink with the trial interval");
    const auto& r = w.lastStepResult();
    CHECK(r.rejectedTrials > 0 && r.acceptedSubsteps > 1, "recovery fixture did not require subdivision");
    CHECK(w.bodies()[moving].pos.x > 0.1f && w.bodies()[moving].pos.x < 0.13f, "wrong recovery displacement");
    CHECK(w.kineticEnergy() == 0, "split recovery added kinetic energy");
}

// Galilean invariance and the analytic central impact of equal masses: vA'=U+(1-e)c/2,
// vB'=U+(1+e)c/2. A common translation U must not change the transferred impulse.
void testRigidStepRelativeImpact() {
    for (float relative : {1.0f, 100.0f}) for (float boost : {0.0f, 99.0f}) {
        RigidWorld w;
        w.params.gravity = Vector3(0);
        w.params.sleeping = w.params.collideWithDomain = w.params.shockPropagation = false;
        w.params.linearDamping = w.params.angularDamping = w.params.restitutionThreshold = 0;
        w.params.ccdAllMoving = true;
        const int a = w.addSphere({-0.2f, 0, 0}, 0.1f, 1000, Vector3(1));
        const int b = w.addSphere({0.2f, 0, 0}, 0.1f, 1000, Vector3(1));
        w.bodies()[a].vel.x = boost + relative; w.bodies()[b].vel.x = boost;
        w.bodies()[a].restitution = w.bodies()[b].restitution = 0.5f;
        const float dt = 0.05f / relative;
        for (int i = 0; i < 10; ++i) {
            const bool complete = w.tryStep(dt);
            CHECK(complete, "relative impact rejected: c=%g U=%g", relative, boost);
            if (!complete) break;
        }
        const float va = w.bodies()[a].vel.x - boost, vb = w.bodies()[b].vel.x - boost;
        CHECK(std::fabs(va - 0.25f * relative) < 1e-3f * relative, "A impact velocity %g vs %g", va, 0.25f * relative);
        CHECK(std::fabs(vb - 0.75f * relative) < 1e-3f * relative, "B impact velocity %g vs %g", vb, 0.75f * relative);
        CHECK(std::fabs(va + vb - relative) < 2e-5f * relative, "impact changed total momentum");
        std::printf("  relative=%g boost=%g: outgoing relative to boost %g/%g m/s\n", relative, boost, va, vb);
    }
}
