// Discrete actuator work: analytical kinetic energy, warm starts, braking and CCD rollback.
#include "TestRunner.h"
#include "Tests.h"

namespace {

void motorWorld(RigidWorld& w, bool measured) {
    w.params.gravity = Vector3(0);
    w.params.sleeping = w.params.collideWithDomain = false;
    w.params.linearDamping = w.params.angularDamping = 0;
    w.params.measureMotorWork = measured;
    const int wheel = w.addBox({0, 3, 0}, {0.1f, 0.1f, 0.05f}, Quaternion(), 1000, Vector3(1));
    auto& motor = w.addHingeJoint(wheel, -1, {0, 3, 0}, {0, 0, 1});
    motor.motorEnabled = true;
    motor.motorSpeed = 2;
    motor.maxMotorTorque = 0.01f;
}

uint64_t workHash(const RigidWorld& w) {
    StateHash hash;
    for (const auto& b : w.bodies()) { hash.add(b.pos); hash.add(b.rot); hash.add(b.vel); hash.add(b.angVel); }
    return hash.h;
}

} // namespace

void testMotorWork() {
    RigidWorld observed, reference;
    motorWorld(observed, true); motorWorld(reference, false);
    double previous = 0, largestError = 0;
    for (int step = 0; step < 600; ++step) {
        const float target = step < 200 ? 2.0f : step < 400 ? 0.0f : -2.0f;
        for (auto* w : {&observed, &reference}) {
            static_cast<HingeJoint&>(*w->joints().front()).motorSpeed = target;
            w->step(0.001f);
        }
        CHECK(workHash(observed) == workHash(reference), "work observation changed dynamics at %d", step);
        largestError = std::max(largestError, std::fabs(observed.kineticEnergy() - observed.motorWork()));
        CHECK(std::fabs(observed.motorWork() - previous - observed.lastStepResult().motorWork) < 1e-12,
              "accepted work is not additive");
        previous = observed.motorWork();
        if (step == 399) CHECK(std::fabs(previous) < 1e-8, "braking did not recover the rotor energy");
    }
    CHECK(largestError < 2e-8, "analytical W=delta T error %.9g J", largestError);
    CHECK(reference.motorWork() == 0, "disabled accounting reports work");
    observed.clear();
    CHECK(observed.motorWork() == 0, "clear must reset cumulative work");
    motorWorld(observed, true);
    observed.bodies()[0].angVel.z = 2;
    for (int step = 0; step < 200; ++step) {
        observed.bodies()[0].torque.z = -0.005f;
        observed.step(0.001f);
    }
    // Constant speed against a known torque: W = tau*omega*time, up to trapezoidal kick error.
    CHECK(std::fabs(observed.motorWork() - 0.005 * 2 * 0.2) < 2e-6,
          "loaded constant-speed motor work %.9g J", observed.motorWork());
}

void testMotorWorkRollback() {
    RigidWorld w;
    motorWorld(w, true);
    w.addBox({-3, 0, 0}, {0.01f, 1, 1}, Quaternion(), 0, Vector3(1));
    const int bullet = w.addSphere({-4.5f, 0, 0}, 0.02f, 1000, Vector3(1));
    w.bodies()[bullet].vel = {100, 0, 0};
    w.params.ccdMaxSubdivisions = 1;
    CHECK(!w.tryStep(0.02f), "fixture must reject its second half");
    CHECK(w.lastStepResult().motorWork == 0 && w.motorWork() == 0, "rejected trial consumed actuator work");
    CHECK(w.bodies()[0].angVel.z == 0, "rollback did not restore the motor body");
    w.params.ccdMaxSubdivisions = 16;
    CHECK(w.tryStep(0.02f) && w.lastStepResult().acceptedSubsteps > 1, "fixture must retry successfully");
    const auto& rotor = w.bodies()[0];
    const double expected = 0.5 * rotor.angVel.z * rotor.angVel.z / rotor.invInertiaLocal.z;
    CHECK(std::fabs(w.motorWork() - expected) < 2e-8, "accepted retries double-count work");
}
