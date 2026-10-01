// Analytic flight/impact oracles and failure contracts for the restricted variational path.
// These test the real RigidWorld API, constant and variable caller steps, and curved CCD.
#include "TestRunner.h"
#include "Tests.h"
#include "rigid/VariationalFlight.h"
#include "rigid/ConservativeAdvancement.h"
#include <limits>

namespace {
void configure(RigidWorld& w) {
    w.params.sleeping = w.params.collideWithDomain = false;
    w.params.linearDamping = w.params.angularDamping = 0;
    w.params.rollingResistance = 0;
    w.params.iterations = 40;
}

int ball(RigidWorld& w, const Vector3& p, float mass = 1) {
    const int i = w.addSphere(p, 0.1f, 1, Vector3(1));
    auto& b = w.bodies()[i]; b.mass = mass; b.invMass = 1 / mass;
    b.invInertiaLocal = Vector3(250 / mass); b.updateInertia();
    b.restitution = 1; b.friction = b.staticFriction = 0;
    return i;
}

double energy(const RigidWorld& w) {
    double e = 0;
    for (const auto& b : w.bodies()) if (b.invMass > 0)
        for (int k = 0; k < 3; ++k) e += b.mass * (0.5 * double(b.vel[k]) * b.vel[k] - double(w.params.gravity[k]) * b.pos[k]);
    return e;
}

uint64_t state(const RigidWorld& w) {
    StateHash h;
    for (const auto& b : w.bodies()) {
        h.add(b.pos); h.add(b.vel); h.add(b.rot); h.add(b.angVel); h.add(b.force); h.add(b.torque);
        h.add(b.prevPos); h.add(b.prevRot); h.add(b.biasVel); h.add(b.biasAngVel);
    }
    return h.h;
}
} // namespace

void testVariationalFlight() {
    RigidWorld w; configure(w);
    const int id = ball(w, {2, 3, -1}, 2);
    w.bodies()[id].vel = {1, -2, 3};
    const Vector3 initial = w.bodies()[id].pos, velocity = w.bodies()[id].vel;
    const Vector3 force{2, 4, -1}, a = w.params.gravity + force * 0.5f;
    double time = 0;
    for (int i = 0; i < 1000; ++i) {
        const float h = i % 3 == 0 ? 0.003f : 0.0005f;
        w.bodies()[id].force = force;
        CHECK(w.tryVariationalStep(h), "free flight rejected"); time += h;
    }
    const auto& b = w.bodies()[id];
    for (int k = 0; k < 3; ++k) {
        const double q = initial[k] + velocity[k] * time + 0.5 * a[k] * time * time;
        CHECK(std::fabs(b.pos[k] - q) < 1e-4, "variable-step ballistic position error");
        CHECK(std::fabs(b.vel[k] - (velocity[k] + a[k] * time)) < 1e-4, "load consumed incorrectly");
    }
    CHECK(length2(b.force) == 0, "completed flight must consume the load");
    // Legendre transform checked by independent central differences in configuration.
    const Vector3 q0{1, 2, 3}, q1{1.5f, 1.75f, 3.25f};
    const double h = 0.125, eps = 1.0 / 1024;
    for (int k = 0; k < 3; ++k) {
        Vector3 lo = q0, hi = q0; lo[k] -= float(eps); hi[k] += float(eps);
        const double momentum = -(VariationalFlight::action(hi, q1, a, 2, h) - VariationalFlight::action(lo, q1, a, 2, h)) / (2 * eps);
        CHECK(std::fabs(momentum - 2 * ((q1[k] - q0[k]) / h - 0.5 * a[k] * h)) < 1e-8,
              "discrete action has the wrong initial momentum");
    }
}

void testVariationalBounce() {
    for (float step : {0.01f, 0.005f, 0.002f, 0.001f, 0.0005f, -0.01f}) {
        RigidWorld w; configure(w);
        w.addBox({0,-5,0}, {25,5,25}, Quaternion(), 0, Vector3(1));
        const int id = ball(w, {0,5,0});
        const double e0 = energy(w), g = -double(w.params.gravity.y), first = std::sqrt(2 * (5 - 0.1) / g);
        double time = 0, worst = 0, positionSq = 0; int frames = 0; size_t events = 0;
        while (time < 20) {
            float h = std::fabs(step);
            if (step < 0 && w.bodies()[id].pos.y < 0.2f) h /= 16;
            h = float(std::min(double(h), 20 - time));
            const bool done = w.tryVariationalStep(h);
            CHECK(done, "bounce step rejected: h=%g t=%g status=%d", h, time, int(w.lastStepResult().status));
            if (!done) break;
            time += h; events += w.lastStepResult().impactEvents;
            worst = std::max(worst, std::fabs(energy(w) / e0 - 1));
            const double tau = std::fmod(time - first, 2 * first);
            const double y = time < first ? 5 - 0.5 * g * time * time : 0.1 + g * first * tau - 0.5 * g * tau * tau;
            positionSq += std::pow(w.bodies()[id].pos.y - y, 2); ++frames;
            CHECK(w.bodies()[id].pos.y >= 0.1f - 1e-5f, "sphere penetrated the floor");
        }
        CHECK(time >= 20 && events >= 9, "incomplete bouncing trajectory");
        CHECK(worst < 2e-4 && std::sqrt(positionSq / frames) < 0.025, "energy or phase error too large");
        std::printf("  h=%g variable=%d: max dE/E=%g, position RMSE=%g m, events=%zu\n",
                    std::fabs(step), int(step < 0), worst, std::sqrt(positionSq / frames), events);
    }
}

void testVariationalImpacts() {
    for (float boost : {0.0f, 10.0f}) {
        RigidWorld w; configure(w); w.params.gravity = Vector3(0);
        const int a = ball(w, {-0.4f,0,0}, 2), b = ball(w, {0.4f,0,0}, 3);
        w.bodies()[a].vel.x = boost + 4; w.bodies()[b].vel.x = boost - 1;
        const double e0 = energy(w);
        CHECK(w.tryVariationalStep(0.2f), "unequal-mass collision rejected");
        CHECK(std::fabs(w.bodies()[a].vel.x - (boost - 2)) < 2e-5f && std::fabs(w.bodies()[b].vel.x - (boost + 3)) < 2e-5f,
              "wrong analytic outgoing velocities");
        CHECK(std::fabs(energy(w) / e0 - 1) < 2e-5, "elastic collision lost energy");
        CHECK(std::fabs(2 * w.bodies()[a].vel.x + 3 * w.bodies()[b].vel.x - (5 * boost + 5)) < 1e-4,
              "pair impulse changed total momentum");
    }
    RigidWorld w; configure(w); w.params.gravity = Vector3(0);
    for (int k = 0; k < 3; ++k) ball(w, {0.2f * k,0,0});
    w.bodies()[0].vel.x = 1;
    CHECK(w.tryVariationalStep(0.01f), "coupled three-body elastic impact rejected");
    CHECK(std::fabs(w.bodies()[0].vel.x + 1.0f / 3) < 1e-4 && std::fabs(w.bodies()[2].vel.x - 2.0f / 3) < 1e-4,
          "simultaneous Newton impact law not satisfied");
}

void testVariationalRejections() {
    RigidWorld w; configure(w);
    w.addBox({0,-5,0}, {25,5,25}, Quaternion(), 0, Vector3(1));
    const int id = ball(w, {0,0.1f,0});
    w.bodies()[id].force = {0,-1,0};
    uint64_t before = state(w);
    CHECK(!w.tryVariationalStep(0.01f) && w.lastStepResult().acceptedDt == 0, "persistent support silently advanced");
    CHECK(state(w) == before, "failed support changed state or consumed loads");
    w.bodies()[id].pos.y = 5;
    w.bodies()[id].friction = 0.5f; before = state(w);
    CHECK(!w.tryVariationalStep(0.01f) && w.lastStepResult().status == RigidStepStatus::Unsupported, "friction silently discarded");
    CHECK(state(w) == before, "unsupported model changed state");
    w.bodies()[id].friction = 0;
    CHECK(!w.tryVariationalStep(std::numeric_limits<float>::quiet_NaN()), "NaN time accepted");
    w.addBox({10,10,10}, Vector3(1), Quaternion(), 1, Vector3(1));
    before = state(w);
    CHECK(!w.tryVariationalStep(0.01f) && w.lastStepResult().status == RigidStepStatus::Unsupported,
          "anisotropic rigid body silently used spherical dynamics");
    CHECK(state(w) == before, "unsupported body changed state");
    RigidWorld bay; configure(bay); bay.params.gravity = Vector3(0);
    bay.params.collideWithDomain = true; bay.setDomain(AABB({-0.2f,-1,-1}, {0.2f,1,1}));
    ball(bay, Vector3(0)); bay.bodies()[0].vel = {1000,0,0};
    const uint64_t original = state(bay);
    CHECK(!bay.tryVariationalStep(1) && bay.lastStepResult().status == RigidStepStatus::EventLimit,
          "unbounded impact sequence was accepted");
    CHECK(state(bay) == original && bay.lastStepResult().acceptedDt == 0, "partial event sequence escaped rollback");
    RigidWorld spin; configure(spin); spin.params.gravity = Vector3(0);
    ball(spin, Vector3(0)); spin.bodies()[0].angVel.x = 1e30f;
    const uint64_t spinBefore = state(spin);
    CHECK(!spin.tryVariationalStep(100) && spin.lastStepResult().status == RigidStepStatus::NonFinite,
          "non-finite orientation escaped the trial guard");
    CHECK(state(spin) == spinBefore, "failed free turn did not restore orientation");
}

void testVariationalCurvedCcd() {
    SphereShape sphere(0.1f);
    SweptPose sweep; sweep.shape = &sphere; sweep.p0 = sweep.p1 = {0,1,0}; sweep.translationCurve = {0,8,0};
    const auto bounds = sweptBounds(sweep);
    CHECK(bounds.lo.y <= -1.1f, "swept bounds missed an interior extremum");
    const auto result = ConservativeAdvancement(1e-4f).againstPlane(sweep, {0,1,0}, 0);
    const double exact = 0.5 * (1 - std::sqrt(1 - 4 * 0.9 / 8));
    CHECK(result.status == ToiStatus::Impact && result.s <= exact && exact - result.s < 1e-4,
          "equal endpoint positions hid an accelerated impact: s=%g vs %g", result.s, exact);
    sweep.p0 = {0,0.101f,0}; sweep.p1 = {0,0.3f,0}; sweep.translationCurve = Vector3(0);
    CHECK(ConservativeAdvancement(0.002f,64,true).againstPlane(sweep,{0,1,0},0).status == ToiStatus::Separated,
          "a separating contact cannot leave the contact band");
    sweep.p1 = {0,-0.1f,0}; sweep.translationCurve = {0,-1,0};
    CHECK(ConservativeAdvancement(0.002f,64,true).againstPlane(sweep,{0,1,0},0).status != ToiStatus::Separated,
          "departure incorrectly ignored return to the same plane");
}
