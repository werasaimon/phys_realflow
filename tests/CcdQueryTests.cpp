// Independent analytic CCD oracles and explicit failure paths: translating spheres/planes,
// an extremely short time to impact, query budget exhaustion, and world-level propagation.
#include "TestRunner.h"
#include "Tests.h"
#include "rigid/ConservativeAdvancement.h"

#include <limits>

void testCcdQueryAnalytic() {
    SphereShape sphere(0.5f);
    ConservativeAdvancement query(1e-4f);
    for (int k = 0; k < 100; ++k) {
        const float offset = 0.008f * float(k);
        SweptPose a{&sphere, {-3, offset, 0}, {3, offset, 0}, Quaternion(), Vector3(0)};
        SweptPose b{&sphere, {2, 0, 0}, {2, 0, 0}, Quaternion(), Vector3(0)};
        const double exact = (5.0 - std::sqrt(1.0 - double(offset) * offset)) / 6.0;
        const ToiResult result = query.between(a, b);
        CHECK(result.status == ToiStatus::Impact && result.hit, "sphere case %d was not an impact", k);
        CHECK(result.s <= exact && exact - result.s < 5e-4, "sphere TOI %g vs exact %.12g", result.s, exact);
        const Vector3 delta = a.at(result.s).p - b.at(result.s).p;
        CHECK(double(length(delta)) >= 1.0, "returned sphere configuration intersects");
        const ToiResult reverse = query.between(b, a);
        CHECK(std::fabs(reverse.s - result.s) < 1e-6f, "swapping pair changes TOI");
    }
    SweptPose a{&sphere, {0, 3, 0}, {0, -3, 0}, Quaternion(), Vector3(0)};
    const ToiResult plane = query.againstPlane(a, {0, 1, 0}, 0);
    CHECK(plane.hit && plane.s <= 2.5 / 6.0 && 2.5 / 6.0 - plane.s < 1e-4, "plane analytic TOI %g", plane.s);
    CHECK(query.againstPlane(a, {0, 7, 0}, 0).s == plane.s, "plane normal scaling changes TOI");
    SphereShape tiny(0.01f);
    SweptPose fast{&tiny, {0, 0, 0}, {1e6f, 0, 0}, Quaternion(), Vector3(0)};
    SweptPose target{&tiny, {0.04f, 0, 0}, {0.04f, 0, 0}, Quaternion(), Vector3(0)};
    const ToiResult shortTime = query.between(fast, target);
    CHECK(shortTime.hit && shortTime.s > 0 && shortTime.s < 2e-8f,
          "minimum time advance skipped a thin obstacle: %.12g status %d", shortTime.s, int(shortTime.status));
    std::printf("  100 analytic sphere impacts, plane and sub-1e-6 impact fractions checked\n");
}

void testCcdQueryStatuses() {
    SphereShape sphere(0.5f);
    SweptPose a{&sphere, {-3, 0, 0}, {3, 0, 0}, Quaternion(), Vector3(0)};
    SweptPose b{&sphere, {2, 0, 0}, {2, 0, 0}, Quaternion(), Vector3(0)};
    for (int budget : {0, 1}) {
        const ToiResult r = ConservativeAdvancement(1e-4f, budget).between(a, b);
        CHECK(r.status == ToiStatus::Unresolved && r.reason == ToiReason::IterationLimit && !r.hit,
              "exhausted query incorrectly reports a conclusion");
        CHECK(r.needsClamping() && r.s >= 0 && r.s < 4.0f / 6,
              "budget failure advanced beyond the analytic collision");
        CHECK((budget == 0 && r.s == 0) || (budget == 1 && r.s > 0.6f),
              "the bounded prefix was lost on exhaustion");
        const ToiResult p = ConservativeAdvancement(1e-4f, budget).againstPlane(a, {-1, 0, 0}, 0);
        CHECK(p.status == ToiStatus::Unresolved && p.needsClamping(), "plane budget failure ignored");
    }
    ConservativeAdvancement query(1e-4f);
    a.p1 = a.p0;
    CHECK(query.between(a, b).status == ToiStatus::Separated, "stationary separated pair");
    a.p0 = a.p1 = b.p0;
    const ToiResult contact = query.between(a, b);
    CHECK(contact.status == ToiStatus::InitialContact && !contact.needsClamping(), "initial overlap hidden as separation");
    CHECK(query.againstPlane(a, {-1, 0, 0}, 0).status == ToiStatus::InitialContact, "initial plane overlap");
    a.p1.x = std::numeric_limits<float>::quiet_NaN();
    CHECK(query.between(a, b).reason == ToiReason::InvalidInput, "NaN motion accepted");
    a.shape = nullptr;
    CHECK(query.between(a, b).reason == ToiReason::InvalidInput, "null shape accepted");
    a = b;
    CHECK(query.againstPlane(a, Vector3(0), 0).reason == ToiReason::InvalidInput, "zero plane normal accepted");
    CHECK(ConservativeAdvancement(0).between(a, b).reason == ToiReason::InvalidInput, "zero tolerance accepted");
    CHECK(ConservativeAdvancement(std::numeric_limits<float>::quiet_NaN()).between(a, b).reason == ToiReason::InvalidInput,
          "NaN tolerance accepted");
    a.p0 = {1000002, 0, 0}; a.p1 = {999998, 0, 0};
    const ToiResult distantPlane = query.againstPlane(a, {1, 0, 0}, 1000000);
    CHECK(distantPlane.status == ToiStatus::Unresolved && distantPlane.reason == ToiReason::NoSeparatingPlane,
          "rounding uncertainty must not masquerade as an initial contact");
    TriangleShape triangle({-2, 0, 0}, {2, 0, 0}, {0, 0.1f, 0});
    SweptPose rotating;
    rotating.shape = &triangle; rotating.dTheta = {0, 0, 1};
    CHECK(rotating.angularReach() >= 2, "moving triangle cannot use its static zero bounding radius");
}

void testCcdQuerySliding() {
    BoxShape box(Vector3(0.05f));
    for (float angle : {0.0f, 0.4f, 1.1f}) {
        const Quaternion rotation = Quaternion::fromAxisAngle({1, 2, -1}, angle);
        const Matrix3x3 r = rotation.toMatrix3x3();
        TriangleShape triangle(r * Vector3(-5, 0, -5), r * Vector3(5, 0, -5), r * Vector3(0, 0, 5));
        SweptPose a{&box, r * Vector3(0, 0.061f, 0), r * Vector3(0.8f, 0.061f, 0), rotation, Vector3(0)};
        SweptPose b; b.shape = &triangle;
        ConservativeAdvancement query(1e-4f, 2);
        // The exact gap remains 11 mm throughout this translation along the triangle face.
        const ToiResult result = query.between(a, b), reverse = query.between(b, a);
        CHECK(result.status == ToiStatus::Separated && reverse.status == ToiStatus::Separated,
              "parallel sliding exhausted the CCD budget, angle %g", angle);
        CHECK(result.s == 1 && !result.needsClamping(), "separated sliding was clamped");
        CHECK(query.againstPlane(a, r * Vector3(0, 1, 0), 0).status == ToiStatus::Separated,
              "parallel plane query must finish with the same small budget");
    }
}

static void checkSphereBoxSeparation() {
    // Captured bullet/plate pair: the supplied supporting plane proves a 2.78 mm gap.
    BoxShape plate({0.02f, 0.02f, 0.5f});
    SphereShape sphere(0.03f);
    SweptPose a, b;
    a.shape = &plate;
    a.p0 = a.p1 = {-0.999770105f, 2.19977975f, 0.599928379f};
    a.q0 = {0.955300987f, 0.00143690326f, 0.295611888f, 0.00341342925f};
    b.shape = &sphere;
    b.p0 = b.p1 = {-1.04958105f, 2.22494078f, 0.620101392f};
    const Vector3 axis = normalize(Vector3(0.814298809f, -0.161686227f, -0.55747205f));
    const float gap = dot(axis, a.at(0).support(-axis) - b.at(0).support(axis));
    CHECK(gap > 0.0027f, "reference sphere/box plane does not separate: %g", gap);
    ConservativeAdvancement query(0.002f);
    for (bool movingAway : {false, true}) {
        if (movingAway) b.p1 -= axis * 0.1f;
        CHECK(query.between(a, b).status == ToiStatus::Separated && query.between(b, a).status == ToiStatus::Separated,
              "sphere/box separating plane missed, moving away %d", int(movingAway));
    }
}

static void checkCcdBulletScene() {
    Simulation sim;
    loadSample(sim, Preset::RigidCcd);
    // This mixed scene couples spheres, boxes and hinged/spinning plates. A frozen anisotropic
    // plate used to diverge to NaN in three frames despite isolated CCD tests passing.
    for (int frame = 0; frame < 6; ++frame) {
        sim.stepFrame();
        for (const RigidBody& b : sim.rigid.bodies())
            CHECK(std::isfinite(length2(b.pos) + length2(b.vel) + length2(b.angVel)),
                  "bullet scene diverged at frame %d", frame);
    }
}

void testCcdQueryBoxSeparation() {
    // Captured before the beam regression froze. An explicit supporting plane proves a gap;
    // GJK's finite simplex alone supplied a direction with only an 8 micrometre lower bound.
    BoxShape beam({1, 0.03f, 0.03f}), cube(Vector3(0.1f));
    SweptPose a, b;
    a.shape = &beam;
    a.p0 = a.p1 = {0.000313842174f, 0.244616494f, -0.000116238785f};
    a.q0 = {0.25076285f, 0.0325911157f, 0.960484922f, -0.116295069f};
    b.shape = &cube;
    b.p0 = b.p1 = {-0.154171601f, 0.108316913f, 0.00634917663f};
    b.q0 = {0.998862445f, -4.84169868e-05f, 0.021784164f, -0.0424174853f};
    const Vector3 axis = -normalize(cross(a.at(0).R.col(0), b.at(0).R.col(0)));
    const float gap = dot(axis, a.at(0).support(-axis) - b.at(0).support(axis));
    CHECK(gap > 0.002f, "reference support plane must separate boxes, gap %g", gap);
    ConservativeAdvancement query(0.002f);
    for (bool movingAway : {false, true}) {
        if (movingAway) a.p1 += axis * 0.1f;
        const ToiResult r = query.between(a, b), reverse = query.between(b, a);
        CHECK(r.status == ToiStatus::Separated && reverse.status == ToiStatus::Separated,
              "a proved separating plane was missed, moving away %d", int(movingAway));
    }
    checkSphereBoxSeparation();
    checkCcdBulletScene();
}

void testCcdWorldUncertainty() {
    for (int environment = 0; environment < 4; ++environment) {
        RigidWorld world;
        world.params.gravity = Vector3(0);
        world.params.sleeping = false;
        world.params.linearDamping = world.params.angularDamping = 0;
        world.params.ccdMaxIterations = 0;
        world.params.ccdMaxPasses = 1;
        world.params.collideWithDomain = environment == 1;
        world.setDomain(AABB({-10, -10, -10}, {0, 10, 10}));
        MeshBVH mesh;
        if (environment == 0) world.addBox(Vector3(0), {0.01f, 1, 1}, Quaternion(), 0, Vector3(1));
        if (environment == 2) {
            mesh.build(primitives::box({0.01f, 1, 1}));
            world.setStaticMesh(&mesh);
        }
        if (environment == 3) {
            const TriMesh box = primitives::box({0.01f, 1, 1});
            world.addCompound(std::make_shared<const CompoundShape>(std::vector<TriMesh>{box}, box),
                              Vector3(0), Quaternion(), 0, Vector3(1));
        }
        const int body = world.addSphere({-2.23f, 0, 0}, 0.02f, 1000, Vector3(1));
        world.bodies()[body].vel = {300, 0, 0};
        CHECK(!world.tryStep(1.0f / 60), "zero query budget must reject the step");
        CHECK(world.bodies()[body].pos.x == -2.23f, "unresolved environment %d advanced the body", environment);
        CHECK(world.ccdDiagnostics().unresolved > 0 && world.ccdHits() == 1, "uncertainty lost by world");
        CHECK(world.lastStepResult().status == RigidStepStatus::Unresolved && world.lastStepResult().acceptedDt == 0,
              "rejection must not consume model time");
        world.params.ccd = false;
        world.step(1.0f / 60);
        CHECK(world.ccdDiagnostics().queries == 0 && world.ccdHits() == 0, "disabled CCD retains stale statistics");
        world.clear();
        CHECK(world.ccdDiagnostics().unresolved == 0, "clear retains uncertainty");
    }
}

void testCcdAllMoving() {
    for (bool allMoving : {false, true}) {
        RigidWorld world;
        world.params.gravity = Vector3(0);
        world.params.sleeping = world.params.collideWithDomain = false;
        world.params.linearDamping = world.params.angularDamping = 0;
        world.params.ccdAllMoving = allMoving;
        world.params.contactMargin = 0;
        world.addBox(Vector3(0), {0.01f, 2, 2}, Quaternion(), 0, Vector3(1));
        const int body = world.addSphere({-1.09f, 0, 0}, 1, 1000, Vector3(1));
        world.bodies()[body].vel = {12, 0, 0};
        world.step(1.0f / 60);
        if (allMoving) {
            CHECK(world.ccdDiagnostics().queries > 0 && world.bodies()[body].pos.x < -1.01f, "slow sweep escaped expanded CCD");
        } else {
            CHECK(world.ccdDiagnostics().queries == 0 && world.bodies()[body].pos.x > -1.01f, "test is not below speed threshold");
        }
    }
}
