// CCD bounds must include the rotation arc, not just the two endpoint boxes. Check both
// sampled material points (including an offset compound child) and an actual world step.
#include "TestRunner.h"
#include "Tests.h"

void testCcdSweptBounds() {
    BoxShape shape({0.5f, 0.01f, 0.02f});
    for (float angle : {0.0f, 0.1f, 0.7f, 2.9f}) {
        SweptPose sweep;
        sweep.shape = &shape;
        sweep.p0 = {0.2f, -0.7f, 1.0f}; sweep.p1 = {-0.4f, 0.3f, 0.8f};
        sweep.partT = {0.1f, 0.3f, -0.1f};
        sweep.q0 = Quaternion::fromAxisAngle(normalize(Vector3(1, -2, 3)), 0.6f);
        sweep.dTheta = normalize(Vector3(-2, 1, 3)) * angle;
        sweep.partR = Quaternion::fromAxisAngle({0, 1, 0}, 0.4f).toMatrix3x3();
        const AABB bounds = sweptBounds(sweep, 2e-6f);
        for (int step = 0; step <= 100; ++step) {
            const PosedShape p = sweep.at(float(step) / 100);
            for (int k = 0; k < 8; ++k) {
                const Vector3 corner((k & 1) ? 0.5f : -0.5f, (k & 2) ? 0.01f : -0.01f, (k & 4) ? 0.02f : -0.02f);
                CHECK(bounds.contains(p.p + p.R * corner), "rotation arc escaped CCD bounds, angle %.2f step %d", angle, step);
            }
        }
    }
    RigidWorld world;
    world.params.gravity = Vector3(0);
    world.params.collideWithDomain = false;
    world.params.sleeping = false;
    world.params.linearDamping = world.params.angularDamping = 0;
    world.addBox({0, 0.4f, 0}, {0.02f, 0.02f, 0.02f}, Quaternion(), 0, Vector3(1));
    const int plate = world.addBox(Vector3(0), {0.5f, 0.01f, 0.01f}, Quaternion(), 800, Vector3(1));
    world.bodies()[size_t(plate)].angVel = {0, 0, 0.95f * kPi * 60};
    world.step(1.0f / 60);
    const float turn = length(world.bodies()[size_t(plate)].rot.log());
    CHECK(world.ccdHits() > 0, "171-degree sweep missed a post outside both endpoint boxes");
    CHECK(turn > 1.3f && turn < 1.55f, "plate must stop just before crossing the post, angle %.6f", turn);
    CHECK(maxOverlap(world) < 1e-4f, "CCD stopped the plate inside the post");
}
