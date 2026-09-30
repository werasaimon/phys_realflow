// Plane contacts must retain the actual supporting region of a compound, including an extreme
// vertex absent from its most aligned face. These fixtures fail before the wall-patch repair:
// a symmetric four-foot support tips on one point, and a tilted cup loses 0.476 mm of true depth.
#include "TestRunner.h"
#include "Tests.h"

namespace {

std::shared_ptr<const CompoundShape> shallowTable() {
    std::vector<TriMesh> parts;
    auto box = [&](const Vector3& half, const Vector3& at) {
        TriMesh mesh = primitives::box(half);
        for (Vector3& p : mesh.positions) p += at;
        parts.push_back(std::move(mesh));
    };
    // The wide top's underside is 6 mm above the feet: inside the 10 mm speculative margin,
    // outside the 4 mm touching zone. It cannot support the body at this pose.
    box({0.15f, 0.01f, 0.15f}, {0, 0.016f, 0});
    for (int i = 0; i < 4; ++i)
        box({0.01f, 0.005f, 0.01f}, {i & 1 ? 0.05f : -0.05f, 0.005f, i & 2 ? 0.05f : -0.05f});
    return std::make_shared<CompoundShape>(parts, parts.front());
}

} // namespace

void testWallSupportRegion() {
    const auto shape = shallowTable();
    RigidWorld world;
    world.setDomain(AABB({-2, 0, -2}, {2, 2, 2}));
    world.params.sleeping = false;
    const int id = world.addCompound(shape, shape->centerOfMass() + Vector3(0, -0.0001f, 0), Quaternion(), 1000, {1, 1, 1});
    constexpr float dt = 1.0f / 600;
    world.step(dt);
    int loaded = 0;
    float impulse = 0;
    for (const auto& c : world.debugContacts()) {
        if (c.a != id || c.b != -3) continue;
        loaded += c.impulse > 1e-7f;
        impulse += c.impulse;
    }
    const RigidBody& body = world.bodies()[size_t(id)];
    const float ratio = impulse / (body.mass * 9.81f * dt);
    std::printf("  four feet: %d loaded points, support / (m g dt) %.6f, spin %.6f rad/s\n", loaded, ratio, length(body.angVel));
    CHECK(loaded >= 3, "lost the load-bearing footprint: %d points", loaded);
    CHECK(std::fabs(ratio - 1) < 0.002f, "the support does not carry the weight: %.6f", ratio);
    CHECK(length(body.angVel) < 1e-3f && length(body.vel) < 1e-4f, "symmetric supported body acquired motion");
}

void testWallExtremeWitness() {
    RigidWorld world;
    world.params.gravity = {0, 0, 0};
    world.setDomain(AABB({-2, 0, -2}, {2, 2, 2}));
    // Captured principal-frame orientation of cup 19 in the mixed sandbox. Geometry alone,
    // without gravity or penetration recovery, must agree with the support-map plane distance.
    const Quaternion q{-0.7252967f, 0.647397f, -0.15878864f, 0.1720703f};
    const int id = world.addBody(cupShape(), {0, 0.1f, 0}, q.normalized(), 1000, {1, 1, 1});
    RigidBody& body = world.bodies()[size_t(id)];
    body.pos.y -= body.posed().support({0, -1, 0}).y + 0.001f;
    const float exactDepth = -body.posed().support({0, -1, 0}).y;
    world.step(1.0f / 600);
    const float measured = world.deepestPenetration();
    std::printf("  tilted cup / plane: support depth %.6f mm, contact depth %.6f mm\n", 1000 * exactDepth, 1000 * measured);
    CHECK(std::fabs(exactDepth - measured) < 5e-6f, "the extreme vertex is absent: %.6f vs %.6f m", exactDepth, measured);
}
