// A large static triangle mesh as the ground (as Bullet's btBvhTriangleMeshShape): a heightfield
// of 51 200 triangles with a BVH, 150 bodies dropped on it. Nothing may fall through, everything
// must come to rest, and the step must stay cheap - the tree makes each body's contact search
// independent of the mesh size.
#include "TestRunner.h"

#include "core/Mesh.h"
#include "core/Probe.h"
#include "rigid/RigidWorld.h"
#include "spatial/BVH.h"

#include <chrono>
#include <cmath>

using namespace rf;

void testTerrain() {
    auto height = [](float x, float z) {
        return 0.4f * std::sin(0.7f * x) * std::cos(0.5f * z) + 0.25f * std::sin(1.3f * z + 0.4f * x) + 0.3f * (x * x + z * z) / 36.0f;
    };
    TriMesh mesh = primitives::heightfield(12.0f, 12.0f, 160, 160, height);
    MeshBVH bvh;
    bvh.build(mesh);
    std::printf("  terrain: %zu triangles\n", mesh.triangles.size());

    RigidWorld w;
    w.setDomain(AABB({-6, -2, -6}, {6, 8, 6}));
    w.setStaticMesh(&bvh);
    uint32_t seed = 11;
    auto rnd = [&] { seed = seed * 1664525u + 1013904223u; return (seed >> 8) * (1.0f / 16777216.0f); };
    const int n = 150;
    for (int i = 0; i < n; ++i) {
        const Vector3 p(-5.0f + 10.0f * rnd(), 3.0f + 3.0f * rnd(), -5.0f + 10.0f * rnd());
        if (i % 2 == 0) w.addBox(p, Vector3(0.06f + 0.1f * rnd(), 0.06f + 0.1f * rnd(), 0.06f + 0.1f * rnd()),
                                 Quaternion::fromAxisAngle({rnd(), rnd(), rnd() + 0.1f}, 6.28f * rnd()), 600.0f, Vector3(1));
        else w.addSphere(p, 0.06f + 0.1f * rnd(), 600.0f, Vector3(1));
    }
    const float dt = 1.0f / 60 / w.params.substeps;
    float worstMs = 0, totalMs = 0;
    int steps = 0;
    for (int f = 0; f < 60 * 8; ++f) { // 8 s: the last ones roll to the valley floor and sleep
        auto t0 = std::chrono::steady_clock::now();
        for (int k = 0; k < w.params.substeps; ++k) w.step(dt);
        const float ms = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
        worstMs = std::max(worstMs, ms);
        totalMs += ms;
        ++steps;
    }
    // Nothing fell through: every body's centre is above the surface under it (a body resting on a
    // slope sits at least its inradius above it; 2 cm covers the contact slop).
    int below = 0, asleep = 0;
    float maxSpeed = 0;
    for (const RigidBody& b : w.bodies()) {
        if (b.pos.y < height(b.pos.x, b.pos.z) - 0.02f) ++below;
        if (b.sleeping) ++asleep;
        maxSpeed = std::max(maxSpeed, length(b.vel));
    }
    std::printf("  terrain: %d bodies, %d below the surface, %d asleep after 8 s, max speed %.3f m/s; %.2f ms/frame (worst %.2f)\n",
                n, below, asleep, maxSpeed, totalMs / steps, worstMs);
    CHECK(below == 0, "%d bodies fell through the terrain", below);
    CHECK(maxSpeed < 0.3f, "bodies still moving at %.3f m/s after 8 s", maxSpeed);
    CHECK(asleep >= n * 9 / 10, "only %d of %d bodies asleep", asleep, n);
    CHECK(totalMs / steps < 40.0f, "a frame of 150 bodies on 51 200 triangles costs %.1f ms", totalMs / steps);
}
