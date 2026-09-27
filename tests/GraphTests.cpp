// The scene graph of the editor: the text file comes back exactly, the magnet formula matches the
// textbook numbers, two free magnets pull together keeping the total momentum, and a scene with
// every kind of entity runs.
#include "TestRunner.h"

#include "scene/Magnets.h"
#include "scene/SceneGraph.h"

#include <cmath>

using namespace rf;

namespace {

SceneGraph everyRoleGraph() {
    SceneGraph g;
    g.world.gravity = {0, -9.81f, 0};
    g.world.size = {3, 2, 2};
    Entity floor;
    floor.name = "Floor";
    floor.shape = ShapeKind::Plane;
    floor.size = {3, 0.02f, 2};
    floor.rigid.enabled = true;
    floor.rigid.fixed = true;
    Entity magnet;
    magnet.name = "Magnet \"A\""; // quotes in a name become apostrophes
    magnet.size = {0.1f, 0.05f, 0.05f};
    magnet.position = {0.1f, 0.5f, -0.2f};
    magnet.rotationDeg = {10, 20, 30};
    magnet.rigid.enabled = true;
    magnet.rigid.density = 7800;
    magnet.rigid.velocity = {0.25f, 0, 0};
    magnet.magnet.enabled = true;
    magnet.magnet.moment = {1, 0.5f, 0};
    Entity jelly;
    jelly.name = "Jelly";
    jelly.shape = ShapeKind::Sphere;
    jelly.size = Vector3(0.16f);
    jelly.position = {-0.5f, 0.6f, 0};
    jelly.soft.enabled = true;
    Entity pool;
    pool.name = "Pool";
    pool.shape = ShapeKind::Box;
    pool.size = {0.4f, 0.12f, 0.3f};
    pool.position = {0.6f, 0.07f, 0};
    pool.liquid.enabled = true;
    Entity stove;
    stove.name = "Stove";
    stove.shape = ShapeKind::Cylinder;
    stove.heat.enabled = true;
    stove.heat.temperature = 412.5f;
    Entity cone;
    cone.name = "Cone";
    cone.shape = ShapeKind::Cone;
    cone.position = {0, 1, 0.3f};
    cone.rigid.enabled = true;
    g.entities = {floor, magnet, jelly, pool, stove, cone};
    return g;
}

} // namespace

void testSceneGraphRoundTrip() {
    const std::string text = everyRoleGraph().save();
    SceneGraph back;
    std::string error;
    CHECK(back.load(text, error), "the saved scene does not load: %s", error.c_str());
    const std::string again = back.save();
    std::printf("  scene file: %zu bytes, %zu entities; save -> load -> save %s\n", text.size(), back.entities.size(),
                again == text ? "identical" : "DIFFERENT");
    CHECK(again == text, "save -> load -> save changed the text:\n%s\n---\n%s", text.c_str(), again.c_str());
    CHECK(back.entities.size() == 6 && back.entities[1].name == "Magnet 'A'", "entities not restored");
    CHECK(back.entities[1].rotationDeg.z == 30.0f && back.entities[4].heat.temperature == 412.5f, "numbers not restored");

    SceneGraph bad;
    const bool loaded = bad.load("world gravity 0 -9.81 0\nentity \"x\"\n  shape box size 1 1 1\n  rigid bouncy 3\nend\n", error);
    std::printf("  a bad line: %s\n", error.c_str());
    CHECK(!loaded && error.find("line 4") != std::string::npos && error.find("bouncy") != std::string::npos,
          "the error must name line 4 and the word 'bouncy': %s", error.c_str());
}

void testMagnetForce() {
    const float mu0 = 4e-7f * kPi, d = 0.1f;
    // Head to tail along x: they attract with 3 mu0 m1 m2 / (2 pi d^4) = 6 mN.
    const Vector3 onSecondCoaxial = dipoleForce({1, 0, 0}, {1, 0, 0}, {d, 0, 0});
    const float coaxial = 3 * mu0 / (2 * kPi * std::pow(d, 4.0f));
    // Side by side (both along y, apart along x): they repel with 3 mu0 m^2 / (4 pi d^4) = 3 mN.
    const Vector3 onSecondSide = dipoleForce({0, 1, 0}, {0, 1, 0}, {d, 0, 0});
    const float side = 3 * mu0 / (4 * kPi * std::pow(d, 4.0f));
    std::printf("  coaxial: F = %.6f N (theory -%.6f, pulls); side by side: F = %.6f N (theory +%.6f, pushes)\n",
                onSecondCoaxial.x, coaxial, onSecondSide.x, side);
    CHECK(std::fabs(onSecondCoaxial.x + coaxial) < 1e-5f * coaxial && std::fabs(onSecondCoaxial.y) < 1e-12f, "coaxial force wrong");
    CHECK(std::fabs(onSecondSide.x - side) < 1e-5f * side, "side-by-side force wrong");
    // Newton's third law for an arbitrary pair: the force on the first is exactly minus that on the second.
    const Vector3 m1(0.3f, -0.7f, 0.2f), m2(-0.5f, 0.1f, 0.9f), r(0.07f, 0.05f, -0.04f);
    const Vector3 f12 = dipoleForce(m1, m2, r), f21 = dipoleForce(m2, m1, -r);
    CHECK(length(f12 + f21) <= 1e-6f * length(f12), "Newton's third law: %e", double(length(f12 + f21)));
    // A needle across the field turns to it: tau = m2 x B1. Field on the axis of m1 = (1,0,0) at d:
    // B = mu0 / (2 pi d^3) along x; a needle m2 = (0,1,0) there feels tau = m2 x B = -B z.
    const Vector3 B = dipoleField({1, 0, 0}, {d, 0, 0});
    const Vector3 tau = cross(Vector3(0, 1, 0), B);
    const float Baxis = mu0 / (2 * kPi * d * d * d);
    std::printf("  field on the axis %.4e T (theory %.4e), torque on a crossed needle %.4e N m\n", B.x, Baxis, tau.z);
    CHECK(std::fabs(B.x - Baxis) < 1e-5f * Baxis && std::fabs(tau.z + Baxis) < 1e-5f * Baxis, "field or torque wrong");
}

void testMagnetsAttract() {
    SceneGraph g;
    g.world.gravity = Vector3(0.0f);
    g.world.size = {2, 3, 2};
    for (int i = 0; i < 2; ++i) {
        Entity e;
        e.name = i ? "right" : "left";
        e.size = Vector3(0.02f);
        e.position = {i ? 0.05f : -0.05f, 1.5f, 0};
        e.rigid.enabled = true;
        e.rigid.density = 7800;
        e.magnet.enabled = true;
        e.magnet.moment = {0.5f, 0, 0};
        g.entities.push_back(e);
    }
    Simulation sim;
    sim.load(std::make_unique<GraphScene>(g));
    auto momentum = [&] { return sim.rigid.bodies()[0].vel * sim.rigid.bodies()[0].mass + sim.rigid.bodies()[1].vel * sim.rigid.bodies()[1].mass; };
    auto gap = [&] { return sim.rigid.bodies()[1].pos.x - sim.rigid.bodies()[0].pos.x; };
    float worstP = 0, maxSpeed = 0, lastGaps[2] = {0, 0};
    for (int f = 0; f < 180; ++f) { // 3 s
        sim.stepFrame();
        worstP = std::max(worstP, length(momentum()));
        maxSpeed = std::max(maxSpeed, length(sim.rigid.bodies()[0].vel));
        if (f == 150) lastGaps[0] = gap();
    }
    lastGaps[1] = gap();
    // Touching: the centres 2 cm apart (the boxes' size), give or take the contact slop; still:
    // the gap did not change over the last half second (the magnets press into each other every
    // frame, so the speed right after the frame is that of the new push, not of motion).
    std::printf("  two magnets 10 cm apart: max speed %.3f m/s, final centre gap %.4f m (boxes 0.02), gap change over the "
                "last 0.5 s %.2e m, worst |P| %.2e kg m/s\n", maxSpeed, lastGaps[1], std::fabs(lastGaps[1] - lastGaps[0]), worstP);
    CHECK(maxSpeed > 0.05f, "the magnets did not move toward each other");
    CHECK(worstP < 1e-6f, "total momentum not conserved: %e", double(worstP));
    CHECK(std::fabs(lastGaps[1] - 0.02f) < 0.005f, "the magnets are not touching: gap %f", double(lastGaps[1]));
    CHECK(std::fabs(lastGaps[1] - lastGaps[0]) < 1e-3f, "the magnets are not at rest together");
}

void testGraphSceneBuilds() {
    SceneGraph g = everyRoleGraph();
    Simulation sim;
    sim.load(std::make_unique<GraphScene>(g));
    for (int f = 0; f < 60; ++f) sim.stepFrame();
    bool finite = true;
    for (const RigidBody& b : sim.rigid.bodies()) finite &= std::isfinite(b.pos.x + b.pos.y + b.pos.z + b.vel.x + b.vel.y + b.vel.z);
    for (const Vector3& p : sim.particles.positions()) finite &= std::isfinite(p.x + p.y + p.z);
    std::printf("  graph scene: %zu rigid bodies, %zu soft bodies, %zu liquid particles, %zu particles in all after 60 frames\n",
                sim.rigid.bodies().size(), sim.particles.softBodies().size(), sim.particles.fluidCount(), sim.particles.size());
    // floor, magnet, cone rigid; the stove has only heat (no gas: skipped with an info line).
    CHECK(sim.rigid.bodies().size() == 3, "rigid bodies: %zu (expected 3)", sim.rigid.bodies().size());
    CHECK(sim.particles.softBodies().size() == 1 && sim.particles.fluidCount() > 100, "soft body or liquid missing");
    CHECK(finite, "NaN in the graph scene");
}
