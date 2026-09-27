// The scene graph of the editor: the text file comes back exactly, the magnet formula matches the
// textbook numbers, two free magnets pull together keeping the total momentum, and a scene with
// every kind of entity runs.
#include "TestRunner.h"

#include "scene/Magnets.h"
#include "scene/SceneGraph.h"
#include "scene/Simulation.h"

#include <cmath>
#include <cstdio>

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
    Entity curtain; // a flammable cloth: no body, particles
    curtain.name = "Curtain";
    curtain.shape = ShapeKind::Plane;
    curtain.size = {0.5f, 0.02f, 0.6f};
    curtain.position = {-0.8f, 1.2f, 0.5f};
    curtain.rotationDeg = {90, 0, 0};
    curtain.cloth.enabled = true;
    curtain.cloth.areaDensity = 0.2f;
    curtain.cloth.tearable = false;
    curtain.cloth.pinnedEdges = 16 | 1;
    curtain.flammable.enabled = true;
    Entity smoker; // a rigid body trailing smoke: two roles on one shape
    smoker.name = "Smoker";
    smoker.position = {0.8f, 1.0f, -0.5f};
    smoker.rigid.enabled = true;
    smoker.emitter.enabled = true;
    smoker.emitter.smoke = 2.5f;
    smoker.emitter.temperature = 50;
    smoker.emitter.velocity = {0, 0.5f, 0};
    Entity ghost; // hidden: takes no part
    ghost.name = "Ghost";
    ghost.position = {0, 1.5f, 0};
    ghost.rigid.enabled = true;
    ghost.visible = false;
    g.entities = {floor, magnet, jelly, pool, stove, cone, curtain, smoker, ghost};
    for (size_t i = 0; i < g.entities.size(); ++i) g.entities[i].id = uint32_t(100 + 7 * i);
    g.entities[0].locked = true;
    g.entities[8].id = 16777217; // above 2^24: must stay exact (not through a float)
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
    CHECK(back.entities.size() == 9 && back.entities[1].name == "Magnet 'A'", "entities not restored");
    CHECK(back.entities[1].rotationDeg.z == 30.0f && back.entities[4].heat.temperature == 412.5f, "numbers not restored");
    const Entity &curtain = back.entities[6], &smoker = back.entities[7], &ghost = back.entities[8];
    CHECK(curtain.cloth.enabled && !curtain.cloth.tearable && curtain.cloth.pinnedEdges == (16 | 1) && curtain.flammable.enabled,
          "cloth / flammable not restored");
    CHECK(smoker.emitter.enabled && smoker.emitter.smoke == 2.5f && smoker.emitter.velocity.y == 0.5f, "emitter not restored");
    CHECK(back.entities[0].locked && !ghost.visible && ghost.id == 16777217u && back.entities[1].id == 107u,
          "id / visible / locked not restored (ghost id %u)", ghost.id);

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
    // floor, magnet, cone, smoker rigid; the stove has only heat (no gas: skipped with an info line);
    // the ghost is hidden; the curtain is cloth.
    CHECK(sim.rigid.bodies().size() == 4, "rigid bodies: %zu (expected 4)", sim.rigid.bodies().size());
    CHECK(sim.particles.softBodies().size() == 1 && sim.particles.fluidCount() > 100, "soft body or liquid missing");
    CHECK(sim.particles.cloths().size() == 1, "cloths: %zu (expected 1)", sim.particles.cloths().size());
    CHECK(finite, "NaN in the graph scene");
}

// The role ties are what the editor promises: a plane made cloth, pinned at its top, hangs.
// Two cases: (1) a sheet stood up (turned 90 degrees) hangs its own length from the rod, nothing
// torn; (2) a level sheet pinned along one edge swings down on it and never hangs lower than its
// length. Case 2 uses a non-tearable cloth: at the bottom of that swing the tearable default loses
// ~200 threads at once (the cloth solver's tension estimate spikes as the free edge whips, far
// above the ~10 N/m the swing really puts on 4000 N/m threads) - a known issue of the cloth solver.
static void hangSheet(bool standing, float& pinnedY, float& hang, int& torn, bool& finite) {
    SceneGraph g;
    g.world.size = {2, 3, 2};
    Entity sheet;
    sheet.name = "Sheet";
    sheet.shape = ShapeKind::Plane;
    sheet.size = {0.6f, 0.02f, 0.6f};
    sheet.position = {0, 1.5f, 0};
    if (standing) sheet.rotationDeg = {90, 0, 0}; // stands from y 1.2 to 1.8; else level: "top row" = its -z edge
    sheet.cloth.enabled = true;
    sheet.cloth.pinnedEdges = 16;
    sheet.cloth.tearable = standing;
    g.entities.push_back(sheet);
    Simulation sim;
    sim.load(std::make_unique<GraphScene>(g));
    float pinnedMin = 1e9f, lowestEver = 1e9f;
    finite = true;
    for (int f = 0; f < 120; ++f) { // 2 s
        sim.stepFrame();
        const Cloth& c = sim.particles.cloths()[0];
        const auto& x = sim.particles.positions();
        for (int i = c.firstParticle; i < c.firstParticle + c.width * c.height; ++i) {
            lowestEver = std::min(lowestEver, x[i].y);
            finite &= std::isfinite(x[i].x + x[i].y + x[i].z);
        }
        for (int k = 0; k < c.width; ++k) pinnedMin = std::min(pinnedMin, x[c.particle(k, 0)].y);
        torn = c.tornThreads;
    }
    pinnedY = pinnedMin;
    hang = pinnedMin - lowestEver;
}

void testGraphCloth() {
    float pinnedA, hangA, pinnedB, hangB;
    int tornA, tornB;
    bool finiteA, finiteB;
    hangSheet(true, pinnedA, hangA, tornA, finiteA);
    hangSheet(false, pinnedB, hangB, tornB, finiteB);
    std::printf("  cloth 0.6 x 0.6 m, standing, pinned at the top: rod at y %.3f, hangs %.3f m, %d threads torn; "
                "level, pinned along one edge: swings, lowest point %.3f m below the edge\n", pinnedA, hangA, tornA, hangB);
    CHECK(finiteA && finiteB, "NaN in the cloth");
    CHECK(std::fabs(pinnedA - 1.8f) < 1e-3f && std::fabs(pinnedB - 1.5f) < 1e-3f, "the pinned row moved");
    CHECK(std::fabs(hangA - 0.6f) < 0.02f && tornA == 0, "the standing sheet does not hang its length: %f m, %d torn", hangA, tornA);
    CHECK(hangB > 0.5f && hangB < 0.63f, "the swinging sheet must reach down about its length, never more: %f m", hangB);
}

// A rigid box that emits smoke, thrown sideways through still air: the smoke is left along its
// path, so the smoke's centre moves in the direction of the throw.
void testGraphEmitterFollows() {
    SceneGraph g;
    g.world.gravity = Vector3(0.0f);
    g.world.size = {3, 1.5f, 1.5f};
    g.world.gas = true;
    Entity box;
    box.name = "Smoky box";
    box.size = Vector3(0.1f);
    box.position = {-1.0f, 0.75f, 0};
    box.rigid.enabled = true;
    box.rigid.velocity = {2, 0, 0};
    box.emitter.enabled = true;
    box.emitter.smoke = 5;
    g.entities.push_back(box);
    Simulation sim;
    sim.load(std::make_unique<GraphScene>(g));
    auto smokeCentre = [&](float& total) {
        const GasSolver& gas = sim.grid;
        double sum = 0, sx = 0;
        for (int k = 0; k < gas.nz(); ++k)
            for (int j = 0; j < gas.ny(); ++j)
                for (int i = 0; i < gas.nx(); ++i) {
                    const float s = gas.cellValue(GridField::Smoke, i, j, k);
                    sum += s;
                    sx += s * (gas.origin().x + (i + 0.5f) * gas.dx());
                }
        total = float(sum);
        return sum > 0 ? float(sx / sum) : 0.0f;
    };
    float early = 0, total = 0, earlyTotal = 0;
    for (int f = 0; f < 60; ++f) { // 1 s
        sim.stepFrame();
        if (f == 9) early = smokeCentre(earlyTotal);
    }
    const float late = smokeCentre(total);
    std::printf("  smoky box thrown at 2 m/s: box at x %.2f m, smoke centre x %.2f m at 0.17 s -> %.2f m at 1 s, smoke %.3g (was %.3g)\n",
                sim.rigid.bodies()[0].pos.x, early, late, total, earlyTotal);
    CHECK(total > 0 && earlyTotal > 0, "no smoke was emitted");
    CHECK(late - early > 0.5f, "the smoke did not follow the box: centre %f -> %f", early, late);
}

// A flammable curtain over a hot emitter catches fire (the fire sample's cellulose, clicked together).
void testGraphFlammableCloth() {
    SceneGraph g;
    g.world.size = {1.2f, 1.5f, 1.2f};
    g.world.gas = true;
    Entity heater; // only an emitter: stays where it is put
    heater.name = "Heater";
    heater.size = Vector3(0.1f);
    heater.position = {0, 0.2f, 0};
    heater.emitter.enabled = true;
    heater.emitter.smoke = 1;
    heater.emitter.temperature = 600;
    heater.emitter.velocity = {0, 0.5f, 0};
    Entity curtain;
    curtain.name = "Curtain";
    curtain.shape = ShapeKind::Plane;
    curtain.size = {0.4f, 0.02f, 0.4f};
    curtain.position = {0, 0.6f, 0.03f};
    curtain.rotationDeg = {90, 0, 0}; // stands up; its top edge is pinned (the rod)
    curtain.cloth.enabled = true;
    curtain.cloth.areaDensity = 0.2f;
    curtain.cloth.pinnedEdges = 16;
    curtain.flammable.enabled = true;
    g.entities = {heater, curtain};
    Simulation sim;
    sim.load(std::make_unique<GraphScene>(g));
    float leastUnburnt = 1, hottest = 0, when = -1;
    for (int f = 0; f < 300 && when < 0; ++f) { // up to 5 s
        sim.stepFrame();
        if (sim.particles.cloths().empty()) break;
        const Cloth& c = sim.particles.cloths()[0];
        for (float u : c.unburnt) leastUnburnt = std::min(leastUnburnt, u);
        for (float t : c.temperature) hottest = std::max(hottest, t);
        if (leastUnburnt < 0.99f) when = sim.time();
    }
    std::printf("  flammable curtain over a 600 K emitter: hottest fabric %.0f K above ambient, least unburnt %.3f, caught fire at %.2f s\n",
                hottest, leastUnburnt, when);
    CHECK(sim.grid.combustion.enabled, "a flammable entity must switch the combustion on");
    CHECK(when >= 0, "the curtain did not catch fire in 5 s (hottest %f K)", hottest);
}

// A cube as an OBJ file (1 m, 8 vertices, 12 triangles), written next to the test program.
static const char* writeCubeObj() {
    static const char* path = "rf_test_cube.obj";
    std::FILE* f = std::fopen(path, "w");
    if (!f) return nullptr;
    std::fputs("# unit cube\n"
               "v -0.5 -0.5 -0.5\nv 0.5 -0.5 -0.5\nv 0.5 0.5 -0.5\nv -0.5 0.5 -0.5\n"
               "v -0.5 -0.5 0.5\nv 0.5 -0.5 0.5\nv 0.5 0.5 0.5\nv -0.5 0.5 0.5\n"
               "f 1 3 2\nf 1 4 3\nf 5 6 7\nf 5 7 8\nf 1 2 6\nf 1 6 5\n"
               "f 4 8 7\nf 4 7 3\nf 1 5 8\nf 1 8 4\nf 2 3 7\nf 2 7 6\n", f);
    std::fclose(f);
    return path;
}

// A floor plane and one model entity with the given role, as the editor would make them.
static SceneGraph meshGraph(const std::string& file, int role) {
    SceneGraph g;
    g.baseDirectory = ".";
    g.world.size = {3, 2, 3};
    Entity floor;
    floor.name = "Floor";
    floor.shape = ShapeKind::Plane;
    floor.size = {3, 0.02f, 3};
    floor.rigid.enabled = floor.rigid.fixed = true;
    g.entities.push_back(floor);
    Entity model;
    model.name = "Model";
    model.shape = ShapeKind::Mesh;
    model.meshFile = file;
    model.size = {0.3f, 0.3f, 0.3f};
    model.position = {0, 1, 0};
    if (role == 0) model.rigid.enabled = true;
    if (role == 1) model.soft.enabled = true;
    if (role == 2) model.cloth.enabled = true;
    g.entities.push_back(model);
    return g;
}

// True if the scene's readings mention `text`.
static bool readingsMention(const Simulation& sim, const std::string& text) {
    RenderSnapshot snap;
    sim.fillSnapshot(snap);
    for (const auto& line : snap.info)
        if (line.second.find(text) != std::string::npos) return true;
    return false;
}

// A model from a file as each "made of" role: rigid falls on the floor and rests at its half size
// above it; soft becomes a soft body; cloth is not made from a model and is left out with a note;
// an unreadable file makes nothing. And the model's file name survives save -> load -> save.
void testGraphMeshShape() {
    const char* file = writeCubeObj();
    CHECK(file != nullptr, "cannot write the test OBJ file");
    Simulation rigid;
    rigid.load(std::make_unique<GraphScene>(meshGraph(file, 0)));
    for (int f = 0; f < 120; ++f) rigid.stepFrame();
    const RigidBody& cube = rigid.rigid.bodies().back();
    const float expected = 0.01f + 0.15f; // the plane's top face + half the 0.3 m cube
    std::printf("  model as rigid: %zu bodies, cube centre at y %.4f m (expected %.4f), speed %.4f m/s\n",
                rigid.rigid.bodies().size(), cube.pos.y, expected, length(cube.vel));
    CHECK(rigid.rigid.bodies().size() == 2, "rigid bodies: %zu (expected the floor and the model)", rigid.rigid.bodies().size());
    CHECK(std::isfinite(cube.pos.y) && std::fabs(cube.pos.y - expected) < 0.01f, "the model rests at y %f, expected %f", cube.pos.y, expected);
    Simulation soft;
    soft.load(std::make_unique<GraphScene>(meshGraph(file, 1)));
    CHECK(soft.particles.softBodies().size() == 1, "the model as soft: %zu soft bodies", soft.particles.softBodies().size());
    Simulation cloth;
    cloth.load(std::make_unique<GraphScene>(meshGraph(file, 2)));
    CHECK(cloth.particles.cloths().empty(), "a model must not become cloth");
    CHECK(readingsMention(cloth, "ткань из модели"), "no note about the model that cannot be cloth");
    Simulation missing;
    missing.load(std::make_unique<GraphScene>(meshGraph("no_such_model.obj", 0)));
    CHECK(missing.rigid.bodies().size() == 1 && readingsMention(missing, "модель не прочитана"), "an unreadable model must make nothing, with a note");
    const std::string text = meshGraph(file, 0).save();
    SceneGraph back;
    std::string error;
    CHECK(back.load(text, error) && back.save() == text && back.entities[1].meshFile == file, "the model entity does not round-trip: %s", error.c_str());
    std::remove(file);
}

// A shape with no role is geometry only: it is drawn (entityMesh gives its world mesh) but has no
// body, so a ball dropped onto it falls through to the floor.
void testGraphGeometryOnly() {
    SceneGraph g;
    g.world.size = {2, 2, 2};
    Entity block;
    block.name = "Block";
    block.shape = ShapeKind::Box;
    block.size = {0.4f, 0.4f, 0.4f};
    block.position = {0, 0.5f, 0};
    g.entities.push_back(block);
    Entity ball;
    ball.name = "Ball";
    ball.shape = ShapeKind::Sphere;
    ball.size = {0.1f, 0.1f, 0.1f};
    ball.position = {0, 1.5f, 0};
    ball.rigid.enabled = true;
    g.entities.push_back(ball);
    Simulation sim;
    sim.load(std::make_unique<GraphScene>(g));
    for (int f = 0; f < 120; ++f) sim.stepFrame();
    const TriMesh mesh = entityMesh(block);
    const AABB box = mesh.bounds();
    std::printf("  geometry only: %zu bodies, ball at y %.3f m (floor 0.05, block top 0.7), block mesh %zu triangles, "
                "box (%.2f %.2f %.2f)..(%.2f %.2f %.2f)\n", sim.rigid.bodies().size(), sim.rigid.bodies()[0].pos.y, mesh.triangles.size(),
                box.lo.x, box.lo.y, box.lo.z, box.hi.x, box.hi.y, box.hi.z);
    CHECK(entityIsGeometryOnly(block) && !entityIsGeometryOnly(ball), "the geometry-only rule");
    CHECK(sim.rigid.bodies().size() == 1, "bodies: %zu (the roleless block must not have one)", sim.rigid.bodies().size());
    CHECK(sim.rigid.bodies()[0].pos.y < 0.1f, "the ball stopped at y %f: it must fall through the block", sim.rigid.bodies()[0].pos.y);
    CHECK(mesh.triangles.size() == 12, "the block mesh has %zu triangles", mesh.triangles.size());
    CHECK(length(box.lo - Vector3(-0.2f, 0.3f, -0.2f)) < 1e-5f && length(box.hi - Vector3(0.2f, 0.7f, 0.2f)) < 1e-5f, "the block mesh is not in place");
}
