// The scene graph of the editor: the text file comes back exactly, the magnet formula matches the
// textbook numbers, two free magnets pull together keeping the total momentum, and a scene with
// every kind of entity runs.
#include "TestRunner.h"

#include <clocale>

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
    floor.collider.enabled = true; // the editor adds a collider with the rigid role
    floor.rigid.fixed = true;
    Entity magnet;
    magnet.name = "Magnet \"A\""; // quotes in a name become apostrophes
    magnet.size = {0.1f, 0.05f, 0.05f};
    magnet.position = {0.1f, 0.5f, -0.2f};
    magnet.rotationDeg = {10, 20, 30};
    magnet.rigid.enabled = true;
    magnet.collider.enabled = true; // the editor adds a collider with the rigid role
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
    cone.collider.enabled = true; // the editor adds a collider with the rigid role
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
    smoker.collider.enabled = true; // the editor adds a collider with the rigid role
    smoker.emitter.enabled = true;
    smoker.emitter.smoke = 2.5f;
    smoker.emitter.temperature = 50;
    smoker.emitter.velocity = {0, 0.5f, 0};
    Entity ghost; // hidden: takes no part
    ghost.name = "Ghost";
    ghost.position = {0, 1.5f, 0};
    ghost.rigid.enabled = true;
    ghost.collider.enabled = true; // the editor adds a collider with the rigid role
    ghost.visible = false;
    g.entities = {floor, magnet, jelly, pool, stove, cone, curtain, smoker, ghost};
    for (size_t i = 0; i < g.entities.size(); ++i) g.entities[i].id = uint32_t(100 + 7 * i);
    g.entities[0].locked = true;
    g.entities[8].id = 16777217; // above 2^24: must stay exact (not through a float)
    return g;
}

} // namespace

// A program may switch the C library to its user's locale (Qt does on Linux): in Russian or German
// the decimal point is a comma, and printf / strtof follow it. The scene file must not: the same
// text is written, and a file written elsewhere still loads. Skipped where no such locale is installed.
void testSceneGraphLocale() {
    const std::string text = everyRoleGraph().save();
    const char* chosen = nullptr;
    for (const char* name : {"ru_RU.UTF-8", "ru_UA.UTF-8", "de_DE.UTF-8", "fr_FR.UTF-8", "Russian_Russia.1251", "German_Germany.1252", "de-DE"})
        if (std::setlocale(LC_NUMERIC, name) && std::localeconv()->decimal_point[0] == ',') {
            chosen = name;
            break;
        }
    if (!chosen) {
        std::setlocale(LC_NUMERIC, "C");
        std::printf("  no locale with a decimal comma installed here: skipped\n");
        return;
    }
    const std::string there = everyRoleGraph().save();
    SceneGraph back;
    std::string error;
    const bool loaded = back.load(text, error);
    std::setlocale(LC_NUMERIC, "C");
    std::printf("  C library in %s: the file %s, it %s\n", chosen, there == text ? "is the same" : "CHANGED",
                loaded ? "loads" : "does not load");
    CHECK(there == text, "the scene file depends on the locale %s:\n%s", chosen, there.c_str());
    CHECK(loaded && back.save() == text, "a file does not load in the locale %s: %s", chosen, error.c_str());
}

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
        e.collider.enabled = true; // the editor adds a collider with the rigid role
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

// A new scene of the editor: the world's 4 x 3 x 4 m box with its floor, a fixed plane.
static SceneGraph floorOnlyGraph() {
    SceneGraph g;
    Entity floor;
    floor.name = "Пол";
    floor.shape = ShapeKind::Plane;
    floor.size = {4, 0.02f, 4};
    floor.position = {0, 0.01f, 0};
    floor.rigid.enabled = true;
    floor.rigid.fixed = true;
    floor.collider.enabled = true;
    g.entities.push_back(floor);
    return g;
}

static Entity ball(const Vector3& position) {
    Entity s;
    s.name = "Сфера";
    s.shape = ShapeKind::Sphere;
    s.size = Vector3(0.2f);
    s.position = position;
    s.rigid.enabled = true; // 500 kg/m^3, friction 0.5, restitution 0.2: the editor's defaults
    s.collider.enabled = true;
    return s;
}

static Entity jelly(const Vector3& position) {
    Entity s;
    s.name = "Желе";
    s.shape = ShapeKind::Box;
    s.size = Vector3(0.3f);
    s.position = position;
    s.soft.enabled = true; // 150 kg/m^3, stiffness 0.3: the editor's defaults
    return s;
}

// A soft body's centre of mass and its speed, and the particles' kinetic and potential energy.
static void softMotion(const Simulation& sim, int body, Vector3& centre, Vector3& velocity) {
    const ParticleSystem& ps = sim.particles;
    double m = 0;
    Vector3 x(0.0f), v(0.0f);
    for (int i : ps.softBodies()[size_t(body)].particles) {
        const float mi = ps.invMasses()[size_t(i)] > 0 ? 1.0f / ps.invMasses()[size_t(i)] : 0.0f;
        x += ps.positions()[size_t(i)] * mi;
        v += ps.velocities()[size_t(i)] * mi;
        m += mi;
    }
    centre = x * float(1.0 / m);
    velocity = v * float(1.0 / m);
}

static double particleEnergy(const Simulation& sim) {
    const ParticleSystem& ps = sim.particles;
    const float g = -sim.particles.params.gravity.y;
    double e = 0;
    for (size_t i = 0; i < ps.size(); ++i) {
        if (ps.invMasses()[i] == 0) continue;
        const double m = 1.0 / ps.invMasses()[i];
        e += 0.5 * m * length2(ps.velocities()[i]) + m * g * ps.positions()[i].y;
    }
    return e;
}

// Two soft bodies made inside each other are pushed apart, not thrown apart (found by the user:
// "soft bodies that went into each other began to fly"). The overlap is removed from where the
// particles start and where they go alike (pre-stabilization, Macklin et al. 2014, sec. 4.4), so
// it never becomes a speed. (1) Weightless, a third of each jelly inside the other: after 2 s
// neither moves faster than 5 cm/s. (2) A jelly 10 cm inside one lying on the floor: it rises by
// no more than the overlap, and the energy of the particles never exceeds its start by more than
// the weight lifted out of the overlap.
void testGraphSoftOverlapNoFlight() {
    SceneGraph weightless;
    weightless.world.gravity = Vector3(0.0f);
    weightless.entities.push_back(jelly({-0.1f, 1.5f, 0}));
    weightless.entities.push_back(jelly({0.1f, 1.5f, 0}));
    Simulation apart;
    apart.load(std::make_unique<GraphScene>(weightless));
    float fastest = 0;
    for (int f = 0; f < 120; ++f) { // 2 s
        apart.stepFrame();
        for (int b = 0; b < 2; ++b) {
            Vector3 x, v;
            softMotion(apart, b, x, v);
            fastest = std::max(fastest, length(v));
        }
    }
    SceneGraph pile = floorOnlyGraph();
    pile.entities.push_back(jelly({0, 0.17f, 0}));
    pile.entities.push_back(jelly({0, 0.37f, 0})); // 10 cm inside the lower one
    Simulation stacked;
    stacked.load(std::make_unique<GraphScene>(pile));
    Vector3 top0, v;
    softMotion(stacked, 1, top0, v);
    const double e0 = particleEnergy(stacked);
    double mass = 0;
    for (float w : stacked.particles.invMasses()) mass += w > 0 ? 1.0 / w : 0.0;
    float highest = top0.y;
    double worstGain = 0;
    for (int f = 0; f < 120; ++f) {
        stacked.stepFrame();
        Vector3 top;
        softMotion(stacked, 1, top, v);
        highest = std::max(highest, top.y);
        worstGain = std::max(worstGain, particleEnergy(stacked) - e0);
    }
    const double allowed = 0.5 * mass * 9.81 * 0.1; // half of the particles lifted by the overlap
    std::printf("  jellies a third inside each other, weightless: fastest centre %.4f m/s in 2 s; a jelly 10 cm inside another: "
                "rose %.4f m, energy gain %.3f J (allowed %.3f J)\n", fastest, highest - top0.y, worstGain, allowed);
    CHECK(fastest < 0.05f, "overlapping jellies were thrown apart at %.3f m/s", fastest);
    CHECK(highest - top0.y < 0.11f, "the upper jelly flew up %.3f m", highest - top0.y);
    CHECK(worstGain < allowed, "the overlap made energy: +%.3f J", worstGain);
}

// The world's box is no wall for rigid bodies - Box2D, Jolt and PhysX have no world walls at all.
// Found by the user: a column of 100 spheres 28 m high in the 3 m room was pushed down through the
// box's lid and flew apart. (1) The column falls as nature has it: the gaps close from the bottom
// up, so for the first second nothing below reaches the top sphere - it falls freely,
// y = y0 - g t^2 / 2 - and no sphere is ever lower than a free fall would take it (the only things
// acting on it are gravity and supports from below). (2) A ball rolled off the floor's edge rolls on
// past where the old side wall stood (x = 2 m).
void testGraphNoInvisibleWalls() {
    SceneGraph tower = floorOnlyGraph();
    for (int i = 0; i < 100; ++i) tower.entities.push_back(ball({-0.9f, 0.12f + 0.28f * float(i), -0.9f})); // 8 cm gaps
    Simulation sim;
    sim.load(std::make_unique<GraphScene>(tower));
    const float g = 9.81f, y0Top = 0.12f + 0.28f * 99;
    const int frames = int(std::lround(1.0f / sim.frameDt)); // 1 s
    float worstBelowFreeFall = 0;
    for (int f = 1; f <= frames; ++f) {
        sim.stepFrame();
        const float t = float(f) * sim.frameDt;
        for (int i = 0; i < 100; ++i) { // body 0 is the floor, the spheres follow in their order
            const float freeFall = 0.12f + 0.28f * float(i) - 0.5f * g * t * t;
            worstBelowFreeFall = std::max(worstBelowFreeFall, freeFall - sim.rigid.bodies()[size_t(1 + i)].pos.y);
        }
    }
    const float t = float(frames) * sim.frameDt;
    const float topFall = y0Top - sim.rigid.bodies()[100].pos.y, freeFall = 0.5f * g * t * t;

    SceneGraph table = floorOnlyGraph();
    table.entities.push_back(ball({1.5f, 0.12f, 0}));
    table.entities.back().rigid.velocity = {3, 0, 0};
    Simulation rolled;
    rolled.load(std::make_unique<GraphScene>(table));
    for (int f = 0; f < 2 * frames; ++f) rolled.stepFrame(); // 2 s
    const Vector3 end = rolled.rigid.bodies()[1].pos;

    std::printf("  100 spheres, 28 m, in a 3 m room: the top one fell %.3f m in %.2f s (free fall %.3f); lowest below free fall "
                "%.4f m. A ball rolled off the floor at 3 m/s: x = %.2f m, y = %.3f m after 2 s\n",
                topFall, t, freeFall, worstBelowFreeFall, end.x, end.y);
    CHECK(std::fabs(topFall - freeFall) < 0.02f * freeFall, "the top sphere does not fall freely: %.3f m instead of %.3f", topFall, freeFall);
    CHECK(worstBelowFreeFall < 0.01f, "a sphere was pushed down, %.3f m below a free fall", worstBelowFreeFall);
    CHECK(end.x > 2.5f, "the ball stopped at an invisible wall: x = %.2f m", end.x);
    CHECK(end.y > 0.09f && end.y < 0.13f, "the ball is not on the ground: y = %.3f m", end.y);
}

// The role ties are what the editor promises: a plane made cloth, pinned at its top, hangs.
// Two cases: (1) a sheet stood up (turned 90 degrees) hangs its own length from the rod, nothing
// torn; (2) a level sheet pinned along one edge swings down on it and never hangs lower than its
// length. Both are tearable and lose no thread (the swing once tore ~200 at its bottom: the cloth
// solver's lag read as tension - fixed, see testClothSwingNoFalseTears).
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
    sheet.cloth.tearable = true;
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
                "level, pinned along one edge: swings, lowest point %.3f m below the edge, %d threads torn\n", pinnedA, hangA, tornA,
                hangB, tornB);
    CHECK(finiteA && finiteB, "NaN in the cloth");
    CHECK(std::fabs(pinnedA - 1.8f) < 1e-3f && std::fabs(pinnedB - 1.5f) < 1e-3f, "the pinned row moved");
    CHECK(std::fabs(hangA - 0.6f) < 0.02f && tornA == 0, "the standing sheet does not hang its length: %f m, %d torn", hangA, tornA);
    CHECK(hangB > 0.5f && hangB < 0.63f, "the swinging sheet must reach down about its length, never more: %f m", hangB);
    CHECK(tornB == 0, "the swinging sheet tore %d threads", tornB);
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
    box.collider.enabled = true; // the editor adds a collider with the rigid role
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
    floor.collider.enabled = true; // the editor adds a collider with the rigid role
    g.entities.push_back(floor);
    Entity model;
    model.name = "Model";
    model.shape = ShapeKind::Mesh;
    model.meshFile = file;
    model.size = {0.3f, 0.3f, 0.3f};
    model.position = {0, 1, 0};
    if (role == 0) model.rigid.enabled = model.collider.enabled = true;
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
    ball.collider.enabled = true; // the editor adds a collider with the rigid role
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

// ---------------------------------------------------------------------------
// Colliders apart from the geometry
// ---------------------------------------------------------------------------
// A floor and two bodies thrown along x at 1.5 m/s: a sphere that collides as a box (slides and
// stops) and a box that collides as a sphere (rolls on).
static SceneGraph lookVsColliderGraph() {
    SceneGraph g;
    g.world.size = {6, 2, 3};
    Entity floor;
    floor.name = "Floor";
    floor.shape = ShapeKind::Plane;
    floor.size = {6, 0.02f, 3};
    floor.rigid.enabled = floor.rigid.fixed = floor.collider.enabled = true;
    g.entities.push_back(floor);
    Entity ball;
    ball.name = "Sphere colliding as a box";
    ball.shape = ShapeKind::Sphere;
    ball.size = Vector3(0.2f);
    ball.position = {-2, 0.12f, -0.5f};
    ball.rigid.enabled = ball.collider.enabled = true;
    ball.rigid.velocity = {1.5f, 0, 0};
    ball.collider.kind = ColliderKind::Box;
    g.entities.push_back(ball);
    Entity block = ball;
    block.name = "Box colliding as a sphere";
    block.shape = ShapeKind::Box;
    block.position = {-2, 0.12f, 0.5f};
    block.collider.kind = ColliderKind::Sphere;
    g.entities.push_back(block);
    return g;
}

// The rigid body an entity made (-1: none).
static int bodyOfEntity(const GraphScene& scene, uint32_t id) {
    for (const MetaObject& m : scene.metaObjects(id))
        if (m.kind == MetaObject::Kind::RigidBody) return m.handle;
    return -1;
}

// What a body looks like and what it collides with are chosen apart: a sphere on a box collider
// slides to a stop (no rolling), a box on a sphere collider rolls on; the viewer gets the sphere's
// own mesh to draw, while the collider (for picking) is the box.
void testGraphColliderApart() {
    auto owned = std::make_unique<GraphScene>(lookVsColliderGraph());
    GraphScene* scene = owned.get();
    Simulation sim;
    sim.load(std::move(owned));
    const int ball = bodyOfEntity(*scene, scene->graph().entities[1].id), block = bodyOfEntity(*scene, scene->graph().entities[2].id);
    RenderSnapshot snap;
    sim.fillSnapshot(snap);
    const RenderSnapshot::Body& drawn = snap.bodies[size_t(ball)];
    const size_t sphereTriangles = entityLocalMesh(scene->graph().entities[1]).triangles.size();
    for (int f = 0; f < 180; ++f) sim.stepFrame();
    const float slid = sim.rigid.bodies()[ball].pos.x + 2, rolled = sim.rigid.bodies()[block].pos.x + 2;
    std::printf("  a sphere on a box collider slid %.3f m (v^2 / 2 mu g = %.3f); a box on a sphere collider rolled %.3f m; "
                "the sphere is drawn with %zu triangles, collides as shape type %d\n", slid, 1.5 * 1.5 / (2 * 0.5 * 9.81), rolled,
                drawn.mesh ? drawn.mesh->triangles.size() : 0, drawn.collisionShape ? int(drawn.collisionShape->type()) : -1);
    CHECK(ball >= 0 && block >= 0, "the two bodies were not made");
    CHECK(drawn.mesh && drawn.mesh->triangles.size() == sphereTriangles && drawn.shape == ShapeType::ConvexHull,
          "the sphere on a box collider is not drawn as the sphere");
    CHECK(drawn.collisionShape && drawn.collisionShape->type() == ShapeType::Box, "the sphere's collider is not the box");
    CHECK(slid > 0.1f && slid < 0.5f && length(sim.rigid.bodies()[ball].angVel) < 0.1f, "the box collider must slide to a stop (slid %f)", slid);
    CHECK(rolled > 3 * slid, "the sphere collider must roll farther (rolled %f, slid %f)", rolled, slid);
}

// A tall model on a capsule collider fitted to it stays on the floor (its lowest collider point on
// the floor's top); the collider fields survive save -> load -> save; a body whose collider sits off
// the object's centre is rebuilt without the object jumping.
void testGraphColliderFit() {
    static const char* path = "rf_test_tall.obj"; // a box 0.2 x 0.8 x 0.2 m
    std::FILE* f = std::fopen(path, "w");
    CHECK(f != nullptr, "cannot write the test OBJ file");
    std::fputs("v -0.1 -0.4 -0.1\nv 0.1 -0.4 -0.1\nv 0.1 0.4 -0.1\nv -0.1 0.4 -0.1\n"
               "v -0.1 -0.4 0.1\nv 0.1 -0.4 0.1\nv 0.1 0.4 0.1\nv -0.1 0.4 0.1\n"
               "f 1 3 2\nf 1 4 3\nf 5 6 7\nf 5 7 8\nf 1 2 6\nf 1 6 5\nf 4 8 7\nf 4 7 3\nf 1 5 8\nf 1 8 4\nf 2 3 7\nf 2 7 6\n", f);
    std::fclose(f);
    SceneGraph g = meshGraph(path, 0);
    Entity& model = g.entities[1];
    model.size = Vector3(0.4f); // scaled uniformly: 0.1 x 0.4 x 0.1 m
    model.collider.kind = ColliderKind::Capsule;
    Entity offset;
    offset.name = "Box with its collider off-centre";
    offset.position = {0.8f, 0.3f, 0};
    offset.rigid.enabled = offset.collider.enabled = true;
    offset.collider = {true, ColliderKind::Box, false, Vector3(0.2f), {0.3f, 0, 0}, {0, 30, 0}};
    g.entities.push_back(offset);
    auto owned = std::make_unique<GraphScene>(g);
    GraphScene* scene = owned.get();
    Simulation sim;
    sim.load(std::move(owned));
    for (int k = 0; k < 240; ++k) sim.stepFrame();
    const RigidBody& capsule = sim.rigid.bodies()[size_t(bodyOfEntity(*scene, scene->graph().entities[1].id))];
    const float lowest = capsule.worldBounds().lo.y;
    const int moved = bodyOfEntity(*scene, scene->graph().entities[2].id);
    const Vector3 before = sim.rigid.bodies()[size_t(moved)].pos;
    scene->rebuildEntity(sim, scene->graph().entities[2]);
    const Vector3 after = sim.rigid.bodies()[size_t(bodyOfEntity(*scene, scene->graph().entities[2].id))].pos;
    std::string error;
    SceneGraph back;
    const std::string text = g.save();
    const bool roundTrip = back.load(text, error) && back.save() == text && back.entities[2].collider.kind == ColliderKind::Box &&
                           length(back.entities[2].collider.offset - Vector3(0.3f, 0, 0)) < 1e-7f && back.entities[1].collider.kind == ColliderKind::Capsule;
    std::printf("  model on a fitted capsule: lowest collider point %.4f m (floor top 0.0100), speed %.4f m/s; off-centre collider "
                "body before / after a rebuild (%.4f %.4f %.4f) / (%.4f %.4f %.4f); round trip %d\n", lowest, length(capsule.vel),
                before.x, before.y, before.z, after.x, after.y, after.z, int(roundTrip));
    CHECK(capsule.shape->type() == ShapeType::Capsule && std::fabs(lowest - 0.01f) < 0.005f && length(capsule.vel) < 0.05f,
          "the model on its capsule does not rest on the floor (lowest %f)", lowest);
    CHECK(length(after - before) < 1e-3f, "the off-centre collider's body jumped by %f m in a rebuild", length(after - before));
    CHECK(roundTrip, "the collider fields do not round-trip: %s", error.c_str());
    std::remove(path);
}

// A collider without the rigid role is a static obstacle: a ball dropped on it bounces off, and the
// obstacle never moves.
void testGraphColliderOnly() {
    SceneGraph g;
    g.world.size = {2, 2, 2};
    Entity wall;
    wall.name = "Collider only";
    wall.size = {0.8f, 0.4f, 0.8f};
    wall.position = {0, 0.2f, 0};
    wall.collider.enabled = true;
    g.entities.push_back(wall);
    Entity ball;
    ball.name = "Ball";
    ball.shape = ShapeKind::Sphere;
    ball.size = Vector3(0.1f);
    ball.position = {0, 1.2f, 0};
    ball.rigid.enabled = ball.collider.enabled = true;
    ball.rigid.restitution = 0.6f;
    g.entities.push_back(ball);
    auto owned = std::make_unique<GraphScene>(g);
    GraphScene* scene = owned.get();
    Simulation sim;
    sim.load(std::move(owned));
    const int obstacle = bodyOfEntity(*scene, scene->graph().entities[0].id), b = bodyOfEntity(*scene, scene->graph().entities[1].id);
    const Vector3 start = sim.rigid.bodies()[size_t(obstacle)].pos;
    float rise = 0;
    for (int k = 0; k < 150; ++k) {
        sim.stepFrame();
        rise = std::max(rise, sim.rigid.bodies()[size_t(b)].vel.y);
    }
    const RigidBody& o = sim.rigid.bodies()[size_t(obstacle)];
    std::printf("  collider only: a static body (inverse mass %.1f) moved %.6f m; the ball bounced up at %.2f m/s, rests at y %.4f "
                "(the obstacle's top 0.40 + 0.05)\n", o.invMass, length(o.pos - start), rise, sim.rigid.bodies()[size_t(b)].pos.y);
    CHECK(obstacle >= 0 && o.invMass == 0 && length(o.pos - start) < 1e-6f, "the collider-only entity must be a static body");
    CHECK(rise > 0.5f, "the ball did not bounce off the obstacle (rose at %f m/s)", rise);
    CHECK(std::fabs(sim.rigid.bodies()[size_t(b)].pos.y - 0.45f) < 0.01f, "the ball rests at y %f, not on the obstacle", sim.rigid.bodies()[size_t(b)].pos.y);
}

// Lights and cameras: the text file keeps them exactly, a spot turned 90 degrees about x shines
// along -z (its -y axis turned), a camera's frame is orthonormal and looks along its -z, and the
// snapshot carries the visible lights (the hidden one stays out) posed in the world.
void testGraphLightsCameras() {
    SceneGraph g;
    Light sun;
    sun.id = 1, sun.name = "Солнце", sun.kind = LightKind::Sun, sun.intensity = 1.5f, sun.shadows = true;
    sun.rotationDeg = {-50, 30, 0};
    Light spot;
    spot.id = 2, spot.name = "Прожектор", spot.kind = LightKind::Spot, spot.coneDeg = 35, spot.softnessDeg = 7.5f;
    spot.position = {0, 2, 0}, spot.rotationDeg = {90, 0, 0};
    Light hidden;
    hidden.id = 3, hidden.name = "Лампа", hidden.visible = false, hidden.range = 3.25f;
    Camera cam;
    cam.id = 4, cam.name = "Камера 1", cam.fovDeg = 42, cam.nearClip = 0.1f, cam.farClip = 80, cam.active = true;
    cam.position = {1, 2, 3}, cam.rotationDeg = {-20, 35, 10};
    g.lights = {sun, spot, hidden};
    g.cameras = {cam};

    const std::string text = g.save();
    SceneGraph back;
    std::string error;
    CHECK(back.load(text, error), "lights and cameras do not load: %s", error.c_str());
    CHECK(back.save() == text, "save -> load -> save changed the lights / cameras:\n%s\n---\n%s", text.c_str(), back.save().c_str());
    CHECK(back.lights.size() == 3 && back.cameras.size() == 1 && back.lights[1].kind == LightKind::Spot &&
              back.lights[1].softnessDeg == 7.5f && back.lights[0].shadows && back.cameras[0].active && back.cameras[0].fovDeg == 42,
          "light / camera fields not restored");
    CHECK(findObject(back, 4) == &back.cameras[0] && nextId(back) == 5, "cameras must be found by id and counted in nextId");

    const Vector3 d = lightDirection(g, spot);
    Vector3 eye, forward, up;
    cameraFrame(g, cam, eye, forward, up);
    const float ortho = std::fabs(dot(forward, up)), unit = std::fabs(length(forward) - 1) + std::fabs(length(up) - 1);
    std::printf("  lights: spot turned 90 deg about x shines (%.3f %.3f %.3f); camera forward.up %.1e, |f|,|u| off by %.1e\n",
                d.x, d.y, d.z, ortho, unit);
    CHECK(length(d - Vector3(0, 0, -1)) < 1e-5f, "the spot's direction is (%f %f %f), expected (0 0 -1)", d.x, d.y, d.z);
    CHECK(ortho < 1e-5f && unit < 1e-5f && length(eye - cam.position) < 1e-6f, "the camera frame is not orthonormal");
    Camera straight;
    cameraFrame(g, straight, eye, forward, up);
    CHECK(length(forward - Vector3(0, 0, -1)) < 1e-6f && length(up - Vector3(0, 1, 0)) < 1e-6f, "an unturned camera looks along -z");

    Simulation sim;
    sim.load(std::make_unique<GraphScene>(g));
    RenderSnapshot s;
    sim.fillSnapshot(s);
    CHECK(s.lights.size() == 2 && s.lights[1].kind == int(LightKind::Spot) && length(s.lights[1].direction - d) < 1e-5f &&
              length(s.lights[1].position - spot.position) < 1e-6f,
          "the snapshot must carry the two visible lights posed in the world (got %zu)", s.lights.size());
    // A light edited while the scene runs: the next snapshot shows it, the bodies do not notice.
    auto* scene = dynamic_cast<GraphScene*>(sim.scene());
    CHECK(scene != nullptr, "the simulation must hold the graph scene");
    if (!scene) return;
    scene->setLightsAndCameras({sun}, {});
    sim.fillSnapshot(s);
    CHECK(s.lights.size() == 1 && s.lights[0].kind == int(LightKind::Sun) && scene->graph().cameras.empty(),
          "lights set between frames must reach the next snapshot (got %zu)", s.lights.size());
}
