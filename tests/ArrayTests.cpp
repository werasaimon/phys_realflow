// Many objects at once: hierarchy (a child of a turned group sits where the turn puts it), instances
// (ten boxes sharing one master's roles all get heavier when the master does), arrays (a line, a
// grid and a circle of copies exactly where the pattern says; fifty cubes from one array fall into
// a pile and sleep; the array grows from 50 to 80 copies while the scene plays) and a glued group
// (three boxes that fall and tumble as one body). The text file keeps all of it.
#include "TestRunner.h"

#include "scene/SceneGraph.h"
#include "scene/Simulation.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <utility>
#include <vector>

using namespace rf;

namespace {

Entity floorEntity(float width) {
    Entity e;
    e.name = "Floor";
    e.shape = ShapeKind::Plane;
    e.size = {width, 0.02f, width};
    e.rigid.enabled = e.collider.enabled = e.rigid.fixed = true;
    return e;
}

Entity boxEntity(const char* name, const Vector3& position, float size) {
    Entity e;
    e.name = name;
    e.size = Vector3(size);
    e.position = position;
    e.rigid.enabled = e.collider.enabled = true;
    return e;
}

void steps(Simulation& sim, int frames) {
    for (int i = 0; i < frames; ++i) sim.stepFrame();
}

int aliveBodies(const Simulation& sim) {
    int n = 0;
    for (int i = 0; i < int(sim.rigid.bodies().size()); ++i) n += sim.rigid.isAlive(i) ? 1 : 0;
    return n;
}

int bodyOf(const GraphScene& scene, uint32_t id) {
    for (const MetaObject& m : scene.metaObjects(id))
        if (m.kind == MetaObject::Kind::RigidBody) return m.handle;
    return -1;
}

// A graph with a turned group holding a child, a master with a chain of instances and an array.
SceneGraph hierarchyGraph() {
    SceneGraph g;
    Group grp;
    grp.id = 10;
    grp.name = "Group";
    grp.position = {1, 0, 0};
    grp.rotationDeg = {0, 90, 0};
    grp.glued = true;
    g.groups.push_back(grp);
    Entity child = boxEntity("Child", {1, 0, 0}, 0.1f);
    child.id = 1;
    child.parent = 10;
    g.entities.push_back(child);
    Entity master = boxEntity("Master", {0, 1, 0}, 0.2f);
    master.id = 2;
    master.rigid.density = 700;
    master.shape = ShapeKind::Sphere;
    g.entities.push_back(master);
    Entity a = boxEntity("A", {0, 2, 0}, 0.3f);
    a.id = 3;
    a.instanceOf = 2;
    g.entities.push_back(a);
    Entity b = boxEntity("B", {0, 3, 0}, 0.4f);
    b.id = 4;
    b.instanceOf = 3; // a chain: B -> A -> Master
    g.entities.push_back(b);
    ArrayObject arr;
    arr.id = 20;
    arr.name = "Row";
    arr.templateId = 2;
    arr.count[0] = 5;
    arr.step = {0.25f, 0, 0};
    arr.jitter = 0.01f;
    arr.seed = 7;
    g.arrays.push_back(arr);
    return g;
}

// Where every array copy's body is now, by the copy's stable id.
std::vector<std::pair<uint32_t, Vector3>> copyPositions(const Simulation& sim, const GraphScene& scene) {
    std::vector<std::pair<uint32_t, Vector3>> out;
    for (int b = 0; b < int(sim.rigid.bodies().size()); ++b) {
        const uint32_t id = scene.entityIdOfBody(b);
        if (sim.rigid.isAlive(b) && isArrayCopyId(id)) out.push_back({id, sim.rigid.bodies()[size_t(b)].pos});
    }
    return out;
}

// The farthest any of those copies is from where it was (a copy without a body counts as far).
float worstCopyMove(const Simulation& sim, const GraphScene& scene, const std::vector<std::pair<uint32_t, Vector3>>& before) {
    float worst = 0;
    for (const auto& [id, pos] : before) {
        const int b = bodyOf(scene, id);
        worst = std::max(worst, b < 0 ? 1e9f : length(sim.rigid.bodies()[size_t(b)].pos - pos));
    }
    return worst;
}

} // namespace

void testHierarchyInstancesRoundTrip() {
    const SceneGraph g = hierarchyGraph();
    const std::string text = g.save();
    SceneGraph back;
    std::string error;
    const bool loaded = back.load(text, error);
    const std::string again = back.save();
    std::printf("  file: %zu bytes, %zu entities, %zu groups, %zu arrays; round trip %s\n", text.size(), back.entities.size(),
                back.groups.size(), back.arrays.size(), text == again ? "identical" : "DIFFERENT");
    CHECK(loaded, "load failed: %s", error.c_str());
    CHECK(text == again, "save -> load -> save changed the text");
    Vector3 p;
    Quaternion q;
    worldPose(g, g.entities[0], p, q); // the child at x = 1 of a group at x = 1 turned 90 deg about y
    const Vector3 e = eulerDegrees(q);
    std::printf("  child of the turned group: world (%.4f %.4f %.4f), yaw %.2f deg (expected (1 0 -1), 90)\n", p.x, p.y, p.z, e.y);
    CHECK(length(p - Vector3(1, 0, -1)) < 1e-5f && std::fabs(e.y - 90) < 1e-3f, "wrong world pose of a child");
    SceneObject moved = g.entities[0];
    setWorldPose(g, moved, {2, 0.5f, -1}, q);
    Vector3 p2;
    Quaternion q2;
    worldPose(g, moved, p2, q2);
    CHECK(length(p2 - Vector3(2, 0.5f, -1)) < 1e-5f, "setWorldPose does not invert worldPose");
    const Entity b = resolveInstance(g, g.entities[3]);
    std::printf("  instance chain B -> A -> Master: shape %s, density %.0f, own position y %.1f, own id %u\n",
                b.shape == ShapeKind::Sphere ? "sphere" : "other", b.rigid.density, b.position.y, b.id);
    CHECK(b.shape == ShapeKind::Sphere && b.rigid.density == 700 && b.position.y == 3 && b.id == 4, "instance not resolved to its master");
    CHECK(gluedGroupOf(g, g.entities[0]) == 10 && nextId(g) == 21, "glued group or next id wrong");
}

void testExpandArrayPatterns() {
    SceneGraph g;
    g.entities.push_back(boxEntity("Template", {5, 5, 5}, 0.1f));
    g.entities[0].id = 1;
    ArrayObject a;
    a.id = 2;
    a.templateId = 1;
    a.position = {0, 1, 0};
    a.count[0] = 5;
    a.step = {0.3f, 0, 0.1f};
    const std::vector<Entity> line = expandArray(g, a);
    float lineErr = 0;
    for (size_t i = 0; i < line.size(); ++i) lineErr = std::max(lineErr, length(line[i].position - (a.position + a.step * float(i))));
    a.pattern = ArrayPattern::Grid;
    a.count[0] = 3, a.count[1] = 2, a.count[2] = 2;
    a.step = {0.2f, 0.3f, 0.4f};
    const std::vector<Entity> grid = expandArray(g, a);
    const Vector3 last = grid.empty() ? Vector3(0.0f) : grid.back().position - a.position; // (2, 1, 1) * spacing
    a.pattern = ArrayPattern::Circle;
    a.count[0] = 8;
    a.radius = 1.5f;
    const std::vector<Entity> circle = expandArray(g, a);
    float radiusErr = 0;
    for (const Entity& c : circle) radiusErr = std::max(radiusErr, std::fabs(length(c.position - a.position) - a.radius));
    a.pattern = ArrayPattern::Line;
    a.count[0] = 6;
    a.jitter = 0.05f;
    const std::vector<Entity> j1 = expandArray(g, a), j2 = expandArray(g, a);
    a.seed = 99;
    const std::vector<Entity> j3 = expandArray(g, a);
    bool same = j1.size() == j2.size(), differs = false;
    for (size_t i = 0; i < j1.size() && same; ++i) same = length(j1[i].position - j2[i].position) == 0;
    for (size_t i = 0; i < j1.size() && i < j3.size(); ++i) differs = differs || length(j1[i].position - j3[i].position) > 1e-4f;
    std::printf("  line: %zu copies, worst spacing error %.1e m; grid: %zu copies, last at (%.2f %.2f %.2f); circle: %zu, "
                "worst radius error %.1e m; jitter same seed %s, other seed %s\n", line.size(), lineErr, grid.size(), last.x, last.y,
                last.z, circle.size(), radiusErr, same ? "same" : "DIFFERENT", differs ? "differs" : "SAME");
    CHECK(line.size() == 5 && lineErr < 1e-6f, "line copies misplaced");
    CHECK(grid.size() == 12 && length(last - Vector3(0.4f, 0.3f, 0.4f)) < 1e-6f, "grid copies misplaced");
    CHECK(circle.size() == 8 && radiusErr < 1e-6f, "circle copies off the radius");
    CHECK(same && differs, "the jitter must repeat with the seed and change with another");
    CHECK(isArrayCopyId(line[3].id) && arrayOfCopyId(line[3].id) == 2 && line[3].id == arrayCopyId(2, 3), "copy ids not stable");
}

// Fifty cubes from one 5 x 2 x 5 array fall onto the floor, settle into a pile and sleep; then the
// array grows to 80 copies while the scene plays, and a resting box elsewhere does not notice.
void testArrayPileAndGrow() {
    SceneGraph g;
    g.world.size = {4, 3, 4};
    g.entities.push_back(floorEntity(4));
    Entity tmpl = boxEntity("Cube", {0, 0, 0}, 0.1f);
    tmpl.visible = false; // the template: the array builds the copies, not the template itself
    g.entities.push_back(tmpl);
    g.entities.push_back(boxEntity("Resting box", {1.7f, 0.05f, 1.7f}, 0.1f));
    ArrayObject a;
    a.name = "Pile";
    a.pattern = ArrayPattern::Grid;
    a.count[0] = 5, a.count[1] = 2, a.count[2] = 5;
    a.step = {0.13f, 0.12f, 0.13f};
    a.position = {-0.3f, 0.3f, -0.3f};
    a.jitter = 0.005f;
    g.arrays.push_back(a);
    auto owned = std::make_unique<GraphScene>(g);
    GraphScene* scene = owned.get();
    Simulation sim;
    sim.load(std::move(owned));
    a = scene->graph().arrays[0];
    a.templateId = scene->graph().entities[1].id;
    scene->rebuildArray(sim, a); // the template's id is known only now: its copies appear
    const int built = aliveBodies(sim);
    steps(sim, 300);
    int asleep = 0;
    bool finite = true;
    for (const RigidBody& b : sim.rigid.bodies()) {
        asleep += b.alive && b.invMass > 0 && b.sleeping ? 1 : 0;
        finite = finite && std::isfinite(b.pos.x + b.pos.y + b.pos.z);
    }
    const uint32_t restingId = scene->graph().entities[2].id;
    const Vector3 restBefore = sim.rigid.bodies()[size_t(bodyOf(*scene, restingId))].pos;
    const auto copiesBefore = copyPositions(sim, *scene);
    a.count[0] = 8; // 8 x 2 x 5 = 80
    scene->rebuildArray(sim, a);
    const int grown = aliveBodies(sim);
    const Vector3 restAfter = sim.rigid.bodies()[size_t(bodyOf(*scene, restingId))].pos;
    // The 50 old copies are rebuilt where they lie (their live state is read by flat index: this
    // read past graph().entities before, and put the copies wherever that memory said).
    const float copyMoved = worstCopyMove(sim, *scene, copiesBefore);
    steps(sim, 60);
    std::printf("  array 5x2x5: %d bodies (floor, box, 50 copies), %d of 51 asleep after 5 s; grown to 8x2x5 during play: %d bodies, "
                "resting box moved %.1e m, the %zu old copies at most %.1e m\n", built, asleep, grown, length(restAfter - restBefore),
                copiesBefore.size(), copyMoved);
    CHECK(copiesBefore.size() == 50 && copyMoved < 1e-5f, "the old copies moved when the array grew (%.3g m)", copyMoved);
    CHECK(built == 52 && finite, "the array did not build 50 bodies (%d)", built);
    CHECK(asleep == 51, "only %d of 51 dynamic bodies asleep after 5 s", asleep);
    CHECK(grown == 82, "the grown array has %d bodies, expected 82", grown);
    CHECK(length(restAfter - restBefore) < 1e-6f, "growing the array moved the resting box");
}

// Ten boxes are instances of one master: making the master twice as dense makes all ten heavier.
void testInstancesShareRoles() {
    SceneGraph g;
    g.world.size = {4, 3, 4};
    g.entities.push_back(floorEntity(4));
    Entity master = boxEntity("Master", {-1.5f, 0.05f, 0}, 0.1f);
    master.id = 50;
    g.entities.push_back(master);
    for (int i = 0; i < 9; ++i) {
        Entity e = boxEntity("Instance", {-1.2f + 0.3f * float(i), 0.05f, 0}, 0.3f); // its own size is ignored
        e.instanceOf = 50;
        g.entities.push_back(e);
    }
    auto owned = std::make_unique<GraphScene>(g);
    GraphScene* scene = owned.get();
    Simulation sim;
    sim.load(std::move(owned));
    steps(sim, 10);
    auto masses = [&](float& lo, float& hi) {
        lo = 1e9f, hi = 0;
        for (const Entity& e : scene->graph().entities) {
            if (e.id != 50 && e.instanceOf != 50) continue;
            const float m = sim.rigid.bodies()[size_t(bodyOf(*scene, e.id))].mass;
            lo = std::min(lo, m), hi = std::max(hi, m);
        }
    };
    float lo0, hi0, lo1, hi1;
    masses(lo0, hi0);
    Entity heavier = scene->graph().entities[1];
    heavier.rigid.density *= 2;
    scene->rebuildEntity(sim, heavier);
    masses(lo1, hi1);
    std::printf("  10 bodies of one master: mass %.3f..%.3f kg, after doubling the master's density %.3f..%.3f kg\n", lo0, hi0, lo1, hi1);
    CHECK(std::fabs(hi0 - lo0) < 1e-6f && std::fabs(lo0 - 0.5f) < 1e-3f, "instances do not share the master's size and density");
    CHECK(std::fabs(lo1 - 2 * lo0) < 1e-4f && std::fabs(hi1 - 2 * hi0) < 1e-4f, "not every instance got heavier");
}

// Three boxes in an L glued into one body: it falls and tumbles, the boxes never move apart.
void testGluedGroupTumbles() {
    SceneGraph g;
    g.world.size = {4, 4, 4};
    g.entities.push_back(floorEntity(4));
    Group grp;
    grp.id = 100;
    grp.glued = true;
    grp.position = {0, 1.5f, 0};
    grp.rotationDeg = {20, 0, 30};
    g.groups.push_back(grp);
    const Vector3 local[3] = {{0, 0, 0}, {0.2f, 0, 0}, {0, 0.2f, 0}};
    for (int i = 0; i < 3; ++i) {
        Entity e = boxEntity("Part", local[i], 0.2f);
        e.id = uint32_t(101 + i);
        e.parent = 100;
        e.rigid.angularVelocity = {0, 0, 3};
        g.entities.push_back(e);
    }
    auto owned = std::make_unique<GraphScene>(g);
    GraphScene* scene = owned.get();
    Simulation sim;
    sim.load(std::move(owned));
    const int body = bodyOf(*scene, 101);
    const bool oneBody = body >= 0 && bodyOf(*scene, 102) == body && bodyOf(*scene, 103) == body && aliveBodies(sim) == 2;
    auto partWorld = [&](uint32_t id) { // a part's centre from the shared body (bodyOffset in the part's frame)
        const MetaObject& m = scene->metaObjects(id)[0];
        const RigidBody& b = sim.rigid.bodies()[size_t(m.handle)];
        const Quaternion er = b.rot * m.bodyToEntity;
        return b.pos - er.rotate(m.bodyOffset);
    };
    const float d12 = length(partWorld(101) - partWorld(102)), d13 = length(partWorld(101) - partWorld(103));
    float worst = 0, maxSpin = 0;
    const float y0 = sim.rigid.bodies()[size_t(body)].pos.y;
    for (int f = 0; f < 120; ++f) {
        steps(sim, 1);
        worst = std::max({worst, std::fabs(length(partWorld(101) - partWorld(102)) - d12), std::fabs(length(partWorld(101) - partWorld(103)) - d13)});
        maxSpin = std::max(maxSpin, length(sim.rigid.bodies()[size_t(body)].angVel));
    }
    const float y1 = sim.rigid.bodies()[size_t(body)].pos.y;
    std::printf("  glued L of 3 boxes: one body %s, fell %.2f m, spin up to %.2f rad/s, parts moved apart by at most %.1e m\n",
                oneBody ? "yes" : "NO", y0 - y1, maxSpin, worst);
    CHECK(oneBody, "the glued group is not one body");
    CHECK(y0 - y1 > 0.5f && maxSpin > 1.0f, "the glued group did not fall and tumble");
    CHECK(worst < 1e-3f, "the glued parts moved apart by %e m", double(worst));
}
