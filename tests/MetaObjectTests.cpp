// Meta-objects of the scene graph: one entity changes what it is made of between frames and the
// rest of the scene does not notice. Water turns into jelly and back where it is, moving as it
// moved; a plane turns into cloth and back, the freed body slot is reused and the editor's
// body -> entity map stays right; a magnet's role goes off and on without a reload; fifty changes
// leave no trace in memory.
#include "TestRunner.h"

#include "core/Probe.h"
#include "scene/SceneGraph.h"
#include "scene/Simulation.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>

using namespace rf;

namespace {

Entity floorEntity(float width) {
    Entity e;
    e.name = "Floor";
    e.shape = ShapeKind::Plane;
    e.size = {width, 0.02f, width};
    e.rigid.enabled = true;
    e.rigid.fixed = true;
    e.locked = true;
    return e;
}

Entity boxEntity(const char* name, const Vector3& position, float size) {
    Entity e;
    e.name = name;
    e.size = Vector3(size);
    e.position = position;
    e.rigid.enabled = true;
    return e;
}

// The particles of one group: count, centre, and whether all are finite.
struct GroupState {
    int count = 0;
    Vector3 centre{0.0f};
    bool finite = true;
};

GroupState groupState(const Simulation& sim, int group) {
    GroupState s;
    for (size_t i = 0; i < sim.particles.size(); ++i) {
        if (sim.particles.groupOf(int(i)) != group) continue;
        const Vector3& x = sim.particles.positions()[i];
        s.finite &= std::isfinite(x.x + x.y + x.z);
        s.centre += x;
        ++s.count;
    }
    if (s.count > 0) s.centre /= float(s.count);
    return s;
}

// The entity turned into liquid (true) or jelly (false), everything else kept.
Entity madeOfWater(Entity e, bool water) {
    e.liquid.enabled = water;
    e.soft.enabled = !water;
    return e;
}

int particleGroupOf(const GraphScene& scene, uint32_t id) {
    for (const MetaObject& m : scene.metaObjects(id))
        if (m.kind == MetaObject::Kind::Liquid || m.kind == MetaObject::Kind::SoftBody || m.kind == MetaObject::Kind::Cloth) return m.handle;
    return -1;
}

int bodyOf(const GraphScene& scene, uint32_t id) {
    for (const MetaObject& m : scene.metaObjects(id))
        if (m.kind == MetaObject::Kind::RigidBody) return m.handle;
    return -1;
}

void stepFrames(Simulation& sim, int n) {
    for (int f = 0; f < n; ++f) sim.stepFrame();
}

} // namespace

void testMetaWaterSoftCycle() {
    SceneGraph g;
    g.world.size = {3.2f, 1.5f, 1.0f};
    g.entities.push_back(floorEntity(3.2f));
    g.entities.push_back(boxEntity("Resting box", {-1.3f, 0.06f, 0}, 0.1f)); // far from the water
    Entity cube = boxEntity("Cube", {1.0f, 0.4f, 0}, 0.2f);
    cube.rigid.enabled = false;
    cube.liquid.enabled = true;
    g.entities.push_back(cube);
    auto owned = std::make_unique<GraphScene>(g);
    GraphScene* scene = owned.get();
    Simulation sim;
    sim.load(std::move(owned));
    const uint32_t boxId = scene->graph().entities[1].id, cubeId = scene->graph().entities[2].id;
    auto boxPos = [&] { return sim.rigid.bodies()[size_t(bodyOf(*scene, boxId))].pos; };
    stepFrames(sim, 30);
    const Vector3 boxStart = boxPos();
    float worstJump = 0, worstOffset = 0;
    bool countsRight = true, finite = true;
    for (int k = 0; k < 3; ++k) { // water -> jelly -> water -> jelly
        const bool toWater = k % 2 == 1;
        const GroupState before = groupState(sim, particleGroupOf(*scene, cubeId));
        const size_t total = sim.particles.size();
        const Vector3 boxBefore = boxPos();
        scene->rebuildEntity(sim, madeOfWater(scene->graph().entities[2], toWater));
        const GroupState after = groupState(sim, particleGroupOf(*scene, cubeId));
        worstJump = std::max(worstJump, length(boxPos() - boxBefore));
        worstOffset = std::max(worstOffset, length(after.centre - before.centre));
        countsRight &= sim.particles.size() == total - size_t(before.count) + size_t(after.count);
        std::printf("  switch %d -> %s: %d particles out, %d in, centre moved %.3f m\n", k + 1, toWater ? "water" : "jelly",
                    before.count, after.count, length(after.centre - before.centre));
        stepFrames(sim, 30);
        finite &= groupState(sim, particleGroupOf(*scene, cubeId)).finite;
    }
    const float boxDrift = length(boxPos() - boxStart);
    std::printf("  the resting box: jump at the switches %.1e m, drift over the run %.1e m\n", worstJump, boxDrift);
    CHECK(countsRight, "the particle count did not change by exactly the removed and added groups");
    CHECK(worstJump < 1e-6f, "switching the cube moved the resting box by %e m", double(worstJump));
    CHECK(boxDrift < 2e-3f, "the resting box drifted %f m", double(boxDrift));
    CHECK(worstOffset < 0.2f, "the new shape appeared %f m away from the old one", double(worstOffset));
    CHECK(finite, "NaN in the cube's particles");
}

// The rotation read back from a body gives the entity's own Euler angles again.
static void checkPoseRoundTrip(Simulation& sim, GraphScene& scene, uint32_t id) {
    const int index = int(std::find_if(scene.graph().entities.begin(), scene.graph().entities.end(),
                                       [id](const Entity& e) { return e.id == id; }) - scene.graph().entities.begin());
    const Entity before = scene.graph().entities[size_t(index)];
    scene.rebuildEntity(sim, before); // same roles, same pose: the live pose is read back
    const Entity& after = scene.graph().entities[size_t(index)];
    const float angle = length(after.rotationDeg - before.rotationDeg), shift = length(after.position - before.position);
    std::printf("  pose read back from the body: rotation off by %.1e deg, position by %.1e m\n", angle, shift);
    CHECK(angle < 1e-2f && shift < 1e-5f, "the entity's pose did not survive a rebuild: %f deg, %f m", double(angle), double(shift));
}

void testMetaRigidToClothAndBack() {
    SceneGraph g;
    g.world.size = {3, 2, 3};
    g.entities.push_back(floorEntity(3));
    Entity sheet = boxEntity("Sheet", {0, 1.2f, 0}, 0.6f);
    sheet.shape = ShapeKind::Plane;
    sheet.size = {0.6f, 0.02f, 0.6f};
    sheet.rigid.fixed = true;
    g.entities.push_back(sheet);
    Entity turned = boxEntity("Turned", {1.0f, 0.5f, 1.0f}, 0.2f);
    turned.rigid.fixed = true;
    turned.rotationDeg = {10, 20, 30};
    g.entities.push_back(turned);
    auto owned = std::make_unique<GraphScene>(g);
    GraphScene* scene = owned.get();
    Simulation sim;
    sim.load(std::move(owned));
    const uint32_t sheetId = scene->graph().entities[1].id;
    checkPoseRoundTrip(sim, *scene, scene->graph().entities[2].id);
    const int sheetSlot = bodyOf(*scene, sheetId);
    Entity cloth = scene->graph().entities[1];
    cloth.rigid.enabled = false;
    cloth.cloth.enabled = true;
    cloth.cloth.pinnedEdges = 16;
    scene->rebuildEntity(sim, cloth);
    const int bodiesAsCloth = sim.rigid.bodyCount(), clothsAsCloth = int(sim.particles.cloths().size());
    stepFrames(sim, 60);
    // A new ball takes the sheet's freed body slot; the sheet, rigid again, gets another one.
    scene->addEntity(sim, boxEntity("Ball", {-1.0f, 0.5f, -1.0f}, 0.1f));
    const uint32_t ballId = scene->graph().entities.back().id;
    Entity rigidAgain = scene->graph().entities[1];
    rigidAgain.cloth.enabled = false;
    rigidAgain.rigid.enabled = true;
    scene->rebuildEntity(sim, rigidAgain);
    const int ballSlot = bodyOf(*scene, ballId), sheetSlotAgain = bodyOf(*scene, sheetId);
    std::printf("  sheet: body slot %d -> cloth (%d bodies, %d cloths) -> the ball took slot %d, the sheet is body %d again "
                "(%d bodies, %d cloths)\n", sheetSlot, bodiesAsCloth, clothsAsCloth, ballSlot, sheetSlotAgain, sim.rigid.bodyCount(),
                int(sim.particles.cloths().size()));
    CHECK(bodiesAsCloth == 2 && clothsAsCloth == 1, "as cloth: %d bodies, %d cloths", bodiesAsCloth, clothsAsCloth);
    CHECK(ballSlot == sheetSlot, "the freed slot %d was not reused (ball in %d)", sheetSlot, ballSlot);
    CHECK(sim.rigid.bodyCount() == 4 && sim.particles.cloths().empty(), "rigid again: %d bodies, %zu cloths", sim.rigid.bodyCount(),
          sim.particles.cloths().size());
    CHECK(scene->entityIdOfBody(ballSlot) == ballId && scene->entityIdOfBody(sheetSlotAgain) == sheetId,
          "body -> entity map wrong after the slot was reused");
}

void testMetaMagnetToggle() {
    SceneGraph g;
    g.world.gravity = Vector3(0.0f);
    g.world.size = {2, 3, 2};
    for (int i = 0; i < 2; ++i) {
        Entity e = boxEntity(i ? "Right" : "Left", {i ? 0.08f : -0.08f, 1.5f, 0}, 0.02f);
        e.rigid.density = 7800;
        e.magnet.enabled = true;
        e.magnet.moment = {0.5f, 0, 0};
        g.entities.push_back(e);
    }
    auto owned = std::make_unique<GraphScene>(g);
    GraphScene* scene = owned.get();
    Simulation sim;
    sim.load(std::move(owned));
    sim.rigid.params.linearDamping = sim.rigid.params.angularDamping = 0; // only the magnets change the speed
    const uint32_t leftId = scene->graph().entities[0].id, rightId = scene->graph().entities[1].id;
    auto rightSpeed = [&] { return sim.rigid.bodies()[size_t(bodyOf(*scene, rightId))].vel.x; };
    stepFrames(sim, 5);
    const uint64_t frameBefore = sim.frame();
    Entity plain = scene->graph().entities[0];
    plain.magnet.enabled = false;
    scene->rebuildEntity(sim, plain);
    const float v0 = rightSpeed();
    sim.stepFrame();
    const float changeWithout = std::fabs(rightSpeed() - v0);
    Entity magnetAgain = scene->graph().entities[0];
    magnetAgain.magnet.enabled = true;
    scene->rebuildEntity(sim, magnetAgain);
    const float v1 = rightSpeed();
    // Coaxial dipoles: F = 3 mu0 m1 m2 / (2 pi d^4) = 6e-7 m1 m2 / d^4 [N], acting for one frame.
    const RigidBody& right = sim.rigid.bodies()[size_t(bodyOf(*scene, rightId))];
    const float d = length(right.pos - sim.rigid.bodies()[size_t(bodyOf(*scene, leftId))].pos);
    const float expected = 6e-7f * 0.25f / (d * d * d * d) / right.mass * sim.frameDt;
    sim.stepFrame();
    const float changeWith = std::fabs(rightSpeed() - v1);
    std::printf("  right magnet's speed change in one frame: %.2e m/s with the left one plain, %.3e m/s with it a magnet again "
                "(theory %.3e at %.3f m); frame %llu -> %llu (no reload)\n", changeWithout, changeWith, expected, d,
                (unsigned long long)frameBefore, (unsigned long long)sim.frame());
    CHECK(changeWithout < 1e-7f, "a plain body still pulled the magnet: %e m/s", double(changeWithout));
    CHECK(std::fabs(changeWith / expected - 1.0f) < 0.05f, "the magnet role did not come back: %e m/s (theory %e)",
          double(changeWith), double(expected));
    CHECK(sim.frame() == frameBefore + 2, "the scene was reloaded (frame %llu)", (unsigned long long)sim.frame());
    (void)leftId;
}

void testMetaNoGrowth() {
    SceneGraph g;
    g.world.size = {2, 1.5f, 1.0f};
    g.entities.push_back(floorEntity(2));
    Entity water = boxEntity("Water", {0.5f, 0.3f, 0}, 0.16f);
    water.rigid.enabled = false;
    water.liquid.enabled = true;
    g.entities.push_back(water);
    g.entities.push_back(boxEntity("Block", {-0.5f, 0.3f, 0}, 0.12f));
    auto owned = std::make_unique<GraphScene>(g);
    GraphScene* scene = owned.get();
    Simulation sim;
    sim.load(std::move(owned));
    auto allocationsPerFrame = [&] {
        double sum = 0;
        for (int f = 0; f < 10; ++f) {
            sim.stepFrame();
            sum += Probe::snapshot().value("memory/allocations per frame", 0.0);
        }
        return sum / 10.0;
    };
    const double before = allocationsPerFrame();
    const size_t slotsBefore = sim.rigid.bodies().size(), particlesBefore = sim.particles.size();
    for (int k = 0; k < 50; ++k) {
        scene->rebuildEntity(sim, madeOfWater(scene->graph().entities[1], k % 2 == 1)); // water <-> jelly
        Entity block = scene->graph().entities[2];
        block.rigid.enabled = k % 2 == 1; // rigid <-> jelly
        block.soft.enabled = k % 2 == 0;
        scene->rebuildEntity(sim, block);
        sim.stepFrame();
    }
    const double after = allocationsPerFrame();
    std::printf("  50 changes: allocations per frame %.1f -> %.1f, body slots %zu -> %zu, particles %zu -> %zu\n", before, after,
                slotsBefore, sim.rigid.bodies().size(), particlesBefore, sim.particles.size());
    CHECK(after <= 2.0 * before + 20.0, "allocations per frame grew: %f -> %f", before, after);
    CHECK(sim.rigid.bodies().size() <= slotsBefore + 1, "body slots grew: %zu -> %zu", slotsBefore, sim.rigid.bodies().size());
}
