// The standard scenes for non-convex rigid bodies, each with a number to check it against (tests
// in NonConvexTests.cpp): a race down a slope (sphere, solid cylinder and a wheel of boxes roll with
// a = g sin(t) / (1 + I / mR^2), so they arrive in that order), a bowl with 60 things poured in, a
// stack of eight tables built of boxes, nesting cups (exact ones would sit 20.4 mm apart, their
// convex parts 30 mm), and teapots, bunnies and rings on the 51 200 triangles of the terrain.
#include "samples/Models.h"
#include "samples/Samples.h"
#include "samples/SamplesInternal.h"

#include <cmath>
#include <memory>

namespace rf {

namespace {

// A static body put down on the arena floor (y = 0) by its lowest point.
void putOnFloor(Simulation& sim, int body) {
    RigidBody& b = sim.rigid.bodies()[size_t(body)];
    b.pos.y -= b.posed().support(Vector3(0, -1, 0)).y;
}

class RaceScene : public Scene {
public:
    void build(Simulation& sim) override {
        sim.useRigidArena(AABB({-6.0f, 0.0f, -2.0f}, {6.0f, 4.0f, 2.0f}));
        const float theta = 15.0f * kPi / 180.0f, R = 0.3f;
        const Quaternion tilt = Quaternion::fromAxisAngle({0, 0, 1}, -theta);
        const Vector3 down(std::cos(theta), -std::sin(theta), 0), up(std::sin(theta), std::cos(theta), 0);
        const Vector3 centre(-1.2f, 4.5f * std::sin(theta) - 0.1f * std::cos(theta), 0); // the low end's top at the floor
        sim.rigid.addBox(centre, {4.5f, 0.1f, 1.0f}, tilt, 0.0f, {0.55f, 0.5f, 0.45f});
        const Vector3 start = centre - down * 4.0f + up * (0.1f + R + 0.001f);
        sim.rigid.addSphere(start + Vector3(0, 0, -0.6f), R, 1000.0f, {0.8f, 0.3f, 0.25f});
        sim.rigid.addConvex(primitives::cylinder(R, 0.25f, 48), start, Quaternion(), 1000.0f, {0.3f, 0.65f, 0.35f});
        sim.rigid.addCompound(wheelShape(), start + Vector3(0, 0, 0.6f), Quaternion(), 1000.0f, {0.75f, 0.75f, 0.78f});
    }
};

class BowlScene : public Scene {
public:
    void build(Simulation& sim) override {
        sim.useRigidArena(Simulation::kDefaultArena);
        putOnFloor(sim, sim.rigid.addCompound(bowlShape(), {0, 1, 0}, Quaternion(), 0.0f, {0.85f, 0.8f, 0.7f}));
        // Poured in a loose spiral, a box, a ball and a capsule in turn. The radius does not follow the
        // kind: 20 boxes dropped in one column at the centre landed on one another, the tower toppled
        // and three boxes went over the rim.
        for (int i = 0; i < 60; ++i) {
            const float r = 0.05f + 0.04f * float(i % 5), a = 2.4f * float(i);
            const Vector3 p(r * std::cos(a), 0.8f + 0.05f * float(i), r * std::sin(a));
            const Quaternion q = Quaternion::fromAxisAngle(normalize(Vector3(1, 0.3f * float(i), 0.7f)), 0.37f * float(i));
            const Vector3 colour(0.3f + 0.6f * float(i % 5) / 4, 0.4f + 0.1f * float(i % 3), 0.9f - 0.6f * float(i % 5) / 4);
            if (i % 3 == 0) sim.rigid.addBox(p, Vector3(0.05f), q, 600.0f, colour);
            else if (i % 3 == 1) sim.rigid.addSphere(p, 0.05f, 600.0f, colour);
            else sim.rigid.addBody(std::make_shared<CapsuleShape>(0.04f, 0.05f), p, q, 600.0f, colour);
        }
    }
};

// Eight tables of five boxes stacked 5 cm apart, each shifted up to 1 cm and turned up to 1.6
// degrees: every leg lands on the top below (the legs stand 2 cm in from its edges). Shifted by 3 cm
// and turned by 5 degrees, a leg lands past the edge, the table stands on three legs with its
// centre of mass on their diagonal and the stack topples - as it should.
class TablesScene : public Scene {
public:
    void build(Simulation& sim) override {
        sim.useRigidArena(Simulation::kDefaultArena);
        const auto table = tableShape();
        const float height = 0.44f, com = table->centerOfMass().y;
        for (int k = 0; k < 8; ++k) {
            const Vector3 at(0.01f * std::sin(1.7f * float(k)), com + float(k) * (height + 0.05f) + 0.05f, 0.01f * std::cos(2.3f * float(k)));
            sim.rigid.addCompound(table, at, Quaternion::fromAxisAngle({0, 1, 0}, 0.029f * std::sin(3.1f * float(k))), 700.0f,
                                  {0.7f - 0.04f * float(k), 0.5f, 0.3f});
        }
    }
};

class CupsScene : public Scene {
public:
    void build(Simulation& sim) override {
        sim.useRigidArena(AABB({-1.0f, 0.0f, -1.0f}, {1.0f, 2.0f, 1.0f}));
        for (int i = 0; i < 6; ++i)
            sim.rigid.addCompound(cupShape(), {0.002f * float(i % 2), 0.15f + 0.12f * float(i), 0}, Quaternion(), 1000.0f,
                                  {0.95f, 0.9f - 0.1f * float(i % 4), 0.4f + 0.1f * float(i)});
    }
};

class TerrainModelsScene : public Scene {
public:
    void configure(Simulation& sim) override { configureTerrain(sim); }
    void build(Simulation& sim) override {
        sim.useRigidArena(AABB({-6.0f, -1.0f, -6.0f}, {6.0f, 8.0f, 6.0f}));
        const std::shared_ptr<const CompoundShape> kinds[3] = {teapotShape(), bunnyShape(), ringShape()};
        const Vector3 colours[3] = {{0.95f, 0.95f, 0.92f}, {0.8f, 0.65f, 0.5f}, {0.62f, 0.64f, 0.68f}};
        for (int i = 0; i < 45; ++i) {
            const Vector3 p(-4.5f + 0.2f * float((i * 37) % 45), 2.5f + 0.1f * float(i), -4.5f + 0.2f * float((i * 17) % 45));
            const Quaternion q = Quaternion::fromAxisAngle(normalize(Vector3(1, 0.5f * float(i), 0.3f)), 0.7f * float(i));
            sim.rigid.addCompound(kinds[i % 3], p, q, i % 3 == 2 ? 7800.0f : 600.0f, colours[i % 3]);
        }
    }
};

} // namespace

void addNonConvexSamples(std::vector<SampleEntry>& out) {
    out.push_back({Preset::RigidRace, "Твёрдые тела", "Скатывание: шар, цилиндр и колесо из 55 ящиков", [] { return std::unique_ptr<Scene>(new RaceScene); }});
    out.push_back({Preset::RigidBowl, "Твёрдые тела", "Невыпуклые: чаша и 60 тел", [] { return std::unique_ptr<Scene>(new BowlScene); }});
    out.push_back({Preset::RigidTables, "Твёрдые тела", "Невыпуклые: стопка из 8 столов", [] { return std::unique_ptr<Scene>(new TablesScene); }});
    out.push_back({Preset::RigidCups, "Твёрдые тела", "Невыпуклые: чашки вкладываются одна в другую", [] { return std::unique_ptr<Scene>(new CupsScene); }});
    out.push_back({Preset::TerrainModels, "Твёрдые тела", "Невыпуклые: чайники, кролики и кольца на рельефе", [] { return std::unique_ptr<Scene>(new TerrainModelsScene); }});
}

} // namespace rf
