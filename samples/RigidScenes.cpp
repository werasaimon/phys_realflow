// Rigid-body scenes: bodies dropped on a mesh, granular spheres, a pyramid hit by a ball, mixed
// polyhedra, the 100-cube tower, the five joint types, CCD bullets, 100 teapots, and a terrain of
// 50 000 static triangles. Each is a Scene on the Simulation facade (see Scene.h).
#include "samples/Samples.h"
#include "samples/Models.h"

#include <cmath>
#include <memory>

namespace rf {

namespace {

// The same numbers every run: a linear congruential generator with a fixed seed, so a scene
// looks the same each time it is loaded (and the tests see the same bodies).
struct Random {
    uint32_t seed = 7;
    float next() {
        seed = seed * 1664525u + 1013904223u;
        return (seed >> 8) * (1.0f / 16777216.0f);
    }
    Vector3 color() { return Vector3(0.35f + 0.6f * next(), 0.35f + 0.6f * next(), 0.35f + 0.6f * next()); }
};

// Bodies dropped on a big sphere: boxes and spheres of random sizes and orientations.
class RigidFallingScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.obstacle.shape = ObstacleShape::Sphere;
        sim.obstacle.size = 1.6f;
        sim.obstacle.position = {0.0f, 0.0f, 0.0f};
    }
    void build(Simulation& sim) override {
        sim.useRigidArena(Simulation::kDefaultArena);
        Random rnd;
        for (int i = 0; i < 60; ++i) {
            Vector3 p(-1.2f + 2.4f * rnd.next(), 2.0f + 2.5f * rnd.next(), -1.2f + 2.4f * rnd.next());
            if (i % 2 == 0)
                sim.rigid.addBox(p, Vector3(0.08f + 0.12f * rnd.next(), 0.08f + 0.12f * rnd.next(), 0.08f + 0.12f * rnd.next()),
                                 Quaternion::fromAxisAngle({rnd.next(), rnd.next(), rnd.next() + 0.1f}, 6.28f * rnd.next()), 800.0f, rnd.color());
            else
                sim.rigid.addSphere(p, 0.08f + 0.12f * rnd.next(), 800.0f, rnd.color());
        }
    }
};

// 2016 spheres poured over a cone: a granular medium.
class RigidGranularScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.obstacle.shape = ObstacleShape::Cone;
        sim.obstacle.size = 1.2f;
        sim.obstacle.rollDeg = 0.0f;
        sim.obstacle.angleOfAttackDeg = 90.0f; // apex up
        sim.obstacle.position = {0.0f, 0.6f, 0.0f};
        sim.rigid.params.iterations = 8;
        sim.rigid.params.substeps = 4; // loose spheres: no tall stacks, fewer substeps suffice
    }
    void build(Simulation& sim) override {
        sim.useRigidArena(Simulation::kDefaultArena);
        Random rnd;
        const float r = 0.05f;
        for (int k = 0; k < 12; ++k)
            for (int j = 0; j < 14; ++j)
                for (int i = 0; i < 12; ++i) {
                    Vector3 p(-0.6f + i * 2.1f * r + 0.01f * rnd.next(), 2.0f + j * 2.1f * r, -0.6f + k * 2.1f * r + 0.01f * rnd.next());
                    float t = float(j) / 14.0f;
                    sim.rigid.addSphere(p, r * (0.85f + 0.3f * rnd.next()), 2000.0f, {0.85f - 0.3f * t, 0.65f, 0.35f + 0.4f * t});
                }
    }
};

// A pyramid of cubes and a heavy ball thrown at it.
class RigidPyramidScene : public Scene {
public:
    void build(Simulation& sim) override {
        sim.useRigidArena(Simulation::kDefaultArena);
        const float h = 0.15f;
        int levels = 7;
        for (int l = 0; l < levels; ++l)
            for (int i = 0; i < levels - l; ++i) {
                float x = (i - (levels - l - 1) * 0.5f) * (2 * h + 0.005f);
                sim.rigid.addBox({x, h + l * 2 * h, 0.0f}, Vector3(h), Quaternion(), 500.0f,
                                 {0.9f - 0.08f * l, 0.5f + 0.05f * l, 0.3f});
            }
        int ball = sim.rigid.addSphere({-1.8f, 0.6f, 0.0f}, 0.2f, 3000.0f, {0.2f, 0.3f, 0.9f});
        sim.rigid.bodies()[ball].vel = {9.0f, 1.5f, 0.0f};
    }
};

// Mixed convex polyhedra dropped on a slope: exercises SAT (box-box) and GJK/EPA with
// perturbation manifolds (hull-hull, hull-box, hull-mesh).
class RigidConvexScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.obstacle.shape = ObstacleShape::Cube;
        sim.obstacle.size = 1.4f;
        sim.obstacle.position = {0.0f, 0.1f, 0.0f};
        sim.obstacle.angleOfAttackDeg = 20.0f;
        sim.obstacle.yawDeg = 30.0f;
    }
    void build(Simulation& sim) override {
        sim.useRigidArena(Simulation::kDefaultArena);
        Random rnd;
        std::vector<TriMesh> kinds = {primitives::cylinder(0.14f, 0.3f, 12), primitives::cone(0.16f, 0.35f, 12),
                                      primitives::sphere(0.16f, 8, 5), primitives::ellipsoid({0.22f, 0.1f, 0.14f}, 10, 6)};
        TriMesh tet;
        tet.positions = {{0.2f, 0.2f, 0.2f}, {-0.2f, -0.2f, 0.2f}, {-0.2f, 0.2f, -0.2f}, {0.2f, -0.2f, -0.2f}};
        tet.triangles = {{0, 1, 2}, {0, 3, 1}, {0, 2, 3}, {1, 3, 2}};
        tet.orientOutward();
        kinds.push_back(tet);
        for (int i = 0; i < 45; ++i) {
            Vector3 p(-1.3f + 2.6f * rnd.next(), 1.5f + 3.0f * rnd.next(), -1.3f + 2.6f * rnd.next());
            Quaternion q = Quaternion::fromAxisAngle({rnd.next() - 0.5f, rnd.next() - 0.5f, rnd.next() - 0.5f + 1e-3f}, 6.28f * rnd.next());
            if (i % 6 == 5)
                sim.rigid.addBox(p, Vector3(0.1f + 0.1f * rnd.next(), 0.08f + 0.08f * rnd.next(), 0.1f + 0.1f * rnd.next()), q, 700.0f, rnd.color());
            else
                sim.rigid.addConvex(kinds[i % kinds.size()], p, q, 700.0f, rnd.color());
        }
    }
};

// Standard stability test: 100 cubes, each released 1 cm above the one below.
class RigidTowerScene : public Scene {
public:
    void build(Simulation& sim) override {
        sim.useRigidArena(AABB({-3.0f, 0.0f, -3.0f}, {3.0f, 24.0f, 3.0f}));
        const float h = 0.1f, gap = 0.01f;
        for (int i = 0; i < 100; ++i) {
            float t = i / 99.0f;
            sim.rigid.addBox({0.0f, h + i * (2 * h + gap), 0.0f}, Vector3(h), Quaternion(), 500.0f,
                             {0.25f + 0.65f * t, 0.55f + 0.2f * std::sin(6.28f * t), 0.9f - 0.6f * t});
        }
    }
};

// The five joint types, one machine each.
class RigidJointsScene : public Scene {
public:
    void build(Simulation& sim) override {
        sim.useRigidArena(Simulation::kDefaultArena);
        RigidWorld& rigid = sim.rigid;
        const Vector3 grey(0.7f, 0.72f, 0.76f);
        // 1) Ball joints: a chain of 10 links hanging from the ceiling, released horizontally.
        {
            Vector3 top(-1.6f, 4.4f, -1.2f);
            int prev = -1;
            for (int i = 0; i < 10; ++i) {
                Vector3 c = top + Vector3(0.12f + 0.24f * i, 0, 0);
                int id = rigid.addBox(c, {0.1f, 0.03f, 0.03f}, Quaternion(), 800.0f, {0.9f, 0.45f, 0.3f});
                rigid.addBallJoint(id, prev, c - Vector3(0.12f, 0, 0));
                prev = id;
            }
        }
        // 2) Hinges: a motor-driven paddle wheel and a swinging door with angle limits.
        {
            int paddle = rigid.addBox({0.0f, 1.0f, -1.2f}, {0.6f, 0.05f, 0.15f}, Quaternion(), 600.0f, {0.3f, 0.6f, 0.95f});
            HingeJoint& h = rigid.addHingeJoint(paddle, -1, {0.0f, 1.0f, -1.2f}, {0, 0, 1});
            h.motorEnabled = true;
            h.motorSpeed = 1.5f;
            h.maxMotorTorque = 400.0f;
            int door = rigid.addBox({-1.2f, 1.0f, 0.3f}, {0.35f, 0.6f, 0.03f}, Quaternion(), 400.0f, {0.95f, 0.8f, 0.3f});
            HingeJoint& d = rigid.addHingeJoint(door, -1, {-1.55f, 1.0f, 0.3f}, {0, 1, 0});
            d.limitEnabled = true;
            d.lower = -1.2f;
            d.upper = 1.2f;
            rigid.bodies()[door].angVel = {0, 4.0f, 0};
        }
        // 3) Slider: a carriage on a tilted rail with travel limits.
        {
            Vector3 axis = normalize(Vector3(1.0f, -0.4f, 0.0f));
            int cart = rigid.addBox({1.0f, 2.2f, 1.2f}, {0.15f, 0.1f, 0.1f}, Quaternion::fromAxisAngle({0, 0, 1}, std::atan2(axis.y, axis.x)),
                                    700.0f, {0.4f, 0.85f, 0.4f});
            SliderJoint& sl = rigid.addSliderJoint(cart, -1, axis);
            sl.limitEnabled = true;
            sl.lower = -0.9f;
            sl.upper = 0.9f;
        }
        // 4) Fixed: an L-shaped body welded from two boxes, dropped with a spin.
        {
            int a = rigid.addBox({1.2f, 3.0f, -0.3f}, {0.3f, 0.07f, 0.07f}, Quaternion(), 600.0f, {0.8f, 0.4f, 0.9f});
            int b = rigid.addBox({0.97f, 3.23f, -0.3f}, {0.07f, 0.3f, 0.07f}, Quaternion(), 600.0f, {0.8f, 0.4f, 0.9f});
            rigid.addFixedJoint(a, b);
            rigid.bodies()[a].angVel = {1.0f, 0.5f, 2.0f};
            rigid.bodies()[b].angVel = {1.0f, 0.5f, 2.0f};
        }
        // 5) Distance: a rope pendulum and a soft spring.
        {
            int bob = rigid.addSphere({1.5f, 3.4f, -1.4f}, 0.12f, 2000.0f, {0.95f, 0.3f, 0.35f});
            DistanceJoint& rope = rigid.addDistanceJoint(bob, -1, {1.5f, 3.4f, -1.4f}, {0.6f, 4.4f, -1.4f});
            rope.rope = true;
            int weight = rigid.addBox({0.4f, 3.0f, 1.3f}, Vector3(0.12f), Quaternion(), 900.0f, {0.3f, 0.8f, 0.85f});
            DistanceJoint& spring = rigid.addDistanceJoint(weight, -1, {0.4f, 3.12f, 1.3f}, {0.4f, 4.4f, 1.3f});
            spring.frequency = 1.2f;
            spring.dampingRatio = 0.1f;
        }
        // A few free bodies to hit with the paddle or grab with the mouse.
        for (int i = 0; i < 6; ++i)
            rigid.addBox({-0.3f + 0.12f * i, 2.0f + 0.3f * i, -1.2f}, Vector3(0.08f), Quaternion(), 500.0f, grey);
    }
};

// A 2 cm thick wall and bullets that travel 25 cm per substep: without CCD they tunnel.
class RigidCcdScene : public Scene {
public:
    void build(Simulation& sim) override {
        sim.useRigidArena(Simulation::kDefaultArena);
        RigidWorld& rigid = sim.rigid;
        rigid.addBox({0.5f, 1.5f, 0.0f}, {0.01f, 1.5f, 1.5f}, Quaternion(), 0.0f, {0.6f, 0.62f, 0.66f});
        for (int i = 0; i < 12; ++i) {
            int b = rigid.addSphere({-1.8f, 0.4f + 0.2f * i, -1.0f + 0.18f * i}, 0.03f, 8000.0f, {1.0f, 0.9f, 0.3f});
            rigid.bodies()[b].vel = {150.0f, 5.0f, 0.0f};
        }
        // Moving targets: a free stack of boxes and thin plates hanging on hinges - bullets must hit
        // them (momentum transfer, bounce), never pass through or overlap.
        for (int i = 0; i < 6; ++i)
            rigid.addBox({-0.3f, 0.1f + 0.2f * i, -0.6f}, Vector3(0.1f), Quaternion(), 400.0f, {0.35f + 0.1f * i, 0.55f, 0.9f - 0.1f * i});
        for (int i = 0; i < 3; ++i) {
            float z = 0.2f + 0.3f * i;
            int plate = rigid.addBox({0.0f, 1.9f, z}, {0.01f, 0.4f, 0.12f}, Quaternion(), 1500.0f, {0.95f, 0.5f, 0.3f});
            rigid.addHingeJoint(plate, -1, {0.0f, 2.3f, z}, {0, 0, 1});
            int b = rigid.addSphere({-1.8f, 1.75f, z}, 0.03f, 8000.0f, {1.0f, 0.95f, 0.4f});
            rigid.bodies()[b].vel = {120.0f, 0.0f, 0.0f};
        }
        for (int i = 0; i < 6; ++i) {
            int b = rigid.addSphere({-1.8f, 0.1f + 0.2f * i, -0.6f}, 0.03f, 8000.0f, {1.0f, 0.95f, 0.4f});
            rigid.bodies()[b].vel = {150.0f, 0.0f, 0.0f};
        }
        // A thin fast-spinning plate next to a thin post.
        rigid.addBox({-1.0f, 1.5f, 1.2f}, {0.01f, 1.5f, 0.01f}, Quaternion(), 0.0f, {0.6f, 0.62f, 0.66f});
        int plate = rigid.addBox({-1.0f, 2.2f, 0.6f}, {0.02f, 0.02f, 0.5f}, Quaternion(), 800.0f, {0.4f, 0.8f, 1.0f});
        rigid.bodies()[plate].angVel = {0.0f, 120.0f, 0.0f};
    }
};

// Non-convex bodies: all 100 share one compound shape (see teapotShape()).
class RigidTeapotsScene : public Scene {
public:
    void build(Simulation& sim) override {
        sim.useRigidArena(Simulation::kDefaultArena);
        Random rnd;
        const auto teapot = teapotShape();
        for (int i = 0; i < 100; ++i) {
            int layer = i / 25, cell = i % 25;
            Vector3 p(-1.2f + 0.6f * (cell % 5) + 0.08f * (rnd.next() - 0.5f), 0.4f + 0.5f * layer,
                      -1.2f + 0.6f * (cell / 5) + 0.08f * (rnd.next() - 0.5f));
            Quaternion q = Quaternion::fromAxisAngle({rnd.next() - 0.5f, rnd.next() - 0.5f, rnd.next() - 0.5f + 1e-3f}, 6.28f * rnd.next());
            float t = i / 99.0f;
            sim.rigid.addCompound(teapot, p, q, 600.0f, {0.85f - 0.4f * t, 0.45f + 0.35f * std::sin(3.14f * t), 0.3f + 0.6f * t});
        }
    }
};

// A terrain: the ground is one big static triangle mesh (a 12 x 12 m height field of 160 x 160
// cells, 51 200 triangles) with a BVH over it - what Bullet does with btBvhTriangleMeshShape.
// Every body asks the BVH for the triangles under its bounding box and collides with just those
// (ContactSolver::collideStatic), CCD sweeps against them too; the particles use the same
// path, so water would run down these hills as well. 150 mixed bodies dropped from 3-6 m roll
// down the slopes and collect in the bowl. The measured ms/frame goes into the docs table.
class TerrainScene : public Scene {
public:
    void configure(Simulation& sim) override {
        // Hills of 0.3-0.8 m from a few waves, in a shallow bowl so that the bodies stay on the
        // terrain instead of leaving over its edge.
        TriMesh ground = primitives::heightfield(12.0f, 12.0f, 160, 160, [](float x, float z) {
            return 0.4f * std::sin(0.8f * x) * std::cos(0.6f * z) + 0.2f * std::sin(1.5f * x + 0.4f * z) +
                   0.1f * std::cos(2.1f * z - 0.9f * x) + 0.3f * (x * x + z * z) / 36.0f;
        });
        // Obstacle.cpp fits a custom mesh to `size` about the origin: the size is the mesh's own
        // largest extent (scale 1), and the position lifts it so that its lowest point is at y = 0.
        const AABB b = ground.bounds();
        sim.obstacle.shape = ObstacleShape::Custom;
        sim.obstacle.size = maxComp(b.extent());
        sim.obstacle.position = {0.0f, 0.5f * b.extent().y, 0.0f};
        sim.obstacle.customName = "рельеф";
        sim.obstacle.customMesh = std::make_shared<TriMesh>(std::move(ground));
    }
    void build(Simulation& sim) override {
        sim.useRigidArena(AABB({-6.0f, -1.0f, -6.0f}, {6.0f, 8.0f, 6.0f}));
        Random rnd;
        const auto teapot = teapotShape();
        const std::vector<TriMesh> kinds = {primitives::cylinder(0.12f, 0.3f, 12), primitives::cone(0.14f, 0.3f, 12)};
        for (int i = 0; i < 150; ++i) {
            Vector3 p(-4.5f + 9.0f * rnd.next(), 3.0f + 3.0f * rnd.next(), -4.5f + 9.0f * rnd.next());
            Quaternion q = Quaternion::fromAxisAngle({rnd.next() - 0.5f, rnd.next() - 0.5f, rnd.next() - 0.5f + 1e-3f}, 6.28f * rnd.next());
            if (i % 15 == 14)
                sim.rigid.addCompound(teapot, p, q, 600.0f, {0.95f, 0.95f, 0.92f}); // 10 teapots
            else if (i % 5 == 4)
                sim.rigid.addConvex(kinds[i % kinds.size()], p, q, 700.0f, rnd.color());
            else if (i % 2 == 0)
                sim.rigid.addBox(p, Vector3(0.08f + 0.12f * rnd.next(), 0.08f + 0.12f * rnd.next(), 0.08f + 0.12f * rnd.next()), q, 800.0f, rnd.color());
            else
                sim.rigid.addSphere(p, 0.08f + 0.12f * rnd.next(), 800.0f, rnd.color());
        }
    }
};

} // namespace

void addRigidSamples(std::vector<SampleEntry>& out) {
    out.push_back({Preset::RigidFalling, "Твёрдые тела", "Твёрдые тела: падение на меш", [] { return std::unique_ptr<Scene>(new RigidFallingScene); }});
    out.push_back({Preset::RigidGranular, "Твёрдые тела", "Твёрдые частицы: сыпучая среда", [] { return std::unique_ptr<Scene>(new RigidGranularScene); }});
    out.push_back({Preset::RigidPyramid, "Твёрдые тела", "Твёрдые тела: пирамида и снаряд", [] { return std::unique_ptr<Scene>(new RigidPyramidScene); }});
    out.push_back({Preset::RigidConvex, "Твёрдые тела", "Твёрдые тела: многогранники (SAT, GJK-EPA)", [] { return std::unique_ptr<Scene>(new RigidConvexScene); }});
    out.push_back({Preset::RigidTower, "Твёрдые тела", "Твёрдые тела: башня из 100 кубиков", [] { return std::unique_ptr<Scene>(new RigidTowerScene); }});
    out.push_back({Preset::RigidJoints, "Твёрдые тела", "Сочленения: 5 типов", [] { return std::unique_ptr<Scene>(new RigidJointsScene); }});
    out.push_back({Preset::RigidCcd, "Твёрдые тела", "CCD: пули и тонкая стена", [] { return std::unique_ptr<Scene>(new RigidCcdScene); }});
    out.push_back({Preset::RigidTeapots, "Твёрдые тела", "Невыпуклые: 100 чайников (выпуклая декомпозиция)", [] { return std::unique_ptr<Scene>(new RigidTeapotsScene); }});
    out.push_back({Preset::Terrain, "Твёрдые тела", "Твёрдые тела: рельеф из 50 000 треугольников (статичный меш + BVH)", [] { return std::unique_ptr<Scene>(new TerrainScene); }});
}

} // namespace rf
