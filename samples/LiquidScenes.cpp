// Liquid scenes: a tank of particles (PBF) with whatever bodies, obstacles and emitters go into
// it. Each scene is a Scene: configure() sets the parameters once when it is loaded, build()
// fills the tank at every reset.
#include "samples/Models.h"
#include "samples/Samples.h"

namespace rf {

// A column of water at the left end of the tank collapses (the classic dam break).
class DamBreakScene : public Scene {
public:
    void build(Simulation& sim) override {
        sim.useLiquidTank(Simulation::kDefaultTank);
        const AABB& d = sim.particles.domain();
        sim.particles.addBlock(AABB(d.lo, {d.lo.x + 0.6f, 0.7f, d.hi.z}));
    }
};

// The same column, with a sphere in the way of the wave.
class FluidObstacleScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.obstacle.shape = ObstacleShape::Sphere;
        sim.obstacle.size = 0.35f;
        sim.obstacle.position = {0.3f, 0.12f, 0.0f};
    }
    void build(Simulation& sim) override {
        sim.useLiquidTank(Simulation::kDefaultTank);
        const AABB& d = sim.particles.domain();
        sim.particles.addBlock(AABB(d.lo, {d.lo.x + 0.6f, 0.7f, d.hi.z}));
    }
};

// Bodies of different densities dropped into a pool: wood floats, the steel-like cube sinks.
class FloatingBodiesScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.particles.params.particleRadius = 0.016f;
    }
    void build(Simulation& sim) override {
        sim.useLiquidTank(Simulation::kDefaultTank);
        const AABB& d = sim.particles.domain();
        sim.rigid.addBox({-0.5f, 0.75f, 0.0f}, {0.12f, 0.06f, 0.1f}, Quaternion::fromAxisAngle({0, 0, 1}, 0.4f), 400.0f,
                         {0.85f, 0.55f, 0.25f});
        sim.rigid.addBox({0.0f, 0.9f, 0.05f}, {0.1f, 0.1f, 0.1f}, Quaternion::fromAxisAngle({1, 1, 0}, 0.7f), 600.0f,
                         {0.3f, 0.75f, 0.35f});
        sim.rigid.addSphere({0.5f, 0.8f, -0.1f}, 0.09f, 300.0f, {0.9f, 0.3f, 0.3f});
        sim.rigid.addBox({0.4f, 1.05f, 0.1f}, {0.06f, 0.06f, 0.06f}, Quaternion(), 2500.0f, {0.5f, 0.5f, 0.55f});
        sim.particles.addBlock(AABB(d.lo, {d.hi.x, 0.32f, d.hi.z}));
    }
};

// A jet of water from an emitter hits a wing set at an angle.
class JetOnObjectScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.obstacle.shape = ObstacleShape::Wing;
        sim.obstacle.size = 0.5f;
        sim.obstacle.span = 0.6f;
        sim.obstacle.angleOfAttackDeg = 25.0f;
        sim.obstacle.position = {0.1f, 0.35f, 0.0f};
        sim.particles.emitter.enabled = true;
        sim.particles.emitter.position = {-0.9f, 0.75f, 0.0f};
        sim.particles.emitter.direction = normalize(Vector3(1.0f, -0.25f, 0.0f));
        sim.particles.emitter.radius = 0.07f;
        sim.particles.emitter.speed = 3.0f;
        sim.particles.params.maxParticles = 80000;
    }
    void build(Simulation& sim) override {
        sim.useLiquidTank(Simulation::kDefaultTank);
        const AABB& d = sim.particles.domain();
        sim.particles.addBlock(AABB(d.lo, {d.hi.x, 0.08f, d.hi.z}));
    }
};

// A pool with floating things; a column of water at the left end collapses into a wave that runs
// through them. Drawn as a water surface (screen-space) instead of spheres.
class WaterScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.particles.params.particleRadius = 0.012f; // ~40 000 particles: a smooth enough surface
        sim.vis.liquidSurface = true;
    }
    void build(Simulation& sim) override {
        sim.useLiquidTank(Simulation::kDefaultTank);
        const AABB& d = sim.particles.domain();
        // Densities: pine 500, plank 600, beach ball 80, teapot (hollow in reality, solid here) 700,
        // steel-like ball 3000 (sinks), soft foam 300.
        sim.rigid.addBox({-0.25f, 0.3f, -0.15f}, Vector3(0.07f), Quaternion::fromAxisAngle({0, 1, 0}, 0.4f), 500.0f,
                         {0.85f, 0.6f, 0.3f});
        sim.rigid.addBox({0.1f, 0.3f, 0.15f}, {0.16f, 0.025f, 0.06f}, Quaternion::fromAxisAngle({0, 1, 0}, -0.3f), 600.0f,
                         {0.7f, 0.45f, 0.25f});
        sim.rigid.addSphere({0.45f, 0.35f, -0.1f}, 0.08f, 80.0f, {0.95f, 0.3f, 0.25f});
        sim.rigid.addSphere({0.7f, 0.5f, 0.2f}, 0.05f, 3000.0f, {0.6f, 0.62f, 0.66f});
        sim.rigid.addCompound(teapotShape(), {0.3f, 0.4f, 0.2f}, Quaternion::fromAxisAngle({0, 1, 0}, 2.5f), 700.0f,
                              {0.95f, 0.95f, 0.92f});
        TriMesh foam = primitives::box(Vector3(0.06f));
        foam.translate({0.75f, 0.35f, -0.2f});
        sim.particles.addSoftBody(foam, 300.0f, 0.4f, {0.4f, 0.85f, 0.4f});
        sim.particles.addBlock(AABB(d.lo, {d.hi.x, 0.2f, d.hi.z}));                      // the pool
        sim.particles.addBlock(AABB({d.lo.x, 0.2f, d.lo.z}, {d.lo.x + 0.45f, 0.75f, d.hi.z})); // the column
    }
};

void addLiquidSamples(std::vector<SampleEntry>& out) {
    out.push_back({Preset::DamBreak, "Жидкость", "Жидкость: разрушение плотины", [] { return std::unique_ptr<Scene>(new DamBreakScene); }});
    out.push_back({Preset::FluidObstacle, "Жидкость", "Жидкость: поток на препятствие", [] { return std::unique_ptr<Scene>(new FluidObstacleScene); }});
    out.push_back({Preset::FloatingBodies, "Жидкость", "Жидкость: плавающие тела", [] { return std::unique_ptr<Scene>(new FloatingBodiesScene); }});
    out.push_back({Preset::JetOnObject, "Жидкость", "Жидкость: струя на объект", [] { return std::unique_ptr<Scene>(new JetOnObjectScene); }});
    out.push_back({Preset::Water, "Жидкость", "Вода: волна в бассейне, плавающие тела (шейдер)", [] { return std::unique_ptr<Scene>(new WaterScene); }});
}

} // namespace rf
