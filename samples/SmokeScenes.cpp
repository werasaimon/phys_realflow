// Smoke scenes: a hot source in a box of gas - a plume in a tunnel, a closed box that fills up,
// and bodies dropped through the rising smoke (two-way coupling).
#include "samples/Samples.h"
#include "samples/SamplesInternal.h"
#include "samples/Models.h"

namespace rf {

void configureClosedSmokeBox(Simulation& sim, bool bodiesInside) {
    // Closed box of voxels, hot smoky sphere near the bottom, gas rises by buoyancy.
    sim.grid.params.domainSize = {1.6f, 2.4f, 1.6f};
    sim.grid.params.resolutionX = 32;               // dx = 5 cm -> 32 x 48 x 32 cells
    sim.grid.params.inflowSpeed = 1.0f;             // only a reference scale for Cp here
    for (auto& b : sim.grid.params.bc) b = BoundaryType::Wall;
    sim.grid.params.smokeRake = false;
    sim.grid.params.heatBuoyancy = 3.0f;
    sim.grid.params.smokeBuoyancy = 0.2f;
    sim.grid.params.temperatureDissipation = 0.4f;
    sim.grid.params.smokeDissipation = 0.08f;       // tracer fades slowly so the closed box does not fog up
    sim.grid.params.vorticityConfinement = 1.0f;
    sim.grid.source.enabled = true;
    sim.grid.source.center = {0.0f, -0.8f, 0.0f};
    sim.grid.source.radius = 0.15f;
    sim.vis.sliceField = GridField::Speed;
    sim.vis.showSlice = false;
    sim.vis.showStreamlines = false;
    // Scenes with bodies: smoke only (no voxel grid / vectors), 3/4 view, dense smoke that
    // leaves through the open ceiling.
    sim.vis.gridDisplay = bodiesInside ? 0 : 1;
    sim.vis.vectorDisplay = bodiesInside ? 0 : 1;
    if (bodiesInside) { // dense smoke, open ceiling: the plume leaves, the box does not fill up
        sim.grid.source.smoke = 2.0f;
        sim.grid.source.radius = 0.2f;
        sim.grid.params.smokeDissipation = 0.3f; // smoke lives ~3 s: the plume stays dense, the box clears
        sim.grid.params.bc[3] = BoundaryType::Outflow;
    }
}

class SmokePlumeScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.obstacle.shape = ObstacleShape::Sphere;
        sim.obstacle.size = 0.45f;
        sim.obstacle.position = {0.0f, 0.1f, 0.0f};
        sim.grid.params.domainSize = {2.0f, 3.0f, 2.0f};
        sim.grid.params.resolutionX = 48;
        sim.grid.params.inflowSpeed = 1.0f;
        sim.grid.params.bc[0] = BoundaryType::Wall;
        sim.grid.params.bc[1] = BoundaryType::Wall;
        sim.grid.params.bc[3] = BoundaryType::Outflow;
        sim.grid.params.smokeRake = false;
        sim.grid.params.heatBuoyancy = 4.0f;
        sim.grid.params.smokeBuoyancy = 0.3f;
        sim.grid.params.vorticityConfinement = 1.5f;
        sim.grid.params.smokeDissipation = 0.05f;
        sim.grid.source.enabled = true;
        sim.grid.source.center = {0.0f, -1.25f, 0.0f};
        sim.grid.source.radius = 0.2f;
        sim.vis.sliceField = GridField::Temperature;
        sim.vis.showStreamlines = false;
    }
    void build(Simulation& sim) override { sim.useGasBox({0.5f, 0.5f, 0.5f}); }
};

class SmokeSphereScene : public Scene {
public:
    void configure(Simulation& sim) override { configureClosedSmokeBox(sim, false); }
    void build(Simulation& sim) override { sim.useGasBox({0.5f, 0.5f, 0.5f}); }
};

class SmokeBodiesScene : public Scene {
public:
    void configure(Simulation& sim) override { configureClosedSmokeBox(sim, true); }
    void build(Simulation& sim) override {
        sim.useGasBox({0.5f, 0.5f, 0.5f});
        RigidWorld& rigid = sim.rigid;
        // A few bodies released above the hot plume: they fall through the smoke, push it aside
        // and shed vortices; the rising gas acts back on them.
        rigid.addBox({-0.35f, 0.5f, 0.1f}, Vector3(0.12f), Quaternion::fromAxisAngle({1, 0, 1}, 0.5f), 400.0f, {0.3f, 0.7f, 0.9f});
        rigid.addBox({0.3f, 0.9f, -0.15f}, {0.18f, 0.05f, 0.14f}, Quaternion::fromAxisAngle({0, 0, 1}, 0.3f), 300.0f,
                     {0.9f, 0.6f, 0.25f});
        rigid.addSphere({0.05f, 0.2f, 0.0f}, 0.12f, 500.0f, {0.9f, 0.35f, 0.3f});
        rigid.addSphere({-0.2f, 1.0f, -0.3f}, 0.09f, 500.0f, {0.95f, 0.85f, 0.3f});
        rigid.addConvex(primitives::cylinder(0.1f, 0.3f, 16), {0.35f, 0.3f, 0.35f}, Quaternion::fromAxisAngle({1, 0, 0}, 1.2f),
                        500.0f, {0.6f, 0.85f, 0.4f});
        rigid.addCompound(teapotShape(), {-0.4f, 0.85f, 0.35f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.8f), 500.0f,
                          {0.8f, 0.55f, 0.85f});
        rigid.addCompound(teapotShape(), {0.4f, 0.55f, -0.35f}, Quaternion::fromAxisAngle({1, 0, 0.3f}, 2.2f), 500.0f,
                          {0.95f, 0.75f, 0.45f});
        rigid.addCompound(bunnyShape(), {0.0f, 0.95f, 0.0f}, Quaternion::fromAxisAngle({0, 1, 0}, -0.6f), 500.0f,
                          {0.92f, 0.9f, 0.86f}); // right above the plume: falls through it
        // They hang still (the smoke flows around them) until the plume has risen, then drop through it.
        sim.releaseTime = 1.5f;
        for (int i = 0; i < int(rigid.bodies().size()); ++i) rigid.hold(i);
    }
};

void addSmokeSamples(std::vector<SampleEntry>& out) {
    out.push_back({Preset::SmokePlume, "Газ и дым", "Газ: тепловой шлейф дыма", [] { return std::unique_ptr<Scene>(new SmokePlumeScene); }});
    out.push_back({Preset::SmokeSphere, "Газ и дым", "Дым: сферический источник (закрытый объём)", [] { return std::unique_ptr<Scene>(new SmokeSphereScene); }});
    out.push_back({Preset::SmokeBodies, "Газ и дым", "Дым + твёрдые тела (двусторонняя связь)", [] { return std::unique_ptr<Scene>(new SmokeBodiesScene); }});
}

} // namespace rf
