// Hydrodynamics: water and air in one box - wind from the left over a pool, a water column that
// collapses into a wave, floating and sinking bodies, a flag in the wind.
#include "samples/Samples.h"
#include "samples/Models.h"

namespace rf {

class HydroScene : public Scene {
public:
    void configure(Simulation& sim) override {
        // Water and air in one box: wind from the left over a pool, a water column that collapses
        // into a wave, floating and sinking bodies, a flag in the wind.
        sim.grid.params.domainSize = {2.0f, 1.2f, 0.8f};
        // The air only has to carry the smoke streaks and push spray and flag: a coarser grid,
        // a looser pressure tolerance and larger (still stable, semi-Lagrangian) steps.
        sim.grid.params.resolutionX = 48;                  // dx ~4 cm
        sim.grid.params.pressureTolerance = 1e-3f;
        sim.grid.params.cfl = 4.0f;
        sim.grid.params.inflowSpeed = 3.0f;
        for (auto& b : sim.grid.params.bc) b = BoundaryType::Wall;
        sim.grid.params.bc[0] = BoundaryType::Inflow;       // wind
        sim.grid.params.bc[1] = BoundaryType::Outflow;
        sim.grid.params.bc[3] = BoundaryType::Outflow;      // open top
        sim.grid.params.smokeRake = true;                   // streaks show the air flow
        sim.grid.params.smokeDissipation = 0.15f;
        sim.grid.params.vorticityConfinement = 0.5f;
        sim.particles.params.particleRadius = 0.018f;
        sim.particles.params.clothSpacing = 2.0f;
        sim.vis.showSlice = false;
        sim.vis.showStreamlines = false;
        sim.vis.gridDisplay = 0;
        sim.vis.vectorDisplay = 0;
    }

    void build(Simulation& sim) override {
        sim.useGasBox({0.5f, 0.0f, 0.5f}); // floor at y = 0
        RigidWorld& rigid = sim.rigid;
        ParticleSystem& particles = sim.particles;
        const AABB d = sim.grid.domain();
        // Bodies first (the water is not placed inside them).
        rigid.addBox({0.05f, 0.3f, 0.1f}, Vector3(0.08f), Quaternion::fromAxisAngle({0, 1, 0}, 0.3f), 500.0f, {0.9f, 0.6f, 0.25f});
        rigid.addBox({0.4f, 0.3f, -0.15f}, {0.12f, 0.03f, 0.08f}, Quaternion(), 600.0f, {0.75f, 0.5f, 0.3f}); // plank
        rigid.addSphere({0.7f, 0.45f, 0.15f}, 0.06f, 3000.0f, {0.6f, 0.6f, 0.65f});                          // sinks
        rigid.addCompound(teapotShape(), {-0.15f, 0.35f, -0.15f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.5f), 400.0f,
                          {0.8f, 0.55f, 0.85f});
        TriMesh foam = primitives::box(Vector3(0.06f));
        foam.translate({0.3f, 0.45f, 0.2f});
        particles.addSoftBody(foam, 150.0f, 0.4f, {0.3f, 0.75f, 0.95f});
        // A flag on a pole downwind: its left edge is fixed, the wind makes it flutter.
        ClothMaterial flag;
        flag.areaDensity = 0.15f;
        flag.bendCompliance = 1e-2f;
        particles.addCloth({0.65f, 1.0f, 0.0f}, {0.3f, 0, 0}, {0, -0.22f, 0}, flag, 64, {0.9f, 0.25f, 0.25f});
        // Water: a pool over the whole floor and a column at the upwind end that collapses.
        particles.addBlock(AABB(d.lo, {d.hi.x, 0.16f, d.hi.z}));
        particles.addBlock(AABB({d.lo.x, 0.16f, d.lo.z}, {d.lo.x + 0.4f, 0.6f, d.hi.z}));
    }
};

void addHydroSamples(std::vector<SampleEntry>& out) {
    out.push_back({Preset::Hydro, "Гидродинамика", "Гидродинамика: вода + воздух + тела", [] { return std::unique_ptr<Scene>(new HydroScene); }});
}

} // namespace rf
