// Fire: a gas burner on the floor of a room-corner box open at the top, a cotton curtain that
// catches fire from it, and bodies beside the flame.
#include "samples/Samples.h"
#include "samples/Models.h"

namespace rf {

class FireScene : public Scene {
public:
    void configure(Simulation& sim) override {
        // A gas burner on the floor of a room-corner box open at the top; a cotton curtain hangs
        // next to it. The flame heats the curtain's lower edge until it ignites; the burning fabric
        // gives off fuel gas that burns in the air, so the fire climbs the curtain, which chars
        // through and falls apart. Temperatures in kelvin above ambient (Combustion).
        sim.grid.params.domainSize = {1.2f, 1.6f, 1.0f};
        sim.grid.params.resolutionX = 40;                // dx = 3 cm -> 40 x 53 x 33 cells
        sim.grid.params.inflowSpeed = 1.0f;              // reference scale only
        for (auto& b : sim.grid.params.bc) b = BoundaryType::Wall;
        sim.grid.params.bc[3] = BoundaryType::Outflow;   // open top: the hot gas and the smoke leave
        sim.grid.params.smokeRake = false;
        sim.grid.params.smokeDissipation = 0.15f;
        sim.grid.params.vorticityConfinement = 2.0f;     // small eddies of the flame lost on the grid
        sim.grid.params.pressureTolerance = 1e-3f;
        sim.grid.combustion.enabled = true;
        sim.grid.source.enabled = true;                  // the burner: fuel gas, lit by a pilot flame
        sim.grid.source.center = {0.0f, -0.75f, 0.03f}; // its flame licks the curtain's lower edge
        sim.grid.source.radius = 0.06f;
        sim.grid.source.fuel = 1.0f;
        sim.grid.source.temperature = 400.0f;            // above ignition: burns as it leaves the burner
        sim.grid.source.smoke = 0.0f;
        sim.grid.source.velocity = {0.0f, 0.5f, 0.0f};
        sim.particles.params.clothSpacing = 2.0f;        // 3 cm - the grid spacing
        sim.vis.sliceField = GridField::Temperature;
        sim.vis.showSlice = false;
        sim.vis.showStreamlines = false;
        sim.vis.gridDisplay = 0;
        sim.vis.vectorDisplay = 0;
    }

    void build(Simulation& sim) override {
        sim.useGasBox({0.5f, 0.5f, 0.5f});
        RigidWorld& rigid = sim.rigid;
        ParticleSystem& particles = sim.particles;
        const float floor = sim.grid.domain().lo.y;
        // A cotton curtain on a rod, its lower edge a hand above the burner flame.
        ClothMaterial cotton;
        cotton.areaDensity = 0.2f;
        cotton.flammable = true; // cellulose: Arrhenius pyrolysis (ClothMaterial defaults)
        particles.addCloth({-0.3f, 0.55f, 0.08f}, {0.6f, 0, 0}, {0, -1.1f, 0}, cotton, 16, {0.85f, 0.8f, 0.7f});
        // Bodies beside the fire: a crate and a teapot on the floor, a ball thrown into the plume,
        // a jelly cube falling next to the curtain.
        rigid.addBox({0.38f, floor + 0.1f, -0.25f}, Vector3(0.1f), Quaternion::fromAxisAngle({0, 1, 0}, 0.4f), 500.0f,
                     {0.65f, 0.45f, 0.25f});
        rigid.addCompound(teapotShape(), {-0.35f, floor + 0.14f, -0.25f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.6f), 500.0f,
                          {0.8f, 0.55f, 0.85f});
        const int ball = rigid.addSphere({0.45f, 0.2f, 0.0f}, 0.06f, 300.0f, {0.9f, 0.35f, 0.3f});
        rigid.bodies()[ball].vel = {-1.2f, 1.0f, 0.0f};
        TriMesh jelly = primitives::box(Vector3(0.06f));
        jelly.translate({-0.4f, 0.3f, 0.2f});
        particles.addSoftBody(jelly, SoftMaterial{150.0f, 4e4f, 0.45f}, {0.55f, 0.9f, 0.35f}); // jelly, 40 kPa
    }
};

void addFireSamples(std::vector<SampleEntry>& out) {
    out.push_back({Preset::Fire, "Огонь", "Огонь: горелка, горящая штора, тела", [] { return std::unique_ptr<Scene>(new FireScene); }});
}

} // namespace rf
