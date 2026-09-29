// Soft bodies and cloth in the unified particle solver: a tank scene with a hammock, a curtain
// and jelly bodies, and the same kind of things in a hot gas box.
#include "samples/Models.h"
#include "samples/Samples.h"
#include "samples/SamplesInternal.h"

namespace rf {

// Left: a canvas hammock pinned at its corners with a foam cube and a jelly ball dropped on it.
// Right: a pool with a curtain hanging above it and a floating soft cube.
class SoftClothScene : public Scene {
public:
    void build(Simulation& sim) override {
        sim.useLiquidTank(Simulation::kDefaultTank);
        const AABB& d = sim.particles.domain();
        sim.particles.addBlock(AABB({0.1f, d.lo.y, d.lo.z}, {d.hi.x, 0.22f, d.hi.z}));
        ClothMaterial canvas;
        canvas.areaDensity = 1.5f;
        canvas.tensileStiffness = 0.0f; // inextensible canvas ...
        canvas.bendCompliance = 1e-4f;
        canvas.strengthWarp = canvas.strengthWeft = 0.0f; // ... that does not tear
        sim.particles.addCloth({-0.9f, 0.55f, -0.3f}, {0.6f, 0, 0}, {0, 0, 0.6f}, canvas, 1 | 2 | 4 | 8, {0.85f, 0.72f, 0.4f});
        // Cotton curtain on a rod, sewn from two panels: the vertical seam in the middle holds 30 %
        // of the fabric strength - pull a panel hard sideways and it opens along the seam.
        ClothMaterial cotton;
        cotton.areaDensity = 0.3f;
        cotton.seamColumns = {20}; // 41 columns: the middle
        sim.particles.addCloth({0.4f, 1.1f, -0.3f}, {0, 0, 0.6f}, {0, -0.55f, 0}, cotton, 16 /* on a rod */, {0.75f, 0.3f, 0.35f});
        // Soft bodies: a jelly cube and a softer ball dropped on the trampoline, a firmer block that
        // floats, and a jelly beam laid across two static blocks - it sags between them.
        TriMesh cube = primitives::box(Vector3(0.09f));
        cube.translate({-0.65f, 0.95f, 0.0f});
        sim.particles.addSoftBody(cube, SoftMaterial{150.0f, 3e4f, 0.3f}, {0.3f, 0.75f, 0.95f}); // foam
        TriMesh ball = primitives::sphere(0.08f, 16, 8);
        ball.translate({-0.45f, 1.1f, 0.1f});
        sim.particles.addSoftBody(ball, SoftMaterial{150.0f, 1.5e4f, 0.45f}, {0.55f, 0.9f, 0.35f}); // jelly
        TriMesh floater = primitives::box(Vector3(0.07f));
        floater.translate({0.75f, 0.6f, 0.0f});
        sim.particles.addSoftBody(floater, SoftMaterial{500.0f, 6e4f, 0.3f}, {0.95f, 0.6f, 0.2f});
        // No slender beams here: shape matching on clusters has no bending stiffness across the
        // clusters, so a long thin jelly bar droops like a rope whatever its stiffness (the FEM
        // soft body with a Young's modulus is the next step; a clamped beam is its test).
    }
};

// The closed smoke box with a hot source, and rigid bodies, soft bodies and cloth in the rising gas.
class GasSoftClothScene : public Scene {
public:
    void configure(Simulation& sim) override {
        configureClosedSmokeBox(sim, true);
        // Cloth in the gas scene only flutters in the flow (no bodies thrown at it): the coarser
        // cloth spacing 2r - four times fewer particles - is enough.
        sim.particles.params.clothSpacing = 2.0f;
    }
    void build(Simulation& sim) override {
        sim.useGasBox({0.5f, 0.5f, 0.5f});
        const float floor = sim.grid.domain().lo.y;
        // Rigid: a teapot, a box and a ball on the floor beside the hot source.
        sim.rigid.addCompound(teapotShape(), {-0.45f, floor + 0.16f, 0.3f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.6f), 500.0f,
                              {0.8f, 0.55f, 0.85f});
        sim.rigid.addBox({0.45f, floor + 0.12f, 0.35f}, Vector3(0.12f), Quaternion::fromAxisAngle({0, 1, 0}, 0.4f), 400.0f,
                         {0.3f, 0.7f, 0.9f});
        sim.rigid.addSphere({0.4f, floor + 0.1f, -0.3f}, 0.1f, 500.0f, {0.9f, 0.35f, 0.3f});
        // Soft: a foam cube and a jelly ball dropped from above.
        TriMesh cube = primitives::box(Vector3(0.08f));
        cube.translate({-0.2f, 0.7f, 0.35f});
        sim.particles.addSoftBody(cube, SoftMaterial{150.0f, 6e4f, 0.3f}, {0.3f, 0.75f, 0.95f});
        TriMesh ball = primitives::sphere(0.08f, 16, 8);
        ball.translate({0.35f, 0.9f, -0.25f});
        sim.particles.addSoftBody(ball, SoftMaterial{150.0f, 2e4f, 0.45f}, {0.55f, 0.9f, 0.35f});
        // Cloth: a silk handkerchief high above the hot source: it floats down into the plume, which
        // holds it up (terminal speed ~0.7 m/s, the plume rises at 2-3 m/s) ...
        ClothMaterial silk;
        silk.areaDensity = 0.04f;
        silk.bendCompliance = 1e-2f;
        sim.particles.addCloth({-0.2f, 0.45f, -0.2f}, {0.4f, 0, 0}, {0, 0, 0.4f}, silk, 0, {0.95f, 0.85f, 0.5f});
        // ... and a cotton curtain on a rod beside it, sewn from two panels.
        ClothMaterial cotton;
        cotton.seamColumns = {13}; // 0.8 m at 3 cm: 28 columns, the middle
        sim.particles.addCloth({-0.45f, 1.05f, -0.45f}, {0.8f, 0, 0}, {0, -0.9f, 0}, cotton, 16, {0.75f, 0.3f, 0.35f});
    }
};

void addSoftBodySamples(std::vector<SampleEntry>& out) {
    out.push_back({Preset::SoftCloth, "Мягкие тела", "Мягкие тела и ткань (единый решатель частиц)", [] { return std::unique_ptr<Scene>(new SoftClothScene); }});
    out.push_back({Preset::GasSoftCloth, "Мягкие тела", "Газ + мягкие тела + ткань + твёрдые тела", [] { return std::unique_ptr<Scene>(new GasSoftClothScene); }});
}

} // namespace rf
