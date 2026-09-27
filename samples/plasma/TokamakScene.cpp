// The tokamak scene (TokamakScene.h): a ring of plasma in a toroidal vessel, built on the MHD
// solver of the engine.
#include "samples/plasma/TokamakScene.h"

#include "core/Format.h"
#include "scene/Simulation.h"

#include <cmath>

namespace rf {

void TokamakScene::configure(Simulation& sim) {
    // A ring of plasma in a toroidal vessel (Tokamak.h): the coils' toroidal field, the plasma
    // current's poloidal field and Shafranov's vertical field hold it; the field lines wind
    // with the safety factor q. With q_a = 0.7 (default) the column kinks into a helix within
    // a couple of seconds; set tokamak.safetyFactorEdge above 1 (or below the wall limit 0.4)
    // and reset to see it stay a ring. Scaled to millitesla: the Alfven speed is a few m/s, so
    // the kink takes seconds instead of microseconds.
    const Tokamak& t = tokamak; // R0 = 0.4 m, a = 0.12 m, b = 0.24 m, B0 = 2.5 mT, q_a = 0.7
    sim.obstacle.shape = ObstacleShape::Custom;
    sim.obstacle.customMesh = std::make_shared<TriMesh>(primitives::torus(t.majorRadius, t.vesselRadius));
    sim.obstacle.customName = "торовый сосуд";
    sim.obstacle.size = 2.0f * (t.majorRadius + t.vesselRadius); // the mesh as built
    sim.grid.params.domainSize = {1.4f, 0.6f, 1.4f};
    sim.grid.params.resolutionX = 56; // dx = 2.5 cm: 5 cells across the current channel's radius
    sim.grid.params.inflowSpeed = 0.0f;
    for (BoundaryType& b : sim.grid.params.bc) b = BoundaryType::Wall;
    sim.grid.params.fluidDensity = 1.0f;
    sim.grid.params.pressureTolerance = 1e-3f;
    sim.grid.params.smokeRake = false;
    sim.grid.params.wallFriction = false; // no boundary layer of a gas: the plasma slips along the wall
    sim.grid.magnetic.enabled = true;
    sim.grid.magnetic.conductivity = 1e10f;      // the current outlives the scene: tau = a^2 / (5.8 eta) ~ 30 s
    sim.grid.magnetic.numericalDissipation = 0.0f; // only the Lax-Wendroff amount: f |u| dx would eat the current
    sim.vis.sliceField = GridField::CurrentDensity;
    sim.vis.sliceAxis = 2;
    sim.vis.slicePosition = 0.5f; // the poloidal cross-section at z = 0
    sim.vis.showSlice = false;
    sim.vis.showStreamlines = false;
    sim.vis.surfacePressure = false;
    sim.vis.gridDisplay = 0;
    sim.vis.vesselGlass = true;
}

void TokamakScene::build(Simulation& sim) {
    // The torus is a vessel: the plasma fills it and everything outside is its conducting
    // wall. The mesh is only drawn (as glass), not voxelised as a body.
    sim.grid.vessel = [t = tokamak](const Vector3& x) { return t.inside(x); };
    sim.useGasBox({0.5f, 0.5f, 0.5f});
    // The coils' fields (toroidal B0 R0 / R and the vertical field) are current-free: the
    // background B0. The plasma current's poloidal field is the evolving part B1. The glowing
    // plasma sits where the current flows.
    const Tokamak t = tokamak;
    sim.grid.magnetic.setBackgroundFromPotential([t](const Vector3& x) { return t.coilPotential(x); });
    sim.grid.magnetic.addFromPotential([t](const Vector3& x) { return t.plasmaPotential(x); });
    // The resistive "vacuum" between the channel and the wall (see Tokamak.h).
    const float etaPlasma = sim.grid.magnetic.resistivity();
    sim.grid.magnetic.setResistivityMap([t, etaPlasma](const Vector3& x) { return t.resistivityAt(x, etaPlasma); });
    sim.grid.setTracer([t](const Vector3& x) { return t.tracer(x); });
    bv_ = t.verticalFieldStrength();
    shift_ = 0;
    controlTime_ = 0;
}

void TokamakScene::afterStep(Simulation& sim) {
    // Radial position control, as the vertical-field coils of a real machine: a PD law on the
    // measured outward shift of the ring (the n = 0 part of its current centroid) around the
    // vertical field of the equilibrium in the shell (Tokamak.h). Without it the ring, never
    // quite in the equilibrium of the formulas on a grid, swings in and out for seconds, and the
    // swing's flow eats the current through the dissipation of the induction step.
    const Tokamak& t = tokamak;
    if (!t.positionControl || !t.verticalField) return;
    GasSolver& grid = sim.grid;
    const float shift = t.measuredShift(grid.magnetic, grid.dx());
    const float dt = grid.time() - controlTime_;
    const float rate = controlTime_ > 0 && dt > 1e-5f ? (shift - shift_) / dt : 0.0f;
    shift_ = shift;
    controlTime_ = grid.time();
    float gain, damping;
    t.controlGains(grid.params.fluidDensity, gain, damping);
    const float bv = t.verticalFieldStrength() + gain * shift + damping * rate;
    if (std::fabs(bv - bv_) < 0.002f * std::fabs(t.verticalFieldStrength())) return; // unchanged: keep the field
    bv_ = bv;
    const Tokamak tc = t;
    grid.magnetic.setBackgroundFromPotential([tc, bv](const Vector3& x) { return tc.coilPotential(x, bv); });
}

void TokamakScene::describe(const Simulation& sim, RenderSnapshot& s) const {
    const Tokamak& t = tokamak;
    const MagneticField& m = sim.grid.magnetic;
    const float dx = sim.grid.dx();
    s.info.push_back({"Ток плазмы I_p", format("%.0f А (задано %.0f А)", t.measuredCurrent(m, dx), t.plasmaCurrent())});
    s.info.push_back({"Запас устойчивости q(0) / q(a)", format("%.2f / %.2f", t.safetyFactor(0), t.safetyFactor(t.minorRadius))});
    s.info.push_back({"Кинк m = 1, n = 1",
                      t.kinkUnstable() ? format("растёт: %.2f < q(a) < 1, γ = %.2f 1/с", t.wallLimit(), t.kinkGrowthRate(sim.grid.params.fluidDensity))
                      : t.safetyFactorEdge >= 1 ? std::string("устойчив: q(a) ≥ 1 (предел Крускала–Шафранова)")
                                                : format("устойчив: стенка держит шнур при q(a) < %.2f", t.wallLimit())});
    s.info.push_back({"Вертикальное поле B_v", format("%.3f мТл (равновесие в оболочке %.3f)", bv_ * 1000,
                                                       t.verticalFieldStrength() * 1000)});
    s.info.push_back({"Сдвиг кольца наружу", format("%.1f мм (без B_v по Шафранову %.1f мм)", t.measuredShift(m, dx) * 1000,
                                                     t.equilibriumShift() * 1000)});
    s.plots.push_back({"Амплитуда кинка, мм", t.kinkAmplitude(m, dx) * 1000});
}

void TokamakScene::fieldLineSeeds(const Simulation&, std::vector<Vector3>& seeds) const {
    // Seeds in the poloidal plane phi = 0 at a few minor radii: the lines wind around the
    // torus and show the twist, q(r) toroidal turns per poloidal turn.
    const Tokamak& t = tokamak;
    for (float f : {0.35f, 0.7f, 1.0f, 1.5f})
        for (int q = 0; q < 4; ++q) {
            const float th = 0.5f * kPi * float(q);
            seeds.push_back(t.centre + Vector3(t.majorRadius + f * t.minorRadius * std::cos(th), f * t.minorRadius * std::sin(th), 0.0f));
        }
}

std::vector<SceneParam> TokamakScene::params() const {
    std::vector<SceneParam> p(4);
    p[0] = {"Запас устойчивости q(a)", tokamak.safetyFactorEdge, 0.2f, 8.0f, 0.1f, 2, false};
    p[1] = {"Тороидальное поле B0, мТл", tokamak.toroidalField * 1000.0f, 0.5f, 20.0f, 0.5f, 1, false};
    p[2] = {"Вертикальное поле Шафранова", tokamak.verticalField ? 1.0f : 0.0f, 0.0f, 1.0f, 1.0f, 0, true};
    p[3] = {"Затравка кинка", tokamak.seedDisplacement, 0.0f, 0.2f, 0.01f, 2, false};
    p[0].tip = "Ток плазмы I_p = 2π a² B0 / (μ0 R0 q(a)). Между 2a²/(a²+b²) = 0.4 и 1 шнур скручивается в винт — "
               "кинк-неустойчивость (Крускал–Шафранов); выше 1 держит натяжение линий, ниже 0.4 — стенка";
    p[1].tip = "Поле катушек на магнитной оси; B_φ = B0 R0 / R. Скорость Альфвена растёт с ним, шаг по времени падает";
    p[2].tip = "Держит кольцо от расширения по большому радиусу; без него его держат только токи изображения в стенке";
    p[3].tip = "Винтовое смещение шнура в начале (m = 1, n = 1), доля малого радиуса a";
    return p;
}

void TokamakScene::setParam(int index, float value) {
    switch (index) {
    case 0: tokamak.safetyFactorEdge = value; break;
    case 1: tokamak.toroidalField = value / 1000.0f; break; // millitesla -> tesla
    case 2: tokamak.verticalField = value > 0.5f; break;
    case 3: tokamak.seedDisplacement = value; break;
    default: break;
    }
}

} // namespace rf
