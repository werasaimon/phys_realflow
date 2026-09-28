// Gas on the MAC grid and fire: uniform flow stays uniform, drag and lift against the
// correlations, smoke mass in a closed box, two-way coupling with bodies, cloth and liquid, and
// the combustion model against its exact solutions (T^4 cooling, Fourier conduction, Arrhenius
// pyrolysis).
#include "TestRunner.h"
#include "Tests.h"

#include "core/Probe.h"

void testGasBodies() {
    // 1) Preset: bodies fall through the hot plume into the closed box - the gas stays
    //    divergence-free around the moving boundaries, the bodies come to rest on the floor.
    {
        Simulation sim;
        loadSample(sim, Preset::SmokeBodies);
        float worstRel = 0, gasMax = 0;
        const float y0 = sim.rigid.bodies()[0].pos.y;
        for (int f = 0; f < 60; ++f) sim.stepFrame();
        CHECK(std::fabs(sim.rigid.bodies()[0].pos.y - y0) < 1e-6f, "bodies must hang still until the plume has risen");
        for (int f = 60; f < 240; ++f) {
            sim.stepFrame();
            if (f < 5) continue; // the gas is still at rest: the relative measure is meaningless
            worstRel = std::max(worstRel, sim.grid.maxDivergence() * sim.grid.dx() / std::max(sim.grid.maxVelocity(), 1e-3f));
            gasMax = std::max(gasMax, sim.grid.maxVelocity());
        }
        float vmax = 0, ylo = 1e9f;
        for (const RigidBody& b : sim.rigid.bodies()) {
            vmax = std::max(vmax, length(b.vel));
            ylo = std::min(ylo, b.pos.y);
        }
        std::printf("  smoke + bodies: max gas speed %.2f m/s, worst div*dx/U %.1e, bodies at rest (max |v| %.3f), %d solid cells\n",
                    gasMax, worstRel, vmax, sim.grid.movingSolidCells());
        CHECK(worstRel < 5e-3f, "moving solids break incompressibility: %e", worstRel);
        CHECK(gasMax > 1.0f, "falling bodies must stir the gas (max %f m/s)", gasMax);
        CHECK(vmax < 0.1f && ylo > sim.grid.domain().lo.y, "bodies not at rest in the box (v %f, y %f)", vmax, ylo);
        CHECK(sim.grid.movingSolidCells() > 100, "bodies not voxelised (%d cells)", sim.grid.movingSolidCells());
    }
    // 2) A heavy block pushed through still gas: the gas ahead is carried along, flows back around
    //    the sides, and the pressure force opposes the motion.
    {
        Simulation sim;
        loadSample(sim, Preset::SmokeSphere);
        sim.grid.source.enabled = false;
        sim.grid.params.heatBuoyancy = sim.grid.params.smokeBuoyancy = 0;
        sim.reset();
        sim.rigid.params.gravity = Vector3(0.0f);
        int b = sim.rigid.addBox({-0.4f, 0.0f, 0.0f}, Vector3(0.15f), Quaternion(), 5000.0f, Vector3(1));
        sim.rigid.bodies()[b].vel = {1.0f, 0, 0};
        for (int f = 0; f < 20; ++f) sim.stepFrame();
        const RigidBody& B = sim.rigid.bodies()[b];
        float ahead = sim.grid.velocityAt(B.pos + Vector3(0.2f, 0, 0)).x;
        float beside = sim.grid.velocityAt(B.pos + Vector3(0, 0.25f, 0)).x;
        Vector3 F = sim.grid.movingForces()[0];
        std::printf("  block at 1 m/s: gas ahead %.2f m/s, beside %.2f m/s, drag (%.3f, %.3f, %.3f) N\n", ahead, beside, F.x, F.y, F.z);
        CHECK(ahead > 0.5f, "gas ahead of the block must be pushed along (%f)", ahead);
        CHECK(beside < 0.0f, "gas must flow back around the block (%f)", beside);
        CHECK(F.x < 0 && std::fabs(F.y) < std::fabs(F.x) && std::fabs(F.z) < std::fabs(F.x), "drag must oppose the motion");
    }
    // 3) Two-way: in a dense gas a light sphere falls clearly slower than with one-way coupling
    //    (drag + buoyancy of the displaced gas); one-way it falls freely.
    float vy[2];
    for (int on = 0; on < 2; ++on) {
        Simulation sim;
        loadSample(sim, Preset::SmokeSphere);
        sim.grid.source.enabled = false;
        sim.grid.params.heatBuoyancy = sim.grid.params.smokeBuoyancy = 0;
        sim.grid.params.fluidDensity = 50.0f;
        sim.reset();
        sim.gasPushesBodies = on != 0;
        int b = sim.rigid.addSphere({0, 0.9f, 0}, 0.12f, 200.0f, Vector3(1));
        for (int f = 0; f < 30; ++f) sim.stepFrame();
        vy[on] = sim.rigid.bodies()[b].vel.y;
    }
    float freeFall = -9.81f * 0.5f * std::exp(-0.02f * 0.5f); // with the bodies' linear damping
    std::printf("  sphere rho 200 in gas rho 50 after 0.5 s: vy %.2f (one-way) vs %.2f m/s (two-way), free fall %.2f\n", vy[0], vy[1],
                freeFall);
    CHECK(std::fabs(vy[0] - freeFall) < 0.05f, "one-way coupling must not slow the body (%f vs %f)", vy[0], freeFall);
    CHECK(vy[1] > vy[0] + 1.0f, "gas must slow the sphere: %f vs %f", vy[1], vy[0]);
}

void testSurfaceLoads() {
    // Sphere in the tunnel: loads on every triangle of the real mesh.
    Simulation sim;
    loadSample(sim, Preset::TunnelSphere);
    sim.grid.params.resolutionX = 64;
    sim.reset();
    while (sim.grid.time() < 1.2f) sim.stepFrame();
    const SurfaceLoads& L = sim.surfaceLoads();
    CHECK(L.triangles.size() == sim.obstacleMesh().triangles.size(), "one load per triangle (%zu)", L.triangles.size());
    // Stagnation point: the triangle facing the flow must see Cp ~ +1 (Bernoulli: p0 - p = q).
    const TriangleLoad* front = &L.triangles[0];
    for (const TriangleLoad& t : L.triangles)
        if (t.normal.x < front->normal.x) front = &t;
    // Pressure force: triangles vs the solver's own integral over the voxel faces.
    Vector3 voxel = sim.grid.bodyForce() - sim.grid.frictionForce();
    float rel = std::fabs(L.pressureForce.x - voxel.x) / std::max(std::fabs(voxel.x), 1e-6f);
    float r = 0.25f, area = kPi * r * r;
    std::printf("  sphere, %zu triangles (wetted %.4f m2, exact %.4f): stagnation Cp %.2f; Cd %.3f = pressure %.3f + friction %.3f; "
                "pressure force triangles %.2f N vs voxel faces %.2f N; Cl %.3f\n",
                L.triangles.size(), L.wettedArea, 4 * kPi * r * r, front->cp, L.cd, L.cdPressure, L.cdFriction, L.pressureForce.x,
                voxel.x, L.cl);
    (void)area;
    CHECK(front->cp > 0.8f && front->cp < 1.2f, "stagnation Cp %f (expected ~1)", front->cp);
    CHECK(rel < 0.3f, "triangle and voxel pressure forces disagree by %.0f%%", rel * 100);
    CHECK(std::fabs(L.cl) < 0.05f, "sphere must have no lift: Cl %f", L.cl);
    CHECK(L.cdFriction > 0 && L.cdFriction < 0.1f * L.cd, "friction must be a small part of the sphere drag");
    // CSV export: header + one row per triangle.
    std::string path = "surface_loads_test.csv";
    CHECK(saveSurfaceLoadsCsv(path, L), "CSV not written");
    FILE* f = std::fopen(path.c_str(), "r");
    int lines = 0;
    for (int ch; f && (ch = std::fgetc(f)) != EOF;) lines += ch == '\n';
    if (f) std::fclose(f);
    std::remove(path.c_str());
    CHECK(lines == int(L.triangles.size()) + 2, "CSV rows %d", lines);
}

void testGasParticles() {
    // Gas + soft bodies + cloth + rigid bodies: the hot plume holds up a silk handkerchief (with the
    // aerodynamic coupling; without it the handkerchief falls to the floor), the gas stays
    // divergence-free around everything, the curtain is not torn by the flow, the soft bodies land.
    float handkerchiefY[2] = {0, 0};
    for (int on = 0; on < 2; ++on) {
        Simulation sim;
        loadSample(sim, Preset::GasSoftCloth);
        sim.gasPushesBodies = on != 0;
        float worstRel = 0;
        for (int f = 1; f <= 120; ++f) {
            sim.stepFrame();
            if (f > 5) worstRel = std::max(worstRel, sim.grid.maxDivergence() * sim.grid.dx() / std::max(sim.grid.maxVelocity(), 1e-3f));
        }
        const ParticleSystem& P = sim.particles;
        const Cloth& h = P.cloths()[0];
        Vector3 c(0.0f);
        for (int i = 0; i < h.width * h.height; ++i) c += P.positions()[h.firstParticle + i];
        handkerchiefY[on] = c.y / float(h.width * h.height);
        bool finite = true;
        for (const Vector3& p : P.positions()) finite &= std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
        Vector3 soft(0.0f);
        for (int i : P.softBodies()[0].particles) soft += P.positions()[i];
        soft /= float(P.softBodies()[0].particles.size());
        if (on) {
            std::printf("  gas + particles: worst div*dx/U %.1e, curtain torn %d, foam cube y %.3f (floor %.2f)%s", worstRel,
                        P.cloths()[1].tornThreads, soft.y, sim.grid.domain().lo.y, "\n");
            CHECK(worstRel < 5e-3f && finite, "gas / particles broke (div %e, finite %d)", worstRel, int(finite));
            CHECK(P.cloths()[1].tornThreads == 0, "the flow tore the curtain (%d threads)", P.cloths()[1].tornThreads);
            CHECK(soft.y < sim.grid.domain().lo.y + 0.15f, "soft body did not land (y %f)", soft.y);
        }
    }
    std::printf("  handkerchief after 2 s: centre y %.2f with the gas drag, %.2f without%s", handkerchiefY[1], handkerchiefY[0], "\n");
    CHECK(handkerchiefY[1] > handkerchiefY[0] + 0.5f, "the plume must hold the handkerchief up (%f vs %f)", handkerchiefY[1],
          handkerchiefY[0]);
}

namespace {

// The x of the flag's centre now (the mean of its particles).
float flagCentreX(const ParticleSystem& P) {
    const Cloth& flag = P.cloths()[0];
    float x = 0;
    for (int i = 0; i < flag.width * flag.height; ++i) x += P.positions()[size_t(flag.firstParticle + i)].x;
    return x / float(flag.width * flag.height);
}

// Takes the collapsing water column out of the hydro scene: the fluid group whose particles start
// above the pool (y > 0.17 m). Without it nothing sets the air in motion - the still reference.
void removeWaterColumn(ParticleSystem& P) {
    for (size_t i = 0; i < P.size(); ++i)
        if (P.phases()[i] == uint8_t(ParticlePhase::Fluid) && P.positions()[i].y > 0.17f) return P.removeGroup(P.groupOf(int(i)));
}

} // namespace

void testHydro() {
    // Water + air + bodies: light bodies float, the heavy ball sinks, the flag streams downwind in
    // the wind; the air stays divergence-free around the moving water.
    // The flag is compared with the flag in STILL air: the same scene without wind and without the
    // water column. The column's collapse stirs the air of the closed box, and that motion is
    // physical - the gas without viscous damping keeps it (advection-reflection loses almost no
    // energy), so "no wind" with the column is not still air; the old check passed only because the
    // numerical dissipation of the old step killed that motion (2026-09-28). The flag swings by a
    // few cm, so its centre is averaged over the last second (frames 120-180) in both runs.
    float flagX[2] = {0, 0};
    for (int windOn = 0; windOn < 2; ++windOn) {
        Simulation sim;
        loadSample(sim, Preset::Hydro);
        sim.grid.params.inflowSpeed = windOn ? 3.0f : 0.0f;
        sim.reset();
        if (!windOn) removeWaterColumn(sim.particles);
        float worstRel = 0;
        for (int f = 1; f <= 180; ++f) {
            sim.stepFrame();
            if (f > 5) worstRel = std::max(worstRel, sim.grid.maxDivergence() * sim.grid.dx() / std::max(sim.grid.maxVelocity(), 1e-3f));
            if (f > 120) flagX[windOn] += flagCentreX(sim.particles) / 60.0f;
        }
        const ParticleSystem& P = sim.particles;
        bool finite = true;
        for (const Vector3& p : P.positions()) finite &= std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
        if (!windOn) continue;
        const auto& B = sim.rigid.bodies();
        std::printf("  hydro: %zu water particles; box y %.2f, plank y %.2f, teapot y %.2f, heavy ball y %.2f; div*dx/U %.1e%s",
                    P.fluidCount(), B[0].pos.y, B[1].pos.y, B[3].pos.y, B[2].pos.y, worstRel, "\n");
        CHECK(finite && worstRel < 5e-3f, "hydro scene broke (div %e, finite %d)", worstRel, int(finite));
        // Floating = clearly above where each would lie on the bottom (box 0.08, plank 0.03, teapot
        // ~0.15), even while the wave of the collapsed column still rocks them.
        CHECK(B[0].pos.y > 0.12f && B[1].pos.y > 0.08f && B[3].pos.y > 0.18f, "light bodies must float (%f %f %f)", B[0].pos.y,
              B[1].pos.y, B[3].pos.y);
        CHECK(B[2].pos.y < 0.1f, "the heavy ball must sink (y %f)", B[2].pos.y);
    }
    std::printf("  flag (pole at x 0.65), centre x over the last second: %.4f in the wind, %.4f in still air%s", flagX[1], flagX[0], "\n");
    // The drag acts on the flag's real area (each particle its share of the sheet, as its mass - the
    // grid spacing squared overstated it by 30 %).
    CHECK(flagX[1] > 0.75f && flagX[1] > flagX[0] + 0.03f, "the wind must stream the flag out (%f vs %f)", flagX[1], flagX[0]);
}

void testGridUniform() {
    GasSolver g;
    g.params.domainSize = {2, 1, 1};
    g.params.resolutionX = 32;
    g.reset({0, 0, 0}, nullptr);
    for (int i = 0; i < 20; ++i) g.step(0.05f);
    Vector3 v = g.velocityAt({1.0f, 0.5f, 0.5f});
    CHECK(std::fabs(v.x - g.params.inflowSpeed) < 0.05f && std::fabs(v.y) < 0.05f, "uniform flow %f %f", v.x, v.y);
}

void testGridSphere() {
    GasSolver g;
    g.params.domainSize = {4, 2, 2};
    g.params.resolutionX = 64;
    TriMesh sp = primitives::sphere(0.25f);
    MeshBVH bvh;
    bvh.build(sp);
    g.reset({-1.2f, -1, -1}, &bvh);
    CHECK(g.hasObstacle(), "sphere voxelised");
    float area = g.frontalArea(), ae = kPi * 0.25f * 0.25f;
    CHECK(std::fabs(area - ae) / ae < 0.25f, "frontal area %f vs %f", area, ae);
    double t = 0;
    while (t < 1.2) t += g.step(0.05f);
    float cd = g.dragCoefficientAvg();
    std::printf("  sphere: Cd(avg)=%.3f  Cl=%.3f  iters=%d  residual=%.1e\n", cd, g.liftCoefficientAvg(),
                g.lastPressureIterations(), g.lastResidual());
    // Coarse inviscid solver: expect the right order of magnitude (experiment: ~0.4-0.5 subcritical).
    CHECK(cd > 0.1f && cd < 1.5f, "sphere Cd %f", cd);
    CHECK(std::fabs(g.liftCoefficientAvg()) < 0.3f, "sphere Cl %f", g.liftCoefficientAvg());
    CHECK(g.lastResidual() < 1e-3f, "pressure residual %f", g.lastResidual());
}

void testGridWingLift() {
    auto liftAt = [](float aoa) {
        Simulation sim;
        loadSample(sim, Preset::TunnelWing);
        sim.grid.params.resolutionX = 64;
        sim.obstacle.angleOfAttackDeg = aoa;
        sim.rebuildObstacle();
        while (sim.grid.time() < 1.0f) sim.stepFrame();
        return sim.grid.liftCoefficientAvg();
    };
    float c0 = liftAt(0.0f), c8 = liftAt(8.0f);
    std::printf("  wing NACA2412: Cl(0)=%.3f  Cl(8)=%.3f\n", c0, c8);
    CHECK(c8 > c0 + 0.1f, "lift must grow with angle of attack: %f -> %f", c0, c8);
}

void testSmokeClosedBox() {
    // Closed box (all walls), hot smoky sphere: gas must rise, stay divergence-free and bounded.
    Simulation sim;
    loadSample(sim, Preset::SmokeSphere);
    GasSolver& g = sim.grid;
    CHECK(g.nx() == 32 && g.ny() == 48 && g.nz() == 32, "grid %d %d %d", g.nx(), g.ny(), g.nz());
    float smoke1 = 0;
    for (int f = 0; f < 120; ++f) {
        sim.stepFrame();
        if (f == 20) smoke1 = g.totalSmoke();
    }
    float smoke2 = g.totalSmoke();
    Vector3 above = g.velocityAt(g.source.center + Vector3(0, 0.4f, 0));
    std::printf("  smoke box: t=%.2f s  max|div|=%.2e  smoke %.4f -> %.4f m^3  v_above=%.3f m/s  iters=%d\n",
                g.time(), g.maxDivergence(), smoke1, smoke2, above.y, g.lastPressureIterations());
    CHECK(above.y > 0.05f, "hot gas above the source must rise, vy=%f", above.y);
    CHECK(smoke2 > smoke1, "source keeps adding smoke: %f -> %f", smoke1, smoke2);
    // |div| relative to the flow scale U/dx must be small after the projection.
    float rel = g.maxDivergence() * g.dx() / std::max(g.maxVelocity(), 1e-3f);
    CHECK(rel < 5e-3f, "divergence not removed: max|div| dx / Umax = %e", rel);
    // Walls are impermeable: normal velocity on the boundary faces is zero.
    Vector3 wallTop = g.velocityAt(g.domain().center() + Vector3(0, 0.5f * g.domain().extent().y, 0));
    CHECK(std::fabs(wallTop.y) < 1e-4f, "top wall normal velocity %f", wallTop.y);
    float smin = 1e9f, smax = -1e9f;
    for (float s : g.smoke().d) { smin = std::min(smin, s); smax = std::max(smax, s); }
    CHECK(smin > -1e-4f && smax < 1.0f + 1e-4f, "smoke out of range [%f, %f]", smin, smax);
}

void testDisturbance() {
    GasSolver g;
    g.params.domainSize = {1, 1, 1};
    g.params.resolutionX = 20;
    for (auto& b : g.params.bc) b = BoundaryType::Wall;
    g.params.inflowSpeed = 1.0f;
    g.reset({0, 0, 0}, nullptr);
    Disturbance d;
    d.center = {0.5f, 0.5f, 0.5f};
    d.radius = 0.2f;
    d.velocity = {2.0f, 0, 0};
    d.smoke = 1.0f;
    g.applyDisturbance(d);
    CHECK(g.velocityAt(d.center).x > 1.5f, "velocity at the brush centre %f", g.velocityAt(d.center).x);
    CHECK(g.smoke().sample((d.center / g.dx())) > 0.5f, "smoke at the brush centre");
    g.step(0.02f);
    CHECK(g.maxDivergence() * g.dx() / std::max(g.maxVelocity(), 1e-3f) < 5e-3f, "projected after disturbance");
    // velocityBlend = 0 leaves the flow untouched.
    GasSolver h;
    h.params = g.params;
    h.reset({0, 0, 0}, nullptr);
    d.velocityBlend = 0;
    h.applyDisturbance(d);
    CHECK(length(h.velocityAt(d.center)) < 1e-6f, "blend 0 must not change velocity");
}

void testCombustion() {
    // Radiative cooling is integrated exactly: one big step == many small ones.
    Combustion comb;
    std::vector<uint8_t> solid(1, 0);
    std::vector<float> fuel(1, 0.0f), products(1, 0.0f), T(1, 1500.0f), smoke(1, 0.0f), expansion(1, 0.0f);
    comb.react(fuel, products, T, smoke, expansion, solid, 0.5f);
    const float oneStep = T[0];
    T[0] = 1500.0f;
    for (int i = 0; i < 5000; ++i) comb.react(fuel, products, T, smoke, expansion, solid, 0.5f / 5000);
    std::printf("  radiative cooling from 1500 K over 0.5 s: %.2f K (1 step) vs %.2f K (5000 steps)\n", oneStep, T[0]);
    CHECK(std::fabs(oneStep - T[0]) < 0.5f && oneStep < 1500.0f, "cooling not exact (%f vs %f)", oneStep, T[0]);

    // Oxygen limit: a cell of pure fuel (10 x stoichiometric) burns only what its air allows, and
    // every unit burnt heats the gas by exactly heatRelease (energy conservation).
    comb.radiativeCooling = 0;
    fuel[0] = 10.0f;
    products[0] = 0.0f;
    T[0] = 500.0f;
    double burnt = 0;
    for (int i = 0; i < 600; ++i) burnt += comb.react(fuel, products, T, smoke, expansion, solid, 0.01f);
    std::printf("  rich cell: burnt %.4f of 10 fuel units, products %.4f, T rise %.1f K (heatRelease %.0f)\n", burnt, products[0],
                T[0] - 500.0f, comb.heatRelease);
    CHECK(burnt <= 1.0 + 1e-4 && burnt > 0.99, "oxygen must limit the burning (burnt %f)", burnt);
    CHECK(std::fabs(T[0] - 500.0f - comb.heatRelease * float(burnt)) < 1.0f, "heat release must match the fuel burnt");
    CHECK(std::fabs(fuel[0] - (10.0f - float(burnt))) < 1e-3f, "fuel must be conserved");
    // Cold fuel does not burn.
    fuel[0] = 1.0f; products[0] = 0.0f; T[0] = 0.0f;
    CHECK(comb.react(fuel, products, T, smoke, expansion, solid, 0.1f) == 0.0, "cold fuel must not burn");
}

void testHeatConduction() {
    // Fourier's law in the gas at rest: the heat is conserved and a spot spreads so that its
    // second moment grows as d<r^2>/dt = 6 alpha (3D diffusion), whatever its shape.
    GasSolver g;
    g.params.domainSize = {1, 1, 1};
    g.params.resolutionX = 32;
    g.params.inflowSpeed = 0;
    g.params.smokeRake = false;
    for (auto& b : g.params.bc) b = BoundaryType::Wall;
    g.combustion.enabled = true;
    g.combustion.gravity = 0;
    g.combustion.radiativeCooling = 0;
    g.combustion.thermalDiffusivity = 1e-3f;
    g.reset({-0.5f, -0.5f, -0.5f}, nullptr);
    Disturbance spot;
    spot.radius = 0.15f;
    spot.velocityBlend = 0;
    spot.heat = 1.0f; // 1 K: the diffusivity stays constant to 0.6 %
    g.applyDisturbance(spot);
    auto moments = [&](double& heat, double& r2) {
        heat = r2 = 0;
        const Field3& T = g.temperature();
        for (int k = 0; k < g.nz(); ++k)
            for (int j = 0; j < g.ny(); ++j)
                for (int i = 0; i < g.nx(); ++i) {
                    const Vector3 x = g.origin() + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * g.dx();
                    heat += T.at(i, j, k);
                    r2 += T.at(i, j, k) * length2(x);
                }
        r2 /= heat;
    };
    double h0, m0, h1, m1;
    moments(h0, m0);
    float t = 0;
    while (t < 1.0f - 1e-5f) t += g.step(std::min(0.05f, 1.0f - t));
    moments(h1, m1);
    const double rate = (m1 - m0) / t, expected = 6.0 * g.combustion.thermalDiffusivity;
    std::printf("  gas: heat %.6f -> %.6f, d<r^2>/dt %.3e vs 6 alpha %.3e\n", h0, h1, rate, expected);
    CHECK(std::fabs(h1 / h0 - 1) < 1e-4, "heat must be conserved (%f -> %f)", h0, h1);
    CHECK(std::fabs(rate / expected - 1) < 0.05, "diffusion rate %e, expected %e", rate, expected);

    // Along a cloth: conduction through the threads alone (no gas exchange) conserves the heat and
    // warms the neighbours of a hot patch; nothing decomposes at 100 K above ambient.
    ParticleSystem ps;
    ps.params.clothSpacing = 2.0f;
    ps.reset(AABB({-1, -1, -1}, {1, 1, 1}));
    ClothMaterial cotton;
    cotton.flammable = true;
    cotton.heatTransfer = 0;
    cotton.emissivity = 0;
    ps.addCloth({-0.15f, 0, 0}, {0.3f, 0, 0}, {0, 0, 0.3f}, cotton, 0, Vector3(1));
    Cloth c = ps.cloths()[0];
    const int mid = c.width / 2 + c.width * (c.height / 2);
    c.temperature[mid] = 100.0f;
    const size_t n = c.temperature.size();
    std::vector<float> zero(n, 0.0f), heat(n, 0.0f), fuel(n, 0.0f);
    for (int s = 0; s < 600; ++s) burnCloth(c, zero, zero, 293.0f, 1.0f / 60, heat, fuel);
    double total = 0;
    for (float T : c.temperature) total += T;
    std::printf("  cloth: heat %.4f K-patches (100 at start), hot patch now %.2f K, neighbour %.3f K\n", total, c.temperature[mid],
                c.temperature[mid + 1]);
    CHECK(std::fabs(total - 100.0) < 1e-2, "cloth conduction must conserve heat (%f)", total);
    CHECK(c.temperature[mid] < 100.0f && c.temperature[mid + 1] > 0.0f, "heat must flow along the fabric");
    CHECK(c.unburnt[mid] == 1.0f && fuel[mid] == 0.0f, "no pyrolysis at 100 K above ambient");
}

void testPyrolysis() {
    // A cotton patch under a constant radiant flux (as in a cone calorimeter): Arrhenius pyrolysis
    // with the energy balance. A strong flux chars it through within seconds; a weak one, below the
    // critical flux of cotton (~10-15 kW/m^2), leaves it intact - there is no ignition temperature
    // in the model, the threshold comes out of the kinetics and the heat balance.
    auto timeToChar = [](float flux, float tMax) {
        ParticleSystem ps;
        ps.params.clothSpacing = 2.0f;
        ps.reset(AABB({-1, -1, -1}, {1, 1, 1}));
        ClothMaterial cotton;
        cotton.flammable = true;
        cotton.areaDensity = 0.2f;
        ps.addCloth({0, 0, 0}, {0.03f, 0, 0}, {0, 0, 0.03f}, cotton, 0, Vector3(1));
        Cloth c = ps.cloths()[0];
        const size_t n = c.temperature.size();
        std::vector<float> gas(n, 0.0f), q(n, flux), heat(n, 0.0f), fuel(n, 0.0f);
        const float dt = 1.0f / 60;
        for (float t = 0; t < tMax; t += dt) {
            burnCloth(c, gas, q, 293.0f, dt, heat, fuel);
            if (c.unburnt[0] < 0.5f) return std::make_pair(t, c.temperature[0]);
        }
        return std::make_pair(-1.0f, c.temperature[0]);
    };
    const auto strong = timeToChar(100e3f, 20.0f), weak = timeToChar(8e3f, 30.0f);
    std::printf("  100 kW/m^2: half decomposed after %.2f s (at %.0f K); 8 kW/m^2: %s, patch at %.0f K\n", strong.first,
                strong.second + 293, weak.first < 0 ? "intact after 30 s" : "decomposed", weak.second + 293);
    CHECK(strong.first > 0.5f && strong.first < 10.0f, "cotton under 100 kW/m^2 must char within seconds (%f)", strong.first);
    CHECK(weak.first < 0, "cotton under 8 kW/m^2 must not ignite");
}

void testFireScene() {
    // Burner + cotton curtain: the flame ignites the curtain, the fire climbs it and burns through
    // threads; the flame stays at physical temperatures (oxygen-limited) and the scene stays finite.
    Simulation sim;
    loadSample(sim, Preset::Fire);
    const ParticleSystem& P = sim.particles;
    float ignition = -1, Tmax = 0, hrrMax = 0, topFire = -1e9f;
    for (int f = 1; f <= 300; ++f) {
        sim.stepFrame();
        const Cloth& c = P.cloths()[0];
        for (int k = 0; k < c.width * c.height; ++k)
            if (c.unburnt[k] < 0.9f) {
                if (ignition < 0) ignition = f / 60.0f;
                topFire = std::max(topFire, P.positions()[c.firstParticle + k].y);
            }
        for (float T : sim.grid.temperature().d) Tmax = std::max(Tmax, T);
        hrrMax = std::max(hrrMax, sim.grid.heatReleaseRate());
    }
    bool finite = true;
    for (const Vector3& x : P.positions()) finite &= std::isfinite(x.x + x.y + x.z);
    const Cloth& c = P.cloths()[0];
    std::printf("  curtain ignites at %.2f s, fire reaches y %.2f (rod at 0.55); %d threads burnt through; flame Tmax %.0f K, "
                "peak power %.0f kW\n",
                ignition, topFire, c.burntThreads, Tmax + sim.grid.combustion.ambientTemperature, hrrMax / 1000);
    CHECK(finite, "fire scene produced NaN");
    CHECK(ignition > 0 && ignition < 4.0f, "the burner flame must ignite the curtain (%f)", ignition);
    CHECK(topFire > 0.3f, "the fire must climb the curtain (%f)", topFire);
    CHECK(c.burntThreads > 20, "threads must burn through (%d)", c.burntThreads);
    CHECK(Tmax + 293 > 1100 && Tmax + 293 < 2300, "flame temperature out of the physical range (%f K)", Tmax + 293);
    CHECK(hrrMax < 2e6f, "heat release ran away (%f W)", hrrMax);
}

// Validation against a published benchmark: the vortex street behind a circular cylinder at
// Re = 100 sheds at the Strouhal number St = f D / U = 0.164 in an unbounded flow (Williamson
// 1996; Roshko 1954: 0.198 (1 - 19.7 / Re) = 0.159). Walls raise it: at 12.5 % blockage (D / H)
// about 0.17, at 25 % about 0.20 (Sahin & Owens 2004) - the tunnel here is 8 D high. The
// cylinder spans the whole depth (a quasi-2D flow, half a diameter deep), the lift coefficient
// oscillates at the shedding frequency; a small nudge in the wake at the start saves the seconds
// the instability would take to grow out of round-off.
void testCylinderStrouhal() {
    GasSolver g;
    g.params.domainSize = {4, 4, 0.5f};
    g.params.resolutionX = 96; // dx = 4.2 cm, 12 cells across the cylinder
    const float D = 0.5f, U = 1.0f, Re = 100.0f;
    g.params.inflowSpeed = U;
    g.params.kinematicViscosity = U * D / Re;
    g.params.smokeRake = false;
    TriMesh cyl = primitives::cylinder(0.5f * D, 0.6f);
    MeshBVH bvh;
    bvh.build(cyl);
    g.reset({-1.2f, -2, -0.25f}, &bvh);
    Disturbance nudge;
    nudge.center = {0.6f, 0.15f, 0.0f};
    nudge.radius = 0.2f;
    nudge.velocity = {U, 0.3f * U, 0.0f};
    g.applyDisturbance(nudge);
    std::vector<std::pair<float, float>> cl; // (t, Cl)
    double t = 0, cdSum = 0;
    int cdCount = 0;
    while (t < 40.0) {
        t += g.step(0.1f);
        cl.push_back({float(t), g.liftCoefficient()});
        if (t > 20.0) { cdSum += g.dragCoefficient(); ++cdCount; }
    }
    // Shedding frequency from the zero crossings of Cl over the last 20 s, amplitude from its extremes.
    int crossings = 0;
    float first = -1, last = -1, clMax = -1e9f, clMin = 1e9f;
    for (size_t i = 1; i < cl.size(); ++i) {
        if (cl[i].first < 20.0f) continue;
        clMax = std::max(clMax, cl[i].second);
        clMin = std::min(clMin, cl[i].second);
        if ((cl[i - 1].second < 0) != (cl[i].second < 0)) {
            if (first < 0) first = cl[i].first;
            last = cl[i].first;
            ++crossings;
        }
    }
    const float f = crossings > 2 ? float(crossings - 1) / 2.0f / (last - first) : 0.0f, St = f * D / U;
    std::printf("  cylinder Re 100, 12.5%% blockage: St = %.3f (Williamson 0.164 unbounded, ~0.17 confined), Cl amplitude %.2f, Cd mean %.2f "
                "(2D reference ~1.3), %d zero crossings in %.0f s, residual %.1e\n",
                St, 0.5f * (clMax - clMin), cdSum / std::max(1, cdCount), crossings, last - first, g.lastResidual());
    CHECK(crossings >= 6, "no vortex shedding: %d zero crossings of Cl", crossings);
    CHECK(St > 0.145f && St < 0.195f, "Strouhal number %f vs 0.164 (0.17 confined)", St);
    CHECK(0.5f * (clMax - clMin) > 0.05f, "lift oscillation too weak: %f", 0.5f * (clMax - clMin));
}

// Solution verification by grid refinement (Roache 1998; ASME V&V 20): the same problem on three
// grids h, h/2, h/4 against an exact solution, the observed order p = log2(e_h / e_h/2) and the
// grid convergence index GCI = 1.25 |f_fine - f_medium| / (|f_fine| (2^p - 1)) of a scalar of the
// solution (Richardson). Two problems with exact answers:
//  1. heat conduction of a smooth spot in still gas - T(r, t) from the radial heat kernel;
//  2. a smooth smoke blob carried by a uniform flow - the blob shifted by U t.
void testGridConvergence() {
    auto rms = [](double sum2, size_t n) { return std::sqrt(sum2 / double(std::max<size_t>(n, 1))); };
    auto order = [](double eCoarse, double eFine) { return std::log(eCoarse / eFine) / std::log(2.0); };

    // 1) Heat: T0(r) = 0.1 K (1 - r^2/R^2)^2 inside R (the disturbance's shape), alpha = 1e-3 m^2/s,
    //    1 s: the exact solution of a radial profile is the integral against the 3D radial kernel
    //    G(r, r') = r' / (r sqrt(4 pi a t)) [exp(-(r-r')^2/4at) - exp(-(r+r')^2/4at)]. The time step
    //    is the same small one on every grid (4 ms, 250 steps), so the spatial order is what is
    //    measured; 0.1 K keeps the T^1.75 diffusivity constant to 0.06 %.
    {
        const float R = 0.15f, alpha = 1e-3f, tEnd = 1.0f, heat = 0.1f;
        auto exact = [&](double r) {
            const double s = 4.0 * alpha * tEnd, n = 400;
            double sum = 0;
            for (int m = 0; m <= n; ++m) { // Simpson over r' in [0, R]
                const double rp = R * m / n, w = (m == 0 || m == n) ? 1 : (m % 2 ? 4 : 2);
                const double T0 = heat * std::pow(1.0 - rp * rp / (double(R) * R), 2.0);
                sum += w * T0 * rp * (std::exp(-(r - rp) * (r - rp) / s) - std::exp(-(r + rp) * (r + rp) / s));
            }
            return sum * (R / n / 3.0) / (r * std::sqrt(kPi * s));
        };
        double err[3], peak[3], dt[3];
        int res[3] = {16, 32, 64};
        for (int lvl = 0; lvl < 3; ++lvl) {
            GasSolver g;
            g.params.domainSize = {1, 1, 1};
            g.params.resolutionX = res[lvl];
            g.params.inflowSpeed = 0;
            g.params.smokeRake = false;
            for (auto& b : g.params.bc) b = BoundaryType::Wall;
            g.combustion.enabled = true;
            g.combustion.gravity = 0;
            g.combustion.radiativeCooling = 0;
            g.combustion.thermalDiffusivity = alpha;
            g.reset({-0.5f, -0.5f, -0.5f}, nullptr);
            Disturbance spot;
            spot.radius = R;
            spot.velocityBlend = 0;
            spot.heat = heat;
            g.applyDisturbance(spot);
            float t = 0;
            dt[lvl] = 0;
            while (t < tEnd - 1e-5f) { const float h = g.step(std::min(0.004f, tEnd - t)); t += h; dt[lvl] = std::max(dt[lvl], double(h)); }
            const Field3& T = g.temperature();
            double sum2 = 0, pk = 0;
            for (int k = 0; k < g.nz(); ++k)
                for (int j = 0; j < g.ny(); ++j)
                    for (int i = 0; i < g.nx(); ++i) {
                        const Vector3 x = g.origin() + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * g.dx();
                        const double e = T.at(i, j, k) - exact(length(x));
                        sum2 += e * e;
                        pk = std::max(pk, double(T.at(i, j, k)));
                    }
            err[lvl] = rms(sum2, T.d.size());
            peak[lvl] = pk;
        }
        const double p1 = order(err[0], err[1]), p2 = order(err[1], err[2]);
        const double pPeak = std::log(std::fabs((peak[0] - peak[1]) / (peak[1] - peak[2]))) / std::log(2.0);
        const double gci = 1.25 * std::fabs(peak[2] - peak[1]) / (std::fabs(peak[2]) * (std::pow(2.0, pPeak) - 1.0));
        std::printf("  heat spot, exact peak %.4f K: grids 16/32/64, RMS error %.2e / %.2e / %.2e, order %.2f then %.2f; peak %.4f / %.4f / %.4f "
                    "(order %.2f, GCI fine %.2f%%), max dt %.3f s\n",
                    exact(1e-6), err[0], err[1], err[2], p1, p2, peak[0], peak[1], peak[2], pPeak, 100 * gci, dt[2]);
        CHECK(p2 > 1.5, "heat conduction converges at order %.2f, expected ~2", p2);
        CHECK(gci < 0.05, "heat peak GCI %.1f%%", 100 * gci);
    }

    // 2) Advection: smoke = exp(-|x - c|^2 / 2 s^2), s = 0.12 m, in a uniform flow U = 1 m/s for 0.5 s
    //    (no obstacle, so the projection leaves the flow uniform): exact = the blob 0.5 m downstream.
    //    The Courant number is the same on every grid, 0.7, and not an integer: with U dt = 2 dx the
    //    back-traces land on cell centres, the interpolation is exact and the error vanishes - a
    //    property of the semi-Lagrangian step, not of the grid.
    {
        const float U = 1.0f, tEnd = 0.5f, sigma = 0.12f;
        const Vector3 c0(0.5f, 0.5f, 0.5f);
        double err[3], peak[3], dt[3];
        int res[3] = {32, 64, 128};
        for (int lvl = 0; lvl < 3; ++lvl) {
            GasSolver g;
            g.params.domainSize = {2, 1, 1};
            g.params.resolutionX = res[lvl];
            g.params.inflowSpeed = U;
            g.params.smokeRake = false;
            g.reset({0, 0, 0}, nullptr);
            g.setTracer([&](const Vector3& x) { return std::exp(-length2(x - c0) / (2 * sigma * sigma)); });
            float t = 0;
            dt[lvl] = 0;
            const float courant = 0.7f * g.dx() / U;
            while (t < tEnd - 1e-5f) { const float h = g.step(std::min(courant, tEnd - t)); t += h; dt[lvl] = std::max(dt[lvl], double(h)); }
            const Vector3 c1 = c0 + Vector3(U * tEnd, 0, 0);
            const Field3& S = g.smoke();
            double sum2 = 0, pk = 0;
            for (int k = 0; k < g.nz(); ++k)
                for (int j = 0; j < g.ny(); ++j)
                    for (int i = 0; i < g.nx(); ++i) {
                        const Vector3 x = g.origin() + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * g.dx();
                        const double e = S.at(i, j, k) - std::exp(-length2(x - c1) / (2 * sigma * sigma));
                        sum2 += e * e;
                        pk = std::max(pk, double(S.at(i, j, k)));
                    }
            err[lvl] = rms(sum2, S.d.size());
            peak[lvl] = pk;
        }
        const double p1 = order(err[0], err[1]), p2 = order(err[1], err[2]);
        const double pPeak = std::log(std::fabs((peak[0] - peak[1]) / (peak[1] - peak[2]))) / std::log(2.0);
        const double gci = 1.25 * std::fabs(peak[2] - peak[1]) / (std::fabs(peak[2]) * (std::pow(2.0, pPeak) - 1.0));
        std::printf("  smoke blob advected 0.5 m: grids 32/64/128, RMS error %.2e / %.2e / %.2e, order %.2f then %.2f; peak %.4f / %.4f / %.4f "
                    "(exact 1, order %.2f, GCI fine %.2f%%), max dt %.4f / %.4f / %.4f s\n",
                    err[0], err[1], err[2], p1, p2, peak[0], peak[1], peak[2], pPeak, 100 * gci, dt[0], dt[1], dt[2]);
        CHECK(p2 > 0.8, "advection converges at order %.2f", p2);
        CHECK(peak[2] > 0.9, "the blob's peak after 0.5 m is %.3f of 1 on the fine grid", peak[2]);
    }
}

namespace {

// A cube of side 1 m with n cells per side: a sphere of radius 0.2 m in the middle, gas coming in
// at 1 m/s through -x and leaving through +x, walls elsewhere. The first projection from rest is
// the potential flow around the sphere - a full, cold-started pressure solve.
struct PressureRun {
    int iterations = 0;
    double ms = 0;
    GasSolver gas;
};

void firstProjection(PressureRun& run, int n, PressurePreconditioner pre, float tolerance) {
    static MeshBVH sphere; // built once: the same obstacle for every grid
    if (sphere.empty()) sphere.build(primitives::sphere(0.2f));
    GasSolver& g = run.gas;
    g.params.domainSize = {1, 1, 1};
    g.params.resolutionX = n;
    g.params.inflowSpeed = 1.0f;
    g.params.smokeRake = false;
    g.params.wallFriction = false;
    g.params.pressureTolerance = tolerance;
    g.params.maxPressureIterations = 5000;
    g.params.pressurePreconditioner = pre;
    g.reset({-0.5f, -0.5f, -0.5f}, &sphere);
    g.step(0.5f / n);
    run.ms = g.lastPressureMs(); // the projection alone (no magnetic field here)
    run.iterations = g.lastPressureIterations();
}

// The largest difference of the cell velocities of two solvers on the same grid [m/s].
float maxVelocityDifference(const GasSolver& a, const GasSolver& b) {
    float worst = 0;
    for (int k = 0; k < a.nz(); ++k)
        for (int j = 0; j < a.ny(); ++j)
            for (int i = 0; i < a.nx(); ++i)
                worst = std::max(worst, length(a.cellVelocity(i, j, k) - b.cellVelocity(i, j, k)));
    return worst;
}

} // namespace

// The multigrid preconditioner (Multigrid.h, McAdams et al. 2010) against the Jacobi one on the
// same equation: (1) with a tight tolerance both give the same flow; (2) the iteration count of the
// multigrid PCG hardly grows with the grid, while the Jacobi PCG's grows about like n.
void testMultigridPressure() {
    Probe::beginFrame(); // counters to zero: "gas/multigrid fallbacks" must stay 0
    {
        PressureRun jacobi, multigrid;
        firstProjection(jacobi, 48, PressurePreconditioner::Jacobi, 1e-7f);
        firstProjection(multigrid, 48, PressurePreconditioner::Multigrid, 1e-7f);
        const float diff = maxVelocityDifference(jacobi.gas, multigrid.gas);
        std::printf("  48^3, tolerance 1e-7: Jacobi %d iterations, multigrid %d; the flows differ by at most %.1e m/s (inflow 1 m/s), "
                    "max |div| %.1e / %.1e 1/s\n",
                    jacobi.iterations, multigrid.iterations, diff, jacobi.gas.maxDivergence(), multigrid.gas.maxDivergence());
        CHECK(diff < 1e-4f, "multigrid and Jacobi PCG give different flows: %e m/s", diff);
    }
    const int sizes[3] = {32, 64, 128};
    int itJ[3] = {}, itM[3] = {};
    double msJ[3] = {}, msM[3] = {};
    for (int s = 0; s < 3; ++s) {
        PressureRun jacobi, multigrid;
        firstProjection(jacobi, sizes[s], PressurePreconditioner::Jacobi, 1e-4f);
        firstProjection(multigrid, sizes[s], PressurePreconditioner::Multigrid, 1e-4f);
        itJ[s] = jacobi.iterations;
        itM[s] = multigrid.iterations;
        msJ[s] = jacobi.ms;
        msM[s] = multigrid.ms;
        std::printf("  %3d^3 cold start to 1e-4, the projection: Jacobi PCG %4d iterations %7.1f ms | multigrid PCG %3d iterations %6.1f ms "
                    "(%.1fx faster)\n",
                    sizes[s], itJ[s], msJ[s], itM[s], msM[s], msJ[s] / std::max(msM[s], 1e-3));
    }
    const double fallbacks = Probe::snapshot().value("gas/multigrid fallbacks");
    CHECK(fallbacks == 0, "the multigrid PCG broke down %.0f times and handed over to Jacobi", fallbacks);
    CHECK(itM[2] <= 2 * itM[0] + 2, "multigrid iterations grow with the grid: %d on 32^3, %d on 128^3", itM[0], itM[2]);
    CHECK(itM[2] * 5 < itJ[2], "multigrid should need far fewer iterations on 128^3: %d vs %d", itM[2], itJ[2]);
}

namespace {

// The pressure's cost in a scene: frames timed after a warm-up, mean milliseconds of the pressure
// solve and mean iterations per frame (a frame may hold several gas steps).
struct SceneCost { double ms = 0, iterations = 0; };

SceneCost timeScene(const std::function<void(PressurePreconditioner)>& load, const std::function<void()>& frame,
                    const std::function<double()>& pressureMs, const std::function<int()>& iterations,
                    PressurePreconditioner pre, int warmup, int frames) {
    load(pre);
    for (int f = 0; f < warmup; ++f) frame();
    SceneCost c;
    for (int f = 0; f < frames; ++f) {
        frame();
        c.ms += pressureMs();
        c.iterations += iterations();
    }
    c.ms /= frames;
    c.iterations /= frames;
    return c;
}

double median3(double a, double b, double c) { return std::max(std::min(a, b), std::min(std::max(a, b), c)); }

} // namespace

// Speed of the two preconditioners on the scenes that use the gas: Jacobi then multigrid, three
// rounds back to back (other programs share the processor), the median of each. Run on request:
// RF_TEST=benchmark: gas pressure
void testPressureBenchmark() {
    struct Row { const char* name; int warmup, frames; std::function<void(PressurePreconditioner)> load; std::function<void()> frame; };
    Simulation sim;
    GasSolver box;
    auto simLoad = [&](Preset p) {
        return [&sim, p](PressurePreconditioner pre) { loadSample(sim, p); sim.grid.params.pressurePreconditioner = pre; };
    };
    auto simFrame = [&] { sim.stepFrame(); };
    auto boxLoad = [&](PressurePreconditioner pre) {
        box = GasSolver();
        box.params.domainSize = {1, 1, 1};
        box.params.resolutionX = 64;
        for (BoundaryType& b : box.params.bc) b = BoundaryType::Wall;
        box.params.smokeRake = false;
        box.params.heatBuoyancy = 4.0f;
        box.source.enabled = true;
        box.source.center = {0.5f, 0.2f, 0.5f};
        box.source.radius = 0.1f;
        box.params.pressurePreconditioner = pre;
        box.reset({0, 0, 0}, nullptr);
    };
    auto boxFrame = [&] { box.step(1.0f / 60); };
    std::vector<Row> rows = {
        {"smoke box 64^3, closed", 30, 60, boxLoad, boxFrame},
        {"fire (burner and curtain)", 30, 60, simLoad(Preset::Fire), simFrame},
        {"wind tunnel: sphere", 30, 60, simLoad(Preset::TunnelSphere), simFrame},
        {"wind tunnel: cylinder", 30, 60, simLoad(Preset::TunnelCylinder), simFrame},
        {"plasma: magnetosphere (Boris weights)", 30, 60, simLoad(Preset::Magnetosphere), simFrame},
    };
    std::printf("  %-40s %12s %12s %10s %10s %8s\n", "scene (grid)", "Jacobi ms", "multigrid ms", "Jacobi it", "MG it", "speedup");
    for (Row& row : rows) {
        const bool isBox = row.name[0] == 's';
        auto pressureMs = [&] { return isBox ? double(box.lastPressureMs()) : Probe::snapshot().value("gas/pressure ms"); };
        auto iterations = [&] { return isBox ? box.lastPressureIterations() : sim.grid.lastPressureIterations(); };
        SceneCost j[3], m[3];
        for (int round = 0; round < 3; ++round) {
            j[round] = timeScene(row.load, row.frame, pressureMs, iterations, PressurePreconditioner::Jacobi, row.warmup, row.frames);
            m[round] = timeScene(row.load, row.frame, pressureMs, iterations, PressurePreconditioner::Multigrid, row.warmup, row.frames);
        }
        const double jm = median3(j[0].ms, j[1].ms, j[2].ms), mm = median3(m[0].ms, m[1].ms, m[2].ms);
        const GasSolver& g = isBox ? box : sim.grid;
        char name[80];
        std::snprintf(name, sizeof name, "%s %dx%dx%d", row.name, g.nx(), g.ny(), g.nz());
        std::printf("  %-40s %12.2f %12.2f %10.0f %10.0f %7.1fx\n", name, jm, mm, j[1].iterations, m[1].iterations, jm / std::max(mm, 1e-6));
    }
}

namespace {

// Kinetic energy 1/2 rho |u|^2 of the gas, summed over the cells from the cell-centred velocity [J].
double gasKineticEnergy(const GasSolver& g) {
    const double volume = double(g.dx()) * g.dx() * g.dx();
    double e = 0;
    for (int k = 0; k < g.nz(); ++k)
        for (int j = 0; j < g.ny(); ++j)
            for (int i = 0; i < g.nx(); ++i) e += 0.5 * g.params.fluidDensity * double(length2(g.cellVelocity(i, j, k))) * volume;
    return e;
}

// The Taylor-Green vortex u = sin(pi x) cos(pi y), v = -cos(pi x) sin(pi y) in the unit box with
// free-slip walls, WITHOUT viscosity: an exact steady flow of the Euler equations, so every joule
// it loses is lost by the numerics.
void startInviscidTaylorGreen(GasSolver& g, bool reflection) {
    const float dx = 1.0f / 32;
    g.params.domainSize = {1, 1, 4 * dx};
    g.params.resolutionX = 32;
    g.params.inflowSpeed = 0;
    g.params.smokeRake = false;
    g.params.kinematicViscosity = 0;
    g.params.wallFriction = false;
    g.params.advectionReflection = reflection;
    for (auto& bc : g.params.bc) bc = BoundaryType::Wall;
    g.reset({0, 0, 0}, nullptr);
    g.setVelocity([](const Vector3& x) {
        return Vector3(std::sin(kPi * x.x) * std::cos(kPi * x.y), -std::cos(kPi * x.x) * std::sin(kPi * x.y), 0);
    });
}

} // namespace

// Advection-reflection against the classic projection (AdvectionReflection.cpp). The inviscid
// Taylor-Green vortex runs 2 s in both solvers side by side, the same step dt = dx / 2: the
// classic step loses energy to its projection every step, the reflected one keeps nearly all of it
// (Zehnder, Narain & Thomaszewski 2018, fig. 3: "two orders of magnitude" less loss). Neither may
// make energy. Printed: the energy kept by each and the milliseconds per step (the price of the
// second pressure solve); with RF_PLOT_DIR, the energy over time for the chapter's graph.
void testAdvectionReflection() {
    GasSolver classic, mirrored;
    startInviscidTaylorGreen(classic, false);
    startInviscidTaylorGreen(mirrored, true);
    const double e0 = gasKineticEnergy(classic), dt = 0.5 / 32;
    double classicMs = 0, mirroredMs = 0, largest = 0;
    const char* dir = std::getenv("RF_PLOT_DIR");
    FILE* csv = dir && *dir ? std::fopen((std::string(dir) + "/advection_reflection.csv").c_str(), "w") : nullptr;
    if (csv) std::fprintf(csv, "t_s,projection,reflection\n0,1,1\n");
    const int steps = int(2.0 / dt + 0.5);
    for (int s = 1; s <= steps; ++s) {
        auto t0 = std::chrono::steady_clock::now();
        classic.step(float(dt));
        auto t1 = std::chrono::steady_clock::now();
        mirrored.step(float(dt));
        auto t2 = std::chrono::steady_clock::now();
        classicMs += std::chrono::duration<double, std::milli>(t1 - t0).count();
        mirroredMs += std::chrono::duration<double, std::milli>(t2 - t1).count();
        const double kept = gasKineticEnergy(mirrored) / e0;
        largest = std::max(largest, kept);
        if (csv) std::fprintf(csv, "%.6g,%.9g,%.9g\n", s * dt, gasKineticEnergy(classic) / e0, kept);
    }
    if (csv) std::fclose(csv);
    const double keptClassic = gasKineticEnergy(classic) / e0, keptMirrored = gasKineticEnergy(mirrored) / e0;
    std::printf("  inviscid Taylor-Green, 32 cells, 2 s: energy kept %.4f with the projection, %.4f with the reflection "
                "(lost %.2f%% vs %.2f%%); %.2f vs %.2f ms per step (x%.2f)\n",
                keptClassic, keptMirrored, 100 * (1 - keptClassic), 100 * (1 - keptMirrored), classicMs / steps, mirroredMs / steps,
                mirroredMs / std::max(classicMs, 1e-9));
    CHECK(1 - keptMirrored < 0.25 * (1 - keptClassic), "the reflection must lose far less energy: %f vs %f", 1 - keptMirrored, 1 - keptClassic);
    CHECK(largest < 1.0 + 1e-3, "the reflected flow made energy: %f of the start", largest);

    // A steady inflow stays what it is set to: an empty wind tunnel at 2 m/s, reflection on, 1 s.
    // The reflection and the half-step advections must not change a uniform flow in from the inlet
    // and out of the outlet - neither at the inlet nor anywhere downstream.
    GasSolver tunnel;
    tunnel.params.domainSize = {2, 1, 1};
    tunnel.params.resolutionX = 32;
    tunnel.params.inflowSpeed = 2.0f;
    tunnel.params.smokeRake = false;
    tunnel.params.advectionReflection = true;
    tunnel.reset({0, 0, 0}, nullptr);
    for (double t = 0; t < 1.0;) t += tunnel.step(0.05f);
    double inlet = 0, everywhere = 0;
    int inletCells = 0, cells = 0;
    for (int k = 0; k < tunnel.nz(); ++k)
        for (int j = 0; j < tunnel.ny(); ++j) {
            inlet += tunnel.cellVelocity(0, j, k).x, ++inletCells;
            for (int i = 0; i < tunnel.nx(); ++i) everywhere += tunnel.cellVelocity(i, j, k).x, ++cells;
        }
    inlet /= inletCells;
    everywhere /= cells;
    std::printf("  empty tunnel at 2 m/s with the reflection: mean u %.5f m/s at the inlet, %.5f m/s everywhere\n", inlet, everywhere);
    CHECK(std::fabs(inlet / 2.0 - 1) < 1e-3 && std::fabs(everywhere / 2.0 - 1) < 1e-3, "the reflection changed a steady inflow: %f / %f m/s",
          inlet, everywhere);
}
