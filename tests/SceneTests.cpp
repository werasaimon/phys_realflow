// The scenes as a whole: every preset runs without NaN, scene switches carry nothing over, and
// the solvers agree on the shared quantities (gravity, Coulomb friction).
#include "TestRunner.h"
#include "Tests.h"

void testSimulationPresets() {
    // The registry of samples/ and enum Preset agree: one scene per value, in order, named.
    const std::vector<SampleEntry>& all = samples();
    CHECK(int(all.size()) == int(Preset::Count), "samples(): %zu scenes for %d presets", all.size(), int(Preset::Count));
    for (int p = 0; p < int(all.size()); ++p) {
        const SampleEntry& e = all[p];
        CHECK(int(e.id) == p, "sample %d registered as preset %d", p, int(e.id));
        CHECK(e.name && e.name[0] && e.category && e.category[0], "sample %d has no name or category", p);
        CHECK(e.create != nullptr, "sample %d cannot be created", p);
    }
    Simulation sim;
    for (int p = 0; p < int(Preset::Count); ++p) {
        loadSample(sim, Preset(p));
        CHECK(sim.scene() && sim.scene()->name == all[p].name, "preset %d: the loaded scene is not '%s'", p, all[p].name);
        if (sim.mode() == SimMode::WindTunnel) sim.grid.params.resolutionX = 32, sim.reset();
        const int frames = Preset(p) == Preset::Terrain ? 1 : 3; // the terrain's 150 bodies on 50 000 triangles: heavier
        for (int f = 0; f < frames; ++f) sim.stepFrame();
        RenderSnapshot snap;
        sim.fillSnapshot(snap);
        bool finite = true;
        for (const Vector3& q : snap.particles) finite &= std::isfinite(q.x + q.y + q.z);
        for (const auto& b : snap.bodies) finite &= std::isfinite(b.pos.x + b.pos.y + b.pos.z);
        CHECK(finite, "preset %s produced NaN", all[p].name);
        CHECK(snap.sceneName == all[p].name, "snapshot names the scene '%s', expected '%s'", snap.sceneName.c_str(), all[p].name);
    }
}

// Scene switches and resets must not carry state of the old scene over; the solvers must agree on
// the shared quantities (gravity, friction) whatever the order or the mode.
void testCoherence() {
    // 1) The water of the Hydro scene does not stay in the gas grid of the next (gas-only) scene.
    {
        Simulation sim;
        loadSample(sim, Preset::Hydro);
        sim.grid.params.resolutionX = 32;
        sim.reset();
        for (int f = 0; f < 3; ++f) sim.stepFrame();
        loadSample(sim, Preset::TunnelSphere);
        sim.grid.params.resolutionX = 32;
        sim.reset();
        for (int f = 0; f < 2; ++f) sim.stepFrame();
        CHECK(sim.grid.liquidCellCount() == 0, "phantom water in the tunnel: %d cells", sim.grid.liquidCellCount());
        CHECK(!sim.surfaceLoads().triangles.empty(), "tunnel loads missing");
        loadSample(sim, Preset::RigidPyramid);
        CHECK(sim.surfaceLoads().triangles.empty(), "surface loads of the old scene kept");
    }
    // 2) XPBD after a reset to fewer bodies: fresh contacts (the old ones index bodies that are gone).
    {
        Simulation sim;
        loadSample(sim, Preset::RigidPyramid);
        sim.rigid.params.solver = RigidSolver::XPBD;
        sim.rigid.params.substeps = 31; // leaves the collision counter unaligned
        for (int i = 0; i < 60; ++i) sim.rigid.addSphere({-1.5f + 0.05f * i, 0.2f, 1.0f}, 0.1f, 500.0f, Vector3(1));
        for (int f = 0; f < 11; ++f) sim.stepFrame();
        sim.reset();
        sim.stepFrame();
        float vmax = 0;
        for (const RigidBody& b : sim.rigid.bodies()) vmax = std::max(vmax, length(b.vel));
        CHECK(vmax < 20.0f, "XPBD after reset: bodies at %f m/s", vmax);
    }
    // 3) Coulomb friction: a box sliding at 3 m/s stops after v^2 / (2 mu g) (the shock pass adds none).
    {
        RigidWorld w;
        w.setDomain(AABB({-10, 0, -10}, {10, 10, 10}));
        w.params.sleeping = false;
        int b = w.addBox({-5.0f, 0.05f, 0.0f}, Vector3(0.05f), Quaternion(), 500.0f, Vector3(1));
        for (int k = 0; k < 300; ++k) w.step(1.0f / 600);
        const float x0 = w.bodies()[b].pos.x;
        w.bodies()[b].vel = {3.0f, 0, 0};
        for (int k = 0; k < 2400; ++k) w.step(1.0f / 600);
        const float mu = std::sqrt(0.5f * 0.6f), coulomb = 9.0f / (2 * mu * 9.81f), d = w.bodies()[b].pos.x - x0;
        std::printf("  sliding box stops after %.3f m (Coulomb %.3f m)\n", d, coulomb);
        CHECK(std::fabs(d - coulomb) < 0.1f * coulomb, "friction is not Coulomb: %f vs %f m", d, coulomb);
    }
    // 4) Rolling resistance of a small ball on a big static platform: independent of which body was
    //    created first.
    {
        float speed[2];
        for (int order = 0; order < 2; ++order) {
            RigidWorld w;
            w.setDomain(AABB({-10, -1, -10}, {10, 10, 10}));
            w.params.sleeping = false;
            if (order == 0) w.addBox({0, 0.05f, 0}, {4, 0.05f, 4}, Quaternion(), 0.0f, Vector3(1));
            int ball = w.addSphere({-2, 0.15f, 0}, 0.05f, 500, Vector3(1));
            if (order == 1) w.addBox({0, 0.05f, 0}, {4, 0.05f, 4}, Quaternion(), 0.0f, Vector3(1));
            for (int k = 0; k < 60; ++k) w.step(1.0f / 600);
            w.bodies()[ball].vel = {2, 0, 0};
            w.bodies()[ball].angVel = {0, 0, -40};
            for (int k = 0; k < 600; ++k) w.step(1.0f / 600);
            speed[order] = length(w.bodies()[ball].vel);
        }
        CHECK(std::fabs(speed[0] - speed[1]) < 0.05f && speed[0] > 1.0f, "rolling ball: %f vs %f m/s", speed[0], speed[1]);
    }
    // 5) The mouse joint also works with the gravity switched off.
    {
        RigidWorld w;
        w.setDomain(AABB({-5, 0, -5}, {5, 5, 5}));
        w.params.gravity = Vector3(0.0f);
        int b = w.addBox({0, 2, 0}, Vector3(0.1f), Quaternion(), 500, Vector3(1));
        w.grab(b, {0, 2, 0});
        w.setGrabTarget({1, 2, 0});
        for (int k = 0; k < 600; ++k) w.step(1.0f / 600);
        CHECK(w.bodies()[b].pos.x > 0.8f, "grab at g = 0: body at x %f", w.bodies()[b].pos.x);
    }
    // 6) One gravity for the scene; a preset restores the default everywhere.
    {
        Simulation sim;
        loadSample(sim, Preset::Fire);
        sim.setGravity({0.0f, -1.62f, 0.0f});
        CHECK(sim.particles.params.gravity.y == -1.62f && sim.grid.combustion.gravity == 1.62f, "gravity not shared");
        loadSample(sim, Preset::Fire);
        CHECK(sim.gravity().y == -9.81f && sim.particles.params.gravity.y == -9.81f, "preset gravity not restored");
    }
    // 7) Fabric that cannot tear, burnt through along a band: the part below falls off the rod.
    {
        ParticleSystem ps;
        ps.params.clothSpacing = 2.0f;
        ps.reset(AABB({-1, -1, -1}, {1, 1, 1}));
        ClothMaterial cotton;
        cotton.flammable = true;
        cotton.areaDensity = 0.2f;
        cotton.strengthWarp = cotton.strengthWeft = 0.0f;
        ps.addCloth({-0.2f, 0.5f, 0.0f}, {0.4f, 0, 0}, {0, -0.6f, 0}, cotton, 16, Vector3(1));
        auto gas = [](const Vector3& x) {
            const bool hot = std::fabs(x.y - 0.185f) < 0.025f;
            return GasHeat{hot ? 1500.0f : 0.0f, hot ? 150e3f : 0.0f};
        };
        for (int f = 0; f < 180; ++f) {
            for (int s = 0; s < 3; ++s) ps.step(1.0f / 180);
            std::vector<FireOutput> out;
            ps.burnCloths(1.0f / 60, gas, 293.0f, out);
        }
        const Cloth& c = ps.cloths()[0];
        float lowest = 1e9f;
        for (int k = 0; k < c.width * c.height; ++k) lowest = std::min(lowest, ps.positions()[c.firstParticle + k].y);
        CHECK(c.burntThreads > 0 && lowest < -0.5f, "burnt-off cloth still held (lowest %f)", lowest);
    }
    // 8) No liquid inside a static obstacle; a non-numeric NACA code builds the default wing.
    {
        Simulation sim;
        sim.obstacle.shape = ObstacleShape::Sphere;
        sim.obstacle.size = 0.4f;
        sim.obstacle.position = {-0.7f, 0.3f, 0.0f};
        sim.rebuildObstacle();
        int inside = 0;
        for (const Vector3& x : sim.particles.positions()) inside += length(x - Vector3(-0.7f, 0.3f, 0.0f)) < 0.2f;
        CHECK(inside == 0, "%d liquid particles inside the obstacle", inside);
        TriMesh wing = primitives::nacaWing("NACA", 1.0f, 1.0f);
        bool finite = !wing.empty();
        for (const Vector3& p : wing.positions) finite &= std::isfinite(p.x + p.y + p.z);
        CHECK(finite, "NACA code 'NACA' gave no wing");
    }
}

// Two runs of the same scene must agree to the bit: a result nobody can reproduce cannot be
// debugged, and a race in the thread pool would show up here first. Rigid bodies (contacts,
// islands, shock passes), particles (a dam break) and the gas (smoke in a closed box, PCG with
// parallel reductions), each fingerprinted after a few frames.
void testDeterminism() {
    auto mix = [](uint64_t& h, float v) { // FNV-1a over the float bits
        uint32_t bits;
        std::memcpy(&bits, &v, 4);
        h ^= bits;
        h *= 1099511628211ull;
    };
    auto fingerprint = [&](Preset p, int frames, int resolution) {
        Simulation sim;
        loadSample(sim, p);
        if (resolution > 0) { sim.grid.params.resolutionX = resolution; sim.reset(); }
        for (int f = 0; f < frames; ++f) sim.stepFrame();
        uint64_t h = 14695981039346656037ull;
        for (const RigidBody& b : sim.rigid.bodies()) for (float v : {b.pos.x, b.pos.y, b.pos.z, b.vel.x, b.vel.y, b.vel.z, b.angVel.x, b.angVel.y, b.angVel.z}) mix(h, v);
        for (const Vector3& x : sim.particles.positions()) for (float v : {x.x, x.y, x.z}) mix(h, v);
        if (sim.mode() == SimMode::WindTunnel)
            for (int k = 0; k < sim.grid.nz(); ++k)
                for (int j = 0; j < sim.grid.ny(); ++j)
                    for (int i = 0; i < sim.grid.nx(); ++i) { const Vector3 u = sim.grid.cellVelocity(i, j, k); mix(h, u.x); mix(h, u.y); mix(h, u.z); }
        return h;
    };
    struct Case { const char* name; Preset preset; int frames, resolution; };
    const Case cases[] = {{"rigid: pyramid + projectile", Preset::RigidPyramid, 90, 0}, {"particles: dam break", Preset::DamBreak, 30, 0},
                          {"gas: smoke in a closed box", Preset::SmokeSphere, 8, 32}};
    for (const Case& c : cases) {
        const uint64_t a = fingerprint(c.preset, c.frames, c.resolution), b = fingerprint(c.preset, c.frames, c.resolution);
        std::printf("  %s: %d frames, fingerprints %016llx / %016llx %s\n", c.name, c.frames, (unsigned long long)a, (unsigned long long)b,
                    a == b ? "equal" : "DIFFER");
        CHECK(a == b, "%s is not reproducible run to run", c.name);
    }
    // The probe is wired into stepFrame: the last run left the frame's channels behind.
    const std::vector<std::string> names = Probe::channels();
    CHECK(std::find(names.begin(), names.end(), "frame/step ms") != names.end(), "stepFrame reports frame/step ms to the Probe");
    CHECK(Probe::snapshot().value("frame/step ms") > 0, "frame/step ms is %g", Probe::snapshot().value("frame/step ms"));
}
