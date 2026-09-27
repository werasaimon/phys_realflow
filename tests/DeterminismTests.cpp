// Determinism across machines: the same scene gives the same bits on 1, 2 and all threads of the
// pool - the rule of Box2D v3 (E. Catto, "Determinism", box2d.org, 2024). A result that changes
// with the number of cores cannot be reproduced by a reader on another computer. Each scene is run
// from the same start once per thread count (ThreadPool::setActiveThreads) and its whole state is
// hashed: FNV-1a, 64 bits, over the bits of every position, velocity, orientation and field value.
// How the pool keeps the rule is written at the head of src/core/Parallel.h.
#include "TestRunner.h"
#include "Tests.h"

#include "core/Parallel.h"

namespace {

// FNV-1a (Fowler, Noll, Vo): each byte is xor-ed in, then the hash is multiplied by the FNV prime.
// Components are fed one by one - never a struct's raw bytes, whose padding is not part of the state.
struct StateHash {
    uint64_t h = 14695981039346656037ull; // the FNV-1a offset basis
    void add(float v) {
        unsigned char b[4];
        std::memcpy(b, &v, 4);
        for (unsigned char c : b) {
            h ^= c;
            h *= 1099511628211ull; // the FNV prime for 64 bits
        }
    }
    void add(const Vector3& v) { add(v.x); add(v.y); add(v.z); }
    void add(const Quaternion& q) { add(q.w); add(q.x); add(q.y); add(q.z); }
    void add(const std::vector<float>& d) { for (float v : d) add(v); }
};

// Everything that moves in a simulation: bodies, particles, the gas and its magnetic field.
uint64_t hashSimulation(const Simulation& sim) {
    StateHash s;
    for (const RigidBody& b : sim.rigid.bodies()) { s.add(b.pos); s.add(b.rot); s.add(b.vel); s.add(b.angVel); }
    for (const Vector3& x : sim.particles.positions()) s.add(x);
    for (const Vector3& v : sim.particles.velocities()) s.add(v);
    const GasSolver& g = sim.grid;
    for (int k = 0; k < g.nz(); ++k)
        for (int j = 0; j < g.ny(); ++j)
            for (int i = 0; i < g.nx(); ++i) s.add(g.cellVelocity(i, j, k));
    s.add(g.smoke().d);
    s.add(g.temperature().d);
    s.add(g.fuel().d);
    if (g.magnetic.enabled) { s.add(g.magnetic.bx.d); s.add(g.magnetic.by.d); s.add(g.magnetic.bz.d); }
    return s.h;
}

// A sample scene stepped `frames` frames on `threads` threads (0: all), hashed.
uint64_t runSample(Preset p, int frames, int resolution, int threads) {
    ThreadPool::instance().setActiveThreads(threads);
    Simulation sim;
    loadSample(sim, p);
    if (resolution > 0) { sim.grid.params.resolutionX = resolution; sim.reset(); }
    for (int f = 0; f < frames; ++f) sim.stepFrame();
    ThreadPool::instance().setActiveThreads(0);
    return hashSimulation(sim);
}

// A torsional Alfven wave on its own (the set-up of PlasmaTests, smaller): a swirl in a uniform
// field, the Lorentz force and the induction in turn, 60 steps, hashed.
uint64_t runAlfvenWave(int threads) {
    ThreadPool::instance().setActiveThreads(threads);
    const int nx = 24, ny = 12, nz = 12;
    const float dx = 1.0f / nx, a = 0.15f, yc = 0.25f, zc = 0.25f;
    MagneticField m;
    m.applied = {0.01f, 0, 0};
    m.conductivity = 1e12f;
    m.numericalDissipation = 0;
    m.reset(nx, ny, nz, dx, Vector3(0.0f));
    Field3 u, v, w;
    u.init(nx + 1, ny, nz, {0, 0.5f, 0.5f});
    v.init(nx, ny + 1, nz, {0.5f, 0, 0.5f});
    w.init(nx, ny, nz + 1, {0.5f, 0.5f, 0});
    auto swirl = [&](float y, float z) { const float r2 = (sqr(y - yc) + sqr(z - zc)) / (a * a); return r2 < 1 ? 0.3f * sqr(1 - r2) : 0.0f; };
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i < nx; ++i) v.at(i, j, k) = -swirl(j * dx, (k + 0.5f) * dx) * ((k + 0.5f) * dx - zc) * std::sin(kPi * (i + 0.5f) * dx);
    for (int k = 0; k <= nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) w.at(i, j, k) = swirl((j + 0.5f) * dx, k * dx) * ((j + 0.5f) * dx - yc) * std::sin(kPi * (i + 0.5f) * dx);
    const std::vector<uint8_t> none;
    for (int step = 0; step < 60; ++step) {
        m.applyLorentzForce(u, v, w, none, 1.0f, 2e-3f);
        m.induce(u, v, w, 1.0f, 2e-3f);
    }
    ThreadPool::instance().setActiveThreads(0);
    StateHash s;
    s.add(m.bx.d); s.add(m.by.d); s.add(m.bz.d);
    s.add(u.d); s.add(v.d); s.add(w.d);
    return s.h;
}

// One line of the report and its check. knownSite: a place in src/rigid or src/particles that is
// known to depend on the thread count and waits for its fix - then a difference is reported, not
// failed.
void report(const char* name, const uint64_t h[3], const char* knownSite) {
    const bool same = h[0] == h[1] && h[1] == h[2];
    std::printf("  %s: 1 thread %016llx, 2 threads %016llx, all %016llx %s\n", name, (unsigned long long)h[0],
                (unsigned long long)h[1], (unsigned long long)h[2], same ? "equal" : "DIFFER");
    if (!same && knownSite) { std::printf("  EXPECTED-NONDETERMINISTIC: %s\n", knownSite); return; }
    CHECK(same, "%s depends on the number of threads", name);
}

} // namespace

void testDeterminismAcrossThreads() {
    // A few frames are enough: a sum rounded differently changes the bits in the first frame that
    // makes it, no chaos is needed to show it. Each scene runs its parallel passes every frame. The
    // frame counts are long enough to catch the old fault: with reductions cut by the thread count,
    // the smoke box differed at 8 frames and the fire from 30 frames on (its gas is nearly still
    // before), while the plasma box is too small to be cut differently at all.
    struct Case { const char* name; Preset preset; int frames, resolution; const char* knownSite; };
    const Case cases[] = {
        {"gas: smoke in a closed box", Preset::SmokeSphere, 8, 32, nullptr},
        {"fire: burner and cotton curtain", Preset::Fire, 30, 0, nullptr},
        {"plasma: solar wind against a magnet", Preset::Magnetosphere, 3, 0, nullptr},
        {"rigid: pyramid and projectile", Preset::RigidPyramid, 90, 0, nullptr},
        {"rigid: a hundred teapots (parallel narrow phase)", Preset::RigidTeapots, 30, 0, nullptr},
        {"rigid: 150 bodies on a terrain mesh (BVH)", Preset::Terrain, 60, 0, nullptr},
        {"particles: dam break", Preset::DamBreak, 10, 0, nullptr},
        {"particles: cloth and soft bodies", Preset::SoftCloth, 5, 0, nullptr},
    };
    std::printf("  threads in the pool: %d\n", ThreadPool::instance().threadCount());
    for (const Case& c : cases) {
        const auto t0 = std::chrono::steady_clock::now();
        uint64_t h[3];
        const int threads[3] = {1, 2, 0};
        for (int t = 0; t < 3; ++t) h[t] = runSample(c.preset, c.frames, c.resolution, threads[t]);
        report(c.name, h, c.knownSite);
        std::printf("    (%d frames, %.1f s for the three runs)\n", c.frames,
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }
    const uint64_t alfven[3] = {runAlfvenWave(1), runAlfvenWave(2), runAlfvenWave(0)};
    report("MHD: torsional Alfven wave, 60 steps", alfven, nullptr);
}

namespace {

// The rigid benchmark scene of the golden hash: a pyramid of six boxes, a stack of five, a ball
// rolling across the floor and a chain of four links hinged to the world, released horizontal.
// Only exactly representable numbers and no trigonometry in the set-up, so the start is the same
// bits on every compiler.
void buildRigidBenchmark(RigidWorld& world) {
    world.setDomain(AABB({-3, 0, -3}, {3, 4, 3}));
    const Vector3 half(0.1f), grey(0.5f);
    for (int row = 0; row < 3; ++row)
        for (int k = 0; k < 3 - row; ++k)
            world.addBox({-1.5f + 0.1f * float(row) + 0.2f * float(k), 0.1f + 0.2f * float(row), 0.0f}, half, Quaternion(), 500.0f, grey);
    for (int k = 0; k < 5; ++k) world.addBox({0.0f, 0.1f + 0.2f * float(k), -1.0f}, half, Quaternion(), 500.0f, grey);
    const int ball = world.addSphere({1.0f, 0.1f, 0.0f}, 0.1f, 1000.0f, grey);
    world.bodies()[size_t(ball)].vel = {-1.0f, 0.0f, 0.5f};
    world.bodies()[size_t(ball)].angVel = {5.0f, 0.0f, 10.0f}; // rolling: w = n x v / r
    int previous = -1;
    for (int k = 0; k < 4; ++k) {
        const int link = world.addBox({1.6f + 0.2f * float(k), 3.0f, 1.0f}, {0.1f, 0.025f, 0.025f}, Quaternion(), 500.0f, grey);
        world.addHingeJoint(link, previous, {1.5f + 0.2f * float(k), 3.0f, 1.0f}, {0, 0, 1});
        previous = link;
    }
}

} // namespace

// The golden hash: the benchmark scene run for 10 s (600 frames of 1/60 s) and hashed must give
// one known number. It is kept for the RF_STRICT_FP=ON build (no fused multiply-add, no
// -march=native) with MinGW's C runtime: std::sin and friends round their last bit differently in
// glibc or MSVC, so another compiler needs step 2 of src/math/ElementaryFunctions.h first. The
// default build (-march=native: FMA wherever the compiler likes) gives another number, printed only.
// A change of the constant is a change of the rigid solver's results: explain it in the commit.
void testRigidGoldenHash() {
    constexpr uint64_t kGolden = 0x912fa823d3669448ull; // RF_STRICT_FP=ON, MinGW (GCC 11.2), Release
    RigidWorld world;
    buildRigidBenchmark(world);
    for (int f = 0; f < 600; ++f) world.step(1.0f / 60.0f);
    StateHash s;
    for (const RigidBody& b : world.bodies()) { s.add(b.pos); s.add(b.rot); s.add(b.vel); s.add(b.angVel); }
    const RigidBody& lastLink = world.bodies().back();
    std::printf("  rigid benchmark scene after 10 s: hash %016llx; %zu bodies, %zu asleep, chain end at (%.4f, %.4f, %.4f)\n",
                (unsigned long long)s.h, world.bodies().size(), world.sleepingCount(), lastLink.pos.x, lastLink.pos.y, lastLink.pos.z);
#ifdef RF_STRICT_FP
    std::printf("  strict floating point: the hash must be %016llx\n", (unsigned long long)kGolden);
    CHECK(s.h == kGolden, "the golden hash changed: %016llx (was %016llx) - the rigid results changed", (unsigned long long)s.h,
          (unsigned long long)kGolden);
#else
    std::printf("  default build (-march=native, FMA): printed only; the constant is checked in the RF_STRICT_FP build\n");
#endif
}

// src/math/ElementaryFunctions.h, step 1 of 2: rf::sin and the others still call std::, so they
// must return exactly its bits. Step 2 (our own series) will turn this into an accuracy test.
void testElementaryFunctionsAreStd() {
    auto sameBits = [](float a, float b) { return std::memcmp(&a, &b, sizeof a) == 0; };
    int differ = 0, samples = 0;
    for (int k = -20000; k <= 20000; ++k) {
        const float x = 0.00123f * float(k), u = float(k) / 20000.0f; // x in [-24.6, 24.6], u in [-1, 1]
        differ += !sameBits(rf::sin(x), std::sin(x)) + !sameBits(rf::cos(x), std::cos(x)) + !sameBits(rf::atan(x), std::atan(x));
        differ += !sameBits(rf::atan2(x, 1.5f - u), std::atan2(x, 1.5f - u)) + !sameBits(rf::exp(u), std::exp(u));
        differ += !sameBits(rf::acos(u), std::acos(u)) + !sameBits(rf::asin(u), std::asin(u)) + !sameBits(rf::pow(1.5f + u, x), std::pow(1.5f + u, x));
        samples += 8;
    }
    std::printf("  rf::sin, cos, atan, atan2, exp, acos, asin, pow: %d of %d results differ from std:: in any bit\n", differ, samples);
    CHECK(differ == 0, "%d results of rf:: differ from std::", differ);
}

// The whole fire test on 1, 2 and all threads (minutes; on request): the number of burnt threads
// was 65 / 34 / 20 before the reductions were cut by blocks independent of the thread count.
void testFireAcrossThreads() {
    const int threads[3] = {1, 2, 0};
    int burnt[3];
    for (int t = 0; t < 3; ++t) {
        ThreadPool::instance().setActiveThreads(threads[t]);
        Simulation sim;
        loadSample(sim, Preset::Fire);
        float Tmax = 0;
        for (int f = 1; f <= 300; ++f) {
            sim.stepFrame();
            for (float T : sim.grid.temperature().d) Tmax = std::max(Tmax, T);
        }
        ThreadPool::instance().setActiveThreads(0);
        burnt[t] = sim.particles.cloths()[0].burntThreads;
        std::printf("  %s: %d threads burnt through, flame Tmax %.0f K\n", t == 0 ? "1 thread" : t == 1 ? "2 threads" : "all threads",
                    burnt[t], Tmax + sim.grid.combustion.ambientTemperature);
    }
    CHECK(burnt[0] == burnt[1] && burnt[1] == burnt[2], "the fire burns %d / %d / %d threads on 1 / 2 / all threads", burnt[0],
          burnt[1], burnt[2]);
}
