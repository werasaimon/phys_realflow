// Standard rigid-body stress scenes, shared with the editor's samples: columns, interlocked
// compound rings and a mixed compound sandbox. Inspired by Jolt StackTest, Box3D WavePileTest
// and Bullet BenchmarkDemo; geometry, durations and acceptance limits here are our own.
// Timings exclude the measurements below. Sleeping-off runs expose motion that sleep can hide.
#include "TestRunner.h"
#include "Tests.h"
#include "core/Format.h"
#include "core/Parallel.h"

#include <filesystem>
#include <fstream>

namespace {

struct SceneReport {
    bool finite = true;
    int dynamic = 0, asleep = 0;
    float sleepAt = -1, floorDepth = 0, overlap = 0, restPath = 0, restSpeed = 0, restSpin = 0, linkDistance = 0;
    std::vector<double> milliseconds, activeMilliseconds;
};

bool finiteVector(const Vector3& v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

// Check every component separately: summing components could conceal a non-finite state.
bool finiteBody(const RigidBody& b) {
    return finiteVector(b.pos) && finiteVector(b.vel) && finiteVector(b.angVel) &&
           std::isfinite(b.rot.x) && std::isfinite(b.rot.y) && std::isfinite(b.rot.z) && std::isfinite(b.rot.w);
}

std::ofstream frameCsv(const char* id) {
    if (const char* dir = std::getenv("RF_PLOT_DIR"); dir && *dir) {
        std::filesystem::create_directories(dir);
        std::ofstream out(std::filesystem::path(dir) / (std::string(id) + ".csv"));
        CHECK(bool(out), "cannot open frame CSV for %s", id);
        out << "frame,time_s,step_ms,asleep,dynamic,max_floor_penetration_m,late_max_path_m,late_max_speed_m_s,link_distance_m\n";
        return out;
    }
    return {};
}

// The floor is a plane: support(-y) is exact even for a compound. Body-body penetration is
// measured part by part below; calling GJK on a whole compound would close the ring/cup holes.
void observeBodies(const Simulation& sim, SceneReport& r, std::vector<Vector3>& last, std::vector<float>& paths, bool late) {
    r.dynamic = r.asleep = 0;
    const auto& bodies = sim.rigid.bodies();
    for (size_t i = 0; i < bodies.size(); ++i) {
        const RigidBody& b = bodies[i];
        r.finite &= finiteBody(b);
        if (!b.alive || b.mass <= 0) continue;
        ++r.dynamic;
        r.asleep += b.sleeping;
        r.floorDepth = std::max(r.floorDepth, -b.posed().support(Vector3(0, -1, 0)).y);
        if (late) {
            paths[i] += length(b.pos - last[i]);
            r.restPath = std::max(r.restPath, paths[i]);
            r.restSpeed = std::max(r.restSpeed, length(b.vel));
            r.restSpin = std::max(r.restSpin, length(b.angVel));
        }
        last[i] = b.pos;
    }
}

// The first 28 bodies are two chains of 14 rings, with fixed roots and no joints. A separated
// neighbouring pair is a real loss of interlocking, not a constraint's positional error.
float chainDistance(const RigidWorld& w) {
    float worst = 0;
    for (int root : {0, 14})
        for (int k = 1; k < 14; ++k)
            worst = std::max(worst, length(w.bodies()[size_t(root + k)].pos - w.bodies()[size_t(root + k - 1)].pos));
    return worst;
}

SceneReport measureScene(Simulation& sim, const char* id, bool sleep, bool chains = false) {
    constexpr int frames = 600, lateStart = 480; // 10 s, with a final 2 s observation window
    sim.rigid.params.sleeping = sleep;
    SceneReport r;
    std::vector<Vector3> last;
    for (const RigidBody& b : sim.rigid.bodies()) last.push_back(b.pos);
    std::vector<float> paths(last.size(), 0);
    std::ofstream csv = frameCsv(id);
    for (int f = 0; f < frames; ++f) {
        const auto start = std::chrono::steady_clock::now();
        sim.stepFrame();
        r.milliseconds.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count());
        observeBodies(sim, r, last, paths, f >= lateStart);
        if (!r.finite) break;
        if (r.asleep < r.dynamic) r.activeMilliseconds.push_back(r.milliseconds.back());
        if (r.dynamic > 0 && r.asleep == r.dynamic && r.sleepAt < 0) r.sleepAt = float(f + 1) / 60;
        if (chains) r.linkDistance = std::max(r.linkDistance, chainDistance(sim.rigid));
        if (f >= lateStart && f % 60 == 0) r.overlap = std::max(r.overlap, maxPartOverlap(sim.rigid));
        if (csv) csv << f + 1 << ',' << numberText(double(f + 1) / 60) << ',' << numberText(r.milliseconds.back()) << ','
                     << r.asleep << ',' << r.dynamic << ',' << numberText(r.floorDepth) << ',' << numberText(r.restPath)
                     << ',' << numberText(r.restSpeed) << ',' << numberText(r.linkDistance) << '\n';
    }
    if (r.finite) r.overlap = std::max(r.overlap, maxPartOverlap(sim.rigid));
    return r;
}

double percentile(std::vector<double> values, double fraction) {
    if (values.empty()) return 0;
    std::sort(values.begin(), values.end());
    // Nearest-rank percentile, rank counted from one.
    const size_t rank = size_t(std::ceil(fraction * double(values.size())));
    return values[std::min(values.size() - 1, std::max(size_t(1), rank) - 1)];
}

void printReport(const char* id, const SceneReport& r) {
    std::printf("  %s: median %.3f ms, active median %.3f ms, p95 %.3f ms; asleep %d/%d (first all %.3f s); floor %.3f mm, late part overlap %.3f mm; "
                "late path %.3f mm, speed %.5f m/s, spin %.5f rad/s; chain separation %.3f mm\n",
                id, percentile(r.milliseconds, 0.5), percentile(r.activeMilliseconds, 0.5), percentile(r.milliseconds, 0.95), r.asleep, r.dynamic, r.sleepAt,
                1000 * r.floorDepth, 1000 * r.overlap, 1000 * r.restPath, r.restSpeed, r.restSpin, 1000 * r.linkDistance);
    CHECK(r.finite && r.milliseconds.size() == 600, "%s did not complete 600 finite frames", id);
    CHECK(r.floorDepth < 0.01f, "%s penetrated the floor by %.3f mm", id, 1000 * r.floorDepth);
}

} // namespace

void testStandardSampleColumns() {
    for (const auto preset : {Preset::RigidTower, Preset::RigidTower200}) {
        Simulation sim;
        loadSample(sim, preset);
        const int count = preset == Preset::RigidTower ? 100 : 200;
        const char* id = count == 100 ? "column-100" : "column-200";
        CHECK(int(sim.rigid.bodies().size()) == count, "%s has the wrong number of bodies", id);
        if (int(sim.rigid.bodies().size()) != count) continue;
        const SceneReport r = measureScene(sim, id, true);
        printReport(id, r);
        const Vector3 top = sim.rigid.bodies().back().pos;
        const float expected = 0.1f + 0.2f * float(count - 1);
        const float drift = std::sqrt(top.x * top.x + top.z * top.z);
        std::printf("  %s: height error %.3f mm, top drift %.3f mm\n", id, 1000 * (top.y - expected), 1000 * drift);
        CHECK(std::fabs(top.y - expected) < 0.01f, "%s lost its height", id);
        CHECK(drift < 0.01f, "%s drift %.4f m (target 1 cm)", id, drift);
        CHECK(r.asleep == count, "%s has %d awake bodies", id, count - r.asleep);
        CHECK(r.overlap < 0.01f, "%s late overlap %.3f mm", id, 1000 * r.overlap);
    }
}

void testStandardSampleChains() {
    Simulation sim;
    loadSample(sim, Preset::RigidChains);
    CHECK(sim.rigid.bodies().size() == 35 && sim.rigid.joints().empty(), "the chain sample must have 35 bodies and no joints");
    if (sim.rigid.bodies().size() != 35) return;
    const SceneReport r = measureScene(sim, "compound-chains-14", false, true);
    printReport("compound-chains-14", r);
    // Ring geometry: R=50 mm, tube=12 mm. Taut centres are 2(R-r)=76 mm apart; 7 mm
    // accommodates decomposition plus contact slop, as in testChainOfRings (8 rings).
    CHECK(r.linkDistance < 0.083f, "chain lost interlocking: neighbours %.3f mm apart", 1000 * r.linkDistance);
    CHECK(r.overlap < 0.02f, "chain late part penetration %.3f mm", 1000 * r.overlap);
}

void testStandardSandbox() {
    for (bool sleep : {true, false}) {
        Simulation sim;
        loadSample(sim, Preset::RigidSandbox);
        const char* id = sleep ? "compound-sandbox-sleep" : "compound-sandbox-awake";
        CHECK(sim.rigid.bodies().size() == 24, "%s must contain 24 bodies", id);
        const SceneReport r = measureScene(sim, id, sleep);
        printReport(id, r);
        if (sleep && r.asleep != 24)
            for (size_t i = 0; i < sim.rigid.bodies().size(); ++i) {
                const RigidBody& b = sim.rigid.bodies()[i];
                if (!b.sleeping) std::printf("  awake body %zu (%s): speed %.6f m/s, spin %.6f rad/s, sleep timer %.4f s\n",
                    i, i % 4 == 0 ? "table" : i % 4 == 1 ? "teapot" : i % 4 == 2 ? "ring" : "cup",
                    length(b.vel), length(b.angVel), b.sleepTimer);
            }
        CHECK(r.overlap < 0.02f, "%s late part penetration %.3f mm", id, 1000 * r.overlap);
        if (sleep) CHECK(r.asleep == 24, "%s: %d bodies did not sleep within 10 s", id, 24 - r.asleep);
        else CHECK(r.restPath < 0.01f, "%s: residual path %.3f mm in the last 2 s", id, 1000 * r.restPath);
    }
}

void testStandardRigidDeterminism() {
    ThreadPool& pool = ThreadPool::instance();
    for (const auto preset : {Preset::RigidTower200, Preset::RigidChains, Preset::RigidSandbox}) {
        uint64_t hashes[3];
        const int threads[] = {1, 2, 0}; // 0 = every worker in this process's pool
        for (int t = 0; t < 3; ++t) {
            pool.setActiveThreads(threads[t]);
            Simulation sim;
            loadSample(sim, preset);
            for (int frame = 0; frame < 180; ++frame) sim.stepFrame();
            StateHash hash;
            for (const RigidBody& b : sim.rigid.bodies()) {
                hash.add(b.pos); hash.add(b.rot); hash.add(b.vel); hash.add(b.angVel);
                hash.add(b.sleepTimer); hash.add(float(b.sleeping));
            }
            hashes[t] = hash.h;
        }
        std::printf("  preset %d, 180 frames, 1/2/%d threads: %016llx / %016llx / %016llx\n", int(preset), pool.threadCount(),
                    (unsigned long long)hashes[0], (unsigned long long)hashes[1], (unsigned long long)hashes[2]);
        CHECK(hashes[0] == hashes[1] && hashes[1] == hashes[2], "preset %d changes across worker counts", int(preset));
    }
    pool.setActiveThreads(0);
}
