// Long contact-only torus chains: verify the topology meter against known circle pairs, then
// observe the same loaded/swinging/bridge scene as the editor, with sleep disabled.
#include "TestRunner.h"
#include "Tests.h"
#include "samples/Models.h"
#include "samples/TorusChainsScene.h"
#include "core/Format.h"
#include "core/Parallel.h"

#include <filesystem>
#include <fstream>

namespace {

double totalEnergy(const RigidWorld& world) {
    double energy = world.kineticEnergy();
    for (const RigidBody& body : world.bodies())
        if (body.invMass > 0) energy -= double(body.mass) * dot(world.params.gravity, body.pos);
    return energy;
}

bool finiteTorus(const RigidBody& body) {
    const float values[] = {body.pos.x, body.pos.y, body.pos.z, body.vel.x, body.vel.y, body.vel.z,
        body.angVel.x, body.angVel.y, body.angVel.z, body.rot.x, body.rot.y, body.rot.z, body.rot.w};
    for (float value : values) if (!std::isfinite(value)) return false;
    return true;
}

std::ofstream torusCsv() {
    if (const char* dir = std::getenv("RF_PLOT_DIR"); dir && *dir) {
        std::filesystem::create_directories(dir);
        std::ofstream file(std::filesystem::path(dir) / "torus-chains-40.csv");
        CHECK(bool(file), "cannot write torus chain CSV");
        file << "frame,time_s,step_ms,broken_pairs,max_distance_m,energy_j,energy_gain_j,max_floor_depth_m\n";
        return file;
    }
    return {};
}

double torusPercentile(std::vector<double> samples, double fraction) {
    std::sort(samples.begin(), samples.end());
    return samples[std::min(samples.size() - 1, size_t(std::ceil(fraction * samples.size())) - 1)];
}

} // namespace

void testTorusLinkingMeter() {
    RigidWorld world;
    const auto ring = ringShape();
    const Quaternion crossing = Quaternion::fromAxisAngle({1, 0, 0}, 0.5f * kPi);
    world.addCompound(ring, {0, 0, 0}, Quaternion(), 0, Vector3(1));
    world.addCompound(ring, {0.07f, 0, 0}, crossing, 0, Vector3(1));
    world.addCompound(ring, {0, 0.07f, 0}, crossing, 0, Vector3(1));
    world.addCompound(ring, {0, 0.04f, 0}, crossing, 0, Vector3(1)); // both disk crossings cancel
    world.addCompound(ring, {0.14f, 0, 0}, crossing, 0, Vector3(1)); // both outside the disk
    world.addCompound(ring, {0, 0.04f, 0}, Quaternion(), 0, Vector3(1)); // parallel planes
    auto& b = world.bodies();
    CHECK(std::abs(TorusChainsScene::linkingNumber(b[0], b[1])) == 1, "linked circles must have |Lk|=1");
    CHECK(std::abs(TorusChainsScene::linkingNumber(b[1], b[0])) == 1, "swapping linked circles loses topology");
    for (int k = 2; k < 6; ++k)
        CHECK(TorusChainsScene::linkingNumber(b[0], b[size_t(k)]) == 0, "unlinked fixture %d reports linking", k);
    const Quaternion turn = Quaternion::fromAxisAngle(normalize(Vector3(1, 2, 3)), 1.234f);
    for (auto& body : b) {
        body.pos = turn.rotate(body.pos) + Vector3(1.3f, 2.7f, -0.8f);
        body.rot = (turn * body.rot).normalized();
    }
    CHECK(std::abs(TorusChainsScene::linkingNumber(b[0], b[1])) == 1, "linking changes with a global rigid transform");
    for (int k = 2; k < 6; ++k)
        CHECK(TorusChainsScene::linkingNumber(b[0], b[size_t(k)]) == 0, "rotated unlinked fixture %d reports linking", k);
    Simulation sim;
    loadSample(sim, Preset::TorusChains);
    for (int count = 8; count <= 100; ++count) {
        sim.scene()->setParam(0, float(count));
        sim.reset();
        CHECK(sim.rigid.bodies().size() == size_t(3 * count), "reset ignores link count %d", count);
        for (int root : {0, count, 2 * count})
            for (int k = 1; k < count; ++k)
                CHECK(std::abs(TorusChainsScene::linkingNumber(sim.rigid.bodies()[size_t(root + k - 1)],
                    sim.rigid.bodies()[size_t(root + k)])) == 1, "count %d: unlinked initial pair %d/%d", count, root, k);
        CHECK(maxPartOverlap(sim.rigid) < 0.001f, "count %d: release geometry intersects", count);
    }
}

// Body-frame poses captured after substep 492 of sample 39, bodies 67/68. Cold restart:
// no gravity, velocity, angular velocity, domain walls or joints. A genuine overlap must
// produce contacts and relax to the configured slop. This guards the compound seam filter;
// no topology claim is made for this already crossed pair.
void testTorusContactStarvation() {
    RigidWorld world;
    world.params.gravity = Vector3(0);
    world.params.collideWithDomain = false;
    world.params.sleeping = false;
    const auto ring = ringShape();
    world.addBody(ring, {0.49678934f, 1.689123f, 0.86407435f},
        {0.87695223f, 0.1506506f, 0.07180688f, -0.45066956f}, 7800, Vector3(1));
    world.addBody(ring, {0.5551778f, 1.5995147f, 0.8649265f},
        {0.5809721f, 0.6733655f, -0.35156918f, -0.2923175f}, 7800, Vector3(1));
    const float initial = maxPartOverlap(world);
    CHECK(initial > 0.013f && initial < 0.015f, "replay fixture has wrong penetration %.3f mm", 1000 * initial);
    world.step(1.0f / 600);
    const size_t firstContacts = world.contactCount();
    CHECK(firstContacts > 0, "stationary intersecting tori have no solver contacts");
    for (int step = 1; step < 120; ++step) world.step(1.0f / 600);
    const float final = maxPartOverlap(world);
    std::printf("  cold two-torus replay: initial/final penetration %.3f/%.3f mm, first contacts %zu, slop %.3f mm\n",
        1000 * initial, 1000 * final, firstContacts, 1000 * world.params.slop);
    CHECK(final < world.params.slop + 0.0001f, "unloaded tori remain intersecting %.3f mm", 1000 * final);
    for (const auto& body : world.bodies()) CHECK(finiteTorus(body), "cold torus replay became non-finite");
}

void testTorusChainStress() {
    Simulation sim;
    loadSample(sim, Preset::TorusChains);
    const auto* scene = static_cast<const TorusChainsScene*>(sim.scene());
    const int count = scene->linksPerChain();
    CHECK(count == 40 && sim.rigid.bodies().size() == 120, "stress preset must have three chains of 40 tori");
    CHECK(sim.rigid.joints().empty() && !sim.rigid.params.sleeping, "stress scene conceals contacts with joints or sleep");
    int fixed = 0, initialBroken = 0;
    for (const auto& b : sim.rigid.bodies()) fixed += b.invMass == 0;
    CHECK(fixed == 4, "only four endpoint rings may be fixed, found %d", fixed);
    for (int root : {0, count, 2 * count})
        for (int k = 1; k < count; ++k)
            initialBroken += std::abs(TorusChainsScene::linkingNumber(sim.rigid.bodies()[size_t(root + k - 1)],
                                                                     sim.rigid.bodies()[size_t(root + k)])) != 1;
    CHECK(initialBroken == 0, "%d initial torus pairs are not interlocked", initialBroken);
    const float initialOverlap = maxPartOverlap(sim.rigid);
    CHECK(initialOverlap < 0.001f, "invalid release: initial part penetration %.3f mm", 1000 * initialOverlap);
    const double initialEnergy = totalEnergy(sim.rigid);
    double energyGain = 0;
    float floorDepth = 0, lateOverlap = 0, maxDistance = 0;
    int broken = 0, frames = 0;
    bool finite = true;
    std::vector<double> timings;
    std::ofstream csv = torusCsv();
    for (int frame = 1; frame <= 600; ++frame) {
        const auto start = std::chrono::steady_clock::now();
        sim.stepFrame();
        timings.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count());
        const double energy = totalEnergy(sim.rigid);
        energyGain = std::max(energyGain, energy - initialEnergy);
        broken = std::max(broken, scene->brokenPairs());
        for (const auto& b : sim.rigid.bodies()) {
            finite &= finiteTorus(b);
            if (b.invMass > 0) floorDepth = std::max(floorDepth, -b.posed().support({0, -1, 0}).y);
        }
        for (int root : {0, count, 2 * count})
            for (int k = 1; k < count; ++k)
                maxDistance = std::max(maxDistance, length(sim.rigid.bodies()[size_t(root + k)].pos - sim.rigid.bodies()[size_t(root + k - 1)].pos));
        if (frame >= 480 && frame % 60 == 0) lateOverlap = std::max(lateOverlap, maxPartOverlap(sim.rigid));
        if (csv) csv << frame << ',' << numberText(double(frame) / 60) << ',' << numberText(timings.back()) << ',' << scene->brokenPairs()
                     << ',' << numberText(maxDistance) << ',' << numberText(energy) << ',' << numberText(energyGain) << ',' << numberText(floorDepth) << '\n';
        frames = frame;
        if (!finite) break;
    }
    std::printf("  120 tori / 117 links / 4 fixed / 0 joints: %d frames, first break %d, worst broken %d; initial/late overlap %.3f/%.3f mm; max separation %.3f mm; floor %.3f mm; energy gain %.6f/%.3f J; step median/p95 %.3f/%.3f ms\n",
        frames, scene->firstBrokenFrame(), broken, 1000 * initialOverlap, 1000 * lateOverlap, 1000 * maxDistance, 1000 * floorDepth,
        energyGain, initialEnergy, torusPercentile(timings, 0.5), torusPercentile(timings, 0.95));
    CHECK(finite && frames == 600, "torus scene did not finish 10 finite seconds");
    CHECK(broken == 0 && scene->firstBrokenFrame() < 0, "tori passed through each other; first frame %d, worst %d broken pairs", scene->firstBrokenFrame(), broken);
    CHECK(floorDepth < 0.01f, "torus floor penetration %.3f mm", 1000 * floorDepth);
    CHECK(lateOverlap < 0.0005f, "late torus part penetration %.3f mm exceeds 0.5 mm", 1000 * lateOverlap);
    CHECK(energyGain < 0.01 * initialEnergy, "torus chains gained more than 1%% of initial mechanical energy");
}

void testTorusChainDeterminism() {
    ThreadPool& pool = ThreadPool::instance();
    uint64_t hashes[3];
    const int workers[] = {1, 2, 0};
    for (int run = 0; run < 3; ++run) {
        pool.setActiveThreads(workers[run]);
        Simulation sim;
        loadSample(sim, Preset::TorusChains);
        for (int frame = 0; frame < 120; ++frame) sim.stepFrame();
        StateHash hash;
        for (const auto& body : sim.rigid.bodies()) {
            hash.add(body.pos); hash.add(body.rot); hash.add(body.vel); hash.add(body.angVel);
        }
        hashes[run] = hash.h;
    }
    pool.setActiveThreads(0);
    std::printf("  120 tori, 120 frames, 1/2/%d threads: %016llx / %016llx / %016llx\n", pool.threadCount(),
        (unsigned long long)hashes[0], (unsigned long long)hashes[1], (unsigned long long)hashes[2]);
    CHECK(hashes[0] == hashes[1] && hashes[1] == hashes[2], "torus chain state changes across worker counts");
}
