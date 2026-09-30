// The shock pass can switch from pushing only an upper body to pushing both bodies. Contacts
// sharing a moving support must therefore never write it concurrently, even on different levels.
// Compare full states across workers and repeats; the dense samples also check every frame.
#include "TestRunner.h"
#include "Tests.h"

namespace {

uint64_t rigidStateHash(const RigidWorld& world) {
    StateHash hash;
    for (const auto& body : world.bodies()) {
        hash.add(body.pos); hash.add(body.rot); hash.add(body.vel); hash.add(body.angVel);
        hash.add(body.sleepTimer); hash.add(float(body.sleeping));
    }
    return hash.h;
}

uint64_t movingSupport(int threads) {
    ThreadPool::instance().setActiveThreads(threads);
    RigidWorld world;
    world.setDomain(AABB({-10, 0, -10}, {10, 10, 10}));
    world.params.sleeping = false;
    const int platform = world.addBox({0, 0.1f, 0}, {3, 0.1f, 3}, Quaternion(), 1000, {1, 1, 1});
    world.bodies()[size_t(platform)].vel = {0.5f, 0, 0}; // exceeds the resting-support threshold
    for (int z = 0; z < 8; ++z)
        for (int x = 0; x < 8; ++x)
            world.addBox({(float(x) - 3.5f) * 0.5f, 0.25f, (float(z) - 3.5f) * 0.5f}, Vector3(0.05f), Quaternion(), 1000, {1, 1, 1});
    // 64 upper contacts exceed the shock parallel threshold. The regular contact pass, with
    // fewer than 256 manifolds, is serial: this fixture isolates the shock pass's shared writes.
    for (int step = 0; step < 20; ++step) world.step(1.0f / 600);
    return rigidStateHash(world);
}

std::vector<uint64_t> sampleStates(Preset preset, int frames, int threads) {
    ThreadPool::instance().setActiveThreads(threads);
    Simulation sim;
    loadSample(sim, preset);
    std::vector<uint64_t> states;
    states.reserve(size_t(frames));
    for (int frame = 0; frame < frames; ++frame) {
        sim.stepFrame();
        states.push_back(rigidStateHash(sim.rigid));
    }
    return states;
}

} // namespace

void testShockSharedSupport() {
    const uint64_t reference = movingSupport(1);
    for (int threads : {1, 2, 0, 0, 2, 0, 1, 0}) {
        const uint64_t actual = movingSupport(threads);
        std::printf("  moving support, %d threads: %016llx, serial %016llx\n", ThreadPool::instance().activeThreads(),
                    (unsigned long long)actual, (unsigned long long)reference);
        CHECK(actual == reference, "contacts sharing a moving support depend on worker scheduling");
    }
    ThreadPool::instance().setActiveThreads(0);
}

void testDensePileDeterminism() {
    for (Preset preset : {Preset::RigidGranular, Preset::RigidTeapots}) {
        const int frames = preset == Preset::RigidGranular ? 480 : 420;
        const auto reference = sampleStates(preset, frames, 1);
        for (int threads : {2, 0, 0}) {
            const auto actual = sampleStates(preset, frames, threads);
            const auto difference = std::mismatch(reference.begin(), reference.end(), actual.begin());
            const int first = difference.first == reference.end() ? -1 : int(difference.first - reference.begin()) + 1;
            std::printf("  dense preset %d, %d frames, %d threads: first differing frame %d, hash %016llx\n", int(preset), frames,
                        ThreadPool::instance().activeThreads(), first, (unsigned long long)actual.back());
            CHECK(first < 0, "dense preset %d differs from serial at frame %d", int(preset), first);
        }
    }
    ThreadPool::instance().setActiveThreads(0);
}
