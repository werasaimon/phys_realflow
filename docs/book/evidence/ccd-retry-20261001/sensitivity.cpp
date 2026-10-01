// Bounded numerical profile of the beam and the mixed CCD sample at three time steps.
// Compile with C++17, -O3 -march=native, -Isrc -I., librfsamples.a, librfcore.a and -lpthread.
#include "scene/Simulation.h"
#include "samples/Samples.h"
#include "core/Format.h"
#include <chrono>
#include <iostream>
#include <cstdlib>
using namespace rf;

static double energy(const RigidWorld& world) {
    double total = world.kineticEnergy();
    for (const RigidBody& body : world.bodies())
        if (body.mass > 0) total -= body.mass * dot(world.params.gravity, body.pos);
    return total;
}

static void beam(RigidWorld& world) {
    world.clear();
    world.setStaticMesh(nullptr);
    world.setDomain(AABB({-5, 0, -5}, {5, 10, 5}));
    world.params.sleeping = false;
    for (int i = 0; i < 6; ++i)
        world.addBox({-0.75f + 0.3f * i, 0.1f, 0}, Vector3(0.1f), Quaternion(), 500, Vector3(1));
    const int id = world.addBox({0, 1.5f, 0}, {1, 0.03f, 0.03f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.3f), 800, Vector3(1));
    world.bodies()[id].vel = {0, -25, 0};
    world.bodies()[id].angVel = {0, 40, 0};
}

static void profile(bool mixed, int factor) {
    Simulation sim;
    if (mixed) loadSample(sim, Preset::RigidCcd); else beam(sim.rigid);
    RigidWorld& world = sim.rigid;
    world.params.ccdAllMoving = std::getenv("RF_PROFILE_ALL_MOVING") != nullptr;
    const float dt = 1.0f / (600 * factor);
    const int steps = 1200 * factor;
    const double initial = energy(world);
    double peak = initial, final = initial, ccdMs = 0, solveMs = 0;
    size_t queries = 0, unresolved = 0, clamped = 0;
    float depth = 0;
    int completed = 0, longestStall = 0;
    std::vector<int> stalled(world.bodies().size(), 0);
    const auto start = std::chrono::steady_clock::now();
    for (int step = 0; step < steps; ++step) {
        const auto before = world.bodies();
        world.step(dt);
        ++completed;
        final = energy(world);
        if (!std::isfinite(final)) break;
        peak = std::max(peak, final);
        depth = std::max(depth, world.deepestPenetration());
        const auto& diag = world.ccdDiagnostics();
        queries += diag.queries; unresolved += diag.unresolved; clamped += world.ccdHits();
        ccdMs += world.timings().ccd; solveMs += world.timings().solve;
        for (size_t i = 0; i < before.size(); ++i) {
            const auto& body = world.bodies()[i]; const auto& old = before[i];
            const bool fixedPose = length2(body.pos - old.pos) == 0 && body.rot.w == old.rot.w
                && body.rot.x == old.rot.x && body.rot.y == old.rot.y && body.rot.z == old.rot.z;
            const bool fast = length(body.vel) + body.radius() * length(body.angVel) > 1;
            stalled[i] = body.invMass > 0 && fixedPose && fast ? stalled[i] + 1 : 0;
            longestStall = std::max(longestStall, stalled[i]);
        }
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "{\"scene\":\"" << (mixed ? "mixed_bullets" : "beam") << "\",\"dt_s\":" << numberText(dt)
        << ",\"completed_steps\":" << completed << ",\"expected_steps\":" << steps << ",\"finite\":" << (std::isfinite(final) ? "true" : "false")
        << ",\"initial_energy_j\":" << numberText(initial) << ",\"peak_energy_j\":" << numberText(peak)
        << ",\"final_energy_j\":" << (std::isfinite(final) ? numberText(final) : "null")
        << ",\"max_manifold_depth_m\":" << numberText(depth) << ",\"longest_fast_pose_stall_steps\":" << longestStall
        << ",\"queries\":" << queries << ",\"unresolved\":" << unresolved << ",\"clamped\":" << clamped
        << ",\"ccd_ms\":" << numberText(ccdMs) << ",\"solve_ms\":" << numberText(solveMs)
        << ",\"wall_seconds\":" << numberText(seconds) << "}\n";
}

int main() {
    for (bool mixed : {false, true}) for (int factor : {1, 2, 4}) profile(mixed, factor);
}
