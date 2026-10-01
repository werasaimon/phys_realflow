// Bounded, headless refinement/load study of the stationary chain-wheel sample (preset 41).
// Build with rfsamples/rfcore; arguments: SUBSTEPS FRAMES LOAD_RATIO. JSON goes to stdout.
// This reports sampled poses and energies, not a certificate of continuous non-intersection.
#include "samples/Samples.h"
#include "samples/ChainWheelScene.h"
#include "core/Format.h"
#include "core/Parallel.h"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <iostream>
#include <stdexcept>

namespace {

double energy(const rf::Simulation& sim) {
    double result = sim.rigid.kineticEnergy();
    for (const auto& b : sim.rigid.bodies())
        if (b.invMass > 0) result -= double(b.mass) * dot(sim.rigid.params.gravity, b.pos);
    return result;
}

void hashState(uint64_t& hash, const rf::RigidWorld& world) {
    for (const auto& b : world.bodies())
        for (float v : {b.pos.x, b.pos.y, b.pos.z, b.rot.w, b.rot.x, b.rot.y, b.rot.z,
                        b.vel.x, b.vel.y, b.vel.z, b.angVel.x, b.angVel.y, b.angVel.z}) {
            uint32_t bits;
            std::memcpy(&bits, &v, sizeof(bits));
            hash = (hash ^ bits) * 1099511628211ull;
        }
}

void finalPoses(const rf::RigidWorld& world) {
    std::cout << ",\"final_poses\":[";
    bool first = true;
    for (const auto& b : world.bodies()) {
        if (!first) std::cout << ',';
        first = false;
        std::cout << '[';
        const float values[] = {b.pos.x, b.pos.y, b.pos.z, b.rot.w, b.rot.x, b.rot.y, b.rot.z};
        for (int k = 0; k < 7; ++k) std::cout << (k ? "," : "") << rf::numberText(values[k]);
        std::cout << ']';
    }
    std::cout << ']';
}

void study(int substeps, int frames, float load) {
    rf::Simulation sim;
    rf::loadSample(sim, rf::Preset::ChainWheel);
    sim.scene()->setParam(0, load);
    sim.rigid.params.substeps = substeps;
    sim.reset();
    auto* scene = static_cast<const rf::ChainWheelScene*>(sim.scene());
    const double initial = energy(sim);
    double peak = initial, lateT = 0, sumMs = 0;
    float contactDepth = 0;
    std::vector<double> ms;
    uint64_t hash = 14695981039346656037ull;
    for (int f = 0; f < frames; ++f) {
        const auto start = std::chrono::steady_clock::now();
        sim.stepFrame();
        ms.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count());
        sumMs += ms.back();
        hashState(hash, sim.rigid);
        const double e = energy(sim);
        if (!std::isfinite(e) || scene->firstBrokenFrame() >= 0) throw std::runtime_error("non-finite energy or chain break");
        peak = std::max(peak, e);
        contactDepth = std::max(contactDepth, sim.rigid.deepestPenetration());
        if (f >= frames - 60) lateT += sim.rigid.kineticEnergy() / 60.0;
    }
    std::sort(ms.begin(), ms.end());
    std::cout << "{\"substeps\":" << substeps << ",\"frames\":" << frames << ",\"load_ratio\":" << rf::numberText(load)
        << ",\"threads\":" << rf::ThreadPool::instance().threadCount() << ",\"broken_pairs\":" << scene->brokenPairs()
        << ",\"initial_energy_J\":" << rf::numberText(initial) << ",\"peak_energy_J\":" << rf::numberText(peak)
        << ",\"final_energy_J\":" << rf::numberText(energy(sim)) << ",\"late_kinetic_J\":" << rf::numberText(lateT)
        << ",\"solver_contact_depth_m\":" << rf::numberText(contactDepth) << ",\"mean_step_ms\":" << rf::numberText(sumMs / frames)
        << ",\"p95_step_ms\":" << rf::numberText(ms[size_t(frames * 0.95)])
        << ",\"trajectory_hash\":\"" << std::hex << hash << std::dec << '"';
    finalPoses(sim.rigid);
    std::cout << "}\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        float substeps = 0, frames = 0, load = 0;
        if (argc != 4 || !rf::parseNumber(argv[1], substeps) || !rf::parseNumber(argv[2], frames)
            || !rf::parseNumber(argv[3], load) || substeps < 1 || substeps > 160 || frames < 120 || frames > 3600
            || !std::isfinite(substeps) || !std::isfinite(frames) || substeps != std::floor(substeps) || frames != std::floor(frames)
            || !std::isfinite(load) || load < 1 || load > 3)
            throw std::runtime_error("SUBSTEPS (1..160) FRAMES (120..3600) LOAD_RATIO (1..3)");
        study(int(substeps), int(frames), load);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 2;
    }
}
