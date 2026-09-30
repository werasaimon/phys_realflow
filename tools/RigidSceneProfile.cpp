// Profile a rigid sample through the same Simulation path as the editor. CSV stage timers are
// separate from geometry observations and file output; final state hashes compare exact runs.
// Build: c++ -std=c++17 -O3 -DNDEBUG -march=native -Isrc -I. tools/RigidSceneProfile.cpp
//        build-core/librfsamples.a build-core/librfcore.a -pthread -o /tmp/rf_scene_profile
// Run: RF_THREADS=4 /tmp/rf_scene_profile 38 600 -1 /tmp/sandbox-profile
// Arguments: preset, frames, sleeping (-1 = scene default, 0 = off, 1 = on), new output directory.
#include "samples/Samples.h"
#include "core/Format.h"
#include "core/Probe.h"

#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace rf;
namespace {

int integer(const char* text, int lo, int hi) {
    float value = 0;
    if (!parseNumber(text, value) || !std::isfinite(value) || value < lo || value > hi || std::floor(value) != value)
        throw std::runtime_error("invalid integer argument");
    return int(value);
}

void row(std::ofstream& stream, std::initializer_list<double> values) {
    bool first = true;
    for (double value : values) {
        if (!first) stream << ',';
        stream << numberText(value);
        first = false;
    }
    stream << '\n';
}

bool finite(const RigidBody& body) {
    for (const Vector3& v : {body.pos, body.vel, body.angVel})
        if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z)) return false;
    return std::isfinite(body.rot.w) && std::isfinite(body.rot.x) && std::isfinite(body.rot.y) && std::isfinite(body.rot.z);
}

uint64_t stateHash(const RigidWorld& world) {
    uint64_t hash = 14695981039346656037ull;
    auto add = [&](float value) {
        unsigned char bytes[sizeof(float)];
        std::memcpy(bytes, &value, sizeof(value));
        for (unsigned char b : bytes) hash = (hash ^ b) * 1099511628211ull;
    };
    for (const RigidBody& b : world.bodies()) {
        for (const Vector3& v : {b.pos, b.vel, b.angVel}) { add(v.x); add(v.y); add(v.z); }
        add(b.rot.w); add(b.rot.x); add(b.rot.y); add(b.rot.z);
        add(b.sleepTimer); add(float(b.sleeping));
    }
    return hash;
}

void profile(Simulation& sim, int frames, const std::filesystem::path& output) {
    std::ofstream csv(output / "frames.csv"), summary(output / "summary.txt");
    if (!csv || !summary) throw std::runtime_error("cannot open profile output");
    csv << "frame,step_ms,collide_ms,solve_ms,ccd_ms,integrate_ms,islands_ms,awake,contacts,manifolds,energy_j,max_floor_depth_m\n";
    float floor = 0;
    for (int f = 0; f < frames; ++f) {
        sim.stepFrame();
        const auto probe = Probe::snapshot();
        double energy = sim.rigid.kineticEnergy();
        int awake = 0;
        for (const auto& body : sim.rigid.bodies()) {
            if (!body.alive || body.mass <= 0) continue;
            if (!finite(body)) throw std::runtime_error("non-finite rigid state");
            awake += !body.sleeping;
            energy -= double(body.mass) * dot(sim.rigid.params.gravity, body.pos);
            floor = std::max(floor, sim.rigid.domain().lo.y - body.posed().support({0, -1, 0}).y);
        }
        row(csv, {double(f + 1), probe.value("frame/step ms"), probe.value("rigid/collide ms"), probe.value("rigid/solve ms"),
            probe.value("rigid/ccd ms"), probe.value("rigid/integrate ms"), probe.value("rigid/islands ms"),
            double(awake), double(sim.rigid.contactCount()), double(sim.rigid.manifoldCount()), energy, floor});
    }
    summary << "frames=" << numberText(double(frames)) << "\nsubsteps=" << numberText(double(sim.rigid.params.substeps))
            << "\niterations=" << numberText(double(sim.rigid.params.iterations)) << "\nslop_m=" << numberText(sim.rigid.params.slop)
            << "\nsleeping=" << (sim.rigid.params.sleeping ? "1" : "0") << "\nstate_hash=" << std::hex << stateHash(sim.rigid) << '\n';
    if (!csv || !summary) throw std::runtime_error("failed writing profile output");
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 5) throw std::runtime_error("usage: rf_scene_profile PRESET FRAMES SLEEP(-1/0/1) NEW_OUTPUT_DIR");
        const int preset = integer(argv[1], 0, int(Preset::Count) - 1), frames = integer(argv[2], 1, 1000000);
        const int sleep = integer(argv[3], -1, 1);
        const std::filesystem::path output(argv[4]);
        Simulation sim;
        loadSample(sim, Preset(preset));
        if (sim.mode() != SimMode::Rigid || sim.particles.size() != 0) throw std::runtime_error("sample must contain only rigid bodies");
        if (sleep >= 0) sim.rigid.params.sleeping = sleep != 0;
        if (std::filesystem::exists(output) || !std::filesystem::create_directories(output))
            throw std::runtime_error("choose a new output directory");
        profile(sim, frames, output);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
