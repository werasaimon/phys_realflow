// Reproduce strict-FP residual cup motion; vary one solver mechanism per run.
#include "tests/TestRunner.h"
#include "core/Format.h"
#include <iostream>

void probe(const std::string& variant, bool sleep) {
    Simulation sim;
    loadSample(sim, Preset::RigidSandbox);
    auto& p = sim.rigid.params;
    p.substeps = 10;
    p.iterations = 6;
    p.slop = 0.004f;
    p.sleeping = sleep;
    if (variant == "iterations12") p.iterations = 12;
    if (variant == "iterations24") p.iterations = 24;
    if (variant == "substeps20") p.substeps = 20;
    if (variant == "substeps40") p.substeps = 40;
    if (variant == "shock-off") p.shockPropagation = false;
    if (variant == "shock-friction-off") p.shockFriction = false;
    if (variant == "lock-off") p.rotationalLock = false;
    if (variant == "rolling-off") p.rollingResistance = 0;
    if (variant == "slop1mm") p.slop = 0.001f;
    if (variant == "slop0.1mm" || variant == "substeps20-slop0.1mm") p.slop = 0.0001f;
    if (variant == "substeps20-slop0.1mm") p.substeps = 20;
    StateHash hash;
    double lateKinetic = 0, lastEnergy = 0, lateGain = 0;
    float spin = 0, speed = 0, lateFloor = 0;
    int allSleep = -1;
    std::vector<Vector3> last;
    std::vector<float> path(sim.rigid.bodies().size(), 0);
    for (const auto& b : sim.rigid.bodies()) last.push_back(b.pos);
    for (int frame = 0; frame < 600; ++frame) {
        sim.stepFrame();
        double energy = sim.rigid.kineticEnergy();
        if (frame >= 480) lateKinetic += energy / 120;
        for (size_t i = 0; i < sim.rigid.bodies().size(); ++i) {
            const auto& b = sim.rigid.bodies()[i];
            energy -= double(b.mass) * dot(p.gravity, b.pos);
            hash.add(b.pos); hash.add(b.rot); hash.add(b.vel); hash.add(b.angVel);
            if (frame >= 480) {
                path[i] += length(b.pos - last[i]);
                spin = std::max(spin, length(b.angVel)); speed = std::max(speed, length(b.vel));
                lateFloor = std::max(lateFloor, -b.posed().support(Vector3(0,-1,0)).y);
            }
            last[i] = b.pos;
        }
        if (frame >= 480) lateGain = std::max(lateGain, energy - lastEnergy);
        lastEnergy = energy;
        if (sim.rigid.sleepingCount() == 24 && allSleep < 0) allSleep = frame + 1;
    }
    const auto& cup = sim.rigid.bodies()[23];
    const float overlap = maxPartOverlap(sim.rigid);
    std::cout << "{\"variant\":\"" << variant << "\",\"sleep_enabled\":" << (sleep ? "true" : "false")
              << ",\"asleep\":" << sim.rigid.sleepingCount() << ",\"all_sleep_frame\":" << allSleep
              << ",\"trajectory_hash\":\"" << std::hex << hash.h << std::dec << "\""
              << ",\"late_mean_kinetic_J\":" << numberText(lateKinetic)
              << ",\"late_max_frame_energy_gain_J\":" << numberText(lateGain)
              << ",\"late_peak_spin_rad_s\":" << numberText(spin)
              << ",\"late_peak_speed_m_s\":" << numberText(speed)
              << ",\"late_max_path_m\":" << numberText(*std::max_element(path.begin(), path.end()))
              << ",\"final_overlap_m\":" << numberText(overlap)
              << ",\"late_floor_depth_m\":" << numberText(lateFloor)
              << ",\"cup_spin_rad_s\":" << numberText(length(cup.angVel))
              << ",\"cup_speed_m_s\":" << numberText(length(cup.vel)) << "}\n" << std::flush;
}

int main(int argc, char** argv) {
    const std::string selected = argc > 1 ? argv[1] : "default";
    if (selected == "geometry") {
        for (const char* variant : {"default", "slop1mm", "slop0.1mm", "substeps20-slop0.1mm"})
            for (bool sleep : {true, false}) probe(variant, sleep);
        return 0;
    }
    if (selected != "all") { probe(selected, true); probe(selected, false); return 0; }
    for (const char* variant : {"default", "iterations12", "iterations24", "substeps20", "substeps40",
                                "shock-off", "shock-friction-off", "lock-off", "rolling-off"})
        for (bool sleep : {true, false}) probe(variant, sleep);
}
