// Controlled shock-propagation A/B: identical hanging-torus geometry and initial tip kick.
// Measure topology, sampled part overlap and energy without changing solver tolerances.
#include "tests/TestRunner.h"
#include "samples/HangingTorusScene.h"
#include "core/Format.h"
#include <iostream>

double totalEnergy(const RigidWorld& world) {
    double result = world.kineticEnergy();
    for (const auto& body : world.bodies())
        if (body.alive && body.invMass > 0)
            result -= double(body.mass) * dot(world.params.gravity, body.pos);
    return result;
}

void runProbe(bool shock, int substeps) {
    Simulation sim;
    loadSample(sim, Preset::HangingTorus);
    sim.rigid.params.shockPropagation = shock;
    sim.rigid.params.substeps = substeps;
    sim.rigid.bodies().back().vel = {0.75f, 0, 0.3f};
    const auto* scene = static_cast<const HangingTorusScene*>(sim.scene());
    const double initial = totalEnergy(sim.rigid);
    double peak = initial, peakStepGain = 0, previous = initial, lateKinetic = 0;
    float overlap = maxPartOverlap(sim.rigid), lateOverlap = 0;
    int broken = 0;
    const auto begin = std::chrono::steady_clock::now();
    for (int frame = 0; frame < 600; ++frame) {
        sim.stepFrame();
        const double energy = totalEnergy(sim.rigid);
        peak = std::max(peak, energy);
        peakStepGain = std::max(peakStepGain, energy - previous);
        previous = energy;
        broken = std::max(broken, scene->brokenPairs());
        if (frame >= 480) lateKinetic += sim.rigid.kineticEnergy() / 120;
        if (frame % 10 == 0) {
            const float depth = maxPartOverlap(sim.rigid);
            overlap = std::max(overlap, depth);
            if (frame >= 480) lateOverlap = std::max(lateOverlap, depth);
        }
    }
    const Vector3 tip = sim.rigid.bodies().back().pos;
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    std::cout << "{\"shock\":" << (shock ? "true" : "false")
              << ",\"substeps\":" << substeps << ",\"frames\":600,\"broken_pairs\":" << broken
              << ",\"first_broken_frame\":" << scene->firstBrokenFrame()
              << ",\"initial_energy_J\":" << numberText(initial)
              << ",\"final_energy_J\":" << numberText(previous)
              << ",\"peak_energy_gain_J\":" << numberText(peak - initial)
              << ",\"peak_frame_energy_gain_J\":" << numberText(peakStepGain)
              << ",\"late_mean_kinetic_J\":" << numberText(lateKinetic)
              << ",\"sampled_peak_overlap_m\":" << numberText(overlap)
              << ",\"sampled_late_overlap_m\":" << numberText(lateOverlap)
              << ",\"final_tip_m\":[" << numberText(tip.x) << ',' << numberText(tip.y) << ',' << numberText(tip.z)
              << "],\"wall_seconds\":" << numberText(seconds) << "}\n" << std::flush;
}

int main() {
    for (int substeps : {40, 80})
        for (bool shock : {true, false}) runProbe(shock, substeps);
}
