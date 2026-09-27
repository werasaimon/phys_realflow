// The registry of the sample scenes: every file of samples/ adds its scenes, the list is sorted
// by Preset so that the menu number is the enum value.
#include "samples/Samples.h"

#include <algorithm>
#include <stdexcept>

namespace rf {

const std::vector<SampleEntry>& samples() {
    static const std::vector<SampleEntry> all = [] {
        std::vector<SampleEntry> v;
        addLiquidSamples(v);
        addWindTunnelSamples(v);
        addSmokeSamples(v);
        addRigidSamples(v);
        addSoftBodySamples(v);
        addHydroSamples(v);
        addFireSamples(v);
        addPlasmaSamples(v);
        std::sort(v.begin(), v.end(), [](const SampleEntry& a, const SampleEntry& b) { return a.id < b.id; });
        for (size_t i = 0; i < v.size(); ++i)
            if (int(v[i].id) != int(i)) throw std::logic_error("samples(): a Preset without a scene, or registered twice");
        if (v.size() != size_t(Preset::Count)) throw std::logic_error("samples(): the list and enum Preset disagree");
        return v;
    }();
    return all;
}

void loadSample(Simulation& sim, Preset p) {
    const SampleEntry& e = samples().at(size_t(p));
    std::unique_ptr<Scene> scene = e.create();
    scene->name = e.name;
    sim.load(std::move(scene));
}

} // namespace rf
