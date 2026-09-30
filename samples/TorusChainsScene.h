#pragma once
// Contact-only chains of circular tori: a loaded hanging chain, a sideways release and a bridge.
// Their centreline linking number is observed without applying any forces or constraints.
#include "scene/Scene.h"
#include "rigid/RigidWorld.h"

namespace rf {

class TorusChainsScene final : public Scene {
public:
    static constexpr float majorRadius = 0.05f, tubeRadius = 0.012f, spacing = 0.07f;
    void configure(Simulation& sim) override;
    void build(Simulation& sim) override;
    void afterStep(Simulation& sim) override;
    void describe(const Simulation& sim, RenderSnapshot& out) const override;
    std::vector<SceneParam> params() const override;
    void setParam(int index, float value) override;
    // Signed disk crossings for two ringShape() bodies; not an arbitrary compound-shape query.
    static int linkingNumber(const RigidBody& a, const RigidBody& b);
    int linksPerChain() const { return count_; }
    int brokenPairs() const { return broken_; }
    int firstBrokenFrame() const { return firstBrokenFrame_; }
private:
    void addChain(Simulation& sim, int kind, float length, float height);
    int count_ = 40, broken_ = 0, firstBrokenFrame_ = -1;
    float loadRatio_ = 10, kick_ = 1, maxSeparation_ = 0;
    std::vector<int> initialLinking_;
};

void addTorusChainSamples(std::vector<struct SampleEntry>& out);

} // namespace rf
