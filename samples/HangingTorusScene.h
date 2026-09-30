#pragma once
// A small contact-only hanging chain for interactive mouse pulls. The topology meter observes
// the original links; it never moves bodies or adds constraints between rings.
#include "scene/Scene.h"

#include <array>

namespace rf {

class HangingTorusScene final : public Scene {
public:
    static constexpr int linkCount = 8;
    void configure(Simulation& sim) override;
    void build(Simulation& sim) override;
    void afterStep(Simulation& sim) override;
    void describe(const Simulation& sim, RenderSnapshot& out) const override;
    int brokenPairs() const { return broken_; }
    int firstBrokenFrame() const { return firstBrokenFrame_; }
private:
    std::array<int, linkCount - 1> initialLinking_{};
    int broken_ = 0, firstBrokenFrame_ = -1;
    float maxSeparation_ = 0;
};

void addHangingTorusSample(std::vector<struct SampleEntry>& out);

} // namespace rf
