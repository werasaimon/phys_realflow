#pragma once
// A short, loaded torus chain over a toothed wheel: fixed support or a torque-limited hinge drive.
// The driven mode frees both chain ends. Diagnostics observe, never repair, ring linking.
#include "scene/Scene.h"

#include <array>

namespace rf {

class ChainWheelScene final : public Scene {
public:
    enum class Mode { Stationary, Driven };
    explicit ChainWheelScene(Mode mode = Mode::Stationary) : mode_(mode) {}
    static constexpr int toothCount = 10, linkCount = toothCount + 6;
    static constexpr int wheelBody = 0, firstLink = 1, lastLink = linkCount;
    static float pitchRadius();
    void configure(Simulation& sim) override;
    void build(Simulation& sim) override;
    void afterStep(Simulation& sim) override;
    void describe(const Simulation& sim, RenderSnapshot& out) const override;
    std::vector<SceneParam> params() const override;
    void setParam(int index, float value) override;
    int brokenPairs() const { return broken_; }
    int firstBrokenFrame() const { return firstBrokenFrame_; }
    bool driven() const { return mode_ == Mode::Driven; }
private:
    void addChain(Simulation& sim);
    void updateMotor(Simulation& sim);
    std::array<int, linkCount - 1> initialLinking_{};
    int broken_ = 0, firstBrokenFrame_ = -1;
    float loadRatio_ = 2;
    Mode mode_;
    float motorSpeed_ = 0.15f, maxMotorTorque_ = 20;
    int phase_ = 0;
};

void addChainWheelSample(std::vector<struct SampleEntry>& out);

} // namespace rf
