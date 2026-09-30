// Eight actual non-convex tori, one fixed at the top, for a responsive mouse-grab experiment.
// Geometry and contact settings match the long-chain stress sample; fewer bodies make it cheap.
#include "samples/HangingTorusScene.h"
#include "samples/Models.h"
#include "samples/Samples.h"
#include "samples/TorusChainsScene.h"
#include "core/Format.h"
#include "core/Probe.h"

#include <algorithm>
#include <cmath>

namespace rf {

void HangingTorusScene::configure(Simulation& sim) {
    sim.rigid.params.sleeping = false;
    sim.rigid.params.substeps = 40;
    sim.rigid.params.iterations = 12;
    sim.rigid.params.slop = 0.0001f;
}

void HangingTorusScene::build(Simulation& sim) {
    sim.useRigidArena(AABB({-0.45f, 0, -0.45f}, {0.45f, 0.9f, 0.45f}));
    const auto ring = ringShape();
    const Quaternion xy = Quaternion::fromAxisAngle({1, 0, 0}, 0.5f * kPi);
    const Quaternion yz = Quaternion::fromAxisAngle({0, 0, 1}, -0.5f * kPi);
    for (int k = 0; k < linkCount; ++k) {
        const Vector3 color = k == 0 ? Vector3(0.95f, 0.25f, 0.15f)
            : k == linkCount - 1 ? Vector3(1, 0.7f, 0.15f) : Vector3(0.4f, 0.7f, 0.95f);
        sim.rigid.addCompound(ring, {0, 0.72f - TorusChainsScene::spacing * float(k), 0},
                              k % 2 ? yz : xy, k == 0 ? 0.0f : 7800.0f, color);
    }
    broken_ = 0;
    firstBrokenFrame_ = -1;
    maxSeparation_ = 0;
    const auto& bodies = sim.rigid.bodies();
    for (int k = 1; k < linkCount; ++k)
        initialLinking_[size_t(k - 1)] = TorusChainsScene::linkingNumber(bodies[size_t(k - 1)], bodies[size_t(k)]);
}

void HangingTorusScene::afterStep(Simulation& sim) {
    broken_ = 0;
    const auto& bodies = sim.rigid.bodies();
    for (int k = 1; k < linkCount; ++k) {
        const auto& a = bodies[size_t(k - 1)];
        const auto& b = bodies[size_t(k)];
        const int initial = initialLinking_[size_t(k - 1)];
        broken_ += std::abs(initial) != 1 || TorusChainsScene::linkingNumber(a, b) != initial;
        maxSeparation_ = std::max(maxSeparation_, length(a.pos - b.pos));
    }
    if (broken_ && firstBrokenFrame_ < 0) firstBrokenFrame_ = int(sim.frame()) + 1;
    Probe::set("chains/broken pairs", broken_);
    Probe::set("chains/first broken frame", firstBrokenFrame_);
    Probe::set("chains/max separation m", maxSeparation_);
}

void HangingTorusScene::describe(const Simulation&, RenderSnapshot& out) const {
    out.info.push_back({"Мышь", "ЛКМ по ободу звена и тянуть; отпустить — свободное движение"});
    out.info.push_back({"Торы / шарниры", format("%d / 0", linkCount)});
    out.info.push_back({"Зацепления сейчас", format("%d / %d", linkCount - 1 - broken_, linkCount - 1)});
    out.info.push_back({"Первый разрыв, кадр", format("%d", firstBrokenFrame_)});
    out.info.push_back({"Макс. расстояние соседей", format("%.2f мм", 1000 * maxSeparation_)});
}

void addHangingTorusSample(std::vector<SampleEntry>& out) {
    out.push_back({Preset::HangingTorus, "Твёрдые тела", "Одна цепь: 8 торов — потяните звено мышью",
                   [] { return std::unique_ptr<Scene>(new HangingTorusScene); }});
}

} // namespace rf
