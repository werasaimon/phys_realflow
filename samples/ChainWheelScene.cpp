// Chain-wheel fixtures: one faceted hub, ten slender teeth and sixteen existing tori.
// This is an original sample, not an import of the Rigid IPC bicycle-chain benchmark.
#include "samples/ChainWheelScene.h"
#include "samples/Models.h"
#include "samples/Samples.h"
#include "samples/TorusChainsScene.h"
#include "core/Format.h"
#include "core/Probe.h"

#include <algorithm>
#include <cmath>

namespace rf {
namespace {

const Vector3 wheelCentre(0, 0.55f, 0);

std::shared_ptr<const CompoundShape> chainWheelShape() {
    static const auto shape = [] {
        const float radius = ChainWheelScene::pitchRadius(), root = radius - 0.066f, tip = radius + 0.005f;
        std::vector<TriMesh> parts{primitives::cylinder(root, 0.018f, 40)};
        for (int k = 0; k < ChainWheelScene::toothCount; ++k) {
            const float a = kPi * float(2 * k + 1) / ChainWheelScene::toothCount, mid = 0.5f * (root + tip);
            TriMesh tooth = primitives::box({0.5f * (tip - root), 0.004f, 0.008f});
            tooth.transform(Quaternion::fromAxisAngle({0, 0, 1}, a).toMatrix3x3(), Vector3(1),
                            {mid * rf::cos(a), mid * rf::sin(a), 0});
            parts.push_back(std::move(tooth));
        }
        return std::make_shared<const CompoundShape>(parts, primitives::merge(parts));
    }();
    return shape;
}

} // namespace

float ChainWheelScene::pitchRadius() {
    // Adjacent link centres form a chord s = 2 R sin(alpha/2); teeth engage alternate rings.
    return TorusChainsScene::spacing / (2 * rf::sin(kPi / (2 * toothCount)));
}

void ChainWheelScene::configure(Simulation& sim) {
    sim.rigid.params.sleeping = false;
    sim.rigid.params.substeps = 40;
    sim.rigid.params.iterations = 12;
    sim.rigid.params.slop = 0.0001f;
}

void ChainWheelScene::addChain(Simulation& sim) {
    const float radius = pitchRadius(), spacing = TorusChainsScene::spacing;
    const auto ring = ringShape();
    // Two links precede the semicircle, three follow it. Only the stationary mode anchors an end.
    for (int k = -2; k <= toothCount + 3; ++k) {
        const float a = k < 0 ? kPi : k > toothCount ? 0 : kPi - float(k) * kPi / toothCount;
        Vector3 at = wheelCentre + Vector3(radius * rf::cos(a), radius * rf::sin(a), 0);
        if (k < 0) at.y += float(k) * spacing;
        if (k > toothCount) at.y -= float(k - toothCount) * spacing;
        const Quaternion pose = k % 2 ? Quaternion::fromAxisAngle({0, 0, 1}, a - 0.5f * kPi)
            : Quaternion::fromAxisAngle({1, 0, 0}, 0.5f * kPi);
        const bool fixed = !driven() && k == -2, loaded = k == toothCount + 3;
        const Vector3 color = fixed ? Vector3(0.95f, 0.25f, 0.15f) : loaded ? Vector3(1, 0.7f, 0.15f)
            : Vector3(0.4f, 0.7f, 0.95f);
        sim.rigid.addCompound(ring, at, pose, fixed ? 0 : 7800.0f * (loaded ? loadRatio_ : 1), color);
    }
}

void ChainWheelScene::build(Simulation& sim) {
    sim.useRigidArena(AABB({-0.65f, 0, -0.4f}, {0.65f, 1.05f, 0.4f}));
    sim.rigid.addCompound(chainWheelShape(), wheelCentre, Quaternion(), driven() ? 7800.0f : 0, {0.6f, 0.62f, 0.66f});
    addChain(sim);
    phase_ = 0;
    if (driven()) {
        auto& motor = sim.rigid.addHingeJoint(wheelBody, -1, wheelCentre, {0, 0, 1});
        motor.motorEnabled = true;
        motor.motorSpeed = 0;
        motor.maxMotorTorque = maxMotorTorque_;
    }
    broken_ = 0;
    firstBrokenFrame_ = -1;
    const auto& bodies = sim.rigid.bodies();
    for (int i = 1; i < linkCount; ++i)
        initialLinking_[size_t(i - 1)] = TorusChainsScene::linkingNumber(bodies[size_t(i)], bodies[size_t(i + 1)]);
}

void ChainWheelScene::afterStep(Simulation& sim) {
    broken_ = 0;
    const auto& bodies = sim.rigid.bodies();
    for (int i = 1; i < linkCount; ++i) {
        const int initial = initialLinking_[size_t(i - 1)];
        broken_ += std::abs(initial) != 1 || TorusChainsScene::linkingNumber(bodies[size_t(i)], bodies[size_t(i + 1)]) != initial;
    }
    if (broken_ && firstBrokenFrame_ < 0) firstBrokenFrame_ = int(sim.frame()) + 1;
    Probe::set("chain-wheel/broken pairs", broken_);
    Probe::set("chain-wheel/first broken frame", firstBrokenFrame_);
    if (driven()) updateMotor(sim);
}

void ChainWheelScene::updateMotor(Simulation& sim) {
    auto& motor = static_cast<HingeJoint&>(*sim.rigid.joints().front());
    // A 10 s cycle: settle/brake, lift, hold, lower, hold. Speed changes are commands, not poses.
    phase_ = int(std::fmod(sim.time(), 10.0) / 2.0);
    motor.motorSpeed = phase_ == 1 ? motorSpeed_ : phase_ == 3 ? -motorSpeed_ : 0;
    Probe::set("chain-wheel/phase", phase_);
    Probe::set("chain-wheel/angle rad", motor.angle(sim.rigid.bodies()));
    Probe::set("chain-wheel/angular speed rad_s", sim.rigid.bodies()[wheelBody].angVel.z);
}

void ChainWheelScene::describe(const Simulation& sim, RenderSnapshot& out) const {
    out.info.push_back({"Опора", driven() ? "Колесо на осевом шарнире; оба конца цепи свободны"
                                        : "Звёздочка и красное возвратное звено неподвижны; привод выключен"});
    out.info.push_back({"Торы / шарниры", format("%d / %d", linkCount, driven() ? 1 : 0)});
    out.info.push_back({"Зацепления сейчас", format("%d / %d", linkCount - 1 - broken_, linkCount - 1)});
    out.info.push_back({"Первый разрыв, кадр", format("%d", firstBrokenFrame_)});
    out.info.push_back({"Груз", format("Масса жёлтого звена: %.1f массы обычного", loadRatio_)});
    out.info.push_back({"Мышь", "ЛКМ по ободу жёлтого звена — тянуть; отпустить — свободное движение"});
    if (driven()) {
        const char* phases[] = {"Пауза с торможением", "Подъём", "Удержание", "Обратный ход", "Удержание"};
        out.info.push_back({"Привод", phases[phase_]});
        out.info.push_back({"Скорость / предел момента", format("%.2f рад/с / %.1f Н·м", motorSpeed_, maxMotorTorque_)});
        out.info.push_back({"Вращение сейчас", format("%.3f рад/с", sim.rigid.bodies()[wheelBody].angVel.z)});
    }
}

std::vector<SceneParam> ChainWheelScene::params() const {
    std::vector<SceneParam> out{{"Масса груза / звена", loadRatio_, 1, 3, 0.5f, 1, false, "Масса жёлтого кольца при прежней геометрии."}};
    if (driven()) {
        out.push_back({"Скорость привода", motorSpeed_, 0.05f, 0.15f, 0.05f, 2, false, "рад/с; подъём и обратный ход по 2 с, между ними удержание."});
        out.push_back({"Предел момента", maxMotorTorque_, 5, 30, 5, 0, false, "Н·м; при недостаточном моменте привод может остановиться под нагрузкой."});
    }
    return out;
}

void ChainWheelScene::setParam(int index, float value) {
    if (!std::isfinite(value)) return;
    if (index == 0) loadRatio_ = std::clamp(value, 1.0f, 3.0f);
    if (driven() && index == 1) motorSpeed_ = std::clamp(value, 0.05f, 0.15f);
    if (driven() && index == 2) maxMotorTorque_ = std::clamp(value, 5.0f, 30.0f);
}

void addChainWheelSample(std::vector<SampleEntry>& out) {
    out.push_back({Preset::ChainWheel, "Твёрдые тела", "Цепь на неподвижной звёздочке: 16 торов и груз",
                   [] { return std::unique_ptr<Scene>(new ChainWheelScene); }});
    out.push_back({Preset::ChainDrive, "Твёрдые тела", "Цепной привод: подъём, удержание и обратный ход",
                   [] { return std::unique_ptr<Scene>(new ChainWheelScene(ChainWheelScene::Mode::Driven)); }});
}

} // namespace rf
