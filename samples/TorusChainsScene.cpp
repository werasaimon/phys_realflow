// Three long chains of true non-convex tori. Only the endpoint rings are fixed; contact with
// convex parts holds every other ring. Diagnostics measure linking; they never correct a pose.
#include "samples/TorusChainsScene.h"
#include "samples/Models.h"
#include "samples/Samples.h"
#include "core/Format.h"
#include "core/Probe.h"

#include <algorithm>
#include <cmath>

namespace rf {
namespace {

struct Circle {
    Vector3 centre, x, z, normal;
};

Circle circleOf(const RigidBody& body) {
    const auto& shape = static_cast<const CompoundShape&>(*body.shape);
    const Matrix3x3 rotation = body.rotation() * shape.principalRotation().transposed();
    return {body.pos - rotation * shape.centerOfMass(), rotation * Vector3(1, 0, 0),
            rotation * Vector3(0, 0, 1), rotation * Vector3(0, 1, 0)};
}

} // namespace

// Linking is the signed intersection count of B's centreline with A's spanning disk.
// B(t)=centre+R(x cos(t)+z sin(t)); its plane crossings solve a cos(t)+b sin(t)+c=0.
// Both crossings inside the disk cancel. Distance alone cannot distinguish these unlinked rings.
// Oriented intersections: G. Moore, Physics 511 (2014), topology notes, pp. 8-9:
// https://www.physics.rutgers.edu/~gmoore/511Fall2014/Physics511-2014-Ch2-Topology.pdf
int TorusChainsScene::linkingNumber(const RigidBody& a, const RigidBody& b) {
    const Circle first = circleOf(a), second = circleOf(b);
    const double ca = majorRadius * double(dot(first.normal, second.x));
    const double cb = majorRadius * double(dot(first.normal, second.z));
    const double cc = dot(first.normal, second.centre - first.centre);
    const double amplitude = std::hypot(ca, cb);
    if (amplitude < 1e-9 || std::fabs(cc) >= amplitude) return 0;
    const double phase = std::atan2(cb, ca), delta = std::acos(-cc / amplitude);
    int crossings = 0;
    for (double t : {phase - delta, phase + delta}) {
        const Vector3 at = second.centre + majorRadius * (second.x * float(std::cos(t)) + second.z * float(std::sin(t)));
        if (double(length2(at - first.centre)) < double(majorRadius) * majorRadius)
            crossings += -ca * std::sin(t) + cb * std::cos(t) > 0 ? 1 : -1;
    }
    return crossings;
}

void TorusChainsScene::configure(Simulation& sim) {
    sim.rigid.params.sleeping = false; // do not freeze a failure or conceal residual motion
    // The tube is only 24 mm thick. Refresh geometry at 2400 Hz and keep the allowed
    // penetration small relative to that thickness; the generic 4 mm slop loses links.
    sim.rigid.params.substeps = 40;
    sim.rigid.params.iterations = 12;
    sim.rigid.params.slop = 0.0001f;
}

void TorusChainsScene::addChain(Simulation& sim, int kind, float length, float height) {
    const auto ring = ringShape();
    // Short bridges start with a shallow arc: bending adjacent thick tori too far intersects
    // their tubes before release. The default 40-link bridge still spans a semicircle.
    const float arc = std::min(kPi, 0.1f * float(count_ - 1));
    const float arcStep = arc / float(count_ - 1), radius = spacing / (2 * rf::sin(0.5f * arcStep));
    const float start = arc == kPi ? kPi : 1.5f * kPi - 0.5f * arc;
    const float lift = arc == kPi ? 0 : radius * rf::cos(0.5f * arc);
    const Quaternion inXY = Quaternion::fromAxisAngle({1, 0, 0}, 0.5f * kPi);
    for (int k = 0; k < count_; ++k) {
        Vector3 at;
        float angle;
        if (kind == 0) {
            at = {-0.45f * length, height - spacing * float(k), 0.65f};
            angle = -0.5f * kPi;
        } else if (kind == 1) {
            at = {-0.1f * length + spacing * float(k), height, 0.65f};
            angle = 0;
        } else {
            const float t = start + arcStep * float(k);
            at = {radius * rf::cos(t), height - 0.4f + lift + radius * rf::sin(t), -0.65f};
            angle = t + 0.5f * kPi;
        }
        const bool fixed = k == 0 || (kind == 2 && k == count_ - 1);
        const bool loaded = kind == 2 ? k == count_ / 2 : k == count_ - 1;
        const float density = fixed ? 0 : 7800.0f * (loaded ? loadRatio_ : 1);
        const Vector3 color = fixed ? Vector3(0.9f, 0.25f, 0.15f) : loaded ? Vector3(1, 0.65f, 0.15f)
            : kind == 0 ? Vector3(0.5f, 0.7f, 0.9f) : kind == 1 ? Vector3(0.85f, 0.6f, 0.35f) : Vector3(0.6f, 0.8f, 0.55f);
        const Quaternion pose = k % 2 ? Quaternion::fromAxisAngle({0, 0, 1}, angle) : inXY;
        const int id = sim.rigid.addCompound(ring, at, pose, density, color);
        if (loaded && !fixed) sim.rigid.bodies()[size_t(id)].vel.z = kick_;
    }
}

void TorusChainsScene::build(Simulation& sim) {
    const float length = spacing * float(count_ - 1), height = length + 1;
    sim.useRigidArena(AABB({-length - 0.5f, 0, -1.5f}, {length + 0.5f, height + 0.4f, 1.5f}));
    for (int kind = 0; kind < 3; ++kind) addChain(sim, kind, length, height);
    broken_ = 0; firstBrokenFrame_ = -1; maxSeparation_ = 0;
    initialLinking_.clear();
    const auto& bodies = sim.rigid.bodies();
    for (int root : {0, count_, 2 * count_})
        for (int k = 1; k < count_; ++k)
            initialLinking_.push_back(linkingNumber(bodies[size_t(root + k - 1)], bodies[size_t(root + k)]));
}

void TorusChainsScene::afterStep(Simulation& sim) {
    const auto& bodies = sim.rigid.bodies();
    broken_ = 0;
    size_t pair = 0;
    for (int root : {0, count_, 2 * count_})
        for (int k = 1; k < count_; ++k, ++pair) {
            const auto& a = bodies[size_t(root + k - 1)];
            const auto& b = bodies[size_t(root + k)];
            broken_ += std::abs(initialLinking_[pair]) != 1 || linkingNumber(a, b) != initialLinking_[pair];
            maxSeparation_ = std::max(maxSeparation_, length(a.pos - b.pos));
        }
    if (broken_ && firstBrokenFrame_ < 0) firstBrokenFrame_ = int(sim.frame()) + 1;
    Probe::set("chains/broken pairs", broken_);
    Probe::set("chains/first broken frame", firstBrokenFrame_);
    Probe::set("chains/max separation m", maxSeparation_);
}

void TorusChainsScene::describe(const Simulation&, RenderSnapshot& out) const {
    out.info.push_back({"Торы / шарниры", format("%d / 0", 3 * count_)});
    out.info.push_back({"Зацепления сейчас", format("%d / %d", 3 * (count_ - 1) - broken_, 3 * (count_ - 1))});
    out.info.push_back({"Первый разрыв, кадр", format("%d", firstBrokenFrame_)});
    out.info.push_back({"Макс. расстояние соседей", format("%.2f мм", 1000 * maxSeparation_)});
}

std::vector<SceneParam> TorusChainsScene::params() const {
    return {{"Звеньев в каждой цепи", float(count_), 8, 100, 1, 0, false, "Три цепи; максимум 300 невыпуклых торов."},
            {"Масса груза / звена", loadRatio_, 1, 50, 1, 0, false, "Жёлтый тор тяжелее остальных, но имеет ту же форму."},
            {"Поперечный толчок груза", kick_, 0, 5, 0.1f, 1, false, "Начальная скорость жёлтых торов вдоль z, м/с."}};
}

void TorusChainsScene::setParam(int index, float value) {
    if (!std::isfinite(value)) return;
    if (index == 0) count_ = int(std::clamp(value, 8.0f, 100.0f));
    if (index == 1) loadRatio_ = std::clamp(value, 1.0f, 50.0f);
    if (index == 2) kick_ = std::clamp(value, 0.0f, 5.0f);
}

void addTorusChainSamples(std::vector<SampleEntry>& out) {
    out.push_back({Preset::TorusChains, "Твёрдые тела", "Стресс-тест: 120 торов — груз, падение, подвес",
                   [] { return std::unique_ptr<Scene>(new TorusChainsScene); }});
}

} // namespace rf
