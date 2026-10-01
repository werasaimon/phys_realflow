// Piecewise-convex square helices. Pitch is deliberately coarse and clearance explicit.
// Ideal screw kinematics y-y0 = pitch*(theta-theta0)/(2*pi) is a test oracle, not a constraint.
#include "samples/NutBoltScene.h"
#include "rigid/ConvexDecomposition.h"
#include "core/Format.h"
#include "core/Probe.h"

#include <algorithm>
#include <cmath>

namespace rf {
namespace {

constexpr float rootRadius = 0.030f, crestRadius = 0.038f, boreRadius = 0.0395f;
constexpr float toothWidth = 0.0036f, nutHeight = 0.032f;
const Vector3 boltOrigin(0, 0.14f, 0);

Vector3 radial(float radius, float angle, float y) {
    return {radius * rf::cos(angle), y, -radius * rf::sin(angle)};
}

TriMesh ridge(int index, float inner, float outer, float phase) {
    std::vector<Vector3> points;
    for (int end : {0, 1}) {
        const float turn = float(index + end) / NutBoltScene::segments, angle = 2 * kPi * turn;
        for (float radius : {inner, outer}) for (float side : {-0.5f, 0.5f})
            points.push_back(radial(radius, angle, NutBoltScene::pitch * turn + phase + side * toothWidth));
    }
    return buildConvexHull(points);
}

void addRidges(std::vector<TriMesh>& parts, bool internal) {
    const float phase = internal ? NutBoltScene::pitch * 0.5f : 0;
    const float limit = internal ? nutHeight * 0.5f : 0.060f;
    const float inner = internal ? rootRadius + 0.0008f : rootRadius;
    const float outer = internal ? boreRadius : crestRadius;
    for (int k = -6 * NutBoltScene::segments; k < 6 * NutBoltScene::segments; ++k) {
        const float lo = NutBoltScene::pitch * float(k) / NutBoltScene::segments + phase - 0.5f * toothWidth;
        const float hi = lo + NutBoltScene::pitch / NutBoltScene::segments + toothWidth;
        if (lo >= -limit && hi <= limit) parts.push_back(ridge(k, inner, outer, phase));
    }
}

float hexRadius(float a) {
    const float sector = kPi / 3, middle = (std::floor(a / sector) + 0.5f) * sector;
    return 0.055f * rf::cos(sector * 0.5f) / rf::cos(a - middle);
}

TriMesh nutWall(int index) {
    std::vector<Vector3> points;
    for (int end : {0, 1}) {
        const float angle = 2 * kPi * float(index + end) / NutBoltScene::segments;
        for (float radius : {boreRadius, hexRadius(angle)}) for (float y : {-nutHeight * 0.5f, nutHeight * 0.5f})
            points.push_back(radial(radius, angle, y));
    }
    return buildConvexHull(points);
}

} // namespace

std::shared_ptr<const CompoundShape> NutBoltScene::boltShape() {
    static const auto shape = [] {
        std::vector<TriMesh> parts;
        for (const auto& spec : std::vector<Vector3>{{rootRadius, 0.16f, 0}, {0.050f, 0.018f, -0.080f}}) {
            auto mesh = primitives::cylinder(spec.x, spec.y, spec.z == 0 ? segments : 6);
            mesh.transform(Quaternion::fromAxisAngle({1, 0, 0}, -0.5f * kPi).toMatrix3x3(), Vector3(1), {0, spec.z, 0});
            parts.push_back(std::move(mesh));
        }
        addRidges(parts, false);
        return std::make_shared<const CompoundShape>(parts, primitives::merge(parts));
    }();
    return shape;
}

std::shared_ptr<const CompoundShape> NutBoltScene::nutShape() {
    static const auto shape = [] {
        std::vector<TriMesh> parts;
        for (int k = 0; k < segments; ++k) parts.push_back(nutWall(k));
        addRidges(parts, true);
        return std::make_shared<const CompoundShape>(parts, primitives::merge(parts));
    }();
    return shape;
}

Vector3 NutBoltScene::modelOrigin(const RigidBody& body) {
    const auto& shape = static_cast<const CompoundShape&>(*body.shape);
    return body.pos - body.rotation() * shape.principalRotation().transposed() * shape.centerOfMass();
}

float NutBoltScene::modelAngle(const RigidBody& body) {
    const auto& shape = static_cast<const CompoundShape&>(*body.shape);
    const Vector3 x = body.rotation() * shape.principalRotation().transposed() * Vector3(1, 0, 0);
    return rf::atan2(-x.z, x.x);
}

void NutBoltScene::configure(Simulation& sim) {
    sim.rigid.params.substeps = 40;
    sim.rigid.params.iterations = 16;
    sim.rigid.params.slop = 0.0001f;
    sim.rigid.params.sleeping = false;
    sim.rigid.params.rollingResistance = 0;
    sim.rigid.params.rotationalLock = false;
    sim.rigid.params.shockPropagation = false;
}

void NutBoltScene::build(Simulation& sim) {
    sim.useRigidArena(AABB({-0.14f, 0, -0.14f}, {0.14f, 0.30f, 0.14f}));
    for (int i = 0; i < 2; ++i) {
        const auto shape = i == boltBody ? boltShape() : nutShape();
        const Vector3 origin = boltOrigin + Vector3(0, i == nutBody ? pitch : 0, 0);
        const int body = sim.rigid.addCompound(shape, origin + shape->centerOfMass(), Quaternion(),
                                               i == boltBody ? 0 : 7800, i == boltBody ? Vector3(0.6f) : Vector3(1, 0.7f, 0.15f));
        auto& b = sim.rigid.bodies()[size_t(body)];
        b.friction = 0.08f; b.staticFriction = 0.10f; b.restitution = 0;
    }
    driveWork_ = 0; commandedSpeed_ = 0;
}

void NutBoltScene::drive(Simulation& sim) {
    auto& b = sim.rigid.bodies()[nutBody];
    const float time = float(std::fmod(sim.time(), 12.0));
    const int phase = int(time / 3);
    const float target = phase == 1 ? 0.8f * (time - 3) : phase == 2 ? 0.8f * (9 - time) : 0;
    const float error = rf::atan2(rf::sin(target - modelAngle(b)), rf::cos(target - modelAngle(b)));
    const float feedforward = phase == 1 ? 0.8f : phase == 2 ? -0.8f : 0;
    commandedSpeed_ = clampv(feedforward + 8 * error, -1.4f, 1.4f);
    const Vector3 axis(0, 1, 0);
    const float k = dot(axis, b.applyInvInertiaWorld(axis));
    // Angle feedback compensates load-induced lag; only a bounded angular impulse is applied.
    // tau_max=0.12 Nm; neither the angle nor the axial travel is overwritten.
    const float impulse = clampv((commandedSpeed_ - dot(axis, b.angVel)) / k,
                                 -0.12f * sim.frameDt, 0.12f * sim.frameDt);
    driveWork_ += double(impulse) * dot(axis, b.angVel) + 0.5 * double(impulse) * impulse * k;
    sim.rigid.applyExternalWrench(nutBody, Vector3(0), axis * impulse);
}

void NutBoltScene::afterStep(Simulation& sim) {
    drive(sim);
    Probe::set("nut-bolt/height m", modelOrigin(sim.rigid.bodies()[nutBody]).y);
    Probe::set("nut-bolt/angle rad", modelAngle(sim.rigid.bodies()[nutBody]));
    Probe::set("nut-bolt/drive work J", driveWork_);
}

void NutBoltScene::describe(const Simulation& sim, RenderSnapshot& out) const {
    out.info.push_back({"Резьба", "Грубая квадратная, шаг 12 мм; полигональная геометрия"});
    out.info.push_back({"Опора", "Болт закреплён; гайка свободна, винтового шарнира нет"});
    out.info.push_back({"Привод", "3 с удержание, 3 с подъём, 3 с обратный ход, 3 с удержание"});
    out.info.push_back({"Высота / угол", format("%.4f м / %.3f рад", modelOrigin(sim.rigid.bodies()[nutBody]).y,
                                             modelAngle(sim.rigid.bodies()[nutBody]))});
    out.info.push_back({"Работа привода", format("%.6f Дж", driveWork_)});
}

void addNutBoltSample(std::vector<SampleEntry>& out) {
    out.push_back({Preset::NutBolt, "Твёрдые тела", "Гайка и болт: контакт винтовой резьбы",
                   [] { return std::unique_ptr<Scene>(new NutBoltScene); }});
}

} // namespace rf
