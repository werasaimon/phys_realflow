// Calls the compiled solver block: -fno-access-control is for this diagnostic only.
// W has eigenvalues 1.4, 0.6, 1; its inverse is a realizable rigid-body inertia.
#include "rigid/RigidWorld.h"
#include "core/Format.h"
#include <iostream>

using namespace rf;

RigidBody body(Vector3 omega) {
    RigidBody b;
    b.shape = std::make_shared<BoxShape>(Vector3(std::sqrt(41.0f / 14),
        std::sqrt(1.0f / 14), std::sqrt(29.0f / 14)));
    b.invInertiaLocal = {1.4f, 0.6f, 1};
    b.rot = Quaternion::fromAxisAngle({0, 0, 1}, kPi / 4);
    b.angVel = omega;
    b.updateInertia();
    return b;
}

double energy(const RigidWorld& world) {
    double result = 0;
    for (const auto& b : world.bodies())
        result += 0.5 * dot(b.angVel, b.invInertiaWorld.inverse() * b.angVel);
    return result;
}

Vector3 momentum(const RigidWorld& world) {
    Vector3 result;
    for (const auto& b : world.bodies()) result += b.invInertiaWorld.inverse() * b.angVel;
    return result;
}

void run(const char* name, Vector3 omega, bool pair, float limit) {
    RigidWorld world;
    world.bodies().push_back(body(omega));
    if (pair) world.bodies().push_back(body({0.03f, -0.2f, 0.1f}));
    RigidWorld::Manifold m;
    m.a = 0;
    m.b = pair ? 1 : -1;
    m.normal = {0, 1, 0};
    m.patchRadius = 1;
    Matrix3x3 W = world.bodies()[0].invInertiaWorld;
    if (pair) W += world.bodies()[1].invInertiaWorld;
    m.rollMass = W.inverse();
    m.massTwist = 1 / dot(m.normal, W * m.normal);
    const double before = energy(world);
    const Vector3 initialL = momentum(world);
    const float speed = omega.y - (pair ? world.bodies()[1].angVel.y : 0);
    const float expected = clampv(-m.massTwist * speed, -limit, limit);
    world.solveRotationalLock(m, limit);
    const float residual = world.bodies()[0].angVel.y - (pair ? world.bodies()[1].angVel.y : 0);
    std::cout << "{\"case\":\"" << name << "\",\"delta_energy_J\":" << numberText(energy(world) - before)
        << ",\"impulse_Nms\":" << numberText(m.jlock.y) << ",\"expected_Nms\":" << numberText(expected)
        << ",\"normal_speed_after\":" << numberText(residual)
        << ",\"delta_angular_momentum\":" << numberText(length(momentum(world) - initialL)) << "}\n";
}

int main() {
    run("zero normal speed", {0.1f, 0, 0}, false, 1);
    run("normal spin", {0.1f, 0.4f, 0}, false, 1);
    run("dynamic pair", {0.1f, 0.4f, 0}, true, 1);
    run("clamped pair", {0.1f, 0.4f, 0}, true, 0.01f);
}
