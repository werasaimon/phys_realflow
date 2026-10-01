// Reproduce the new book examples with the SDK's C++17 math types, without Python.
// These are local algebraic checks, not a test of the complete contact solver or CCD.
#include "core/Format.h"
#include "rigid/RigidBody.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace rf;

static void closeTo(float actual, float expected, const char* what) {
    if (!std::isfinite(actual) || std::fabs(actual - expected) > 2e-6f)
        throw std::runtime_error(what);
}

static void poseCorrectionEnergy() {
    const Vector3 gravity{0, -9.81f, 0}, displacement{0, 0.001f, 0};
    const float potentialGain = -dot(gravity, displacement); // mass = 1 kg
    const auto inertia = Matrix3x3::diag({1, 2, 2});
    const auto rotation = Quaternion::fromAxisAngle({0, 0, 1}, kPi / 4).toMatrix3x3();
    const auto rotated = rotation * inertia * rotation.transposed();
    const Vector3 spin{1, 0, 0}; // Held fixed in world coordinates during this pose-only example.
    const float before = 0.5f * dot(spin, inertia * spin);
    const float after = 0.5f * dot(spin, rotated * spin);
    const auto angularChange = (rotated - inertia) * spin;
    closeTo(potentialGain, 0.00981f, "pose correction: potential energy");
    closeTo(before, 0.5f, "pose correction: initial kinetic energy");
    closeTo(after, 0.75f, "pose correction: rotated kinetic energy");
    closeTo(length(angularChange - Vector3{0.5f, -0.5f, 0}), 0, "pose correction: angular momentum");
    std::cout << "pose correction: delta_U_J=" << numberText(potentialGain)
              << ", T_before_J=" << numberText(before) << ", T_after_J=" << numberText(after) << '\n';
}

static void zeroTargetSplit() {
    RigidBody a, b; // Unit masses and identity inverse inertias.
    b.pos = {1, 0, 0};
    a.updateInertia();
    b.updateInertia();
    a.biasVel = {0, -1, 0};
    const Vector3 point{0.5f, 0, 0}, normal{0, 1, 0};
    const auto ra = point - a.pos, rb = point - b.pos;
    const auto ta = cross(ra, normal), tb = cross(rb, normal);
    const float inverseMass = a.invMass + b.invMass
        + dot(ta, a.applyInvInertiaWorld(ta)) + dot(tb, b.applyInvInertiaWorld(tb));
    const float impulse = 1 / inverseMass; // Closing speed -1 m/s, zero target.
    const auto pairImpulse = normal * impulse;
    a.biasVel += pairImpulse * a.invMass;
    b.biasVel -= pairImpulse * b.invMass;
    a.biasAngVel += a.applyInvInertiaWorld(cross(ra, pairImpulse));
    b.biasAngVel -= b.applyInvInertiaWorld(cross(rb, pairImpulse));
    const float after = dot(a.biasVel + cross(a.biasAngVel, ra)
        - b.biasVel - cross(b.biasAngVel, rb), normal);
    const auto angular = cross(a.pos, a.biasVel) + cross(b.pos, b.biasVel)
        + a.biasAngVel + b.biasAngVel;
    closeTo(inverseMass, 2.5f, "split: effective inverse mass");
    closeTo(impulse, 0.4f, "split: impulse");
    closeTo(after, 0, "split: contact velocity");
    closeTo(length(a.biasVel + b.biasVel - Vector3{0, -1, 0}), 0, "split: linear pseudo momentum");
    closeTo(length(angular), 0, "split: angular pseudo momentum");
    closeTo(length(a.vel) + length(b.vel) + length(a.angVel) + length(b.angVel), 0, "split: real velocities");
    std::cout << "zero-target split: impulse_Ns=" << numberText(impulse)
              << ", relative_pseudo_velocity_m_s=" << numberText(after) << '\n';
}

static void rayleighPower() {
    const float speed = 2, viscousCoefficient = 3, friction = 0.5f, load = 4;
    const float viscousR = 0.5f * viscousCoefficient * speed * speed;
    const float viscousPower = (-viscousCoefficient * speed) * speed;
    const float dryR = friction * load * std::fabs(speed);
    const float dryPower = (-friction * load) * speed;
    closeTo(viscousPower, -2 * viscousR, "Rayleigh: quadratic dissipation");
    closeTo(dryPower, -dryR, "Rayleigh: dry sliding dissipation");
    closeTo((-friction * load) * 0.0f, 0, "Rayleigh: static friction power");
    std::cout << "Rayleigh: viscous_power_W=" << numberText(viscousPower)
              << ", dry_power_W=" << numberText(dryPower) << '\n';
}

int main() {
    try {
        poseCorrectionEnergy();
        zeroTargetSplit();
        rayleighPower();
        std::cout << "PASS: 3 book examples; scope excludes complete simulation\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
