// Dipole-dipole forces between magnets on rigid bodies (formulas and sources in Magnets.h).
#include "scene/Magnets.h"

#include "rigid/RigidWorld.h"

namespace rf {

namespace {
const float kMu0Over4Pi = 1e-7f; // mu0 / (4 pi) [T m / A]
} // namespace

Vector3 dipoleField(const Vector3& m, const Vector3& r) {
    const float d = length(r);
    if (d < 1e-9f) return Vector3(0.0f); // at the dipole itself the point model has no answer
    const Vector3 n = r / d;
    return (n * (3.0f * dot(m, n)) - m) * (kMu0Over4Pi / (d * d * d));
}

Vector3 dipoleForce(const Vector3& m1, const Vector3& m2, const Vector3& r) {
    const float d = length(r);
    if (d < 1e-9f) return Vector3(0.0f);
    const Vector3 n = r / d;
    const float a = dot(m1, n), b = dot(m2, n);
    const Vector3 bracket = m2 * a + m1 * b + n * dot(m1, m2) - n * (5.0f * a * b);
    return bracket * (3.0f * kMu0Over4Pi / (d * d * d * d));
}

float applyMagnetForces(RigidWorld& world, const std::vector<int>& bodies, const std::vector<Vector3>& momentsBody, float dt) {
    const int n = int(bodies.size());
    auto& all = world.bodies();
    // The moments in the world frame: each turns with its body.
    std::vector<Vector3> m(n), x(n);
    std::vector<float> radius(n);
    for (int i = 0; i < n; ++i) {
        const RigidBody& b = all[bodies[i]];
        m[i] = b.rotation() * momentsBody[i];
        x[i] = b.pos;
        radius[i] = b.boundingRadius();
    }
    float maxForce = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            Vector3 r = x[j] - x[i];
            // Softening: two touching magnets are not points - their centres cannot come closer
            // than about the sum of their sizes, and the 1/d^4 force of the point model would
            // blow up as d -> 0 (a numerical explosion, not physics). Keep d at least half the sum
            // of the bounding radii.
            const float dMin = 0.5f * (radius[i] + radius[j]), d = length(r);
            if (d < dMin) r = d > 1e-9f ? r * (dMin / d) : Vector3(dMin, 0, 0);
            const Vector3 F = dipoleForce(m[i], m[j], r);     // on j; on i it is -F
            const Vector3 tauJ = cross(m[j], dipoleField(m[i], r));
            const Vector3 tauI = cross(m[i], dipoleField(m[j], -r));
            world.applyExternalWrench(bodies[j], F * dt, tauJ * dt);
            world.applyExternalWrench(bodies[i], -F * dt, tauI * dt);
            maxForce = std::max(maxForce, length(F));
        }
    return maxForce;
}

} // namespace rf
