// Exact Ld = m|q1-q0|^2/(2h) + mh*a.(q0+q1)/2 - m|a|^2*h^3/24.
// p0=-D1 Ld, p1=D2 Ld: q1=q0+h*p0/m+a*h^2/2, p1=p0+m*a*h.
// Fetecau et al. (2003), sections 3.1-3.2, eqs. (47)-(56), exact-action specialization.
#include "rigid/VariationalFlight.h"

namespace rf {

double VariationalFlight::action(const Vector3& q0, const Vector3& q1, const Vector3& a, double mass, double h) {
    double value = 0;
    for (int k = 0; k < 3; ++k) {
        const double delta = double(q1[k]) - q0[k];
        value += mass * (delta * delta / (2 * h) + h * a[k] * (double(q0[k]) + q1[k]) / 2
                        - double(a[k]) * a[k] * h * h * h / 24);
    }
    return value;
}

void VariationalFlight::advance(RigidBody& b, const Vector3& a, float h) {
    b.prevPos = b.pos; b.prevRot = b.rot;
    for (int k = 0; k < 3; ++k) {
        b.pos[k] = float(double(b.pos[k]) + double(h) * b.vel[k] + 0.5 * double(h) * h * a[k]);
        b.vel[k] = float(double(b.vel[k]) + double(h) * a[k]);
    }
    b.rot = b.rot.integrated(b.angVel, h); // spherical inertia and zero torque in this mode
    b.updateInertia();
}

} // namespace rf
