#pragma once
// Exact discrete action for a translational body under a constant acceleration. Its discrete
// Legendre transforms give the ballistic flow; see docs/book/14-variational-impact.md.
#include "rigid/RigidBody.h"

namespace rf {

class VariationalFlight {
public:
    static double action(const Vector3& q0, const Vector3& q1, const Vector3& acceleration, double mass, double h);
    static void advance(RigidBody& body, const Vector3& acceleration, float h);
};

} // namespace rf
