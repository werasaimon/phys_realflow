#pragma once
// Helpers shared by the files of GasSolver (not part of the SDK's interface).

#include "gas/Field3.h"
#include "math/Math.h"

#include <chrono>

namespace rf {

// The three face fields sampled at a point in cell units: the velocity there.
inline Vector3 sampleVel(const Field3& u, const Field3& v, const Field3& w, const Vector3& gp) {
    return {u.sample(gp), v.sample(gp), w.sample(gp)};
}

// Milliseconds since `start`: the stages' times for lastPressureMs / lastSolidCouplingMs.
inline float millisecondsSince(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - start).count();
}

} // namespace rf
