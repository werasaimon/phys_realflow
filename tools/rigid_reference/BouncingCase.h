#pragma once
// Shared SI parameters and observed state for the adapted SimBenchmark Bouncing scene.
// Author scenario: https://leggedrobotics.github.io/SimBenchmark/bouncing/
#include "math/Math.h"
#include <cstddef>

namespace rf::benchmark {
constexpr int side = 7, count = side * side;
constexpr float radius = 0.1f, mass = 10, height = 5, gravity = 9.81f, duration = 20;
constexpr double referenceEnergy = 24034.5, sphereInertia = 0.04;

struct State {
    rf::Vector3 position, velocity, angularVelocity;
};

struct Work {
    size_t accepted = 0, rejected = 0;
    double time = 0;
    size_t sdkImpactEvents = 0, sdkCcdQueries = 0;
    double sdkBroadMs = 0, sdkNarrowMs = 0, sdkSolveMs = 0, sdkCcdMs = 0;
};

inline rf::Vector3 startPosition(int i) {
    return {float(i / side) * 2 - 10, height, float(i % side) * 2 - 10};
}

} // namespace rf::benchmark
