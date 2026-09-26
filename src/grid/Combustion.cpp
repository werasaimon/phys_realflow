#include "grid/Combustion.h"

#include "core/Parallel.h"

#include <algorithm>
#include <cmath>

namespace rf {

double Combustion::react(std::vector<float>& fuel, std::vector<float>& products, std::vector<float>& temperature,
                         std::vector<float>& smoke, std::vector<float>& expansionRate, const std::vector<uint8_t>& solid,
                         float dt) const {
    const int n = int(fuel.size());
    const float burntFraction = 1.0f - std::exp(-burnRate * dt); // exact for a first-order reaction
    // Radiative cooling dT/dt = -c (T/1000)^4 has the exact solution T / cbrt(1 + 3 c (T/1000)^3 dt / 1000).
    const float k = 3.0f * radiativeCooling * 1e-3f * dt;
    return parallelSum<double>(n, [&](int b, int e) {
        double burntSum = 0;
        for (int c = b; c < e; ++c) {
            expansionRate[c] = 0;
            if (solid[c]) continue;
            float& T = temperature[c];
            const float air = 1.0f - products[c]; // what can still burn here (oxygen)
            if (fuel[c] > 0 && air > 0 && T >= ignitionTemperature) {
                const float burnt = std::min(fuel[c], air) * burntFraction;
                fuel[c] -= burnt;
                products[c] += burnt;
                T += heatRelease * burnt;
                smoke[c] += sootYield * burnt;
                expansionRate[c] = expansion * burnt / dt;
                burntSum += burnt;
            }
            if (T > 0) {
                const float t = T * 1e-3f;
                T /= std::cbrt(1.0f + k * t * t * t);
            }
        }
        return burntSum;
    }, 4096);
}

} // namespace rf
