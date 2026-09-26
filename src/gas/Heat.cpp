// Heat in the gas of GasSolver (fire): Fourier conduction with the conductivity of hot air, the
// thermal radiation of the flame's hot cells, and the fuel, heat and smoke given off by burning
// things in the flow.
#include "gas/GasSolver.h"

#include "core/Parallel.h"

namespace rf {

void GasSolver::conductHeat(float dt) {
    // Fourier's law in the gas: dT/dt = div(alpha grad T). The diffusivity of air grows with
    // temperature (~T^1.75); on a face it is the mean of its two cells. No flux into solids or
    // through the walls. Explicit, in as many sub-steps as stability needs (6 alpha h / dx^2 <= 1).
    const float Tair = combustion.ambientTemperature;
    const size_t n = temp_.d.size();
    std::vector<float> alpha(n);
    const float alphaMax = parallelMax<float>(int(n), 0.0f, [&](int b, int e) {
        float m = 0;
        for (int c = b; c < e; ++c) {
            alpha[c] = combustion.thermalDiffusivity * std::pow((temp_.d[c] + Tair) / Tair, 1.75f);
            m = std::max(m, alpha[c]);
        }
        return m;
    }, 4096);
    const int steps = std::max(1, int(std::ceil(6.0f * alphaMax * dt / (dx_ * dx_) / 0.9f)));
    const float h = dt / float(steps) / (dx_ * dx_);
    const long sx = 1, sy = nx_, sz = long(nx_) * ny_;
    std::vector<float> T;
    for (int s = 0; s < steps; ++s) {
        T = temp_.d;
        parallelFor(nz_, [&](int k) {
            for (int j = 0; j < ny_; ++j)
                for (int i = 0; i < nx_; ++i) {
                    const size_t c = cidx(i, j, k);
                    if (solid_[c]) continue;
                    float flux = 0;
                    auto face = [&](bool inside, long step) {
                        if (!inside || solid_[c + step]) return;
                        flux += 0.5f * (alpha[c] + alpha[c + step]) * (T[c + step] - T[c]);
                    };
                    face(i > 0, -sx); face(i < nx_ - 1, sx);
                    face(j > 0, -sy); face(j < ny_ - 1, sy);
                    face(k > 0, -sz); face(k < nz_ - 1, sz);
                    temp_.d[c] = T[c] + h * flux;
                }
        }, 1);
    }
}

void GasSolver::collectRadiators() {
    const float sigma = 5.670e-8f, Tair = combustion.ambientTemperature;
    const float kappaV = combustion.absorptionCoefficient * dx_ * dx_ * dx_;
    radiators_.clear();
    for (int k = 0; k < nz_; ++k)
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                const float T = temp_.at(i, j, k) + Tair;
                if (T < 700.0f) continue; // below that the flame's radiation is negligible
                const float power = 4.0f * kappaV * sigma * (T * T * T * T - Tair * Tair * Tair * Tair);
                radiators_.push_back({origin_ + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * dx_, power});
            }
}

float GasSolver::irradianceAt(const Vector3& world) const {
    const float minR2 = dx_ * dx_; // a point source is only fair from about a cell away
    float q = 0;
    for (const Radiator& r : radiators_) q += r.power / (4.0f * kPi * std::max(length2(r.position - world), minR2));
    return q;
}

void GasSolver::applyPendingEmissions() {
    const float cellVolume = dx_ * dx_ * dx_;
    const float heatCapacity = params.fluidDensity * combustion.specificHeat * cellVolume; // [J/K] of a cell
    for (const GasEmission& e : pendingEmissions_) {
        const Vector3 g = (e.position - origin_) / dx_;
        const int i = std::clamp(int(std::floor(g.x)), 0, nx_ - 1);
        const int j = std::clamp(int(std::floor(g.y)), 0, ny_ - 1);
        const int k = std::clamp(int(std::floor(g.z)), 0, nz_ - 1);
        if (solid(i, j, k)) continue;
        fuel_.at(i, j, k) += e.fuel / cellVolume;
        smoke_.at(i, j, k) += e.smoke / cellVolume;
        temp_.at(i, j, k) = std::max(0.0f, temp_.at(i, j, k) + e.heat / heatCapacity);
    }
    pendingEmissions_.clear();
}

} // namespace rf
