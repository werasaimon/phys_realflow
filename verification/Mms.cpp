// The manufactured heat solution on the gas solver (see Mms.h).
#include "Mms.h"

#include "gas/GasSolver.h"

#include <algorithm>
#include <cmath>

namespace rf::verify {

static const double kPiD = 3.14159265358979323846;

double ManufacturedHeat::diffusionTime() const { return 1.0 / (alpha0 * kPiD * kPiD * (1.0 / (Lx * Lx) + 1.0 / (Ly * Ly))); }

double ManufacturedHeat::omega() const { return 2 * kPiD / diffusionTime(); }

double ManufacturedHeat::temperature(double x, double y, double t) const {
    const double phi = std::cos(kPiD * x / Lx) * std::cos(kPiD * y / Ly), g = 1 + 0.5 * std::sin(omega() * t);
    return amplitude * (2 + phi * g);
}

double ManufacturedHeat::source(double x, double y, double t) const {
    const double cx = std::cos(kPiD * x / Lx), sx = std::sin(kPiD * x / Lx);
    const double cy = std::cos(kPiD * y / Ly), sy = std::sin(kPiD * y / Ly);
    const double phi = cx * cy, w = omega(), g = 1 + 0.5 * std::sin(w * t);
    const double k2 = kPiD * kPiD * (1 / (Lx * Lx) + 1 / (Ly * Ly));
    const double grad2 = std::pow(kPiD / Lx * sx * cy, 2) + std::pow(kPiD / Ly * cx * sy, 2); // |grad phi|^2
    const double T = amplitude * (2 + phi * g), ratio = (T + ambient) / ambient;
    const double alpha = alpha0 * std::pow(ratio, 1.75), dAlpha = 1.75 * alpha0 / ambient * std::pow(ratio, 0.75);
    const double dTdt = amplitude * phi * 0.5 * w * std::cos(w * t);
    return dTdt + alpha * amplitude * g * k2 * phi - dAlpha * amplitude * amplitude * g * g * grad2;
}

// A still gas in a closed box with conduction only: fire on (the conduction lives there), no
// gravity, no radiation, no fuel.
static void setUpStillHeatBox(GasSolver& g, const ManufacturedHeat& m, int cells) {
    const float dx = float(m.Lx) / float(cells);
    g.params.domainSize = {float(m.Lx), float(m.Ly), 4 * dx};
    g.params.resolutionX = cells;
    g.params.inflowSpeed = 0;
    g.params.smokeRake = false;
    for (auto& b : g.params.bc) b = BoundaryType::Wall;
    g.combustion.enabled = true;
    g.combustion.gravity = 0;
    g.combustion.radiativeCooling = 0;
    g.combustion.thermalDiffusivity = float(m.alpha0);
    g.combustion.ambientTemperature = float(m.ambient);
    g.reset({0, 0, 0}, nullptr);
}

MmsLevel runManufacturedHeat(const ManufacturedHeat& m, int cells) {
    GasSolver g;
    setUpStillHeatBox(g, m, cells);
    g.setTemperature([&](const Vector3& x) { return float(m.temperature(x.x, x.y, 0)); });
    // The source runs on the solver's own clock, which it sums in double (in float, ~10^3 steps of
    // 10^-4 s drifted by ~10^-4 s and dropped the order on the finest grid from 2 to 1.5).
    g.heatSource = [&](const Vector3& x, double t) { return float(m.source(x.x, x.y, t)); };
    // One conduction substep per step at the hottest point (6 alpha_max dt / dx^2 = 0.8 < 0.9):
    // the splitting error of the source then falls as dt ~ dx^2, like the spatial error.
    const double alphaMax = m.alpha0 * std::pow((m.temperature(0, 0, 0.25 * 2 * kPiD / m.omega()) + m.ambient) / m.ambient, 1.75);
    const float dx = g.dx(), dt = float(0.8 * dx * dx / (6.0 * alphaMax));
    const double tEnd = m.diffusionTime();
    MmsLevel out;
    out.cells = cells;
    out.dx = dx;
    while (g.time() < tEnd - 1e-9) {
        g.step(float(std::min<double>(dt, tEnd - g.time()))); // taken whole: the CFL limit (no flow) is far larger
        ++out.steps;
    }
    const double t = g.time();
    double sum2 = 0;
    const Field3& T = g.temperature();
    for (int k = 0; k < g.nz(); ++k)
        for (int j = 0; j < g.ny(); ++j)
            for (int i = 0; i < g.nx(); ++i) {
                const Vector3 x = g.origin() + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * dx;
                const double e = double(T.at(i, j, k)) - m.temperature(x.x, x.y, t);
                sum2 += e * e;
                out.errorMax = std::max(out.errorMax, std::fabs(e));
            }
    out.errorRms = std::sqrt(sum2 / double(T.d.size()));
    return out;
}

} // namespace rf::verify
