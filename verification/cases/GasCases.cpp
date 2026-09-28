// Verification and validation cases of the gas solver (GasSolver, incompressible Navier-Stokes on
// the MAC grid):
//   manufactured heat   - the heat equation with a manufactured solution: order 2 expected
//   Taylor-Green vortex - an exact decaying flow: the amplitude falls as exp(-2 nu k^2 t)
//   Poiseuille, Couette-Poiseuille - the exact laminar profiles between two plates
//   heat spot, advected blob - the refinement studies of tests/GasTests.cpp (testGridConvergence)
//   cylinder at Re 100  - the Strouhal number against Williamson's experiments
#include "../Benchmark.h"
#include "../Convergence.h"
#include "../Mms.h"

#include "core/Format.h"
#include "core/Mesh.h"
#include "gas/GasSolver.h"
#include "spatial/BVH.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace rf::verify {

static const double kPiD = 3.14159265358979323846;

// Order from the errors of a study; u_num of the order = the standard error of the fit.
static void orderFromStudy(Result& r) {
    const OrderFit fit = fitOrder(r.convergence);
    r.value = r.observedOrder = fit.order;
    r.numericalUncertainty = std::isfinite(fit.stdError) ? fit.stdError : 0;
}

// ---------------------------------------------------------------------------------------------
static Result runManufacturedHeatCase(const RunOptions& o) {
    Result r;
    r.hLabel = "dx, м";
    r.theoreticalOrder = 2;
    const ManufacturedHeat m;
    const int levels = 3; // 16/32/64; a 128 grid costs 16 times the 64 one (dt ~ dx^2)
    for (int l = 0; l < levels; ++l) {
        const MmsLevel level = runManufacturedHeat(m, 16 << l);
        r.convergence.push_back({level.dx, level.errorRms, level.errorRms});
    }
    orderFromStudy(r);
    r.detail = format("RMS ошибка T через одно время диффузии: %.2e K (16 ячеек) … %.2e K (%d ячеек), амплитуда %.0f K, α(T) ~ T^1.75",
                      r.convergence.front().error, r.convergence.back().error, 16 << (levels - 1), m.amplitude);
    return r;
}

// ---------------------------------------------------------------------------------------------
// Taylor-Green vortex in the unit square, u = U sin(pi x) cos(pi y), v = -U cos(pi x) sin(pi y):
// the box walls of the solver are free-slip (no normal flow, no tangential stress), and this flow
// has exactly that on x, y = 0 and 1, so the half-period box holds it without error. The
// nonlinear term is a pure gradient, absorbed by the pressure; the amplitude decays as
// exp(-2 nu pi^2 t) (Taylor & Green 1937). Thin in z (4 cells), w = 0.
static double taylorGreenDecay(int cells, double nu, double tEnd) {
    GasSolver g;
    const float dx = 1.0f / float(cells), U = 1.0f;
    g.params.domainSize = {1, 1, 4 * dx};
    g.params.resolutionX = cells;
    g.params.inflowSpeed = 0;
    g.params.smokeRake = false;
    g.params.kinematicViscosity = float(nu);
    g.params.wallFriction = false;
    for (auto& b : g.params.bc) b = BoundaryType::Wall;
    g.reset({0, 0, 0}, nullptr);
    const float k = float(kPiD);
    g.setVelocity([&](const Vector3& x) {
        return Vector3(U * std::sin(k * x.x) * std::cos(k * x.y), -U * std::cos(k * x.x) * std::sin(k * x.y), 0);
    });
    // Amplitude: the projection of the cell velocity on the mode (the same measure at t = 0 and
    // at the end, so the averaging of the faces into cells cancels in the ratio).
    auto amplitude = [&] {
        double num = 0, den = 0;
        for (int kk = 0; kk < g.nz(); ++kk)
            for (int j = 0; j < g.ny(); ++j)
                for (int i = 0; i < g.nx(); ++i) {
                    const double x = (i + 0.5) * dx, y = (j + 0.5) * dx;
                    const double phi = std::sin(kPiD * x) * std::cos(kPiD * y);
                    num += g.cellVelocity(i, j, kk).x * phi;
                    den += phi * phi;
                }
        return num / den;
    };
    const double a0 = amplitude();
    double t = 0;
    const float dt = 0.5f * dx / U; // Courant number 0.5 on every grid
    while (t < tEnd - 1e-9) t += g.step(float(std::min<double>(dt, tEnd - t)));
    return amplitude() / a0;
}

static Result runTaylorGreen(const RunOptions& o) {
    Result r;
    r.hLabel = "dx, м";
    // The order in time of the scheme that runs: 2 for advection-reflection (the implicit midpoint
    // rule), 1 for the classic projection step (semi-Lagrangian + implicit Euler viscosity, dt ~ dx).
    const bool reflection = GasParams().advectionReflection;
    r.theoreticalOrder = reflection ? 2 : 1;
    const double nu = 0.01, tEnd = 1.0, exact = std::exp(-2 * nu * kPiD * kPiD * tEnd);
    const int levels = o.full ? 4 : 3;
    for (int l = 0; l < levels; ++l) {
        const int cells = 16 << l;
        const double ratio = taylorGreenDecay(cells, nu, tEnd);
        r.convergence.push_back({1.0 / cells, ratio, std::fabs(ratio - exact)});
    }
    orderFromStudy(r);
    r.detail = format("схема: %s; A(1 с)/A(0): %.4f (16 ячеек) … %.4f (%d ячеек), точное exp(-2νπ²t) = %.4f, ν = 0.01",
                      reflection ? "перенос с отражением (advectionReflection)" : "проекция в конце шага", r.convergence.front().value,
                      r.convergence.back().value, 16 << (levels - 1), exact);
    return r;
}

// ---------------------------------------------------------------------------------------------
// Laminar flow between two plates H = 1 m apart, 4 m long, at Re = U H / nu = 2. The plates are
// wall cells (the solver's vessel: every cell outside is static solid, 2 cells thick), the z
// walls are free-slip, so the flow is plane. The wall-function friction is off: at Re = 2 the
// flow is resolved, the no-slip condition alone holds the profile. Uniform inflow U (a flux U H):
//   Poiseuille:          still plates, u = 6 U eta (1 - eta), eta = y / H;
//   Couette-Poiseuille:  the top plate (a moving solid) at U, u = U eta + 3 U eta (1 - eta)
//                        (White, "Viscous Fluid Flow", 3rd ed., sec. 3-2).
// Pure Couette (the lid at 2U, the flux of the linear profile) has no curvature: a second-order
// scheme holds a linear profile exactly, so its error is the solvers' tolerance on every grid and
// has no order to measure - it is reported beside the Couette-Poiseuille order.
// The profile is measured at x = 3 m (developed; entrance length ~0.6 H at Re 2), started from
// the exact profile, after 2 viscous times H^2 / nu. The coarsest grid runs again at Courant 0.1:
// a steady flow must not depend on the time step.
//
// All the profiles are one formula: plates at `bottom` and `bottom + width`, the top one moving at
// Ul, a flux U H through the channel: u = Ul eta + 6 (U H / width - Ul / 2) eta (1 - eta).
// The diagnostic `shifted` measures against the same formula with each plate half a cell farther
// out. Until 2026-09-27 the no-slip speed sat on the first velocity sample inside the solid, half
// a cell behind the wall: the error against the shifted profile converged and the true one stayed
// first order (0.97 Poiseuille, 0.95 Couette), and the steady profile moved with the time step
// (0.165 at Courant 0.5, 0.135 at 0.1: the viscosity had a fixed 20 Jacobi sweeps, and the wall
// slipped by dt |grad p| / rho after the projection). src/gas/Viscosity.cpp now puts the no-slip
// speed on the wall itself, solves to a tolerance and carries the old pressure through the solve.
static double channelProfile(double y, double bottom, double width, double U, double H, double Ul) {
    const double eta = (y - bottom) / width;
    return Ul * eta + 6 * (U * H / width - 0.5 * Ul) * eta * (1 - eta);
}

struct ChannelErrors {
    double error = 0;   // RMS / U against the exact profile
    double shifted = 0; // RMS / U against the profile with the plates half a cell out
};

// The errors of the profile at x = 3 m; the profile itself (u of every gas cell there) if asked.
static ChannelErrors measureProfile(const GasSolver& g, float y0, float H, float U, float Ul, std::vector<double>* profile) {
    const float dx = g.dx();
    const int i = int(std::lround(3.0f / dx)); // x = 3 m
    double sum2 = 0, sum2s = 0;
    int n = 0;
    for (int k = 0; k < g.nz(); ++k)
        for (int j = 0; j < g.ny(); ++j) {
            const float y = (j + 0.5f) * dx;
            if (y < y0 || y > y0 + H) continue;
            const double u = g.cellVelocity(i, j, k).x;
            sum2 += std::pow(u - channelProfile(y, y0, H, U, H, Ul), 2);
            sum2s += std::pow(u - channelProfile(y, y0 - 0.5 * dx, H + dx, U, H, Ul), 2);
            if (profile) profile->push_back(u);
            ++n;
        }
    return {std::sqrt(sum2 / std::max(n, 1)) / U, std::sqrt(sum2s / std::max(n, 1)) / U};
}

// One channel, `cellsAcross` cells between the plates, the top plate moving at lidSpeed x U (0:
// Poiseuille), stepped at the Courant number `courant` of the profile's fastest speed.
static ChannelErrors channelProfileError(int cellsAcross, float lidSpeed, float courant, std::vector<double>* profile = nullptr) {
    const float H = 1, U = 1, dx = H / float(cellsAcross), y0 = 2 * dx, Ul = lidSpeed * U;
    const bool lid = lidSpeed != 0;
    GasSolver g;
    g.params.domainSize = {4 * H, H + 4 * dx, 4 * dx};
    g.params.resolutionX = 4 * cellsAcross;
    g.params.inflowSpeed = U;
    g.params.kinematicViscosity = 0.5f;
    g.params.smokeRake = false;
    g.params.wallFriction = false;
    g.vessel = [&](const Vector3& x) { return lid ? x.y > y0 : (x.y > y0 && x.y < y0 + H); };
    g.reset({0, 0, 0}, nullptr);
    auto exact = [&](float y) { return float(channelProfile(y, y0, H, U, H, Ul)); };
    if (lid) {
        MovingSolid top;
        top.bounds = AABB({-1, y0 + H, -1}, {5 * H, 2 * H, 1});
        top.position = Vector3(2, y0 + H, 0);
        top.velocity = Vector3(Ul, 0, 0);
        top.inside = [&](const Vector3& x) { return x.y > y0 + H; };
        g.setMovingSolids({top});
    }
    g.setVelocity([&](const Vector3& x) { return Vector3(x.y > y0 && x.y < y0 + H ? exact(x.y) : 0.0f, 0, 0); });
    const double tEnd = 2 * H * H / g.params.kinematicViscosity;
    const float dt = courant * dx / std::max(1.5f * U, Ul);
    for (double t = 0; t < tEnd - 1e-9;) t += g.step(float(std::min<double>(dt, tEnd - t)));
    return measureProfile(g, y0, H, U, Ul, profile);
}

// RMS difference of two profiles of the same grid [U].
static double profileDifference(const std::vector<double>& a, const std::vector<double>& b) {
    double sum2 = 0;
    for (size_t m = 0; m < a.size() && m < b.size(); ++m) sum2 += (a[m] - b[m]) * (a[m] - b[m]);
    return std::sqrt(sum2 / double(std::max<size_t>(a.size(), 1)));
}

static Result runChannel(const RunOptions& o, float lidSpeed) {
    Result r;
    r.hLabel = "dx, м";
    r.theoreticalOrder = 2;
    const int levels = o.full ? 4 : 3;
    std::vector<ConvergencePoint> shifted;
    std::vector<double> atHalf, atTenth;
    for (int l = 0; l < levels; ++l) {
        const int cells = 8 << l;
        const ChannelErrors e = channelProfileError(cells, lidSpeed, 0.5f, l == 0 ? &atHalf : nullptr);
        r.convergence.push_back({1.0 / cells, e.error, e.error});
        shifted.push_back({1.0 / cells, e.shifted, e.shifted});
    }
    orderFromStudy(r);
    channelProfileError(8, lidSpeed, 0.1f, &atTenth);
    r.detail = format("RMS ошибка профиля / U: %.2e (8 ячеек на H) … %.2e (%d ячеек), Re = 2; ", r.convergence.front().error,
                      r.convergence.back().error, 8 << (levels - 1)) +
               format("шаг Куранта 0.5 против 0.1 (8 ячеек): профили различаются на %.1e U; ", profileDifference(atHalf, atTenth)) +
               format("против профиля со стенками на полъячейки дальше: %.2e … %.2e, порядок %.2f", shifted.front().error,
                      shifted.back().error, fitOrder(shifted).order);
    if (lidSpeed != 0) { // pure Couette: a linear profile, exact for a second-order scheme
        const ChannelErrors pure = channelProfileError(8, 2.0f, 0.5f);
        r.detail += format("; чистый Куэтт (крышка 2U, линейный профиль): ошибка %.1e U — точно, до допусков решателей", pure.error);
    }
    return r;
}

static Result runPoiseuille(const RunOptions& o) { return runChannel(o, 0.0f); }
static Result runCouette(const RunOptions& o) { return runChannel(o, 1.0f); }

// ---------------------------------------------------------------------------------------------
// The heat spot of testGridConvergence: T0(r) = 0.1 K (1 - r^2/R^2)^2 in still gas, alpha = 1e-3,
// 1 s; the peak temperature on grids 16/32/64 -> Richardson, GCI.
static double heatSpotPeak(int cells) {
    GasSolver g;
    g.params.domainSize = {1, 1, 1};
    g.params.resolutionX = cells;
    g.params.inflowSpeed = 0;
    g.params.smokeRake = false;
    for (auto& b : g.params.bc) b = BoundaryType::Wall;
    g.combustion.enabled = true;
    g.combustion.gravity = 0;
    g.combustion.radiativeCooling = 0;
    g.combustion.thermalDiffusivity = 1e-3f;
    g.reset({-0.5f, -0.5f, -0.5f}, nullptr);
    Disturbance spot;
    spot.radius = 0.15f;
    spot.velocityBlend = 0;
    spot.heat = 0.1f;
    g.applyDisturbance(spot);
    for (float t = 0; t < 1.0f - 1e-5f;) t += g.step(std::min(0.004f, 1.0f - t));
    double peak = 0;
    for (float T : g.temperature().d) peak = std::max(peak, double(T));
    return peak;
}

// The exact peak: the centre of the spot after 1 s, the 3D heat kernel against T0 (radial).
static double heatSpotExactPeak() {
    const double R = 0.15, a = 1e-3, t = 1.0, s = 4 * a * t;
    const int n = 2000;
    double sum = 0;
    for (int m = 0; m <= n; ++m) { // T(0) = int T0(r') 4 pi r'^2 G(r') dr', G = exp(-r'^2/s) / (pi s)^1.5
        const double rp = R * m / n, w = (m == 0 || m == n) ? 1 : (m % 2 ? 4 : 2);
        const double T0 = 0.1 * std::pow(1 - rp * rp / (R * R), 2.0);
        sum += w * T0 * 4 * kPiD * rp * rp * std::exp(-rp * rp / s) / std::pow(kPiD * s, 1.5);
    }
    return sum * (R / n) / 3;
}

static Result runHeatSpot(const RunOptions&) {
    Result r;
    r.unit = "K";
    r.hLabel = "dx, м";
    r.theoreticalOrder = 2;
    const double exact = heatSpotExactPeak();
    for (int cells : {16, 32, 64}) {
        const double peak = heatSpotPeak(cells);
        r.convergence.push_back({1.0 / cells, peak, std::fabs(peak - exact)});
    }
    const RichardsonEstimate re = richardson(r.convergence[2].value, r.convergence[1].value, r.convergence[0].value, 2.0);
    r.value = r.convergence[2].value;
    r.numericalUncertainty = numericalUncertaintyFromGci(re.gciFine);
    r.observedOrder = re.order;
    r.detail = format("пик %.5f / %.5f / %.5f K (16/32/64), Ричардсон %.5f K, точный %.5f K, p = %.2f, GCI %.2f%%",
                      r.convergence[0].value, r.convergence[1].value, r.convergence[2].value, re.extrapolated, exact, re.order,
                      100 * re.gciFine / r.value);
    return r;
}

// The smoke blob of testGridConvergence: exp(-|x - c|^2 / 2 s^2), s = 0.12 m, carried 0.5 m by a
// uniform flow at Courant 0.7 (MacCormack): the RMS error against the shifted blob.
static double advectedBlobError(int cells) {
    const float U = 1.0f, tEnd = 0.5f, sigma = 0.12f;
    const Vector3 c0(0.5f, 0.5f, 0.5f), c1 = c0 + Vector3(U * tEnd, 0, 0);
    GasSolver g;
    g.params.domainSize = {2, 1, 1};
    g.params.resolutionX = cells;
    g.params.inflowSpeed = U;
    g.params.smokeRake = false;
    g.reset({0, 0, 0}, nullptr);
    g.setTracer([&](const Vector3& x) { return std::exp(-length2(x - c0) / (2 * sigma * sigma)); });
    const float dt = 0.7f * g.dx() / U;
    for (float t = 0; t < tEnd - 1e-5f;) t += g.step(std::min(dt, tEnd - t));
    double sum2 = 0;
    const Field3& S = g.smoke();
    for (int k = 0; k < g.nz(); ++k)
        for (int j = 0; j < g.ny(); ++j)
            for (int i = 0; i < g.nx(); ++i) {
                const Vector3 x = g.origin() + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * g.dx();
                const double e = S.at(i, j, k) - std::exp(-length2(x - c1) / (2 * sigma * sigma));
                sum2 += e * e;
            }
    return std::sqrt(sum2 / double(S.d.size()));
}

static Result runAdvection(const RunOptions&) {
    Result r;
    r.hLabel = "dx, м";
    r.theoreticalOrder = 2;
    for (int cells : {32, 64, 128}) {
        const double e = advectedBlobError(cells);
        r.convergence.push_back({2.0 / cells, e, e});
    }
    orderFromStudy(r);
    r.detail = format("RMS ошибка дыма: %.2e / %.2e / %.2e (32/64/128), шаги %.2f и %.2f",
                      r.convergence[0].error, r.convergence[1].error, r.convergence[2].error,
                      orderFromErrors(r.convergence[0].error, r.convergence[1].error, 2),
                      orderFromErrors(r.convergence[1].error, r.convergence[2].error, 2));
    return r;
}

// ---------------------------------------------------------------------------------------------
// The cylinder of testCylinderStrouhal: D = 0.5 m in a 4 x 4 m channel (12.5 % blockage), Re 100,
// 40 s; St = f D / U from the zero crossings of the lift over the last 20 s.
static double cylinderStrouhal(int resolution) {
    GasSolver g;
    g.params.domainSize = {4, 4, 0.5f};
    g.params.resolutionX = resolution;
    const float D = 0.5f, U = 1.0f;
    g.params.inflowSpeed = U;
    g.params.kinematicViscosity = U * D / 100.0f;
    g.params.smokeRake = false;
    const TriMesh cylinder = primitives::cylinder(0.5f * D, 0.6f);
    MeshBVH bvh;
    bvh.build(cylinder);
    g.reset({-1.2f, -2, -0.25f}, &bvh);
    Disturbance nudge;
    nudge.center = {0.6f, 0.15f, 0.0f};
    nudge.radius = 0.2f;
    nudge.velocity = {U, 0.3f * U, 0.0f};
    g.applyDisturbance(nudge);
    double t = 0, prev = 0, first = -1, last = -1;
    int crossings = 0;
    while (t < 40.0) {
        t += g.step(0.1f);
        const double cl = g.liftCoefficient();
        if (t > 20.0 && (prev < 0) != (cl < 0)) {
            if (first < 0) first = t;
            last = t;
            ++crossings;
        }
        prev = cl;
    }
    return crossings > 2 ? (crossings - 1) / 2.0 / (last - first) * D / U : kNaN;
}

static Result runStrouhal(const RunOptions&) {
    Result r;
    r.hLabel = "dx, м";
    const int res[3] = {48, 72, 108}; // r = 1.5
    for (int n : res) {
        const double st = cylinderStrouhal(n);
        r.convergence.push_back({4.0 / n, st, kNaN});
    }
    const RichardsonEstimate re = richardson(r.convergence[2].value, r.convergence[1].value, r.convergence[0].value, 1.5);
    r.value = r.convergence[2].value;
    r.numericalUncertainty = numericalUncertaintyFromGci(re.gciFine);
    r.observedOrder = re.order;
    r.detail = format("St = %.4f / %.4f / %.4f (сетки 48/72/108), Ричардсон %.4f, p = %.2f%s", r.convergence[0].value,
                      r.convergence[1].value, r.convergence[2].value, re.extrapolated, re.order, re.oscillating ? " (колебательная)" : "");
    return r;
}

// ---------------------------------------------------------------------------------------------
void addGasCases(std::vector<Case>& cases) {
    cases.push_back({"gas-mms-heat", "Теплопроводность: построенное решение (MMS)", "code-verification", "",
                     {"Roache 2002, MMS; ожидаемый порядок схемы 2", "https://doi.org/10.1115/1.1436090", 2, 0, "analytic",
                      "T = A (2 + cos πx cos πy (1 + sin ωt / 2)) с α(T) ~ T^1.75, источник — из уравнения (Mms.h)"},
                     runManufacturedHeatCase, {0.3, 0.6, false}, false});
    // The expected order depends on the scheme that runs (a change logged 2026-09-28): the classic
    // projection step is first order in time (dt ~ dx, so in dx); advection-reflection is the
    // implicit midpoint rule, second order (Narain, Zehnder & Thomaszewski 2019, sec. 3.1).
    const bool reflection = GasParams().advectionReflection;
    cases.push_back({"gas-taylor-green", "Вихрь Тейлора–Грина: затухание", "code-verification", "",
                     {reflection ? "Taylor & Green 1937: A = e^{−2νk²t}; порядок 2 (перенос с отражением — неявная средняя точка, "
                                   "Narain et al. 2019)"
                                 : "Taylor & Green 1937: A = e^{−2νk²t}; порядок 1 (проекция в конце шага + неявная вязкость, dt ~ dx)",
                      "https://doi.org/10.1098/rspa.1937.0036", reflection ? 2.0 : 1.0, 0, "analytic",
                      "полупериод в ящике со скольжением — точное решение; перенос с отражением (advectionReflection, выключен по умолчанию): порядок 2.268 при ожидаемом 2"},
                     runTaylorGreen, {0.3, 0.6, false}, false});
    cases.push_back({"gas-poiseuille", "Течение Пуазейля между пластинами", "code-verification", "",
                     {"Пуазейль: u = 6U η(1−η); порядок MAC-схемы 2", "", 2, 0, "analytic", "Re = 2, стенки — ячейки сосуда; перенос с отражением (выключен по умолчанию): порядок 1.986"},
                     runPoiseuille, {0.3, 0.6, false}, false});
    cases.push_back({"gas-couette", "Течение Куэтта–Пуазейля: подвижная пластина", "code-verification", "",
                     {"Куэтт–Пуазейль: u = U_w η + 6(U − U_w/2) η(1−η); порядок MAC-схемы 2 (White, Viscous Fluid Flow, §3-2)", "", 2, 0,
                      "analytic", "верхняя пластина — движущееся тело со скоростью U; чистый Куэтт (линейный) схема передаёт точно; перенос с отражением (выключен по умолчанию): порядок 1.988"},
                     runCouette, {0.3, 0.6, false}, false});
    cases.push_back({"gas-heat-spot", "Остывание пятна тепла: GCI пика", "solution-verification", "grid convergence",
                     {"Точное решение уравнения теплопроводности (ядро)", "", heatSpotExactPeak(), 0, "analytic", "пик в центре через 1 с"},
                     runHeatSpot, {0.01, 0.03, true}, false});
    cases.push_back({"gas-advection", "Перенос MacCormack: порядок", "solution-verification", "grid convergence",
                     {"Selle et al. 2008, MacCormack: порядок 2", "https://doi.org/10.1007/s10915-007-9166-4", 2, 0, "analytic", ""},
                     runAdvection, {0.3, 0.6, false}, true});
    // The reference matches our channel, not an unbounded flow (a change logged 2026-09-28 in
    // verification/unblinding-log.md). Behr, Hastreiter, Mittal & Tezduyar 1995 computed exactly our
    // set-up - a uniform inflow, slip (symmetry) lateral walls, Re = 100 - for walls at A = 9 ... 32
    // radii from the cylinder (their table 4, two formulations). Our walls are 4 D = 8 radii away
    // (blockage 12.5 %), just below their closest; St is linear in the blockage there, so their two
    // closest points are extended to it: 0.1735 (space-time) and 0.1761 (velocity-pressure-stress),
    // mean 0.175. u_D = 0.004: half their spread (0.0013), the extension (~0.002) and their own
    // discretisation (~0.002), rounded up. Williamson's 0.164 is the unbounded flow of experiment.
    cases.push_back({"gas-cylinder-strouhal", "Цилиндр, Re = 100, загромождение 12.5 %: число Струхаля", "validation",
                     "benchmark: cylinder vortex street",
                     {"Behr, Hastreiter, Mittal, Tezduyar 1995: равномерный приток, скользящие стенки, Re = 100 — продолжение их табл. 4 до "
                      "загромождения 12.5 %",
                      "https://doi.org/10.1016/0045-7825(94)00736-7", 0.175, 0.004, "code",
                      "их A = 9 и 12.5 радиуса: St 0.1711 / 0.1658 и 0.1739 / 0.1690, продолжено до A = 8; безграничный поток "
                      "(эксперимент Williamson 1996, doi 10.1146/annurev.fl.28.010196.002401): 0.164; перенос с отражением (выключен по умолчанию): 0.1745"},
                     runStrouhal, {0.05, 0.10, true}, true});
}

} // namespace rf::verify
