// Verification cases of magnetohydrodynamics.
//
// mhd-alfvenic-state - code verification with an exact NONLINEAR solution of incompressible MHD.
//   In Alfvén units, b = B / sqrt(mu0 rho), the equations in Elsässer's variables z± = v ± b
//   (Elsasser 1950, Phys. Rev. 79, 183) are
//       dz±/dt + (z∓ . grad) z± = -grad P + (nu + eta)/2 lap z± + (nu - eta)/2 lap z∓,
//   each field carried only by the other one. Where v = b everywhere, z- = 0 carries nothing, the
//   nonlinear terms cancel for any flow however tangled (Dobrowolny, Mangeney & Veltri 1980,
//   Phys. Rev. Lett. 45, 144), and with nu = eta every eigenmode of the Laplacian just decays:
//       A_k(t) = A_k(0) exp(-nu k^2 t).
//   The flow: two convection cells of different size in the unit square, psi = U1/pi sin(pi x)
//   sin(pi y) + U2/(2 pi) sin(2 pi x) sin(2 pi y), v = (dpsi/dy, -dpsi/dx), the field the same flow
//   (B = sqrt(mu0 rho) v), free-slip perfectly conducting walls that both cells fit. Without the
//   field the cells stir each other and the (2,2) cell even turns over (tests/PlasmaTests.cpp,
//   testAlfvenicState, runs that control); with it they must decay exactly as above.
//   The value is the observed order of the worst amplitude error at t = 1 s over three grids.
//   Its theoretical order is 1: the induction step adds the stabilising resistivity |u|^2 h of a
//   centred, forward-stepped advection of B (the Lax-Wendroff amount, MagneticField.cpp,
//   computeElectricField), which is first order in the step h ~ dx; the velocity side is second
//   order with advection-reflection.
//
// Why not the Orszag-Tang vortex (Orszag & Tang 1979, J. Fluid Mech. 90, 129; in its standard
// compressible form, Stone et al. 2008, ApJS 178, 137, sec. 4.2.3): it is compressible MHD with
// gamma = 5/3 and shocks, in a doubly periodic box. Our gas is incompressible (a pressure
// projection at constant density) and has walls, no periodic boundaries. docs/06, "Границы".
//
// The same setup as the test (tests/PlasmaTests.cpp): the test can run without this library.
#include "../Benchmark.h"
#include "../Convergence.h"

#include "core/Format.h"
#include "gas/GasSolver.h"

#include <algorithm>
#include <cmath>

namespace rf::verify {

namespace {

constexpr double kPiD = 3.14159265358979323846;
constexpr double kCellSpeed1 = 1.0, kCellSpeed2 = 0.5; // U1, U2 [m/s]
constexpr double kNu = 0.01, kSeconds = 1.0;            // nu = eta [m^2/s], the run [s]

double psi(double x, double y) {
    return kCellSpeed1 / kPiD * std::sin(kPiD * x) * std::sin(kPiD * y) +
           kCellSpeed2 / (2 * kPiD) * std::sin(2 * kPiD * x) * std::sin(2 * kPiD * y);
}

// The two cells as the flow and as the field: both the discrete curl of psi on the same faces,
// so v = b holds on the grid from the start. Unit square, 4 cells thick, eta = nu.
void startState(GasSolver& g, int cells) {
    const double dx = 1.0 / cells;
    g.params.domainSize = {1, 1, float(4 * dx)};
    g.params.resolutionX = cells;
    g.params.inflowSpeed = 0;
    g.params.fluidDensity = 1;
    g.params.kinematicViscosity = float(kNu);
    g.params.smokeRake = false;
    g.params.wallFriction = false;
    for (auto& b : g.params.bc) b = BoundaryType::Wall;
    g.magnetic.enabled = true;
    g.magnetic.conductivity = float(1.0 / (MagneticField::kMu0 * kNu));
    g.magnetic.numericalDissipation = 0;
    g.reset({0, 0, 0}, nullptr);
    g.setVelocity([dx](const Vector3& p) {
        return Vector3(float((psi(p.x, p.y + dx / 2) - psi(p.x, p.y - dx / 2)) / dx),
                       float(-(psi(p.x + dx / 2, p.y) - psi(p.x - dx / 2, p.y)) / dx), 0.0f);
    });
    const double alfven = std::sqrt(MagneticField::kMu0 * 1.0);
    g.magnetic.addFromPotential([alfven](const Vector3& p) { return Vector3(0, 0, float(alfven * psi(p.x, p.y))); });
}

// The amplitude of the cell (m, m): the projection of the cell velocities on its shape.
double mode(const GasSolver& g, int m) {
    const double dx = 1.0 / g.nx(), k = m * kPiD;
    double num = 0, den = 0;
    for (int kk = 0; kk < g.nz(); ++kk)
        for (int j = 0; j < g.ny(); ++j)
            for (int i = 0; i < g.nx(); ++i) {
                const double x = (i + 0.5) * dx, y = (j + 0.5) * dx;
                const double sx = std::sin(k * x) * std::cos(k * y), sy = -std::cos(k * x) * std::sin(k * y);
                const Vector3 v = g.cellVelocity(i, j, kk);
                num += v.x * sx + v.y * sy;
                den += sx * sx + sy * sy;
            }
    return num / den;
}

// |z-| / |z+| = |v - b| / |v + b| over the grid: 0 for an Alfvénic state.
double imbalance(const GasSolver& g) {
    const double alfven = std::sqrt(MagneticField::kMu0 * g.params.fluidDensity);
    double minus = 0, plus = 0;
    for (int kk = 0; kk < g.nz(); ++kk)
        for (int j = 0; j < g.ny(); ++j)
            for (int i = 0; i < g.nx(); ++i) {
                const Vector3 v = g.cellVelocity(i, j, kk), b = g.magnetic.cellField(i, j, kk) * float(1.0 / alfven);
                minus += length2(v - b);
                plus += length2(v + b);
            }
    return std::sqrt(minus / plus);
}

// One grid: the worst error of the two amplitude ratios at t = 1 s against exp(-nu k^2 t).
struct Level {
    double error = 0, ratio11 = 0, ratio22 = 0, imbalance = 0;
};

Level runLevel(int cells) {
    GasSolver g;
    startState(g, cells);
    const double a1 = mode(g, 1), a2 = mode(g, 2);
    const double dtMax = 0.25 / (cells * (kCellSpeed1 + kCellSpeed2)); // a quarter cell of the fastest speed
    while (g.time() < kSeconds - 1e-9) g.step(float(std::min(dtMax, kSeconds - g.time())));
    Level l;
    l.ratio11 = mode(g, 1) / a1;
    l.ratio22 = mode(g, 2) / a2;
    const double e1 = std::exp(-kNu * 2 * kPiD * kPiD * kSeconds), e2 = std::exp(-kNu * 8 * kPiD * kPiD * kSeconds);
    l.error = std::max(std::fabs(l.ratio11 - e1), std::fabs(l.ratio22 - e2));
    l.imbalance = imbalance(g);
    return l;
}

Result runAlfvenicState(const RunOptions& o) {
    Result r;
    r.hLabel = "dx, м";
    r.theoreticalOrder = 1;
    // 1. Three (four in --full) grids: the worst amplitude error of each.
    const int levels = o.full ? 4 : 3;
    std::vector<Level> runs;
    for (int l = 0; l < levels; ++l) {
        const int cells = 16 << l;
        runs.push_back(runLevel(cells));
        r.convergence.push_back({1.0 / cells, runs.back().error, runs.back().error});
    }
    // 2. The observed order and its standard error from the least-squares fit.
    const OrderFit fit = fitOrder(r.convergence);
    r.value = r.observedOrder = fit.order;
    r.numericalUncertainty = std::isfinite(fit.stdError) ? fit.stdError : 0;
    const Level& fine = runs.back();
    r.detail = format("A(1 с)/A(0) ячеек (1,1) и (2,2): %.4f / %.4f на сетке %d, точно %.4f / %.4f; худшая ошибка %.3f (16) … %.4f (%d); "
                      "|v − b| / |v + b|: %.1e (16) … %.1e (%d)",
                      fine.ratio11, fine.ratio22, 16 << (levels - 1), std::exp(-kNu * 2 * kPiD * kPiD), std::exp(-kNu * 8 * kPiD * kPiD),
                      runs.front().error, fine.error, 16 << (levels - 1), runs.front().imbalance, fine.imbalance, 16 << (levels - 1));
    return r;
}

} // namespace

void addPlasmaCases(std::vector<Case>& cases) {
    cases.push_back({"mhd-alfvenic-state", "Альфвеновское состояние v = b: точное нелинейное решение МГД", "code-verification",
                     "MHD: an Alfvénic state v = b is an exact nonlinear solution (Elsässer)",
                     {"Elsasser 1950; Dobrowolny, Mangeney & Veltri 1980: при z⁻ = 0 нелинейные члены исчезают, A_k = A_k(0) e^{−νk²t}; "
                      "порядок 1 (стабилизирующее сопротивление |u|²h в индукции)",
                      "https://doi.org/10.1103/PhysRevLett.45.144", 1, 0, "analytic",
                      "две ячейки (1,1) и (2,2) в квадрате со скользящими проводящими стенками, ν = η = 0,01, 1 с; без поля ячейки "
                      "перемешивают друг друга (контроль в тесте)"},
                     runAlfvenicState, {0.3, 0.6, false}, false});
}

} // namespace rf::verify
