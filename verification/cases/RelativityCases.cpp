// Verification cases of general relativity: the symplectic geodesic integrator of
// src/relativity/SymplecticGeodesic.h on a Kerr orbit.
//
// gr-symplectic-drift - code verification of the property a symplectic method is chosen for: its
//   energy error only wobbles, it does not drift. The value is |least-squares trend of H - H0 over
//   the run| / (largest |H - H0|): 0 for a pure wobble, 1 for a pure drift (RK4 gives 1.00). The
//   reference is 0 (Hairer, Lubich & Wanner 2006, ch. IX: a symplectic step is the exact flow of a
//   nearby Hamiltonian, so H stays within O(h^p) of its start for exponentially long times). The
//   refinement study: the largest |H - H0| over 20 orbits at h = 2, 1, 1/2 M, of order 4 for Tao-4.
// The orbit: a unit mass around a Kerr hole a = 0.9 M, turning points 8 M and 15 M, Q = 4 M^2.
#include "../Benchmark.h"
#include "../Convergence.h"

#include "core/Format.h"
#include "relativity/KerrHamiltonian.h"
#include "relativity/SymplecticGeodesic.h"

#include <algorithm>
#include <cmath>

namespace rf::verify {

namespace {

constexpr double kSpin = 0.9, kCarterQ = 4.0, kOmega = 0.1;

// The particle at the apoapsis (15 M) on the equator of the orbit turning at 8 M and 15 M.
GeodesicState kerrOrbit() {
    const KerrHamiltonian H(1.0, kSpin);
    double E, L;
    H.boundOrbit(8.0, 15.0, kCarterQ, E, L);
    GeodesicState s;
    s.x[1] = 15.0;
    s.x[2] = 0.5 * 3.14159265358979323846;
    s.p[0] = -E;
    s.p[2] = std::sqrt(kCarterQ);
    s.p[3] = L;
    return s;
}

// What one run of a method measured: the largest excursion of H and the trend/excursion ratio.
struct EnergyError {
    double largest = 0;
    double trendRatio = 0;
};

// `orbits` radial periods with Tao-4 (tao = true) or RK4, step h: H - H0 at every step.
EnergyError energyError(bool tao, double h, double orbits) {
    const KerrHamiltonian H(1.0, kSpin);
    const SymplecticGeodesic sg(H);
    const GeodesicState start = kerrOrbit();
    const long long steps = std::llround(orbits * sg.radialPeriod(start) / h);
    GeodesicState s = start;
    TaoState z = SymplecticGeodesic::taoStart(start);
    const double h0 = H.value(start.x, start.p);
    EnergyError e;
    double st = 0, sy = 0, stt = 0, sty = 0;
    for (long long k = 1; k <= steps; ++k) {
        if (tao) {
            sg.stepTao(z, h, 4, kOmega);
            s = z.s;
        } else {
            sg.stepRK4(s, h);
        }
        const double t = double(k) * h, dH = H.value(s.x, s.p) - h0;
        e.largest = std::max(e.largest, std::fabs(dH));
        st += t, sy += dH, stt += t * t, sty += t * dH;
    }
    const double n = double(steps), slope = (n * sty - st * sy) / (n * stt - st * st);
    e.trendRatio = std::fabs(slope * n * h) / std::max(e.largest, 1e-300);
    return e;
}

Result runSymplecticDrift(const RunOptions& o) {
    Result r;
    r.hLabel = "h, M";
    r.theoreticalOrder = 4;
    // 1. The drift over a long run: Tao-4 against RK4 at the same step.
    const double orbits = o.full ? 3000 : 300;
    const EnergyError tao = energyError(true, 1.0, orbits), rk4 = energyError(false, 1.0, orbits);
    r.value = tao.trendRatio;
    // 2. The order of the excursion over 20 orbits: h = 2, 1, 1/2 M.
    for (double h : {2.0, 1.0, 0.5}) {
        const double largest = energyError(true, h, 20).largest;
        r.convergence.push_back({h, largest, largest});
    }
    const OrderFit fit = fitOrder(r.convergence);
    r.observedOrder = fit.order;
    r.detail = format("Тао-4 (h = M, ω = 0,1): max|ΔH| %.2e, тренд/размах %.3f за %.0f витков; RK4 с тем же шагом: max|ΔH| %.2e, "
                      "тренд/размах %.2f — чистый уход",
                      tao.largest, tao.trendRatio, orbits, rk4.largest, rk4.trendRatio);
    return r;
}

} // namespace

void addRelativityCases(std::vector<Case>& cases) {
    cases.push_back({"gr-symplectic-drift", "Геодезические Керра: энергия без векового дрейфа (Тао-4)", "code-verification",
                     "symplectic geodesics: energy bounded over thousands of orbits, RK4 drifts",
                     {"Hairer, Lubich, Wanner 2006: симплектический шаг — точный поток близкого гамильтониана, ошибка энергии ограничена",
                      "https://doi.org/10.1007/3-540-30666-8", 0, 0, "analytic",
                      "значение — |тренд H за прогон| / наибольшее |H − H₀|: 0 — колебание, 1 — чистый уход (у RK4 1,00)"},
                     runSymplecticDrift, {0.2, 0.5, false}, false});
}

} // namespace rf::verify
