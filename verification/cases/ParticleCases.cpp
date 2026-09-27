// Verification and validation cases of the particle system (ParticleSystem: liquid, cloth):
//   elastic catenary - a cloth strip pinned at both ends sags under its weight; the exact shape of
//                      an elastic cable (Irvine 1981) gives the sag
//   dam break        - the front of a collapsing water column against Martin & Moyce 1952, with
//                      the experiment's uncertainty, a particle-size refinement study, the
//                      input uncertainty by Monte Carlo and several seeds; a blinded holdout case
#include "../Benchmark.h"
#include "../Convergence.h"
#include "../Uncertainty.h"

#include "core/Format.h"
#include "particles/ParticleSystem.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace rf::verify {

// ---------------------------------------------------------------------------------------------
// The elastic catenary (H. M. Irvine, "Cable Structures", MIT Press 1981, ch. 2): a cable of
// unstretched length L0, weight W in all, axial stiffness EA, hung between two supports at the
// same height a span l apart. With H the horizontal tension, the span and the mid-span sag are
//     l = H L0 / EA + 2 (H L0 / W) asinh(W / 2H),
//     d = W L0 / (8 EA) + (H L0 / W) (sqrt(1 + (W / 2H)^2) - 1).
// Our strip is pinned at its rest length (l = L0), so it sags only because it stretches; the
// inextensible catenary y = a cosh(x / a) is the limit EA -> infinity with slack. For a strip
// every quantity is per unit width: W = rho_A g L0, EA = the tensile stiffness [N/m].
struct CatenaryInput {
    double L0 = 1.0;          // unstretched length = span [m]
    double areaDensity = 0.3; // rho_A [kg/m^2]
    double stiffness = 300.0; // EA per unit width [N/m]: soft, so the sag is large enough to measure
    double g = 9.81;          // [m/s^2]
};

// Irvine's mid-span sag in two steps: 1) the horizontal tension H that makes the span equal L0
// (the span grows monotonically with H, so bisection on a log scale finds it); 2) H into the
// sag formula.
static double irvineSag(const CatenaryInput& in) {
    const double W = in.areaDensity * in.g * in.L0, EA = in.stiffness, L0 = in.L0;
    auto spanExcess = [&](double H) { return H * L0 / EA + 2 * (H * L0 / W) * std::asinh(W / (2 * H)) - L0; };
    double lo = 1e-6, hi = 1e6; // spanExcess < 0 at lo (too slack), > 0 at hi
    for (int k = 0; k < 200; ++k) {
        const double mid = std::sqrt(lo * hi);
        (spanExcess(mid) < 0 ? lo : hi) = mid;
    }
    const double H = std::sqrt(lo * hi);
    return W * L0 / (8 * EA) + (H * L0 / W) * (std::sqrt(1 + W * W / (4 * H * H)) - 1);
}

// How the strip is built, beyond Irvine's inputs:
//   bending  - negligible (compliance 1000 m/N): a cable carries no moment;
//   shear    - negligible (1 N/m) unless defaultShear: the diagonal links also resist a stretch
//              along the strip, and that is not the tensile stiffness being verified;
//   strength - the default (4000 N/m, far above the ~5 N/m carried): a cloth that can tear gets
//              tethers with 15 % slack, so its threads carry the load. A cloth that cannot tear
//              (tearable = false) keeps its tethers tight at the rest distance and cannot stretch
//              at all - measured too, as a warning about that setting.
struct StripSetup {
    float spacing = 0.01f; // particle radius r [m]; the strip is 4 spacings wide
    bool tearable = true;
    bool defaultShear = false;
};

// The sag of the middle column right now: minus the mean height of its particles.
static double midspanSag(const ParticleSystem& ps, const CatenaryInput& in, float r) {
    double height = 0;
    int n = 0;
    for (const Vector3& p : ps.positions())
        if (std::fabs(p.x - 0.5f * float(in.L0)) < 0.51f * r * ps.params.clothSpacing) { height += p.y; ++n; }
    return n > 0 ? -height / n : kNaN;
}

// The strip of `setup`, 1 m long, both short edges pinned, gravity still off.
static void hangStrip(ParticleSystem& ps, const CatenaryInput& in, const StripSetup& setup) {
    const float r = setup.spacing;
    ps.params.particleRadius = r;
    ps.reset(AABB({-1, -2, -1}, {2, 1, 1}));
    ClothMaterial m;
    m.areaDensity = float(in.areaDensity);
    m.tensileStiffness = float(in.stiffness);
    if (!setup.defaultShear) m.shearStiffness = 1.0f;
    m.bendCompliance = 1e3f;
    if (!setup.tearable) m.strengthWarp = m.strengthWeft = 0;
    ps.params.gravity = Vector3(0.0f);
    ps.addCloth({0, 0, 0}, {float(in.L0), 0, 0}, {0, 0, 4 * r * ps.params.clothSpacing}, m, 64 | 128, Vector3(1));
}

struct SagMeasurement {
    double sag = kNaN;   // equilibrium sag of the middle [m]
    double swing = kNaN; // max - min of the sag over the measuring window [m]: how still the strip is
};

// Irvine's sag is a static equilibrium, so the equilibrium is what is measured (since 2026-09-28,
// see verification/unblinding-log.md; before, the mean of a strip swinging freely after a sudden
// release). The recipe, at dt = 1/600 s:
//   1) 0..3 s   gravity ramps in smoothly, g(t) = g s(t / 3), s(x) = 3x^2 - 2x^3: no sudden load;
//   2) 3..6 s   gravity held; light damping dv = -c v dt, c = 2 1/s, on the free particles takes
//               the remaining swing out (e^-6 over the hold). A force proportional to velocity is
//               zero at rest, so it does not move the equilibrium - it only shortens the wait;
//   3) 5.5..6 s the mid-span sag is averaged, and its spread in the window is the residual swing.
static SagMeasurement equilibriumSag(const CatenaryInput& in, const StripSetup& setup) {
    ParticleSystem ps;
    hangStrip(ps, in, setup);
    const float dt = 1.0f / 600, rampTime = 3.0f, damping = 2.0f;
    const int steps = 3600, measureFrom = 3300;
    double sum = 0, lowest = 1e30, highest = -1e30;
    int samples = 0;
    for (int k = 0; k < steps; ++k) {
        const float x = std::min(1.0f, k * dt / rampTime); // 1) the ramp
        ps.params.gravity = Vector3(0, -float(in.g) * x * x * (3 - 2 * x), 0);
        if (k * dt >= rampTime)                            // 2) the damped hold; pinned edges untouched
            for (size_t i = 0; i < ps.size(); ++i)
                if (ps.invMasses()[i] > 0) ps.addVelocity(int(i), ps.velocities()[i] * (-damping * dt));
        ps.step(dt);
        if (k < measureFrom) continue;                     // 3) the measurement
        const double sag = midspanSag(ps, in, setup.spacing);
        sum += sag;
        lowest = std::min(lowest, sag);
        highest = std::max(highest, sag);
        ++samples;
    }
    SagMeasurement out;
    if (samples) { out.sag = sum / samples; out.swing = highest - lowest; }
    return out;
}

// The case: 1) the equilibrium sag on particle spacings 4, 2, 1 (and 0.5 in --full) cm;
// 2) Richardson and the GCI of the three finest; 3) two settings a user may meet, measured at
// 2 cm and reported, not judged: the default shear links, and a cloth that cannot tear.
static Result runCatenary(const RunOptions& o) {
    Result r;
    r.unit = "м";
    r.hLabel = "шаг частиц, м";
    const CatenaryInput in;
    const double exact = irvineSag(in);
    const float spacings[4] = {0.04f, 0.02f, 0.01f, 0.005f};
    const int levels = o.full ? 4 : 3;
    double swing = 0; // the largest residual swing over the levels
    for (int l = 0; l < levels; ++l) {
        const SagMeasurement m = equilibriumSag(in, StripSetup{spacings[l]});
        r.convergence.push_back({spacings[l], m.sag, std::fabs(m.sag - exact)});
        swing = std::max(swing, m.swing);
    }
    const size_t n = r.convergence.size();
    const RichardsonEstimate re = richardson(r.convergence[n - 1].value, r.convergence[n - 2].value, r.convergence[n - 3].value, 2.0);
    r.value = r.convergence[n - 1].value;
    r.numericalUncertainty = numericalUncertaintyFromGci(re.gciFine);
    r.observedOrder = fitOrder(r.convergence).order;
    const double withShear = equilibriumSag(in, StripSetup{0.02f, true, true}).sag;
    const double untearable = equilibriumSag(in, StripSetup{0.02f, false, false}).sag;
    r.detail = format("провис %.4f / %.4f / %.4f м (шаг %.3f … %.3f м), Ирвин %.4f м, Ричардсон %.4f м; ", r.convergence[n - 3].value,
                      r.convergence[n - 2].value, r.convergence[n - 1].value, r.convergence[n - 3].h, r.convergence[n - 1].h, exact,
                      re.extrapolated) +
               format("равновесие: остаточное качание ≤ %.1e м за последние 0.5 с; ", swing) +
               format("со сдвигом по умолчанию %.0f Н/м: %.4f м; ", double(ClothMaterial().shearStiffness), withShear) +
               format("неразрываемая (тетеры без слабины): %.4f м", untearable);
    return r;
}

// ---------------------------------------------------------------------------------------------
// Dam break (J. C. Martin & W. J. Moyce 1952, Phil. Trans. R. Soc. A 244, 312-324): a water
// column a wide and n^2 a = 2a high collapses along a floor; the front Z = z / a at the
// dimensionless time T = t sqrt(2 g / a) = 2. A holdout case of the blind analysis: its reference
// (their curve at T = 2) is sealed in verification/sealed/liquid-dam-break.ref, so it is not in
// this file, and the run never compares with the experiment while the case is blinded.
// (Honest note: the curve was seen before the sealing - it is in tests/ParticleTests.cpp - so the
// blind holds from the sealing on: tuning must not look at that test's deviation either.)
// SHA-256("value|uncertainty|salt") of the sealed file, printed by rf_verify --seal.
static const char* kDamBreakCommitment = "d3ddb03c0e8242c84924f041d01fa6860b11973d265dcad5060ebac194c26fe5";

// The bulk front: the last bin (one spacing wide) along the floor holding at least 5 particles of
// the floor layer (a few splashed particles ahead do not count), in units of a.
static double frontPosition(const ParticleSystem& s, float a, float r) {
    const float bin = 2 * r;
    std::vector<int> bins(size_t(6 * a / bin) + 1, 0);
    for (const Vector3& p : s.positions())
        if (p.y < 2.5f * r) ++bins[std::min(size_t((p.x + r) / bin), bins.size() - 1)];
    double front = 0;
    for (size_t b = 0; b < bins.size(); ++b)
        if (bins[b] >= 5) front = bin * double(b + 1) / a;
    return front;
}

// One run: particle radius r (the time step shrinks with r, so the speed limit of the particles,
// 2 r / dt, is the same on every level), column width a, floor friction. The seed perturbs the
// start: every particle gets a random velocity of at most 1 mm/s - nothing physical, but a
// splash amplifies it, and the spread over seeds is the run-to-run scatter of the front.
static double damBreakFront(float r, float a, float floorFriction, uint32_t seed) {
    const float g = 9.81f, depth = 0.1f;
    ParticleSystem s;
    s.params.particleRadius = r;
    s.params.gravity = {0, -g, 0};
    s.params.wallFriction = floorFriction;
    s.reset(AABB({0, 0, 0}, {6 * 0.2f, 3 * 0.2f, depth}));
    s.addBlock(AABB({0, 0, 0}, {a, 2 * a, depth}));
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> jitter(-0.577e-3f, 0.577e-3f); // |dv| <= 1 mm/s
    for (size_t i = 0; i < s.size(); ++i) s.addVelocity(int(i), Vector3(jitter(rng), jitter(rng), jitter(rng)));
    const float dt = (1.0f / 600) * (r / 0.005f), scale = std::sqrt(2 * g / a);
    for (int step = 0;; ++step) {
        if (step * dt * scale >= 2.0) return frontPosition(s, a, r);
        s.step(dt);
    }
}

static Result runDamBreak(const RunOptions& o) {
    Result r;
    r.hLabel = "радиус частицы, м";
    const float radii[3] = {0.01f, 0.00707f, 0.005f};
    for (float rad : radii) r.convergence.push_back({rad, damBreakFront(rad, 0.2f, 0.1f, o.seed), kNaN});
    const RichardsonEstimate re = richardson(r.convergence[2].value, r.convergence[1].value, r.convergence[0].value, std::sqrt(2.0));
    r.value = r.convergence[2].value;
    r.numericalUncertainty = numericalUncertaintyFromGci(re.gciFine);
    r.observedOrder = re.order;
    std::string uq = "UQ не считалась (--quick)";
    if (o.uqSamples > 1) { // on the coarse level: the sensitivities, not the resolution, matter here
        const std::vector<UncertainInput> inputs = {
            {"ширина столба a, м", 0.2, 0.2 * 0.015, UncertainInput::Distribution::Uniform},
            {"трение о дно", 0.1, 0.05, UncertainInput::Distribution::Uniform}};
        const UqSummary uqs = propagate(inputs, [&](const std::vector<double>& x) {
            return damBreakFront(radii[0], float(x[0]), float(x[1]), o.seed);
        }, o.uqSamples, o.seed);
        r.inputUncertainty = uqs.stddev;
        uq = format("UQ %d прогонов: σ %.3f, 95%% [%.2f, %.2f]", uqs.samples, uqs.stddev, uqs.lo95, uqs.hi95);
    }
    r.detail = format("Z(T=2) = %.3f / %.3f / %.3f (r = 10/7.1/5 мм); ", r.convergence[0].value, r.convergence[1].value,
                      r.convergence[2].value) + uq;
    return r;
}

// ---------------------------------------------------------------------------------------------
void addParticleCases(std::vector<Case>& cases) {
    cases.push_back({"cloth-elastic-catenary", "Ткань между двумя опорами: провис", "code-verification", "",
                     {"Irvine 1981, упругая цепная линия", "https://mitpress.mit.edu/9780262090230/cable-structures/", irvineSag(CatenaryInput()),
                      0, "analytic", "полоса 1 м натянута без слабины, 300 Н/м, 0.3 кг/м²"},
                     runCatenary, {0.02, 0.05, true}, false});
    // Holdout, blinded: value and uncertainty sealed (Blinding.h); only their commitment is here.
    cases.push_back({"liquid-dam-break", "Обрушение плотины: фронт Z при T = 2", "validation", "benchmark: dam break front",
                     {"Martin & Moyce 1952, эксперимент, n² = 2", "https://doi.org/10.1098/rsta.1952.0006", kNaN, kNaN, "experiment",
                      "u_D — наша оценка, не из статьи: чтение графика и момент старта, в квадратуре (в запечатанном файле)"},
                     runDamBreak, {0.08, 0.15, true}, false, "holdout", kDamBreakCommitment, true});
}

} // namespace rf::verify
