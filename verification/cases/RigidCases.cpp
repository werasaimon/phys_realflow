// Verification cases of the rigid-body solver (RigidWorld): exact answers of mechanics that the
// integrator, the joints, the friction and the contacts must reproduce.
//   projectile      - free flight, x(t) = x0 + v0 t + g t^2 / 2; the order of the integrator in dt
//   incline         - a block on a slope: sliding at a = g (sin theta - mu_k cos theta), sticking
//                     below tan theta = mu_s (Coulomb)
//   pendulum        - the period of a compound pendulum at 10 degrees (elliptic integral), and
//                     its convergence as the time step shrinks
//   Noether         - energy, momentum and angular momentum stay what they were
//   Newton's cradle - a hit passes down a row of elastic balls
// The scenes are built the way tests/RigidTests.cpp builds them; damping is off, because what is
// verified is the solver of the equations, not the damping model.
#include "../Benchmark.h"
#include "../Convergence.h"

#include "core/Format.h"
#include "rigid/Joints.h"
#include "rigid/RigidWorld.h"

#include <cmath>

namespace rf::verify {

static const double kPiD = 3.14159265358979323846;

// A world without damping and sleeping: what is left is the integrator and the solver.
static void undamped(RigidWorld& w) {
    w.params.sleeping = false;
    w.params.linearDamping = w.params.angularDamping = 0;
    w.params.rollingResistance = 0;
}

// ---------------------------------------------------------------------------------------------
// Projectile: the semi-implicit (symplectic) Euler step v += g dt, x += v dt gives
// x_n = x0 + v0 t + g t^2 / 2 + g t dt / 2 - first order in dt, whatever the forces.
static double projectileError(float dt, double tEnd) {
    RigidWorld w;
    w.setDomain(AABB({-50, -50, -50}, {50, 50, 50}));
    undamped(w);
    w.params.collideWithDomain = false;
    w.params.ccd = false;
    const Vector3 v0(3, 4, 0);
    const int s = w.addSphere({0, 0, 0}, 0.1f, 1000, Vector3(1));
    w.bodies()[s].vel = v0;
    const int steps = int(std::lround(tEnd / dt));
    for (int k = 0; k < steps; ++k) w.step(dt);
    const double t = steps * double(dt), g = w.params.gravity.y;
    const double ex = v0.x * t, ey = v0.y * t + 0.5 * g * t * t;
    const Vector3 p = w.bodies()[s].pos;
    return std::hypot(p.x - ex, p.y - ey);
}

static Result runProjectile(const RunOptions& o) {
    Result r;
    r.unit = "";
    r.hLabel = "dt, с";
    r.theoreticalOrder = 1;
    const int levels = o.full ? 5 : 3;
    for (int l = 0; l < levels; ++l) {
        const float dt = 1.0f / 60 / float(1 << l);
        const double e = projectileError(dt, 0.8);
        r.convergence.push_back({dt, e, e});
    }
    const OrderFit fit = fitOrder(r.convergence);
    r.value = r.observedOrder = fit.order;
    r.numericalUncertainty = std::isfinite(fit.stdError) ? fit.stdError : 0;
    const ConvergencePoint& fine = r.convergence.back();
    r.detail = format("ошибка положения за 0.8 с: %.2e м при dt = 1/60, %.2e м при dt = %.1e с; ожидание g t dt / 2",
                      r.convergence.front().error, fine.error, fine.h);
    return r;
}

// ---------------------------------------------------------------------------------------------
// A block on an incline: the slope is the horizontal floor with gravity turned by theta (the same
// physics, no placement error). Both bodies carry the same coefficients: the solver mixes them as
// sqrt(mu_a mu_b), so the pair has exactly mu_k and mu_s.
struct InclineRun {
    double acceleration = 0; // along the slope [m/s^2], from the velocity change
    double slide = 0;        // distance slid [m]
};

static InclineRun runIncline(double thetaDeg, float muK, float muS, double seconds) {
    RigidWorld w;
    w.setDomain(AABB({-20, -5, -5}, {20, 5, 5}));
    undamped(w);
    const double th = thetaDeg * kPiD / 180.0, g = 9.81;
    w.params.gravity = Vector3(float(g * std::sin(th)), float(-g * std::cos(th)), 0);
    const int floor = w.addBox({0, -0.5f, 0}, {15, 0.5f, 2}, Quaternion(), 0, Vector3(1));
    const int block = w.addBox({-10, 0.1f, 0}, Vector3(0.1f), Quaternion(), 1000, Vector3(1));
    for (int b : {floor, block}) {
        w.bodies()[b].friction = muK;
        w.bodies()[b].staticFriction = muS;
        w.bodies()[b].restitution = 0;
    }
    const float dt = 1.0f / 600;
    auto advance = [&](double t) { for (int k = 0; k < int(std::lround(t / dt)); ++k) w.step(dt); };
    advance(0.1); // contact established
    const double x0 = w.bodies()[block].pos.x, v0 = w.bodies()[block].vel.x;
    advance(seconds);
    InclineRun out;
    out.acceleration = (w.bodies()[block].vel.x - v0) / seconds;
    out.slide = w.bodies()[block].pos.x - x0;
    return out;
}

static Result runInclineSlide(const RunOptions&) {
    Result r;
    r.unit = "м/с²";
    const InclineRun run = runIncline(30, 0.3f, 0.4f, 1.0);
    r.value = run.acceleration;
    r.detail = format("склон 30°, μk = 0.3, μs = 0.4: a = %.4f м/с², путь за 1 с %.3f м", run.acceleration, run.slide);
    return r;
}

// The steepest slope on which the block stays: bisection between 15 and 30 degrees; "stays" =
// slides less than 5 mm in 1 s from rest.
static Result runInclineStick(const RunOptions& o) {
    Result r;
    r.unit = "°";
    double lo = 15, hi = 30;
    const int halvings = o.full ? 10 : 7;
    for (int k = 0; k < halvings; ++k) {
        const double mid = 0.5 * (lo + hi);
        (runIncline(mid, 0.3f, 0.4f, 1.0).slide < 0.005 ? lo : hi) = mid;
    }
    r.value = 0.5 * (lo + hi);
    r.numericalUncertainty = 0.5 * (hi - lo) / std::sqrt(3.0); // the bracket as a uniform distribution
    r.detail = format("граница покоя между %.3f° и %.3f° (μs = 0.4, μk = 0.3)", lo, hi);
    return r;
}

// ---------------------------------------------------------------------------------------------
// Pendulum: a steel sphere (r = 5 cm) on a ball joint 1 m below the pivot, released at 10 degrees.
// A compound pendulum: I = m (L^2 + 2 r^2 / 5) about the pivot, small-angle period
// T0 = 2 pi sqrt(I / (m g L)); at amplitude theta0 the exact period is T = T0 (2 / pi) K(sin(theta0/2)),
// K the complete elliptic integral of the first kind (from the arithmetic-geometric mean).
static double ellipticK(double k) {
    double a = 1, b = std::sqrt(1 - k * k);
    for (int i = 0; i < 30; ++i) { const double an = 0.5 * (a + b); b = std::sqrt(a * b); a = an; }
    return kPiD / (2 * a);
}

static double exactPendulumPeriod() {
    const double L = 1, r = 0.05, g = 9.81, theta0 = 10 * kPiD / 180;
    const double T0 = 2 * kPiD * std::sqrt((L * L + 0.4 * r * r) / (g * L));
    return T0 * (2 / kPiD) * ellipticK(std::sin(theta0 / 2));
}

// The period from the upward zero crossings of the bob's x velocity, interpolated inside the step.
static double pendulumPeriod(int substeps) {
    RigidWorld w;
    w.setDomain(AABB({-10, -10, -10}, {10, 10, 10}));
    undamped(w);
    const double theta0 = 10 * kPiD / 180;
    const int bob = w.addSphere({float(std::sin(theta0)), float(5 - std::cos(theta0)), 0}, 0.05f, 7800, Vector3(1));
    w.addBallJoint(bob, -1, {0, 5, 0});
    const float dt = 1.0f / 60 / float(substeps);
    double t = 0, prev = 0, first = -1, last = -1;
    int crossings = 0;
    for (int k = 0; t < 7.0; ++k) {
        w.step(dt);
        t += dt;
        const double vx = w.bodies()[bob].vel.x;
        if (k > 0 && prev < 0 && vx >= 0) {
            const double tc = t - dt * vx / (vx - prev); // linear between the two steps
            if (first < 0) first = tc;
            last = tc;
            ++crossings;
        }
        prev = vx;
    }
    return crossings > 1 ? (last - first) / (crossings - 1) : kNaN;
}

static Result runPendulum(const RunOptions& o) {
    Result r;
    r.unit = "с";
    r.hLabel = "dt, с";
    const double exact = exactPendulumPeriod();
    const int levels = o.full ? 4 : 3;
    for (int l = 0; l < levels; ++l) {
        const int sub = 1 << l;
        const double T = pendulumPeriod(sub);
        r.convergence.push_back({1.0 / 60 / sub, T, std::fabs(T - exact)});
    }
    const size_t n = r.convergence.size();
    const RichardsonEstimate re = richardson(r.convergence[n - 1].value, r.convergence[n - 2].value, r.convergence[n - 3].value, 2.0);
    r.value = r.convergence.back().value;
    r.numericalUncertainty = numericalUncertaintyFromGci(re.gciFine);
    r.observedOrder = fitOrder(r.convergence).order;
    r.detail = format("T = %.5f / %.5f / %.5f с (подшагов %d..%d), точное %.5f с; Ричардсон %.5f, p = %.2f",
                      r.convergence[n - 3].value, r.convergence[n - 2].value, r.convergence[n - 1].value, 1 << (levels - 3), 1 << (levels - 1),
                      exact, re.extrapolated, re.order);
    return r;
}

// ---------------------------------------------------------------------------------------------
// Noether: an elastic ball (e = 1) dropped from 1 m keeps its energy over 20 bounces.
static Result runNoetherEnergy(const RunOptions&) {
    RigidWorld w;
    w.setDomain(AABB({-2, 0, -2}, {2, 5, 2}));
    undamped(w);
    const float r = 0.1f, dt = 1.0f / 600;
    const int s = w.addSphere({0, 1.0f + r, 0}, r, 1000, Vector3(1));
    RigidBody& b = w.bodies()[s];
    b.restitution = 1.0f;
    b.friction = b.staticFriction = 0;
    const double m = b.mass, g = -w.params.gravity.y;
    auto energy = [&] { return m * g * (b.pos.y - r) + 0.5 * m * length2(b.vel); };
    const double E0 = energy();
    int bounces = 0;
    double worst = 0, prevVy = 0;
    for (int k = 0; k < 600 * 20 && bounces < 20; ++k) {
        w.step(dt);
        if (prevVy < 0 && b.vel.y > 0) {
            ++bounces;
            worst = std::max(worst, std::fabs(energy() / E0 - 1));
        }
        prevVy = b.vel.y;
    }
    Result out;
    out.value = bounces >= 20 ? worst : kNaN;
    out.detail = format("худший дрейф энергии %.2f%% за %d отскоков (dt = 1/600)", 100 * worst, bounces);
    return out;
}

// Two boxes in zero gravity, a glancing hit that spins them: the pair impulses cancel, P stays.
static Result runNoetherMomentum(const RunOptions&) {
    RigidWorld w;
    w.setDomain(AABB({-10, -10, -10}, {10, 10, 10}));
    undamped(w);
    w.params.gravity = Vector3(0.0f);
    w.params.collideWithDomain = false;
    const int a = w.addBox({-1, 0, 0}, {0.15f, 0.1f, 0.2f}, Quaternion(), 800, Vector3(1));
    const int c = w.addBox({1, 0.12f, 0.05f}, {0.2f, 0.15f, 0.1f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.4f), 500, Vector3(1));
    w.bodies()[a].vel = {3, 0, 0};
    w.bodies()[a].angVel = {0, 0, 2};
    w.bodies()[c].vel = {-1, 0, 0};
    auto momentum = [&] { return w.bodies()[a].vel * w.bodies()[a].mass + w.bodies()[c].vel * w.bodies()[c].mass; };
    const Vector3 P0 = momentum();
    double worst = 0;
    for (int k = 0; k < 600; ++k) {
        w.step(1.0f / 600);
        worst = std::max(worst, double(length(momentum() - P0) / length(P0)));
    }
    Result out;
    out.value = worst;
    out.detail = format("худший |ΔP|/|P| = %.1e за 1 с, спин второго ящика %.2f рад/с", worst, length(w.bodies()[c].angVel));
    return out;
}

// A free box spun about its middle axis flips over and over (Dzhanibekov); |L| stays.
static Result runNoetherAngular(const RunOptions&) {
    RigidWorld w;
    w.setDomain(AABB({-10, -10, -10}, {10, 10, 10}));
    undamped(w);
    w.params.gravity = Vector3(0.0f);
    w.params.collideWithDomain = false;
    const int i = w.addBox({0, 0, 0}, {0.05f, 0.2f, 0.1f}, Quaternion(), 1000, Vector3(1));
    RigidBody& b = w.bodies()[i];
    b.angVel = {0.2f, 0, 10.0f};
    auto angularMomentum = [&] {
        const Matrix3x3 R = b.rotation();
        const Matrix3x3 I = R * Matrix3x3::diag({1 / b.invInertiaLocal.x, 1 / b.invInertiaLocal.y, 1 / b.invInertiaLocal.z}) * R.transposed();
        return length(I * b.angVel);
    };
    const double L0 = angularMomentum();
    double worst = 0, prevWz = b.angVel.z;
    int flips = 0;
    for (int k = 0; k < 6000; ++k) {
        w.step(1.0f / 600);
        worst = std::max(worst, std::fabs(angularMomentum() / L0 - 1));
        const double wz = dot(b.rotation() * Vector3(0, 0, 1), b.angVel);
        flips += prevWz * wz < 0;
        prevWz = wz;
    }
    Result out;
    out.value = flips >= 2 ? worst : kNaN;
    out.detail = format("|L| дрейф %.1e за 10 с, %d переворотов (эффект Джанибекова)", worst, flips);
    return out;
}

// Six steel balls 1 cm apart, e = 1, no friction; the first comes in at 2 m/s.
static Result runNewtonCradle(const RunOptions&) {
    RigidWorld w;
    w.setDomain(AABB({-2, 0, -1}, {2, 2, 1}));
    undamped(w);
    const float r = 0.05f, gap = 0.01f;
    for (int i = 0; i < 6; ++i) {
        const int b = w.addSphere({-0.5f + i * (2 * r + gap), r, 0}, r, 7800, Vector3(1));
        w.bodies()[b].restitution = 1.0f;
        w.bodies()[b].friction = w.bodies()[b].staticFriction = 0;
    }
    w.bodies()[0].vel = {2, 0, 0};
    for (int k = 0; k < 30 * w.params.substeps; ++k) w.step(1.0f / 60 / w.params.substeps);
    double others = 0;
    for (int i = 0; i < 5; ++i) others = std::max(others, double(std::fabs(w.bodies()[i].vel.x)));
    Result out;
    out.unit = "м/с";
    out.value = w.bodies()[5].vel.x;
    out.detail = format("последний шар %.4f м/с, остальные не быстрее %.4f м/с", out.value, others);
    return out;
}

// ---------------------------------------------------------------------------------------------
void addRigidCases(std::vector<Case>& cases) {
    const Reference noether{"Нётер (1918): сохранение при симметрии", "https://doi.org/10.1080/00411457108231446", 0, 0, "analytic",
                            "точный закон сохранения; допуск — наш"};
    cases.push_back({"rigid-projectile", "Полёт тела: порядок интегратора по dt", "code-verification", "",
                     {"Полусимплектический Эйлер: ошибка g t dt / 2 (Hairer, Lubich, Wanner 2006)", "https://doi.org/10.1007/3-540-30666-8", 1, 0,
                      "analytic", "теоретический порядок схемы"},
                     runProjectile, {0.1, 0.3, false}, false});
    cases.push_back({"rigid-incline-slide", "Брусок скользит по склону 30°", "code-verification", "",
                     {"Кулон: a = g (sin θ − μk cos θ)", "", 9.81 * (0.5 - 0.3 * std::sqrt(3.0) / 2), 0, "analytic", "μk = 0.3"},
                     runInclineSlide, {0.01, 0.03, true}, false});
    cases.push_back({"rigid-incline-stick", "Брусок на склоне: угол срыва", "code-verification", "",
                     {"Кулон: tg θc = μs", "", std::atan(0.4) * 180 / kPiD, 0, "analytic", "μs = 0.4"},
                     runInclineStick, {0.25, 1.0, false}, false});
    cases.push_back({"rigid-noether-energy", "Энергия упругого мяча, 20 отскоков", "code-verification", "Noether", noether,
                     runNoetherEnergy, {0.02, 0.05, false}, false});
    cases.push_back({"rigid-noether-momentum", "Импульс в косом ударе", "code-verification", "Noether", noether,
                     runNoetherMomentum, {1e-4, 1e-3, false}, false});
    cases.push_back({"rigid-noether-angular", "|L| свободного тела (Джанибеков)", "code-verification", "Noether", noether,
                     runNoetherAngular, {1e-3, 1e-2, false}, false});
    cases.push_back({"rigid-newton-cradle", "Колыбель Ньютона на полу, e = 1", "code-verification", "Newton's cradle on the floor",
                     {"Сохранение импульса и энергии: последний шар уходит с 2 м/с", "", 2.0, 0, "analytic", ""},
                     runNewtonCradle, {0.01, 0.05, true}, false});
    cases.push_back({"rigid-pendulum-dt", "Маятник 10°: период и сходимость по dt", "solution-verification", "joints",
                     {"Физический маятник, эллиптический интеграл K", "", exactPendulumPeriod(), 0, "analytic",
                      "T = T0 (2/π) K(sin θ0/2), I = m (L² + 2r²/5)"},
                     runPendulum, {1e-3, 5e-3, true}, false});
}

} // namespace rf::verify
