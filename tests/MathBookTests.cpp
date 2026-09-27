// The pictures of the book "Язык природы" (docs/math/), drawn by the engine itself.
//
// Every chapter of the book shows one piece of mathematics at work, and every picture there is
// made from a real run of one of the tests below - never drawn by hand. Each test also checks a
// number, so the book cannot quietly go out of date: if the engine changes, the test says so.
//
//   rounding     - float and double: how many digits survive a long sum and a small difference
//   derivative   - the slope of the height of a falling ball is its speed; the engine's step
//                  makes an error proportional to the step (first order)
//   schemes      - one pendulum, three ways to step it: Euler gains energy, symplectic Euler
//                  keeps it bounded, Runge-Kutta 4 keeps it almost exactly
//   vectors      - the cross product at work: a charge circles in a magnetic field, a push off
//                  the centre makes a body spin
//   linear       - a big linear system (the pressure's kind) by Jacobi and by conjugate
//                  gradients; the eigenvalues of an inertia tensor are the principal moments
//   flip         - the tennis-racket flip of a box spun about its middle axis
//   divergence   - a gas field before and after the pressure projection removes div u
//   sigma        - the spread of an average falls as 1 / sqrt(N): the root of every error bar
//
// With RF_PLOT_DIR=<folder> the tests write their curves there as CSV; tools/plot_mathbook.py
// turns them into the SVGs of docs/img/math/. Without it nothing is written.
#include "TestRunner.h"
#include "Tests.h"

#include "math/Matrix3x3.h"

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <vector>

namespace {

// ------------------------------------------------------------------------------------------------
//  Helpers
// ------------------------------------------------------------------------------------------------

// One CSV of a picture, written only when RF_PLOT_DIR names a folder.
class BookCsv {
public:
    BookCsv(const char* name, const char* header) {
        const char* folder = std::getenv("RF_PLOT_DIR");
        if (!folder || !*folder) return;
        const std::string path = std::string(folder) + "/" + name + ".csv";
        file_ = std::fopen(path.c_str(), "w");
        if (file_) std::fprintf(file_, "%s\n", header);
    }
    ~BookCsv() {
        if (file_) std::fclose(file_);
    }
    BookCsv(const BookCsv&) = delete;
    BookCsv& operator=(const BookCsv&) = delete;

    void row(std::initializer_list<double> values) {
        if (!file_) return;
        const char* separator = "";
        for (double v : values) {
            std::fprintf(file_, "%s%.9g", separator, v);
            separator = ",";
        }
        std::fprintf(file_, "\n");
    }

private:
    FILE* file_ = nullptr;
};

// A rigid world where nothing dissipates and nothing falls asleep: only what a test adds acts.
void frictionlessWorld(RigidWorld& world, const AABB& domain, const Vector3& gravity) {
    world.setDomain(domain);
    world.params.gravity = gravity;
    world.params.linearDamping = 0;
    world.params.angularDamping = 0;
    world.params.rollingResistance = 0;
    world.params.sleeping = false;
}

// A small random-number generator with a fixed seed (the linear congruential generator of
// Numerical Recipes): the same numbers on every machine, so the book prints the same figures.
struct Dice {
    uint32_t state = 12345;
    double roll() { // uniform in [0, 1)
        state = state * 1664525u + 1013904223u;
        return (state >> 8) * (1.0 / 16777216.0);
    }
};

} // namespace

// ------------------------------------------------------------------------------------------------
//  Chapter 1. Numbers and precision
// ------------------------------------------------------------------------------------------------

// A float keeps about 7 significant digits, a double about 16. Two experiments show what that means:
//   1. add 0.1 ten million times: the float total drifts, because every addition rounds to the
//      total's own spacing, which grows with the total; the double total stays right to 1e-10;
//   2. the story of the Alfven wave (docs/06): a wave of 2.618e-6 J riding on a field of 9.947 J.
//      In float the wave's energy is the difference of two nearly equal numbers, and that
//      difference is mostly rounding; in double it survives.
void testMathBookRounding() {
    static const long long kCheckpoints[] = {10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000,
                                             300000, 1000000, 3000000, 10000000};
    BookCsv csv("rounding", "n,float_relative_error,double_relative_error");
    float floatTotal = 0;
    double doubleTotal = 0, floatError = 0, doubleError = 0;
    long long n = 0;
    for (long long checkpoint : kCheckpoints) {
        for (; n < checkpoint; ++n) {
            floatTotal += 0.1f;
            doubleTotal += 0.1;
        }
        const double exact = 0.1 * double(n);
        floatError = std::fabs(floatTotal - exact) / exact;
        doubleError = std::fabs(doubleTotal - exact) / exact;
        csv.row({double(n), floatError, doubleError});
    }
    const double fieldEnergy = 9.947, waveEnergy = 2.618e-6;
    const float waveInFloat = (float(fieldEnergy) + float(waveEnergy)) - float(fieldEnergy);
    const double waveInDouble = (fieldEnergy + waveEnergy) - fieldEnergy;
    const double floatWaveError = std::fabs(waveInFloat - waveEnergy) / waveEnergy;
    const double doubleWaveError = std::fabs(waveInDouble - waveEnergy) / waveEnergy;
    std::printf("  0.1 added 1e7 times: float %.1f (error %.2f%%), double %.9f (error %.1e); a wave of %.3e J on %.3f J: "
                "float sees %.3e J (error %.0f%%), double %.6e J (error %.1e)\n",
                floatTotal, 100 * floatError, doubleTotal, doubleError, waveEnergy, fieldEnergy, waveInFloat,
                100 * floatWaveError, waveInDouble, doubleWaveError);
    CHECK(floatError > 1e-3, "the float total should drift visibly (%e)", floatError);
    CHECK(doubleError < 1e-9, "the double total drifted %e", doubleError);
    CHECK(floatWaveError > 0.01 && doubleWaveError < 1e-6, "wave energy: float %e, double %e", floatWaveError, doubleWaveError);
}

// ------------------------------------------------------------------------------------------------
//  Chapter 2. Vectors
// ------------------------------------------------------------------------------------------------

// Two cross products, both at work in the engine:
//   1. the Lorentz force F = q v x B turns a charge in a circle of radius r = m v / (q B). The
//      step is Boris's rotation (Boris 1970; Birdsall & Langdon, "Plasma Physics via Computer
//      Simulation" 4-3): it turns the velocity without changing its length;
//   2. an impulse J pushed at a point r from the centre of mass spins a body up by
//      delta omega = I^-1 (r x J) - the engine's applyImpulse, checked against the formula.
void testMathBookVectors() {
    const double charge = 1, mass = 1, field = 1, speed = 1, h = 0.05;
    const Vector3 B(0, 0, float(field));
    // Leapfrog: the position lives at whole steps, the velocity half a step between them - each
    // move x += v h uses the velocity of the middle of the step. So the velocity starts half a step
    // back: turned back by half of Boris's turn per step, 2 atan(q B h / 2m). Started at the same
    // moment as the position instead, the whole circle slides sideways by half a turn (the radius
    // measured from (0, -1) then swings 0.975 .. 1.025 m - how this test first failed).
    const double halfTurn = std::atan(0.5 * charge * field * h / mass);
    Vector3 x(0, 0, 0), v(float(speed * std::cos(halfTurn)), float(speed * std::sin(halfTurn)), 0);
    BookCsv csv("lorentz", "x_m,y_m");
    csv.row({0, 0});
    double farthest = 0, nearest = 1e9;
    const Vector3 centre(0, -1, 0); // v x B points to -y at the start: the circle's centre
    for (int k = 1; k <= 200; ++k) { // about 1.6 turns (period 2 pi m / (q B))
        const Vector3 t = B * float(0.5 * charge * h / mass);  // 1. half the turn, as a vector
        const Vector3 s = t * (2.0f / (1.0f + length2(t)));
        const Vector3 vPrime = v + cross(v, t);              // 2. turn once by t ...
        v = v + cross(vPrime, s);                            // 3. ... and finish the turn by s
        x = x + v * float(h);                                // 4. move
        const double r = length(x - centre);
        farthest = std::max(farthest, r);
        nearest = std::min(nearest, r);
        csv.row({x.x, x.y});
    }
    const double expected = mass * speed / (charge * field);

    RigidWorld world;
    frictionlessWorld(world, AABB({-5, -5, -5}, {5, 5, 5}), Vector3(0.0f));
    const int box = world.addBox({0, 0, 0}, {0.2f, 0.1f, 0.4f}, Quaternion(), 1000, Vector3(1));
    RigidBody& b = world.bodies()[size_t(box)];
    const Vector3 impulse(0, 3, 0), arm(0.4f, 0, 0);
    world.applyImpulse(box, impulse, b.pos + arm);
    const Vector3 armCrossImpulse = cross(arm, impulse);
    const Vector3 expectedSpin(armCrossImpulse.x * b.invInertiaLocal.x, armCrossImpulse.y * b.invInertiaLocal.y,
                               armCrossImpulse.z * b.invInertiaLocal.z);
    const double spinError = length(b.angVel - expectedSpin) / length(expectedSpin);
    std::printf("  Lorentz circle: radius %.5f ... %.5f m, theory m v / (q B) = %.5f m, speed kept %.7f; "
                "push off centre: omega (%.4f, %.4f, %.4f) rad/s vs I^-1 (r x J), error %.1e\n",
                nearest, farthest, expected, double(length(v)), double(b.angVel.x), double(b.angVel.y),
                double(b.angVel.z), spinError);
    CHECK(std::fabs(farthest / expected - 1) < 1e-3 && std::fabs(nearest / expected - 1) < 1e-3, "circle radius %f..%f", nearest, farthest);
    CHECK(std::fabs(length(v) - speed) < 1e-5, "the magnetic force did work: |v| = %f", double(length(v)));
    CHECK(spinError < 1e-5, "angular velocity from the push is off by %e", spinError);
}

// ------------------------------------------------------------------------------------------------
//  Chapter 3. The derivative and the integral
// ------------------------------------------------------------------------------------------------

namespace {

// Drops a ball in empty space for `seconds` at frame step dt and returns how far its height is
// from the exact y0 - g t^2 / 2 at the end. The engine steps v += g h, then y += v h (symplectic
// Euler): the ball runs ahead of the parabola by g t h / 2, an error proportional to the step.
double freeFallError(float dt, float seconds, BookCsv* csv) {
    RigidWorld world;
    frictionlessWorld(world, AABB({-50, -100, -50}, {50, 100, 50}), {0, -9.81f, 0});
    const float y0 = 50;
    const int ball = world.addSphere({0, y0, 0}, 0.1f, 1000, Vector3(1));
    const int frames = int(std::lround(seconds / dt));
    double previousHeight = y0;
    for (int k = 1; k <= frames; ++k) {
        world.step(dt);
        const RigidBody& b = world.bodies()[size_t(ball)];
        const double slope = (b.pos.y - previousHeight) / dt; // the speedometer from the height alone
        if (csv) csv->row({k * dt, b.pos.y, b.vel.y, slope});
        previousHeight = b.pos.y;
    }
    const double t = frames * double(dt);
    return world.bodies()[size_t(ball)].pos.y - (y0 - 0.5 * 9.81 * t * t);
}

} // namespace

// The derivative is a speedometer: the slope of the height is the velocity. The integral is the
// odometer: the sum of velocity times time is the distance. The test
//   1. drops a ball for 2 s and compares the slope of its height with its velocity;
//   2. halves the step three times: the error at t = 2 s halves too - the step is first order.
void testMathBookDerivative() {
    BookCsv fall("derivative", "t_s,height_m,velocity_mps,slope_mps");
    const double error = freeFallError(1.0f / 60, 2.0f, &fall);
    BookCsv orders("step_order", "dt_s,height_error_m");
    std::vector<double> errors;
    for (float dt : {1.0f / 15, 1.0f / 30, 1.0f / 60, 1.0f / 120}) {
        errors.push_back(std::fabs(freeFallError(dt, 2.0f, nullptr)));
        orders.row({dt, errors.back()});
    }
    const double order = std::log2(errors[2] / errors[3]);
    std::printf("  free fall 2 s at dt = 1/60 s: height error %.4f m; errors at dt = 1/15 ... 1/120: %.4f %.4f %.4f %.4f m, "
                "observed order %.3f\n", error, errors[0], errors[1], errors[2], errors[3], order);
    CHECK(std::fabs(order - 1) < 0.05, "the engine's step should be first order, got %f", order);
}

// ------------------------------------------------------------------------------------------------
//  Chapter 4. Equations of motion and the time step
// ------------------------------------------------------------------------------------------------

namespace {

// The pendulum theta'' = -(g / L) sin theta as two first-order equations for (theta, omega), and
// its energy per unit mass E = L^2 omega^2 / 2 + g L (1 - cos theta).
struct Pendulum {
    double g = 9.81, L = 1.0;
    double acceleration(double theta) const { return -g / L * std::sin(theta); }
    double energy(double theta, double omega) const { return 0.5 * L * L * omega * omega + g * L * (1 - std::cos(theta)); }
};

// One step of each scheme (h = the step):
//   Euler            - both from the old state:         theta += h omega,     omega += h a(theta)
//   symplectic Euler - the speed first, then the angle: omega += h a(theta),  theta += h omega
//   Runge-Kutta 4    - four slopes, averaged 1 : 2 : 2 : 1 (Kutta 1901)
void eulerStep(const Pendulum& p, double& theta, double& omega, double h) {
    const double a = p.acceleration(theta);
    theta += h * omega;
    omega += h * a;
}
void symplecticEulerStep(const Pendulum& p, double& theta, double& omega, double h) {
    omega += h * p.acceleration(theta);
    theta += h * omega;
}
void rungeKutta4Step(const Pendulum& p, double& theta, double& omega, double h) {
    const double k1t = omega, k1w = p.acceleration(theta);
    const double k2t = omega + 0.5 * h * k1w, k2w = p.acceleration(theta + 0.5 * h * k1t);
    const double k3t = omega + 0.5 * h * k2w, k3w = p.acceleration(theta + 0.5 * h * k2t);
    const double k4t = omega + h * k3w, k4w = p.acceleration(theta + h * k3t);
    theta += h / 6 * (k1t + 2 * k2t + 2 * k3t + k4t);
    omega += h / 6 * (k1w + 2 * k2w + 2 * k3w + k4w);
}

} // namespace

// The same pendulum (1 m, released at 60 degrees) stepped by three schemes at h = 0.05 s for 20 s.
// Nature keeps the energy constant; the schemes do not all agree. The test checks that
//   1. Euler gains energy without bound (it spirals outwards);
//   2. symplectic Euler never drifts: its error stays within a few percent, forever;
//   3. Runge-Kutta 4 keeps the energy to 1e-3 over the run.
void testMathBookSchemes() {
    const Pendulum p;
    const double h = 0.05, start = 60 * kPi / 180, E0 = p.energy(start, 0);
    double eulerTheta = start, eulerOmega = 0, symplecticTheta = start, symplecticOmega = 0, rkTheta = start, rkOmega = 0;
    double symplecticWorst = 0, rkWorst = 0;
    BookCsv csv("schemes", "t_s,euler,symplectic_euler,runge_kutta_4");
    csv.row({0, 1, 1, 1});
    for (int k = 1; k <= 400; ++k) {
        eulerStep(p, eulerTheta, eulerOmega, h);
        symplecticEulerStep(p, symplecticTheta, symplecticOmega, h);
        rungeKutta4Step(p, rkTheta, rkOmega, h);
        const double eulerRatio = p.energy(eulerTheta, eulerOmega) / E0;
        const double symplecticRatio = p.energy(symplecticTheta, symplecticOmega) / E0;
        const double rkRatio = p.energy(rkTheta, rkOmega) / E0;
        symplecticWorst = std::max(symplecticWorst, std::fabs(symplecticRatio - 1));
        rkWorst = std::max(rkWorst, std::fabs(rkRatio - 1));
        csv.row({k * h, eulerRatio, symplecticRatio, rkRatio});
    }
    const double eulerFinal = p.energy(eulerTheta, eulerOmega) / E0;
    std::printf("  pendulum 60 deg, h = 0.05 s, 20 s: energy E/E0 at the end - Euler %.3f; worst |E/E0 - 1| - symplectic "
                "Euler %.4f, Runge-Kutta 4 %.1e\n", eulerFinal, symplecticWorst, rkWorst);
    CHECK(eulerFinal > 1.5, "Euler should gain energy, E/E0 = %f", eulerFinal);
    CHECK(symplecticWorst < 0.1, "symplectic Euler drifted %f", symplecticWorst);
    CHECK(rkWorst < 1e-3, "Runge-Kutta 4 drifted %e", rkWorst);
}

// ------------------------------------------------------------------------------------------------
//  Chapter 5. Matrices and systems of equations
// ------------------------------------------------------------------------------------------------

namespace {

// The Poisson equation -laplace p = f on an n x n grid with p = 0 on the border, the kind of system
// the gas solves for its pressure every step. apply() is the matrix A times a vector: the
// five-point stencil 4 p_ij - p_(i-1)j - p_(i+1)j - p_i(j-1) - p_i(j+1) (h^2 folded into f).
struct PoissonGrid {
    int n;
    double at(const std::vector<double>& p, int i, int j) const {
        return (i < 0 || j < 0 || i >= n || j >= n) ? 0.0 : p[size_t(j) * n + i];
    }
    void apply(const std::vector<double>& p, std::vector<double>& out) const {
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                out[size_t(j) * n + i] = 4 * at(p, i, j) - at(p, i - 1, j) - at(p, i + 1, j) - at(p, i, j - 1) - at(p, i, j + 1);
    }
};

double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0;
    for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

// Jacobi: every unknown takes the value its own equation wants, all at once - p = (f + neighbours) / 4.
// Returns the relative residual |f - A p| / |f| after each sweep.
std::vector<double> solveJacobi(const PoissonGrid& grid, const std::vector<double>& f, int sweeps) {
    std::vector<double> p(f.size(), 0.0), next(f.size()), Ap(f.size()), history;
    const double fNorm = std::sqrt(dotProduct(f, f));
    for (int s = 0; s < sweeps; ++s) {
        for (int j = 0; j < grid.n; ++j)
            for (int i = 0; i < grid.n; ++i)
                next[size_t(j) * grid.n + i] = (f[size_t(j) * grid.n + i] + grid.at(p, i - 1, j) + grid.at(p, i + 1, j) +
                                                grid.at(p, i, j - 1) + grid.at(p, i, j + 1)) / 4;
        p.swap(next);
        grid.apply(p, Ap);
        double r2 = 0;
        for (size_t k = 0; k < f.size(); ++k) r2 += (f[k] - Ap[k]) * (f[k] - Ap[k]);
        history.push_back(std::sqrt(r2) / fNorm);
    }
    return history;
}

// Conjugate gradients (Hestenes & Stiefel 1952): each step goes down the hill of
// E(p) = p.A p / 2 - f.p along a direction that does not spoil the earlier ones.
//   1. r = f - A p, the residual; d = r, the first direction
//   2. alpha = r.r / d.A d: how far to go along d
//   3. p += alpha d, r -= alpha A d
//   4. beta = r.r(new) / r.r(old); d = r + beta d: the next direction
std::vector<double> solveConjugateGradients(const PoissonGrid& grid, const std::vector<double>& f, int steps) {
    std::vector<double> p(f.size(), 0.0), r = f, d = f, Ad(f.size()), history;
    const double fNorm = std::sqrt(dotProduct(f, f));
    double rr = dotProduct(r, r);
    for (int s = 0; s < steps && rr > 0; ++s) {
        grid.apply(d, Ad);
        const double alpha = rr / dotProduct(d, Ad);
        for (size_t k = 0; k < p.size(); ++k) {
            p[k] += alpha * d[k];
            r[k] -= alpha * Ad[k];
        }
        const double rrNew = dotProduct(r, r);
        for (size_t k = 0; k < d.size(); ++k) d[k] = r[k] + rrNew / rr * d[k];
        rr = rrNew;
        history.push_back(std::sqrt(rr) / fNorm);
    }
    return history;
}

} // namespace

// Two faces of linear algebra:
//   1. a system of 1024 equations (a 32 x 32 pressure grid) solved by Jacobi and by conjugate
//      gradients: the second reaches a millionth of the residual in tens of steps, the first is
//      still far after hundreds - why the engine's pressure solver is a (multigrid) CG;
//   2. the inertia tensor of a turned box: its eigenvalues are the box's principal moments
//      m (b^2 + c^2) / 12 whatever the turn - the axes a body likes to spin about.
void testMathBookLinearSystems() {
    const PoissonGrid grid{32};
    std::vector<double> f(size_t(grid.n) * grid.n);
    for (int j = 0; j < grid.n; ++j)
        for (int i = 0; i < grid.n; ++i) // a smooth load with a sharp spot: both kinds of error
            f[size_t(j) * grid.n + i] = std::sin(kPi * (i + 1) / (grid.n + 1)) * std::sin(kPi * (j + 1) / (grid.n + 1)) +
                                        ((i == 8 && j == 20) ? 1.0 : 0.0);
    const std::vector<double> jacobi = solveJacobi(grid, f, 300);
    const std::vector<double> cg = solveConjugateGradients(grid, f, 300);
    BookCsv csv("linear_systems", "iteration,jacobi,conjugate_gradients");
    int cgSteps = 0;
    for (size_t k = 0; k < jacobi.size(); ++k) {
        const double cgResidual = k < cg.size() ? std::max(cg[k], 1e-16) : 1e-16;
        if (cgSteps == 0 && cgResidual < 1e-6) cgSteps = int(k) + 1;
        csv.row({double(k + 1), jacobi[k], cgResidual});
    }

    const Vector3 half(0.1f, 0.3f, 0.5f);
    const float mass = 8 * half.x * half.y * half.z * 1000;
    const Vector3 principal(mass * (4 * half.y * half.y + 4 * half.z * half.z) / 12,
                            mass * (4 * half.x * half.x + 4 * half.z * half.z) / 12,
                            mass * (4 * half.x * half.x + 4 * half.y * half.y) / 12);
    const Matrix3x3 R = Quaternion::fromAxisAngle(normalize(Vector3(1, 2, 3)), 0.7f).toMatrix3x3();
    const Matrix3x3 inertiaTurned = R * Matrix3x3::diag(principal) * R.transposed();
    Vector3 eigenvalues;
    Matrix3x3 axes;
    symmetricEigen(inertiaTurned, eigenvalues, axes);
    float worst = 0;
    for (int a = 0; a < 3; ++a) { // the eigenvalues come in any order: match each to the nearest moment
        float best = 1e30f;
        for (int b = 0; b < 3; ++b) best = std::min(best, std::fabs(eigenvalues[a] - principal[b]) / principal[b]);
        worst = std::max(worst, best);
    }
    std::printf("  32 x 32 pressure system: CG reaches 1e-6 in %d steps, Jacobi after 300 sweeps is at %.2e; "
                "inertia of a turned box: eigenvalues %.4f %.4f %.4f vs moments %.4f %.4f %.4f kg m^2 (worst %.1e)\n",
                cgSteps, jacobi.back(), double(eigenvalues.x), double(eigenvalues.y), double(eigenvalues.z),
                double(principal.x), double(principal.y), double(principal.z), double(worst));
    CHECK(cgSteps > 0 && cgSteps < 100, "conjugate gradients took %d steps", cgSteps);
    CHECK(jacobi.back() > 1e-2, "Jacobi should still be far after 300 sweeps (%e)", jacobi.back());
    CHECK(worst < 1e-5f, "the eigenvalues missed the principal moments by %e", double(worst));
}

// ------------------------------------------------------------------------------------------------
//  Chapter 6. Quaternions and rotations
// ------------------------------------------------------------------------------------------------

namespace {

// Euler's equations of a free body in its own principal axes (Euler 1765):
//   I1 w1' = (I2 - I3) w2 w3,   I2 w2' = (I3 - I1) w3 w1,   I3 w3' = (I1 - I2) w1 w2.
// Their right-hand side for the spin w, the principal moments I.
Vector3 eulerTorqueFree(const Vector3& w, const Vector3& I) {
    return {(I.y - I.z) * w.y * w.z / I.x, (I.z - I.x) * w.z * w.x / I.y, (I.x - I.y) * w.x * w.y / I.z};
}

// The exact reference: Euler's equations stepped by Runge-Kutta 4 at a step a thousand times finer
// than the engine's, so fine that its own error is far below what the picture can show. Returns
// the spin about the middle axis at the frame times k * frame, k = 0 ... frames.
std::vector<double> exactMiddleSpin(Vector3 w, const Vector3& I, float frame, int frames) {
    const int fine = 1000;
    const float h = frame / fine;
    std::vector<double> spin{w.y};
    for (int k = 0; k < frames; ++k) {
        for (int s = 0; s < fine; ++s) {
            const Vector3 k1 = eulerTorqueFree(w, I);
            const Vector3 k2 = eulerTorqueFree(w + k1 * (0.5f * h), I);
            const Vector3 k3 = eulerTorqueFree(w + k2 * (0.5f * h), I);
            const Vector3 k4 = eulerTorqueFree(w + k3 * h, I);
            w += (k1 + k2 * 2.0f + k3 * 2.0f + k4) * (h / 6);
        }
        spin.push_back(w.y);
    }
    return spin;
}

int signChanges(const std::vector<double>& values) {
    int changes = 0;
    for (size_t k = 1; k < values.size(); ++k) changes += (values[k] > 0) != (values[k - 1] > 0);
    return changes;
}

} // namespace

// A box with half extents 0.1 x 0.3 x 0.5 m spun about its middle axis (y) with a tiny nudge about
// the others turns over again and again (Dzhanibekov 1985; Ashbaugh, Chicone & Cushman 1991). The
// orientation is a unit quaternion. The test
//   1. spins the box in the engine for 20 s and records the spin in the box's own axes;
//   2. steps Euler's equations exactly (fine Runge-Kutta 4) from the same start;
//   3. counts the turn-overs (sign changes of the middle-axis spin) of both, and checks that the
//      quaternion stays of length 1 - a rotation, never a stretch.
// The engine turns over more often than the exact motion: the known defect of its rotation step
// (docs/12, docs/13: a symplectic splitting is the planned cure). Printed, not hidden.
void testMathBookFlip() {
    RigidWorld world;
    frictionlessWorld(world, AABB({-50, -50, -50}, {50, 50, 50}), Vector3(0.0f));
    const int box = world.addBox({0, 0, 0}, {0.1f, 0.3f, 0.5f}, Quaternion(), 1000, Vector3(1));
    RigidBody& b = world.bodies()[size_t(box)];
    const Vector3 startSpin(0.01f, 4.0f, 0.01f);
    b.angVel = startSpin;
    const Vector3 moments(1 / b.invInertiaLocal.x, 1 / b.invInertiaLocal.y, 1 / b.invInertiaLocal.z);
    const float frame = 1.0f / 60;
    const int frames = 20 * 60;
    const std::vector<double> exact = exactMiddleSpin(startSpin, moments, frame, frames);
    std::vector<double> engine{startSpin.y};
    double worstLength = 0;
    BookCsv csv("flip", "t_s,engine_omega_x,engine_omega_y,engine_omega_z,exact_omega_y");
    for (int k = 1; k <= frames; ++k) {
        world.step(frame);
        const Vector3 spinInBody = b.rotation().transposed() * b.angVel;
        engine.push_back(spinInBody.y);
        const double length = std::sqrt(double(b.rot.w) * b.rot.w + double(b.rot.x) * b.rot.x + double(b.rot.y) * b.rot.y +
                                        double(b.rot.z) * b.rot.z);
        worstLength = std::max(worstLength, std::fabs(length - 1));
        if (k % 2 == 0) csv.row({k * frame, spinInBody.x, spinInBody.y, spinInBody.z, exact[size_t(k)]});
    }
    const int engineTurns = signChanges(engine), exactTurns = signChanges(exact);
    std::printf("  tennis-racket box, 20 s: turned over %d times in the engine, %d times by Euler's equations "
                "(exact); |q| - 1 at worst %.1e\n", engineTurns, exactTurns, worstLength);
    CHECK(exactTurns >= 2 && engineTurns >= 2, "the box should turn over (engine %d, exact %d)", engineTurns, exactTurns);
    CHECK(worstLength < 1e-5, "the quaternion lost its unit length (%e)", worstLength);
}

// ------------------------------------------------------------------------------------------------
//  Chapter 7. Fields: gradient, divergence, curl
// ------------------------------------------------------------------------------------------------

// A gas that cannot be squeezed must have div u = 0 everywhere: what flows into a cell flows out.
// The test
//   1. sets a velocity u = (sin(pi x) sin(pi y), 0, 0) whose divergence pi cos(pi x) sin(pi y) is
//      far from zero;
//   2. takes one tiny step of the gas: its pressure projection subtracts grad p, the part that
//      squeezes;
//   3. measures div u along the middle row before and after.
void testMathBookDivergence() {
    GasSolver gas;
    const int cells = 32;
    const float dx = 1.0f / cells;
    gas.params.domainSize = {1, 1, 4 * dx};
    gas.params.resolutionX = cells;
    gas.params.inflowSpeed = 0;
    gas.params.smokeRake = false;
    gas.params.kinematicViscosity = 0;
    gas.params.wallFriction = false;
    for (auto& bc : gas.params.bc) bc = BoundaryType::Wall;
    gas.reset({0, 0, 0}, nullptr);
    gas.setVelocity([](const Vector3& x) { return Vector3(std::sin(kPi * x.x) * std::sin(kPi * x.y), 0, 0); });
    const int row = cells / 2, slab = 1;
    std::vector<float> before(static_cast<size_t>(cells)), after(static_cast<size_t>(cells));
    float largestBefore = 0, largestAfter = 0;
    for (int i = 0; i < cells; ++i) before[size_t(i)] = gas.cellDivergence(i, row, slab);
    gas.step(1e-4f);
    for (int i = 0; i < cells; ++i) after[size_t(i)] = gas.cellDivergence(i, row, slab);
    BookCsv csv("divergence", "x_m,divergence_before,divergence_after");
    for (int i = 0; i < cells; ++i) {
        largestBefore = std::max(largestBefore, std::fabs(before[size_t(i)]));
        largestAfter = std::max(largestAfter, std::fabs(after[size_t(i)]));
        csv.row({(i + 0.5) * dx, before[size_t(i)], after[size_t(i)]});
    }
    std::printf("  divergence along the middle row: largest %.3f 1/s before the projection, %.1e 1/s after (%.0e of it)\n",
                double(largestBefore), double(largestAfter), double(largestAfter / largestBefore));
    CHECK(largestBefore > 1.0f, "the start field should squeeze (%f)", double(largestBefore));
    CHECK(largestAfter < 1e-3f * largestBefore, "the projection left div u = %e", double(largestAfter));
}

// ------------------------------------------------------------------------------------------------
//  Chapter 10. Statistics and sigma
// ------------------------------------------------------------------------------------------------

// Throw N dice (uniform numbers in [0, 1)) and average them; repeat the experiment 400 times. The
// averages scatter around 1/2 with the spread sigma = 1 / sqrt(12 N) (the central limit theorem;
// 1/12 is the variance of one die). The test checks
//   1. the measured spread against 1 / sqrt(12 N) for N = 4 ... 4096 (15 %);
//   2. for N = 1024: about 68 % of the averages lie within one sigma, 95 % within two.
void testMathBookSigma() {
    Dice dice;
    BookCsv csv("sigma", "n,measured_spread,theory_spread");
    const int repeats = 400;
    double worst = 0, withinOne = 0, withinTwo = 0;
    for (int n : {4, 16, 64, 256, 1024, 4096}) {
        std::vector<double> averages;
        for (int r = 0; r < repeats; ++r) {
            double sum = 0;
            for (int k = 0; k < n; ++k) sum += dice.roll();
            averages.push_back(sum / n);
        }
        double mean = 0, variance = 0;
        for (double a : averages) mean += a / repeats;
        for (double a : averages) variance += (a - mean) * (a - mean) / (repeats - 1);
        const double spread = std::sqrt(variance), theory = 1 / std::sqrt(12.0 * n);
        worst = std::max(worst, std::fabs(spread / theory - 1));
        if (n == 1024)
            for (double a : averages) {
                withinOne += std::fabs(a - 0.5) < theory ? 1.0 / repeats : 0.0;
                withinTwo += std::fabs(a - 0.5) < 2 * theory ? 1.0 / repeats : 0.0;
            }
        csv.row({double(n), spread, theory});
    }
    std::printf("  averages of N dice: spread vs 1/sqrt(12 N) at worst %.1f%% off; N = 1024: %.1f%% within 1 sigma, "
                "%.1f%% within 2 sigma (theory 68.3%%, 95.4%%)\n", 100 * worst, 100 * withinOne, 100 * withinTwo);
    CHECK(worst < 0.15, "the spread missed 1/sqrt(12 N) by %f", worst);
    CHECK(std::fabs(withinOne - 0.683) < 0.06 && std::fabs(withinTwo - 0.954) < 0.04, "within 1 and 2 sigma: %f %f", withinOne, withinTwo);
}
