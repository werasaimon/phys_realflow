// The action principle, checked by energy balances (docs/11-action.md). What has an action keeps
// what its symmetries say (Noether): a free body keeps its energy and its angular momentum. What
// dissipates has no action, only the Lagrange-d'Alembert principle with a Rayleigh function R: the
// energy then falls by exactly the work of the dissipative forces, dE/dt = -2R. Each test runs one
// physics, adds up the dissipation the way the theory counts it, and compares the two ledgers:
//   friction   - a block sliding to rest: lost energy = mu m g * sliding distance
//   damping    - a pendulum on a ball joint with linear damping: lost energy = integral c m v^2 dt
//                (the joint's multiplier force does no work)
//   viscosity  - a Taylor-Green vortex: lost kinetic energy vs 2 nu rho integral |S|^2 dV dt
//   resistance - a decaying magnetic field: lost magnetic energy = Joule heat integral J^2 / sigma
//   free spin  - the tennis-racket box: energy and the vector L constant (time and rotation symmetry)
// With RF_PLOT_DIR=<folder> every test also writes its curves as CSV there for tools/plot_action.py
// (the plots of docs/11-action.md); without it nothing is written.
#include "TestRunner.h"
#include "Tests.h"

#include "plasma/MagneticField.h"

#include <cmath>
#include <initializer_list>

namespace {

// One CSV of a plot, written only when RF_PLOT_DIR is set.
class PlotCsv {
public:
    PlotCsv(const char* name, const char* header) {
        const char* dir = std::getenv("RF_PLOT_DIR");
        if (!dir || !*dir) return;
        const std::string path = std::string(dir) + "/" + name + ".csv";
        file_ = std::fopen(path.c_str(), "w");
        if (file_) std::fprintf(file_, "%s\n", header);
    }
    ~PlotCsv() {
        if (file_) std::fclose(file_);
    }
    PlotCsv(const PlotCsv&) = delete;
    PlotCsv& operator=(const PlotCsv&) = delete;
    void row(std::initializer_list<double> values) {
        if (!file_) return;
        bool first = true;
        for (double v : values) {
            std::fprintf(file_, first ? "%.9g" : ",%.9g", v);
            first = false;
        }
        std::fprintf(file_, "\n");
    }

private:
    FILE* file_ = nullptr;
};

// The inertia tensor of a body in its principal frame (diagonal), from the inverse the body keeps.
Vector3 principalInertia(const RigidBody& b) {
    return {1.0f / b.invInertiaLocal.x, 1.0f / b.invInertiaLocal.y, 1.0f / b.invInertiaLocal.z};
}

// omega . I omega: twice the rotational energy (the body frame, where I is diagonal).
double spinSquare(const RigidBody& b) {
    const Vector3 w = b.rotation().transposed() * b.angVel;
    const Vector3 I = principalInertia(b);
    return double(I.x) * w.x * w.x + double(I.y) * w.y * w.y + double(I.z) * w.z * w.z;
}

// Kinetic plus potential energy in gravity g (heights from y = 0).
double mechanicalEnergy(const RigidBody& b, float g) {
    return 0.5 * b.mass * double(length2(b.vel)) + 0.5 * spinSquare(b) + double(b.mass) * g * b.pos.y;
}

// A rigid world without the global damping and without sleep: only what a test adds dissipates.
void quietWorld(RigidWorld& w, const AABB& domain) {
    w.setDomain(domain);
    w.params.linearDamping = w.params.angularDamping = 0;
    w.params.rollingResistance = 0;
    w.params.sleeping = false;
}

} // namespace

// Friction. A slab 0.4 x 0.1 x 0.4 m (32 kg) slides on the floor at 3 m/s and stops. Coulomb's
// friction is the Rayleigh-type force F = -mu N v/|v| (not linear in v, but it still only takes
// energy): the energy lost must equal mu m g times the distance slid, and the stop distance is
// v0^2 / (2 mu g) = 0.765 m. Tolerance 2 %.
void testActionFriction() {
    RigidWorld w;
    quietWorld(w, AABB({-5, 0, -5}, {5, 2, 5}));
    const int s = w.addBox({-3, 0.05f, 0}, {0.2f, 0.05f, 0.2f}, Quaternion(), 1000, Vector3(1));
    RigidBody& b = w.bodies()[s];
    b.friction = 0.6f; // with the floor's 0.6: mu = sqrt(0.6 * 0.6) = 0.6 (ContactSolver's mix)
    b.staticFriction = 0.8f;
    const float dt = 1.0f / 600, g = -w.params.gravity.y, mu = 0.6f, v0 = 3.0f;
    for (int k = 0; k < 300; ++k) w.step(dt); // settle on the floor
    b.vel = {v0, 0, 0};
    const double E0 = mechanicalEnergy(b, g);
    const Vector3 start = b.pos;
    double work = 0; // mu m g * distance, summed step by step
    PlotCsv csv("friction", "t_s,energy_lost_J,friction_work_J");
    csv.row({0, 0, 0});
    for (int k = 1; k <= 1200; ++k) {
        const Vector3 before = b.pos;
        w.step(dt);
        const Vector3 d = b.pos - before;
        work += double(mu) * b.mass * g * std::sqrt(double(d.x) * d.x + double(d.z) * d.z);
        if (k % 6 == 0) csv.row({k * dt, E0 - mechanicalEnergy(b, g), work});
    }
    const double lost = E0 - mechanicalEnergy(b, g), slid = length(b.pos - start);
    const double stop = double(v0) * v0 / (2.0 * mu * g), balance = std::fabs(lost - work) / lost;
    std::printf("  friction: lost %.3f J, mu m g s = %.3f J (balance %.2f%%); slid %.4f m, v0^2 / 2 mu g = %.4f m\n", lost,
                work, 100 * balance, slid, stop);
    CHECK(balance < 0.02, "energy lost %.4f J vs friction work %.4f J", lost, work);
    CHECK(std::fabs(slid / stop - 1) < 0.02, "slid %.4f m, theory %.4f m", slid, stop);
}

// Linear damping as a Rayleigh function. The solver's global damping v <- v (1 - c dt) is the
// force -c m v (and the torque -c I w): R = c (m v^2 + w.Iw) / 2, dE/dt = -2R. A 1 m pendulum on a
// ball joint, released at 30 degrees, damped at c = 0.3 1/s for 8 s. The joint's impulse is a
// Lagrange multiplier: it holds the bob on the sphere and does no work, so the whole loss must be
// the damping's. Tolerance 2 %.
void testActionDampedPendulum() {
    RigidWorld w;
    quietWorld(w, AABB({-5, 0, -5}, {5, 10, 5}));
    const float c = 0.3f, angle = 30 * kPi / 180, dt = 1.0f / 600;
    w.params.linearDamping = w.params.angularDamping = c;
    const Vector3 pivot(0, 5, 0);
    const int bob = w.addSphere(pivot + Vector3(std::sin(angle), -std::cos(angle), 0), 0.05f, 5000, Vector3(1));
    w.addBallJoint(bob, -1, pivot);
    RigidBody& b = w.bodies()[bob];
    w.joints().back()->localAnchorA = b.rot.conjugate().rotate(pivot - b.pos);
    const float g = -w.params.gravity.y;
    const double E0 = mechanicalEnergy(b, g);
    double dissipated = 0, worstDrift = 0;
    PlotCsv csv("damped_pendulum", "t_s,energy_lost_J,damping_work_J");
    csv.row({0, 0, 0});
    for (int k = 1; k <= 8 * 600; ++k) {
        const double before = 2 * c * (0.5 * b.mass * double(length2(b.vel)) + 0.5 * spinSquare(b)); // 2R
        w.step(dt);
        const double after = 2 * c * (0.5 * b.mass * double(length2(b.vel)) + 0.5 * spinSquare(b));
        dissipated += 0.5 * (before + after) * dt; // trapezoid
        worstDrift = std::max(worstDrift, double(length(w.joints().back()->worldAnchorA(w.bodies()) - pivot)));
        if (k % 20 == 0) csv.row({k * dt, E0 - mechanicalEnergy(b, g), dissipated});
    }
    const double lost = E0 - mechanicalEnergy(b, g), balance = std::fabs(lost - dissipated) / lost;
    std::printf("  damped pendulum: E0 %.4f J, lost %.4f J, integral 2R dt = %.4f J (balance %.2f%%), joint drift %.1e m\n",
                E0 - double(b.mass) * g * (pivot.y - 1), lost, dissipated, 100 * balance, worstDrift);
    CHECK(balance < 0.02, "energy lost %.5f J vs damping work %.5f J", lost, dissipated);
}

namespace {

// Kinetic energy 1/2 rho |u|^2 and viscous dissipation 2 nu rho S:S of the gas, summed over the
// cells from the cell-centred velocity (central differences inside, one-sided at the walls).
struct GasLedger {
    double kinetic = 0, dissipationRate = 0;
};

GasLedger gasLedger(const GasSolver& g) {
    const int nx = g.nx(), ny = g.ny(), nz = g.nz();
    const double dx = g.dx(), vol = dx * dx * dx, rho = g.params.fluidDensity, nu = g.params.kinematicViscosity;
    auto diff = [&](int i, int j, int k, int axis) { // d u / d x_axis at a cell
        int lo[3] = {i, j, k}, hi[3] = {i, j, k};
        const int n[3] = {nx, ny, nz};
        lo[axis] = std::max(0, lo[axis] - 1);
        hi[axis] = std::min(n[axis] - 1, hi[axis] + 1);
        if (hi[axis] == lo[axis]) return Vector3(0.0f);
        return (g.cellVelocity(hi[0], hi[1], hi[2]) - g.cellVelocity(lo[0], lo[1], lo[2])) / float((hi[axis] - lo[axis]) * dx);
    };
    GasLedger L;
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                L.kinetic += 0.5 * rho * double(length2(g.cellVelocity(i, j, k))) * vol;
                const Vector3 grad[3] = {diff(i, j, k, 0), diff(i, j, k, 1), diff(i, j, k, 2)}; // grad[a][b] = d u_b / d x_a
                double SS = 0;
                for (int a = 0; a < 3; ++a)
                    for (int c = 0; c < 3; ++c) {
                        const double s = 0.5 * (double(grad[a][c]) + double(grad[c][a]));
                        SS += s * s;
                    }
                L.dissipationRate += 2 * nu * rho * SS * vol;
            }
    return L;
}

} // namespace

// Viscosity. The Taylor-Green vortex u = sin(pi x) cos(pi y), v = -cos(pi x) sin(pi y) in the unit
// box with free-slip walls, nu = 0.01 m^2/s, 32 cells, Courant 0.5, 1 s: its energy decays as
// exp(-4 nu pi^2 t) (Taylor & Green 1937), all of it through the Rayleigh function of viscosity,
// integral 2 nu rho |S|^2 dV. The test checks that ledger against the exact decay (5 %), and that
// the scheme never makes energy (loss >= 98 % of the viscous part). Whatever the solver loses
// beyond the viscous part is the numerical dissipation of the advection: printed and plotted.
void testActionViscousBalance() {
    GasSolver g;
    const int cells = 32;
    const float dx = 1.0f / cells, U = 1.0f;
    g.params.domainSize = {1, 1, 4 * dx};
    g.params.resolutionX = cells;
    g.params.inflowSpeed = 0;
    g.params.smokeRake = false;
    g.params.kinematicViscosity = 0.01f;
    g.params.wallFriction = false;
    for (auto& bc : g.params.bc) bc = BoundaryType::Wall;
    g.reset({0, 0, 0}, nullptr);
    g.setVelocity([&](const Vector3& x) {
        return Vector3(U * std::sin(kPi * x.x) * std::cos(kPi * x.y), -U * std::cos(kPi * x.x) * std::sin(kPi * x.y), 0);
    });
    const double nu = g.params.kinematicViscosity, rate = 4 * nu * kPi * kPi;
    GasLedger L = gasLedger(g);
    const double K0 = L.kinetic;
    double t = 0, dissipated = 0, prevRate = L.dissipationRate;
    PlotCsv csv("viscous", "t_s,energy_lost_J,viscous_work_J,exact_loss_J");
    csv.row({0, 0, 0, 0});
    while (t < 1.0 - 1e-9) {
        const double h = g.step(float(std::min(0.5 * dx / U, 1.0 - t)));
        t += h;
        L = gasLedger(g);
        dissipated += 0.5 * (prevRate + L.dissipationRate) * h; // trapezoid
        prevRate = L.dissipationRate;
        csv.row({t, K0 - L.kinetic, dissipated, K0 * (1 - std::exp(-rate * t))});
    }
    const double lost = K0 - L.kinetic, exact = K0 * (1 - std::exp(-rate * t));
    std::printf("  viscous: K0 %.5f J; after %.2f s lost %.5f J, integral 2 nu rho |S|^2 = %.5f J, exact loss %.5f J; "
                "viscous ledger vs exact %+.1f%%, numerical dissipation %.1f%% of the loss\n",
                K0, t, lost, dissipated, exact, 100 * (dissipated / exact - 1), 100 * (1 - dissipated / lost));
    CHECK(std::fabs(dissipated / exact - 1) < 0.05, "viscous dissipation %.5f J vs exact %.5f J", dissipated, exact);
    CHECK(lost >= 0.98 * dissipated, "the scheme made energy: lost %.5f J < viscous %.5f J", lost, dissipated);
}

namespace {

// The ledger of the constrained-transport scheme itself: the magnetic energy of every face,
// sum B^2 dx^3 / 2mu0 (the faces in the walls are frozen: they add a constant), and the Joule
// power of every FREE edge, sum mu0 eta J^2 dx^3, with J = curl B / mu0 on the edge as
// MagneticField::computeCurrent takes it. The edges on a perfectly conducting wall carry E = 0:
// no work is done there. Summation by parts makes the two exact partners: dW/dt = -sum E.J dV.
double faceEnergy(const MagneticField& m, float dx) {
    double e = 0;
    for (const Field3* f : {&m.bx, &m.by, &m.bz})
        for (float b : f->d) e += double(b) * b;
    return e * double(dx) * dx * dx / (2.0 * MagneticField::kMu0);
}

double freeEdgeJoulePower(const MagneticField& m, int nx, int ny, int nz, float dx) {
    const double s = 1.0 / (double(dx) * MagneticField::kMu0);
    const Field3 &bx = m.bx, &by = m.by, &bz = m.bz;
    double J2 = 0;
    for (int k = 0; k <= nz; ++k)
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i <= nx; ++i) {
                const bool iIn = i > 0 && i < nx, jIn = j > 0 && j < ny, kIn = k > 0 && k < nz;
                if (i < nx && jIn && kIn) J2 += sqr(s * ((bz.at(i, j, k) - bz.at(i, j - 1, k)) - (by.at(i, j, k) - by.at(i, j, k - 1))));
                if (j < ny && iIn && kIn) J2 += sqr(s * ((bx.at(i, j, k) - bx.at(i, j, k - 1)) - (bz.at(i, j, k) - bz.at(i - 1, j, k))));
                if (k < nz && iIn && jIn) J2 += sqr(s * ((by.at(i, j, k) - by.at(i - 1, j, k)) - (bx.at(i, j, k) - bx.at(i, j - 1, k))));
            }
    return MagneticField::kMu0 * m.resistivity() * J2 * double(dx) * dx * dx;
}

} // namespace

// Resistance. A field Bz = B1 cos(pi x / L) in a resting conductor (sigma = 1e5 S/m, eta = 7.96
// m^2/s) between perfectly conducting walls: E = eta J vanishes on the walls, no Poynting flux
// leaves, so the magnetic energy integral B^2 / 2mu0 can only turn into Joule heat integral
// J^2 / sigma (the Rayleigh function of Ohm's law). 30 steps of 0.1 ms. Tolerance 1 %, on the
// scheme's own ledger (above). Checked beside it: MagneticField::jouleHeating, the heat the gas is
// given. It once averaged the edge currents into cells, the wall edges included, and so counted
// current where E = 0 does no work (+3.2 %); now it heats from the free edges, as the ledger does.
void testActionJouleBalance() {
    const int nx = 32, ny = 4, nz = 32;
    const float dx = 1.0f / nx, dt = 1e-4f;
    MagneticField m;
    m.conductivity = 1e5f;
    m.numericalDissipation = 0;
    m.reset(nx, ny, nz, dx, Vector3(0.0f));
    for (int k = 0; k <= nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) m.bz.at(i, j, k) = 1e-3f * std::cos(kPi * (i + 0.5f) * dx);
    Field3 u, v, w;
    u.init(nx + 1, ny, nz, {0, 0.5f, 0.5f});
    v.init(nx, ny + 1, nz, {0.5f, 0, 0.5f});
    w.init(nx, ny, nz + 1, {0.5f, 0.5f, 0});
    std::vector<float> heat;
    auto joulePower = [&] { // [W]: the heat of a unit time
        m.jouleHeating(heat, 1.0f);
        double sum = 0;
        for (float q : heat) sum += q;
        return sum * double(dx) * dx * dx;
    };
    m.induce(u, v, w, 1.0f, 0.0f); // computes J of the start field, B unchanged
    const double W0 = faceEnergy(m, dx);
    double joule = 0, power = freeEdgeJoulePower(m, nx, ny, nz, dx);
    double diagnostic = 0, diagPower = joulePower();
    PlotCsv csv("joule", "t_ms,magnetic_energy_lost_uJ,joule_heat_uJ,cell_diagnostic_uJ");
    csv.row({0, 0, 0, 0});
    for (int s = 1; s <= 30; ++s) {
        m.induce(u, v, w, 1.0f, dt);
        const double next = freeEdgeJoulePower(m, nx, ny, nz, dx), nextDiag = joulePower();
        joule += 0.5 * (power + next) * dt; // trapezoid
        diagnostic += 0.5 * (diagPower + nextDiag) * dt;
        power = next;
        diagPower = nextDiag;
        csv.row({s * dt * 1e3, (W0 - faceEnergy(m, dx)) * 1e6, joule * 1e6, diagnostic * 1e6});
    }
    const double lost = W0 - faceEnergy(m, dx), balance = std::fabs(lost - joule) / lost;
    std::printf("  Joule: W0 %.4e J, lost %.4e J (%.1f%%), integral J^2 / sigma over the free edges %.4e J (balance %.3f%%); "
                "jouleHeating (the heat the gas gets) %.4e J (%+.2f%%)\n",
                W0, lost, 100 * lost / W0, joule, 100 * balance, diagnostic, 100 * (diagnostic / lost - 1));
    CHECK(balance < 0.01, "magnetic energy lost %e J vs Joule heat %e J", lost, joule);
    CHECK(std::fabs(diagnostic / lost - 1) < 0.01, "the gas got %e J of Joule heat for %e J of field energy", diagnostic, lost);
}

namespace {

// A free box spun about its middle axis (the tennis-racket / Dzhanibekov flip): the largest
// relative drift of its energy and of the vector L (world) over `seconds` at the step dt.
struct SpinDrift {
    double energy = 0, momentum = 0;
    int flips = 0;
};

SpinDrift freeSpin(float dt, float seconds, PlotCsv* csv) {
    RigidWorld w;
    quietWorld(w, AABB({-10, -10, -10}, {10, 10, 10}));
    w.params.gravity = Vector3(0.0f);
    w.params.collideWithDomain = false;
    const int id = w.addBox({0, 0, 0}, {0.05f, 0.1f, 0.2f}, Quaternion(), 1000, Vector3(1)); // I: x > y > z
    RigidBody& b = w.bodies()[id];
    b.angVel = {0.01f, 5.0f, 0.01f}; // the middle axis, slightly off
    auto angularMomentum = [&] { return b.rotation() * (principalInertia(b) * (b.rotation().transposed() * b.angVel)); };
    const double E0 = 0.5 * spinSquare(b);
    const Vector3 L0 = angularMomentum();
    SpinDrift d;
    float prevSign = 1;
    const int steps = int(std::lround(seconds / dt)), every = std::max(1, steps / 600);
    for (int k = 1; k <= steps; ++k) {
        w.step(dt);
        d.energy = std::max(d.energy, std::fabs(0.5 * spinSquare(b) / E0 - 1));
        d.momentum = std::max(d.momentum, double(length(angularMomentum() - L0) / length(L0)));
        const float wy = (b.rotation().transposed() * b.angVel).y;
        if (wy * prevSign < 0) {
            ++d.flips;
            prevSign = wy;
        }
        if (csv && k % every == 0) {
            const Vector3 wb = b.rotation().transposed() * b.angVel;
            csv->row({k * dt, 0.5 * spinSquare(b) / E0 - 1, length(angularMomentum() - L0) / length(L0), wb.x, wb.y, wb.z});
        }
    }
    return d;
}

} // namespace

// Free rotation. No force, no torque: the action is unchanged by a shift in time (energy kept) and
// by a turn of the whole world (the vector L kept) - Noether. A box spun about its middle axis
// flips over and over while both stay put. The free turn is a splitting into exact turns about the
// principal axes (FreeRotation.cpp): every sub-step is a rotation, so the vector L is kept to
// rounding at any step, and the energy error is that of a second-order symplectic method - it
// shrinks four times when the step halves, and stays bounded, with no steps at the flips.
// (The implicit midpoint turn it replaced drifted L at first order - 2.8e-3 over 20 s - and the
// energy rose in steps at every flip, more with a smaller step: 3e-3 over 20 s.)
void testActionFreeRotation() {
    SpinDrift d;
    {
        PlotCsv csv("free_rotation", "t_s,energy_drift,momentum_drift,wx,wy,wz");
        d = freeSpin(1.0f / 600, 20, &csv);
    }
    std::printf("  free spin 20 s at 1/600 s: %d flips, energy drift %.1e, |L - L0| / |L0| %.1e (target 1e-3%s)\n", d.flips, d.energy,
                d.momentum, d.energy < 1e-3 && d.momentum < 1e-3 ? ": met" : ": NOT met, see docs/11-action.md");
    CHECK(d.flips >= 4, "the box must flip over (it flipped %d times)", d.flips);
    CHECK(d.energy < 1e-4 && d.momentum < 1e-5, "energy drift %e, momentum drift %e", d.energy, d.momentum);
    // The order of the energy error, where the truncation - not the rounding - is what is measured:
    // steps of 1/60, 1/120 and 1/240 s over 5 s. The vector L at every step: rounding only.
    PlotCsv conv("free_rotation_convergence", "dt_s,momentum_drift,energy_drift");
    const float steps[3] = {1.0f / 60, 1.0f / 120, 1.0f / 240};
    double prevEnergy = 0;
    std::printf("  drift over 5 s:");
    for (float dt : steps) {
        const SpinDrift s = freeSpin(dt, 5, nullptr);
        conv.row({dt, s.momentum, s.energy});
        std::printf("  dt 1/%.0f: L %.1e, E %.1e", 1 / dt, s.momentum, s.energy);
        if (prevEnergy > 0) {
            const double order = std::log2(prevEnergy / std::max(s.energy, 1e-30));
            std::printf(" (energy order %.2f)", order);
            CHECK(order > 1.8, "the energy error must fall as dt^2: order %f", order);
        }
        CHECK(s.momentum < 1e-5, "the vector L drifted %e at dt 1/%.0f", s.momentum, 1 / dt);
        prevEnergy = s.energy;
    }
    std::printf("\n");
}
