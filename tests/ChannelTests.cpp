// The channels (scene/Channels.h) checked against physics, not against themselves: a chart is only
// worth what the number under it is worth. Each test lets a law of nature read its own channel:
//   free fall   - energy kept up to the known offset of the integrator, momentum m g t
//   collision   - momentum and energy of an elastic collision, the group's centre of mass moving
//                 uniformly
//   friction    - the lost energy equals mu m g x distance (docs/11-action.md)
//   stack       - Newton's third law: the ground feels the whole stack, each box its own weight
//   Taylor-Green - the gas's kinetic energy decays as exp(-4 nu pi^2 t)
//   Poiseuille  - a probe line across a channel reads the parabola 6 U eta (1 - eta)
//   hydrostatics - a probe line down a heavy gas reads p = rho g h
//   cost        - nothing asked: nothing measured, no allocation; measuring changes no bit of the
//                 physics, and the numbers are the same bits on 1 thread and on all
#include "TestRunner.h"
#include "Tests.h"

#include "core/Parallel.h"

#include <cmath>
#include <cstring>
#include <limits>

namespace {

// The value of channel `id` among the measurements (NaN: not measured).
double valueOf(const std::vector<Measurement>& ms, const std::string& id) {
    for (const Measurement& m : ms)
        if (m.info.id == id) return m.value;
    return std::numeric_limits<double>::quiet_NaN();
}

// A rigid arena without the global damping and without sleep: only what a test adds dissipates.
void quietArena(Simulation& sim, const AABB& box) {
    sim.useRigidArena(box);
    sim.rigid.params.linearDamping = sim.rigid.params.angularDamping = 0;
    sim.rigid.params.rollingResistance = 0;
    sim.rigid.params.sleeping = false;
}

// FNV-1a over the raw bits of a list of vectors (tests/DeterminismTests.cpp does the same).
void hashVectors(uint64_t& h, const Vector3* v, size_t n) {
    const unsigned char* bytes = reinterpret_cast<const unsigned char*>(v);
    for (size_t k = 0; k < n * sizeof(Vector3); ++k) {
        h ^= bytes[k];
        h *= 1099511628211ull;
    }
}

uint64_t stateHash(const Simulation& sim) {
    uint64_t h = 14695981039346656037ull;
    for (const RigidBody& b : sim.rigid.bodies()) {
        const Vector3 v[3] = {b.pos, b.vel, b.angVel};
        hashVectors(h, v, 3);
    }
    hashVectors(h, sim.particles.positions().data(), sim.particles.size());
    hashVectors(h, sim.particles.velocities().data(), sim.particles.size());
    return h;
}

} // namespace

// Free fall. A ball dropped from 5 m falls 0.5 s in an empty arena. The integrator (symplectic
// Euler: v += g h, then x += v h) puts the ball a known bit lower than the exact parabola,
// x_n = x0 - g t^2 / 2 - g h t / 2, so the mechanical energy falls by exactly m g^2 h t / 2 while
// the momentum grows exactly as m g t. The channels must read both to 1e-4.
void testChannelsFreeFall() {
    Simulation sim;
    quietArena(sim, AABB({-2, 0, -2}, {2, 6, 2}));
    const int ball = sim.rigid.addSphere({0, 5, 0}, 0.1f, 1000, Vector3(1));
    sim.channels.scene = true;
    sim.channels.objects.push_back({ObjectRef::Kind::Body, ball, "ball"});
    const double E0 = valueOf(measureScene(sim), "rigid/mechanical energy");
    for (int f = 0; f < 30; ++f) sim.stepFrame();
    const double m = sim.rigid.bodies()[size_t(ball)].mass, g = -sim.gravity().y, t = sim.time();
    const double h = sim.rigid.lastStepDt(), predicted = m * g * g * h * t / 2;
    const double E = valueOf(sim.measurements(), "rigid/mechanical energy");
    const double P = valueOf(sim.measurements(), "rigid/momentum"), speed = valueOf(sim.measurements(), "ball/speed");
    std::printf("  free fall %.2f s: energy %.4f -> %.4f J, drop %.5f J vs integrator's m g^2 h t / 2 = %.5f J; "
                "|P| %.5f vs m g t %.5f N s; ball speed %.5f vs g t %.5f m/s\n",
                t, E0, E, E0 - E, predicted, P, m * g * t, speed, g * t);
    CHECK(std::fabs((E0 - E) - predicted) < 1e-4 * E0, "energy drop %f J, the integrator predicts %f J", E0 - E, predicted);
    CHECK(std::fabs(P - m * g * t) < 1e-4 * m * g * t, "momentum %f vs m g t %f", P, m * g * t);
    CHECK(std::fabs(speed - g * t) < 1e-4 * g * t, "speed %f vs g t %f", speed, g * t);
}

// An elastic collision in weightless space: a steel ball at 2 m/s hits one at rest (e = 1, no
// friction). The scene's momentum and energy stay, the struck ball leaves with 2 m/s, and the
// pair's centre of mass (a group) moves on uniformly at P / M through the collision.
void testChannelsCollision() {
    Simulation sim;
    quietArena(sim, AABB({-2, 0, -2}, {2, 4, 2}));
    sim.setGravity({0, 0, 0});
    const int a = sim.rigid.addSphere({-0.5f, 2, 0}, 0.1f, 7800, Vector3(1));
    const int b = sim.rigid.addSphere({0.5f, 2, 0}, 0.1f, 7800, Vector3(1));
    for (int i : {a, b}) {
        RigidBody& body = sim.rigid.bodies()[size_t(i)];
        body.restitution = 1, body.friction = body.staticFriction = 0;
    }
    sim.rigid.bodies()[size_t(a)].vel = {2, 0, 0};
    sim.channels.scene = true;
    sim.channels.objects.push_back({ObjectRef::Kind::Body, b, "struck"});
    sim.channels.groups.push_back({"pair", {{ObjectRef::Kind::Body, a, "a"}, {ObjectRef::Kind::Body, b, "b"}}});
    const std::vector<Measurement> before = measureAll(sim, sim.channels);
    for (int f = 0; f < 60; ++f) sim.stepFrame();
    const std::vector<Measurement>& after = sim.measurements();
    const double P0 = valueOf(before, "rigid/momentum x"), P1 = valueOf(after, "rigid/momentum x");
    const double E0 = valueOf(before, "rigid/kinetic energy"), E1 = valueOf(after, "rigid/kinetic energy");
    const double x0 = valueOf(before, "pair/centre x"), x1 = valueOf(after, "pair/centre x");
    const double drift = P0 / valueOf(before, "pair/mass") * sim.time();
    std::printf("  collision: momentum %.5f -> %.5f N s, energy %.5f -> %.5f J, struck ball %.4f m/s, "
                "pair's centre moved %.5f m (P/M t = %.5f m), group momentum %.5f\n",
                P0, P1, E0, E1, valueOf(after, "struck/speed"), x1 - x0, drift, valueOf(after, "pair/momentum"));
    CHECK(std::fabs(P1 / P0 - 1) < 1e-3, "momentum %f -> %f", P0, P1);
    CHECK(std::fabs(E1 / E0 - 1) < 1e-3, "energy %f -> %f", E0, E1);
    CHECK(std::fabs(valueOf(after, "struck/speed") - 2) < 0.02, "the struck ball leaves with %f m/s", valueOf(after, "struck/speed"));
    CHECK(std::fabs((x1 - x0) - drift) < 1e-3, "the pair's centre moved %f m, P/M t = %f m", x1 - x0, drift);
    CHECK(std::fabs(valueOf(after, "pair/momentum") - P1) < 1e-6 * P0, "the group's momentum differs from the scene's");
}

// Friction as the channels see it: a slab slides on the floor at 3 m/s and stops (the setup of
// testActionFriction). Its lost mechanical energy, read from its own channel, equals the friction
// work mu m g x the distance its centre slid (read from its channels too). Tolerance 2 %.
void testChannelsFriction() {
    Simulation sim;
    quietArena(sim, AABB({-5, 0, -5}, {5, 2, 5}));
    const int s = sim.rigid.addBox({-3, 0.05f, 0}, {0.2f, 0.05f, 0.2f}, Quaternion(), 1000, Vector3(1));
    RigidBody& slab = sim.rigid.bodies()[size_t(s)];
    slab.friction = 0.6f, slab.staticFriction = 0.8f;
    sim.channels.objects.push_back({ObjectRef::Kind::Body, s, "slab"});
    for (int f = 0; f < 30; ++f) sim.stepFrame(); // settle on the floor
    slab.vel = {3, 0, 0};
    std::vector<Measurement> now = measureObject(sim, sim.channels.objects[0]);
    const double E0 = valueOf(now, "slab/energy"), m = valueOf(now, "slab/mass"), g = -sim.gravity().y, mu = 0.6;
    double work = 0, x = valueOf(now, "slab/centre x"), z = valueOf(now, "slab/centre z");
    for (int f = 0; f < 90; ++f) {
        sim.stepFrame();
        const double nx = valueOf(sim.measurements(), "slab/centre x"), nz = valueOf(sim.measurements(), "slab/centre z");
        work += mu * m * g * std::hypot(nx - x, nz - z);
        x = nx, z = nz;
    }
    const double lost = E0 - valueOf(sim.measurements(), "slab/energy");
    std::printf("  friction: lost %.3f J (its energy channel), mu m g * distance = %.3f J (its centre channels), balance %.2f%%\n",
                lost, work, 100 * std::fabs(lost - work) / lost);
    CHECK(std::fabs(lost - work) < 0.02 * lost, "energy lost %f J vs friction work %f J", lost, work);
}

// Newton's third law on a chart. Four 4 kg boxes stand on a static slab (the ground, a body of its
// own). Each box's contacts hold it up with its own weight, m g, upwards; the ground's contacts
// carry the whole stack, 4 m g, downwards. Read from the contact-force channels, within 2 %.
void testChannelsStackForces() {
    Simulation sim;
    quietArena(sim, AABB({-2, 0, -2}, {2, 3, 2}));
    const int ground = sim.rigid.addBox({0, 0.1f, 0}, {1, 0.1f, 1}, Quaternion(), 0, Vector3(1));
    std::vector<int> boxes;
    for (int k = 0; k < 4; ++k) boxes.push_back(sim.rigid.addBox({0, 0.3f + 0.2f * k, 0}, {0.1f, 0.1f, 0.1f}, Quaternion(), 500, Vector3(1)));
    sim.channels.objects.push_back({ObjectRef::Kind::Body, ground, "ground"});
    sim.channels.objects.push_back({ObjectRef::Kind::Body, boxes.front(), "bottom"});
    sim.channels.objects.push_back({ObjectRef::Kind::Body, boxes.back(), "top"});
    for (int f = 0; f < 120; ++f) sim.stepFrame();
    const std::vector<Measurement>& ms = sim.measurements();
    const double w = sim.rigid.bodies()[size_t(boxes[0])].mass * -sim.gravity().y;
    const double ground_ = valueOf(ms, "ground/contact force y"), bottom = valueOf(ms, "bottom/contact force y"),
                 top = valueOf(ms, "top/contact force y");
    std::printf("  stack of 4 x %.2f N: ground feels %.2f N (-4 m g = %.2f), bottom box %.2f N, top box %.2f N (m g = %.2f); "
                "ground touches %g points\n",
                w, ground_, -4 * w, bottom, top, w, valueOf(ms, "ground/contacts"));
    CHECK(std::fabs(ground_ / (-4 * w) - 1) < 0.02, "the ground feels %f N, the stack weighs %f N", ground_, 4 * w);
    CHECK(std::fabs(bottom / w - 1) < 0.02, "the bottom box is held with %f N, its weight is %f N", bottom, w);
    CHECK(std::fabs(top / w - 1) < 0.02, "the top box is held with %f N, its weight is %f N", top, w);
    // Asleep, the solver skips the stack: the boxes, at rest, report the balance of their weight,
    // exactly m g (see contactForce in ObjectChannels.cpp); the ground has no solved contact left.
    sim.rigid.params.sleeping = true;
    for (int f = 0; f < 120; ++f) sim.stepFrame();
    const double asleep = valueOf(sim.measurements(), "bottom/asleep"), bottomAsleep = valueOf(sim.measurements(), "bottom/contact force y");
    std::printf("  asleep (%g): bottom box %.3f N, top box %.3f N (m g = %.3f); the ground %.2f N (no solved contacts)\n", asleep,
                bottomAsleep, valueOf(sim.measurements(), "top/contact force y"), w, valueOf(sim.measurements(), "ground/contact force y"));
    CHECK(asleep == 1 && std::fabs(bottomAsleep / w - 1) < 1e-5, "a sleeping box reads %f N, its weight is %f N", bottomAsleep, w);
}

namespace {

// A gas box of `cells` cells across [0, 1] x [0, 1] x [0, 4 dx], closed on every side, still.
void closedGasBox(Simulation& sim, int cells, float viscosity) {
    GasParams& p = sim.grid.params;
    p.domainSize = {1, 1, 4.0f / cells};
    p.resolutionX = cells;
    p.inflowSpeed = 0;
    p.smokeRake = false;
    p.kinematicViscosity = viscosity;
    p.wallFriction = false;
    for (auto& bc : p.bc) bc = BoundaryType::Wall;
    sim.useGasBox({0, 0, 0});
}

} // namespace

// The Taylor-Green vortex u = sin(pi x) cos(pi y), v = -cos(pi x) sin(pi y) (Taylor & Green 1937)
// in a closed unit box: its kinetic energy decays as exp(-4 nu pi^2 t). The gas channel reads it;
// the scheme's own numerical dissipation adds a little loss on top (testActionViscousBalance
// measures it: 1.6 % of the loss with advection-reflection, 13.8 % without), so the tolerance is
// 8 % of the energy left - and the channel must never read more than the exact value (+1 %).
void testChannelsTaylorGreen() {
    Simulation sim;
    const int cells = 32;
    closedGasBox(sim, cells, 0.01f);
    sim.grid.setVelocity([](const Vector3& x) {
        return Vector3(std::sin(kPi * x.x) * std::cos(kPi * x.y), -std::cos(kPi * x.x) * std::sin(kPi * x.y), 0);
    });
    const double K0 = valueOf(measureScene(sim), "gas/kinetic energy"), nu = sim.grid.params.kinematicViscosity;
    double t = 0;
    while (t < 1.0 - 1e-9) t += sim.grid.step(float(std::min(0.5 / cells, 1.0 - t)));
    const double K = valueOf(measureScene(sim), "gas/kinetic energy"), exact = K0 * std::exp(-4 * nu * kPi * kPi * t);
    std::printf("  Taylor-Green: gas kinetic energy %.5f -> %.5f J after %.2f s, exact %.5f J (%+.2f%%)\n", K0, K, t, exact,
                100 * (K / exact - 1));
    CHECK(K0 > 0 && std::fabs(K / exact - 1) < 0.08, "gas energy %f J vs exact %f J", K, exact);
    CHECK(K < 1.01 * exact, "the gas gained energy: %f J > %f J", K, exact);
}

// Poiseuille flow between two plates H = 1 m apart (the setup of verification's gas-poiseuille:
// uniform inflow U, Re = 2, 16 cells across, developed after two viscous times). A probe line at
// x = 3 m through the centres of the cells across the channel must read the parabola
// u = 6 U eta (1 - eta) (White, "Viscous Fluid Flow", sec. 3-2) - the case reaches ~4e-3 U here.
void testChannelsPoiseuilleProbe() {
    Simulation sim;
    const int cells = 16;
    const float H = 1, U = 1, dx = H / cells, y0 = 2 * dx;
    GasParams& p = sim.grid.params;
    p.domainSize = {4 * H, H + 4 * dx, 4 * dx};
    p.resolutionX = 4 * cells;
    p.inflowSpeed = U;
    p.kinematicViscosity = 0.5f;
    p.smokeRake = false;
    p.wallFriction = false;
    sim.grid.vessel = [=](const Vector3& x) { return x.y > y0 && x.y < y0 + H; };
    sim.useGasBox({0, 0, 0});
    auto exact = [=](double y) { const double eta = (y - y0) / H; return 6 * U * eta * (1 - eta); };
    sim.grid.setVelocity([&](const Vector3& x) { return Vector3(x.y > y0 && x.y < y0 + H ? float(exact(x.y)) : 0.0f, 0, 0); });
    const double tEnd = 2 * H * H / p.kinematicViscosity;
    for (double t = 0; t < tEnd - 1e-9;) t += sim.grid.step(float(std::min<double>(0.5 * dx / (1.5 * U), tEnd - t)));
    const FieldProbe line{"across", {3, y0 + 0.5f * dx, 1.5f * dx}, {3, y0 + H - 0.5f * dx, 1.5f * dx}, cells, FieldQuantity::VelocityX};
    const std::vector<Measurement> ms = measureProbe(sim, line);
    double worst = 0;
    for (size_t k = 0; !ms.empty() && k < ms[0].profile.size(); ++k)
        worst = std::max(worst, std::fabs(ms[0].profile[k] - exact(y0 + 0.5 * dx + ms[0].distance[k])));
    std::printf("  Poiseuille probe: %zu points across, worst |u - 6 U eta (1 - eta)| = %.2e U, centre %.4f U (exact 1.5)\n",
                ms.empty() ? size_t(0) : ms[0].profile.size(), worst / U, ms.empty() ? 0.0 : ms[0].profile[cells / 2] / U);
    CHECK(!ms.empty() && ms[0].profile.size() == size_t(cells), "the probe line gave no profile");
    CHECK(worst < 0.01 * U, "the probe reads the parabola to %e U only", worst / U);
}

// Hydrostatics. A closed box of gas made heavy by smoke (the solver's smoke buoyancy, -9.81 m/s^2
// per unit smoke, is a gravity g on the whole gas): at rest the pressure must balance it,
// dp/dy = -rho g, so between the bottom and the top of a probe line p differs by rho g h. The
// pressure in a closed box is defined up to a constant: the difference is the physics. 2 %.
void testChannelsHydrostaticProbe() {
    Simulation sim;
    const int cells = 16;
    closedGasBox(sim, cells, 0);
    sim.grid.params.smokeBuoyancy = 9.81f;
    sim.grid.setTracer([](const Vector3&) { return 1.0f; });
    for (int s = 0; s < 5; ++s) sim.grid.step(0.01f);
    const float dx = 1.0f / cells;
    const FieldProbe line{"column", {0.5f + 0.5f * dx, 0.5f * dx, 2.5f * dx}, {0.5f + 0.5f * dx, 1 - 0.5f * dx, 2.5f * dx}, cells,
                          FieldQuantity::Pressure};
    const std::vector<Measurement> ms = measureProbe(sim, line);
    const double rho = sim.grid.params.fluidDensity, h = 1 - dx, expected = rho * 9.81 * h;
    const double dp = ms.empty() ? 0 : ms[0].profile.front() - ms[0].profile.back();
    const double speed = valueOf(measureScene(sim), "gas/max speed");
    std::printf("  hydrostatics: p(bottom) - p(top) = %.4f Pa over %.4f m, rho g h = %.4f Pa (%+.2f%%); gas speed %.1e m/s\n", dp, h,
                expected, 100 * (dp / expected - 1), speed);
    CHECK(std::fabs(dp / expected - 1) < 0.02, "pressure difference %f Pa vs rho g h %f Pa", dp, expected);
    CHECK(speed < 1e-3, "the heavy gas should stay at rest, it moves at %f m/s", speed);
}

// The cost and the innocence of measuring. Two identical scenes (twelve falling boxes and a block
// of liquid), one measured every frame (scene, an object, a group), one not: after 40 frames every
// bit of their state is equal - measuring only reads. Nothing asked: no channel appears and the
// measurement allocates nothing. And the numbers are the same bits on 1 thread and on all.
void testChannelsCostAndInnocence() {
    auto build = [](Simulation& sim) {
        sim.useRigidArena(AABB({-1, 0, -1}, {1, 2, 1}));
        for (int k = 0; k < 12; ++k)
            sim.rigid.addBox({-0.5f + 0.09f * k, 0.5f + 0.1f * k, 0.2f}, {0.04f, 0.04f, 0.04f}, Quaternion(), 500, Vector3(1));
        sim.particles.addBlock(AABB({-0.3f, 0.05f, -0.5f}, {0.1f, 0.25f, -0.2f}));
    };
    Simulation plain, watched;
    build(plain);
    build(watched);
    watched.channels.scene = true;
    watched.channels.objects.push_back({ObjectRef::Kind::Body, 3, "box3"});
    watched.channels.groups.push_back({"pair", {{ObjectRef::Kind::Body, 0, "b0"}, {ObjectRef::Kind::Body, 1, "b1"}}});
    for (int f = 0; f < 40; ++f) plain.stepFrame(), watched.stepFrame();
    const uint64_t a = stateHash(plain), b = stateHash(watched);
    Probe::clear();
    plain.stepFrame();
    const bool silent = Probe::snapshot().value("rigid/kinetic energy", -1.0) == -1.0 && plain.measurements().empty();
    const long long before = Probe::allocations.load();
    const std::vector<Measurement> none = measureAll(plain, ChannelRequests());
    const long long allocations = Probe::allocations.load() - before;
    ThreadPool::instance().setActiveThreads(1);
    const std::vector<Measurement> one = measureScene(watched);
    ThreadPool::instance().setActiveThreads(0);
    const std::vector<Measurement> all = measureScene(watched);
    bool sameBits = one.size() == all.size();
    for (size_t k = 0; sameBits && k < one.size(); ++k) sameBits = std::memcmp(&one[k].value, &all[k].value, sizeof(double)) == 0;
    std::printf("  cost: state hash %016llx vs %016llx (measured), %zu measurements a frame; nothing asked: silent %d, "
                "%lld allocations; 1 thread vs all: same bits %d\n",
                (unsigned long long)a, (unsigned long long)b, watched.measurements().size(), int(silent), allocations, int(sameBits));
    CHECK(a == b, "measuring changed the physics");
    CHECK(silent && none.empty() && allocations == 0, "nothing asked, yet something was measured or allocated");
    CHECK(sameBits, "the scene channels differ between 1 thread and all");
    CHECK(!watched.measurements().empty(), "the measured scene reported nothing");
}

// Every name card is complete: an ASCII id, a unit, a name.
void testChannelsCatalog() {
    bool complete = true;
    for (const ChannelInfo& c : sceneChannelCatalog()) {
        for (char ch : c.id) complete &= (unsigned char)ch < 128;
        complete &= !c.id.empty() && !c.unit.empty() && !c.name.empty() && describeChannel(c.id) == &c;
    }
    for (int q = 0; q <= int(FieldQuantity::CurrentDensity); ++q) {
        const ChannelInfo c = fieldQuantityInfo(FieldQuantity(q));
        complete &= !c.id.empty() && !c.unit.empty() && !c.name.empty();
    }
    std::printf("  catalog: %zu scene channels, %d probe quantities, all with id, unit and name: %d\n", sceneChannelCatalog().size(),
                int(FieldQuantity::CurrentDensity) + 1, int(complete));
    CHECK(complete, "a channel's name card is incomplete");
}
