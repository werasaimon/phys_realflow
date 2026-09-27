// Lessons: scenes built to show one law of physics and to put its number next to the theory on
// the panel. The Dzhanibekov effect (a handle spun about its middle axis flips over and over while
// its angular momentum and energy stay put), Newton's cradle (momentum and energy handed down a
// row of balls) and Galileo's tower (every mass falls alike - until the air is switched on).
// Each scene reports "what to check" through describe(): the panel shows the measured value and
// the predicted one side by side. The physics is in src/rigid; here is only the setup.
#include "samples/Samples.h"

#include "core/Format.h"
#include "core/Mesh.h"
#include "core/Probe.h"
#include "rigid/Shapes.h"

#include <cmath>
#include <memory>

namespace rf {

namespace {

// Angular momentum of a body about its centre, |L| = |I w| in the body's principal frame.
float angularMomentum(const RigidBody& b) {
    const Vector3 wl = b.rotation().transposed() * b.angVel;
    return length(Vector3(wl.x / b.invInertiaLocal.x, wl.y / b.invInertiaLocal.y, wl.z / b.invInertiaLocal.z));
}

// --- Dzhanibekov ------------------------------------------------------------------------------
// A free body spun about the axis of its middle moment of inertia is unstable: Euler's equations
// let a tiny wobble grow at the rate sigma = w sqrt((I2 - I1)(I3 - I2) / (I1 I3)) until the body
// flips over, then the wobble dies down, and it flips back - forever (Ashbaugh, Chicone & Cushman
// 1991, "The twisting tennis racket"). What stays constant through every flip: |L| and the
// rotational energy - the two Noether invariants of a free body. The scene is the T-handle of the
// 1985 Salyut-7 footage: a bar with a stub on its middle, made of two boxes as one compound.
class DzhanibekovScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.setGravity(Vector3(0.0f));               // in orbit
        sim.rigid.params.sleeping = false;           // it never comes to rest
        sim.rigid.params.linearDamping = sim.rigid.params.angularDamping = 0;
        sim.rigid.params.collideWithDomain = false;  // nothing to hit
        sim.vis.showSlice = sim.vis.showStreamlines = false;
    }
    void build(Simulation& sim) override {
        sim.useRigidArena(AABB({-2, -2, -2}, {2, 2, 2}));
        // The handle: a bar 0.3 x 0.04 x 0.04 m with a stub 0.04 x 0.15 x 0.04 m on its middle.
        std::vector<TriMesh> parts = {primitives::box({0.15f, 0.02f, 0.02f}), primitives::box({0.02f, 0.075f, 0.02f})};
        parts[1].translate({0.0f, 0.095f, 0.0f});
        auto shape = std::make_shared<CompoundShape>(parts, primitives::merge(parts));
        handle_ = sim.rigid.addCompound(shape, {0, 0, 0}, Quaternion(), 2700.0f, {0.85f, 0.85f, 0.9f}); // aluminium
        RigidBody& b = sim.rigid.bodies()[handle_];
        // The principal moments (the body frame is the principal frame) and the middle one.
        const Vector3 I(1 / b.invInertiaLocal.x, 1 / b.invInertiaLocal.y, 1 / b.invInertiaLocal.z);
        const int axes[3] = {0, 1, 2};
        int lo = 0, hi = 0;
        for (int k : axes) {
            if (I[k] < I[lo]) lo = k;
            if (I[k] > I[hi]) hi = k;
        }
        mid_ = 3 - lo - hi;
        // Spin about the middle axis at 6 rad/s with a 1 % nudge about the smallest one: the seed
        // of the flip (a perfect spin would stay on the axis to rounding for a while).
        const float w = 6.0f;
        Vector3 wl(0.0f);
        wl[mid_] = w;
        wl[lo] = 0.01f * w;
        b.angVel = b.rotation() * wl;
        sigma_ = w * std::sqrt((I[mid_] - I[lo]) * (I[hi] - I[mid_]) / (I[lo] * I[hi]));
        I_ = I;
        L0_ = angularMomentum(b);
        E0_ = sim.rigid.kineticEnergy();
        flips_ = 0;
        prevSpin_ = w;
        worstL_ = worstE_ = 0;
    }
    void afterStep(Simulation& sim) override {
        const RigidBody& b = sim.rigid.bodies()[handle_];
        // The spin about the body's own middle axis changes sign at every flip.
        Vector3 axis(0.0f);
        axis[mid_] = 1;
        const float spin = dot(b.rotation() * axis, b.angVel);
        if (prevSpin_ * spin < 0) ++flips_;
        prevSpin_ = spin;
        worstL_ = std::max(worstL_, std::fabs(angularMomentum(b) / L0_ - 1));
        worstE_ = std::max(worstE_, std::fabs(sim.rigid.kineticEnergy() / E0_ - 1));
        Probe::set("lesson/flips", flips_);
        Probe::set("lesson/L drift", worstL_);
        Probe::set("lesson/energy drift", worstE_);
    }
    void describe(const Simulation& sim, RenderSnapshot& s) const override {
        s.info.push_back({"Моменты инерции I1 < I2 < I3", format("%.4f / %.4f / %.4f кг·м²", minComp(I_), I_[mid_], maxComp(I_))});
        s.info.push_back({"Инкремент неустойчивости σ", format("%.2f 1/с", sigma_)});
        s.info.push_back({"Период переворотов, теория ≈ 20.4/σ", format("%.2f с", 20.4f / sigma_)});
        s.info.push_back({"Переворотов за время сцены", format("%d (ожидается ~%.0f)", flips_, sim.time() * sigma_ / 20.4f * 2)});
        s.info.push_back({"Дрейф |L| / энергии (Нётер)", format("%.1e / %.1e", worstL_, worstE_)});
        s.plots.push_back({"Спин о среднюю ось, рад/с", prevSpin_});
    }

private:
    int handle_ = -1, mid_ = 2, flips_ = 0;
    Vector3 I_;
    float sigma_ = 0, L0_ = 1, E0_ = 1, prevSpin_ = 0, worstL_ = 0, worstE_ = 0;
};

// --- Newton's cradle ----------------------------------------------------------------------------
// Five steel balls on ropes, touching in a row. A ball lifted and released hits the row; momentum
// and energy pass through it, and only the last ball flies out (both m v and m v^2 / 2 must be
// conserved at once, so one ball out at the same speed is the only answer). The balls are a hair
// apart (1 mm) so the hits happen one after another, as in a real cradle, where the tiny elastic
// delay does the same. Restitution 1, no friction, no damping: an ideal cradle.
class NewtonCradleScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.rigid.params.sleeping = false;
        sim.rigid.params.linearDamping = sim.rigid.params.angularDamping = 0;
        sim.rigid.params.rollingResistance = 0;
        // A hit at 0.8 m/s must bounce: below the default threshold (1 m/s) a touch is treated as
        // resting, so the threshold goes down for this lesson.
        sim.rigid.params.restitutionThreshold = 0.2f;
        // Gaps under the slop (4 mm) count as touching, and a row of "touching" balls behaves like
        // one lump. No stacks here, so the slop can be small and the 5 mm gaps stay real gaps.
        sim.rigid.params.slop = 0.0005f;
        sim.vis.showSlice = sim.vis.showStreamlines = false;
    }
    void build(Simulation& sim) override {
        sim.useRigidArena(AABB({-1, 0, -1}, {1, 1.5f, 1}));
        const float r = 0.03f, rope = 0.5f, gap = 0.005f, x0 = -2 * (2 * r + gap);
        const Vector3 pivot0(x0, 1.0f, 0.0f);
        balls_.clear();
        for (int i = 0; i < 5; ++i) {
            const Vector3 pivot = pivot0 + Vector3(i * (2 * r + gap), 0, 0);
            Vector3 p = pivot - Vector3(0, rope, 0);
            if (i == 0) p = pivot + Vector3(-std::sin(angle_), -std::cos(angle_), 0) * rope; // lifted aside
            const int id = sim.rigid.addSphere(p, r, 7800.0f, {0.75f, 0.75f, 0.8f}); // steel
            RigidBody& b = sim.rigid.bodies()[id];
            b.restitution = 1.0f;
            b.friction = b.staticFriction = 0;
            DistanceJoint& j = sim.rigid.addDistanceJoint(id, -1, p, pivot);
            j.rope = true;
            balls_.push_back(id);
        }
        const RigidBody& b0 = sim.rigid.bodies()[balls_[0]];
        restY_ = pivot0.y - rope;
        E0_ = b0.mass * 9.81f * (b0.pos.y - restY_); // the lifted ball's potential energy
    }
    void describe(const Simulation& sim, RenderSnapshot& s) const override {
        float px = 0, kinetic = 0, potential = 0;
        int flying = 0;
        for (int id : balls_) {
            const RigidBody& b = sim.rigid.bodies()[id];
            px += b.mass * b.vel.x;
            kinetic += 0.5f * b.mass * length2(b.vel);
            potential += b.mass * 9.81f * (b.pos.y - restY_);
            if (b.pos.x > 0 && length(b.vel) > 0.1f) ++flying;
        }
        s.info.push_back({"Импульс по x", format("%.4f кг·м/с", px)});
        s.info.push_back({"Энергия: кинетическая + потенциальная", format("%.4f + %.4f = %.4f Дж (в начале %.4f)", kinetic, potential, kinetic + potential, E0_)});
        s.info.push_back({"Шаров летит с правого края", format("%d (теория: 1)", flying)});
        s.plots.push_back({"Энергия, Дж", kinetic + potential});
        Probe::set("lesson/energy J", kinetic + potential);
        Probe::set("lesson/flying balls", flying);
    }

private:
    std::vector<int> balls_;
    float angle_ = 20.0f * 3.14159265f / 180.0f, restY_ = 0.5f, E0_ = 0;
};

// --- Galileo's tower -----------------------------------------------------------------------------
// Two spheres of the same size, one of wood (500 kg/m^3) and one of steel (7800), and a thin plate
// as light as a feather (20 kg/m^3), dropped from 4 m. Without air all three land together at
// t = sqrt(2h/g). With the air switched on, quadratic drag F = 1/2 rho Cd A v^2 slows the light
// ones: the plate settles at its terminal speed v_t = sqrt(2 m g / (rho Cd A)) - the feather of
// the Apollo 15 demonstration, done on the Moon to make the air go away.
class GalileoTowerScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.vis.showSlice = sim.vis.showStreamlines = false;
        sim.rigid.params.sleeping = false;
    }
    void build(Simulation& sim) override {
        sim.useRigidArena(AABB({-2, 0, -2}, {2, 5, 2}));
        const float h = 4.0f, r = 0.05f, t = 0.0025f;
        bodies_[0] = sim.rigid.addSphere({-0.4f, h + r, 0}, r, 500.0f, {0.8f, 0.6f, 0.3f});   // wood
        bodies_[1] = sim.rigid.addSphere({0.0f, h + r, 0}, r, 7800.0f, {0.6f, 0.62f, 0.66f}); // steel
        bodies_[2] = sim.rigid.addBox({0.4f, h + t, 0}, {0.05f, t, 0.05f}, Quaternion(), 20.0f, {0.95f, 0.95f, 0.9f}); // the "feather"
        for (int id : bodies_) sim.rigid.bodies()[id].restitution = 0; // land dead: the first touch is the landing
        bottom_[0] = bottom_[1] = r;
        bottom_[2] = t;
        area_[0] = area_[1] = 3.14159265f * r * r;
        area_[2] = 0.1f * 0.1f; // the face of the plate, falling flat
        drag_[0] = drag_[1] = 0.47f;
        drag_[2] = 1.2f;
        for (float& f : landed_) f = -1;
        for (bool& f : fell_) f = false;
        predicted_ = std::sqrt(2 * h / 9.81f);
    }
    void afterStep(Simulation& sim) override {
        for (int k = 0; k < 3; ++k) {
            const RigidBody& b = sim.rigid.bodies()[bodies_[k]];
            // The landing: the first frame the fall stops (the floor took the speed away).
            if (landed_[k] < 0 && fell_[k] && b.vel.y > -1.0f) landed_[k] = sim.time();
            if (b.vel.y < -1.0f) fell_[k] = true;
            if (!air_) continue;
            // Quadratic drag of the air (rho 1.2 kg/m^3) as an impulse over the frame.
            const float v = length(b.vel);
            const Vector3 F = b.vel * (-0.5f * 1.2f * drag_[k] * area_[k] * v);
            sim.rigid.applyExternalWrench(bodies_[k], F * sim.frameDt, Vector3(0.0f));
        }
        Probe::set("lesson/plate speed", length(sim.rigid.bodies()[bodies_[2]].vel));
    }
    void describe(const Simulation& sim, RenderSnapshot& s) const override {
        auto landed = [&](int k) { return landed_[k] < 0 ? std::string("ещё летит") : format("%.3f с", landed_[k]); };
        s.info.push_back({"Время падения, теория sqrt(2h/g)", format("%.3f с", predicted_)});
        s.info.push_back({"Дерево / сталь / пластина упали за", landed(0) + " / " + landed(1) + " / " + landed(2)});
        const RigidBody& plate = sim.rigid.bodies()[bodies_[2]];
        const float vt = std::sqrt(2 * plate.mass * 9.81f / (1.2f * drag_[2] * area_[2]));
        s.info.push_back({"Воздух", air_ ? format("есть: предельная скорость пластины %.2f м/с, сейчас %.2f", vt, length(plate.vel)) : "нет (вакуум): все падают одинаково"});
        s.plots.push_back({"Скорость пластины, м/с", length(plate.vel)});
    }
    std::vector<SceneParam> params() const override {
        SceneParam air{"Воздух", air_ ? 1.0f : 0.0f, 0, 1, 1, 0, true, "Квадратичное сопротивление F = ½ ρ Cd A v² на каждое тело; без него все три падают за одно время"};
        return {air};
    }
    void setParam(int index, float value) override {
        if (index == 0) air_ = value > 0.5f;
    }

private:
    int bodies_[3] = {-1, -1, -1};
    float bottom_[3] = {0, 0, 0}, area_[3] = {0, 0, 0}, drag_[3] = {0, 0, 0}, landed_[3] = {-1, -1, -1};
    bool fell_[3] = {false, false, false};
    float predicted_ = 0;
    bool air_ = false;
};

} // namespace

void addLessonSamples(std::vector<SampleEntry>& out) {
    out.push_back({Preset::Dzhanibekov, "Уроки", "Урок: эффект Джанибекова — гайка кувыркается в невесомости", [] { return std::unique_ptr<Scene>(new DzhanibekovScene); }});
    out.push_back({Preset::NewtonCradle, "Уроки", "Урок: колыбель Ньютона — импульс и энергия", [] { return std::unique_ptr<Scene>(new NewtonCradleScene); }});
    out.push_back({Preset::GalileoTower, "Уроки", "Урок: башня Галилея — массы падают одинаково", [] { return std::unique_ptr<Scene>(new GalileoTowerScene); }});
}

} // namespace rf
