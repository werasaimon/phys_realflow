// The channels of one object and of a group of objects (see Channels.h).
//
// A rigid body: where it is, how fast it moves and turns, its energy, and the force its contacts
// put on it - the impulses of the last solver step (normal and friction) divided by that step:
// F = J / h. A resting body's contact force is its weight, upwards; the ground under a stack feels
// the weight of the whole stack, downwards - Newton's third law, readable on a chart.
// A particle group (a soft body, a cloth, a block of liquid): its centre of mass, mean velocity and
// energy; a cloth also the highest tension in its threads, per metre of fabric.
// A group: the totals of its members - mass, centre of mass, momentum, energy.
#include "scene/ChannelMath.h"
#include "scene/Simulation.h"

#include "core/Parallel.h"

#include <algorithm>

namespace rf {

namespace {

// What is summed over the members of anything: mass, first moment Σ m x, momentum Σ m v, energy.
struct Totals {
    double mass = 0, kinetic = 0, potential = 0;
    Double3 firstMoment, momentum;
    long particles = 0;
    Totals& operator+=(const Totals& o) {
        mass += o.mass, kinetic += o.kinetic, potential += o.potential, particles += o.particles;
        firstMoment += o.firstMoment, momentum += o.momentum;
        return *this;
    }
    Double3 centre() const { return mass > 0 ? firstMoment * (1.0 / mass) : Double3(); }
};

// One rigid body's totals (a static body has no mass to count: it is the ground).
Totals bodyTotals(const Simulation& sim, const RigidBody& b) {
    Totals t;
    if (b.invMass == 0) return t;
    const double m = b.mass;
    const Double3 x(b.pos), v(b.vel), g(sim.gravity()), floor(potentialFloor(sim));
    const Vector3 wBody = b.rotation().transposed() * b.angVel, invI = b.invInertiaLocal;
    const double spin = (invI.x > 0 ? wBody.x * wBody.x / invI.x : 0.0) + (invI.y > 0 ? wBody.y * wBody.y / invI.y : 0.0) +
                        (invI.z > 0 ? wBody.z * wBody.z / invI.z : 0.0); // ω·Iω in the principal frame
    t.mass = m;
    t.kinetic = 0.5 * m * dot(v, v) + 0.5 * spin;
    t.potential = -m * dot(g, x - floor);
    t.firstMoment = x * m;
    t.momentum = v * m;
    return t;
}

// The totals of one particle group (every particle of it: O(n) over the particles).
Totals particleGroupTotals(const Simulation& sim, int group) {
    const ParticleSystem& ps = sim.particles;
    const Double3 g(sim.gravity()), floor(potentialFloor(sim));
    return parallelSum<Totals>(int(ps.size()), [&](int begin, int end) {
        Totals t;
        for (int i = begin; i < end; ++i) {
            if (ps.groupOf(i) != group) continue;
            ++t.particles;
            const float w = ps.invMasses()[size_t(i)];
            if (w == 0) continue; // pinned: held in place, like the ground
            const double m = 1.0 / w;
            const Double3 x(ps.positions()[size_t(i)]), v(ps.velocities()[size_t(i)]);
            t.mass += m;
            t.kinetic += 0.5 * m * dot(v, v);
            t.potential -= m * dot(g, x - floor);
            t.firstMoment += x * m;
            t.momentum += v * m;
        }
        return t;
    }, 1024);
}

Totals objectTotals(const Simulation& sim, const ObjectRef& o) {
    if (o.kind == ObjectRef::Kind::Body)
        return sim.rigid.isAlive(o.index) ? bodyTotals(sim, sim.rigid.bodies()[size_t(o.index)]) : Totals();
    return particleGroupTotals(sim, o.index);
}

// The measurements every object and group share: mass, centre (and height above the floor of the
// scene box), velocity of the centre of mass, energy.
void addTotals(std::vector<Measurement>& out, const std::string& label, ChannelGroup group, const Totals& t, const Simulation& sim) {
    auto add = [&](const char* quantity, const char* unit, const char* name, double value) {
        out.push_back(measured(label + "/" + quantity, unit, name, group, value));
    };
    const Double3 c = t.centre(), v = t.mass > 0 ? t.momentum * (1.0 / t.mass) : Double3();
    add("mass", "кг", "Масса", t.mass);
    add("centre x", "м", "Центр масс, x", c.x);
    add("centre y", "м", "Центр масс, y", c.y);
    add("centre z", "м", "Центр масс, z", c.z);
    add("height", "м", "Высота центра масс над дном сцены", c.y - potentialFloor(sim).y);
    add("speed", "м/с", "Скорость центра масс", v.length());
    add("velocity x", "м/с", "Скорость, x", v.x);
    add("velocity y", "м/с", "Скорость, y", v.y);
    add("velocity z", "м/с", "Скорость, z", v.z);
    add("momentum", "Н·с", "Импульс |P|", t.momentum.length());
    add("kinetic energy", "Дж", "Энергия движения (с вращением)", t.kinetic);
    add("potential energy", "Дж", "Энергия высоты над дном сцены", t.potential);
    add("energy", "Дж", "Механическая энергия", t.kinetic + t.potential);
}

// The force of a body's contacts: the impulses of the last solver step over its length, F = J / h.
// A sleeping body is the exception: the solver skips sleeping islands, so there are no impulses to
// sum - but sleep means nothing moves the body, the forces on it balance (Newton's first law), and
// its contacts hold exactly its weight, F = -m g. (A static body touching only sleeping bodies has
// no solved contact left: it reads 0 until they wake - docs/10-verification.md, channels.)
Double3 contactForce(const RigidWorld& w, const RigidBody& b, int index, const Vector3& gravity, int& points) {
    const Vector3 J = w.contactImpulseOn(index, points);
    if (b.sleeping && b.invMass > 0) return Double3(gravity) * -double(b.mass);
    const double h = w.lastStepDt();
    return h > 0 ? Double3(J) * (1.0 / h) : Double3();
}

// A rigid body's own channels beside the totals: spin, sleep and the force of its contacts.
void addBodyExtras(std::vector<Measurement>& out, const Simulation& sim, const ObjectRef& o) {
    const RigidWorld& w = sim.rigid;
    const RigidBody& b = w.bodies()[size_t(o.index)];
    int points = 0;
    const Double3 F = contactForce(w, b, o.index, sim.gravity(), points);
    auto add = [&](const char* quantity, const char* unit, const char* name, double value) {
        out.push_back(measured(o.label + "/" + quantity, unit, name, ChannelGroup::Object, value));
    };
    add("angular speed", "рад/с", "Угловая скорость |ω|", length(b.angVel));
    add("asleep", "да/нет", "Спит (1) или движется (0)", b.sleeping ? 1.0 : 0.0);
    add("contact force", "Н", "Сила контактов |F|", F.length());
    add("contact force x", "Н", "Сила контактов, x", F.x);
    add("contact force y", "Н", "Сила контактов, y", F.y);
    add("contact force z", "Н", "Сила контактов, z", F.z);
    add("contacts", "шт", "Точек контакта", points);
}

// A cloth's highest thread tension per metre of fabric (tension / the width the thread stands for).
void addClothTension(std::vector<Measurement>& out, const Simulation& sim, const ObjectRef& o) {
    for (const Cloth& c : sim.particles.cloths()) {
        if (c.group != o.index) continue;
        double worst = 0;
        for (const DistanceConstraint& d : c.constraints) {
            if (d.broken || (d.kind != DistanceConstraint::Warp && d.kind != DistanceConstraint::Weft)) continue;
            const double width = d.kind == DistanceConstraint::Warp ? c.warpWidth : c.weftWidth;
            if (width > 0) worst = std::max(worst, double(d.force) / width);
        }
        out.push_back(measured(o.label + "/max thread tension", "Н/м", "Наибольшее натяжение нитей", ChannelGroup::Object, worst));
    }
}

} // namespace

std::vector<Measurement> measureObject(const Simulation& sim, const ObjectRef& o) {
    std::vector<Measurement> out;
    const bool body = o.kind == ObjectRef::Kind::Body;
    if (body && !sim.rigid.isAlive(o.index)) return out;
    const Totals t = objectTotals(sim, o);
    if (!body && t.particles == 0) return out; // the group is gone
    addTotals(out, o.label, ChannelGroup::Object, t, sim);
    if (body) addBodyExtras(out, sim, o);
    else {
        out.push_back(measured(o.label + "/particles", "шт", "Частиц", ChannelGroup::Object, double(t.particles)));
        addClothTension(out, sim, o);
    }
    return out;
}

std::vector<Measurement> measureGroup(const Simulation& sim, const GroupRef& group) {
    Totals t;
    for (const ObjectRef& o : group.members) t += objectTotals(sim, o);
    std::vector<Measurement> out;
    addTotals(out, group.label, ChannelGroup::Group, t, sim);
    out.push_back(measured(group.label + "/members", "шт", "Объектов в группе", ChannelGroup::Group, double(group.members.size())));
    return out;
}

} // namespace rf
