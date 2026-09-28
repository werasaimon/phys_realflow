// The scene channels (see Channels.h): the name cards of every physical channel, and the sums that
// measure the whole scene - the energy, momentum and angular momentum of the rigid bodies, the
// energy of the particles and the density error of the liquid, the energy of the gas flow, its
// divergence and speed, the magnetic energy. Objects, groups: ObjectChannels.cpp; probes in the
// fields: FieldProbes.cpp.
//
// Every sum runs over the bodies / particles / cells once (O(n)), in double, with parallelSum -
// blocks that depend only on n, so the numbers are the same bits on any number of threads.
#include "scene/ChannelMath.h"
#include "scene/Simulation.h"

#include "core/Parallel.h"

#include <algorithm>

namespace rf {

// ---------------------------------------------------------------------------------------------
// The name cards
// ---------------------------------------------------------------------------------------------
const std::vector<ChannelInfo>& sceneChannelCatalog() {
    static const std::vector<ChannelInfo> cards = {
        {"rigid/kinetic energy", "Дж", "Твёрдые тела: энергия движения, Σ m v²/2", ChannelGroup::Scene},
        {"rigid/rotational energy", "Дж", "Твёрдые тела: энергия вращения, Σ ω·Iω/2", ChannelGroup::Scene},
        {"rigid/potential energy", "Дж", "Твёрдые тела: энергия высоты над дном сцены", ChannelGroup::Scene},
        {"rigid/mechanical energy", "Дж", "Твёрдые тела: полная механическая энергия", ChannelGroup::Scene},
        {"rigid/momentum", "Н·с", "Твёрдые тела: импульс |P|", ChannelGroup::Scene},
        {"rigid/momentum x", "Н·с", "Твёрдые тела: импульс, x", ChannelGroup::Scene},
        {"rigid/momentum y", "Н·с", "Твёрдые тела: импульс, y", ChannelGroup::Scene},
        {"rigid/momentum z", "Н·с", "Твёрдые тела: импульс, z", ChannelGroup::Scene},
        {"rigid/angular momentum", "Н·м·с", "Твёрдые тела: момент импульса |L| около центра масс", ChannelGroup::Scene},
        {"rigid/deepest penetration", "м", "Твёрдые тела: самое глубокое проникновение", ChannelGroup::Scene},
        {"particles/kinetic energy", "Дж", "Частицы: энергия движения", ChannelGroup::Scene},
        {"particles/potential energy", "Дж", "Частицы: энергия высоты над дном сцены", ChannelGroup::Scene},
        {"particles/elastic energy", "Дж", "Мягкие тела: упругая энергия, запасённая в деформации", ChannelGroup::Scene},
        {"liquid/density error mean", "%", "Жидкость: средняя ошибка плотности (сжатие)", ChannelGroup::Scene},
        {"liquid/density error max", "%", "Жидкость: наибольшая ошибка плотности (сжатие)", ChannelGroup::Scene},
        {"gas/kinetic energy", "Дж", "Газ: энергия течения", ChannelGroup::Scene},
        {"gas/max div u", "1/с", "Газ: наибольшая дивергенция скорости", ChannelGroup::Scene},
        {"gas/max speed", "м/с", "Газ: наибольшая скорость", ChannelGroup::Scene},
        {"mhd/energy", "Дж", "Магнитное поле: энергия, ∫B²/2μ₀ dV", ChannelGroup::Scene},
        {"scene/mechanical energy", "Дж", "Сцена: механическая энергия (тела, частицы, газ)", ChannelGroup::Scene},
    };
    return cards;
}

const ChannelInfo* describeChannel(const std::string& id) {
    for (const ChannelInfo& c : sceneChannelCatalog())
        if (c.id == id) return &c;
    return nullptr;
}

// A scene measurement with the card of its id from the catalog (one place for every name and unit).
static void put(std::vector<Measurement>& out, const char* id, double value) {
    Measurement m;
    m.info = *describeChannel(id);
    m.value = value;
    out.push_back(std::move(m));
}

Vector3 potentialFloor(const Simulation& sim) { return sim.particles.domain().lo; }

// ---------------------------------------------------------------------------------------------
// Rigid bodies
// ---------------------------------------------------------------------------------------------
namespace {

// The sums over the moving bodies. The angular momentum is first taken about the origin,
// L_O = Σ (x × m v + I ω), and moved to the centre of mass X of the system afterwards:
// L_X = L_O - X × P (the momentum P "carried" by the centre of mass is taken out).
struct RigidSums {
    double translational = 0, rotational = 0, potential = 0, mass = 0;
    Double3 momentum, firstMoment, angularAboutOrigin; // Σ m v, Σ m x, Σ (x × m v + I ω)
    RigidSums& operator+=(const RigidSums& o) {
        translational += o.translational, rotational += o.rotational, potential += o.potential, mass += o.mass;
        momentum += o.momentum, firstMoment += o.firstMoment, angularAboutOrigin += o.angularAboutOrigin;
        return *this;
    }
};

// One body's share: ½ m v², ½ ω·Iω (in its principal frame, where I is diagonal), -m g·(x - floor),
// m v, m x and x × m v + I ω.
void addBody(RigidSums& s, const RigidBody& b, const Double3& g, const Double3& floor) {
    const double m = b.mass;
    const Double3 x(b.pos), v(b.vel);
    const Vector3 wBody = b.rotation().transposed() * b.angVel;
    const Vector3 invI = b.invInertiaLocal;
    const Double3 Iw(invI.x > 0 ? wBody.x / invI.x : 0.0, invI.y > 0 ? wBody.y / invI.y : 0.0, invI.z > 0 ? wBody.z / invI.z : 0.0);
    const Vector3 IwWorld = b.rotation() * Vector3(float(Iw.x), float(Iw.y), float(Iw.z));
    s.translational += 0.5 * m * dot(v, v);
    s.rotational += 0.5 * dot(Double3(wBody), Iw);
    s.potential -= m * dot(g, x - floor);
    s.mass += m;
    s.momentum += v * m;
    s.firstMoment += x * m;
    s.angularAboutOrigin += cross(x, v * m) + Double3(IwWorld);
}

RigidSums sumRigid(const Simulation& sim) {
    const std::vector<RigidBody>& bodies = sim.rigid.bodies();
    const Double3 g(sim.gravity()), floor(potentialFloor(sim));
    return parallelSum<RigidSums>(int(bodies.size()), [&](int begin, int end) {
        RigidSums s;
        for (int i = begin; i < end; ++i) {
            const RigidBody& b = bodies[size_t(i)];
            if (b.alive && b.invMass > 0) addBody(s, b, g, floor); // a static body neither moves nor stores energy
        }
        return s;
    }, 64);
}

// The rigid channels; returns their mechanical energy.
double addRigidChannels(const Simulation& sim, std::vector<Measurement>& out) {
    const RigidSums s = sumRigid(sim);
    const Double3 centre = s.mass > 0 ? s.firstMoment * (1.0 / s.mass) : Double3();
    const Double3 L = s.angularAboutOrigin - cross(centre, s.momentum);
    const double mechanical = s.translational + s.rotational + s.potential;
    put(out, "rigid/kinetic energy", s.translational);
    put(out, "rigid/rotational energy", s.rotational);
    put(out, "rigid/potential energy", s.potential);
    put(out, "rigid/mechanical energy", mechanical);
    put(out, "rigid/momentum", s.momentum.length());
    put(out, "rigid/momentum x", s.momentum.x);
    put(out, "rigid/momentum y", s.momentum.y);
    put(out, "rigid/momentum z", s.momentum.z);
    put(out, "rigid/angular momentum", L.length());
    put(out, "rigid/deepest penetration", sim.rigid.deepestPenetration());
    return mechanical;
}

// ---------------------------------------------------------------------------------------------
// Particles and liquid
// ---------------------------------------------------------------------------------------------
struct ParticleSums {
    double kinetic = 0, potential = 0, error = 0, errorMax = 0;
    long liquid = 0;
    ParticleSums& operator+=(const ParticleSums& o) {
        kinetic += o.kinetic, potential += o.potential, error += o.error, liquid += o.liquid;
        errorMax = std::max(errorMax, o.errorMax);
        return *this;
    }
};

// ½ m v² and -m g·(x - floor) of every moving particle; the density error of the liquid is the
// compression the solver removes, max(ρ/ρ₀ - 1, 0) - the same measure as its own statistics.
ParticleSums sumParticles(const Simulation& sim) {
    const ParticleSystem& ps = sim.particles;
    const Double3 g(sim.gravity()), floor(potentialFloor(sim));
    const double invRho0 = 1.0 / ps.params.restDensity;
    return parallelSum<ParticleSums>(int(ps.size()), [&](int begin, int end) {
        ParticleSums s;
        for (int i = begin; i < end; ++i) {
            const float w = ps.invMasses()[size_t(i)];
            if (w > 0) {
                const double m = 1.0 / w;
                const Double3 v(ps.velocities()[size_t(i)]);
                s.kinetic += 0.5 * m * dot(v, v);
                s.potential -= m * dot(g, Double3(ps.positions()[size_t(i)]) - floor);
            }
            if (ps.phases()[size_t(i)] != uint8_t(ParticlePhase::Fluid) || ps.densities().size() <= size_t(i)) continue;
            const double e = std::max(ps.densities()[size_t(i)] * invRho0 - 1.0, 0.0);
            s.error += e;
            s.errorMax = std::max(s.errorMax, e);
            ++s.liquid;
        }
        return s;
    }, 1024);
}

// The particle channels; returns their mechanical energy.
double addParticleChannels(const Simulation& sim, std::vector<Measurement>& out) {
    const ParticleSums s = sumParticles(sim);
    put(out, "particles/kinetic energy", s.kinetic);
    put(out, "particles/potential energy", s.potential);
    const double elastic = sim.particles.softElasticEnergy(); // the soft bodies' squashed tetrahedra
    put(out, "particles/elastic energy", elastic);
    if (s.liquid > 0) {
        put(out, "liquid/density error mean", 100.0 * s.error / double(s.liquid));
        put(out, "liquid/density error max", 100.0 * s.errorMax);
    }
    return s.kinetic + s.potential + elastic;
}

// ---------------------------------------------------------------------------------------------
// Gas and magnetic field
// ---------------------------------------------------------------------------------------------
// The kinetic energy of the flow, ½ ρ Σ |u|² dx³ over the gas cells, with u the cell-centred
// velocity (the mean of the two faces per axis) - the measure of tests/ActionTests.cpp.
double gasKineticEnergy(const GasSolver& g) {
    const int nx = g.nx(), ny = g.ny(), nz = g.nz();
    const double cell = 0.5 * g.params.fluidDensity * double(g.dx()) * g.dx() * g.dx();
    return cell * parallelSum<double>(nz, [&](int k0, int k1) {
        double sum = 0;
        for (int k = k0; k < k1; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                    if (!g.solid(i, j, k)) sum += double(length2(g.cellVelocity(i, j, k)));
        return sum;
    }, 1);
}

// The gas channels (only while the scene has a gas grid); returns the flow's kinetic energy.
double addGasChannels(const Simulation& sim, std::vector<Measurement>& out) {
    const double kinetic = gasKineticEnergy(sim.grid);
    put(out, "gas/kinetic energy", kinetic);
    put(out, "gas/max div u", sim.grid.maxDivergence());
    put(out, "gas/max speed", sim.grid.maxVelocity());
    if (sim.grid.magnetic.enabled) put(out, "mhd/energy", sim.grid.magnetic.energy());
    return kinetic;
}

} // namespace

// ---------------------------------------------------------------------------------------------
// The scene, and everything asked for
// ---------------------------------------------------------------------------------------------
std::vector<Measurement> measureScene(const Simulation& sim) {
    std::vector<Measurement> out;
    double mechanical = 0;
    if (!sim.rigid.bodies().empty()) mechanical += addRigidChannels(sim, out);
    if (sim.particles.size() > 0) mechanical += addParticleChannels(sim, out);
    if (sim.mode() == SimMode::WindTunnel) mechanical += addGasChannels(sim, out);
    if (!out.empty()) put(out, "scene/mechanical energy", mechanical);
    return out;
}

std::vector<Measurement> measureAll(const Simulation& sim, const ChannelRequests& requests) {
    std::vector<Measurement> out;
    auto append = [&](std::vector<Measurement>&& part) {
        for (Measurement& m : part) out.push_back(std::move(m));
    };
    if (requests.scene) append(measureScene(sim));
    for (const ObjectRef& o : requests.objects) append(measureObject(sim, o));
    for (const GroupRef& g : requests.groups) append(measureGroup(sim, g));
    for (const FieldProbe& p : requests.probes) append(measureProbe(sim, p));
    return out;
}

} // namespace rf
