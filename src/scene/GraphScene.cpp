// A SceneGraph running as a Scene: every entity goes to the solver its roles name - a rigid body
// to RigidWorld, a soft body and liquid to the particle system, a heat source to the gas - and the
// magnets push and turn each other every frame (Magnets.h). The editor's scenes and the samples go
// through the same Simulation facade, so a scene clicked together in the editor behaves exactly as
// the same scene written in code.
#include "scene/SceneGraph.h"

#include "core/Format.h"
#include "core/Mesh.h"
#include "core/Probe.h"
#include "scene/Magnets.h"
#include "scene/Simulation.h"

namespace rf {

namespace {

// The entity's orientation: yaw about y, then pitch about x, then roll about z (degrees).
Quaternion entityRotation(const Entity& e) {
    const Vector3& d = e.rotationDeg;
    return Quaternion::fromAxisAngle({0, 1, 0}, degToRad(d.y)) * Quaternion::fromAxisAngle({1, 0, 0}, degToRad(d.x)) *
           Quaternion::fromAxisAngle({0, 0, 1}, degToRad(d.z));
}

// The shape as a closed mesh around the origin in the entity's own frame: size is the full size,
// round shapes take the diameter from size.x and the height (along y) from size.y.
TriMesh entityMesh(const Entity& e) {
    const Vector3 s = e.size;
    TriMesh m;
    switch (e.shape) {
    case ShapeKind::Box: m = primitives::box(s * 0.5f); break;
    case ShapeKind::Plane: m = primitives::box(Vector3(0.5f * s.x, 0.01f, 0.5f * s.z)); break;
    case ShapeKind::Sphere: m = primitives::sphere(0.5f * s.x, 24, 12); break;
    case ShapeKind::Cylinder: // the primitive's axis is z: turn it to y
        m = primitives::cylinder(0.5f * s.x, s.y, 24);
        m.transform(Quaternion::fromAxisAngle({1, 0, 0}, -0.5f * kPi).toMatrix3x3(), Vector3(1.0f), Vector3(0.0f));
        break;
    case ShapeKind::Cone: // the primitive's apex points to -x: turn it up (+y)
        m = primitives::cone(0.5f * s.x, s.y, 24);
        m.transform(Quaternion::fromAxisAngle({0, 0, 1}, -0.5f * kPi).toMatrix3x3(), Vector3(1.0f), Vector3(0.0f));
        break;
    }
    m.translate(-m.bounds().center()); // centred on the entity's position
    return m;
}

bool hasLiquid(const SceneGraph& g) {
    for (const Entity& e : g.entities)
        if (e.visible && e.liquid.enabled) return true;
    return false;
}

// What the shape is made of: exactly one of rigid, soft, liquid, cloth - the first in that order
// when several are switched on (describe() names the ones left out).
enum class Matter { None, Rigid, Soft, Liquid, Cloth };

Matter madeOf(const Entity& e) {
    if (e.rigid.enabled) return Matter::Rigid;
    if (e.soft.enabled) return Matter::Soft;
    if (e.liquid.enabled) return Matter::Liquid;
    if (e.cloth.enabled) return Matter::Cloth;
    return Matter::None;
}

// The "made of" roles switched on besides the one used, as words for the readings.
std::string leftOutMatter(const Entity& e, Matter used) {
    std::string out;
    auto add = [&](bool on, Matter m, const char* word) {
        if (on && m != used) out += out.empty() ? word : std::string(", ") + word;
    };
    add(e.rigid.enabled, Matter::Rigid, "твёрдое");
    add(e.soft.enabled, Matter::Soft, "мягкое");
    add(e.liquid.enabled, Matter::Liquid, "жидкость");
    add(e.cloth.enabled, Matter::Cloth, "ткань");
    return out;
}

bool anyVisible(const SceneGraph& g, bool (*has)(const Entity&)) {
    for (const Entity& e : g.entities)
        if (e.visible && has(e)) return true;
    return false;
}

} // namespace

// ---------------------------------------------------------------------------
// Scene
// ---------------------------------------------------------------------------
void GraphScene::configure(Simulation& sim) {
    const WorldSettings& w = graph_.world;
    sim.setGravity(w.gravity);
    if (!w.gas) return;
    // A closed room of air: walls on every side, no wind; heat makes it rise.
    sim.grid.params.domainSize = w.size;
    sim.grid.params.resolutionX = 48;
    sim.grid.params.inflowSpeed = 1.0f; // only the reference scale of the pressure coefficient
    for (BoundaryType& b : sim.grid.params.bc) b = BoundaryType::Wall;
    sim.grid.params.smokeRake = false;
    sim.grid.params.heatBuoyancy = 3.0f;
    sim.grid.magnetic.enabled = w.magneticGas;
    sim.vis.showSlice = false;
    sim.vis.showStreamlines = false;
    // Fire: a flammable cloth takes heat from the gas, decomposes and burns (the fire sample); the
    // combustion model has to run for that. Cloth in the gas uses the coarser 2-radius spacing of
    // the samples: the grid cannot resolve finer folds anyway.
    sim.grid.combustion.enabled = anyVisible(graph_, [](const Entity& e) { return e.flammable.enabled; });
    if (anyVisible(graph_, [](const Entity& e) { return e.cloth.enabled; })) sim.particles.params.clothSpacing = 2.0f;
    for (const Entity& e : graph_.entities)
        if (e.visible && e.heat.enabled) { // the first hot entity is the gas's heat source
            sim.grid.source.enabled = true;
            sim.grid.source.center = e.position;
            sim.grid.source.radius = 0.5f * e.size.x;
            sim.grid.source.temperature = e.heat.temperature;
            sim.grid.source.smoke = e.heat.smoke;
            break;
        }
}

void GraphScene::build(Simulation& sim) {
    const WorldSettings& w = graph_.world;
    const AABB box({-0.5f * w.size.x, 0.0f, -0.5f * w.size.z}, {0.5f * w.size.x, w.size.y, 0.5f * w.size.z});
    if (w.gas) sim.useGasBox({0.5f, 0.0f, 0.5f}); // floor at y = 0
    else if (hasLiquid(graph_)) sim.useLiquidTank(box);
    else sim.useRigidArena(box);
    magnetBody_.clear();
    magnetMoment_.clear();
    entityOfBody_.clear();
    emitters_.clear();
    notes_.clear();
    maxMagnetForce_ = 0;
    sim.particles.emitter.enabled = false; // a liquid emitter switches it on in releaseFromEmitter()
    for (int i = 0; i < int(graph_.entities.size()); ++i) {
        if (!graph_.entities[i].visible) continue; // hidden: no body, no particles, no emission
        const int bodiesBefore = int(sim.rigid.bodies().size());
        addEntity(sim, i);
        // Whatever bodies this entity made belong to it: the editor picks an entity by its body.
        entityOfBody_.resize(sim.rigid.bodies().size(), -1);
        for (int b = bodiesBefore; b < int(sim.rigid.bodies().size()); ++b) entityOfBody_[b] = i;
    }
    if (sim.grid.magnetic.enabled) addMagnetFieldToGas(sim);
}

// One entity into the solvers: first what it is made of, then what it also does.
void GraphScene::addEntity(Simulation& sim, int index) {
    const Entity& e = graph_.entities[size_t(index)];
    const Quaternion q = entityRotation(e);
    const Matter matter = madeOf(e);
    const std::string leftOut = leftOutMatter(e, matter);
    if (!leftOut.empty()) notes_.push_back(e.name + ": сделано из одного — не использовано: " + leftOut);
    int body = -1;
    switch (matter) {
    case Matter::Rigid: body = addRigidBody(sim, e, e.rigid.fixed || e.shape == ShapeKind::Plane); break;
    case Matter::Soft: {
        TriMesh m = entityMesh(e);
        m.transform(q.toMatrix3x3(), Vector3(1.0f), e.position);
        sim.particles.addSoftBody(m, e.soft.density, e.soft.stiffness, e.color);
        break;
    }
    case Matter::Liquid: sim.particles.addBlock(AABB(e.position - e.size * 0.5f, e.position + e.size * 0.5f)); break;
    case Matter::Cloth: addCloth(sim, e); break;
    case Matter::None:
        if (e.magnet.enabled) body = addRigidBody(sim, e, true); // a magnet with nothing else stays put
        break;
    }
    if (e.magnet.enabled && body < 0) notes_.push_back(e.name + ": магнит бывает только твёрдым телом");
    if (e.flammable.enabled && matter != Matter::Cloth) notes_.push_back(e.name + ": горит только ткань (пока)");
    if (e.emitter.enabled) addEmitter(sim, index, body);
    if (e.magnet.enabled && body >= 0) {
        // The moment is given in the entity's frame; the body's own frame may differ (a hull body
        // lives in its principal frame), so turn it into the body frame once here.
        const Vector3 world = q.toMatrix3x3() * e.magnet.moment;
        magnetBody_.push_back(body);
        magnetMoment_.push_back(sim.rigid.bodies()[body].rotation().transposed() * world);
    }
}

int GraphScene::addRigidBody(Simulation& sim, const Entity& e, bool fixed) {
    const Quaternion q = entityRotation(e);
    const float density = fixed ? 0.0f : e.rigid.density; // density 0: a static body
    int id;
    switch (e.shape) {
    case ShapeKind::Box: id = sim.rigid.addBox(e.position, e.size * 0.5f, q, density, e.color); break;
    case ShapeKind::Plane: id = sim.rigid.addBox(e.position, Vector3(0.5f * e.size.x, 0.01f, 0.5f * e.size.z), q, 0.0f, e.color); break;
    case ShapeKind::Sphere: id = sim.rigid.addSphere(e.position, 0.5f * e.size.x, density, e.color); break;
    default: id = sim.rigid.addConvex(entityMesh(e), e.position, q, density, e.color); break;
    }
    RigidBody& b = sim.rigid.bodies()[id];
    if (e.rigid.enabled) {
        b.friction = e.rigid.friction;
        b.staticFriction = std::max(b.staticFriction, e.rigid.friction);
        b.restitution = e.rigid.restitution;
    }
    if (!fixed) {
        b.vel = e.rigid.velocity;
        b.angVel = e.rigid.angularVelocity;
    }
    return id;
}

// A sheet of cloth in place of the shape. A Plane becomes the sheet itself (size.x by size.z in its
// plane); any other shape gives its top face. In the entity's frame the sheet starts at the corner
// (-x, -z) and runs along +x (the cloth's rows) and +z (its columns), so the role's edges map onto
// ParticleSystem::addCloth's pin bits: -z edge = first row (16), +z edge = last row (32),
// -x edge = first column (64), +x edge = last column (128). "The top row" (role bit 16) is the edge
// that is highest in the world after the turn - the rod of a curtain; for a level sheet, the -z edge.
void GraphScene::addCloth(Simulation& sim, const Entity& e) const {
    const Matrix3x3 R = entityRotation(e).toMatrix3x3();
    const float lift = e.shape == ShapeKind::Plane ? 0.0f : 0.5f * e.size.y;
    const Vector3 corner = e.position + R * Vector3(-0.5f * e.size.x, lift, -0.5f * e.size.z);
    const Vector3 u = R * Vector3(e.size.x, 0, 0), v = R * Vector3(0, 0, e.size.z);
    const int edgeBit[4] = {16, 32, 64, 128}; // -z, +z, -x, +x
    const Vector3 edgeMiddle[4] = {corner + u * 0.5f, corner + u * 0.5f + v, corner + v * 0.5f, corner + u + v * 0.5f};
    const int roleEdges = e.cloth.pinnedEdges;
    int pinMask = 0;
    if (roleEdges & 4) pinMask |= edgeBit[0];
    if (roleEdges & 8) pinMask |= edgeBit[1];
    if (roleEdges & 1) pinMask |= edgeBit[2];
    if (roleEdges & 2) pinMask |= edgeBit[3];
    if (roleEdges & 16) {
        int top = 0;
        for (int k = 1; k < 4; ++k)
            if (edgeMiddle[k].y > edgeMiddle[top].y + 1e-4f) top = k;
        pinMask |= edgeBit[top];
    }
    ClothMaterial m;
    m.areaDensity = e.cloth.areaDensity;
    m.bendCompliance = e.cloth.bendCompliance;
    if (!e.cloth.tearable) m.strengthWarp = m.strengthWeft = 0.0f; // 0: the threads never break
    m.flammable = e.flammable.enabled;
    sim.particles.addCloth(corner, u, v, m, pinMask, e.color);
}

void GraphScene::addEmitter(Simulation& sim, int index, int body) {
    EmitterRef ref;
    ref.entity = index;
    ref.body = body;
    if (body >= 0) // entity frame = body frame * (body frame at the start)^T * entity frame at the start
        ref.bodyToEntity = sim.rigid.bodies()[body].rotation().transposed() * entityRotation(graph_.entities[size_t(index)]).toMatrix3x3();
    emitters_.push_back(ref);
}

// What an emitter releases in one frame, at its entity's pose of now. The release point is just
// outside the shape - the gas treats a body as solid and would drop anything put inside it - on
// the side the release velocity points to, else behind a moving body (a trail), else on top.
// Units: smoke 1 = 10 litres of dense smoke per second; hot gas is released at the same rate (at
// least 10 litres a second) with the emitter's temperature above ambient. Liquid: the particle
// system has one nozzle, so the first liquid emitter drives it (others are named in the readings).
void GraphScene::releaseFromEmitter(Simulation& sim, const EmitterRef& ref, bool& liquidDone) const {
    const Entity& e = graph_.entities[size_t(ref.entity)];
    const EmitterRole& em = e.emitter;
    Vector3 centre = e.position, moving(0.0f);
    Matrix3x3 R = entityRotation(e).toMatrix3x3();
    float reach = 0.5f * maxComp(e.size);
    if (ref.body >= 0) {
        const RigidBody& b = sim.rigid.bodies()[ref.body];
        centre = b.pos;
        moving = b.vel;
        R = b.rotation() * ref.bodyToEntity;
        reach = b.boundingRadius();
    }
    const Vector3 release = R * em.velocity;
    Vector3 dir = length(release) > 1e-6f ? normalize(release) : length(moving) > 0.1f ? -normalize(moving) : Vector3(0, 1, 0);
    const bool gas = graph_.world.gas;
    const Vector3 point = centre + dir * (reach + (gas ? sim.grid.dx() : sim.particles.params.particleRadius * 2));
    if (gas && (em.smoke > 0 || em.temperature > 0)) {
        const float dt = sim.frameDt;
        const float smokeVolume = 0.01f * em.smoke * dt;                     // [m^3]
        const float hotVolume = std::max(smokeVolume, 0.01f * dt);           // [m^3]
        const float rho = sim.grid.params.fluidDensity;
        const float heat = rho * sim.grid.combustion.specificHeat * hotVolume * em.temperature; // [J]
        sim.grid.addEmission({point, 0.0f, heat, smokeVolume});
        if (length(release) > 0) sim.grid.addImpulse(point, release * (rho * hotVolume)); // the jet's momentum
    }
    if (em.liquid > 0 && !liquidDone) {
        ParticleEmitter& nozzle = sim.particles.emitter;
        nozzle.enabled = true;
        nozzle.position = point;
        nozzle.direction = dir;
        nozzle.radius = std::max(0.25f * std::min(e.size.x, e.size.z), 2.0f * sim.particles.params.particleRadius);
        nozzle.speed = std::max(length(release), 1.0f);
        liquidDone = true;
    }
}

// A plasma feels the magnets: their dipole fields are the gas's background field, from the vector
// potential A = mu0/4pi (m x r) / |r|^3 so that div B = 0 on the grid (as the magnetosphere sample).
// Set once at build from the magnets' starting poses; magnets that move afterwards do not drag the
// field along yet (future work: rebuild the background when a magnet has moved more than a cell).
void GraphScene::addMagnetFieldToGas(Simulation& sim) const {
    std::vector<Vector3> centre, moment;
    std::vector<float> core;
    for (size_t i = 0; i < magnetBody_.size(); ++i) {
        const RigidBody& b = sim.rigid.bodies()[magnetBody_[i]];
        centre.push_back(b.pos);
        moment.push_back(b.rotation() * magnetMoment_[i]);
        core.push_back(std::max(b.boundingRadius(), 1e-3f)); // inside the magnet: no singularity
    }
    sim.grid.magnetic.setBackgroundFromPotential([centre, moment, core](const Vector3& x) {
        Vector3 A(0.0f);
        for (size_t i = 0; i < centre.size(); ++i) {
            const Vector3 r = x - centre[i];
            const float d = std::max(length(r), core[i]);
            A += cross(moment[i], r) * (1e-7f / (d * d * d));
        }
        return A;
    });
}

void GraphScene::afterStep(Simulation& sim) {
    // The force of this frame acts as an impulse at its end (the next frame's contacts see it).
    maxMagnetForce_ = magnetBody_.size() >= 2 ? applyMagnetForces(sim.rigid, magnetBody_, magnetMoment_, sim.frameDt) : 0.0f;
    Probe::set("magnets/max force N", maxMagnetForce_);
    // Emitters release at the pose their body has now: the trail follows it.
    bool liquidDone = false;
    for (const EmitterRef& ref : emitters_) releaseFromEmitter(sim, ref, liquidDone);
    if (!emitters_.empty() && graph_.world.gas) Probe::set("gas/smoke dm3", sim.grid.totalSmoke() * 1000.0f);
}

void GraphScene::describe(const Simulation& sim, RenderSnapshot& s) const {
    s.info.push_back({"Сущностей", format("%zu", graph_.entities.size())});
    s.info.push_back({"Магнитов", format("%zu", magnetBody_.size())});
    if (!magnetBody_.empty()) s.info.push_back({"Наибольшая магнитная сила", format("%.4g Н", maxMagnetForce_)});
    if (!emitters_.empty()) s.info.push_back({"Излучателей", format("%zu", emitters_.size())});
    // Roles that need the gas, in a world without it: say so instead of doing nothing silently.
    bool needsGas = false;
    for (const Entity& e : graph_.entities)
        needsGas |= e.visible && (e.heat.enabled || e.flammable.enabled || (e.emitter.enabled && (e.emitter.smoke > 0 || e.emitter.temperature > 0)));
    if (needsGas && !graph_.world.gas) s.info.push_back({"Дым, тепло, огонь", "нужен газ: включите «газ» в мире"});
    int liquidEmitters = 0;
    for (const EmitterRef& ref : emitters_) liquidEmitters += graph_.entities[size_t(ref.entity)].emitter.liquid > 0 ? 1 : 0;
    if (liquidEmitters > 1) s.info.push_back({"Струи жидкости", format("работает первая из %d (сопло одно)", liquidEmitters)});
    for (const std::string& note : notes_) s.info.push_back({"Роли", note});
    (void)sim;
}

} // namespace rf
