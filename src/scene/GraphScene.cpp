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

#include <algorithm>

namespace rf {

namespace {

bool hasLiquid(const std::vector<Entity>& entities) {
    for (const Entity& e : entities)
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

bool anyVisible(const std::vector<Entity>& entities, bool (*has)(const Entity&)) {
    for (const Entity& e : entities)
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
    sim.grid.combustion.enabled = anyVisible(flat_, [](const Entity& e) { return e.flammable.enabled; });
    if (anyVisible(flat_, [](const Entity& e) { return e.cloth.enabled; })) sim.particles.params.clothSpacing = 2.0f;
    for (const Entity& e : flat_)
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
    else if (hasLiquid(flat_)) sim.useLiquidTank(box);
    else sim.useRigidArena(box);
    magnetBody_.clear();
    magnetMoment_.clear();
    magnetEntity_.clear();
    entityIdOfBody_.clear();
    emitters_.clear();
    meta_.clear();
    notes_.clear();
    maxMagnetForce_ = 0;
    glueBody_.clear();
    sim.particles.emitter.enabled = false; // a liquid emitter switches it on in releaseFromEmitter()
    // Every simulated entity (flat_: world poses, instances resolved, array copies); the rigid
    // members of a glued group are built afterwards, together, as one body.
    std::vector<uint32_t> glued;
    for (int i = 0; i < int(flat_.size()); ++i) {
        if (!gluedMember(i)) buildEntity(sim, i);
        else if (std::find(glued.begin(), glued.end(), flatGlue_[size_t(i)]) == glued.end()) glued.push_back(flatGlue_[size_t(i)]);
    }
    for (uint32_t g : glued) buildGluedGroup(sim, g);
    if (sim.grid.magnetic.enabled) addMagnetFieldToGas(sim);
}

// Whether the entity's shape can be built: a model file must be readable, and a model cannot be
// made cloth yet. A shape that cannot be built is left out with a note in the readings.
bool GraphScene::shapeUsable(const Entity& e) {
    if (e.shape != ShapeKind::Mesh) return true;
    std::string error;
    if (entityLocalMesh(e, graph_.baseDirectory, &error).empty()) {
        notes_.push_back(e.name + ": модель не прочитана — " + error);
        return false;
    }
    if (madeOf(e) == Matter::Cloth) {
        notes_.push_back(e.name + ": ткань из модели пока не поддерживается - только из плоскости или верхней грани формы");
        return false;
    }
    return true;
}

// One entity into the solvers: first what it is made of, then what it also does; everything made
// is recorded as the entity's meta-objects, so it can be removed and built again alone. An entity
// with no role at all is geometry only: no body, no collision - the editor just draws it; a hidden
// one takes no part at all. A collider without the rigid role is a static obstacle (a wall, a
// floor); the rigid role without a collider falls back to the Auto collider (noted). Each role has its own add* function that builds its part fresh from the
// entity's source (the shape or model, EntityShapes.cpp).
void GraphScene::buildEntity(Simulation& sim, int index) {
    const Entity& e = flat_[size_t(index)];
    std::vector<MetaObject>& meta = meta_[e.id];
    meta.clear();
    if (!e.visible || entityIsGeometryOnly(e) || !shapeUsable(e)) return;
    const Matter matter = madeOf(e);
    const std::string leftOut = leftOutMatter(e, matter);
    if (!leftOut.empty()) notes_.push_back(e.name + ": сделано из одного — не использовано: " + leftOut);
    int body = -1;
    switch (matter) {
    case Matter::Rigid: body = addRigid(sim, index); break;
    case Matter::Soft: meta.push_back({MetaObject::Kind::SoftBody, addSoft(sim, index), Quaternion()}); break;
    case Matter::Liquid: meta.push_back({MetaObject::Kind::Liquid, addLiquid(sim, index), Quaternion()}); break;
    case Matter::Cloth: meta.push_back({MetaObject::Kind::Cloth, addCloth(sim, index), Quaternion()}); break;
    case Matter::None: // a collider alone is a static obstacle; a magnet with nothing else stays put too
        if (e.collider.enabled || e.magnet.enabled) body = addRigid(sim, index);
        break;
    }
    if (matter == Matter::Rigid && !e.collider.enabled) notes_.push_back(e.name + ": у твёрдого тела нет коллайдера — взят «Авто»");
    if (e.collider.enabled && matter != Matter::Rigid && matter != Matter::None)
        notes_.push_back(e.name + ": коллайдер нужен только твёрдому телу — не использован");
    if (body >= 0) { // the body and which entity it belongs to (the editor picks entities by body)
        const RigidBody& b = sim.rigid.bodies()[size_t(body)];
        const Quaternion bodyToEntity = b.rot.conjugate() * entityRotation(e);
        const Vector3 bodyOffset = entityRotation(e).conjugate().rotate(b.pos - e.position);
        meta.push_back({MetaObject::Kind::RigidBody, body, bodyToEntity, bodyOffset});
        if (int(entityIdOfBody_.size()) <= body) entityIdOfBody_.resize(size_t(body) + 1, 0);
        entityIdOfBody_[size_t(body)] = e.id;
    }
    if (e.flammable.enabled && matter != Matter::Cloth) notes_.push_back(e.name + ": горит только ткань (пока)");
    if (e.magnet.enabled) addMagnet(sim, index, body);
    if (e.emitter.enabled) addEmitter(sim, index, body);
}

// Rigid: one body whose collider is chosen apart from its look (entityCollider: by default the
// geometry itself - a box or a sphere exactly, a model's convex decomposition, the hull of anything
// else; or a box / sphere / capsule / hull / decomposition at its own offset). The body is drawn
// with the entity's geometry (RigidBody::visualMesh, in the body's frame), so a horse on a capsule
// still looks like a horse. A fixed body, a plane, and a magnet with no "made of" role do not move
// (density 0). Returns the body index.
int GraphScene::addRigid(Simulation& sim, int index) {
    const Entity& e = flat_[size_t(index)];
    const bool fixed = !e.rigid.enabled || e.rigid.fixed || e.shape == ShapeKind::Plane;
    const Quaternion q = entityRotation(e);
    const float density = fixed ? 0.0f : e.rigid.density; // density 0: a static body
    const EntityCollider c = entityCollider(e, graph_.baseDirectory);
    if (!c.shape) return -1;
    const int id = sim.rigid.addBody(c.shape, e.position + q.rotate(c.position), (q * c.rotation).normalized(), density, e.color);
    RigidBody& b = sim.rigid.bodies()[id];
    // The geometry from the entity's frame into the body's: v_body = R_c^T (v - p_c).
    TriMesh look = entityLocalMesh(e, graph_.baseDirectory);
    const Matrix3x3 toBody = c.rotation.toMatrix3x3().transposed();
    for (Vector3& v : look.positions) v = toBody * (v - c.position);
    b.visualMesh = std::make_shared<const TriMesh>(std::move(look));
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

// Soft: the shape's surface mesh in the world, filled with particles held by shape matching.
// Returns its particle group.
int GraphScene::addSoft(Simulation& sim, int index) {
    const Entity& e = flat_[size_t(index)];
    const int body = sim.particles.addSoftBody(entityMesh(e, graph_.baseDirectory), e.soft.density, e.soft.stiffness, e.color);
    return sim.particles.softBodyGroup(body);
}

// Liquid: the shape's box (turned and moved as the entity) filled with water particles. Returns
// its particle group.
int GraphScene::addLiquid(Simulation& sim, int index) {
    const Entity& e = flat_[size_t(index)];
    return sim.particles.addBlock(entityMesh(e, graph_.baseDirectory).bounds());
}

// Magnet: a dipole riding on the entity's body (-1: none - a magnet has to be a rigid body). The
// moment is given in the entity's frame; the body's own frame may differ (a hull body lives in its
// principal frame), so it is turned into the body frame once here.
void GraphScene::addMagnet(Simulation& sim, int index, int body) {
    const Entity& e = flat_[size_t(index)];
    if (body < 0) {
        notes_.push_back(e.name + ": магнит бывает только твёрдым телом");
        return;
    }
    const Vector3 world = entityRotation(e).toMatrix3x3() * e.magnet.moment;
    magnetBody_.push_back(body);
    magnetMoment_.push_back(sim.rigid.bodies()[size_t(body)].rotation().transposed() * world);
    magnetEntity_.push_back(e.id);
    meta_[e.id].push_back({MetaObject::Kind::Magnet, body, Quaternion()});
}

// A sheet of cloth in place of the shape. A Plane becomes the sheet itself (size.x by size.z in its
// plane); any other shape gives its top face. In the entity's frame the sheet starts at the corner
// (-x, -z) and runs along +x (the cloth's rows) and +z (its columns), so the role's edges map onto
// ParticleSystem::addCloth's pin bits: -z edge = first row (16), +z edge = last row (32),
// -x edge = first column (64), +x edge = last column (128). "The top row" (role bit 16) is the edge
// that is highest in the world after the turn - the rod of a curtain; for a level sheet, the -z edge.
// Returns its particle group.
int GraphScene::addCloth(Simulation& sim, int index) {
    const Entity& e = flat_[size_t(index)];
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
    return sim.particles.clothGroup(sim.particles.addCloth(corner, u, v, m, pinMask, e.color));
}

// Emitter: remembered with the body it rides on; releaseFromEmitter() runs it every frame.
void GraphScene::addEmitter(Simulation& sim, int index, int body) {
    const Entity& e = flat_[size_t(index)];
    EmitterRef ref;
    ref.entity = e.id;
    ref.body = body;
    if (body >= 0) { // entity frame = body frame * (body frame at the start)^T * entity frame at the start
        const RigidBody& b = sim.rigid.bodies()[size_t(body)];
        ref.bodyToEntity = b.rotation().transposed() * entityRotation(e).toMatrix3x3();
        ref.bodyOffset = entityRotation(e).conjugate().rotate(b.pos - e.position);
    }
    emitters_.push_back(ref);
    meta_[e.id].push_back({MetaObject::Kind::Emitter, body, Quaternion()});
}

// What an emitter releases in one frame, at its entity's pose of now. The release point is just
// outside the shape - the gas treats a body as solid and would drop anything put inside it - on
// the side the release velocity points to, else behind a moving body (a trail), else on top.
// Units: smoke 1 = 10 litres of dense smoke per second; hot gas is released at the same rate (at
// least 10 litres a second) with the emitter's temperature above ambient. Liquid: the particle
// system has one nozzle, so the first liquid emitter drives it (others are named in the readings).
void GraphScene::releaseFromEmitter(Simulation& sim, const EmitterRef& ref, bool& liquidDone) const {
    const int index = indexOf(ref.entity);
    if (index < 0) return;
    const Entity& e = flat_[size_t(index)];
    const EmitterRole& em = e.emitter;
    Vector3 centre = e.position, moving(0.0f);
    Matrix3x3 R = entityRotation(e).toMatrix3x3();
    float reach = 0.5f * maxComp(e.size);
    if (ref.body >= 0) {
        const RigidBody& b = sim.rigid.bodies()[ref.body];
        moving = b.vel;
        R = b.rotation() * ref.bodyToEntity;
        centre = b.pos - R * ref.bodyOffset; // the entity's centre, not the collider's
        reach = b.boundingRadius() + length(ref.bodyOffset);
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

// The scene's visible lights, posed in the world, for the viewer's shading (the physics has no use
// for them). None: the snapshot's list stays empty and the viewer keeps its default light.
void describeLights(const SceneGraph& g, RenderSnapshot& s) {
    s.lights.clear();
    for (const Light& l : g.lights) {
        if (!effectivelyVisible(g, l)) continue;
        RenderSnapshot::LightInfo info;
        Quaternion rotation;
        worldPose(g, l, info.position, rotation);
        info.kind = int(l.kind);
        info.direction = lightDirection(g, l);
        info.color = l.color;
        info.intensity = l.intensity;
        info.range = l.range;
        info.coneDeg = l.coneDeg;
        info.softnessDeg = l.softnessDeg;
        info.shadows = l.shadows;
        s.lights.push_back(info);
    }
}

void GraphScene::setLightsAndCameras(const std::vector<Light>& lights, const std::vector<Camera>& cameras) {
    graph_.lights = lights;
    graph_.cameras = cameras;
}

void GraphScene::describe(const Simulation& sim, RenderSnapshot& s) const {
    describeLights(graph_, s);
    s.info.push_back({"Сущностей", format("%zu", flat_.size())});
    s.info.push_back({"Магнитов", format("%zu", magnetBody_.size())});
    if (!magnetBody_.empty()) s.info.push_back({"Наибольшая магнитная сила", format("%.4g Н", maxMagnetForce_)});
    if (!emitters_.empty()) s.info.push_back({"Излучателей", format("%zu", emitters_.size())});
    // Roles that need the gas, in a world without it: say so instead of doing nothing silently.
    bool needsGas = false;
    for (const Entity& e : flat_)
        needsGas |= e.visible && (e.heat.enabled || e.flammable.enabled || (e.emitter.enabled && (e.emitter.smoke > 0 || e.emitter.temperature > 0)));
    if (needsGas && !graph_.world.gas) s.info.push_back({"Дым, тепло, огонь", "нужен газ: включите «газ» в мире"});
    int liquidEmitters = 0;
    for (const EmitterRef& ref : emitters_) {
        const int index = indexOf(ref.entity);
        liquidEmitters += index >= 0 && flat_[size_t(index)].emitter.liquid > 0 ? 1 : 0;
    }
    if (liquidEmitters > 1) s.info.push_back({"Струи жидкости", format("работает первая из %d (сопло одно)", liquidEmitters)});
    for (const std::string& note : notes_) s.info.push_back({"Роли", note});
    (void)sim;
}

} // namespace rf
