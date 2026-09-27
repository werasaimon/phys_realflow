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
        if (e.liquid.enabled) return true;
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
    for (const Entity& e : graph_.entities)
        if (e.heat.enabled) { // the first hot entity is the gas's heat source
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
    maxMagnetForce_ = 0;
    for (int i = 0; i < int(graph_.entities.size()); ++i) {
        const int bodiesBefore = int(sim.rigid.bodies().size());
        addEntity(sim, graph_.entities[i]);
        // Whatever bodies this entity made belong to it: the editor picks an entity by its body.
        entityOfBody_.resize(sim.rigid.bodies().size(), -1);
        for (int b = bodiesBefore; b < int(sim.rigid.bodies().size()); ++b) entityOfBody_[b] = i;
    }
    if (sim.grid.magnetic.enabled) addMagnetFieldToGas(sim);
}

void GraphScene::addEntity(Simulation& sim, const Entity& e) {
    const Quaternion q = entityRotation(e);
    int body = -1;
    if (e.rigid.enabled) body = addRigidBody(sim, e, e.rigid.fixed || e.shape == ShapeKind::Plane);
    else if (e.magnet.enabled) body = addRigidBody(sim, e, true); // a magnet without a body role stays put
    if (e.soft.enabled && !e.rigid.enabled) {
        TriMesh m = entityMesh(e);
        m.transform(q.toMatrix3x3(), Vector3(1.0f), e.position);
        sim.particles.addSoftBody(m, e.soft.density, e.soft.stiffness, e.color);
    }
    if (e.liquid.enabled) sim.particles.addBlock(AABB(e.position - e.size * 0.5f, e.position + e.size * 0.5f));
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
}

void GraphScene::describe(const Simulation& sim, RenderSnapshot& s) const {
    s.info.push_back({"Сущностей", format("%zu", graph_.entities.size())});
    s.info.push_back({"Магнитов", format("%zu", magnetBody_.size())});
    if (!magnetBody_.empty()) s.info.push_back({"Наибольшая магнитная сила", format("%.4g Н", maxMagnetForce_)});
    bool heat = false;
    for (const Entity& e : graph_.entities) heat |= e.heat.enabled;
    if (heat && !graph_.world.gas) s.info.push_back({"Источник тепла", "нужен газ: включите «газ» в мире"});
    (void)sim;
}

} // namespace rf
