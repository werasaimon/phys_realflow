// The meta-objects of a GraphScene: what the roles of one entity made in the solvers, changed one
// entity at a time between frames. Turning water into jelly removes the water's particle group and
// builds a soft body from the entity's SOURCE (its shape, EntityShapes.cpp) where the water is now,
// moving as the water moved; the rest of the scene goes on as if nothing happened. Nothing is ever
// rebuilt from a meta-object's own state, so any role can be turned into any other and back.
// The world (gravity, box, gas, the gas's heat source, combustion, a plasma's field) is set when
// the scene loads; a role that needs it switched on asks for a reload (rebuildEntity returns false).
#include "scene/SceneGraph.h"

#include "scene/Simulation.h"

#include <algorithm>
#include <cmath>

namespace rf {

namespace {

// Euler angles in degrees of a rotation R = Ry(y) Rx(x) Rz(z) - the inverse of entityRotation().
// R = [.. .. sy cx; cx sz cx cz -sx; .. .. cy cx] gives x from R12, y from R02 / R22 and z from
// R10 / R11; at x = +-90 degrees (cx = 0) y and z turn about the same axis: z is set to 0.
Vector3 eulerDegrees(const Quaternion& q) {
    const Matrix3x3 R = q.toMatrix3x3();
    const float sx = clampv(-R.m[1][2], -1.0f, 1.0f);
    const float x = std::asin(sx);
    float y, z;
    if (std::fabs(sx) < 0.99999f) {
        y = std::atan2(R.m[0][2], R.m[2][2]);
        z = std::atan2(R.m[1][0], R.m[1][1]);
    } else {
        y = std::atan2(-R.m[2][0], R.m[0][0]);
        z = 0.0f;
    }
    return {radToDeg(x), radToDeg(y), radToDeg(z)};
}

bool samePose(const Entity& a, const Entity& b) {
    return length(a.position - b.position) < 1e-6f && length(a.rotationDeg - b.rotationDeg) < 1e-4f;
}

} // namespace

// ---------------------------------------------------------------------------
// Ids and lookups
// ---------------------------------------------------------------------------
GraphScene::GraphScene(SceneGraph graph) : graph_(std::move(graph)) {
    uint32_t next = 1;
    for (const Entity& e : graph_.entities) next = std::max(next, e.id + 1);
    for (Entity& e : graph_.entities)
        if (e.id == 0) e.id = next++;
}

int GraphScene::indexOf(uint32_t id) const {
    for (size_t i = 0; i < graph_.entities.size(); ++i)
        if (graph_.entities[i].id == id) return int(i);
    return -1;
}

uint32_t GraphScene::entityIdOfBody(int body) const {
    return body >= 0 && body < int(entityIdOfBody_.size()) ? entityIdOfBody_[size_t(body)] : 0;
}

int GraphScene::entityOfBody(int body) const {
    const uint32_t id = entityIdOfBody(body);
    return id != 0 ? indexOf(id) : -1;
}

const std::vector<MetaObject>& GraphScene::metaObjects(uint32_t id) const {
    static const std::vector<MetaObject> none;
    auto it = meta_.find(id);
    return it != meta_.end() ? it->second : none;
}

// ---------------------------------------------------------------------------
// One entity between frames
// ---------------------------------------------------------------------------
bool GraphScene::rebuildEntity(Simulation& sim, const Entity& updated) {
    const int index = indexOf(updated.id);
    if (index < 0) return addEntity(sim, updated);
    // Where the entity is now: its body or particles have moved since the graph pose was written.
    // If the editor moved it (a different pose than the stored one), the editor's pose wins.
    const bool moved = !samePose(updated, graph_.entities[size_t(index)]);
    const LiveState live = liveState(sim, index);
    Entity next = updated;
    if (!moved) {
        next.position = live.position;
        next.rotationDeg = live.rotationDeg;
        standOnParticles(next, live);
    }
    destroyMetaObjects(sim, updated.id);
    const std::string prefix = graph_.entities[size_t(index)].name + ":";
    notes_.erase(std::remove_if(notes_.begin(), notes_.end(), [&](const std::string& n) { return n.rfind(prefix, 0) == 0; }),
                 notes_.end());
    graph_.entities[size_t(index)] = next;
    buildEntity(sim, index);
    if (!moved) carryVelocity(sim, next.id, live);
    return !needsReload(sim, next);
}

bool GraphScene::addEntity(Simulation& sim, const Entity& entity) {
    Entity e = entity;
    if (e.id == 0 || indexOf(e.id) >= 0) {
        uint32_t next = 1;
        for (const Entity& other : graph_.entities) next = std::max(next, other.id + 1);
        e.id = next;
    }
    graph_.entities.push_back(e);
    buildEntity(sim, int(graph_.entities.size()) - 1);
    return !needsReload(sim, e);
}

void GraphScene::removeEntity(Simulation& sim, uint32_t id) {
    const int index = indexOf(id);
    if (index < 0) return;
    destroyMetaObjects(sim, id);
    meta_.erase(id);
    graph_.entities.erase(graph_.entities.begin() + index);
}

// Everything the entity's roles made goes: its body (the slot becomes free for the next body), its
// particle group, its magnet and its emitter. The other entities' bodies keep their indices and
// their particle groups their numbers, so their meta-objects stay valid.
void GraphScene::destroyMetaObjects(Simulation& sim, uint32_t id) {
    auto it = meta_.find(id);
    if (it == meta_.end()) return;
    for (const MetaObject& m : it->second) {
        switch (m.kind) {
        case MetaObject::Kind::RigidBody:
            if (sim.rigid.isAlive(m.handle)) sim.rigid.destroyBody(m.handle);
            if (m.handle < int(entityIdOfBody_.size())) entityIdOfBody_[size_t(m.handle)] = 0;
            break;
        case MetaObject::Kind::SoftBody:
        case MetaObject::Kind::Liquid:
        case MetaObject::Kind::Cloth: sim.particles.removeGroup(m.handle); break;
        case MetaObject::Kind::Magnet:
        case MetaObject::Kind::Emitter: break; // their lists are cleaned below, by entity
        }
    }
    it->second.clear();
    for (size_t k = magnetEntity_.size(); k-- > 0;)
        if (magnetEntity_[k] == id) {
            magnetEntity_.erase(magnetEntity_.begin() + long(k));
            magnetBody_.erase(magnetBody_.begin() + long(k));
            magnetMoment_.erase(magnetMoment_.begin() + long(k));
        }
    emitters_.erase(std::remove_if(emitters_.begin(), emitters_.end(), [id](const EmitterRef& r) { return r.entity == id; }),
                    emitters_.end());
}

// The live pose and motion of an entity, from what it is made of: a rigid body gives its pose and
// velocities (the entity's orientation is the body's turned back by bodyToEntity); particles give
// their centre and mean velocity (the orientation stays as it was). Nothing made: the graph pose.
GraphScene::LiveState GraphScene::liveState(const Simulation& sim, int index) const {
    const Entity& e = graph_.entities[size_t(index)];
    LiveState s;
    s.position = e.position;
    s.rotationDeg = e.rotationDeg;
    for (const MetaObject& m : metaObjects(e.id)) {
        if (m.kind == MetaObject::Kind::RigidBody && sim.rigid.isAlive(m.handle)) {
            const RigidBody& b = sim.rigid.bodies()[size_t(m.handle)];
            s.position = b.pos;
            s.rotationDeg = eulerDegrees(b.rot * m.bodyToEntity);
            s.velocity = b.vel;
            s.angularVelocity = b.angVel;
            return s;
        }
        if (m.kind == MetaObject::Kind::SoftBody || m.kind == MetaObject::Kind::Liquid || m.kind == MetaObject::Kind::Cloth) {
            Vector3 sumX(0.0f), sumV(0.0f);
            int n = 0;
            float lowest = kInf;
            for (size_t i = 0; i < sim.particles.size(); ++i)
                if (sim.particles.groupOf(int(i)) == m.handle) {
                    sumX += sim.particles.positions()[i];
                    sumV += sim.particles.velocities()[i];
                    lowest = std::min(lowest, sim.particles.positions()[i].y);
                    ++n;
                }
            if (n > 0) {
                s.position = sumX / float(n);
                s.velocity = sumV / float(n);
                s.fromParticles = true;
                s.lowest = lowest;
            }
            return s;
        }
    }
    return s;
}

// A shape formed from particles stands where they lie: a puddle's centre is a few millimetres above
// the floor, and the whole shape put there would be half inside it (and thrown out of it). So the
// shape is lifted until its bottom is at the lowest particle centre - the first layer of new liquid
// is put exactly there, one radius above the floor the puddle lay on - keeping its x and z.
void GraphScene::standOnParticles(Entity& next, const LiveState& live) const {
    if (!live.fromParticles) return;
    const float shapeBottom = entityMesh(next, graph_.baseDirectory).bounds().lo.y;
    if (shapeBottom < live.lowest) next.position.y += live.lowest - shapeBottom;
}

// The new meta-objects move as the old ones did: a dynamic body takes the velocities, every new
// particle the mean velocity (it starts at rest).
void GraphScene::carryVelocity(Simulation& sim, uint32_t id, const LiveState& live) const {
    for (const MetaObject& m : metaObjects(id)) {
        if (m.kind == MetaObject::Kind::RigidBody) {
            RigidBody& b = sim.rigid.bodies()[size_t(m.handle)];
            if (b.invMass == 0) continue;
            b.vel = live.velocity;
            b.angVel = live.angularVelocity;
        } else if (m.kind == MetaObject::Kind::SoftBody || m.kind == MetaObject::Kind::Liquid || m.kind == MetaObject::Kind::Cloth) {
            for (size_t i = 0; i < sim.particles.size(); ++i)
                if (sim.particles.groupOf(int(i)) == m.handle) sim.particles.addVelocity(int(i), live.velocity);
        }
    }
}

// A role that works only with something set when the scene loaded: smoke, heat and fire need the
// gas; fire needs the combustion model; a heat source is the gas's one source, set at load.
bool GraphScene::needsReload(const Simulation& sim, const Entity& e) const {
    if (!e.visible) return false;
    const bool smokeOrHeat = e.heat.enabled || e.flammable.enabled || (e.emitter.enabled && (e.emitter.smoke > 0 || e.emitter.temperature > 0));
    if (smokeOrHeat && !graph_.world.gas) return true;
    if (e.flammable.enabled && !sim.grid.combustion.enabled) return true;
    if (e.heat.enabled && !sim.grid.source.enabled) return true;
    return false;
}

} // namespace rf
