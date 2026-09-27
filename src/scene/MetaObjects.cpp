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

bool samePose(const Entity& a, const Entity& b) {
    return length(a.position - b.position) < 1e-6f && length(a.rotationDeg - b.rotationDeg) < 1e-4f;
}

// Whether e is `target` or an instance whose chain of masters passes through it.
bool sharesRolesWith(const SceneGraph& g, const Entity& e, uint32_t target) {
    const Entity* cur = &e;
    for (size_t guard = 0; cur && guard <= g.entities.size(); ++guard) {
        if (cur->id == target) return true;
        if (cur->instanceOf == 0) return false;
        const Entity* next = nullptr;
        for (const Entity& o : g.entities)
            if (o.id == cur->instanceOf) next = &o;
        cur = next;
    }
    return false;
}

} // namespace

// ---------------------------------------------------------------------------
// Ids and lookups
// ---------------------------------------------------------------------------
// Every object gets a stable id (those without one are numbered after the largest); the flat list
// is what the simulation builds: world poses, instances resolved, array copies appended.
GraphScene::GraphScene(SceneGraph graph) : graph_(std::move(graph)) {
    uint32_t next = nextId(graph_);
    for (Entity& e : graph_.entities)
        if (e.id == 0) e.id = next++;
    for (Group& gr : graph_.groups)
        if (gr.id == 0) gr.id = next++;
    for (ArrayObject& a : graph_.arrays)
        if (a.id == 0) a.id = next++;
    flat_ = worldEntities(graph_);
    flatGlue_.assign(flat_.size(), 0);
    for (size_t i = 0; i < graph_.entities.size(); ++i) flatGlue_[i] = gluedGroupOf(graph_, graph_.entities[i]);
}

int GraphScene::indexOf(uint32_t id) const {
    for (size_t i = 0; i < flat_.size(); ++i)
        if (flat_[i].id == id) return int(i);
    return -1;
}

int GraphScene::authoredIndexOf(uint32_t id) const {
    for (size_t i = 0; i < graph_.entities.size(); ++i)
        if (graph_.entities[i].id == id) return int(i);
    return -1;
}

uint32_t GraphScene::entityIdOfBody(int body) const {
    return body >= 0 && body < int(entityIdOfBody_.size()) ? entityIdOfBody_[size_t(body)] : 0;
}

int GraphScene::entityOfBody(int body) const {
    const uint32_t id = entityIdOfBody(body);
    return id != 0 ? authoredIndexOf(id) : -1;
}

uint32_t GraphScene::objectIdOfBody(int body) const {
    for (const auto& [group, b] : glueBody_)
        if (b == body) return group;
    const uint32_t id = entityIdOfBody(body);
    if (!isArrayCopyId(id)) return id;
    for (const ArrayObject& a : graph_.arrays)
        if ((a.id & 0x7FFFu) == arrayOfCopyId(id)) return a.id;
    return 0;
}

// The graph keeps poses relative to the parents: a world pose read from the simulation is written
// back relative. Array copies have no entity of their own in the graph.
void GraphScene::writeBackPose(int index) {
    const Entity& w = flat_[size_t(index)];
    const int ai = authoredIndexOf(w.id);
    if (ai < 0) return;
    setWorldPose(graph_, graph_.entities[size_t(ai)], w.position, objectRotation(w));
}

const std::vector<MetaObject>& GraphScene::metaObjects(uint32_t id) const {
    static const std::vector<MetaObject> none;
    auto it = meta_.find(id);
    return it != meta_.end() ? it->second : none;
}

// ---------------------------------------------------------------------------
// One entity between frames
// ---------------------------------------------------------------------------
// The edited entity, every instance sharing its roles, and the arrays copying it are rebuilt; each
// keeps where it is now (its live pose) unless the editor moved it (then the editor's pose wins).
bool GraphScene::rebuildEntity(Simulation& sim, const Entity& updated) {
    const int ai = authoredIndexOf(updated.id);
    if (ai < 0) return addEntity(sim, updated);
    const bool moved = !samePose(updated, graph_.entities[size_t(ai)]);
    graph_.entities[size_t(ai)] = updated;
    std::vector<uint32_t> ids;
    for (const Entity& e : graph_.entities)
        if (sharesRolesWith(graph_, e, updated.id)) ids.push_back(e.id);
    bool ok = true;
    for (uint32_t id : ids) ok = rebuildOne(sim, id, !(id == updated.id && moved)) && ok;
    for (const ArrayObject& a : std::vector<ArrayObject>(graph_.arrays))
        if (std::find(ids.begin(), ids.end(), a.templateId) != ids.end()) ok = rebuildArray(sim, a) && ok;
    return ok;
}

// One simulated entity again, from the graph as it is now: its live pose written back first (unless
// the editor's pose should win), its meta-objects removed, built again, its motion carried on.
bool GraphScene::rebuildOne(Simulation& sim, uint32_t id, bool keepLivePose) {
    const int fi = indexOf(id), ai = authoredIndexOf(id);
    if (fi < 0 || ai < 0) return true; // not simulated, or an array's copy (rebuildArray builds those)
    if (gluedMember(fi)) { // the whole group is one body
        rebuildGluedGroup(sim, flatGlue_[size_t(fi)]);
        return !needsReload(sim, flat_[size_t(fi)]);
    }
    const LiveState live = liveState(sim, fi);
    if (keepLivePose) {
        flat_[size_t(fi)].position = live.position;
        flat_[size_t(fi)].rotationDeg = live.rotationDeg;
        writeBackPose(fi);
    }
    destroyMetaObjects(sim, id);
    const std::string prefix = flat_[size_t(fi)].name + ":";
    notes_.erase(std::remove_if(notes_.begin(), notes_.end(), [&](const std::string& n) { return n.rfind(prefix, 0) == 0; }),
                 notes_.end());
    Entity next = resolveInstance(graph_, graph_.entities[size_t(ai)]);
    next.position = flat_[size_t(fi)].position; // the world pose (kept live, or the editor's below)
    next.rotationDeg = flat_[size_t(fi)].rotationDeg;
    if (!keepLivePose) {
        Vector3 p;
        Quaternion q;
        worldPose(graph_, graph_.entities[size_t(ai)], p, q);
        next.position = p;
        next.rotationDeg = graph_.entities[size_t(ai)].parent != 0 ? eulerDegrees(q) : graph_.entities[size_t(ai)].rotationDeg;
    }
    next.parent = 0;
    next.visible = effectivelyVisible(graph_, graph_.entities[size_t(ai)]);
    if (keepLivePose) standOnParticles(next, live);
    flat_[size_t(fi)] = next;
    flatGlue_[size_t(fi)] = gluedGroupOf(graph_, graph_.entities[size_t(ai)]);
    if (keepLivePose) writeBackPose(fi);
    if (gluedMember(fi)) rebuildGluedGroup(sim, flatGlue_[size_t(fi)]); // it joined a glued group
    else buildEntity(sim, fi);
    if (keepLivePose) carryVelocity(sim, id, live);
    return !needsReload(sim, next);
}

bool GraphScene::addEntity(Simulation& sim, const Entity& entity) {
    Entity e = entity;
    if (e.id == 0 || findObject(graph_, e.id)) e.id = nextId(graph_);
    graph_.entities.push_back(e);
    Entity w = resolveInstance(graph_, e);
    Vector3 p;
    Quaternion q;
    worldPose(graph_, e, p, q);
    if (e.parent != 0) w.position = p, w.rotationDeg = eulerDegrees(q);
    w.parent = 0;
    w.visible = effectivelyVisible(graph_, e);
    flat_.push_back(w);
    flatGlue_.push_back(gluedGroupOf(graph_, e));
    const int fi = int(flat_.size()) - 1;
    if (gluedMember(fi)) rebuildGluedGroup(sim, flatGlue_[size_t(fi)]);
    else buildEntity(sim, fi);
    return !needsReload(sim, w);
}

void GraphScene::removeEntity(Simulation& sim, uint32_t id) {
    const int fi = indexOf(id), ai = authoredIndexOf(id);
    if (fi < 0 || ai < 0) return;
    const uint32_t group = gluedMember(fi) ? flatGlue_[size_t(fi)] : 0;
    if (group != 0) rebuildGluedGroup(sim, group); // live poses of the members into the graph first
    destroyMetaObjects(sim, id);
    meta_.erase(id);
    flat_.erase(flat_.begin() + fi);
    flatGlue_.erase(flatGlue_.begin() + fi);
    graph_.entities.erase(graph_.entities.begin() + ai);
    if (group != 0) rebuildGluedGroup(sim, group);
}

// The array's old copies go (each one's live state remembered by its stable id), the new copies are
// built; a copy that existed before stands where it is now and keeps moving as it moved.
bool GraphScene::rebuildArray(Simulation& sim, const ArrayObject& array) {
    ArrayObject a = array;
    auto it = std::find_if(graph_.arrays.begin(), graph_.arrays.end(), [&](const ArrayObject& o) { return o.id == a.id && a.id != 0; });
    if (it == graph_.arrays.end()) {
        if (a.id == 0 || findObject(graph_, a.id)) a.id = nextId(graph_);
        graph_.arrays.push_back(a);
    } else {
        *it = a;
    }
    std::unordered_map<uint32_t, LiveState> before;
    for (int i = int(flat_.size()) - 1; i >= 0; --i) {
        const uint32_t id = flat_[size_t(i)].id;
        if (!isArrayCopyId(id) || arrayOfCopyId(id) != (a.id & 0x7FFFu)) continue;
        before[id] = liveState(sim, i);
        destroyMetaObjects(sim, id);
        meta_.erase(id);
        flat_.erase(flat_.begin() + i);
        flatGlue_.erase(flatGlue_.begin() + i);
    }
    bool ok = true;
    for (Entity& c : expandArray(graph_, a)) {
        auto old = before.find(c.id);
        if (old != before.end()) {
            c.position = old->second.position;
            c.rotationDeg = old->second.rotationDeg;
            standOnParticles(c, old->second);
        }
        flat_.push_back(c);
        flatGlue_.push_back(0);
        buildEntity(sim, int(flat_.size()) - 1);
        if (old != before.end()) carryVelocity(sim, c.id, old->second);
        ok = !needsReload(sim, c) && ok;
    }
    return ok;
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
// their centre and mean velocity (the orientation stays as it was). Nothing made: the world pose.
// `index` is a flat index (every caller passes one). An array's copies live only in flat_, past the
// end of graph_.entities, so reading graph_.entities here read past the vector's end: a sporadic
// crash, caught by _GLIBCXX_ASSERTIONS in "the array grows to eighty during play".
GraphScene::LiveState GraphScene::liveState(const Simulation& sim, int index) const {
    const Entity& e = flat_[size_t(index)];
    LiveState s;
    s.position = e.position;
    s.rotationDeg = e.rotationDeg;
    for (const MetaObject& m : metaObjects(e.id)) {
        if (m.kind == MetaObject::Kind::RigidBody && sim.rigid.isAlive(m.handle)) {
            const RigidBody& b = sim.rigid.bodies()[size_t(m.handle)];
            const Quaternion entityRot = b.rot * m.bodyToEntity;
            s.position = b.pos - entityRot.rotate(m.bodyOffset); // the entity's centre, not the collider's
            s.rotationDeg = eulerDegrees(entityRot);
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
                // A tenth of a radius above the lowest centre: new particles are laid with a jitter
                // of +-2 % of a radius, and one laid lower than its old neighbour could come closer
                // to the floor than a radius - "inside a solid", so it was never made (316 of 343).
                s.lowest = lowest + 0.1f * sim.particles.params.particleRadius;
            }
            return s;
        }
    }
    return s;
}

// A shape formed from particles stands where they lie: a puddle's centre is a few millimetres above
// the floor, and the whole shape put there would be half inside it (and thrown out of it). So the
// shape is lifted until its bottom is at the lowest particle centre - the first layer of new liquid
// is put there (a tenth of a radius higher, see liveState), keeping its x and z.
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
