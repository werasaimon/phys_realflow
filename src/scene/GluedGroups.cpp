// Glued groups of a GraphScene: the rigid members of a group marked `glued` become ONE rigid body.
// Each member's collider (its world mesh, taken into the group's frame) is a part of a compound
// shape; the body's mass is the sum of the dynamic members' masses spread over the compound's
// volume; its look is every member's geometry merged (drawn in the first member's colour). The body
// is fixed only when every glued member is. Every member keeps a meta-object pointing at the shared
// body, so its live pose (liveState) reads back through it; changing one member rebuilds the whole
// group from the graph at the group's live pose. Members that are not rigid (soft, cloth, liquid) are
// built on their own and do not follow the body.
#include "scene/SceneGraph.h"

#include "core/Mesh.h"
#include "scene/Simulation.h"

#include <algorithm>

namespace rf {

// A member of a glued group that goes into the group's body: rigid, or a collider alone (static).
bool GraphScene::gluedMember(int index) const {
    if (flatGlue_[size_t(index)] == 0) return false;
    const Entity& e = flat_[size_t(index)];
    const bool otherMatter = e.soft.enabled || e.liquid.enabled || e.cloth.enabled;
    return e.rigid.enabled || (!otherMatter && e.collider.enabled);
}

void GraphScene::buildGluedGroup(Simulation& sim, uint32_t group) {
    std::vector<int> members;
    for (int i = 0; i < int(flat_.size()); ++i)
        if (flatGlue_[size_t(i)] == group && gluedMember(i) && flat_[size_t(i)].visible && shapeUsable(flat_[size_t(i)]))
            members.push_back(i);
    const SceneObject* g = findObject(graph_, group);
    if (members.empty() || !g) return;
    Vector3 gp;
    Quaternion gq;
    worldPose(graph_, *g, gp, gq);
    const Quaternion toGroup = gq.conjugate();
    std::vector<TriMesh> hulls, looks;
    float mass = 0;
    bool fixed = true;
    for (int i : members) { // colliders and looks in the group's frame
        const Entity& e = flat_[size_t(i)];
        TriMesh hull = colliderMesh(e, graph_.baseDirectory), look = entityMesh(e, graph_.baseDirectory);
        for (Vector3& v : hull.positions) v = toGroup.rotate(v - gp);
        for (Vector3& v : look.positions) v = toGroup.rotate(v - gp);
        hulls.push_back(std::move(hull));
        looks.push_back(std::move(look));
        const bool memberFixed = !e.rigid.enabled || e.rigid.fixed || e.shape == ShapeKind::Plane;
        const EntityCollider c = entityCollider(e, graph_.baseDirectory);
        if (!memberFixed && c.shape) mass += e.rigid.density * c.shape->volume();
        fixed = fixed && memberFixed;
    }
    auto shape = std::make_shared<CompoundShape>(hulls, primitives::merge(looks));
    const float density = fixed ? 0.0f : mass / std::max(shape->volume(), 1e-9f);
    const Entity& first = flat_[size_t(members[0])];
    const Quaternion rot = (gq * Quaternion::fromMatrix3x3(shape->principalRotation())).normalized();
    const int body = sim.rigid.addBody(shape, gp + gq.rotate(shape->centerOfMass()), rot, density, first.color);
    RigidBody& b = sim.rigid.bodies()[size_t(body)];
    b.friction = first.rigid.friction;
    b.staticFriction = std::max(b.staticFriction, first.rigid.friction);
    b.restitution = first.rigid.restitution;
    if (!fixed) b.vel = first.rigid.velocity, b.angVel = first.rigid.angularVelocity;
    glueBody_[group] = body;
    if (int(entityIdOfBody_.size()) <= body) entityIdOfBody_.resize(size_t(body) + 1, 0);
    entityIdOfBody_[size_t(body)] = first.id;
    for (int i : members) {
        const Entity& e = flat_[size_t(i)];
        std::vector<MetaObject>& meta = meta_[e.id];
        meta.clear();
        const Quaternion er = entityRotation(e);
        meta.push_back({MetaObject::Kind::RigidBody, body, b.rot.conjugate() * er, er.conjugate().rotate(b.pos - e.position)});
        if (e.magnet.enabled) addMagnet(sim, i, body);
        if (e.emitter.enabled) addEmitter(sim, i, body);
    }
}

// The group again from the graph: each member's live pose (read through the shared body) written
// back, the body removed, the members' roles taken fresh from the graph, the body built again and
// given the old body's motion.
void GraphScene::rebuildGluedGroup(Simulation& sim, uint32_t group) {
    Vector3 vel(0.0f), angVel(0.0f);
    auto it = glueBody_.find(group);
    const bool had = it != glueBody_.end() && sim.rigid.isAlive(it->second);
    if (had) {
        vel = sim.rigid.bodies()[size_t(it->second)].vel;
        angVel = sim.rigid.bodies()[size_t(it->second)].angVel;
    }
    std::vector<int> leftTheBody; // members edited into something that is not rigid any more
    for (int i = 0; i < int(flat_.size()); ++i) {
        if (flatGlue_[size_t(i)] != group) continue;
        bool inBody = false;
        for (const MetaObject& m : metaObjects(flat_[size_t(i)].id))
            inBody = inBody || (m.kind == MetaObject::Kind::RigidBody && had && m.handle == it->second);
        if (!inBody && !gluedMember(i)) continue; // a soft or cloth member stays as it is
        if (had) {
            const LiveState live = liveState(sim, i);
            flat_[size_t(i)].position = live.position;
            flat_[size_t(i)].rotationDeg = live.rotationDeg;
            writeBackPose(i);
        }
        destroyMetaObjects(sim, flat_[size_t(i)].id);
        const int ai = authoredIndexOf(flat_[size_t(i)].id);
        if (ai < 0) continue;
        Entity fresh = resolveInstance(graph_, graph_.entities[size_t(ai)]); // roles as the graph has them now
        fresh.position = flat_[size_t(i)].position;
        fresh.rotationDeg = flat_[size_t(i)].rotationDeg;
        fresh.parent = 0;
        fresh.visible = effectivelyVisible(graph_, graph_.entities[size_t(ai)]);
        flat_[size_t(i)] = fresh;
        if (!gluedMember(i)) leftTheBody.push_back(i);
    }
    if (it != glueBody_.end()) glueBody_.erase(it);
    buildGluedGroup(sim, group);
    for (int i : leftTheBody) buildEntity(sim, i);
    auto now = glueBody_.find(group);
    if (had && now != glueBody_.end() && sim.rigid.bodies()[size_t(now->second)].invMass > 0) {
        sim.rigid.bodies()[size_t(now->second)].vel = vel;
        sim.rigid.bodies()[size_t(now->second)].angVel = angVel;
    }
}

} // namespace rf
