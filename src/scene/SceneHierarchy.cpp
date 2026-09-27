// Many objects at once, as data: the hierarchy (an object's pose relative to its parent group),
// instances (entities sharing the geometry and roles of a master, 3ds Max's Instance), and arrays
// (one object standing for a pattern of copies, Blender's Array modifier / Houdini's copy-to-points).
// Everything here only reads the SceneGraph and computes poses; GraphScene turns the result into
// meta-objects, and an editor draws it the same way while the scene is being edited.
#include "scene/SceneGraph.h"

#include <algorithm>
#include <cmath>

namespace rf {

namespace {

template <class T>
const SceneObject* findIn(const std::vector<T>& list, uint32_t id) {
    for (const T& o : list)
        if (o.id == id) return &o;
    return nullptr;
}

const Entity* findEntity(const SceneGraph& g, uint32_t id) {
    for (const Entity& e : g.entities)
        if (e.id == id) return &e;
    return nullptr;
}

// The parents of an object, nearest first; stops at a missing parent or a repeated one (a cycle).
std::vector<const SceneObject*> ancestors(const SceneGraph& g, const SceneObject& o) {
    std::vector<const SceneObject*> chain;
    std::vector<uint32_t> seen{o.id};
    for (uint32_t p = o.parent; p != 0;) {
        if (std::find(seen.begin(), seen.end(), p) != seen.end()) break;
        const SceneObject* parent = findObject(g, p);
        if (!parent) break;
        chain.push_back(parent);
        seen.push_back(p);
        p = parent->parent;
    }
    return chain;
}

// A small fixed-seed generator for the jitter of array copies: the same offsets every run.
float jitterNumber(uint32_t& state) {
    state = state * 1664525u + 1013904223u;
    return float(state >> 8) / 16777216.0f * 2.0f - 1.0f; // -1 .. 1
}

// Where copy n of an array sits in the array's own frame (before the jitter).
Vector3 patternOffset(const ArrayObject& a, int n) {
    const int nx = std::max(1, a.count[0]), ny = std::max(1, a.count[1]);
    switch (a.pattern) {
    case ArrayPattern::Line: return a.step * float(n);
    case ArrayPattern::Grid: {
        const int i = n % nx, j = (n / nx) % ny, k = n / (nx * ny);
        return {a.step.x * float(i), a.step.y * float(j), a.step.z * float(k)};
    }
    case ArrayPattern::Circle: {
        const float angle = 2.0f * kPi * float(n) / float(nx);
        return {a.radius * std::cos(angle), 0.0f, a.radius * std::sin(angle)};
    }
    }
    return Vector3(0.0f);
}

int copyCount(const ArrayObject& a) {
    const int nx = std::max(0, a.count[0]);
    if (a.pattern != ArrayPattern::Grid) return nx;
    return nx * std::max(0, a.count[1]) * std::max(0, a.count[2]);
}

} // namespace

const SceneObject* findObject(const SceneGraph& g, uint32_t id) {
    if (id == 0) return nullptr;
    if (const SceneObject* o = findIn(g.entities, id)) return o;
    if (const SceneObject* o = findIn(g.groups, id)) return o;
    if (const SceneObject* o = findIn(g.arrays, id)) return o;
    if (const SceneObject* o = findIn(g.lights, id)) return o;
    return findIn(g.cameras, id);
}

Quaternion objectRotation(const SceneObject& o) {
    const Vector3& d = o.rotationDeg;
    return Quaternion::fromAxisAngle({0, 1, 0}, degToRad(d.y)) * Quaternion::fromAxisAngle({1, 0, 0}, degToRad(d.x)) *
           Quaternion::fromAxisAngle({0, 0, 1}, degToRad(d.z));
}

// Euler angles in degrees of a rotation R = Ry(y) Rx(x) Rz(z) - the inverse of objectRotation().
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

void worldPose(const SceneGraph& g, const SceneObject& o, Vector3& position, Quaternion& rotation) {
    position = o.position;
    rotation = objectRotation(o);
    for (const SceneObject* p : ancestors(g, o)) { // nearest parent first: wrap outward
        const Quaternion pr = objectRotation(*p);
        position = p->position + pr.rotate(position);
        rotation = (pr * rotation).normalized();
    }
}

void setWorldPose(const SceneGraph& g, SceneObject& o, const Vector3& position, const Quaternion& rotation) {
    const std::vector<const SceneObject*> chain = ancestors(g, o);
    if (chain.empty()) { // a root: its pose is the world pose
        o.position = position;
        o.rotationDeg = eulerDegrees(rotation);
        return;
    }
    Vector3 pp;
    Quaternion pq;
    worldPose(g, *chain.front(), pp, pq);
    const Quaternion inv = pq.conjugate();
    o.position = inv.rotate(position - pp);
    o.rotationDeg = eulerDegrees((inv * rotation).normalized());
}

bool effectivelyVisible(const SceneGraph& g, const SceneObject& o) {
    if (!o.visible) return false;
    for (const SceneObject* p : ancestors(g, o))
        if (!p->visible) return false;
    return true;
}

uint32_t gluedGroupOf(const SceneGraph& g, const SceneObject& o) {
    for (const SceneObject* p : ancestors(g, o))
        for (const Group& gr : g.groups)
            if (gr.id == p->id && gr.glued) return gr.id;
    return 0;
}

Entity resolveInstance(const SceneGraph& g, const Entity& e) {
    const Entity* master = &e;
    std::vector<uint32_t> seen{e.id};
    while (master->instanceOf != 0) {
        const Entity* next = findEntity(g, master->instanceOf);
        if (!next || std::find(seen.begin(), seen.end(), next->id) != seen.end()) break;
        seen.push_back(next->id);
        master = next;
    }
    Entity out = *master;                                             // geometry and roles of the master
    static_cast<SceneObject&>(out) = static_cast<const SceneObject&>(e); // name, pose, visibility, id: its own
    out.instanceOf = e.instanceOf;
    return out;
}

uint32_t arrayCopyId(uint32_t arrayId, int copy) {
    return 0x80000000u | ((arrayId & 0x7FFFu) << 16) | (uint32_t(copy) & 0xFFFFu);
}
bool isArrayCopyId(uint32_t id) { return (id & 0x80000000u) != 0; }
uint32_t arrayOfCopyId(uint32_t id) { return isArrayCopyId(id) ? (id >> 16) & 0x7FFFu : 0; }

// Copy n: at the pattern's place (plus its jitter) in the array's frame, turned by n rotation steps
// and then by the template's own rotation; the template's own position is not used.
std::vector<Entity> expandArray(const SceneGraph& g, const ArrayObject& a) {
    std::vector<Entity> out;
    const Entity* tmpl = findEntity(g, a.templateId);
    if (!tmpl) return out;
    const Entity base = resolveInstance(g, *tmpl);
    Vector3 ap;
    Quaternion aq;
    worldPose(g, a, ap, aq);
    const bool visible = effectivelyVisible(g, a);
    const int n = copyCount(a);
    out.reserve(size_t(std::max(0, n)));
    for (int i = 0; i < n; ++i) {
        uint32_t state = a.seed * 2654435761u + uint32_t(i) * 40503u + 1u;
        Vector3 offset = patternOffset(a, i);
        if (a.jitter > 0) offset += Vector3(jitterNumber(state), jitterNumber(state), jitterNumber(state)) * a.jitter;
        SceneObject turn;
        turn.rotationDeg = a.rotationStepDeg * float(i);
        Entity c = base;
        c.id = arrayCopyId(a.id, i);
        c.name = a.name + "[" + std::to_string(i) + "]";
        c.parent = 0;
        c.instanceOf = 0;
        c.visible = visible;
        c.locked = a.locked;
        c.position = ap + aq.rotate(offset);
        c.rotationDeg = eulerDegrees((aq * objectRotation(turn) * objectRotation(base)).normalized());
        out.push_back(std::move(c));
    }
    return out;
}

// A root keeps its pose exactly as written (no round trip through a quaternion), so a flat scene
// builds bit for bit as it did before hierarchies existed.
std::vector<Entity> worldEntities(const SceneGraph& g) {
    std::vector<Entity> out;
    out.reserve(g.entities.size());
    for (const Entity& e : g.entities) {
        Entity w = resolveInstance(g, e);
        if (e.parent != 0 && findObject(g, e.parent)) {
            Vector3 p;
            Quaternion q;
            worldPose(g, e, p, q);
            w.position = p;
            w.rotationDeg = eulerDegrees(q);
        }
        w.parent = 0;
        w.visible = effectivelyVisible(g, e);
        out.push_back(std::move(w));
    }
    for (const ArrayObject& a : g.arrays) {
        std::vector<Entity> copies = expandArray(g, a);
        out.insert(out.end(), copies.begin(), copies.end());
    }
    return out;
}

uint32_t nextId(const SceneGraph& g) {
    uint32_t next = 1;
    for (const Entity& e : g.entities) next = std::max(next, e.id + 1);
    for (const Group& gr : g.groups) next = std::max(next, gr.id + 1);
    for (const ArrayObject& a : g.arrays) next = std::max(next, a.id + 1);
    for (const Light& l : g.lights) next = std::max(next, l.id + 1);
    for (const Camera& c : g.cameras) next = std::max(next, c.id + 1);
    return next;
}

Vector3 lightDirection(const SceneGraph& g, const Light& l) {
    Vector3 position;
    Quaternion rotation;
    worldPose(g, l, position, rotation);
    return normalize(rotation.rotate(Vector3(0, -1, 0)));
}

void cameraFrame(const SceneGraph& g, const Camera& c, Vector3& eye, Vector3& forward, Vector3& up) {
    Quaternion rotation;
    worldPose(g, c, eye, rotation);
    forward = normalize(rotation.rotate(Vector3(0, 0, -1)));
    up = normalize(rotation.rotate(Vector3(0, 1, 0)));
}

} // namespace rf
