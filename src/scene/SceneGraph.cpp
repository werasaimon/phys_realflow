// Saving and loading a SceneGraph as text: one line per fact, so a scene file can be read, edited
// by hand and compared with diff. A sample:
//
//   # PhysRealFlow scene
//   world gravity 0 -9.81 0 size 4 3 4 gas 0 magneticGas 0
//   entity "Magnet A"
//     object id 3 visible 1 locked 0
//     shape box size 0.1 0.05 0.05 position 0 0.5 0 rotation 0 0 0 color 0.8 0.2 0.2
//     rigid density 7800 friction 0.5 restitution 0.2 fixed 0 velocity 0 0 0 spin 0 0 0
//     magnet moment 1 0 0
//     emitter smoke 1 temperature 0 liquid 0 velocity 0 0.5 0
//   end
//   entity "Horse"
//     shape mesh file "models/horse.obj" size 0.5 0.5 0.5 position 1 0.3 0 rotation 0 90 0 color 0.6 0.4 0.3
//     rigid density 600 friction 0.5 restitution 0.2 fixed 0 velocity 0 0 0 spin 0 0 0
//     collider kind capsule fit 1 size 0.2 0.2 0.2 offset 0 0 0 rotation 0 0 90
//   end
//   entity "Curtain"
//     object id 4 visible 1 locked 0
//     shape plane size 1 0.02 1.2 position 0 1.2 0 rotation 90 0 0 color 0.9 0.85 0.7
//     cloth areaDensity 0.2 bendCompliance 0.001 tearable 1 pinned 16
//     flammable
//   end
//
// A role is written only when it is on - the collider too, a role of its own (kinds: auto, box,
// sphere, capsule, hull, decomposition); a rigid entity without a collider line (a file from before
// colliders were a role) loads with the Auto collider. Numbers are written as the shortest decimal that reads
// back as the same float, so save -> load -> save gives the same text.
#include "scene/SceneGraph.h"

#include "core/Format.h"

#include <cstdlib>
#include <sstream>

namespace rf {

namespace {

const char* kShapeNames[] = {"box", "sphere", "cylinder", "cone", "plane", "mesh"};
const int kShapeCount = 6;
const char* kColliderNames[] = {"auto", "box", "sphere", "capsule", "hull", "decomposition"};
const int kColliderCount = 6;

std::string num(float v) { return numberText(v); }
std::string vec(const Vector3& v) { return num(v.x) + " " + num(v.y) + " " + num(v.z); }

// Splits a line into words; a word in double quotes may contain spaces.
std::vector<std::string> words(const std::string& line) {
    std::vector<std::string> out;
    size_t i = 0;
    while (i < line.size()) {
        while (i < line.size() && (line[i] == ' ' || line[i] == '\t' || line[i] == '\r')) ++i;
        if (i >= line.size()) break;
        if (line[i] == '"') {
            const size_t end = line.find('"', i + 1);
            out.push_back(line.substr(i + 1, end == std::string::npos ? std::string::npos : end - i - 1));
            i = end == std::string::npos ? line.size() : end + 1;
        } else {
            const size_t start = i;
            while (i < line.size() && line[i] != ' ' && line[i] != '\t' && line[i] != '\r') ++i;
            out.push_back(line.substr(start, i - start));
        }
    }
    return out;
}

// Reads the words of one line as "key value(s)" pairs. `need(k)` says how many numbers follow key k.
class LineReader {
public:
    LineReader(const std::vector<std::string>& w, size_t start) : w_(w), i_(start) {}
    bool done() const { return i_ >= w_.size(); }
    const std::string& key() { return w_[i_++]; }
    bool number(float& out) {
        if (i_ >= w_.size()) return false;
        if (!parseNumber(w_[i_], out)) return false;
        ++i_;
        return true;
    }
    bool vector(Vector3& v) { return number(v.x) && number(v.y) && number(v.z); }
    // A whole number: ids and bit masks stay exact (a float holds integers only up to 2^24).
    bool integer(uint32_t& out) {
        if (i_ >= w_.size()) return false;
        char* end = nullptr;
        const unsigned long v = std::strtoul(w_[i_].c_str(), &end, 10);
        if (end == w_[i_].c_str() || *end) return false;
        out = uint32_t(v);
        ++i_;
        return true;
    }
    bool integer(int& out) {
        uint32_t u;
        if (!integer(u)) return false;
        out = int(u);
        return true;
    }
    // A word as it is (a quoted file path arrives without its quotes).
    bool text(std::string& out) {
        if (i_ >= w_.size()) return false;
        out = w_[i_++];
        return true;
    }
    bool flag(bool& b) {
        float f;
        if (!number(f)) return false;
        b = f != 0;
        return true;
    }

private:
    const std::vector<std::string>& w_;
    size_t i_;
};

// Remembers what could not be read and says "failed".
bool fail(std::string& bad, const std::string& what) {
    bad = what;
    return false;
}

} // namespace

// ---------------------------------------------------------------------------
// Saving
// ---------------------------------------------------------------------------
// The roles added for the "Midas touch" of the editor: cloth, emitter, flammable.
static void saveMoreRoles(std::ostringstream& o, const Entity& e) {
    if (e.cloth.enabled)
        o << "  cloth areaDensity " << num(e.cloth.areaDensity) << " bendCompliance " << num(e.cloth.bendCompliance)
          << " tearable " << (e.cloth.tearable ? 1 : 0) << " pinned " << e.cloth.pinnedEdges << "\n";
    if (e.emitter.enabled)
        o << "  emitter smoke " << num(e.emitter.smoke) << " temperature " << num(e.emitter.temperature) << " liquid "
          << num(e.emitter.liquid) << " velocity " << vec(e.emitter.velocity) << "\n";
    if (e.flammable.enabled) o << "  flammable\n";
}

// The collider role, a line of its own with every field (written only when the role is on).
static void saveCollider(std::ostringstream& o, const ColliderRole& c) {
    o << "  collider kind " << kColliderNames[int(c.kind)] << " fit " << (c.fitToGeometry ? 1 : 0) << " size " << vec(c.size)
      << " offset " << vec(c.offset) << " rotation " << vec(c.rotationDeg) << "\n";
}

// The part every object has, for groups and arrays: a quoted name, its object line and its pose.
static void saveObjectHead(std::ostringstream& o, const char* kind, const SceneObject& s) {
    std::string name = s.name;
    for (char& c : name)
        if (c == '"') c = '\'';
    o << kind << " \"" << name << "\"\n";
    o << "  object id " << s.id << " visible " << (s.visible ? 1 : 0) << " locked " << (s.locked ? 1 : 0);
    if (s.parent != 0) o << " parent " << s.parent;
    o << "\n  pose position " << vec(s.position) << " rotation " << vec(s.rotationDeg) << " color " << vec(s.color) << "\n";
}

static void saveGroup(std::ostringstream& o, const Group& g) {
    saveObjectHead(o, "group", g);
    if (g.glued) o << "  glued\n";
    o << "end\n";
}

static const char* kPatternNames[] = {"line", "grid", "circle"};

static void saveArray(std::ostringstream& o, const ArrayObject& a) {
    saveObjectHead(o, "array", a);
    o << "  pattern " << kPatternNames[int(a.pattern)] << " template " << a.templateId << " count " << a.count[0] << " "
      << a.count[1] << " " << a.count[2] << " step " << vec(a.step) << " radius " << num(a.radius) << " rotationStep "
      << vec(a.rotationStepDeg) << " jitter " << num(a.jitter) << " seed " << a.seed << "\n";
    o << "end\n";
}

static const char* kLightNames[] = {"sun", "point", "spot"};

// Lights and cameras: only their own settings follow the object head (the physics ignores them).
static void saveLight(std::ostringstream& o, const Light& l) {
    saveObjectHead(o, "light", l);
    o << "  kind " << kLightNames[int(l.kind)] << " intensity " << num(l.intensity) << " range " << num(l.range) << " cone "
      << num(l.coneDeg) << " softness " << num(l.softnessDeg) << " shadows " << (l.shadows ? 1 : 0) << "\n";
    o << "end\n";
}

static void saveCamera(std::ostringstream& o, const Camera& c) {
    saveObjectHead(o, "camera", c);
    o << "  lens fov " << num(c.fovDeg) << " near " << num(c.nearClip) << " far " << num(c.farClip) << " active "
      << (c.active ? 1 : 0) << "\n";
    o << "end\n";
}

static void saveEntity(std::ostringstream& o, const Entity& e) {
    std::string name = e.name;
    for (char& c : name)
        if (c == '"') c = '\''; // the name is written in double quotes
    o << "entity \"" << name << "\"\n";
    o << "  object id " << e.id << " visible " << (e.visible ? 1 : 0) << " locked " << (e.locked ? 1 : 0);
    if (e.parent != 0) o << " parent " << e.parent;
    if (e.instanceOf != 0) o << " instance " << e.instanceOf;
    o << "\n";
    o << "  shape " << kShapeNames[int(e.shape)];
    if (e.shape == ShapeKind::Mesh) o << " file \"" << e.meshFile << "\"";
    o << " size " << vec(e.size) << " position " << vec(e.position)
      << " rotation " << vec(e.rotationDeg) << " color " << vec(e.color) << "\n";
    if (e.rigid.enabled)
        o << "  rigid density " << num(e.rigid.density) << " friction " << num(e.rigid.friction) << " restitution "
          << num(e.rigid.restitution) << " fixed " << (e.rigid.fixed ? 1 : 0) << " velocity " << vec(e.rigid.velocity)
          << " spin " << vec(e.rigid.angularVelocity) << "\n";
    if (e.collider.enabled) saveCollider(o, e.collider);
    if (e.soft.enabled)
        o << "  soft density " << num(e.soft.material.density) << " young " << num(e.soft.material.youngModulus) << " poisson "
          << num(e.soft.material.poissonRatio) << " damping " << num(e.soft.material.damping) << "\n";
    if (e.liquid.enabled) o << "  liquid\n";
    if (e.magnet.enabled) o << "  magnet moment " << vec(e.magnet.moment) << "\n";
    if (e.heat.enabled) o << "  heat temperature " << num(e.heat.temperature) << " smoke " << num(e.heat.smoke) << "\n";
    saveMoreRoles(o, e);
    o << "end\n";
}

std::string SceneGraph::save() const {
    std::ostringstream o;
    o << "# PhysRealFlow scene\n";
    o << "world gravity " << vec(world.gravity) << " size " << vec(world.size) << " gas " << (world.gas ? 1 : 0)
      << " magneticGas " << (world.magneticGas ? 1 : 0) << "\n";
    for (const Entity& e : entities) saveEntity(o, e);
    for (const Group& g : groups) saveGroup(o, g);
    for (const ArrayObject& a : arrays) saveArray(o, a);
    for (const Light& l : lights) saveLight(o, l);
    for (const Camera& c : cameras) saveCamera(o, c);
    return o.str();
}

// ---------------------------------------------------------------------------
// Loading
// ---------------------------------------------------------------------------
static bool readWorld(LineReader& r, WorldSettings& w, std::string& bad) {
    while (!r.done()) {
        const std::string k = r.key();
        const bool ok = k == "gravity" ? r.vector(w.gravity) : k == "size" ? r.vector(w.size) : k == "gas" ? r.flag(w.gas)
                      : k == "magneticGas" ? r.flag(w.magneticGas) : false;
        if (!ok) return fail(bad, k);
    }
    return true;
}

static bool readShape(LineReader& r, Entity& e, const std::string& shapeName, std::string& bad) {
    bool known = false;
    for (int s = 0; s < kShapeCount; ++s)
        if (shapeName == kShapeNames[s]) e.shape = ShapeKind(s), known = true;
    if (!known) return fail(bad, shapeName);
    while (!r.done()) {
        const std::string k = r.key();
        const bool ok = k == "size" ? r.vector(e.size) : k == "position" ? r.vector(e.position)
                      : k == "rotation" ? r.vector(e.rotationDeg) : k == "color" ? r.vector(e.color)
                      : k == "file" ? r.text(e.meshFile) : false;
        if (!ok) return fail(bad, k);
    }
    return true;
}

// One "key value(s)" of a role line; false if the key is not one of that role's.
static bool readRoleKey(const std::string& role, const std::string& k, LineReader& r, Entity& e) {
    if (role == "object") {
        if (k == "instance") return r.integer(e.instanceOf);
        return k == "id" ? r.integer(e.id) : k == "visible" ? r.flag(e.visible) : k == "locked" ? r.flag(e.locked)
             : k == "parent" ? r.integer(e.parent) : false;
    }
    if (role == "rigid")
        return k == "density" ? r.number(e.rigid.density) : k == "friction" ? r.number(e.rigid.friction)
             : k == "restitution" ? r.number(e.rigid.restitution) : k == "fixed" ? r.flag(e.rigid.fixed)
             : k == "velocity" ? r.vector(e.rigid.velocity) : k == "spin" ? r.vector(e.rigid.angularVelocity) : false;
    if (role == "collider") {
        ColliderRole& c = e.collider;
        if (k != "kind")
            return k == "fit" ? r.flag(c.fitToGeometry) : k == "size" ? r.vector(c.size) : k == "offset" ? r.vector(c.offset)
                 : k == "rotation" ? r.vector(c.rotationDeg) : false;
        std::string kind;
        if (!r.text(kind)) return false;
        for (int i = 0; i < kColliderCount; ++i)
            if (kind == kColliderNames[i]) { c.kind = ColliderKind(i); return true; }
        return false;
    }
    if (role == "soft") {
        SoftMaterial& m = e.soft.material;
        if (k == "stiffness") { // a file from the shape-matching soft bodies: its 0..1 as a modulus, 10 kPa .. 1 MPa
            float stiffness = 0;
            if (!r.number(stiffness)) return false;
            m.youngModulus = std::pow(10.0f, 4.0f + 2.0f * std::clamp(stiffness, 0.0f, 1.0f));
            return true;
        }
        return k == "density" ? r.number(m.density) : k == "young" ? r.number(m.youngModulus) : k == "poisson" ? r.number(m.poissonRatio)
             : k == "damping" ? r.number(m.damping) : false;
    }
    if (role == "magnet") return k == "moment" && r.vector(e.magnet.moment);
    if (role == "heat") return k == "temperature" ? r.number(e.heat.temperature) : k == "smoke" ? r.number(e.heat.smoke) : false;
    if (role == "cloth")
        return k == "areaDensity" ? r.number(e.cloth.areaDensity) : k == "bendCompliance" ? r.number(e.cloth.bendCompliance)
             : k == "tearable" ? r.flag(e.cloth.tearable) : k == "pinned" ? r.integer(e.cloth.pinnedEdges) : false;
    if (role == "emitter")
        return k == "smoke" ? r.number(e.emitter.smoke) : k == "temperature" ? r.number(e.emitter.temperature)
             : k == "liquid" ? r.number(e.emitter.liquid) : k == "velocity" ? r.vector(e.emitter.velocity) : false;
    return false; // liquid and flammable have no keys
}

static bool readRole(const std::string& role, LineReader& r, Entity& e, std::string& bad) {
    while (!r.done()) {
        const std::string k = r.key();
        if (!readRoleKey(role, k, r, e)) return fail(bad, k);
    }
    return true;
}

// One line inside "entity ... end". Returns false with `bad` naming what was not understood.
static bool readEntityLine(const std::vector<std::string>& w, Entity& e, std::string& bad) {
    const std::string& what = w[0];
    if (what == "shape") {
        if (w.size() < 2) return fail(bad, what);
        LineReader r(w, 2); // "shape <kind> key values ..."
        return readShape(r, e, w[1], bad);
    }
    LineReader r(w, 1);
    if (what == "object") {} // the object's id and flags, not a role
    else if (what == "collider") e.collider.enabled = true;
    else if (what == "rigid") e.rigid.enabled = true;
    else if (what == "soft") e.soft.enabled = true;
    else if (what == "liquid") e.liquid.enabled = true;
    else if (what == "magnet") e.magnet.enabled = true;
    else if (what == "heat") e.heat.enabled = true;
    else if (what == "cloth") e.cloth.enabled = true;
    else if (what == "emitter") e.emitter.enabled = true;
    else if (what == "flammable") e.flammable.enabled = true;
    else return fail(bad, what);
    return readRole(what, r, e, bad);
}

// The object line and the pose line of a group or an array.
static bool readObjectLine(const std::vector<std::string>& w, SceneObject& s, std::string& bad) {
    LineReader r(w, 1);
    while (!r.done()) {
        const std::string k = r.key();
        const bool ok = w[0] == "object" ? (k == "id" ? r.integer(s.id) : k == "visible" ? r.flag(s.visible)
                                            : k == "locked" ? r.flag(s.locked) : k == "parent" ? r.integer(s.parent) : false)
                                         : (k == "position" ? r.vector(s.position) : k == "rotation" ? r.vector(s.rotationDeg)
                                            : k == "color" ? r.vector(s.color) : false);
        if (!ok) return fail(bad, k);
    }
    return true;
}

static bool readGroupLine(const std::vector<std::string>& w, Group& g, std::string& bad) {
    if (w[0] == "object" || w[0] == "pose") return readObjectLine(w, g, bad);
    if (w[0] == "glued" && w.size() == 1) return g.glued = true;
    return fail(bad, w[0]);
}

static bool readArrayLine(const std::vector<std::string>& w, ArrayObject& a, std::string& bad) {
    if (w[0] == "object" || w[0] == "pose") return readObjectLine(w, a, bad);
    if (w[0] != "pattern" || w.size() < 2) return fail(bad, w[0]);
    bool known = false;
    for (int p = 0; p < 3; ++p)
        if (w[1] == kPatternNames[p]) a.pattern = ArrayPattern(p), known = true;
    if (!known) return fail(bad, w[1]);
    LineReader r(w, 2);
    while (!r.done()) {
        const std::string k = r.key();
        const bool ok = k == "template" ? r.integer(a.templateId)
                      : k == "count" ? r.integer(a.count[0]) && r.integer(a.count[1]) && r.integer(a.count[2])
                      : k == "step" ? r.vector(a.step) : k == "radius" ? r.number(a.radius)
                      : k == "rotationStep" ? r.vector(a.rotationStepDeg) : k == "jitter" ? r.number(a.jitter)
                      : k == "seed" ? r.integer(a.seed) : false;
        if (!ok) return fail(bad, k);
    }
    return true;
}

static bool readLightLine(const std::vector<std::string>& w, Light& l, std::string& bad) {
    if (w[0] == "object" || w[0] == "pose") return readObjectLine(w, l, bad);
    if (w[0] != "kind" || w.size() < 2) return fail(bad, w[0]);
    bool known = false;
    for (int k = 0; k < 3; ++k)
        if (w[1] == kLightNames[k]) l.kind = LightKind(k), known = true;
    if (!known) return fail(bad, w[1]);
    LineReader r(w, 2);
    while (!r.done()) {
        const std::string k = r.key();
        const bool ok = k == "intensity" ? r.number(l.intensity) : k == "range" ? r.number(l.range)
                      : k == "cone" ? r.number(l.coneDeg) : k == "softness" ? r.number(l.softnessDeg)
                      : k == "shadows" ? r.flag(l.shadows) : false;
        if (!ok) return fail(bad, k);
    }
    return true;
}

static bool readCameraLine(const std::vector<std::string>& w, Camera& c, std::string& bad) {
    if (w[0] == "object" || w[0] == "pose") return readObjectLine(w, c, bad);
    if (w[0] != "lens") return fail(bad, w[0]);
    LineReader r(w, 1);
    while (!r.done()) {
        const std::string k = r.key();
        const bool ok = k == "fov" ? r.number(c.fovDeg) : k == "near" ? r.number(c.nearClip)
                      : k == "far" ? r.number(c.farClip) : k == "active" ? r.flag(c.active) : false;
        if (!ok) return fail(bad, k);
    }
    return true;
}

namespace {
enum class Block { None, Entity, Group, Array, Light, Camera };
}

// One line of a scene file in the block it belongs to (or starting / ending one).
static bool readLine(const std::vector<std::string>& w, SceneGraph& g, Block& block, std::string& bad) {
    if (block == Block::None) {
        const std::string name = w.size() > 1 ? w[1] : std::string();
        if (w[0] == "world") {
            LineReader r(w, 1);
            return readWorld(r, g.world, bad);
        }
        if (w[0] == "entity") g.entities.emplace_back(), g.entities.back().name = name, block = Block::Entity;
        else if (w[0] == "group") g.groups.emplace_back(), g.groups.back().name = name, block = Block::Group;
        else if (w[0] == "array") g.arrays.emplace_back(), g.arrays.back().name = name, block = Block::Array;
        else if (w[0] == "light") g.lights.emplace_back(), g.lights.back().name = name, block = Block::Light;
        else if (w[0] == "camera") g.cameras.emplace_back(), g.cameras.back().name = name, block = Block::Camera;
        else return fail(bad, w[0]);
        return true;
    }
    if (w[0] == "end") {
        // A file from before colliders were a role of their own: its rigid bodies collide as their
        // geometry - the Auto collider, now written out.
        if (block == Block::Entity && g.entities.back().rigid.enabled) g.entities.back().collider.enabled = true;
        block = Block::None;
        return true;
    }
    if (block == Block::Entity) return readEntityLine(w, g.entities.back(), bad);
    if (block == Block::Group) return readGroupLine(w, g.groups.back(), bad);
    if (block == Block::Light) return readLightLine(w, g.lights.back(), bad);
    if (block == Block::Camera) return readCameraLine(w, g.cameras.back(), bad);
    return readArrayLine(w, g.arrays.back(), bad);
}

bool SceneGraph::load(const std::string& text, std::string& error) {
    SceneGraph g;
    std::istringstream in(text);
    std::string line, bad;
    Block block = Block::None;
    for (int number = 1; std::getline(in, line); ++number) {
        const std::vector<std::string> w = words(line);
        if (w.empty() || w[0][0] == '#') continue;
        if (!readLine(w, g, block, bad)) {
            error = "line " + std::to_string(number) + ": cannot read '" + bad + "'";
            return false;
        }
    }
    if (block != Block::None) {
        error = "the last block has no 'end'";
        return false;
    }
    g.baseDirectory = baseDirectory; // where this graph's model files are: not in the text
    *this = std::move(g);
    return true;
}

} // namespace rf
