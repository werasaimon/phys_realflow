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
//   entity "Curtain"
//     object id 4 visible 1 locked 0
//     shape plane size 1 0.02 1.2 position 0 1.2 0 rotation 90 0 0 color 0.9 0.85 0.7
//     cloth areaDensity 0.2 bendCompliance 0.001 tearable 1 pinned 16
//     flammable
//   end
//
// A role is written only when it is on. Numbers are written as the shortest decimal that reads
// back as the same float, so save -> load -> save gives the same text.
#include "scene/SceneGraph.h"

#include <cstdio>
#include <cstdlib>
#include <sstream>

namespace rf {

namespace {

const char* kShapeNames[] = {"box", "sphere", "cylinder", "cone", "plane"};

// The shortest decimal that reads back as the same float: 0.02 stays "0.02" (not the exact
// "0.0199999996"), yet save -> load -> save is still exact. Nine digits always suffice for a float.
std::string num(float v) {
    char buf[32];
    for (int digits = 6; digits <= 9; ++digits) {
        std::snprintf(buf, sizeof(buf), "%.*g", digits, double(v));
        if (std::strtof(buf, nullptr) == v) break;
    }
    return buf;
}

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
        char* end = nullptr;
        out = std::strtof(w_[i_].c_str(), &end);
        if (end == w_[i_].c_str() || *end) return false;
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

static void saveEntity(std::ostringstream& o, const Entity& e) {
    std::string name = e.name;
    for (char& c : name)
        if (c == '"') c = '\''; // the name is written in double quotes
    o << "entity \"" << name << "\"\n";
    o << "  object id " << e.id << " visible " << (e.visible ? 1 : 0) << " locked " << (e.locked ? 1 : 0) << "\n";
    o << "  shape " << kShapeNames[int(e.shape)] << " size " << vec(e.size) << " position " << vec(e.position)
      << " rotation " << vec(e.rotationDeg) << " color " << vec(e.color) << "\n";
    if (e.rigid.enabled)
        o << "  rigid density " << num(e.rigid.density) << " friction " << num(e.rigid.friction) << " restitution "
          << num(e.rigid.restitution) << " fixed " << (e.rigid.fixed ? 1 : 0) << " velocity " << vec(e.rigid.velocity)
          << " spin " << vec(e.rigid.angularVelocity) << "\n";
    if (e.soft.enabled) o << "  soft density " << num(e.soft.density) << " stiffness " << num(e.soft.stiffness) << "\n";
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
    for (int s = 0; s < 5; ++s)
        if (shapeName == kShapeNames[s]) e.shape = ShapeKind(s), known = true;
    if (!known) return fail(bad, shapeName);
    while (!r.done()) {
        const std::string k = r.key();
        const bool ok = k == "size" ? r.vector(e.size) : k == "position" ? r.vector(e.position)
                      : k == "rotation" ? r.vector(e.rotationDeg) : k == "color" ? r.vector(e.color) : false;
        if (!ok) return fail(bad, k);
    }
    return true;
}

// One "key value(s)" of a role line; false if the key is not one of that role's.
static bool readRoleKey(const std::string& role, const std::string& k, LineReader& r, Entity& e) {
    if (role == "object")
        return k == "id" ? r.integer(e.id) : k == "visible" ? r.flag(e.visible) : k == "locked" ? r.flag(e.locked) : false;
    if (role == "rigid")
        return k == "density" ? r.number(e.rigid.density) : k == "friction" ? r.number(e.rigid.friction)
             : k == "restitution" ? r.number(e.rigid.restitution) : k == "fixed" ? r.flag(e.rigid.fixed)
             : k == "velocity" ? r.vector(e.rigid.velocity) : k == "spin" ? r.vector(e.rigid.angularVelocity) : false;
    if (role == "soft") return k == "density" ? r.number(e.soft.density) : k == "stiffness" ? r.number(e.soft.stiffness) : false;
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

bool SceneGraph::load(const std::string& text, std::string& error) {
    SceneGraph g;
    std::istringstream in(text);
    std::string line, bad;
    bool inEntity = false;
    for (int number = 1; std::getline(in, line); ++number) {
        const std::vector<std::string> w = words(line);
        if (w.empty() || w[0][0] == '#') continue;
        bool ok = true;
        if (w[0] == "world" && !inEntity) {
            LineReader r(w, 1);
            ok = readWorld(r, g.world, bad);
        } else if (w[0] == "entity" && !inEntity) {
            g.entities.emplace_back();
            g.entities.back().name = w.size() > 1 ? w[1] : std::string();
            inEntity = true;
        } else if (w[0] == "end" && inEntity) {
            inEntity = false;
        } else if (inEntity) {
            ok = readEntityLine(w, g.entities.back(), bad);
        } else {
            ok = fail(bad, w[0]);
        }
        if (!ok) {
            error = "line " + std::to_string(number) + ": cannot read '" + bad + "'";
            return false;
        }
    }
    if (inEntity) {
        error = "the last entity has no 'end'";
        return false;
    }
    *this = std::move(g);
    return true;
}

} // namespace rf
