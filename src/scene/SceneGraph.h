#pragma once
// The scene as data: a list of entities, each a shape with a pose and a set of roles ("you are a
// rigid body", "you are a magnet", "you are liquid"). An editor builds it with buttons and a node
// graph, saves it as a small text file anyone can read, and GraphScene (below) turns it into a
// running Simulation - the same way a hand-written sample scene does, so everything the engine can
// do is reachable from the editor without code.
//
// A role is plain data here; what it means physically is in GraphScene.cpp (which solver it goes
// to) and, for the magnet, in Magnets.cpp (dipole-dipole forces between bodies).
#include "math/Math.h"
#include "scene/Scene.h"

#include <cstdint>
#include <string>
#include <vector>

namespace rf {

enum class ShapeKind { Box, Sphere, Cylinder, Cone, Plane };

// The roles an entity can have - the "Midas touch" of the editor: any shape can be turned into one
// physical model or several interacting ones. Two kinds:
//   * what the shape IS MADE OF - exactly one: rigid, soft, liquid or cloth (none: only drawn);
//   * what it also DOES - any number: magnet, emitter (smoke, heat, liquid), flammable.
// A rigid body that is a magnet and trails smoke is three roles on one entity; the solvers meet in
// the Simulation's frame step (the gas pushes bodies and cloth, bodies are walls for the gas,
// magnets pull each other, a flame heats the cloth it touches).
struct RigidRole {
    bool enabled = false;
    float density = 500;        // kg/m^3
    float friction = 0.5f;      // kinetic friction coefficient
    float restitution = 0.2f;   // bounciness, 0..1
    bool fixed = false;         // does not move (a wall, a table)
    Vector3 velocity{0.0f};     // initial, m/s
    Vector3 angularVelocity{0.0f};
};
struct SoftRole {
    bool enabled = false;
    float density = 150;
    float stiffness = 0.3f;     // shape-matching stiffness per pass, 0..1
};
struct LiquidRole {
    bool enabled = false;       // the shape's box is filled with liquid particles
};
struct MagnetRole {
    bool enabled = false;
    Vector3 moment{0, 1, 0};    // magnetic dipole moment in the body frame, A m^2
};
struct HeatRole {
    bool enabled = false;       // a hot spot in the gas at the shape (needs gas)
    float temperature = 300;    // kelvin above ambient
    float smoke = 1;            // smoke released per second (relative)
};
// A sheet of cloth in place of the shape: a Plane becomes a cloth of size.x x size.z; for another
// shape the cloth is its top face. Pinned edges hang it (a curtain on a rod, a hammock).
struct ClothRole {
    bool enabled = false;
    float areaDensity = 0.3f;   // kg/m^2 (cotton 0.2-0.3, canvas 1.5, silk 0.04)
    float bendCompliance = 1e-3f; // 0 = stiff as card, larger = drapes more
    bool tearable = true;       // tears along the threads when pulled too hard
    int pinnedEdges = 0;        // bit mask: 1 -x edge, 2 +x edge, 4 -z edge, 8 +z edge, 16 the top row (rod)
};
// Something that leaves the shape every second and follows it wherever it moves: a smoke trail
// behind a flying box, a hot jet, a tap of water. Needs the gas for smoke and heat.
struct EmitterRole {
    bool enabled = false;
    float smoke = 1;            // smoke per second (relative)
    float temperature = 0;      // kelvin above ambient of what it releases (0: cold smoke)
    float liquid = 0;           // liquid particles per second (0: none)
    Vector3 velocity{0.0f};     // the speed it is released with, in the shape's frame, m/s
};
// The shape burns: it takes heat from the gas, ignites above its ignition temperature and gives
// fuel and heat back (the cotton of the fire sample). For cloth: charring and burning through.
struct FlammableRole {
    bool enabled = false;
};

// What every thing in a scene has, whatever it is (a shape with roles today; a camera, a light or
// a standalone emitter later): a name, where it is, how it is turned, its colour, and the editor's
// switches. The inspector shows this part first and always the same way; the rest depends on
// what the thing is.
struct SceneObject {
    uint32_t id = 0;            // stable handle (unique in the graph; the editor assigns it)
    std::string name;
    Vector3 position{0.0f};
    Vector3 rotationDeg{0.0f};  // Euler angles in degrees: x pitch, y yaw, z roll (applied Ry Rx Rz)
    Vector3 color{0.8f, 0.8f, 0.8f};
    bool visible = true;        // hidden things take no part in the simulation either
    bool locked = false;        // the editor does not select or move it (a floor, a backdrop)
};

// A shape in the scene with its roles: what it is made of and what else it does.
struct Entity : SceneObject {
    ShapeKind shape = ShapeKind::Box;
    Vector3 size{0.2f};         // full size: box edges; sphere/cylinder/cone: diameter, height in y
    RigidRole rigid;
    SoftRole soft;
    LiquidRole liquid;
    MagnetRole magnet;
    HeatRole heat;
    ClothRole cloth;
    EmitterRole emitter;
    FlammableRole flammable;
};

struct WorldSettings {
    Vector3 gravity{0, -9.81f, 0};
    Vector3 size{4, 3, 4};      // the box everything lives in, m (floor at y = 0)
    bool gas = false;           // air on a grid: smoke, heat, drag on bodies
    bool magneticGas = false;   // the gas is a plasma (MHD): magnets push it and it pushes back
};

struct SceneGraph {
    WorldSettings world;
    std::vector<Entity> entities;

    // Text format, one line per fact, readable and diff-able (see SceneGraph.cpp for a sample).
    std::string save() const;
    bool load(const std::string& text, std::string& error);
};

// A Scene built from a SceneGraph: the editor's scenes run through the same Scene interface as the
// samples (configure the world, build the entities, apply the magnet forces every step).
class GraphScene : public Scene {
public:
    explicit GraphScene(SceneGraph graph) : graph_(std::move(graph)) {}
    void configure(Simulation& sim) override;
    void build(Simulation& sim) override;
    void afterStep(Simulation& sim) override;
    void describe(const Simulation& sim, RenderSnapshot& s) const override;
    const SceneGraph& graph() const { return graph_; }
    // Which entity made rigid body b (-1: none); filled by build(). The editor selects the entity
    // the mouse clicked on through it.
    int entityOfBody(int body) const { return body >= 0 && body < int(entityOfBody_.size()) ? entityOfBody_[body] : -1; }

private:
    // An emitter and what carries it: the rigid body it rides on (-1: it stays where it was put) and
    // the turn from that body's frame to the entity's own frame (a hull body lives in its principal
    // frame), so the release direction turns with the body.
    struct EmitterRef {
        int entity = -1;
        int body = -1;
        Matrix3x3 bodyToEntity = Matrix3x3::identity();
    };

    void addEntity(Simulation& sim, int index);
    int addRigidBody(Simulation& sim, const Entity& e, bool fixed);
    void addCloth(Simulation& sim, const Entity& e) const;
    void addEmitter(Simulation& sim, int index, int body);
    void releaseFromEmitter(Simulation& sim, const EmitterRef& ref, bool& liquidDone) const;
    void addMagnetFieldToGas(Simulation& sim) const;

    SceneGraph graph_;
    std::vector<int> magnetBody_;       // rigid body index of each magnet entity
    std::vector<Vector3> magnetMoment_; // its moment in the body frame
    float maxMagnetForce_ = 0;          // the largest pair force of the last frame [N]
    std::vector<int> entityOfBody_;     // body index -> entity index
    std::vector<EmitterRef> emitters_;  // every visible entity with the emitter role
    std::vector<std::string> notes_;    // roles that could not be honoured, shown in the readings
};

} // namespace rf
