#pragma once
// The scene as data: a list of entities, each a shape with a pose and a set of roles ("you are a
// rigid body", "you are a magnet", "you are liquid"). An editor builds it with buttons and a node
// graph, saves it as a small text file anyone can read, and GraphScene (below) turns it into a
// running Simulation - the same way a hand-written sample scene does, so everything the engine can
// do is reachable from the editor without code.
//
// A role is plain data here; what it means physically is in GraphScene.cpp (which solver it goes
// to) and, for the magnet, in Magnets.cpp (dipole-dipole forces between bodies). The shape of an
// entity as a mesh - for the solvers and for a viewer drawing the scene - is in EntityShapes.cpp.
//
// Edit and play: an editor edits the graph (geometry and roles); Play builds a Simulation from it
// through GraphScene; Stop throws the simulation away and shows the graph again. Creating a shape
// makes geometry only: without roles it is drawn but takes no part in the simulation.
#include "core/Mesh.h"
#include "math/Math.h"
#include "scene/Scene.h"

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace rf {

class CompoundShape;
class ConvexShape;

// Mesh: a model read from an OBJ or STL file (Entity::meshFile).
enum class ShapeKind { Box, Sphere, Cylinder, Cone, Plane, Mesh };

// The roles an entity can have - the "Midas touch" of the editor: any shape can be turned into one
// physical model or several interacting ones. Two kinds:
//   * what the shape IS MADE OF - exactly one: rigid, soft, liquid or cloth (none: only drawn);
//   * what it also DOES - any number: magnet, emitter (smoke, heat, liquid), flammable.
// A rigid body that is a magnet and trails smoke is three roles on one entity; the solvers meet in
// the Simulation's frame step (the gas pushes bodies and cloth, bodies are walls for the gas,
// magnets pull each other, a flame heats the cloth it touches).
// What a body collides with - its own component on the entity, apart from what the entity looks
// like (Unity's Collider next to its MeshRenderer, Unreal's simple collision): a horse model may
// collide as a capsule. The body is always DRAWN with its geometry; the collider is only for
// contacts (the editor shows it as a thin wireframe). A collider alone makes a STATIC obstacle (a
// wall, a floor); with the rigid role it is the collider of a moving body; a rigid role without
// one falls back to Auto (with a note: the editor adds a collider together with the rigid role).
// Auto: the geometry itself - a box for a box, a sphere for a sphere, a thin box for a plane, a
// model's convex decomposition, the convex hull of anything else; offset and rotation are then
// ignored. The other kinds sit at `offset` / `rotationDeg` in the object's frame: Box, Sphere and
// Capsule sized by `size` or fitted to the geometry's extent along the collider's own axes;
// ConvexHull the hull of the geometry, Decomposition its convex parts (a model; a primitive is
// convex already and gets its hull). A body's mass is the collider's volume times the density, its
// centre of mass the collider's.
enum class ColliderKind { Auto, Box, Sphere, Capsule, ConvexHull, Decomposition };
struct ColliderRole {
    bool enabled = false;
    ColliderKind kind = ColliderKind::Auto;
    bool fitToGeometry = true;  // size from the geometry's extent (then `size` is ignored)
    Vector3 size{0.2f};         // box: edges; sphere: diameter = x; capsule: diameter = x, full height = y
    Vector3 offset{0.0f};       // collider centre relative to the object's centre, in the object's frame
    Vector3 rotationDeg{0.0f};  // collider orientation relative to the object (capsule axis = y)
};
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

// A shape in the scene with its roles: what it is made of and what else it does. With no role at
// all it is geometry only: drawn, but no body, no collision, no part in the simulation.
struct Entity : SceneObject {
    ShapeKind shape = ShapeKind::Box;
    // Full size: box edges; sphere/cylinder/cone: diameter (x) and height (y). A mesh keeps the
    // model's proportions: it is scaled uniformly so that its largest extent equals the largest
    // component of size (entityLocalMesh(e).bounds() is its real box).
    Vector3 size{0.2f};
    std::string meshFile;       // ShapeKind::Mesh: the OBJ/STL file as written in the scene file
    RigidRole rigid;
    ColliderRole collider;
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
    // Not saved: the folder relative mesh files are found in (the scene file's own folder). load()
    // keeps it as it was.
    std::string baseDirectory;

    // Text format, one line per fact, readable and diff-able (see SceneGraph.cpp for a sample).
    std::string save() const;
    bool load(const std::string& text, std::string& error);
};

// ---------------------------------------------------------------------------
// The shape of an entity as a mesh (EntityShapes.cpp) - the same mesh the solvers get, so a viewer
// can draw geometry-only entities and the whole scene while it is being edited.
// ---------------------------------------------------------------------------
// The entity's orientation: yaw about y, then pitch about x, then roll about z (degrees).
Quaternion entityRotation(const Entity& e);
// The shape as a closed mesh centred on the origin in the entity's own frame, at its size. A model
// file is read once and cached; an unreadable file gives an empty mesh (and `error`, if asked).
TriMesh entityLocalMesh(const Entity& e, const std::string& baseDirectory = std::string(), std::string* error = nullptr);
// The same mesh posed in the world: turned by the rotation and moved to the position.
TriMesh entityMesh(const Entity& e, const std::string& baseDirectory = std::string());
// No role at all: drawn, but not part of the simulation.
bool entityIsGeometryOnly(const Entity& e);
// A model as a rigid body: its convex decomposition at the entity's size, computed once per file
// and size and shared (decomposing takes seconds). Null for an unreadable file.
std::shared_ptr<const CompoundShape> entityCompound(const Entity& e, const std::string& baseDirectory);

// The collider of an entity (Entity::collider; Auto when it has none) as a shape, and where the shape's
// own frame - its centre of mass and principal axes - sits in the entity's frame. The body is
// created there: pos = entity position + entity rotation * position, rot = entity rotation * rotation.
struct EntityCollider {
    std::shared_ptr<const ConvexShape> shape; // null: the geometry could not be read
    Vector3 position{0.0f};
    Quaternion rotation;
};
EntityCollider entityCollider(const Entity& e, const std::string& baseDirectory = std::string());
// The collider as a closed mesh posed in the world (every part merged), for an editor's wireframe.
TriMesh colliderMesh(const Entity& e, const std::string& baseDirectory = std::string());

// A meta-object: one physical incarnation a role of an entity created in the solvers. The entity
// (the object) keeps its source and its roles; its meta-objects can be removed and built again from
// the source one entity at a time, without touching the rest of the scene.
struct MetaObject {
    enum class Kind { RigidBody, SoftBody, Liquid, Cloth, Magnet, Emitter };
    Kind kind = Kind::RigidBody;
    // RigidBody: body index; SoftBody / Liquid / Cloth: the particle group; Magnet / Emitter: the
    // rigid body it rides on (-1: none).
    int handle = -1;
    // RigidBody: the turn from the body's frame to the entity's frame (a hull body lives in its
    // principal frame), to read the entity's orientation back from the body.
    Quaternion bodyToEntity;
    // RigidBody: where the body's centre sits in the entity's frame (a collider with an offset, a
    // cone's centre of mass): the entity's position = body position - entity rotation * this.
    Vector3 bodyOffset{0.0f};
};

// A Scene built from a SceneGraph: the editor's scenes run through the same Scene interface as the
// samples (configure the world, build the entities, apply the magnet forces every step).
//
// Three layers: the SOURCE (the shape or model, EntityShapes.cpp - physics never changes it), the
// OBJECT (an Entity: pose and roles) and its META-OBJECTS (what the roles made in the solvers).
// rebuildEntity / addEntity / removeEntity change one entity between frames: its meta-objects are
// removed and built again from the source, the rest of the scene goes on untouched. The world
// itself (gravity, box size, gas on or off, the gas's heat source and combustion, a plasma's
// background field) is set when the scene is loaded; changing it needs a full reload.
class GraphScene : public Scene {
public:
    // Every entity gets a stable id: those without one (id 0) are numbered after the largest.
    explicit GraphScene(SceneGraph graph);
    void configure(Simulation& sim) override;
    void build(Simulation& sim) override;
    void afterStep(Simulation& sim) override;
    void describe(const Simulation& sim, RenderSnapshot& s) const override;
    const SceneGraph& graph() const { return graph_; }
    // Which entity made rigid body b: its index in graph().entities (-1: none) or its id (0: none).
    // The editor selects the entity the mouse clicked on through it.
    int entityOfBody(int body) const;
    uint32_t entityIdOfBody(int body) const;

    // One entity between frames, without reloading the scene. rebuildEntity replaces the entity of
    // updated.id: its meta-objects go, new ones are built from the source where the entity is NOW
    // (its body's or particles' live pose, unless `updated` moves it) and carry on its velocity.
    // Returns false when a role of it needs something only a reload can switch on (the gas).
    bool rebuildEntity(Simulation& sim, const Entity& updated);
    bool addEntity(Simulation& sim, const Entity& entity); // an id of 0 gets a new one
    void removeEntity(Simulation& sim, uint32_t id);
    const std::vector<MetaObject>& metaObjects(uint32_t id) const;

private:
    // An emitter and what carries it: the rigid body it rides on (-1: it stays where it was put) and
    // the turn from that body's frame to the entity's own frame (a hull body lives in its principal
    // frame), so the release direction turns with the body.
    struct EmitterRef {
        uint32_t entity = 0; // the entity's id
        int body = -1;
        Matrix3x3 bodyToEntity = Matrix3x3::identity();
        Vector3 bodyOffset{0.0f}; // the body's centre in the entity's frame (see MetaObject)
    };
    // Where an entity is and how it moves, read from its meta-objects (the live state).
    struct LiveState {
        Vector3 position{0.0f}, rotationDeg{0.0f}, velocity{0.0f}, angularVelocity{0.0f};
        bool fromParticles = false; // read from a particle group: the lowest particle centre below
        float lowest = 0;
    };
    void standOnParticles(Entity& next, const LiveState& live) const;

    int indexOf(uint32_t id) const;
    bool shapeUsable(const Entity& e);
    void buildEntity(Simulation& sim, int index);
    void destroyMetaObjects(Simulation& sim, uint32_t id);
    LiveState liveState(const Simulation& sim, int index) const;
    void carryVelocity(Simulation& sim, uint32_t id, const LiveState& live) const;
    bool needsReload(const Simulation& sim, const Entity& e) const;
    // One function per role: builds that role of entity `index` fresh from the source of the entity
    // and returns what it created.
    int addRigid(Simulation& sim, int index);                  // rigid body index
    int addSoft(Simulation& sim, int index);                   // particle group
    int addLiquid(Simulation& sim, int index);                 // particle group
    int addCloth(Simulation& sim, int index);                  // particle group
    void addMagnet(Simulation& sim, int index, int body);
    void addEmitter(Simulation& sim, int index, int body);
    void releaseFromEmitter(Simulation& sim, const EmitterRef& ref, bool& liquidDone) const;
    void addMagnetFieldToGas(Simulation& sim) const;

    SceneGraph graph_;
    std::vector<int> magnetBody_;       // rigid body index of each magnet
    std::vector<Vector3> magnetMoment_; // its moment in the body frame
    std::vector<uint32_t> magnetEntity_; // the id of its entity
    float maxMagnetForce_ = 0;          // the largest pair force of the last frame [N]
    std::vector<uint32_t> entityIdOfBody_; // body index -> entity id (0: none)
    std::vector<EmitterRef> emitters_;  // every visible entity with the emitter role
    std::unordered_map<uint32_t, std::vector<MetaObject>> meta_; // entity id -> its meta-objects
    std::vector<std::string> notes_;    // roles that could not be honoured, shown in the readings
};

} // namespace rf
