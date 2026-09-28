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
// A soft body's material (particles/SoftBody.h, SoftMaterial). The defaults are the editor's
// «Желе» preset. rfscene: `soft density D young E poisson NU friction MU`; an old file's
// `stiffness K` (the shape-matching fraction) still loads and sets E = youngFromStiffness(K);
// `model shape-matching` selects the legacy clusters (for comparison).
struct SoftRole {
    bool enabled = false;
    float density = 1050;        // kg/m^3
    float youngModulus = 1.5e4f; // E [Pa]
    float poissonRatio = 0.45f;  // nu
    float friction = 0.5f;       // Coulomb coefficient
    bool shapeMatching = false;  // the legacy model
    float stiffness = 0.3f;      // shape matching only: fraction back to the rest shape per substep, 0..1
};

// The soft materials the editor offers at one click (handbook values, rounded; Custom = the
// numbers typed in). Jelly: gelatin dessert, E 15 kPa, nu 0.45, 1050 kg/m^3, friction 0.5.
// Rubber: soft natural rubber, E 1 MPa, nu 0.47, 1100 kg/m^3, friction 0.8. Soft plastic: a
// flexible plastic, E 5 MPa, nu 0.4, 950 kg/m^3, friction 0.4 (stiffer plastics are beyond the
// small steps: they would behave like this one, docs/03-particles.md).
enum class SoftPreset { Jelly, Rubber, SoftPlastic, Custom };
SoftRole softPreset(SoftPreset preset);           // the preset's material (the role enabled; Custom: the defaults)
SoftPreset softPresetOf(const SoftRole& role);    // which preset the role's numbers are, Custom if none
const char* softPresetName(SoftPreset preset);    // «Желе», «Резина», «Мягкий пластик», «Своё»
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
    float bendCompliance = 0.1f; // m/N, ClothMaterial's default (cotton); 1e-3 stiff as card, larger drapes more
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
    // Hierarchy: the id of the parent object (0: none). position / rotationDeg are then RELATIVE to
    // the parent, and moving the parent moves the children (a group). worldPose() composes them.
    uint32_t parent = 0;
};

// A group: an object with no geometry of its own that other objects hang under (Ctrl+G in the
// editor). For physics its members stay separate bodies unless `glued`: then they are welded into
// one rigid body (the group's members' colliders become one compound).
struct Group : SceneObject {
    bool glued = false;
};

// A light of the scene, the three kinds every 3D package has. It is for the viewer only (the
// physics does not see it); the direction of Sun and Spot is the object's -y axis turned by its
// rotation (a light "looks down" by default).
enum class LightKind { Sun, Point, Spot };
struct Light : SceneObject {
    LightKind kind = LightKind::Point;
    float intensity = 1.0f;     // relative brightness (1 = the default light of the viewer)
    float range = 5.0f;         // Point / Spot: the distance where the light has faded out, m
    float coneDeg = 40.0f;      // Spot: full cone angle
    float softnessDeg = 10.0f;  // Spot: the soft edge of the cone
    bool shadows = false;       // the viewer may cast shadows from it (Sun first)
};

// A camera of the scene: a viewpoint the viewer can look through, take screenshots and recordings
// from. Looks along the object's -z axis turned by its rotation.
struct Camera : SceneObject {
    float fovDeg = 50.0f;       // vertical field of view
    float nearClip = 0.05f, farClip = 200.0f; // m
    bool active = false;        // the one the viewer looks through when "through the camera" is on
};

// An array: one object standing for many copies of a template entity laid out in a pattern, as the
// Array modifier of Blender or copy-to-points of Houdini. The simulation expands it into one set of
// meta-objects per copy; the editor changes the count or the spacing with one number. The copies
// are placed relative to the array's own pose; the template itself is not simulated.
enum class ArrayPattern { Line, Grid, Circle };
struct ArrayObject : SceneObject {
    uint32_t templateId = 0;          // the entity that is copied (usually hidden, as a template)
    ArrayPattern pattern = ArrayPattern::Line;
    int count[3] = {5, 1, 1};         // Line: count[0]; Grid: nx, ny, nz; Circle: count[0]
    Vector3 step{0.3f, 0.3f, 0.3f};   // Line: step between copies (a vector); Grid: spacing per axis
    float radius = 1.0f;              // Circle: the radius (copies around the array's y axis)
    Vector3 rotationStepDeg{0.0f};    // each next copy turned by this much more (a spiral staircase)
    float jitter = 0;                 // random offset of each copy, m (a fixed seed: the same every run)
    uint32_t seed = 1;
};

// A shape in the scene with its roles: what it is made of and what else it does. With no role at
// all it is geometry only: drawn, but no body, no collision, no part in the simulation.
struct Entity : SceneObject {
    // An instance shares the geometry and the roles of its master (the entity with this id): change
    // the density, the collider or the shape of one and all change (as 3ds Max's Instance). Only the
    // SceneObject part (name, pose, visibility...) stays its own. 0: an ordinary, independent entity.
    // resolveInstance() gives the effective entity.
    uint32_t instanceOf = 0;
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
    std::vector<Group> groups;
    std::vector<ArrayObject> arrays;
    std::vector<Light> lights;
    std::vector<Camera> cameras;
    // Not saved: the folder relative mesh files are found in (the scene file's own folder). load()
    // keeps it as it was.
    std::string baseDirectory;

    // Text format, one line per fact, readable and diff-able (see SceneGraph.cpp for a sample).
    std::string save() const;
    bool load(const std::string& text, std::string& error);
};

// ---------------------------------------------------------------------------
// Many objects at once (SceneHierarchy.cpp): hierarchy, instances, arrays. None of these call a
// solver; GraphScene and an editor use them alike.
// ---------------------------------------------------------------------------
// The object with this id - an entity, a group or an array (null: none).
const SceneObject* findObject(const SceneGraph& g, uint32_t id);
// Any object's orientation from its rotationDeg (Ry Rx Rz), and back: Euler degrees of a rotation.
Quaternion objectRotation(const SceneObject& o);
Vector3 eulerDegrees(const Quaternion& q);
// Where an object is in the world: its pose composed with its parents'. A missing parent or a
// cycle in the parent links ends the chain there (that object counts as a root).
void worldPose(const SceneGraph& g, const SceneObject& o, Vector3& position, Quaternion& rotation);
// The inverse: sets o's relative pose so that its world pose is (position, rotation).
void setWorldPose(const SceneGraph& g, SceneObject& o, const Vector3& position, const Quaternion& rotation);
// Visible, and so are all its parents.
bool effectivelyVisible(const SceneGraph& g, const SceneObject& o);
// The nearest glued group above the object (0: none): its rigid members are one body.
uint32_t gluedGroupOf(const SceneGraph& g, const SceneObject& o);
// An instance as it acts: the geometry and roles of its root master (instanceOf followed to the end)
// with its own SceneObject part (name, pose, visibility, id). Not an instance, or a missing master:
// the entity itself.
Entity resolveInstance(const SceneGraph& g, const Entity& e);
// The copies an array stands for, posed in the WORLD (parent 0), each with the template's geometry
// and roles, visible as the array is. Copy n of the array with id A has the id arrayCopyId(A, n):
// stable across rebuilds as long as the array keeps its id (ids of arrays up to 32767, 65536 copies).
std::vector<Entity> expandArray(const SceneGraph& g, const ArrayObject& a);
uint32_t arrayCopyId(uint32_t arrayId, int copy);
bool isArrayCopyId(uint32_t id);
uint32_t arrayOfCopyId(uint32_t id);
// Every entity the simulation builds, in the world: the graph's entities (instances resolved,
// parents composed, visibility inherited) followed by the copies of every array.
std::vector<Entity> worldEntities(const SceneGraph& g);
// A fresh id: the largest over entities, groups and arrays, plus one.
uint32_t nextId(const SceneGraph& g);

// Lights and cameras in the world (their parents included). A light shines along its -y axis and a
// camera looks along its -z axis, both turned by the object's world rotation; up is the camera's +y.
Vector3 lightDirection(const SceneGraph& g, const Light& l);
void cameraFrame(const SceneGraph& g, const Camera& c, Vector3& eye, Vector3& forward, Vector3& up);

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
    // The object to select for a body: its entity, the array a copy belongs to, or the glued group
    // whose members it is made of (0: none).
    uint32_t objectIdOfBody(int body) const;
    // The entities as the simulation built them: world poses, instances resolved, array copies.
    const std::vector<Entity>& simulatedEntities() const { return flat_; }

    // One entity between frames, without reloading the scene. rebuildEntity replaces the entity of
    // updated.id: its meta-objects go, new ones are built from the source where the entity is NOW
    // (its body's or particles' live pose, unless `updated` moves it) and carry on its velocity.
    // Returns false when a role of it needs something only a reload can switch on (the gas).
    bool rebuildEntity(Simulation& sim, const Entity& updated);
    // An instance master rebuilds its instances too (they share its roles), and arrays copying it.
    // Inside a glued group the whole group's body is rebuilt. The graph keeps poses relative to the
    // parents: a live world pose is written back relative.
    bool addEntity(Simulation& sim, const Entity& entity); // an id of 0 gets a new one
    void removeEntity(Simulation& sim, uint32_t id);
    // An array changed (count, pattern, spacing...) or was added: its old copies go; the new ones are
    // built, and a copy that existed before keeps where it is now and how it moves.
    bool rebuildArray(Simulation& sim, const ArrayObject& array);
    const std::vector<MetaObject>& metaObjects(uint32_t id) const;
    // Lights and cameras between frames: the physics never sees them, the next snapshot shows them.
    void setLightsAndCameras(const std::vector<Light>& lights, const std::vector<Camera>& cameras);

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

    int indexOf(uint32_t id) const;         // in flat_
    int authoredIndexOf(uint32_t id) const; // in graph_.entities
    bool shapeUsable(const Entity& e);
    void buildEntity(Simulation& sim, int index);
    // Glued groups: every rigid member of the group becomes a part of one compound body.
    bool gluedMember(int index) const;
    void buildGluedGroup(Simulation& sim, uint32_t group);
    void rebuildGluedGroup(Simulation& sim, uint32_t group);
    // The graph entity of flat_[index] written with the flat (world) pose, relative to its parent.
    void writeBackPose(int index);
    bool rebuildOne(Simulation& sim, uint32_t id, bool keepLivePose);
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

    SceneGraph graph_;                  // as authored: poses relative to parents, instances, arrays
    std::vector<Entity> flat_;          // as simulated: worldEntities(graph_), kept in step with it
    std::vector<uint32_t> flatGlue_;    // flat_[i]'s glued group (0: none)
    std::unordered_map<uint32_t, int> glueBody_; // glued group id -> its one body
    std::vector<int> magnetBody_;       // rigid body index of each magnet
    std::vector<Vector3> magnetMoment_; // its moment in the body frame
    std::vector<uint32_t> magnetEntity_; // the id of its entity
    float maxMagnetForce_ = 0;          // the largest pair force of the last frame [N]
    std::vector<uint32_t> entityIdOfBody_; // body index -> entity id (0: none)
    std::vector<EmitterRef> emitters_;  // every visible entity with the emitter role
    std::unordered_map<uint32_t, std::vector<MetaObject>> meta_; // entity id -> its meta-objects
    std::vector<std::string> notes_;    // roles that could not be honoured, shown in the readings
    AABB box_;                          // the box the scene lives in: the world's, grown to hold every entity
};

// The visible lights of a graph, posed in the world, into s.lights (GraphScene::describe, and an
// editor drawing the graph without a simulation). None: the list stays empty (the viewer's default).
void describeLights(const SceneGraph& g, RenderSnapshot& s);

} // namespace rf
