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

#include <string>
#include <vector>

namespace rf {

enum class ShapeKind { Box, Sphere, Cylinder, Cone, Plane };

// The roles an entity can have. One entity may have several (a rigid body that is also a magnet);
// a few exclude each other (rigid, soft and liquid - the shape is one of them).
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

struct Entity {
    std::string name;
    ShapeKind shape = ShapeKind::Box;
    Vector3 size{0.2f};         // full size: box edges; sphere/cylinder/cone: diameter, height in y
    Vector3 position{0.0f};
    Vector3 rotationDeg{0.0f};  // Euler angles in degrees: yaw (y), pitch (x), roll (z)
    Vector3 color{0.8f, 0.8f, 0.8f};
    RigidRole rigid;
    SoftRole soft;
    LiquidRole liquid;
    MagnetRole magnet;
    HeatRole heat;
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
    void addEntity(Simulation& sim, const Entity& e);
    int addRigidBody(Simulation& sim, const Entity& e, bool fixed);
    void addMagnetFieldToGas(Simulation& sim) const;

    SceneGraph graph_;
    std::vector<int> magnetBody_;       // rigid body index of each magnet entity
    std::vector<Vector3> magnetMoment_; // its moment in the body frame
    float maxMagnetForce_ = 0;          // the largest pair force of the last frame [N]
    std::vector<int> entityOfBody_;     // body index -> entity index
};

} // namespace rf
