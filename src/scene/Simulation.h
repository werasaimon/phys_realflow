#pragma once
// Scene container tying the solvers together: the boxes they share, the obstacle geometry, time
// stepping, the coupling between them and the extraction of render/analysis data
// (RenderSnapshot). What is in the scene comes from a Scene (Scene.h); the ready-made scenes are
// in samples/. No GUI dependencies.

#include "spatial/BVH.h"
#include "core/Mesh.h"
#include "core/Probe.h"
#include "gas/GasSolver.h"
#include "gas/SurfaceLoads.h"
#include "rigid/RigidWorld.h"
#include "particles/ParticleSystem.h"
#include "scene/Scene.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace rf {

// How a frame is stepped: liquid/soft bodies with rigid bodies (Fluid, Rigid) or the gas with
// everything in it (WindTunnel). Set by the box the scene builds in (useLiquidTank, useGasBox,
// useRigidArena).
enum class SimMode { Fluid, WindTunnel, Rigid };

enum class ObstacleShape { None, Sphere, Cube, Cylinder, Wing, Streamlined, Ellipsoid, Cone, Custom };

struct ObstacleSettings {
    ObstacleShape shape = ObstacleShape::None;
    float size = 0.5f;        // characteristic length (diameter / chord / edge) [m]
    float span = 1.2f;        // wing span / cylinder length [m]
    float thickness = 0.3f;   // thickness ratio for the streamlined body
    std::string nacaCode = "2412";
    Vector3 position{0.0f};
    float yawDeg = 0, angleOfAttackDeg = 0, rollDeg = 0;
    std::shared_ptr<TriMesh> customMesh;
    std::string customName;
};

enum class ParticleColoring { Speed, Density, Uniform };

struct VisSettings {
    GridField sliceField = GridField::Speed;
    int sliceAxis = 2;           // 0 = X, 1 = Y, 2 = Z (plane normal)
    float slicePosition = 0.5f;  // 0..1
    bool showSlice = true;
    bool showSmoke = true;
    bool liquidSurface = false;  // liquid drawn as a water surface (screen-space) instead of spheres
    bool showFieldLines = true;  // magnetic field lines (plasma scenes)
    bool planetSurface = false;  // draw the (spherical) obstacle as an Earth-like planet
    bool vesselGlass = false;    // draw the obstacle as a glass vessel (translucent, after the volume)
    bool showStreamlines = true;
    int streamlineSeeds = 12;    // per side
    bool surfacePressure = true; // color the obstacle by Cp
    ParticleColoring particleColoring = ParticleColoring::Speed;
    bool autoRange = true;
    float rangeMin = 0, rangeMax = 1;

    // Voxel grid / velocity vectors (grid solver)
    int gridDisplay = 0;          // 0 hidden, 1 in the slice plane, 2 whole 3D lattice
    int vectorDisplay = 0;        // 0 hidden, 1 in the slice plane, 2 whole volume
    int vectorStride = 1;         // every n-th cell
    float vectorScale = 1.0f;     // arrow length multiplier (1 = longest arrow spans stride cells)
    bool vectorsWhereSmoke = false;
};

struct RenderSnapshot {
    uint64_t frame = 0;
    SimMode mode = SimMode::Fluid;
    AABB domain;

    // SPH particles
    std::vector<Vector3> particles;
    std::vector<float> particleScalar;
    float particleRadius = 0.01f;
    // Colour of particles that are not coloured by the scalar (soft bodies): used where
    // particleScalar is kOwnColor.
    static constexpr float kOwnColor = -1e9f;
    std::vector<Vector3> particleColor;
    // Liquid particles drawn as a water surface (vis.liquidSurface) - then not in `particles`.
    std::vector<Vector3> liquid;
    float liquidRadius = 0.01f;
    // Cloth sheets (particle grids), drawn as surfaces.
    struct ClothMesh {
        int width = 0, height = 0;
        std::vector<Vector3> positions; // x + width * y
        std::vector<uint8_t> cellIntact; // (w-1)*(h-1); cut cells are not drawn (see Cloth)
        std::vector<float> burnt;        // per particle 0 fresh .. 1 charred (empty: does not burn)
        Vector3 color;
    };
    std::vector<ClothMesh> cloths;
    // Soft bodies drawn as their skinned surfaces.
    struct SoftMesh {
        std::vector<Vector3> positions;
        std::vector<std::array<uint32_t, 3>> triangles;
        Vector3 color;
    };
    std::vector<SoftMesh> softMeshes;

    // Rigid bodies
    struct Body {
        ShapeType shape;
        Vector3 pos;
        Quaternion rot;
        Vector3 halfExtents;
        float radius;
        Vector3 color;
        std::shared_ptr<const TriMesh> mesh; // convex polyhedra (local, principal frame)
        bool sleeping = false;
        std::shared_ptr<const ConvexShape> collisionShape; // for ray picking in the GUI (immutable)
        bool movable = true; // not static (nor held by the scene): the mouse can grab it
    };
    std::vector<Body> bodies;
    // Joints
    struct JointVis {
        JointType type;
        Vector3 anchorA, anchorB, axis;
    };
    std::vector<JointVis> joints;
    // Mouse joint
    bool grabActive = false;
    Vector3 grabAnchor, grabTarget;

    // Static obstacle
    uint64_t obstacleVersion = 0;
    std::shared_ptr<const TriMesh> obstacle;
    std::vector<float> obstacleScalar; // per vertex Cp (wind tunnel)

    // Grid slice
    bool hasSlice = false;
    int sliceW = 0, sliceH = 0;
    std::vector<float> slice;
    std::vector<uint8_t> sliceSolid;
    Vector3 sliceCorner, sliceU, sliceV;

    // Smoke volume
    bool hasVolume = false;
    int volX = 0, volY = 0, volZ = 0;
    std::vector<uint8_t> volume;
    // Fire: 2 channels per cell (smoke, temperature) instead of 1; temperature byte 255 =
    // volumeTemperatureScale kelvin above ambientTemperature.
    int volumeChannels = 1;
    float volumeTemperatureScale = 2000.0f, ambientTemperature = 293.0f;
    bool volumePlasma = false; // the tracer is glowing plasma (emits light) rather than smoke

    // Magnetic field lines (traced along B) and |B| along them [T]
    std::vector<std::vector<Vector3>> fieldLines;
    std::vector<std::vector<float>> fieldLineStrength;
    float fieldLineMax = 1;

    // Voxel lattice of the grid solver
    int gridNx = 0, gridNy = 0, gridNz = 0, sliceLayer = 0;
    float gridDx = 0;
    Vector3 gridOrigin;

    // Velocity vectors at cell centres
    std::vector<Vector3> arrowPos, arrowVel;
    float arrowMax = 1;    // max |v| among arrows [m/s]
    float arrowLength = 0; // world length of the longest arrow [m]

    // Streamlines
    std::vector<std::vector<Vector3>> streamlines;
    std::vector<std::vector<float>> streamlineSpeed;
    float streamlineMax = 1;

    // Colour mapping for the main scalar (slice / particles / surface)
    float colorMin = 0, colorMax = 1;
    std::string colorLabel;

    // Statistics
    float time = 0;
    float stepMs = 0;
    std::vector<std::pair<std::string, std::string>> info;
    std::vector<std::pair<std::string, float>> plots;
    // Everything the engine reported to the Probe this frame: every channel by name (values,
    // counters, timers) and the debug drawing. The plots take any channel from here.
    Probe::Snapshot probe;

    // Current settings, for synchronising the UI
    uint64_t paramsVersion = 0;
    std::string sceneName;              // of the loaded Scene ("" - none)
    std::vector<SceneParam> sceneParams; // its knobs
    ParticleParams particleParams;
    GasParams gasParams;
    RigidParams rigid;
    ObstacleSettings obstacleSettings;
    VisSettings vis;
    ParticleEmitter emitter;
    HeatSource heat;
    bool gasPushesBodies = true;
};

class Simulation {
public:
    ParticleSystem particles;
    GasSolver grid;
    RigidWorld rigid;
    ObstacleSettings obstacle;
    VisSettings vis;
    float frameDt = 1.0f / 60.0f;
    // Gas <-> rigid / soft / cloth coupling (gas mode): bodies are always moving obstacles for the
    // gas; with this flag the gas pressure (and the buoyancy of the displaced gas) also acts back on
    // the bodies, and the aerodynamic drag on the cloth.
    bool gasPushesBodies = true;
    // Bodies held by the scene (RigidWorld::hold) are released at this time.
    float releaseTime = 0;

    Simulation();

    // --- The scene -------------------------------------------------------------------------
    // Takes the scene, puts every parameter back to its default, lets the scene configure and
    // build. Without a scene the simulation is an empty tank.
    void load(std::unique_ptr<Scene> scene);
    Scene* scene() { return scene_.get(); }
    const Scene* scene() const { return scene_.get(); }
    void reset();                    // rebuild the scene with the current parameters
    void rebuildObstacle();          // regenerate geometry from `obstacle`, then reset()
    bool loadCustomMesh(const std::string& path, std::string& error);
    std::vector<SceneParam> sceneParams() const { return scene_ ? scene_->params() : std::vector<SceneParam>(); }
    void setSceneParam(int index, float value); // then reset()

    // --- The boxes a scene builds in (called first in Scene::build) ---------------------------
    // A tank of liquid (and whatever bodies, cloth and soft bodies go into it): the particle
    // solver and the rigid bodies share `box`; the obstacle is their static geometry.
    void useLiquidTank(const AABB& box);
    // The gas box of grid.params.domainSize; its origin is -originFraction * size: {0.5, 0.5, 0.5}
    // centres it, {0.5, 0, 0.5} puts the floor at y = 0, {0.3, 0.5, 0.5} puts the obstacle (at the
    // origin) at 30 % of the length. Bodies, particles and cloth live in the same box. With
    // grid.vessel set the obstacle mesh is only drawn (as the vessel), not voxelised.
    void useGasBox(const Vector3& originFraction);
    // Rigid bodies (and particles, if any) in `box` with the obstacle as static geometry.
    void useRigidArena(const AABB& box);

    // --- Running --------------------------------------------------------------------------
    SimMode mode() const { return mode_; }
    // One gravity for the whole scene: rigid bodies, particles (liquid, soft bodies, cloth) and the
    // buoyancy of the flame. Each solver keeps its own copy (it can run alone); set it here.
    void setGravity(const Vector3& g);
    Vector3 gravity() const { return rigid.params.gravity; }
    void stepFrame();

    float time() const { return time_; }
    uint64_t frame() const { return frame_; }
    void touchParams() { ++paramsVersion_; }
    void fillSnapshot(RenderSnapshot& s) const;
    const TriMesh& obstacleMesh() const { return *obstacleMesh_; }
    const MeshBVH& obstacleBVH() const { return obstacleBVH_; }
    // Wind tunnel: loads on every triangle of the obstacle (updated each frame).
    const SurfaceLoads& surfaceLoads() const { return surfaceLoads_; }

    static const AABB kDefaultTank;  // the liquid tank of the samples: 2 x 1.2 x 0.8 m
    static const AABB kDefaultArena; // the rigid arena: 4 x 5 x 4 m

private:
    void buildObstacleGeometry();
    void stepGasWithBodies();
    // Rigid bodies and particles over one frame, interleaved (floating and every contact between
    // them needs it): rigid.params.substeps rigid steps and particles.params.substeps particle
    // steps, in time order. gasDrag: the gas acts on cloth and liquid before each particle step.
    void stepBodiesAndParticles(bool gasDrag);
    void applyGasOnSoftBodies();
    void applyGasDragOnCloth(float dt);
    void applyGasDragOnLiquid(float dt);
    void passLiquidToGas();
    void fillParticleSolids(RenderSnapshot& s) const;
    std::vector<MovingSolid> movingSolids() const;
    void computeStreamlines(RenderSnapshot& s) const;
    void computeFieldLines(RenderSnapshot& s) const;
    void extractSlice(RenderSnapshot& s) const;
    void extractVectors(RenderSnapshot& s) const;
    int sliceLayer() const;
    void updateSurfaceLoads();

    std::unique_ptr<Scene> scene_;
    SimMode mode_ = SimMode::Fluid;
    std::shared_ptr<TriMesh> obstacleMesh_;
    MeshBVH obstacleBVH_;
    uint64_t obstacleVersion_ = 1;
    uint64_t paramsVersion_ = 1;
    float time_ = 0;
    uint64_t frame_ = 0;
    float lastStepMs_ = 0;
    float lastGridDt_ = 0;
    std::vector<Vector3> gasImpulse_, gasAngularImpulse_; // gas -> bodies, accumulated over a frame
    std::vector<Vector3> softGasImpulse_;                  // gas -> soft bodies (their centre of mass)
    float gasForceMax_ = 0;
    SurfaceLoads surfaceLoads_;
};

} // namespace rf
