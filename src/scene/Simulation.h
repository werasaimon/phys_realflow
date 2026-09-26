#pragma once
// Scene container tying the solvers together: presets, obstacle geometry, time stepping and
// extraction of render/analysis data (RenderSnapshot). No GUI dependencies.

#include "spatial/BVH.h"
#include "core/Mesh.h"
#include "gas/GasSolver.h"
#include "gas/SurfaceLoads.h"
#include "plasma/Tokamak.h"
#include "rigid/RigidWorld.h"
#include "particles/ParticleSystem.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace rf {

enum class SimMode { Fluid, WindTunnel, Rigid };

enum class Preset {
    DamBreak, FluidObstacle, FloatingBodies, JetOnObject,
    TunnelSphere, TunnelCylinder, TunnelWing, TunnelStreamlined, TunnelCube, SmokePlume, SmokeSphere,
    RigidFalling, RigidGranular, RigidPyramid, RigidConvex, RigidTower, RigidJoints, RigidCcd, RigidTeapots, SmokeBodies, SoftCloth, GasSoftCloth, Hydro, Fire, Water, Magnetosphere, Tokamak,
    Count
};
const char* presetName(Preset p);
SimMode presetMode(Preset p);

// Non-convex models decomposed once into convex parts (cached, shared by all bodies).
std::shared_ptr<const CompoundShape> teapotShape();
std::shared_ptr<const CompoundShape> bunnyShape();

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

    // Current settings, for synchronising the UI
    uint64_t paramsVersion = 0;
    Preset preset = Preset::DamBreak;
    ParticleParams particleParams;
    GasParams gasParams;
    Tokamak tokamak;
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
    Tokamak tokamak; // the Tokamak scene: its vessel, coils and plasma current (rebuilt by reset())
    ObstacleSettings obstacle;
    VisSettings vis;
    float frameDt = 1.0f / 60.0f;
    // Gas <-> rigid / soft / cloth coupling (gas mode): bodies are always moving obstacles for the
    // gas; with this flag the gas pressure (and the buoyancy of the displaced gas) also acts back on
    // the bodies, and the aerodynamic drag on the cloth.
    bool gasPushesBodies = true;

    Simulation();

    SimMode mode() const { return mode_; }
    Preset preset() const { return preset_; }
    // One gravity for the whole scene: rigid bodies, particles (liquid, soft bodies, cloth) and the
    // buoyancy of the flame. Each solver keeps its own copy (it can run alone); set it here.
    void setGravity(const Vector3& g);
    Vector3 gravity() const { return rigid.params.gravity; }
    void loadPreset(Preset p);
    void reset();                    // rebuild the current preset with the current parameters
    void rebuildObstacle();          // regenerate geometry from `obstacle`, then reset()
    bool loadCustomMesh(const std::string& path, std::string& error);
    void stepFrame();

    float time() const { return time_; }
    uint64_t frame() const { return frame_; }
    void touchParams() { ++paramsVersion_; }
    void fillSnapshot(RenderSnapshot& s) const;
    const TriMesh& obstacleMesh() const { return *obstacleMesh_; }
    // Wind tunnel: loads on every triangle of the obstacle (updated each frame).
    const SurfaceLoads& surfaceLoads() const { return surfaceLoads_; }

private:
    void buildObstacleGeometry();
    void setupFluidScene();
    void setupTunnelScene();
    void setupRigidScene();
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

    SimMode mode_ = SimMode::Fluid;
    Preset preset_ = Preset::DamBreak;
    std::shared_ptr<TriMesh> obstacleMesh_;
    MeshBVH obstacleBVH_;
    uint64_t obstacleVersion_ = 1;
    uint64_t paramsVersion_ = 1;
    AABB fluidDomain_, rigidDomain_;
    float time_ = 0;
    uint64_t frame_ = 0;
    float lastStepMs_ = 0;
    float lastGridDt_ = 0;
    // Tokamak position control (the vertical field under feedback): its state.
    void controlTokamakPosition();
    float tokamakBv_ = 0, tokamakShift_ = 0, tokamakControlTime_ = 0;
    std::vector<Vector3> gasImpulse_, gasAngularImpulse_; // gas -> bodies, accumulated over a frame
    std::vector<Vector3> softGasImpulse_;                  // gas -> soft bodies (their centre of mass)
    float releaseTime_ = 0; // bodies held by the scene (RigidWorld::hold) fall from this time on
    float gasForceMax_ = 0;
    SurfaceLoads surfaceLoads_;
    void updateSurfaceLoads();
};

} // namespace rf
