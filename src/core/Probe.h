#pragma once
// The probe: one place every part of the engine reports to, and the viewer, the plots, the tests
// and the CSV take from by name - as Bullet's btIDebugDraw and profiler, Box2D's b2Profile and
// b2Counters, PhysX's Visual Debugger. Global on purpose: a solver deep inside a substep, a scene
// and the application all reach it with one line, and the reader asks for any channel by name
// without anybody having wired it up in advance.
//
//   Probe::set("rigid/contacts", n);              a value: the latest number of the channel
//   Probe::add("bvh/queries", 1);                 a counter: summed over the frame, zero at beginFrame
//   { Probe::Timer t("gas/pressure"); ... }       a timer: milliseconds summed over the frame
//   if (Probe::layerOn(DrawLayer::ContactNormals))  debug drawing by LAYER (as PhysX's
//       Probe::arrow(DrawLayer::ContactNormals, p, n, colour);   visualization parameters)
//
// Rules: report per frame or per substep, not per particle (the map is behind a mutex); the
// names are "part/quantity" in ASCII so the channels sort into groups. Test the layer BEFORE a
// loop that draws: with every layer off the drawing costs one atomic load per solver per step.
// Every layer is capped at kLayerCap primitives per frame (a label "(обрезано)" then says so).
// Simulation::stepFrame calls beginFrame() and copies snapshot() into the RenderSnapshot.
#include "math/AABB.h"
#include "math/Vector3.h"

#include <atomic>
#include <chrono>
#include <cstdint>
#include <string>
#include <vector>

namespace rf {

// What the research view can show, one bit each (the editor's panel toggles them). The analogues:
// PhysX PxVisualizationParameter, Bullet btIDebugDraw::DebugDrawModes, Box2D b2Draw flags.
enum class DrawLayer : uint32_t {
    Misc,               // drawings without a layer (the old line/arrow/point/box calls)
    ContactPoints,      // rigid contact points
    ContactNormals,     // their normals (0.1 m)
    ContactImpulses,    // normal force per point: impulse / dt, 1 m of arrow = 100 N (kForceScale)
    PenetrationDepth,   // a segment of the penetration depth along the normal, red
    BodyAabbs,          // every awake body's world box
    WorldTree,          // the rigid world's dynamic AABB tree, nodes coloured by depth
    MeshBvh,            // the static mesh's BVH nodes down to kBvhDepth
    CentreOfMass,       // a small cross at each body's centre of mass
    InertiaAxes,        // principal axes, length proportional to sqrt(I / m)
    Velocities,         // linear (cyan) and angular (magenta) velocity arrows
    Islands,            // bodies' boxes, one colour per contact island
    Sleeping,           // sleeping bodies' boxes, dimmed
    JointFrames,        // joint anchors, the link between them and the joint axis
    GjkSimplex,         // the watched pair's GJK simplices (world witnesses + a Minkowski inset)
    EpaPolytope,        // the watched pair's final EPA polytope (in the inset)
    WitnessPoints,      // the watched pair's witness points and separation / penetration vector
    ParticleNeighbours, // the particle nearest to probePoint() and lines to its neighbours
    DensityError,       // particles coloured by rho / rho0 - 1 (blue low, red high)
    SoftTetrahedra,     // soft-body tetrahedra, edges coloured by the volume ratio det F
    ClothTension,       // cloth edges coloured by tension / strength (blue slack, red tearing)
    GasGrid,            // the grid cells of the slice plane
    GasVelocity,        // gas velocity arrows on the slice plane
    PressureGradient,   // -grad p / rho arrows on the slice plane (the acceleration p gives)
    Divergence,         // div u per cell on the slice, coloured (what the projection removes)
    Vorticity,          // curl u arrows on the slice
    FieldLinesB,        // magnetic field lines, coloured by |B|
    CurrentDensity,     // J = curl B / mu0 arrows on the slice
    Count
};

class Probe {
public:
    enum class Kind { Value, Counter, TimerMs };
    struct Channel {
        std::string name;
        double value = 0;
        Kind kind = Kind::Value;
    };
    struct Line { Vector3 a, b, color; DrawLayer layer = DrawLayer::Misc; };
    struct Point { Vector3 p, color; float size; DrawLayer layer = DrawLayer::Misc; };
    struct Label { Vector3 at; std::string text; DrawLayer layer = DrawLayer::Misc; };
    struct Snapshot {
        std::vector<Channel> channels; // sorted by name
        std::vector<Line> lines;
        std::vector<Point> points;
        std::vector<Label> labels;
        double value(const std::string& name, double fallback = 0) const;
    };

    static void set(const char* name, double value);
    static void add(const char* name, double delta);
    class Timer {
    public:
        explicit Timer(const char* name) : name_(name), t0_(std::chrono::steady_clock::now()) {}
        ~Timer();
        Timer(const Timer&) = delete;
        Timer& operator=(const Timer&) = delete;
    private:
        const char* name_;
        std::chrono::steady_clock::time_point t0_;
    };

    // --- Layers --------------------------------------------------------------------------------
    static constexpr uint32_t bit(DrawLayer l) { return 1u << uint32_t(l); }
    static void setLayers(uint32_t mask) { layers_.store(mask, std::memory_order_relaxed); }
    static uint32_t layers() { return layers_.load(std::memory_order_relaxed); }
    static bool layerOn(DrawLayer l) { return (layers() & bit(l)) != 0; }
    // The bits of the layers first..last (in the order of DrawLayer), e.g. every rigid layer.
    static constexpr uint32_t bits(DrawLayer first, DrawLayer last) {
        return (bit(last) - bit(first)) | bit(last); // first .. last-1, and last
    }
    // Removes this frame's drawing of the layers in mask. A solver stepped several times a frame
    // (substeps) calls it before it draws, so the picture is that of its last step.
    static void clearLayers(uint32_t mask);
    static const char* layerName(DrawLayer l); // "ContactPoints" ... (ASCII, for panels and files)
    static constexpr int kLayerCap = 20000;    // primitives per layer per frame
    static constexpr float kForceScale = 0.01f; // ContactImpulses: metres of arrow per newton
    static constexpr int kBvhDepth = 8;          // MeshBvh: deepest level drawn

    // The old switch, kept for the viewer's "Отладочная отрисовка" box: on = the basic layers.
    static bool drawEnabled() { return layers() != 0; }
    static void enableDraw(bool on);

    // Drawing by layer (nothing is stored when the layer is off). Colours 0..1.
    static void line(DrawLayer l, const Vector3& a, const Vector3& b, const Vector3& color);
    static void arrow(DrawLayer l, const Vector3& p, const Vector3& v, const Vector3& color); // a line with a head
    static void point(DrawLayer l, const Vector3& p, const Vector3& color, float size = 0.01f);
    static void box(DrawLayer l, const AABB& b, const Vector3& color); // its 12 edges
    static void label(DrawLayer l, const Vector3& at, const std::string& text);
    // Without a layer: DrawLayer::Misc.
    static void line(const Vector3& a, const Vector3& b, const Vector3& color) { line(DrawLayer::Misc, a, b, color); }
    static void arrow(const Vector3& p, const Vector3& v, const Vector3& color) { arrow(DrawLayer::Misc, p, v, color); }
    static void point(const Vector3& p, const Vector3& color, float size = 0.01f) { point(DrawLayer::Misc, p, color, size); }
    static void box(const AABB& b, const Vector3& color) { box(DrawLayer::Misc, b, color); }

    // --- What the research view points at ------------------------------------------------------
    // A body pair whose GJK / EPA is drawn (GjkSimplex, EpaPolytope, WitnessPoints); -1 -1: none.
    static void watchPair(int bodyA, int bodyB);
    static void watchedPair(int& bodyA, int& bodyB);
    // A point in the world the viewer's mouse is at (ParticleNeighbours picks the particle nearest to it).
    static void setProbePoint(const Vector3& p);
    static Vector3 probePoint();

    // The frame: counters and timers back to zero, drawings cleared; values stay.
    static void beginFrame();
    // How many times beginFrame() and clear() ran: whoever draws once a frame from outside the
    // step (the viewer's slice) remembers it and draws again only when it changed.
    static uint64_t frameIndex() { return frames_.load(std::memory_order_relaxed); }
    static Snapshot snapshot();
    static std::vector<std::string> channels();
    static void clear(); // everything, including the values (tests)
    static size_t primitiveCount(); // lines + points + labels stored this frame (tests)

    // Memory: whoever replaces the global operator new (the tests, the app) adds to this; the
    // engine itself never touches the allocator. stepFrame reports the difference per frame.
    static std::atomic<long long> allocations;

private:
    static std::atomic<uint32_t> layers_;
    static std::atomic<uint64_t> frames_;
};

// A colour from blue (t = 0) through green to red (t = 1), for the coloured layers.
Vector3 heatColor(float t);

} // namespace rf
