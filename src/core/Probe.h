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
//   Probe::line(a, b, colour);                    debug drawing, only while drawEnabled()
//
// Rules: report per frame or per substep, not per particle (the map is behind a mutex); the
// names are "part/quantity" in ASCII so the channels sort into groups. Simulation::stepFrame
// calls beginFrame() and copies snapshot() into the RenderSnapshot.
#include "math/AABB.h"
#include "math/Vector3.h"

#include <atomic>
#include <chrono>
#include <string>
#include <vector>

namespace rf {

class Probe {
public:
    enum class Kind { Value, Counter, TimerMs };
    struct Channel {
        std::string name;
        double value = 0;
        Kind kind = Kind::Value;
    };
    struct Line { Vector3 a, b, color; };
    struct Point { Vector3 p, color; float size; };
    struct Snapshot {
        std::vector<Channel> channels; // sorted by name
        std::vector<Line> lines;
        std::vector<Point> points;
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

    // Debug drawing. Cheap when off: one flag test. Colours 0..1.
    static bool drawEnabled() { return drawEnabled_.load(std::memory_order_relaxed); }
    static void enableDraw(bool on) { drawEnabled_.store(on, std::memory_order_relaxed); }
    static void line(const Vector3& a, const Vector3& b, const Vector3& color);
    static void arrow(const Vector3& p, const Vector3& v, const Vector3& color); // a line with a head
    static void point(const Vector3& p, const Vector3& color, float size = 0.01f);
    static void box(const AABB& b, const Vector3& color);                       // its 12 edges

    // The frame: counters and timers back to zero, drawings cleared; values stay.
    static void beginFrame();
    static Snapshot snapshot();
    static std::vector<std::string> channels();
    static void clear(); // everything, including the values (tests)

    // Memory: whoever replaces the global operator new (the tests, the app) adds to this; the
    // engine itself never touches the allocator. stepFrame reports the difference per frame.
    static std::atomic<long long> allocations;

private:
    static std::atomic<bool> drawEnabled_;
};

} // namespace rf
