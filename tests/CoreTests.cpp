// The core utilities: the Probe (channels, counters, timers, debug drawing) and the allocation
// census - how many heap allocations a frame of each scene makes, the number that decides whether
// the solvers need a per-step arena allocator (Box2D's b2StackAllocator, PhysX's scratch memory).
#include "TestRunner.h"

#include "core/Format.h"
#include "core/Probe.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <thread>

void testProbe() {
    Probe::clear();
    Probe::enableDraw(false);

    // Values keep the latest number, counters sum within a frame, timers add up milliseconds.
    Probe::set("test/value", 1.0);
    Probe::set("test/value", 2.5);
    Probe::add("test/counter", 3);
    Probe::add("test/counter", 4);
    {
        Probe::Timer t("test/timer ms");
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    {
        Probe::Timer t("test/timer ms"); // a second scope adds to the same channel
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    Probe::Snapshot s = Probe::snapshot();
    CHECK(s.value("test/value") == 2.5, "set keeps the latest value: %g", s.value("test/value"));
    CHECK(s.value("test/counter") == 7.0, "add sums within a frame: %g", s.value("test/counter"));
    CHECK(s.value("test/timer ms") >= 3.0, "timers add up: %g ms", s.value("test/timer ms"));
    CHECK(s.value("test/missing", -1.0) == -1.0, "a missing channel returns the fallback");
    // Sorted by name, and the kinds are remembered.
    CHECK(s.channels.size() == 3, "3 channels, got %zu", s.channels.size());
    for (size_t i = 1; i < s.channels.size(); ++i)
        CHECK(s.channels[i - 1].name < s.channels[i].name, "channels are sorted by name");
    for (const Probe::Channel& c : s.channels) {
        if (c.name == "test/value") CHECK(c.kind == Probe::Kind::Value, "kind of a value");
        if (c.name == "test/counter") CHECK(c.kind == Probe::Kind::Counter, "kind of a counter");
        if (c.name == "test/timer ms") CHECK(c.kind == Probe::Kind::TimerMs, "kind of a timer");
    }
    std::vector<std::string> names = Probe::channels();
    CHECK(names.size() == 3 && names[0] == "test/counter" && names[1] == "test/timer ms" && names[2] == "test/value",
          "channels() lists the sorted names");

    // A new frame: counters and timers back to zero, values stay.
    Probe::beginFrame();
    s = Probe::snapshot();
    CHECK(s.value("test/value") == 2.5, "beginFrame keeps values: %g", s.value("test/value"));
    CHECK(s.value("test/counter") == 0.0, "beginFrame zeroes counters: %g", s.value("test/counter"));
    CHECK(s.value("test/timer ms") == 0.0, "beginFrame zeroes timers: %g", s.value("test/timer ms"));
    CHECK(s.channels.size() == 3, "the channels survive the frame: %zu", s.channels.size());

    // Debug drawing: nothing while disabled, everything while enabled, cleared by the frame.
    Probe::line({0, 0, 0}, {1, 0, 0}, {1, 1, 1});
    Probe::point({0, 0, 0}, {1, 1, 1});
    s = Probe::snapshot();
    CHECK(s.lines.empty() && s.points.empty(), "drawing is off by default: %zu lines, %zu points", s.lines.size(), s.points.size());
    Probe::enableDraw(true);
    CHECK(Probe::drawEnabled(), "enableDraw(true)");
    Probe::line({0, 0, 0}, {1, 0, 0}, {1, 0, 0});
    Probe::arrow({0, 0, 0}, {0, 1, 0}, {0, 1, 0});
    Probe::box(AABB({0, 0, 0}, {1, 1, 1}), {0, 0, 1});
    Probe::point({0.5f, 0.5f, 0.5f}, {1, 1, 0}, 0.02f);
    s = Probe::snapshot();
    CHECK(s.lines.size() == 1 + 3 + 12, "a line, an arrow (3 segments) and a box (12 edges): %zu segments", s.lines.size());
    CHECK(s.points.size() == 1 && s.points[0].size == 0.02f, "one point of size 0.02");
    CHECK(length(s.lines[0].b - Vector3(1, 0, 0)) < 1e-6f && s.lines[0].color.x == 1.0f, "the line as given");
    // The arrow's shaft ends at the tip; its two head strokes start there.
    CHECK(length(s.lines[1].b - Vector3(0, 1, 0)) < 1e-6f, "arrow shaft ends at the tip");
    CHECK(length(s.lines[2].a - Vector3(0, 1, 0)) < 1e-6f && length(s.lines[3].a - Vector3(0, 1, 0)) < 1e-6f, "arrow head starts at the tip");
    // The box's 12 edges all have unit length.
    for (size_t i = 4; i < 16; ++i) CHECK(std::fabs(length(s.lines[i].b - s.lines[i].a) - 1.0f) < 1e-6f, "box edge %zu has unit length", i);
    Probe::beginFrame();
    s = Probe::snapshot();
    CHECK(s.lines.empty() && s.points.empty(), "beginFrame clears the drawing: %zu lines, %zu points", s.lines.size(), s.points.size());
    Probe::enableDraw(false);
    CHECK(!Probe::drawEnabled(), "enableDraw(false)");
    Probe::clear();
    CHECK(Probe::channels().empty(), "clear() removes every channel");
}

void testAllocationsPerFrame() {
    // A census, not a pass/fail: heap allocations per frame of six scenes, counted by the replaced
    // operator new in main.cpp (the SDK never touches the allocator). Thousands per frame would
    // call for a per-step arena; tens mean the vectors already keep their capacity.
    // The limits are twice what the scenes make after the rigid solver was given persistent
    // scratch (manifold points inside the manifold, per-worker narrow-phase and EPA scratch, the
    // contact cache pruned in place): the tower went from 29 943 allocations per frame to 61, the
    // teapots from 140 075 to 85, the terrain from 50 152 to 102. What remains is the probe's own
    // channel lookups (a std::string key per timer) and the particle / gas solvers.
    struct Case { const char* name; Preset preset; int warmup, measure; double limit; };
    const Case cases[] = {{"rigid tower (100 boxes)", Preset::RigidTower, 30, 30, 150},
                          {"100 teapots", Preset::RigidTeapots, 30, 30, 200},
                          {"water: wave in a pool", Preset::Water, 30, 30, 2500},
                          {"fire: burning curtain", Preset::Fire, 10, 10, 1100},
                          {"magnetosphere", Preset::Magnetosphere, 10, 10, 1300},
                          {"terrain: 150 bodies on a mesh", Preset::Terrain, 30, 30, 250}};
    std::printf("  %-32s %12s %12s %12s %12s\n", "scene", "allocs mean", "allocs max", "step ms", "limit");
    for (const Case& c : cases) {
        Simulation sim;
        loadSample(sim, c.preset);
        for (int f = 0; f < c.warmup; ++f) sim.stepFrame();
        double sum = 0, worst = 0, ms = 0;
        for (int f = 0; f < c.measure; ++f) {
            sim.stepFrame();
            const Probe::Snapshot s = Probe::snapshot();
            const double a = s.value("memory/allocations per frame", -1.0);
            CHECK(std::isfinite(a) && a >= 0, "%s: allocations per frame = %g", c.name, a);
            sum += a;
            worst = std::max(worst, a);
            ms += s.value("frame/step ms");
        }
        std::printf("  %-32s %12.0f %12.0f %12.2f %12.0f\n", c.name, sum / c.measure, worst, ms / c.measure, c.limit);
        CHECK(sum / c.measure < c.limit, "%s: %.0f allocations per frame, limit %.0f", c.name, sum / c.measure, c.limit);
    }
}

void testTimersCoverTheStep() {
    // The Laboratory's profiler bar adds the stage timers of a frame up and paints what they miss
    // grey («Вне замеров»). So every stage of the step has a timer, and no timer runs inside another
    // (the bar would count that time twice): the timers add up to the whole step. The rest - the
    // probe's own bookkeeping, the loops between the stages - is allowed 10 %.
    struct Case { const char* name; Preset preset; int warmup, measure; };
    const Case cases[] = {{"rigid tower (100 boxes)", Preset::RigidTower, 30, 30},
                          {"100 teapots", Preset::RigidTeapots, 30, 30},
                          {"joints", Preset::RigidJoints, 30, 30},
                          {"water: wave in a pool", Preset::Water, 20, 20},
                          {"cloth and soft bodies", Preset::SoftCloth, 20, 20},
                          {"cloth in the wind", Preset::GasSoftCloth, 10, 10},
                          {"smoke and bodies", Preset::SmokeBodies, 10, 10},
                          {"wind tunnel: sphere", Preset::TunnelSphere, 10, 10},
                          {"fire: burning curtain", Preset::Fire, 10, 10},
                          {"magnetosphere", Preset::Magnetosphere, 10, 10},
                          {"terrain: 150 bodies on a mesh", Preset::Terrain, 30, 30}};
    std::printf("  %-32s %10s %10s %9s   %s\n", "scene", "step ms", "timed ms", "covered", "the three longest stages");
    for (const Case& c : cases) {
        Simulation sim;
        loadSample(sim, c.preset);
        for (int f = 0; f < c.warmup; ++f) sim.stepFrame();
        double step = 0, timed = 0;
        std::map<std::string, double> stages;
        for (int f = 0; f < c.measure; ++f) {
            sim.stepFrame();
            const Probe::Snapshot s = Probe::snapshot();
            step += s.value("frame/step ms");
            for (const Probe::Channel& ch : s.channels)
                if (ch.kind == Probe::Kind::TimerMs) timed += ch.value, stages[ch.name] += ch.value;
        }
        std::vector<std::pair<double, std::string>> longest;
        for (const auto& [name, ms] : stages) longest.push_back({ms / c.measure, name});
        std::sort(longest.rbegin(), longest.rend());
        std::string top;
        for (size_t i = 0; i < longest.size() && i < 3; ++i) top += format("%s %.2f  ", longest[i].second.c_str(), longest[i].first);
        const double covered = timed / std::max(step, 1e-9);
        std::printf("  %-32s %10.2f %10.2f %8.1f %%   %s\n", c.name, step / c.measure, timed / c.measure, covered * 100, top.c_str());
        CHECK(covered > 0.90, "%s: the timers cover %.1f %% of the step (no grey gap over 10 %%)", c.name, covered * 100);
        CHECK(covered < 1.02, "%s: the timers add up to %.1f %% of the step - one runs inside another", c.name, covered * 100);
    }
}
