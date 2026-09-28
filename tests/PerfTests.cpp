// Speed, measured where the user felt it: a thousand cubes falling onto a floor, exactly as the
// editor builds them (a 10 x 10 x 10 array of 0.2 m cubes from a scene file, the Simulation facade,
// 60 frames a second, the rigid solver's own substeps). The test prints the milliseconds of every
// stage of a frame - collide, solve, integrate, continuous collision - for three phases: the fall
// (everything moving, thousands of contacts), the settling (a pile shaking down) and the rest
// (almost all asleep). Machines differ, so the test only guards against a regression several times
// over the measured value; the numbers themselves are for reading (docs/02-rigid-bodies.md).
#include "TestRunner.h"

#include "core/Parallel.h"
#include "core/Probe.h"
#include "scene/SceneGraph.h"
#include "scene/Simulation.h"

#include <algorithm>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using namespace rf;

namespace {

// The scene file of the editor's benchmark (C:/Users/PC/rf-bench/cubes-1000.rfscene): a floor, one
// template cube and an array of 10 x 10 x 10 of its copies, 0.3 m apart, lightly jittered.
const char* kThousandCubes =
    "world gravity 0 -9.81 0 size 8 12 8 gas 0 magneticGas 0\n"
    "entity \"Floor\"\n"
    "  object id 1 visible 1 locked 1\n"
    "  shape plane size 8 0.02 8 position 0 0.01 0 rotation 0 0 0 color 0.6 0.62 0.66\n"
    "  rigid density 500 friction 0.5 restitution 0.2 fixed 1 velocity 0 0 0 spin 0 0 0\n"
    "end\n"
    "entity \"Cube\"\n"
    "  object id 2 visible 1 locked 0\n"
    "  shape box size 0.2 0.2 0.2 position -1.4 0.3 -1.4 rotation 0 0 0 color 0.9 0.5 0.25\n"
    "  rigid density 500 friction 0.5 restitution 0.2 fixed 0 velocity 0 0 0 spin 0 0 0\n"
    "end\n"
    "array \"Thousand\"\n"
    "  object id 3 visible 1 locked 0\n"
    "  pose position 0 0 0 rotation 0 0 0 color 1 1 1\n"
    "  pattern grid template 2 count 10 10 10 step 0.3 0.3 0.3 radius 0 rotationStep 0 0 0 jitter 0.02 seed 7\n"
    "end\n";

// What one frame cost, stage by stage [ms], and how busy it was.
struct FrameCost {
    double step = 0, collide = 0, solve = 0, integrate = 0, ccd = 0, contacts = 0, awake = 0;
};

// The median of each stage over frames [first, last): a median, not a mean, so one frame the
// operating system took away does not move the number.
FrameCost medianOf(const std::vector<FrameCost>& frames, int first, int last) {
    auto median = [&](double FrameCost::*field) {
        std::vector<double> v;
        for (int i = first; i < last && i < int(frames.size()); ++i) v.push_back(frames[size_t(i)].*field);
        if (v.empty()) return 0.0;
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
    };
    FrameCost m;
    m.step = median(&FrameCost::step), m.collide = median(&FrameCost::collide), m.solve = median(&FrameCost::solve);
    m.integrate = median(&FrameCost::integrate), m.ccd = median(&FrameCost::ccd);
    m.contacts = median(&FrameCost::contacts), m.awake = median(&FrameCost::awake);
    return m;
}

void printPhase(const char* name, const FrameCost& c) {
    std::printf("  %-14s step %6.2f ms = collide %5.2f + solve %5.2f + integrate %5.2f (of it ccd %5.2f); "
                "%5.0f contacts, %4.0f awake\n",
                name, c.step, c.collide, c.solve, c.integrate, c.ccd, c.contacts, c.awake);
}

} // namespace

void testRigidPerfThousandCubes() {
    // 1. The scene, as the editor loads it.
    SceneGraph graph;
    std::string error;
    CHECK(graph.load(kThousandCubes, error), "the benchmark scene does not load: %s", error.c_str());
    Simulation sim;
    sim.load(std::make_unique<GraphScene>(graph));

    // 2. Three hundred frames (5 s), the probe read after each one.
    std::vector<FrameCost> frames;
    for (int f = 0; f < 300; ++f) {
        sim.stepFrame();
        const Probe::Snapshot s = Probe::snapshot();
        FrameCost c;
        c.step = s.value("frame/step ms");
        c.collide = s.value("rigid/collide ms"), c.solve = s.value("rigid/solve ms");
        c.integrate = s.value("rigid/integrate ms"), c.ccd = s.value("rigid/ccd ms");
        c.contacts = s.value("rigid/contacts"), c.awake = s.value("rigid/bodies awake");
        frames.push_back(c);
    }

    // 3. The three phases: the fall, the settling pile, the rest.
    const FrameCost fall = medianOf(frames, 0, 60), settle = medianOf(frames, 60, 180), rest = medianOf(frames, 180, 300);
    std::printf("  1000 cubes, %d threads, medians per 1/60 s frame:\n", ThreadPool::instance().threadCount());
    printPhase("fall 0-60", fall);
    printPhase("settle 60-180", settle);
    printPhase("rest 180-300", rest);

    // 4. A loose guard: three times what the 24-thread machine of docs/02-rigid-bodies.md ("Скорость")
    //    measured - 31.5 ms a frame for the fall, 9.7 ms for the rest (the goals were 16 and 2).
    const double measuredFall = 31.5, measuredRest = 9.7;
    CHECK(fall.step < 3 * measuredFall, "the fall costs %.1f ms a frame", fall.step);
    CHECK(rest.step < 3 * measuredRest, "the rest costs %.1f ms a frame", rest.step);
}
