#pragma once
// A scene is a solution built on the engine, not a part of it: it chooses the parameters, puts
// bodies, particles and fields into the solvers of a Simulation, and may run a controller and
// report readings. The engine ships none; the ready-made ones live in samples/ (as Box2D keeps its
// demos in samples/), and an application writes its own the same way.
//
// Life of a scene:
//   Simulation::load(scene)  -> all parameters back to defaults -> configure() -> build()
//   Simulation::reset()      -> build() again with the parameters as they are now
//   Simulation::stepFrame()  -> the solvers -> afterStep()
//   Simulation::fillSnapshot -> the generic readings -> describe()
#include "math/Math.h"

#include <string>
#include <vector>

namespace rf {

class Simulation;
struct RenderSnapshot;

// A knob of the scene for the viewer's panel (a tokamak's safety factor, a wind speed ...). The
// viewer shows every scene's knobs the same way; setting one is followed by a reset().
struct SceneParam {
    std::string name;
    float value = 0, min = 0, max = 1, step = 0.1f;
    int decimals = 2;
    bool toggle = false; // on/off: value 0 or 1
    std::string tip;     // what the knob does, for the panel's tooltip
};

class Scene {
public:
    virtual ~Scene() = default;

    // Once, when the scene is loaded and every parameter is at its default: the scene's own
    // parameters, the obstacle and the view. The user may change them afterwards; a reset() keeps
    // the changes and only rebuilds the contents.
    virtual void configure(Simulation&) {}
    // At every reset: the contents. Starts with one of the boxes of Simulation (useLiquidTank,
    // useGasBox, useRigidArena), then adds bodies, particles, cloth, fields.
    virtual void build(Simulation&) = 0;
    // Once per frame after the solvers: controllers, timed events.
    virtual void afterStep(Simulation&) {}
    // Readings for the viewer: info lines and plotted values, after the generic ones.
    virtual void describe(const Simulation&, RenderSnapshot&) const {}
    // Where the viewer starts its magnetic field lines (plasma scenes); none: around the obstacle.
    virtual void fieldLineSeeds(const Simulation&, std::vector<Vector3>&) const {}
    // The scene's knobs, and the change of one of them (Simulation resets afterwards).
    virtual std::vector<SceneParam> params() const { return {}; }
    virtual void setParam(int, float) {}

    std::string name; // shown by the viewer; the sample registry sets it
};

} // namespace rf
