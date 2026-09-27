#pragma once
// The tokamak scene: a solution built on the engine, not a part of it. The MHD solver
// (plasma/MagneticField) knows nothing about tokamaks; this scene gives it the coils' potentials,
// the plasma current, the vessel and the resistive vacuum of a Tokamak (Tokamak.h), runs the
// radial position control of a real machine after every frame and reports the readings a
// physicist looks at (current, safety factor, kink growth, Shafranov shift).
#include "samples/plasma/Tokamak.h"
#include "scene/Scene.h"

namespace rf {

class TokamakScene : public Scene {
public:
    Tokamak tokamak; // the device: R0 = 0.4 m, a = 0.12 m, b = 0.24 m, B0 = 2.5 mT, q_a = 0.7

    void configure(Simulation& sim) override;
    void build(Simulation& sim) override;
    void afterStep(Simulation& sim) override;
    void describe(const Simulation& sim, RenderSnapshot& s) const override;
    void fieldLineSeeds(const Simulation& sim, std::vector<Vector3>& seeds) const override;
    std::vector<SceneParam> params() const override;
    void setParam(int index, float value) override;

private:
    // Position control (the vertical field under feedback): its state.
    float bv_ = 0, shift_ = 0;
    double controlTime_ = 0; // the gas clock at the last control step [s]
    float shiftIntegral_ = 0; // time integral of the shift [m s] (the PID's I part)
};

} // namespace rf
