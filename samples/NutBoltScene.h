#pragma once
// Coarse square-thread contact fixture: a fixed bolt and a freely moving driven nut.
// Axial travel is produced by geometry, without a screw joint or overwritten poses.
#include "samples/Samples.h"

namespace rf {

class NutBoltScene final : public Scene {
public:
    static constexpr float pitch = 0.012f;
    static constexpr int segments = 24;
    static constexpr int boltBody = 0, nutBody = 1;
    void configure(Simulation& sim) override;
    void build(Simulation& sim) override;
    void afterStep(Simulation& sim) override;
    void describe(const Simulation& sim, RenderSnapshot& out) const override;
    static std::shared_ptr<const CompoundShape> boltShape();
    static std::shared_ptr<const CompoundShape> nutShape();
    static Vector3 modelOrigin(const RigidBody& body);
    static float modelAngle(const RigidBody& body);
    double driveWork() const { return driveWork_; }
    float commandedSpeed() const { return commandedSpeed_; }

private:
    void drive(Simulation& sim);
    double driveWork_ = 0;
    float commandedSpeed_ = 0;
};

void addNutBoltSample(std::vector<SampleEntry>& out);

} // namespace rf
