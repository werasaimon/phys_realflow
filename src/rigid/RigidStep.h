#pragma once
// Atomic rigid-world stepping: retry a rejected CCD trajectory by reintegrating two half steps.
// Status concerns completion of this numerical step, not a certificate of physical accuracy.
#include <cstddef>

namespace rf {

class RigidWorld;

enum class RigidStepStatus { Completed, InvalidInput, Unresolved, SubdivisionLimit, NonFinite,
                             Unsupported, PersistentContact, EventLimit };

struct RigidStepResult {
    RigidStepStatus status = RigidStepStatus::Completed;
    float requestedDt = 0;
    float acceptedDt = 0;
    size_t acceptedSubsteps = 0;
    size_t rejectedTrials = 0;
    size_t impactEvents = 0; // instantaneous coupled solves in the experimental variational mode
    float maxContactDepth = 0; // maximum manifold depth over accepted subdivisions, not a geometry certificate
    double motorWork = 0; // sum of signed actuator impulse work over accepted subdivisions [J]
    bool completed() const { return status == RigidStepStatus::Completed; }
};

class RigidStep {
public:
    explicit RigidStep(RigidWorld& world) : world_(world) {}
    RigidStepResult advance(float dt);

private:
    RigidWorld& world_;
};

} // namespace rf
