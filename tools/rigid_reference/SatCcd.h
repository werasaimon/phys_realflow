// Experimental GJK candidate ordering, complete polytope SAT and rotational CCD.
// Built by the isolated experiment and optional backend checks; no default solver change.
// SAT: https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf, sections 4/6.
#pragma once

#include "rigid/NarrowPhase.h"
#include "rigid/TimeOfImpact.h"
#include <iosfwd>

namespace rf::experimental {

struct SatResult {
    bool valid = false;
    double gap = -1e100; // positive: separating projection; negative: convex translation depth
    Vector3 normal;     // B -> A
    int axes = 0;
};

enum class SweepStatus { Clear, Hit, InitialContact, Unresolved, Unsupported };
struct SweepResult {
    SweepStatus status = SweepStatus::Clear;
    float s = 1;
    int iterations = 0;
};

bool isPolytope(const ConvexShape* shape);
bool narrowEnabled();
bool ccdEnabled();
SatResult polytopeSat(const PosedShape& a, const PosedShape& b, double rejectGap = 1e100, bool allEdgeAxes = false);
bool collideSat(const PosedShape& a, const PosedShape& b, ContactManifold& manifold);
SweepResult sweptSat(const SweptPose& a, const SweptPose& b, float tolerance, int budget = 64);
ToiResult conservativeToi(const SweptPose& a, const SweptPose& b, float tolerance, int budget = 64);
ToiResult compoundToi(const SweptPose& a, const SweptPose& b, float tolerance);
void writeStats(std::ostream& stream);

} // namespace rf::experimental
