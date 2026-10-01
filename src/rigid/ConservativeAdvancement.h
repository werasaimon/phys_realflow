#pragma once
// Stateful continuous collision query for one pair: GJK selects a direction, supporting planes
// bound separation, and a bound on rigid motion limits each advance (Mirtich, 1996). The object
// owns the last separated fraction and every termination status; it never updates world bodies.
#include "rigid/TimeOfImpact.h"

namespace rf {

class ConservativeAdvancement {
public:
    explicit ConservativeAdvancement(float tolerance, int maxIterations = 64, bool allowDeparture = false)
        : tolerance_(tolerance), maxIterations_(maxIterations), allowDeparture_(allowDeparture) {}

    ToiResult between(const SweptPose& a, const SweptPose& b);
    ToiResult againstPlane(const SweptPose& a, const Vector3& normal, float offset);

private:
    struct PlaneSample { double gap; double guardedGap; };
    struct Separation { double gap; Vector3 normal; };
    float tolerance_;
    int maxIterations_;
    bool allowDeparture_;
    float fraction_ = 0;
    float checkedFraction_ = 0;
    float motionBound_ = 0;
    ToiResult result_;

    bool begin(float motionBound);
    bool advance(double gap);
    void touch(bool separatedSample = false);
    ToiResult unresolved(ToiReason reason);
    bool departure(double gap, const Vector3& normal, const Vector3& translation, const Vector3& curve, float rotation);
    static bool validSweep(const SweptPose& sweep);
    static Separation supportGap(const PosedShape& a, const PosedShape& b, const Vector3& direction);
    static Separation separation(const PosedShape& a, const PosedShape& b, const Vector3& direction);
    static Separation boxSeparation(const PosedShape& a, const PosedShape& b, Separation best);
    static Separation sphereBoxSeparation(const PosedShape& a, const PosedShape& b, Separation best);
    static PlaneSample planeGap(const PosedShape& a, const Vector3& normal, double offset);
};

} // namespace rf
