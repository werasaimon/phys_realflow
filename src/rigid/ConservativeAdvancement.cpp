// Conservative advancement with explicit uncertainty and no forced minimum time increment.
// A separating support plane, unlike a finite GJK simplex distance, is a lower distance bound
// in exact arithmetic. Float geometry still needs a certified replacement for an IPC guarantee.
#include "rigid/ConservativeAdvancement.h"

#include <limits>

namespace rf {

namespace {

bool finiteVector(const Vector3& v) {
    return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

double dotDouble(const Vector3& a, const Vector3& b) {
    return double(a.x) * b.x + double(a.y) * b.y + double(a.z) * b.z;
}

// Engineering guard for float supports/transforms, not a proved outward-rounded interval.
double roundingGuard(const Vector3& a, const Vector3& b) {
    const double scale = double(maxComp(vabs(a))) + double(maxComp(vabs(b)));
    return 8.0 * std::numeric_limits<float>::epsilon() * scale;
}

} // namespace

bool ConservativeAdvancement::validSweep(const SweptPose& sweep) {
    if (!sweep.shape || !finiteVector(sweep.p0) || !finiteVector(sweep.p1)
        || !finiteVector(sweep.dTheta) || !finiteVector(sweep.partT) || !finiteVector(sweep.translationCurve)) return false;
    const Quaternion& q = sweep.q0;
    const double norm2 = double(q.w) * q.w + double(q.x) * q.x + double(q.y) * q.y + double(q.z) * q.z;
    if (!std::isfinite(norm2) || std::fabs(norm2 - 1.0) > 1e-3) return false;
    for (const auto& row : sweep.partR.m)
        for (float value : row)
            if (!std::isfinite(value)) return false;
    return true; // valid support shapes and a rigid partR remain API preconditions
}

bool ConservativeAdvancement::begin(float motionBound) {
    fraction_ = 0;
    checkedFraction_ = 0;
    motionBound_ = motionBound;
    result_ = ToiResult();
    if (!std::isfinite(tolerance_) || tolerance_ <= 0 || !std::isfinite(motionBound_) || motionBound_ < 0) {
        unresolved(ToiReason::InvalidInput);
        return false;
    }
    return true;
}

ToiResult ConservativeAdvancement::unresolved(ToiReason reason) {
    result_.status = ToiStatus::Unresolved;
    result_.reason = reason;
    result_.hit = false;
    return result_;
}

void ConservativeAdvancement::touch(bool separatedSample) {
    result_.s = separatedSample ? fraction_ : checkedFraction_;
    result_.status = fraction_ == 0 ? ToiStatus::InitialContact : ToiStatus::Impact;
    result_.reason = ToiReason::None;
    result_.hit = result_.status == ToiStatus::Impact;
    // Keep the last separated fraction if this sample already overlaps.
}

ConservativeAdvancement::Separation ConservativeAdvancement::supportGap(const PosedShape& a, const PosedShape& b, const Vector3& direction) {
    const double norm = std::sqrt(dotDouble(direction, direction));
    if (!(norm > 0) || !std::isfinite(norm)) return {-1, Vector3(0)};
    const Vector3 n(float(direction.x / norm), float(direction.y / norm), float(direction.z / norm));
    const Vector3 pa = a.support(-n), pb = b.support(n);
    if (!finiteVector(pa) || !finiteVector(pb)) return {std::numeric_limits<double>::quiet_NaN(), n};
    const double gap = dotDouble(pa, n) - dotDouble(pb, n);
    return {gap - roundingGuard(pa, pb), n};
}

ConservativeAdvancement::Separation ConservativeAdvancement::separation(const PosedShape& a, const PosedShape& b, const Vector3& direction) {
    Separation best = supportGap(a, b, direction);
    if (!std::isfinite(best.gap)) return best;
    // A triangle's face normal is an exact SAT candidate in real arithmetic. Near a large face,
    // subtracting GJK witnesses amplifies float noise; the face plane avoids that direction error.
    for (const PosedShape* shape : {&a, &b}) {
        if (shape->shape->type() != ShapeType::Triangle) continue;
        const Vector3 axis = shape->R * static_cast<const TriangleShape*>(shape->shape)->normal();
        for (float sign : {-1.0f, 1.0f}) {
            const Separation candidate = supportGap(a, b, sign * axis);
            if (!std::isfinite(candidate.gap)) return candidate;
            if (candidate.gap > best.gap) best = candidate;
        }
    }
    return best;
}

ConservativeAdvancement::PlaneSample ConservativeAdvancement::planeGap(const PosedShape& a, const Vector3& normal, double offset) {
    const Vector3 p = a.support(-normal);
    if (!finiteVector(p)) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return {nan, nan};
    }
    const double gap = dotDouble(p, normal) - offset;
    return {gap, gap - roundingGuard(p, normal * float(offset))};
}

// OBB SAT candidates: six face normals and nine edge cross products (Gottschalk et al.,
// OBBTree, 1996). A poorly converged GJK direction can hide a real separating plane between
// a thin beam and a cube. Each candidate is verified with supports, never used as an upper bound.
ConservativeAdvancement::Separation ConservativeAdvancement::boxSeparation(const PosedShape& a, const PosedShape& b, Separation best) {
    if (a.shape->type() != ShapeType::Box || b.shape->type() != ShapeType::Box) return best;
    auto consider = [&](const Vector3& axis) {
        for (float sign : {-1.0f, 1.0f}) {
            const Separation candidate = supportGap(a, b, axis * sign);
            if (candidate.gap > best.gap) best = candidate;
        }
    };
    for (int i = 0; i < 3; ++i) {
        consider(a.R.col(i));
        consider(b.R.col(i));
        for (int j = 0; j < 3; ++j) consider(cross(a.R.col(i), b.R.col(j)));
    }
    return best;
}

// For a sphere outside a box, the box's closest-point normal points to the sphere centre.
// Convexity makes it a separating-plane candidate; supports verify the actual surface gap.
// This also covers corner/edge contacts where a box face normal alone cannot separate them.
ConservativeAdvancement::Separation ConservativeAdvancement::sphereBoxSeparation(const PosedShape& a, const PosedShape& b, Separation best) {
    const bool sphereFirst = a.shape->type() == ShapeType::Sphere && b.shape->type() == ShapeType::Box;
    const bool boxFirst = a.shape->type() == ShapeType::Box && b.shape->type() == ShapeType::Sphere;
    if (!sphereFirst && !boxFirst) return best;
    const PosedShape& box = sphereFirst ? b : a;
    const PosedShape& sphere = sphereFirst ? a : b;
    Vector3 normal;
    const float distance = box.shape->signedDistance(box.R.transposed() * (sphere.p - box.p), normal);
    if (!(distance > 0)) return best;
    const Separation candidate = supportGap(a, b, box.R * normal * (sphereFirst ? 1.0f : -1.0f));
    return candidate.gap > best.gap ? candidate : best;
}

bool ConservativeAdvancement::advance(double gap) {
    if (!std::isfinite(gap)) { unresolved(ToiReason::InvalidInput); return true; }
    if (!(gap > 0)) { unresolved(ToiReason::NoSeparatingPlane); return true; }
    checkedFraction_ = fraction_; // a positive supporting-plane gap was found here
    result_.s = fraction_;
    // A small LOWER bound does not prove contact. Only the caller's upper-distance test does.
    if (gap <= 0.5 * tolerance_) { unresolved(ToiReason::NoSeparatingPlane); return true; }
    // Mirtich conservative advancement: the gap cannot shrink faster than motionBound_.
    // Half the contact tolerance is reserved as a positive gap at the next sample.
    const double remaining = 1.0 - fraction_;
    const double distance = gap - 0.5 * tolerance_;
    if (motionBound_ == 0 || distance > double(motionBound_) * remaining) {
        result_.status = ToiStatus::Separated;
        result_.reason = ToiReason::None;
        result_.s = 1;
        return true;
    }
    const double candidate = double(fraction_) + distance / motionBound_;
    // Round the parameter toward the previous sample; never enlarge a tiny safe advance.
    const float next = std::nextafter(float(candidate), fraction_);
    if (!(next > fraction_)) { unresolved(ToiReason::NoProgress); return true; }
    fraction_ = next;
    result_.s = next; // the preceding plane and motion bound cover this entire prefix
    return false;
}

ToiResult ConservativeAdvancement::between(const SweptPose& a, const SweptPose& b) {
    begin(0); // reset even when a reused query receives invalid input
    if (!validSweep(a) || !validSweep(b) || a.shape->type() == ShapeType::Compound
        || b.shape->type() == ShapeType::Compound) return unresolved(ToiReason::InvalidInput);
    const Vector3 translation = (a.p1 - a.p0) - (b.p1 - b.p0);
    const Vector3 curve = a.translationCurve - b.translationCurve;
    const float rotation = a.angularReach() + b.angularReach();
    if (!begin(length(translation) + length(curve) + rotation)) return result_;
    for (int it = 0; it < maxIterations_; ++it) {
        result_.iterations = it + 1;
        const PosedShape pa = a.at(fraction_), pb = b.at(fraction_);
        const GjkResult g = gjk(pa, pb);
        if (!std::isfinite(g.distance)) return unresolved(ToiReason::InvalidInput);
        if (allowDeparture_ && fraction_ == 0) {
            Separation start = separation(pa, pb, pa.p - pb.p);
            start = sphereBoxSeparation(pa, pb, start);
            if (departure(start.gap, start.normal, translation, curve, rotation)) return result_;
        }
        if (g.intersect) { touch(); return result_; }
        Separation plane = separation(pa, pb, g.pointA - g.pointB);
        if (!std::isfinite(plane.gap)) return unresolved(ToiReason::InvalidInput);
        if (g.distance <= tolerance_ && plane.gap <= tolerance_) { touch(true); return result_; }
        if (plane.gap <= 0.5 * tolerance_) plane = boxSeparation(pa, pb, plane);
        if (plane.gap <= 0.5 * tolerance_) plane = sphereBoxSeparation(pa, pb, plane);
        // Only translation toward this fixed separating plane closes its gap; tangential sliding
        // does not. Rotation is bounded by angular speed times reach (Mirtich 1996, docs/book/09).
        motionBound_ = float(std::max(0.0, -dotDouble(translation, plane.normal)
                                      + std::fabs(dotDouble(curve, plane.normal)))) + rotation;
        if (advance(plane.gap)) return result_;
    }
    return unresolved(ToiReason::IterationLimit);
}

ToiResult ConservativeAdvancement::againstPlane(const SweptPose& a, const Vector3& normal, float offset) {
    begin(0);
    const double norm = std::sqrt(dotDouble(normal, normal));
    if (!validSweep(a) || !std::isfinite(offset) || !std::isfinite(norm) || !(norm > 0))
        return unresolved(ToiReason::InvalidInput);
    const Vector3 n(float(normal.x / norm), float(normal.y / norm), float(normal.z / norm));
    const double d = double(offset) / norm;
    const float bound = float(std::max(0.0, -dotDouble(a.p1 - a.p0, n)
                                      + std::fabs(dotDouble(a.translationCurve, n)))) + a.angularReach();
    if (!begin(bound)) return result_;
    for (int it = 0; it < maxIterations_; ++it) {
        result_.iterations = it + 1;
        const PlaneSample sample = planeGap(a.at(fraction_), n, d);
        if (fraction_ == 0 && departure(sample.guardedGap, n, a.p1 - a.p0, a.translationCurve, a.angularReach())) return result_;
        if (std::isfinite(sample.gap) && sample.gap <= tolerance_) { touch(sample.gap > 0); return result_; }
        if (advance(sample.guardedGap)) return result_;
    }
    return unresolved(ToiReason::IterationLimit);
}

// For a fixed support normal, delta gap(s) = s*((d-c).n + s*c.n). Both
// endpoint coefficients nonnegative imply no return below the initial gap.
// The small negative allowance is the existing geometric tolerance, not certification.
bool ConservativeAdvancement::departure(double gap, const Vector3& normal, const Vector3& translation,
                                        const Vector3& curve, float rotation) {
    if (!allowDeparture_ || rotation != 0 || !std::isfinite(gap) || gap < -0.25 * tolerance_
        || length2(normal) < 0.5f || dotDouble(translation - curve, normal) < 0
        || dotDouble(translation, normal) < 0) return false;
    result_.status = ToiStatus::Separated; result_.reason = ToiReason::None; result_.s = 1;
    return true;
}

} // namespace rf
