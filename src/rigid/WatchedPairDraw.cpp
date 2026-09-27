// The research view of one watched body pair (Probe::watchPair): how GJK and EPA see it. The pair
// is run once more through GJK (+ EPA when it overlaps) with a GjkTrace recording, whatever path
// the narrow phase took for it (box-box goes through SAT): the picture shows the Minkowski
// difference A - B that the contact is about.
//   GjkSimplex:    the simplex of every GJK iteration, drawn in an INSET: the Minkowski space
//                  scaled so its largest point is 0.3 m from the inset's origin (a white cross,
//                  the origin of A - B), placed above the pair; early iterations blue, late red.
//                  In the world: the support points a on A (green) and b on B (blue) of the last
//                  simplex.
//   EpaPolytope:   the final EPA polytope's edges in the same inset (orange).
//   WitnessPoints: the closest points and the distance (separated), or the deepest points and
//                  the penetration vector normal * depth (overlapping), with a label.
#include "rigid/RigidWorld.h"

#include "core/Format.h"
#include "core/Probe.h"
#include "rigid/GjkEpa.h"

#include <algorithm>

namespace rf {

namespace {

// Where the Minkowski-space inset sits and how it is scaled.
struct Inset {
    Vector3 origin{0.0f};
    float scale = 1;
    Vector3 map(const Vector3& w) const { return origin + w * scale; }
};

Inset makeInset(const RigidBody& A, const RigidBody& B, const GjkTrace& t) {
    float extent = 1e-6f;
    for (const auto& s : t.simplices)
        for (const auto& v : s) extent = std::max(extent, length(v.w));
    for (const auto& v : t.polytope) extent = std::max(extent, length(v.w));
    Inset inset;
    inset.origin = (A.pos + B.pos) * 0.5f + Vector3(0.0f, std::max(A.boundingRadius(), B.boundingRadius()) * 1.2f + 0.35f, 0.0f);
    inset.scale = 0.3f / extent;
    return inset;
}

void drawSimplices(const GjkTrace& t, const Inset& inset) {
    const float r = 0.02f;
    for (int axis = 0; axis < 3; ++axis) { // the origin of A - B: the shapes touch when it is reached
        Vector3 e(0.0f);
        e[axis] = r;
        Probe::line(DrawLayer::GjkSimplex, inset.origin - e, inset.origin + e, Vector3(1.0f));
    }
    const int n = int(t.simplices.size());
    for (int k = 0; k < n; ++k) {
        const auto& s = t.simplices[size_t(k)];
        const Vector3 col = heatColor(n > 1 ? float(k) / float(n - 1) : 1.0f);
        for (size_t i = 0; i < s.size(); ++i) {
            Probe::point(DrawLayer::GjkSimplex, inset.map(s[i].w), col, 0.01f);
            for (size_t j = i + 1; j < s.size(); ++j) Probe::line(DrawLayer::GjkSimplex, inset.map(s[i].w), inset.map(s[j].w), col);
        }
    }
    if (n == 0) return;
    for (const auto& v : t.simplices.back()) { // the support points in the world that made the last simplex
        Probe::point(DrawLayer::GjkSimplex, v.a, Vector3(0.2f, 1.0f, 0.3f), 0.015f);
        Probe::point(DrawLayer::GjkSimplex, v.b, Vector3(0.2f, 0.5f, 1.0f), 0.015f);
        Probe::line(DrawLayer::GjkSimplex, v.a, v.b, Vector3(0.6f));
    }
    Probe::label(DrawLayer::GjkSimplex, inset.origin + Vector3(0.0f, 0.36f, 0.0f), format("GJK: %d итераций", n));
}

void drawPolytope(const GjkTrace& t, const Inset& inset) {
    const Vector3 col(1.0f, 0.6f, 0.1f);
    for (size_t f = 0; f + 2 < t.faces.size(); f += 3)
        for (int e = 0; e < 3; ++e) {
            const Vector3 a = t.polytope[size_t(t.faces[f + size_t(e)])].w;
            const Vector3 b = t.polytope[size_t(t.faces[f + size_t((e + 1) % 3)])].w;
            Probe::line(DrawLayer::EpaPolytope, inset.map(a), inset.map(b), col);
        }
    if (!t.faces.empty())
        Probe::label(DrawLayer::EpaPolytope, inset.origin - Vector3(0.0f, 0.36f, 0.0f),
                     format("EPA: %zu граней, +%d вершин", t.faces.size() / 3, t.epaIterations));
}

void drawWitness(const GjkResult& g, const PenetrationResult& pen) {
    if (!g.intersect) {
        Probe::point(DrawLayer::WitnessPoints, g.pointA, Vector3(0.2f, 1.0f, 0.3f), 0.02f);
        Probe::point(DrawLayer::WitnessPoints, g.pointB, Vector3(0.2f, 0.5f, 1.0f), 0.02f);
        Probe::line(DrawLayer::WitnessPoints, g.pointA, g.pointB, Vector3(1.0f, 1.0f, 0.3f));
        Probe::label(DrawLayer::WitnessPoints, (g.pointA + g.pointB) * 0.5f, format("зазор %.4f м", double(g.distance)));
        return;
    }
    if (!pen.valid) return;
    Probe::point(DrawLayer::WitnessPoints, pen.pointA, Vector3(0.2f, 1.0f, 0.3f), 0.02f);
    Probe::point(DrawLayer::WitnessPoints, pen.pointB, Vector3(0.2f, 0.5f, 1.0f), 0.02f);
    Probe::arrow(DrawLayer::WitnessPoints, pen.pointA, pen.normal * pen.depth, Vector3(1.0f, 0.15f, 0.15f));
    Probe::label(DrawLayer::WitnessPoints, pen.pointA, format("глубина %.4f м", double(pen.depth)));
}

} // namespace

void RigidWorld::drawWatchedPair() const {
    watchReport_ = WatchReport();
    const bool simplex = Probe::layerOn(DrawLayer::GjkSimplex), poly = Probe::layerOn(DrawLayer::EpaPolytope);
    const bool witness = Probe::layerOn(DrawLayer::WitnessPoints);
    if (!simplex && !poly && !witness) return;
    int ia = -1, ib = -1;
    Probe::watchedPair(ia, ib);
    const int n = int(bodies_.size());
    if (ia < 0 || ib < 0 || ia >= n || ib >= n || ia == ib || !bodies_[size_t(ia)].alive || !bodies_[size_t(ib)].alive) return;
    const RigidBody &A = bodies_[size_t(ia)], &B = bodies_[size_t(ib)];
    GjkTrace trace;
    setGjkTrace(&trace);
    const GjkResult g = gjk(A.posed(), B.posed());
    PenetrationResult pen;
    if (g.intersect) pen = epa(A.posed(), B.posed(), g);
    setGjkTrace(nullptr);
    watchReport_ = {true, g.intersect, g.intersect ? 0.0f : g.distance, pen.valid ? pen.depth : 0.0f,
                    pen.valid ? pen.normal : Vector3(0.0f), int(trace.simplices.size()), int(trace.faces.size() / 3)};
    const Inset inset = makeInset(A, B, trace);
    if (simplex) drawSimplices(trace, inset);
    if (poly) drawPolytope(trace, inset);
    if (witness) drawWitness(g, pen);
}

} // namespace rf
