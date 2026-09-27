// The contact manifolds of the narrow phase: for every pair of shape types the points, normals and
// depths where they touch - analytic for spheres, SAT plus clipping for boxes, GJK/EPA plus face
// clipping for convex hulls and triangles. The dispatch table and the overview are in NarrowPhase.h.
#include "rigid/NarrowPhase.h"

#include "core/Parallel.h"

namespace rf {

namespace {

// A point (the core of a sphere) for GJK.
class PointShape final : public ConvexShape {
public:
    ShapeType type() const override { return ShapeType::Sphere; }
    Vector3 support(const Vector3&) const override { return Vector3(0.0f); }
    AABB localBounds() const override { return AABB(Vector3(0.0f), Vector3(0.0f)); }
    float volume() const override { return 0; }
    Vector3 unitInertia() const override { return Vector3(0.0f); }
    float signedDistance(const Vector3& p, Vector3& n) const override { n = normalize(p); return length(p); }
    float boundingRadius() const override { return 0; }
};
const PointShape kPoint;

// The core of a capsule: its axis segment (0, -h, 0) .. (0, h, 0), for GJK/EPA; the radius is added
// to the result as for a sphere's core point.
class SegmentShape final : public ConvexShape {
public:
    explicit SegmentShape(float h) : h_(h) {}
    ShapeType type() const override { return ShapeType::Capsule; }
    Vector3 support(const Vector3& d) const override { return Vector3(0.0f, d.y >= 0 ? h_ : -h_, 0.0f); }
    AABB localBounds() const override { return AABB(Vector3(0.0f, -h_, 0.0f), Vector3(0.0f, h_, 0.0f)); }
    float volume() const override { return 0; }
    Vector3 unitInertia() const override { return Vector3(0.0f); }
    float signedDistance(const Vector3& p, Vector3& n) const override {
        const Vector3 d = p - Vector3(0.0f, clampv(p.y, -h_, h_), 0.0f);
        n = normalize(d);
        return length(d);
    }
    float boundingRadius() const override { return h_; }

private:
    float h_;
};

// A posed capsule's axis segment in the world and its radius.
void capsuleAxis(const PosedShape& s, Vector3& a, Vector3& b, float& r) {
    const auto* c = static_cast<const CapsuleShape*>(s.shape);
    const Vector3 up = s.R * Vector3(0.0f, c->halfHeight(), 0.0f);
    a = s.p - up;
    b = s.p + up;
    r = c->radius();
}

Vector3 closestOnSegment(const Vector3& a, const Vector3& b, const Vector3& p) {
    const Vector3 ab = b - a;
    const float l2 = dot(ab, ab);
    return l2 < 1e-12f ? a : a + ab * clampv(dot(p - a, ab) / l2, 0.0f, 1.0f);
}

float sphereRadius(const PosedShape& s) { return static_cast<const SphereShape*>(s.shape)->radius(); }

void flipAppend(const ContactManifold& src, ContactManifold& dst) {
    for (const ContactPoint& c : src.points) dst.add(c.position, -c.normal, c.depth);
}

// Closest points between segments p1-q1 and p2-q2 (Ericson 5.1.9).
void closestSegmentSegment(const Vector3& p1, const Vector3& q1, const Vector3& p2, const Vector3& q2, Vector3& c1, Vector3& c2) {
    Vector3 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
    float a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    float s, t;
    if (a <= 1e-12f && e <= 1e-12f) { c1 = p1; c2 = p2; return; }
    if (a <= 1e-12f) { s = 0; t = clampv(f / e, 0.0f, 1.0f); }
    else {
        float c = dot(d1, r);
        if (e <= 1e-12f) { t = 0; s = clampv(-c / a, 0.0f, 1.0f); }
        else {
            float b = dot(d1, d2), den = a * e - b * b;
            s = den != 0 ? clampv((b * f - c * e) / den, 0.0f, 1.0f) : 0.0f;
            t = (b * s + f) / e;
            if (t < 0) { t = 0; s = clampv(-c / a, 0.0f, 1.0f); }
            else if (t > 1) { t = 1; s = clampv((b - c) / a, 0.0f, 1.0f); }
        }
    }
    c1 = p1 + d1 * s;
    c2 = p2 + d2 * t;
}

// Sutherland-Hodgman: keep the part of `poly` with dot(n, x) <= d, into `out` (the caller's
// scratch, cleared here: no allocation once it has grown).
void clipPolygon(const std::vector<Vector3>& poly, const Vector3& n, float d, std::vector<Vector3>& out) {
    out.clear();
    if (poly.empty()) return;
    for (size_t i = 0; i < poly.size(); ++i) {
        const Vector3& a = poly[i];
        const Vector3& b = poly[(i + 1) % poly.size()];
        float da = dot(n, a) - d, db = dot(n, b) - d;
        if (da <= 0) out.push_back(a);
        if ((da <= 0) != (db <= 0)) out.push_back(a + (b - a) * (da / (da - db)));
    }
}

Vector3 perpendicular(const Vector3& n) {
    return normalize(std::fabs(n.x) > 0.57f ? Vector3(n.y, -n.x, 0) : Vector3(0, n.z, -n.y));
}

} // namespace

// ---------------------------------------------------------------------------
void reduceManifold(std::vector<ContactPoint>& pts, size_t maxPoints) {
    // As Jolt's PruneContactPoints / Bullet's sortCachedPoints: 1) the deepest point - among
    // points equally deep (flat contact) the one farthest from the centre of the region, a corner,
    // 2) the point farthest from it, 3) the point farthest from their line on one side, 4) the
    // farthest on the other side (signed area about the contact normal) - the quadrilateral keeps
    // the deepest point and spans the support region whatever order the points came in.
    if (pts.size() <= maxPoints) return;
    Vector3 centre(0.0f);
    float maxDepth = -kInf;
    for (const ContactPoint& p : pts) {
        centre += p.position;
        maxDepth = std::max(maxDepth, p.depth);
    }
    centre /= float(pts.size());
    const float tie = 1e-4f * std::fabs(maxDepth) + 1e-6f; // equal up to rounding: a flat contact
    size_t i0 = 0;
    float far0 = -1;
    for (size_t i = 0; i < pts.size(); ++i) {
        if (pts[i].depth < maxDepth - tie) continue;
        const float d = length2(pts[i].position - centre);
        if (d > far0) { far0 = d; i0 = i; }
    }
    size_t i1 = i0;
    float best = -1;
    for (size_t i = 0; i < pts.size(); ++i) {
        float d = length2(pts[i].position - pts[i0].position);
        if (d > best) { best = d; i1 = i; }
    }
    // The kept points, in this order, written back into `pts` (shrinking never allocates).
    ContactPoint kept[4];
    int nk = 0;
    kept[nk++] = pts[i0];
    if (i1 != i0) kept[nk++] = pts[i1];
    if (maxPoints >= 4 && i1 != i0) {
        const Vector3 n = pts[i0].normal, a = pts[i0].position, b = pts[i1].position;
        size_t left = pts.size(), right = pts.size();
        float maxLeft = 1e-9f, maxRight = -1e-9f;
        for (size_t i = 0; i < pts.size(); ++i) {
            const float s = dot(cross(b - a, pts[i].position - a), n); // twice the signed triangle area
            if (s > maxLeft) { maxLeft = s; left = i; }
            if (s < maxRight) { maxRight = s; right = i; }
        }
        if (left != pts.size()) kept[nk++] = pts[left];
        if (right != pts.size()) kept[nk++] = pts[right];
    }
    pts.resize(size_t(nk));
    for (int i = 0; i < nk; ++i) pts[size_t(i)] = kept[i];
}

// Scratch of the narrow phase, one per thread of the pool: the vectors are cleared before use
// and never freed, so a collision costs no allocation. Indexed by the worker number rather than
// held in thread_local objects, whose destructors MinGW's TLS cleanup double-frees at thread exit.
namespace {
struct NarrowScratch {
    std::vector<Vector3> poly, clipped;    // box - box: the incident face and its clipped copy
    std::vector<ContactPoint> boxPoints;   // box - box: the points below the reference face
    std::vector<Vector3> fa, fb, inc, cut; // faceManifold: the two faces, the incident one clipped
    std::vector<ContactPoint> hullPoints;  // convex - convex: the face manifold or the single point
    std::vector<PosedShape> partsA, partsB; // compound: the convex parts of both sides
    std::vector<AABB> partBounds;
    ContactManifold flipped;               // a pair handled with the shapes swapped
};
NarrowScratch& scratch() {
    static std::vector<NarrowScratch> all(size_t(ThreadPool::instance().threadCount()));
    return all[size_t(ThreadPool::workerIndex())];
}
} // namespace

// ---------------------------------------------------------------------------
NarrowPhase::NarrowPhase() {
    const int S = int(ShapeType::Sphere), B = int(ShapeType::Box), H = int(ShapeType::ConvexHull), T = int(ShapeType::Triangle);
    table_[S][S] = &sphereSphere;
    table_[S][B] = &sphereBox;
    table_[S][H] = &sphereConvex;
    table_[S][T] = &sphereConvex;
    table_[B][B] = &boxBox;
    table_[B][H] = &convexConvex;
    table_[B][T] = &convexConvex;
    table_[H][H] = &convexConvex;
    table_[H][T] = &convexConvex;
    const int C = int(ShapeType::Capsule);
    table_[S][C] = &sphereCapsule;
    table_[C][C] = &capsuleCapsule;
    table_[C][B] = &capsuleConvex;
    table_[C][H] = &capsuleConvex;
    table_[C][T] = &capsuleConvex;
    // Reverse pairs are handled by swapping in collide().
}

bool NarrowPhase::collide(const PosedShape& A, const PosedShape& B, ContactManifold& m) const {
    // Compound (non-convex) bodies: every convex part against the other side, culled by bounds.
    auto partsOf = [](const PosedShape& S, std::vector<PosedShape>& out) {
        if (S.shape->type() != ShapeType::Compound) { out.push_back(S); return; }
        for (const auto& c : static_cast<const CompoundShape*>(S.shape)->children())
            out.push_back({c.shape.get(), S.R * c.R, S.p + S.R * c.t});
    };
    if (A.shape->type() == ShapeType::Compound || B.shape->type() == ShapeType::Compound) {
        // The worker's scratch (the parts are convex, so the recursive calls below never come
        // back here and never touch these): no allocation per compound pair.
        NarrowScratch& S = scratch();
        std::vector<PosedShape>&pa = S.partsA, &pb = S.partsB;
        std::vector<AABB>& bb = S.partBounds;
        pa.clear();
        pb.clear();
        bb.clear();
        partsOf(A, pa);
        partsOf(B, pb);
        auto boundsOf = [](const PosedShape& S) {
            AABB bb = S.shape->type() == ShapeType::Sphere ? S.shape->boundsAt(S.R, S.p)
                                                           : orientedBounds(S.shape->localBounds(), S.R, S.p);
            bb.lo -= Vector3(margin);
            bb.hi += Vector3(margin);
            return bb;
        };
        for (const PosedShape& q : pb) bb.push_back(boundsOf(q));
        bool hit = false;
        for (const PosedShape& p : pa) {
            AABB ab = boundsOf(p);
            for (size_t j = 0; j < pb.size(); ++j)
                if (ab.overlaps(bb[j])) hit |= collide(p, pb[j], m);
        }
        return hit;
    }
    int a = int(A.shape->type()), b = int(B.shape->type());
    if (table_[a][b]) return table_[a][b](A, B, m);
    if (table_[b][a]) {
        ContactManifold& tmp = scratch().flipped;
        tmp.points.clear();
        bool hit = table_[b][a](B, A, tmp);
        flipAppend(tmp, m);
        return hit;
    }
    return false;
}

bool NarrowPhase::sphereSphere(const PosedShape& A, const PosedShape& B, ContactManifold& m) {
    float ra = sphereRadius(A), rb = sphereRadius(B);
    Vector3 d = A.p - B.p;
    float l = length(d);
    float depth = ra + rb - l;
    if (depth <= -margin) return false;
    Vector3 n = l > 1e-9f ? d / l : Vector3(0, 1, 0);
    m.add(B.p + n * (rb - 0.5f * depth), n, depth);
    return true;
}

bool NarrowPhase::sphereBox(const PosedShape& A, const PosedShape& B, ContactManifold& m) {
    float r = sphereRadius(A);
    Vector3 q = B.R.transposed() * (A.p - B.p);
    Vector3 nl;
    float sd = static_cast<const BoxShape*>(B.shape)->signedDistance(q, nl);
    float depth = r - sd;
    if (depth <= -margin) return false;
    Vector3 n = B.R * nl;
    m.add(A.p - n * (0.5f * (r + sd)), n, depth);
    return true;
}

bool NarrowPhase::sphereConvex(const PosedShape& A, const PosedShape& B, ContactManifold& m) {
    float r = sphereRadius(A);
    PosedShape core{&kPoint, A.R, A.p};
    GjkResult g = gjk(core, B, r + margin);
    if (!g.intersect) {
        if (g.distance >= r + margin || g.distance < 1e-9f) return false;
        Vector3 n = (A.p - g.pointB) / g.distance;
        float depth = r - g.distance;
        m.add(g.pointB + n * (0.5f * (g.distance - r)), n, depth); // midway between the two surfaces
        return true;
    }
    PenetrationResult pr = epa(core, B, g);
    if (!pr.valid) return false;
    float depth = pr.depth + r;
    m.add(A.p - pr.normal * (0.5f * (r + pr.depth)), pr.normal, depth);
    return true;
}

// ---------------------------------------------------------------------------
// Box-box: separating axis theorem
// ---------------------------------------------------------------------------
// Box against box: the separating axis theorem over the 15 candidate axes (3 + 3 face normals,
// 9 edge cross products) finds the axis of least overlap; a face axis gives a face contact by
// clipping (Sutherland-Hodgman), an edge axis a single point between the two edges.
bool NarrowPhase::boxBox(const PosedShape& A, const PosedShape& B, ContactManifold& m) {
    BoxAxis best, bestEdge;
    if (!boxSeparatingAxes(A, B, best, bestEdge)) return false;
    if (best.kind == 2) {
        boxEdgeContact(A, B, best, m);
        return true;
    }
    std::vector<ContactPoint>& pts = scratch().boxPoints; // the worker's scratch: no allocation per pair
    pts.clear();
    if (!boxFaceContact(A, B, best, pts)) {
        if (bestEdge.sep != -kInf) { boxEdgeContact(A, B, bestEdge, m); return true; }
        return false;
    }
    reduceManifold(pts, 4);
    for (auto& p : pts) m.points.push_back(p);
    return true;
}

// SAT: project both boxes on every candidate axis; the first axis with a gap larger than the
// margin separates them (false). Otherwise the axis of least penetration wins, with face axes
// preferred - another axis must be clearly better - so the manifold stays the same face from
// step to step (stable stacking). The best edge axis is kept apart as a fallback.
bool NarrowPhase::boxSeparatingAxes(const PosedShape& A, const PosedShape& B, BoxAxis& best, BoxAxis& bestEdge) {
    const Vector3 hA = static_cast<const BoxShape*>(A.shape)->halfExtents();
    const Vector3 hB = static_cast<const BoxShape*>(B.shape)->halfExtents();
    const Vector3 a[3] = {A.R.col(0), A.R.col(1), A.R.col(2)};
    const Vector3 b[3] = {B.R.col(0), B.R.col(1), B.R.col(2)};
    const Vector3 T = B.p - A.p;
    const float scale = std::min(minComp(hA), minComp(hB));
    auto test = [&](Vector3 L, int kind, int i, int j) {
        float len = length(L);
        if (len < 1e-5f) return true; // parallel edges: axis degenerate, covered by face axes
        L /= len;
        float rA = hA.x * std::fabs(dot(a[0], L)) + hA.y * std::fabs(dot(a[1], L)) + hA.z * std::fabs(dot(a[2], L));
        float rB = hB.x * std::fabs(dot(b[0], L)) + hB.y * std::fabs(dot(b[1], L)) + hB.z * std::fabs(dot(b[2], L));
        float dist = dot(T, L);
        float sep = std::fabs(dist) - (rA + rB);
        if (sep > margin) return false; // separating axis found
        Vector3 axis = dist < 0 ? -L : L;
        // Face axes are preferred: another axis must be clearly better (stable manifolds).
        const float tol = 0.95f, abs = 0.005f * scale;
        bool better = best.sep == -kInf || sep > tol * best.sep + abs;
        if (kind == 0 && best.sep != -kInf) better = sep > best.sep;
        if (better) best = {sep, axis, kind, i, j};
        if (kind == 2 && sep > bestEdge.sep) bestEdge = {sep, axis, kind, i, j};
        return true;
    };
    for (int i = 0; i < 3; ++i)
        if (!test(a[i], 0, i, 0)) return false;
    for (int j = 0; j < 3; ++j)
        if (!test(b[j], 1, 0, j)) return false;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (!test(cross(a[i], b[j]), 2, i, j)) return false;
    return true;
}

// Edge-edge contact on axis e: the two edges are the ones whose supporting corners lie along
// the axis; the contact point is midway between their closest points.
void NarrowPhase::boxEdgeContact(const PosedShape& A, const PosedShape& B, const BoxAxis& e, ContactManifold& m) {
    const Vector3 hA = static_cast<const BoxShape*>(A.shape)->halfExtents();
    const Vector3 hB = static_cast<const BoxShape*>(B.shape)->halfExtents();
    const Vector3 a[3] = {A.R.col(0), A.R.col(1), A.R.col(2)};
    const Vector3 b[3] = {B.R.col(0), B.R.col(1), B.R.col(2)};
    const Vector3& L = e.axis;
    Vector3 pa = A.p, pb = B.p;
    for (int k = 0; k < 3; ++k) {
        if (k != e.i) pa += a[k] * (dot(a[k], L) > 0 ? hA[k] : -hA[k]);
        if (k != e.j) pb += b[k] * (dot(b[k], L) < 0 ? hB[k] : -hB[k]);
    }
    Vector3 c1, c2;
    closestSegmentSegment(pa - a[e.i] * hA[e.i], pa + a[e.i] * hA[e.i], pb - b[e.j] * hB[e.j], pb + b[e.j] * hB[e.j], c1, c2);
    m.add((c1 + c2) * 0.5f, -L, -e.sep);
}

// Face contact: the face of the other box most anti-parallel to the reference face normal (the
// incident face) is clipped against the four side planes of the reference face; the clipped
// corners below the reference plane (within the margin) are the contact points. False if none.
bool NarrowPhase::boxFaceContact(const PosedShape& A, const PosedShape& B, const BoxAxis& best, std::vector<ContactPoint>& pts) {
    const Vector3 hA = static_cast<const BoxShape*>(A.shape)->halfExtents();
    const Vector3 hB = static_cast<const BoxShape*>(B.shape)->halfExtents();
    const Vector3 a[3] = {A.R.col(0), A.R.col(1), A.R.col(2)};
    const Vector3 b[3] = {B.R.col(0), B.R.col(1), B.R.col(2)};
    const bool refIsA = best.kind == 0;
    const PosedShape& ref = refIsA ? A : B;
    const Vector3* ra = refIsA ? a : b;
    const Vector3* ia = refIsA ? b : a;
    const Vector3 rh = refIsA ? hA : hB, ih = refIsA ? hB : hA;
    const PosedShape& inc = refIsA ? B : A;
    const int r = refIsA ? best.i : best.j;
    const Vector3 n = refIsA ? best.axis : -best.axis; // reference face normal, towards the incident box

    int k = 0;
    float kd = 0;
    for (int c = 0; c < 3; ++c) {
        float d = std::fabs(dot(ia[c], n));
        if (d > kd) { kd = d; k = c; }
    }
    float s = dot(ia[k], n) > 0 ? -1.0f : 1.0f; // incident face most anti-parallel to n
    Vector3 fc = inc.p + ia[k] * (s * ih[k]);
    int k1 = (k + 1) % 3, k2 = (k + 2) % 3;
    Vector3 e1 = ia[k1] * ih[k1], e2 = ia[k2] * ih[k2];
    NarrowScratch& S = scratch(); // the worker's scratch: no allocation per pair
    std::vector<Vector3>&poly = S.poly, &clipped = S.clipped;
    poly.assign({fc + e1 + e2, fc - e1 + e2, fc - e1 - e2, fc + e1 - e2});

    int u = (r + 1) % 3, v = (r + 2) % 3;
    clipPolygon(poly, ra[u], dot(ra[u], ref.p) + rh[u], clipped);
    poly.swap(clipped);
    clipPolygon(poly, -ra[u], -dot(ra[u], ref.p) + rh[u], clipped);
    poly.swap(clipped);
    clipPolygon(poly, ra[v], dot(ra[v], ref.p) + rh[v], clipped);
    poly.swap(clipped);
    clipPolygon(poly, -ra[v], -dot(ra[v], ref.p) + rh[v], clipped);
    poly.swap(clipped);

    Vector3 refCenter = ref.p + n * rh[r];
    const Vector3 nBA = refIsA ? -n : n;
    for (const Vector3& x : poly) {
        float sep = dot(n, x - refCenter);
        if (sep <= margin) pts.push_back({x - n * (0.5f * sep), nBA, -sep});
    }
    return !pts.empty();
}

// ---------------------------------------------------------------------------
// Generic convex: GJK + EPA + perturbation manifold
// ---------------------------------------------------------------------------
// Contact manifold from the supporting faces of both shapes along the contact normal n (from B to A),
// in the spirit of Jolt's ManifoldBetweenTwoFaces: the incident feature is clipped against the side
// planes of the reference face and only points below the reference plane are kept.
// The incident feature (a face, an edge or a vertex) is cut by the side planes of the reference
// face, one edge of the reference polygon at a time: what sticks out sideways cannot touch.
void NarrowPhase::clipIncidentFace(const std::vector<Vector3>& ref, const Vector3& refNormal, std::vector<Vector3>& inc,
                                   std::vector<Vector3>& clipped) {
    Vector3 c(0.0f);
    for (const Vector3& r : ref) c += r;
    c /= float(ref.size());
    for (size_t i = 0; i < ref.size() && !inc.empty(); ++i) {
        const Vector3& r0 = ref[i];
        const Vector3& r1 = ref[(i + 1) % ref.size()];
        Vector3 side = normalize(cross(r1 - r0, refNormal));
        if (dot(side, c - r0) > 0) side = -side; // outward side plane
        float d = dot(side, r0);
        if (inc.size() >= 3) {
            clipPolygon(inc, side, d, clipped);
            inc.swap(clipped);
        } else if (inc.size() == 2) {
            float d0 = dot(side, inc[0]) - d, d1 = dot(side, inc[1]) - d;
            if (d0 > 0 && d1 > 0) inc.clear();
            else if (d0 > 0) inc[0] = inc[0] + (inc[1] - inc[0]) * (d0 / (d0 - d1));
            else if (d1 > 0) inc[1] = inc[1] + (inc[0] - inc[1]) * (d1 / (d1 - d0));
        } else if (dot(side, inc[0]) > d) {
            inc.clear();
        }
    }
}

bool NarrowPhase::faceManifold(const PosedShape& A, const PosedShape& B, const Vector3& n, std::vector<ContactPoint>& pts) {
    NarrowScratch& S = scratch(); // the worker's scratch: no allocation per pair
    std::vector<Vector3>&fa = S.fa, &fb = S.fb, &inc = S.inc, &clipped = S.cut;
    fa.clear();
    fb.clear();
    A.feature(-n, fa); // A's face towards B
    B.feature(n, fb);  // B's face towards A
    auto faceNormal = [](const std::vector<Vector3>& f) {
        return f.size() >= 3 ? normalize(cross(f[1] - f[0], f[2] - f[0])) : Vector3(0.0f);
    };
    const float cosMax = 0.9f; // reference face must be within ~25 degrees of the contact normal

    if (fa.size() == 2 && fb.size() == 2) { // edge - edge
        Vector3 c1, c2;
        closestSegmentSegment(fa[0], fa[1], fb[0], fb[1], c1, c2);
        float depth = dot(c2 - c1, n);
        if (depth < -margin) return false;
        pts.push_back({(c1 + c2) * 0.5f, n, depth});
        return true;
    }
    float alignA = fa.size() >= 3 ? std::fabs(dot(faceNormal(fa), n)) : -1.0f;
    float alignB = fb.size() >= 3 ? std::fabs(dot(faceNormal(fb), n)) : -1.0f;
    bool refIsB = alignB >= alignA;
    float align = refIsB ? alignB : alignA;
    if (align < cosMax) return false; // no face faces the contact: vertex/edge contact

    const std::vector<Vector3>& ref = refIsB ? fb : fa;
    inc = refIsB ? fa : fb;
    Vector3 nr = faceNormal(ref);
    if (dot(nr, refIsB ? n : -n) < 0) nr = -nr; // reference normal pointing towards the other shape
    clipIncidentFace(ref, nr, inc, clipped);
    // Depth of every clipped point along the contact normal n, from the other shape's supporting
    // plane along n (not along the reference face's own normal, which may differ by up to ~25 deg):
    // then the deepest point is exactly the EPA penetration, and all depths agree with n.
    const float planeB = dot(n, B.support(n)), planeA = dot(n, A.support(-n));
    bool any = false;
    for (const Vector3& x : inc) {
        const float depth = refIsB ? planeB - dot(n, x) : dot(n, x) - planeA; // incident points: A's (ref B) or B's (ref A)
        if (depth < -margin) continue;
        pts.push_back({x + n * (refIsB ? 0.5f * depth : -0.5f * depth), n, depth}); // midway between the surfaces
        any = true;
    }
    return any;
}

bool NarrowPhase::convexConvex(const PosedShape& A, const PosedShape& B, ContactManifold& m) {
    GjkResult g = gjk(A, B, margin); // pairs farther apart than the margin stop early
    Vector3 n;
    ContactPoint single;
    if (!g.intersect) {
        // Separated but within the margin: speculative contacts.
        if (g.distance >= margin || g.distance < 1e-9f) return false;
        n = (g.pointA - g.pointB) / g.distance;
        single = {(g.pointA + g.pointB) * 0.5f, n, -g.distance};
    } else {
        PenetrationResult pr = epa(A, B, g);
        if (!pr.valid) return false;
        n = pr.normal;
        single = {(pr.pointA + pr.pointB) * 0.5f, n, pr.depth};
    }

    // 1) Boundary simplices: clip the supporting faces (exact multi-point manifold).
    std::vector<ContactPoint>& pts = scratch().hullPoints; // the worker's scratch
    pts.clear();
    if (faceManifold(A, B, n, pts)) {
        reduceManifold(pts, 4);
        for (auto& p : pts) m.points.push_back(p);
        return true;
    }
    pts.push_back(single);
    if (!g.intersect) {
        m.points.push_back(single);
        return true;
    }
    // 2) Fallback for vertex contacts: perturbation method (Bullet) around the EPA normal.
    perturbationManifold(A, B, n, pts);
    reduceManifold(pts, 4);
    for (auto& p : pts) m.points.push_back(p);
    return true;
}

// Perturbation manifold (Bullet's btPerturbedContactResult): the smaller shape is tilted by a
// small angle around four directions perpendicular to the normal and the penetration is found
// again each time; the points that come back close to the normal, un-tilted, are the corners of
// the contact patch a single GJK/EPA query cannot see.
void NarrowPhase::perturbationManifold(const PosedShape& A, const PosedShape& B, const Vector3& n, std::vector<ContactPoint>& pts) {
    bool perturbA = B.shape->type() == ShapeType::Triangle ||
                    (A.shape->type() != ShapeType::Triangle && A.shape->boundingRadius() <= B.shape->boundingRadius());
    const PosedShape& P = perturbA ? A : B;
    const float radius = std::max(P.shape->boundingRadius(), 1e-3f);
    const float angle = clampv(0.02f / radius, 0.005f, 0.2f);
    const Vector3 u = perpendicular(n), v = cross(n, u);
    const float mergeDist = 0.05f * radius;
    for (int k = 0; k < 4; ++k) {
        float t = 0.5f * kPi * k;
        Vector3 axis = u * std::cos(t) + v * std::sin(t);
        Matrix3x3 Q = Quaternion::fromAxisAngle(axis, angle).toMatrix3x3();
        PosedShape Pp = P;
        Pp.R = Q * P.R; // rotation about the shape's own centre
        PenetrationResult q;
        bool ok = perturbA ? penetration(Pp, B, q) : penetration(A, Pp, q);
        if (!ok || dot(q.normal, n) < 0.9f) continue;
        Matrix3x3 Qt = Q.transposed();
        Vector3 pa = q.pointA, pb = q.pointB;
        if (perturbA) pa = Qt * (pa - P.p) + P.p;
        else pb = Qt * (pb - P.p) + P.p;
        float depth = dot(pb - pa, n);
        if (depth < -margin) continue;
        Vector3 pos = (pa + pb) * 0.5f;
        bool dup = false;
        for (const ContactPoint& c : pts) dup |= length(c.position - pos) < mergeDist;
        if (!dup) pts.push_back({pos, n, depth});
    }
}

// ---------------------------------------------------------------------------
// Capsules: a core segment plus a radius
// ---------------------------------------------------------------------------
// Sphere against capsule: the sphere against the nearest point of the capsule's axis, as two spheres.
bool NarrowPhase::sphereCapsule(const PosedShape& A, const PosedShape& B, ContactManifold& m) {
    Vector3 a, b;
    float rb;
    capsuleAxis(B, a, b, rb);
    const float ra = sphereRadius(A);
    const Vector3 q = closestOnSegment(a, b, A.p);
    const Vector3 d = A.p - q;
    const float l = length(d), depth = ra + rb - l;
    if (depth <= -margin) return false;
    const Vector3 n = l > 1e-9f ? d / l : B.R * Vector3(1, 0, 0);
    m.add(q + n * (rb - 0.5f * depth), n, depth);
    return true;
}

// Capsule against capsule: the nearest points of the two axes; nearly parallel axes that overlap
// touch along a line, so each axis's ends are also put against the other axis (two points).
bool NarrowPhase::capsuleCapsule(const PosedShape& A, const PosedShape& B, ContactManifold& m) {
    Vector3 a0, a1, b0, b1;
    float ra, rb;
    capsuleAxis(A, a0, a1, ra);
    capsuleAxis(B, b0, b1, rb);
    auto touch = [&](const Vector3& pa, const Vector3& pb) {
        const Vector3 d = pa - pb;
        const float l = length(d), depth = ra + rb - l;
        if (depth <= -margin || l < 1e-9f) return false;
        const Vector3 n = d / l;
        m.add(pb + n * (rb - 0.5f * depth), n, depth);
        return true;
    };
    const Vector3 ua = a1 - a0, ub = b1 - b0;
    const float la = length(ua), lb = length(ub);
    const bool parallel = la > 1e-6f && lb > 1e-6f && std::fabs(dot(ua, ub)) > 0.995f * la * lb;
    if (!parallel) {
        Vector3 c1, c2;
        closestSegmentSegment(a0, a1, b0, b1, c1, c2);
        return touch(c1, c2);
    }
    bool hit = false;
    for (const Vector3& e : {a0, a1}) hit |= touch(e, closestOnSegment(b0, b1, e));
    for (const Vector3& e : {b0, b1}) hit |= touch(closestOnSegment(a0, a1, e), e);
    reduceManifold(m.points, 4);
    return hit;
}

// Capsule against a box, hull or triangle: GJK between the core segment and B gives the normal and
// the gap (EPA the depth when the segment itself is inside B); the capsule's surface is the radius
// further out. Lying on a face it touches at both ends of its side (capsuleOnFace), else at one point.
bool NarrowPhase::capsuleConvex(const PosedShape& A, const PosedShape& B, ContactManifold& m) {
    const auto* cap = static_cast<const CapsuleShape*>(A.shape);
    const float r = cap->radius();
    const SegmentShape segment(cap->halfHeight());
    const PosedShape core{&segment, A.R, A.p};
    GjkResult g = gjk(core, B, r + margin);
    Vector3 n, deepest;
    float gap;
    if (!g.intersect) {
        if (g.distance >= r + margin || g.distance < 1e-9f) return false;
        n = (g.pointA - g.pointB) / g.distance;
        gap = g.distance;
        deepest = g.pointA;
    } else {
        PenetrationResult pr = epa(core, B, g);
        if (!pr.valid) return false;
        n = pr.normal;
        gap = -pr.depth;
        deepest = pr.pointA;
    }
    if (capsuleOnFace(A, B, n, m)) return true;
    m.add(deepest - n * (0.5f * (r + gap)), n, r - gap); // midway between the two surfaces
    return true;
}

// Both ends of the capsule's axis, clipped to B's face along n and measured from B's supporting
// plane: the two contacts of a capsule lying on a table (one would let it rock about that point).
bool NarrowPhase::capsuleOnFace(const PosedShape& A, const PosedShape& B, const Vector3& n, ContactManifold& m) {
    Vector3 a0, a1;
    float r;
    capsuleAxis(A, a0, a1, r);
    const Vector3 axis = a1 - a0;
    const float l = length(axis);
    if (l < 1e-6f || std::fabs(dot(axis, n)) > 0.17f * l) return false; // standing on a cap: one point
    NarrowScratch& S = scratch();
    std::vector<Vector3>&face = S.fb, &inc = S.inc, &clipped = S.cut;
    face.clear();
    B.feature(n, face); // B's face towards the capsule
    if (face.size() < 3) return false;
    Vector3 nr = normalize(cross(face[1] - face[0], face[2] - face[0]));
    if (dot(nr, n) < 0) nr = -nr;
    if (dot(nr, n) < 0.9f) return false; // no face faces the capsule
    inc.assign({a0, a1});
    clipIncidentFace(face, nr, inc, clipped);
    const float plane = dot(n, B.support(n));
    bool any = false;
    for (const Vector3& x : inc) {
        const float gap = dot(n, x) - plane; // the axis point's height above B's surface
        if (r - gap <= -margin) continue;
        m.add(x - n * (0.5f * (r + gap)), n, r - gap);
        any = true;
    }
    return any;
}

} // namespace rf
