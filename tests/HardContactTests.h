#pragma once
// Hard rigid-body contact cases - the ones engines (Bullet, PhysX, Jolt, Box2D) are known to get
// wrong. Every check has an answer computed independently of the engine: the crossing point of two
// edges, the octagon two rotated faces share, the minimum translation vector from a brute-force
// 15-axis separating-axis search, a floor of triangles that must feel like one flat plane, a Jenga
// tower that must stand still. Included by tests.cpp (uses its CHECK and run).

#include "rigid/NarrowPhase.h"
#include "rigid/RigidWorld.h"
#include "spatial/BVH.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace hard_contacts {

using namespace rf;

inline Matrix3x3 rotation(const Vector3& axis, float angle) { return Quaternion::fromAxisAngle(normalize(axis), angle).toMatrix3x3(); }

// The two narrow-phase paths a box pair can take: SAT + clipping (BoxShape) and GJK/EPA + face
// clipping (the same box as a general convex hull). Both must give the same contacts.
struct TwoWays {
    const char* name;
    bool hull;
};
const TwoWays kWays[2] = {{"SAT", false}, {"GJK/EPA", true}};

inline bool collideBoxes(bool hull, const Vector3& ha, const Matrix3x3& Ra, const Vector3& pa, const Vector3& hb, const Matrix3x3& Rb,
                         const Vector3& pb, ContactManifold& m) {
    if (!hull) {
        BoxShape A(ha), B(hb);
        return NarrowPhase::boxBox({&A, Ra, pa}, {&B, Rb, pb}, m);
    }
    ConvexHullShape A(primitives::box(ha)), B(primitives::box(hb));
    return NarrowPhase::convexConvex({&A, Ra, pa}, {&B, Rb, pb}, m);
}

// 1. Crossed beams touching edge on edge (both turned 45 deg about their long axes, like two
//    diamonds): one contact, exactly at the crossing, normal along the cross product of the edges
//    (here vertical), depth as set - even when the edges are nearly parallel.
inline void edgeOnEdge() {
    const Vector3 h(1.0f, 0.05f, 0.05f);
    const float ridge = 0.05f * std::sqrt(2.0f), depth = 0.01f;
    const Matrix3x3 Ra = rotation({1, 0, 0}, 0.25f * kPi);
    float worstPos = 0, worstDepth = 0, worstNormal = 1;
    for (const TwoWays& way : kWays)
        for (float deg : {90.0f, 60.0f, 30.0f, 10.0f}) {
            const Matrix3x3 Rb = rotation({0, 1, 0}, deg * kPi / 180) * rotation({1, 0, 0}, 0.25f * kPi);
            ContactManifold m;
            const bool hit = collideBoxes(way.hull, h, Ra, Vector3(0.0f), h, Rb, {0, 2 * ridge - depth, 0}, m);
            CHECK(hit && !m.points.empty(), "%s: crossed beams at %.0f deg not in contact", way.name, deg);
            if (m.points.empty()) continue;
            const ContactPoint* deepest = &m.points[0];
            for (const ContactPoint& c : m.points)
                if (c.depth > deepest->depth) deepest = &c;
            const Vector3 expected(0.0f, ridge - 0.5f * depth, 0.0f);
            worstPos = std::max(worstPos, length(deepest->position - expected));
            worstDepth = std::max(worstDepth, std::fabs(deepest->depth - depth));
            worstNormal = std::min(worstNormal, -deepest->normal.y); // B above A: the normal points down
            CHECK(length(deepest->position - expected) < 3e-3f && std::fabs(deepest->depth - depth) < 1e-3f && -deepest->normal.y > 0.999f,
                  "%s, beams crossed at %.0f deg: contact (%f %f %f) depth %f normal (%f %f %f)", way.name, deg, deepest->position.x,
                  deepest->position.y, deepest->position.z, deepest->depth, deepest->normal.x, deepest->normal.y, deepest->normal.z);
        }
    std::printf("  edge on edge (beams crossed at 90/60/30/10 deg, SAT and GJK/EPA): point off by %.1e m, depth off by %.1e m, "
                "normal . (0,-1,0) >= %.5f\n",
                worstPos, worstDepth, worstNormal);
}

// 2. A cube resting on a cube turned 10 and 45 deg: the shared area is an octagon. At most 4
//    points are kept - they must lie inside both faces and span most of the octagon (the best 4 of
//    a regular octagon cover 71 % of it), or the top cube rocks.
inline void rotatedFaceOnFace() {
    const Vector3 h(0.5f);
    const float depth = 0.01f;
    for (const TwoWays& way : kWays)
        for (float deg : {10.0f, 45.0f}) {
            const float a = deg * kPi / 180;
            ContactManifold m;
            collideBoxes(way.hull, h, Matrix3x3(), Vector3(0.0f), h, rotation({0, 1, 0}, a), {0, 1 - depth, 0}, m);
            bool inside = !m.points.empty();
            float depthErr = 0;
            std::vector<std::pair<float, float>> xz;
            for (const ContactPoint& c : m.points) {
                const float x = c.position.x, z = c.position.z;
                const float u = std::cos(a) * x - std::sin(a) * z, v = std::sin(a) * x + std::cos(a) * z; // in the top cube's frame
                inside &= std::fabs(x) < 0.5f + 1e-3f && std::fabs(z) < 0.5f + 1e-3f && std::fabs(u) < 0.5f + 1e-3f && std::fabs(v) < 0.5f + 1e-3f;
                inside &= -c.normal.y > 0.999f;
                depthErr = std::max(depthErr, std::fabs(c.depth - depth));
                xz.push_back({x, z});
            }
            // Area of the support polygon (points sorted by angle around their centroid).
            float cx = 0, cz = 0;
            for (auto& p : xz) { cx += p.first; cz += p.second; }
            cx /= std::max<size_t>(1, xz.size());
            cz /= std::max<size_t>(1, xz.size());
            std::sort(xz.begin(), xz.end(), [&](auto& p, auto& q) { return std::atan2(p.second - cz, p.first - cx) < std::atan2(q.second - cz, q.first - cx); });
            float area = 0;
            for (size_t i = 0; i < xz.size(); ++i) {
                const auto& p = xz[i];
                const auto& q = xz[(i + 1) % xz.size()];
                area += 0.5f * (p.first * q.second - q.first * p.second);
            }
            area = std::fabs(area);
            // At 45 deg the shared region is a regular octagon: area 8 a^2 tan(pi/8) = 0.828 m^2
            // (apothem a = 0.5); the best 4 of its corners span 0.586 m^2.
            std::printf("  cube on cube turned %2.0f deg (%s): %zu points, all inside both faces: %s, depth off by %.1e, support area %.3f m^2%s\n",
                        deg, way.name, m.points.size(), inside ? "yes" : "NO", depthErr, area,
                        deg == 45.0f ? " (octagon 0.828, best 4 points 0.586)" : "");
            CHECK(inside && depthErr < 1e-3f, "%s at %.0f deg: points outside the shared face or wrong depth", way.name, deg);
            CHECK(area > (deg == 45.0f ? 0.55f : 0.75f), "%s at %.0f deg: support area only %f", way.name, deg, area);
        }
}

// 3. Scale 1 : 1000 - a 2 cm box on a 20 m block: the four contacts at the small box's corners,
//    to a tenth of a millimetre.
inline void tinyOnHuge() {
    const Vector3 big(10.0f), small(0.01f);
    const float depth = 0.001f, a = 0.35f;
    for (const TwoWays& way : kWays) {
        ContactManifold m;
        collideBoxes(way.hull, small, rotation({0, 1, 0}, a), {0.3f, 0.01f - depth, -0.2f}, big, Matrix3x3(), {0, -10, 0}, m);
        float worst = 0;
        for (int cx = -1; cx <= 1; cx += 2)
            for (int cz = -1; cz <= 1; cz += 2) {
                const Vector3 corner = Vector3(0.3f, 0.0f, -0.2f) + rotation({0, 1, 0}, a) * Vector3(0.01f * cx, 0.0f, 0.01f * cz);
                float best = 1e9f;
                for (const ContactPoint& c : m.points) best = std::min(best, length(Vector3(c.position.x - corner.x, 0, c.position.z - corner.z)));
                worst = std::max(worst, best);
            }
        std::printf("  2 cm box on a 20 m block (%s): %zu points, farthest corner %.1e m from its contact\n", way.name, m.points.size(), worst);
        CHECK(m.points.size() == 4 && worst < 1e-4f, "%s: small box on huge box: %zu points, corner error %f", way.name, m.points.size(), worst);
    }
}

// Minimum translation vector of two overlapping boxes by brute force over the 15 SAT axes:
// the reference the narrow phase must reproduce. Returns the depth, `normal` from B to A;
// `secondBest` is the next smallest overlap on a different axis (ambiguity check).
inline float satReference(const Vector3& ha, const Matrix3x3& Ra, const Vector3& pa, const Vector3& hb, const Matrix3x3& Rb,
                          const Vector3& pb, Vector3& normal, float& secondBest) {
    std::vector<Vector3> axes;
    for (int i = 0; i < 3; ++i) {
        axes.push_back(Ra.col(i));
        axes.push_back(Rb.col(i));
    }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            const Vector3 c = cross(Ra.col(i), Rb.col(j));
            if (length(c) > 1e-4f) axes.push_back(normalize(c));
        }
    float best = 1e9f;
    secondBest = 1e9f;
    for (const Vector3& L : axes) {
        const float ra = ha.x * std::fabs(dot(Ra.col(0), L)) + ha.y * std::fabs(dot(Ra.col(1), L)) + ha.z * std::fabs(dot(Ra.col(2), L));
        const float rb = hb.x * std::fabs(dot(Rb.col(0), L)) + hb.y * std::fabs(dot(Rb.col(1), L)) + hb.z * std::fabs(dot(Rb.col(2), L));
        const float d = dot(pa - pb, L);
        const float overlap = ra + rb - std::fabs(d);
        if (overlap < best - 1e-6f) {
            if (length(cross(L, normal)) > 1e-3f || best == 1e9f) secondBest = best;
            best = overlap;
            normal = d >= 0 ? L : -L;
        } else if (overlap < secondBest && length(cross(L, normal)) > 1e-3f) {
            secondBest = overlap;
        }
    }
    return best;
}

// 4. Deep penetration of randomly posed boxes: the GJK/EPA path must find exactly the minimum
//    translation (depth and direction) of the brute-force reference. SAT deliberately prefers a
//    face axis unless an edge axis is clearly better (5 % + a small absolute tolerance, as ODE and
//    Bullet's box-box do) - it gives a stable 4-point face manifold instead of a single edge point;
//    its depth may exceed the minimum by that tolerance, never more.
inline void deepPenetration() {
    std::mt19937 rng(7);
    std::uniform_real_distribution<float> U(-1, 1), S(0.15f, 0.8f);
    int tested = 0, wrong[2] = {0, 0};
    float worstDepth[2] = {0, 0};
    for (int t = 0; t < 400; ++t) {
        const Vector3 ha(S(rng), S(rng), S(rng)), hb(S(rng), S(rng), S(rng));
        const Matrix3x3 Ra = rotation({U(rng), U(rng), U(rng) + 0.01f}, 3 * U(rng)), Rb = rotation({U(rng), U(rng) + 0.01f, U(rng)}, 3 * U(rng));
        const Vector3 pb = Vector3(U(rng), U(rng), U(rng)) * (0.5f * std::min(maxComp(ha), maxComp(hb)));
        Vector3 n;
        float second;
        const float ref = satReference(ha, Ra, Vector3(0.0f), hb, Rb, pb, n, second);
        if (ref <= 0 || second - ref < 0.02f * ref + 1e-3f) continue; // not overlapping / two axes tie: no unique answer
        ++tested;
        for (int w = 0; w < 2; ++w) {
            ContactManifold m;
            collideBoxes(kWays[w].hull, ha, Ra, Vector3(0.0f), hb, Rb, pb, m);
            float deepest = -1;
            Vector3 dn(0.0f);
            for (const ContactPoint& c : m.points)
                if (c.depth > deepest) { deepest = c.depth; dn = c.normal; }
            worstDepth[w] = std::max(worstDepth[w], std::fabs(deepest - ref));
            const float scale = std::min(minComp(ha), minComp(hb));
            const bool exact = std::fabs(deepest - ref) < 2e-3f + 1e-3f * ref && dot(dn, n) > 0.995f;
            const bool faceBias = !kWays[w].hull && deepest >= ref - 2e-3f && deepest <= ref / 0.95f + 0.005f * scale + 2e-3f;
            if (m.points.empty() || !(exact || faceBias)) ++wrong[w];
        }
    }
    std::printf("  deep penetration, %d random box pairs with a unique answer: wrong SAT %d (depth above the minimum <= %.1e, face bias), "
                "wrong EPA %d (depth off <= %.1e)\n",
                tested, wrong[0], worstDepth[0], wrong[1], worstDepth[1]);
    CHECK(tested > 100 && wrong[0] == 0 && wrong[1] == 0, "minimum translation wrong: SAT %d, EPA %d of %d", wrong[0], wrong[1], tested);
}

// 5. Internal edges ("ghost collisions"): a box sliding without friction over a floor made of 512
//    triangles must feel one flat plane - no bumps at the seams between the triangles (the classic
//    problem that Bullet's internal-edge utility and Box2D's ghost vertices exist for).
inline void slidingOverSeams() {
    TriMesh floor;
    const int n = 16;
    const float size = 8.0f;
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i <= n; ++i) floor.positions.push_back({-0.5f * size + size * i / n, 0.0f, -0.5f * size + size * j / n});
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            const uint32_t a = uint32_t(i + (n + 1) * j), b = a + 1, c = a + uint32_t(n + 1), d = c + 1;
            floor.triangles.push_back({a, c, b});
            floor.triangles.push_back({b, c, d});
        }
    MeshBVH bvh;
    bvh.build(floor);
    RigidWorld w;
    w.setDomain(AABB({-6, -2, -6}, {6, 4, 6}));
    w.setStaticMesh(&bvh);
    w.params.sleeping = false;
    const int box = w.addBox({-3.0f, 0.1f, -1.0f}, Vector3(0.1f), Quaternion(), 500, Vector3(1));
    RigidBody& b = w.bodies()[box];
    b.friction = b.staticFriction = 0;
    b.vel = {4.0f, 0.0f, 1.3f}; // diagonally across the grid: over edges and the triangles' diagonals
    float vyMax = 0, wMax = 0, yMin = 1e9f, yMax = -1e9f;
    for (int f = 0; f < 90; ++f)
        for (int k = 0; k < w.params.substeps; ++k) {
            w.step(1.0f / 60 / w.params.substeps);
            if (f < 5) continue; // settling onto the floor
            vyMax = std::max(vyMax, std::fabs(b.vel.y));
            wMax = std::max(wMax, length(b.angVel));
            yMin = std::min(yMin, b.pos.y);
            yMax = std::max(yMax, b.pos.y);
        }
    std::printf("  box sliding over 512 floor triangles: max |vy| %.4f m/s, max |w| %.4f rad/s, height %.4f..%.4f m (0.1), speed %.2f m/s\n",
                vyMax, wMax, yMin, yMax, length(b.vel));
    CHECK(vyMax < 0.05f && wMax < 0.3f && yMin > 0.095f && yMax < 0.105f, "bumps at triangle seams: vy %f w %f y %f..%f", vyMax, wMax, yMin, yMax);
    CHECK(length(b.vel) > 3.5f, "the frictionless box must keep sliding (%f)", length(b.vel));
}

// 6. A Jenga tower: 16 layers of 3 blocks, crossed layer over layer - hundreds of face and edge
//    contacts. It must stand still.
inline void jengaTower() {
    RigidWorld w;
    w.setDomain(AABB({-3, 0, -3}, {3, 8, 3}));
    const Vector3 h(0.375f, 0.075f, 0.125f); // real Jenga proportions, 5x scale
    const float gap = 0.002f;
    std::vector<int> ids;
    std::vector<Vector3> start;
    const int layers = 16;
    for (int l = 0; l < layers; ++l)
        for (int k = -1; k <= 1; ++k) {
            const bool alongX = l % 2 == 0;
            const float off = k * (2 * h.z + gap);
            const Vector3 p = alongX ? Vector3(0.0f, h.y + 2 * h.y * l, off) : Vector3(off, h.y + 2 * h.y * l, 0.0f);
            const Quaternion q = alongX ? Quaternion() : Quaternion::fromAxisAngle({0, 1, 0}, 0.5f * kPi);
            ids.push_back(w.addBox(p, h, q, 500, Vector3(1)));
            start.push_back(p);
        }
    float worstOverlap = 0;
    for (int f = 0; f < 60 * 5; ++f)
        for (int k = 0; k < w.params.substeps; ++k) {
            w.step(1.0f / 60 / w.params.substeps);
            if (f % 30 == 0 && k == 0) worstOverlap = std::max(worstOverlap, maxOverlap(w));
        }
    float drift = 0, speed = 0;
    for (size_t i = 0; i < ids.size(); ++i) {
        const RigidBody& b = w.bodies()[ids[i]];
        drift = std::max(drift, length(b.pos - start[i]));
        speed = std::max(speed, length(b.vel));
    }
    const float topY = w.bodies()[ids.back()].pos.y, expected = start.back().y;
    std::printf("  Jenga tower %d x 3 blocks after 5 s: max drift %.4f m, top at %.4f (start %.4f), max speed %.4f m/s, max overlap %.4f m, "
                "%zu contacts\n",
                layers, drift, topY, expected, speed, worstOverlap, w.contactCount());
    CHECK(drift < 0.01f && speed < 0.02f && worstOverlap < 0.005f, "Jenga tower moved (drift %f speed %f overlap %f)", drift, speed, worstOverlap);
}

// 7. A diamond beam dropped across a ridge (a box turned 45 deg: edge up) - an edge hitting an
//    edge at speed. No passing through, and no energy from nowhere (the total never rises).
inline void edgeDropOnRidge() {
    RigidWorld w;
    w.setDomain(AABB({-3, 0, -3}, {3, 4, 3}));
    w.params.sleeping = false;
    w.addBox({0.0f, 0.3f, 0.0f}, Vector3(0.3f, 0.3f, 1.0f), Quaternion::fromAxisAngle({0, 0, 1}, 0.25f * kPi), 0, Vector3(1)); // ridge along z
    const Quaternion diamond = Quaternion::fromAxisAngle({1, 0, 0}, 0.25f * kPi);
    const int beam = w.addBox({0.02f, 1.3f, 0.1f}, {0.8f, 0.05f, 0.05f}, diamond, 700, Vector3(1)); // crosses the ridge
    const RigidBody& b = w.bodies()[beam];
    auto energy = [&] {
        const Vector3 L = b.angVel; // rotational part with the world inverse inertia
        const Matrix3x3 invI = b.invInertiaWorld;
        const Vector3 Iw = invI.inverse() * L;
        return b.mass * 9.81f * b.pos.y + 0.5f * b.mass * length2(b.vel) + 0.5f * dot(L, Iw);
    };
    const float e0 = energy();
    float eMax = e0, worst = 0;
    for (int f = 0; f < 180; ++f)
        for (int k = 0; k < w.params.substeps; ++k) {
            w.step(1.0f / 60 / w.params.substeps);
            eMax = std::max(eMax, energy());
            worst = std::max(worst, maxOverlap(w));
        }
    std::printf("  diamond beam dropped across a ridge: max overlap %.4f m, energy max/start %.4f, beam comes to y %.3f m\n", worst, eMax / e0,
                b.pos.y);
    CHECK(worst < 0.01f, "the beam went into the ridge (%f)", worst);
    CHECK(eMax < 1.01f * e0, "energy created in the edge impact (%f of the start)", eMax / e0);
}

} // namespace hard_contacts

void testHardContacts() {
    hard_contacts::edgeOnEdge();
    hard_contacts::rotatedFaceOnFace();
    hard_contacts::tinyOnHuge();
    hard_contacts::deepPenetration();
}

void testHardContactDynamics() {
    hard_contacts::slidingOverSeams();
    hard_contacts::jengaTower();
    hard_contacts::edgeDropOnRidge();
}
