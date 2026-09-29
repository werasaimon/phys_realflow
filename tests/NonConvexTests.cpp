// Non-convex rigid bodies in the scenes the physics engines show them in - Havok's and Bullet's
// chains of rings, a ring tossed over a peg, things dropped into a cup, a ring lying on the floor -
// and the convex decomposition behind them held to the geometry: a hole stays a hole, a hollow
// stays hollow. Every number is checked against the shape: the distance two interlocked rings hang
// apart, the height of the floor a body rests on.
#include "TestRunner.h"
#include "Tests.h"

#include "rigid/ConvexDecomposition.h"

#include <chrono>
#include <cmath>
#include <memory>

namespace {

using namespace rf;

constexpr float kDt = 1.0f / 60.0f;
constexpr float kRingR = 0.05f, kRingTube = 0.012f; // a ring of 10 cm, its tube 24 mm thick

// A cup: outer radius R, height 2 halfHeight, wall and bottom `wall` thick, open at the top (+y).
TriMesh cupMesh(float R, float halfHeight, float wall) {
    const std::vector<std::pair<float, float>> profile = {{-halfHeight, 0}, {-halfHeight, R}, {halfHeight, R},
                                                          {halfHeight, R - wall}, {-halfHeight + wall, R - wall}, {-halfHeight + wall, 0}};
    TriMesh m = primitives::revolve(profile, 48);
    for (Vector3& p : m.positions) p = {p.y, p.x, p.z}; // the axis from x to y
    m.flipWinding();                                       // (a reflection)
    m.orientOutward();
    return m;
}

std::shared_ptr<const CompoundShape> compoundOf(const TriMesh& m) {
    return std::make_shared<const CompoundShape>(convexDecomposition({m}), m);
}

// A world with a floor at y = 0 and room to hang a chain.
struct World {
    RigidWorld w;
    World() { w.setDomain(AABB({-1, 0, -1}, {1, 2, 1})); }
    void run(float seconds) {
        for (int k = 0; k < int(seconds / kDt); ++k) w.step(kDt);
    }
    float lowest(int b) { return w.bodies()[size_t(b)].posed().support(Vector3(0, -1, 0)).y; }
    // A static body put down on the floor by its lowest point (its bounds are looser than that).
    void putOnFloor(int b) { w.bodies()[size_t(b)].pos.y -= lowest(b); }
};

// Is a point of the model's space inside one of the compound's parts?
bool inside(const CompoundShape& c, const Vector3& modelPoint) {
    return c.contains(c.principalRotation().transposed() * (modelPoint - c.centerOfMass()));
}

// How far the parts reach into empty space: the largest distance from their surfaces to the model.
float reach(const TriMesh& model, const std::vector<TriMesh>& parts) {
    MeshBVH bvh;
    bvh.build(model);
    float worst = 0;
    for (const TriMesh& p : parts)
        for (const auto& t : p.triangles) {
            const Vector3 &a = p.positions[t[0]], &b = p.positions[t[1]], &c = p.positions[t[2]];
            for (int u = 0; u <= 8; ++u)
                for (int v = 0; u + v <= 8; ++v) worst = std::max(worst, bvh.signedDistance(a + (b - a) * (u / 8.0f) + (c - a) * (v / 8.0f), 1.0f));
        }
    return worst;
}

} // namespace

// The parts of a ring and of a cup keep their hole and their hollow. A convex hull of an arc of
// a ring spans the chord of the arc, into the hole; judged by volume alone (V-HACD's measure) a
// ring came out in 13 parts that reached 15 mm into the hole, and a cup with slabs across its
// hollow (89 mm into it). Judged by how deep a hull reaches into empty space (CoACD's
// collision-aware concavity), and split through its deepest point, both stay open.
void testNonConvexDecomposition() {
    const TriMesh ring = primitives::torus(kRingR, kRingTube), cup = cupMesh(0.1f, 0.08f, 0.01f);
    auto t0 = std::chrono::steady_clock::now();
    const std::vector<TriMesh> ringParts = convexDecomposition({ring});
    auto t1 = std::chrono::steady_clock::now();
    const std::vector<TriMesh> cupParts = convexDecomposition({cup});
    auto t2 = std::chrono::steady_clock::now();
    const CompoundShape ringShape(ringParts, ring), cupShape(cupParts, cup);
    int hollow = 0, blocked = 0;
    for (float y = -0.065f; y <= 0.08f; y += 0.005f) // the axis of the cup, from 5 mm over its floor to the rim
        for (float r = 0; r <= 0.045f; r += 0.015f, ++hollow) blocked += inside(cupShape, {r, y, 0});
    const bool holeOpen = !inside(ringShape, {0, 0, 0}) && !inside(ringShape, {0.02f, 0, 0}) && !inside(ringShape, {0, 0, 0.02f});
    const float ringReach = reach(ring, ringParts), cupReach = reach(cup, cupParts);
    std::printf("  ring: %zu parts in %.0f ms, reach into the hole %.1f mm; cup: %zu parts in %.0f ms, reach %.1f mm, %d of %d points of its hollow inside a part\n",
                ringParts.size(), std::chrono::duration<double, std::milli>(t1 - t0).count(), 1000 * ringReach, cupParts.size(),
                std::chrono::duration<double, std::milli>(t2 - t1).count(), 1000 * cupReach, blocked, hollow);
    CHECK(holeOpen && ringReach < 0.004f, "the ring's parts close its hole: reach %.1f mm", 1000 * ringReach);
    CHECK(blocked == 0, "%d points of the cup's hollow are inside its parts", blocked);
}

// A ring lying flat on the floor, as one convex hull and as its parts: it rests on its lowest
// ring of vertices. The manifold kept the three points spanning the largest area among all
// within the speculative margin - the wider rings above, not yet touching - and the ring sank
// until they did: 9 mm into the floor, of a 24 mm tube.
void testRingOnFloor() {
    const TriMesh ring = primitives::torus(kRingR, kRingTube);
    World hull, parts;
    const int a = hull.w.addConvex(buildConvexHull(ring.positions), {0, 0.1f, 0}, Quaternion(), 1000.0f, Vector3(1));
    const int b = parts.w.addCompound(compoundOf(ring), {0, 0.1f, 0}, Quaternion(), 1000.0f, Vector3(1));
    hull.run(3.0f), parts.run(3.0f);
    const float sinkHull = -hull.lowest(a), sinkParts = -parts.lowest(b);
    std::printf("  ring on the floor sinks: as one hull %.2f mm, as its parts %.2f mm (asleep %d %d)\n", 1000 * sinkHull, 1000 * sinkParts,
                int(hull.w.bodies()[size_t(a)].sleeping), int(parts.w.bodies()[size_t(b)].sleeping));
    CHECK(sinkHull < 0.001f && sinkParts < 0.001f, "the ring sinks into the floor: %.2f / %.2f mm", 1000 * sinkHull, 1000 * sinkParts);
}

// The energy of the bodies that move: kinetic (of their centres and of rotation about their
// principal axes) and potential in gravity.
struct Energy { double kinetic = 0, potential = 0; };
Energy energyOf(const RigidWorld& w) {
    Energy e;
    for (const RigidBody& b : w.bodies()) {
        if (b.invMass == 0) continue;
        const Vector3 wl = b.rotation().transposed() * b.angVel; // in the principal frame
        const Vector3& ii = b.invInertiaLocal;
        e.kinetic += 0.5 * b.mass * double(length2(b.vel));
        e.kinetic += 0.5 * (double(wl.x * wl.x) / ii.x + double(wl.y * wl.y) / ii.y + double(wl.z * wl.z) / ii.z);
        e.potential -= double(b.mass) * double(dot(w.params.gravity, b.pos));
    }
    return e;
}

// A chain of eight steel rings hanging from a fixed one, each at right angles to the next, let go
// 6 mm short of taut: it drops, jerks taut and swings. Two interlocked rings hang with their
// centres 2 (R - r) apart when the tubes touch: 76 mm; the parts of each ring reach up to 3 mm
// into its hole, a pair up to 6 mm, and a resting contact may sink by the solver's slop (4 mm,
// RigidParams::slop) and a little more under the weight of seven rings - a ring pulled through
// its neighbour would be 20 mm and more past. The chain may only lose energy - the rotational
// lock of a resting face held two rings touching on both sides of a tube as one and let go with
// a kick, and the chain gained up to 18 J in bursts before - and in 10 s it loses nearly all the
// drop gave it.
void testChainOfRings() {
    World s;
    const auto link = compoundOf(primitives::torus(kRingR, kRingTube));
    const Quaternion inXY = Quaternion::fromAxisAngle({1, 0, 0}, 0.5f * kPi), inYZ = Quaternion::fromAxisAngle({0, 0, 1}, 0.5f * kPi);
    const int n = 8;
    const float taut = 2 * (kRingR - kRingTube);
    for (int k = 0; k < n; ++k)
        s.w.addCompound(link, {0, 1.5f - float(k) * (taut - 0.006f), 0}, k % 2 ? inYZ : inXY, k == 0 ? 0.0f : 7800.0f, Vector3(0.7f));
    const Energy start = energyOf(s.w);
    double gained = 0;
    for (int f = 0; f < int(10.0f / kDt); ++f) {
        s.w.step(kDt);
        const Energy e = energyOf(s.w);
        gained = std::max(gained, e.kinetic + e.potential - start.kinetic - start.potential);
    }
    const Energy end = energyOf(s.w);
    const double released = start.potential - end.potential; // what the drop gave
    float nearest = kInf, farthest = 0;
    for (int k = 1; k < n; ++k) {
        const float d = length(s.w.bodies()[size_t(k)].pos - s.w.bodies()[size_t(k - 1)].pos);
        nearest = std::min(nearest, d), farthest = std::max(farthest, d);
    }
    std::printf("  chain of %d rings after 10 s: neighbours %.1f..%.1f mm apart (tubes touching: %.1f mm); energy never above the start by more than %.4f J; "
                "of the %.3f J the drop gave, %.4f J still moving\n", n, 1000 * nearest, 1000 * farthest, 1000 * taut, gained, released, end.kinetic);
    CHECK(farthest < taut + s.w.params.slop + 0.003f, "rings pulled through each other: %.1f mm apart", 1000 * farthest);
    CHECK(nearest > taut - 0.007f, "a ring hangs on another %.1f mm short of the metal", 1000 * (taut - nearest));
    CHECK(gained < 0.01, "the chain gained %.3f J", gained);
    CHECK(end.kinetic < 0.05 * released, "the chain keeps %.3f J of the %.3f J the drop gave", end.kinetic, released);
}

// A ring dropped over a peg thinner than its hole slides down the peg to the floor.
void testRingOverPeg() {
    World s;
    s.w.addBox({0, 0.15f, 0}, {0.012f, 0.15f, 0.012f}, Quaternion(), 0.0f, Vector3(0.5f));
    const int ring = s.w.addCompound(compoundOf(primitives::torus(kRingR, kRingTube)), {0, 0.45f, 0}, Quaternion(), 7800.0f, Vector3(0.7f));
    s.run(3.0f);
    const float y = s.w.bodies()[size_t(ring)].pos.y;
    std::printf("  ring over a peg: rests with its centre at %.1f mm (on the floor: %.1f mm, on the peg's top: %.0f mm)\n", 1000 * y, 1000 * kRingTube,
                1000 * (0.3f + kRingTube));
    CHECK(std::fabs(y - kRingTube) < 0.002f, "the ring did not slide down the peg: centre at %.1f mm", 1000 * y);
}

// A ball and boxes dropped into a cup (walls and bottom 1 cm) land on its bottom.
void testDropIntoCup() {
    const auto cup = compoundOf(cupMesh(0.1f, 0.08f, 0.01f));
    float worst = 0;
    for (int kind = 0; kind < 3; ++kind) {
        World s;
        const int c = s.w.addCompound(cup, {0, 0.3f, 0}, Quaternion(), 0.0f, Vector3(0.6f));
        s.putOnFloor(c);
        const int b = kind == 0 ? s.w.addSphere({0, 0.3f, 0}, 0.03f, 1000.0f, Vector3(1))
                                : s.w.addBox({0.01f, 0.3f, 0.02f}, Vector3(0.03f), Quaternion::fromAxisAngle(normalize(Vector3(1, 1, 0)), kind == 1 ? 0.3f : 0.7f),
                                             1000.0f, Vector3(1));
        s.run(4.0f);
        const float above = s.lowest(b) - 0.01f;
        std::printf("  %s dropped into a cup: rests %.1f mm above its bottom\n", kind == 0 ? "ball" : kind == 1 ? "box tilted 17 deg" : "box tilted 40 deg", 1000 * above);
        worst = std::max(worst, std::fabs(above));
    }
    CHECK(worst < 0.003f, "a body dropped into the cup rests %.1f mm off its bottom", 1000 * worst);
}

// A wheel of boxes - a hub, six spokes and a rim of 48 flat segments - let go on a 15 degree slope
// rolls down without slipping. A round wheel would reach v^2 = 2 s g sin(t) / (1 + I / (m R^2)); a
// 48-gon loses about 1 % of its energy at every corner it rolls over, so it is somewhat slower,
// never faster. The rotational lock of a resting face held it standing before: its centre of mass
// 8 cm past the 4 cm segment it stood on, the lock carried the tipping moment as if friction could.
void testWheelRollsDownSlope() {
    const int n = 48;
    const float R = 0.3f, rim = 0.03f, width = 0.08f, theta = 15.0f * kPi / 180.0f;
    auto boxAt = [](const Vector3& half, float angle, const Vector3& at) {
        TriMesh b = primitives::box(half);
        b.transform(Quaternion::fromAxisAngle({0, 0, 1}, angle).toMatrix3x3(), Vector3(1.0f), at);
        return b;
    };
    std::vector<TriMesh> parts;
    for (int k = 0; k < n; ++k) {
        const float a = 2 * kPi * float(k) / n, r = R - 0.5f * rim;
        parts.push_back(boxAt({R * std::sin(kPi / n) + 0.002f, 0.5f * rim, 0.5f * width}, a + 0.5f * kPi, {r * std::cos(a), r * std::sin(a), 0}));
    }
    for (int k = 0; k < 6; ++k) {
        const float a = 2 * kPi * float(k) / 6, r = 0.5f * (R - rim + 0.05f);
        parts.push_back(boxAt({0.5f * (R - rim - 0.05f), 0.012f, 0.012f}, a, {r * std::cos(a), r * std::sin(a), 0}));
    }
    parts.push_back(primitives::box({0.05f, 0.05f, 0.5f * width}));
    RigidWorld w;
    w.setDomain(AABB({-6, -4, -1}, {6, 4, 1}));
    w.addBox({0, 0, 0}, {5.0f, 0.1f, 0.5f}, Quaternion::fromAxisAngle({0, 0, 1}, -theta), 0.0f, Vector3(0.5f));
    const Vector3 down(std::cos(theta), -std::sin(theta), 0), up(std::sin(theta), std::cos(theta), 0);
    const Vector3 start = down * -3.5f + up * (0.1f + R + 0.001f);
    const int b = w.addCompound(std::make_shared<const CompoundShape>(parts, primitives::merge(parts)), start, Quaternion(), 1000.0f, Vector3(0.8f));
    const RigidBody& wheel = w.bodies()[size_t(b)];
    const Matrix3x3 invI = wheel.rotation() * Matrix3x3::diag(wheel.invInertiaLocal) * wheel.rotation().transposed();
    const float a = 9.81f * std::sin(theta) / (1 + 1 / (invI.m[2][2] * wheel.mass * R * R));
    for (int f = 0; f < 90; ++f) w.step(kDt);
    const float s = dot(wheel.pos - start, down), v = dot(wheel.vel, down), round = std::sqrt(2 * a * std::max(s, 0.0f));
    const float slip = std::fabs(v - std::fabs(wheel.angVel.z) * R) / std::max(v, 1e-3f);
    std::printf("  wheel of %d segments on a 15 deg slope, 1.5 s: rolled %.2f m at %.3f m/s (a round wheel: %.3f), slip %.1f %%\n", n, s, v, round, 100 * slip);
    CHECK(s > 0.8f && v > 0.75f * round && v <= round, "the wheel does not roll: %.2f m, %.3f of %.3f m/s", s, v, round);
    CHECK(slip < 0.03f, "the wheel slides: slip %.1f %%", 100 * slip);
}
