// Non-convex rigid bodies in the scenes the physics engines show them in - Havok's and Bullet's
// chains of rings, a ring tossed over a peg, things dropped into a cup, a ring lying on the floor -
// and the convex decomposition behind them held to the geometry: a hole stays a hole, a hollow
// stays hollow. Every number is checked against the shape: the distance two interlocked rings hang
// apart, the height of the floor a body rests on. The standard scenes of samples/NonConvexScenes.cpp
// (32-36) are held each to its own number: a law of rolling, a height of a stack, the surface.
#include "TestRunner.h"
#include "Tests.h"

#include "rigid/ConvexDecomposition.h"
#include "rigid/GjkEpa.h"

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

// The race of sample 32: a ball, a solid cylinder (a hull of 48 sides) and a wheel of boxes (a rim
// of 48 flat segments) let go side by side on a 15 degree slope. The ball is round: with the
// rolling resistance c (a moment c N rho, rho the lever of RigidParams::rollingResistance, the
// bounding radius) it rolls with v^2 = 2 s g (sin t - c cos t rho / R) / (1 + I / (m R^2)). A 48-gon
// rolls over its corners: between two corners it gains m g L sin t (L the side) less the rolling
// resistance, at each corner it turns about the next one keeping its angular momentum about it,
//     w' = w (k + cos a) / (k + 1),   k = I / (m R^2), a = 2 pi / 48,
// and loses 1 - (w'/w)^2 of its energy - 1.1 % for the cylinder: it runs a third slower than a
// round one would after 2.6 m. Each body is held to its own law, and they must finish in the order
// of their k: ball (0.4), cylinder (0.5), wheel (0.73). The rotational lock of a resting face held
// the wheel standing before: its centre of mass 8 cm past the 4 cm segment it stood on.
void testRollingRace() {
    Simulation sim;
    loadSample(sim, Preset::RigidRace);
    const RigidWorld& w = sim.rigid;
    const float theta = 15.0f * kPi / 180.0f, g = 9.81f, c = w.params.rollingResistance, R = 0.3f;
    const Vector3 down(std::cos(theta), -std::sin(theta), 0);
    std::vector<Vector3> start;
    for (const RigidBody& b : w.bodies()) start.push_back(b.pos);
    for (int f = 0; f < 120; ++f) sim.stepFrame(); // 2 s
    const char* names[] = {"ball", "cylinder of 48 sides", "wheel of 48 segments"};
    float travelled[3] = {};
    for (int i = 0; i < 3; ++i) {
        const RigidBody& b = w.bodies()[size_t(i + 1)]; // body 0 is the slope
        const Matrix3x3 invI = b.rotation() * Matrix3x3::diag(b.invInertiaLocal) * b.rotation().transposed();
        const float corner = i == 2 ? R / std::cos(kPi / 48) : R; // the wheel's segments touch R out, its corners are farther
        const float k = 1 / (invI.m[2][2] * b.mass * corner * corner), rho = b.boundingRadius();
        const float s = dot(b.pos - start[size_t(i + 1)], down), v = dot(b.vel, down);
        const float slope = std::sin(theta) - c * std::cos(theta) * rho / corner;
        float expected;
        if (i == 0) {
            expected = std::sqrt(2 * g * s * slope / (1 + k));
        } else { // corner by corner: the energy per unit mass, E = (1 + k) v^2 / 2
            const float a = 2 * kPi / 48, side = 2 * corner * std::sin(kPi / 48), keep = sqr((k + std::cos(a)) / (k + 1));
            double E = 0;
            for (float x = side; x <= s; x += side) E = (E + g * side * slope) * keep;
            expected = float(std::sqrt(2 * E / (1 + k)));
        }
        const float slip = std::fabs(v - std::fabs(b.angVel.z) * R) / v;
        std::printf("  %-21s k %.3f: 2 s down the slope %.3f m at %.3f m/s, its law %.3f m/s (%+.1f %%), slip %.1f %%\n", names[i], k, s, v, expected,
                    100 * (v / expected - 1), 100 * slip);
        CHECK(std::fabs(v / expected - 1) < (i == 0 ? 0.03f : 0.05f), "%s: %.3f m/s, its law gives %.3f", names[i], v, expected);
        CHECK(slip < 0.03f, "%s slides: slip %.1f %%", names[i], 100 * slip);
        travelled[i] = s;
    }
    CHECK(travelled[0] > travelled[1] && travelled[1] > travelled[2], "the order is not ball, cylinder, wheel: %.3f %.3f %.3f m", travelled[0], travelled[1],
          travelled[2]);
}

// Sample 33: 60 boxes, balls and capsules poured into a bowl (a half shell 1.2 m across, 4 cm
// thick, 48 convex parts). Every one stays inside, none sinks into the shell, the pile falls
// asleep, and the energy never exceeds what the pour started with.
void testBowlHoldsPour() {
    Simulation sim;
    loadSample(sim, Preset::RigidBowl);
    const RigidWorld& w = sim.rigid;
    const Energy e0 = energyOf(w);
    double gain = 0;
    for (int f = 0; f < 360; ++f) {
        sim.stepFrame();
        const Energy e = energyOf(w);
        gain = std::max(gain, e.kinetic + e.potential - e0.kinetic - e0.potential);
    }
    int outside = 0, awake = 0;
    float lowest = kInf;
    for (size_t i = 1; i < w.bodies().size(); ++i) {
        const RigidBody& b = w.bodies()[i];
        outside += std::hypot(b.pos.x, b.pos.z) > 0.6f || b.pos.y > 0.6f;
        awake += !b.sleeping;
        lowest = std::min(lowest, b.posed().support(Vector3(0, -1, 0)).y);
    }
    std::printf("  bowl, 6 s: %d of 60 outside, %d awake, the lowest point %.1f mm above the floor (the inner bottom 40 mm), energy gain %.4f J\n",
                outside, awake, 1000 * lowest, gain);
    CHECK(outside == 0, "%d bodies left the bowl", outside);
    CHECK(awake == 0, "%d bodies still awake after 6 s", awake);
    CHECK(lowest > 0.04f - 0.004f, "a body sank into the bowl's bottom: %.1f mm", 1000 * lowest);
    CHECK(gain < 1e-3, "the pour gained %.4f J", gain);
}

// Sample 34: eight tables, each a top on four legs (five boxes), dropped 5 cm onto one another,
// shifted up to 1 cm and turned up to 1.6 degrees, every leg over the top below. They stack as high
// as eight tables are, 3.52 m, stay where they landed and sleep. Of such stacks (other shifts and
// turns of the same size) 12 of 12 stand with 5 tables and with 8, 9 of 12 with 6: there a table
// rocks on its four legs within the slop and walks off the one below (an open question). With a leg
// past the edge of the top below (shifts of 3 cm, turns of 5 degrees) a table stands on three legs,
// its centre of mass on their diagonal, and the stack topples, as it should.
void testStackOfTables() {
    Simulation sim;
    loadSample(sim, Preset::RigidTables);
    const RigidWorld& w = sim.rigid;
    std::vector<Vector3> start;
    for (const RigidBody& b : w.bodies()) start.push_back(b.pos);
    for (int f = 0; f < 240; ++f) sim.stepFrame();
    float drift = 0;
    int awake = 0;
    for (size_t i = 0; i < w.bodies().size(); ++i) {
        const Vector3 d = w.bodies()[i].pos - start[i];
        drift = std::max(drift, std::hypot(d.x, d.z));
        awake += !w.bodies()[i].sleeping;
    }
    const float top = w.bodies().back().posed().support(Vector3(0, 1, 0)).y, stacked = 0.44f * float(w.bodies().size());
    std::printf("  %zu tables, 4 s: top %.4f m (stacked %.4f m), sideways at most %.1f mm, %d awake\n", w.bodies().size(), top, stacked, 1000 * drift, awake);
    CHECK(w.bodies().size() == 8, "the sample holds %zu tables, not 8", w.bodies().size());
    CHECK(std::fabs(top - stacked) < 0.002f, "the stack is %.4f m tall, %zu tables are %.2f m", top, w.bodies().size(), stacked);
    CHECK(drift < 0.005f, "a table moved %.1f mm sideways", 1000 * drift);
    CHECK(awake == 0, "%d tables awake after 4 s", awake);
}

// Sample 35: six tapered cups (8 to 12 cm across, 10 cm tall, wall 4 mm) dropped into one another.
// Exact cups would nest 20.4 mm apart (the wall's horizontal thickness over the taper, 1:5); the
// convex parts reach up to about 2 mm into the cup (the decomposition's tolerance), five times
// that along the axis, so the parts nest farther apart - how far is found here from the parts
// alone: the upper cup lowered, coaxial, until two parts meet. The stack must come to rest that
// high and sleep. Before, a contact across the seam between two parts of a wall pushed the upper
// cup down into the lower one, the stack sank 38 mm into itself, rocked and never slept (60 ms a
// frame); and a speculative contact through the wall of the cup between made the shock pass lift
// the stack in bursts.
void testNestedCups() {
    const auto cup = cupShape();
    const Vector3 up = cup->principalRotation().transposed() * Vector3(0, 1, 0);
    auto overlap = [&](float d) {
        for (const auto& a : cup->children())
            for (const auto& b : cup->children())
                if (gjk({a.shape.get(), a.R, a.t}, {b.shape.get(), b.R, b.t + up * d}, 0.0f).intersect) return true;
        return false;
    };
    float lo = 0.0f, hi = 0.1f; // overlap at lo, clear at hi
    while (hi - lo > 1e-5f) (overlap(0.5f * (lo + hi)) ? lo : hi) = 0.5f * (lo + hi);
    const float pitch = hi, expected = 0.1f + 5 * pitch;
    Simulation sim;
    loadSample(sim, Preset::RigidCups);
    const RigidWorld& w = sim.rigid;
    double kinetic = 0;
    for (int f = 0; f < 300; ++f) {
        sim.stepFrame();
        if (f >= 90) kinetic = std::max(kinetic, energyOf(w).kinetic);
    }
    float bottom = kInf, top = -kInf;
    int awake = 0;
    for (const RigidBody& b : w.bodies()) {
        bottom = std::min(bottom, b.posed().support(Vector3(0, -1, 0)).y);
        top = std::max(top, b.posed().support(Vector3(0, 1, 0)).y);
        awake += !b.sleeping;
    }
    std::printf("  six cups: the parts nest %.1f mm apart (exact cups 20.4 mm), stack %.1f mm (from the parts %.1f mm), %d awake, "
                "kinetic energy after 1.5 s at most %.2e J\n", 1000 * pitch, 1000 * (top - bottom), 1000 * expected, awake, kinetic);
    CHECK(std::fabs(top - bottom - expected) < 0.005f, "the stack is %.1f mm tall, the parts nest to %.1f mm", 1000 * (top - bottom), 1000 * expected);
    CHECK(awake == 0, "%d cups awake after 5 s", awake);
    CHECK(kinetic < 1e-4, "the nested stack moves again: %.2e J", kinetic);
}

// Sample 36: 45 teapots, bunnies and steel rings dropped on the terrain (a static mesh of 51 200
// triangles). None ends up under the surface - no point of a body lower than the surface under it
// by more than the solver's slop - and all sleep by 10 s. Before, all the contact points of a
// body with the mesh were one manifold: a bunny on two slopes had one normal for both, and 21 of
// the 45 still jittered after 24 s.
void testModelsOnTerrain() {
    Simulation sim;
    loadSample(sim, Preset::TerrainModels);
    const RigidWorld& w = sim.rigid;
    for (int f = 0; f < 600; ++f) sim.stepFrame();
    const MeshBVH* mesh = w.staticMesh();
    CHECK(mesh != nullptr, "the terrain is not the world's static mesh");
    if (!mesh) return;
    int awake = 0, under = 0;
    float deepest = 0;
    for (const RigidBody& b : w.bodies()) {
        awake += !b.sleeping;
        const Vector3 low = b.posed().support(Vector3(0, -1, 0));
        RayHit hit;
        if (!mesh->raycast(Vector3(low.x, 10.0f, low.z), Vector3(0, -1, 0), 20.0f, hit)) continue;
        const float below = (10.0f - hit.t) - low.y;
        deepest = std::max(deepest, below);
        under += below > w.params.slop;
    }
    std::printf("  45 models on the terrain, 10 s: %d awake, %d under the surface, the deepest point %.1f mm below it\n", awake, under, 1000 * deepest);
    CHECK(under == 0, "%d models sank into the terrain (deepest %.1f mm)", under, 1000 * deepest);
    CHECK(awake == 0, "%d models awake after 10 s", awake);
}
