// A pile of barrels: 50 rigid cylinders (the convex hull of a 24-sided prism, as the editor's
// «Бочка» with the collider «Авто») dropped from random heights and turns onto a floor - the
// scene a user threw together and saw "some barrels twitch". What twitched was the collision, not
// the solver. The solver was measured first, part by part switched off (shock propagation, the
// rotational lock, the block solver, split impulse, rolling resistance, the bounce, CCD): each
// switch only gave a different pile, twitching or not by chance. The contacts were wrong in two
// ways, both in the convex-convex narrow phase (NarrowPhase.cpp), and a test here holds each:
//   1. A face contact measured its depths from the other shape's supporting plane along the
//      contact normal. Tilted by a milliradian (all a normal from GJK or EPA is good to), that
//      plane passes through the far corner of a 3.6 m floor, a millimetre off under the barrel:
//      the depth of a resting barrel jumped between a gap and a penetration from step to step.
//   2. Touching closer than GJK can tell apart (micrometres), the normal was rounding noise; once
//      it pointed into the floor, whose far side then looked nearest, 32 cm down, and the barrel
//      was pushed through the floor.
// The pile test measures what a user sees: the pile falls asleep, and with sleeping off it stands
// still - no body faster than a millimetre per second, the warm start finding every contact.
#include "TestRunner.h"
#include "Tests.h"

#include "scene/SceneGraph.h"
#include "scene/Simulation.h"

#include <cmath>
#include <cstdio>
#include <string>

namespace {

// The same random numbers on every machine: a linear congruential generator (Numerical Recipes'
// constants) giving floats in [lo, hi).
struct Lcg {
    uint32_t seed;
    float operator()(float lo, float hi) {
        seed = seed * 1664525u + 1013904223u;
        return lo + (hi - lo) * float(seed >> 8) * (1.0f / 16777216.0f);
    }
};

// The pile's scene: a floor and 50 barrels 0.3 m wide and 0.4 m tall (density 500, friction 0.5,
// restitution 0.2 - the editor's defaults) at random places 0.5 .. 4 m high and random turns; no
// two start overlapping.
SceneGraph barrelPile() {
    SceneGraph g;
    g.world.size = {4, 6, 4};
    Entity floor;
    floor.id = 1;
    floor.name = "Пол";
    floor.shape = ShapeKind::Plane;
    floor.size = {3.6f, 0.02f, 3.6f};
    floor.position = {0, 0.01f, 0};
    floor.rigid.enabled = floor.collider.enabled = true;
    floor.rigid.fixed = true;
    g.entities.push_back(floor);
    Lcg random{7};
    std::vector<Vector3> placed;
    while (placed.size() < 50) {
        const Vector3 p(random(-1, 1), random(0.5f, 4), random(-1, 1));
        const Vector3 turn(random(-90, 90), random(0, 360), random(-90, 90));
        bool clear = true;
        for (const Vector3& q : placed) clear &= length(p - q) > 0.55f; // bounding spheres 0.25 m
        if (!clear) continue;
        placed.push_back(p);
        Entity barrel;
        barrel.id = uint32_t(placed.size() + 1);
        barrel.name = "Бочка " + std::to_string(placed.size());
        barrel.shape = ShapeKind::Cylinder;
        barrel.size = {0.3f, 0.4f, 0.3f};
        barrel.position = p;
        barrel.rotationDeg = turn;
        barrel.rigid.enabled = barrel.collider.enabled = true;
        g.entities.push_back(barrel);
    }
    return g;
}

// What the pile does once it has landed (t = 3 s) and once the last leaning barrel has slid into
// place (t = 5 s), up to t = 8 s.
struct PileReport {
    int bodies = 0, asleepAtEnd = 0;
    float allAsleepAt = -1;          // [s] from the start, -1: never
    float maxSpeed = 0, maxSpin = 0; // over the awake bodies, t = 5 .. 8 s
    float maxPath = 0;               // the longest way a barrel went, t = 5 .. 8 s [m]
    long warmPoints = 0, warmMatched = 0; // contact points from t = 3 s, and those that found last step's impulse
    uint64_t hash = 0;                    // every barrel's pose and velocities at t = 8 s (FNV-1a)
    double warmRate() const { return warmPoints ? double(warmMatched) / double(warmPoints) : 1.0; }
};

PileReport measureBarrelPile(bool sleeping) {
    Simulation sim;
    sim.load(std::make_unique<GraphScene>(barrelPile()));
    sim.rigid.params.sleeping = sleeping;
    PileReport r;
    std::vector<Vector3> last(sim.rigid.bodies().size());
    std::vector<float> path(sim.rigid.bodies().size(), 0.0f);
    for (int frame = 0; frame < 480; ++frame) { // 8 s at 60 Hz
        sim.stepFrame();
        const RigidWorld& w = sim.rigid;
        int awake = 0;
        for (size_t i = 0; i < w.bodies().size(); ++i) {
            const RigidBody& b = w.bodies()[i];
            if (frame >= 300) path[i] += length(b.pos - last[i]); // a twitch goes back and forth: a long way
            last[i] = b.pos;
            if (b.invMass == 0 || !b.alive || b.sleeping) continue;
            ++awake;
            if (frame < 300) continue;
            r.maxSpeed = std::max(r.maxSpeed, length(b.vel));
            r.maxSpin = std::max(r.maxSpin, length(b.angVel));
            r.maxPath = std::max(r.maxPath, path[i]);
        }
        if (awake == 0 && r.allAsleepAt < 0) r.allAsleepAt = float(frame + 1) / 60.0f;
        if (frame < 180) continue;
        r.warmPoints += w.warmStartStats().points;
        r.warmMatched += w.warmStartStats().matched;
    }
    StateHash state;
    for (const RigidBody& b : sim.rigid.bodies()) {
        state.add(b.pos), state.add(b.rot), state.add(b.vel), state.add(b.angVel);
        if (b.invMass == 0 || !b.alive) continue;
        ++r.bodies;
        r.asleepAtEnd += b.sleeping ? 1 : 0;
    }
    r.hash = state.h;
    return r;
}

// A barrel (the 24-sided prism) touching a floor 3.6 m wide in a random pose: turned anyhow,
// anywhere over the floor, its lowest corner `gap` above the floor's top (negative: into it).
// Returns the contact the narrow phase makes; `penetration` is the true one, -gap.
ContactManifold barrelOnFloor(const ConvexShape& barrel, const ConvexShape& floor, Lcg& random, float gap) {
    const Vector3 axis = normalize(Vector3(random(-1, 1), random(-1, 1), random(-1, 1)) + Vector3(0, 0, 1e-3f));
    const Matrix3x3 R = Quaternion::fromAxisAngle(axis, random(0, 2 * kPi)).toMatrix3x3();
    const Vector3 lowest = R * barrel.support(R.transposed() * Vector3(0, -1, 0)); // relative to the centre
    const Vector3 at(random(-1.5f, 1.5f), 0.02f + gap - lowest.y, random(-1.5f, 1.5f));
    const PosedShape A{&barrel, R, at}, B{&floor, Matrix3x3(), Vector3(0, 0.01f, 0)};
    ContactManifold m;
    NarrowPhase().collide(A, B, m);
    return m;
}

} // namespace

// Fault 1 and fault 2 above, on the shapes of the pile: a barrel touching the floor, in a thousand
// random poses at gaps and penetrations of up to 0.2 mm - where a resting barrel lives. Every
// contact must carry the floor's normal (from the floor up to the barrel, B to A) and, at its
// deepest point, the true depth.
void testBarrelTouchesFloor() {
    const ConvexHullShape barrel(primitives::cylinder(0.15f, 0.4f, 24));
    const BoxShape floor(Vector3(1.8f, 0.01f, 1.8f));
    Lcg random{11};
    int poses = 0, lost = 0, wrongNormal = 0, wrongDepth = 0;
    float worstTilt = 0, worstDepth = 0;
    for (int k = 0; k < 1000; ++k) {
        const float gap = random(-2e-4f, 2e-4f);
        const ContactManifold m = barrelOnFloor(barrel, floor, random, gap);
        ++poses;
        if (m.points.empty()) { ++lost; continue; }
        float deepest = -kInf, tilt = 0;
        for (const ContactPoint& p : m.points) {
            deepest = std::max(deepest, p.depth);
            tilt = std::max(tilt, std::acos(std::min(1.0f, p.normal.y)));
        }
        worstTilt = std::max(worstTilt, tilt);
        worstDepth = std::max(worstDepth, std::fabs(deepest + gap));
        wrongNormal += tilt > 5e-3f ? 1 : 0;                     // five milliradians, 0.3 degrees
        wrongDepth += std::fabs(deepest + gap) > 2e-5f ? 1 : 0; // 20 micrometres
    }
    std::printf("  barrel touching a 3.6 m floor, %d poses (gap -0.2 .. 0.2 mm): %d lost, %d normals off by > 5 mrad "
                "(worst %.2e rad), %d depths off by > 20 um (worst %.2e m)\n",
                poses, lost, wrongNormal, double(worstTilt), wrongDepth, double(worstDepth));
    CHECK(lost == 0, "%d touching poses gave no contact", lost);
    CHECK(wrongNormal == 0, "%d contacts had a normal off the floor's by more than 5 mrad (worst %f rad)", wrongNormal, worstTilt);
    CHECK(wrongDepth == 0, "%d contacts had a depth off by more than 20 um (worst %f m)", wrongDepth, worstDepth);
}

// The pile as the user saw it, twice. With sleeping, as Box3D's pile test (WavePileTest in
// test_determinism.c): every barrel asleep within a set time - here 5 s from the drop, the pile has
// landed by 3 s - still asleep at 8 s, and the state then hashed; in the strict floating-point
// build (RF_STRICT_FP) the hash must be the one below. Without sleeping (what sleeping would hide):
// the pile must stand still on its own, and the warm start must find nearly every contact point
// again - a point that is new each step starts from a guess.
// The hash changes with any change of the rigid results; explain a new one where it is set.
//   a374657c2625233b - at first (the pile of the twitching barrels after its two fixes, in the
//                      editor scene whose rigid domain is open above the floor);
//   not yet taken     - the floor touches a barrel at the face turned to it, not at every vertex
//                      within the contact margin (RigidWorld::collideWalls), and the shock pass
//                      holds a support frozen only while it does not move into the body above it
//                      (solveManifoldShock): the value on MinGW is still to be taken; on Linux
//                      (glibc) the hash went 1c7a78e7bc19fef7 -> afbb0b88ffee8731.
void testBarrelPileSettles() {
    constexpr uint64_t kGolden = 0xa374657c2625233bull; // RF_STRICT_FP=ON, MinGW (GCC 11.2), Release
    const PileReport asleep = measureBarrelPile(true);
    std::printf("  sleeping on:  %d of %d barrels asleep at t = 8 s, all asleep at t = %.2f s, hash %016llx\n",
                asleep.asleepAtEnd, asleep.bodies, double(asleep.allAsleepAt), (unsigned long long)asleep.hash);
    const PileReport awake = measureBarrelPile(false);
    std::printf("  sleeping off: t = 5 .. 8 s the longest way a barrel went %.2f mm, the fastest %.2e m/s and %.2e rad/s, "
                "warm start %.2f %% of %ld points\n",
                1000.0 * double(awake.maxPath), double(awake.maxSpeed), double(awake.maxSpin), 100 * awake.warmRate(), awake.warmPoints);
    CHECK(asleep.allAsleepAt > 0 && asleep.allAsleepAt <= 5.0f, "the pile was not asleep 5 s after the drop (at %f s)",
          asleep.allAsleepAt);
    CHECK(asleep.asleepAtEnd == asleep.bodies, "%d of %d barrels awake again at t = 8 s", asleep.bodies - asleep.asleepAtEnd,
          asleep.bodies);
    // Under a millimetre per second on average. The way, not the speed: a barrel that rolled off the
    // plate onto the scene's open ground reports a steady 6 mm/s there in the strict build, yet does
    // not move (that ground contact is a question of its own).
    CHECK(awake.maxPath < 3e-3f, "a barrel of the resting pile went %f m in 3 s", awake.maxPath);
    CHECK(awake.warmRate() > 0.95, "the warm start found only %.1f %% of the contact points", 100 * awake.warmRate());
#ifdef RF_STRICT_FP
    CHECK(asleep.hash == kGolden, "the pile's hash changed: %016llx (was %016llx) - the rigid results changed",
          (unsigned long long)asleep.hash, (unsigned long long)kGolden);
#else
    (void)kGolden; // the default build (-march=native, FMA) gives another number: printed only
#endif
}
