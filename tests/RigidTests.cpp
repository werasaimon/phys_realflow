// Rigid bodies: collision geometry (GJK/EPA, SAT, CCD) against exact distances and times of
// impact, stacks and towers that must stand still, joints against the pendulum period, and the
// industry-hard contact cases (HardContactTests.h). Stability is checked by numbers: heights,
// drifts, velocities, overlaps.
#include "TestRunner.h"
#include "Tests.h"

#include "core/Probe.h"

void testRigid() {
    RigidWorld w;
    w.setDomain(AABB({-5, 0, -5}, {5, 10, 5}));
    int s = w.addSphere({0, 2, 0}, 0.25f, 1000, Vector3(1));
    int b = w.addBox({1.5f, 3, 0}, {0.2f, 0.3f, 0.4f}, Quaternion::fromAxisAngle({1, 1, 0}, 0.6f), 1000, Vector3(1));
    for (int i = 0; i < 60 * 6; ++i)
        for (int k = 0; k < w.params.substeps; ++k) w.step(1.0f / 60 / w.params.substeps);
    const RigidBody& S = w.bodies()[s];
    const RigidBody& B = w.bodies()[b];
    CHECK(std::fabs(S.pos.y - 0.25f) < 0.02f, "sphere rest height %f", S.pos.y);
    CHECK(length(S.vel) < 0.05f, "sphere speed %f", length(S.vel));
    float boxLow = B.worldBounds().lo.y;
    CHECK(std::fabs(boxLow) < 0.03f, "box lowest point %f", boxLow);
    CHECK(length(B.vel) < 0.1f, "box speed %f", length(B.vel));
}

void testGjkEpa() {
    SphereShape s1(1.0f), s2(0.5f);
    PosedShape A{&s1, Matrix3x3(), {0, 0, 0}}, B{&s2, Matrix3x3(), {3, 0, 0}};
    GjkResult g = gjk(A, B);
    CHECK(!g.intersect && std::fabs(g.distance - 1.5f) < 1e-3f, "sphere distance %f", g.distance);

    BoxShape b1({0.5f, 0.5f, 0.5f}), b2({0.5f, 0.5f, 0.5f});
    PosedShape BA{&b1, Matrix3x3(), {0, 0, 0}}, BB{&b2, Matrix3x3(), {0.9f, 0.2f, -0.1f}};
    PenetrationResult pr;
    CHECK(penetration(BA, BB, pr), "overlapping boxes");
    CHECK(std::fabs(pr.depth - 0.1f) < 1e-3f && pr.normal.x < -0.99f, "EPA depth %f normal %f %f %f", pr.depth,
          pr.normal.x, pr.normal.y, pr.normal.z);
    // Same pair through SAT: identical depth, 4-point face manifold.
    ContactManifold m;
    CHECK(NarrowPhase::boxBox(BA, BB, m) && m.points.size() == 4, "SAT face manifold size %zu", m.points.size());
    for (auto& c : m.points) CHECK(std::fabs(c.depth - 0.1f) < 1e-3f, "SAT depth %f", c.depth);

    // Hull vs box through GJK/EPA + perturbation must also produce a multi-point manifold.
    ConvexHullShape hull(primitives::box({0.5f, 0.5f, 0.5f}));
    PosedShape H{&hull, Matrix3x3(), {0, 0.95f, 0}}, F{&b1, Matrix3x3(), {0, 0, 0}};
    ContactManifold mh;
    CHECK(NarrowPhase::convexConvex(H, F, mh), "hull-box contact");
    CHECK(mh.points.size() >= 3, "perturbation manifold has %zu points", mh.points.size());
    for (auto& c : mh.points) CHECK(c.normal.y > 0.99f, "normal must push the hull up");

    // Edge-edge: two cubes rotated about perpendicular axes, crossing edges.
    Matrix3x3 Rz = Quaternion::fromAxisAngle({0, 0, 1}, kPi / 4).toMatrix3x3(), Rx = Quaternion::fromAxisAngle({1, 0, 0}, kPi / 4).toMatrix3x3();
    float reach = 0.5f * std::sqrt(2.0f);
    PosedShape E1{&b1, Rz, {0, 0, 0}}, E2{&b2, Rx, {0, 2 * reach - 0.05f, 0}};
    ContactManifold me;
    CHECK(NarrowPhase::boxBox(E2, E1, me) && me.points.size() == 1, "edge-edge contact points %zu", me.points.size());
    if (!me.points.empty())
        CHECK(std::fabs(me.points[0].depth - 0.05f) < 2e-3f && me.points[0].normal.y > 0.99f, "edge depth %f",
              me.points[0].depth);
}

void testBoxStack() {
    RigidWorld w;
    w.setDomain(AABB({-5, 0, -5}, {5, 10, 5}));
    const float h = 0.25f;
    std::vector<int> ids;
    for (int i = 0; i < 6; ++i) ids.push_back(w.addBox({0.01f * (i % 2), h + i * 2 * h, 0}, Vector3(h), Quaternion(), 500, Vector3(1)));
    for (int f = 0; f < 60 * 5; ++f)
        for (int k = 0; k < w.params.substeps; ++k) w.step(1.0f / 60 / w.params.substeps);
    const RigidBody& top = w.bodies()[ids.back()];
    float expectY = h + 5 * 2 * h;
    std::printf("  stack: top y=%.4f (expect %.4f)  x drift=%.4f  speed=%.4f  contacts=%zu\n", top.pos.y, expectY,
                top.pos.x - 0.01f, length(top.vel), w.contactCount());
    CHECK(std::fabs(top.pos.y - expectY) < 0.03f, "stack collapsed/sank: top y %f", top.pos.y);
    CHECK(std::fabs(top.pos.x - 0.01f) < 0.02f, "stack drifted %f", top.pos.x);
    CHECK(length(top.vel) < 0.05f, "stack not at rest %f", length(top.vel));
}

void testBroadPhase() {
    std::mt19937 rng(3);
    std::uniform_real_distribution<float> U(-3, 3), S(0.05f, 0.6f);
    std::vector<AABB> boxes;
    for (int i = 0; i < 400; ++i) {
        Vector3 c(U(rng), U(rng), U(rng)), h(S(rng), S(rng), S(rng));
        boxes.emplace_back(c - h, c + h);
    }
    BruteForceBroadPhase bf;
    BvhBroadPhase bvh;
    bf.update(boxes);
    bvh.update(boxes);
    std::vector<std::pair<int, int>> a, b;
    bf.findPairs(a);
    bvh.findPairs(b);
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    CHECK(a == b && !a.empty(), "BVH pairs %zu vs brute force %zu", b.size(), a.size());

    // Incremental sweep and prune: exact against brute force over many frames of coherent motion,
    // with touching boxes (equal endpoints), and after a teleport (large reorder).
    SweepAndPruneBroadPhase sap;
    AABBTreeBroadPhase treeBp; // exact boxes: must equal brute force
    int treeMismatches = 0;
    std::vector<Vector3> vel(boxes.size());
    std::uniform_real_distribution<float> V(-0.05f, 0.05f);
    for (auto& v : vel) v = Vector3(V(rng), V(rng), V(rng));
    boxes[1] = AABB(boxes[0].hi, boxes[0].hi + Vector3(0.2f)); // corner contact
    int mismatches = 0;
    size_t swaps = 0;
    for (int f = 0; f < 200; ++f) {
        if (f == 100)
            for (size_t i = 0; i < boxes.size(); i += 7) { // teleport some boxes
                Vector3 c(U(rng), U(rng), U(rng));
                Vector3 h = boxes[i].extent() * 0.5f;
                boxes[i] = AABB(c - h, c + h);
            }
        for (size_t i = 2; i < boxes.size(); ++i) {
            boxes[i].lo += vel[i];
            boxes[i].hi += vel[i];
        }
        bf.update(boxes);
        sap.update(boxes);
        treeBp.update(boxes);
        bf.findPairs(a);
        sap.findPairs(b);
        std::sort(a.begin(), a.end());
        mismatches += a != b;
        treeBp.findPairs(b);
        treeMismatches += a != b;
        if (f > 0 && f != 100) swaps += sap.lastSwaps();
    }
    std::printf("  sweep and prune: 200 frames, %d mismatches, %.0f endpoint swaps per frame (%zu endpoints)\n", mismatches,
                swaps / 198.0, boxes.size() * 6);
    CHECK(mismatches == 0, "sweep and prune differs from brute force in %d frames", mismatches);
    CHECK(treeMismatches == 0, "AABB tree broad phase differs from brute force in %d frames", treeMismatches);
}

void testTallStack() {
    RigidWorld w;
    w.setDomain(AABB({-5, 0, -5}, {5, 20, 5}));
    const float h = 0.25f;
    for (int i = 0; i < 10; ++i) w.addBox({0.01f * (i % 2), h + i * 2 * h, 0}, Vector3(h), Quaternion(), 500, Vector3(1));
    float maxW = 0;
    for (int f = 0; f < 60 * 5; ++f) {
        for (int k = 0; k < w.params.substeps; ++k) w.step(1.0f / 60 / w.params.substeps);
        if (f > 60)
            for (auto& b : w.bodies()) maxW = std::max(maxW, length(b.angVel));
    }
    const RigidBody& top = w.bodies().back();
    std::printf("  10-stack: top y=%.4f (expect %.4f) x=%.4f  max|w|=%.4f  contacts=%zu pairs=%zu\n", top.pos.y,
                h + 18 * h, top.pos.x, maxW, w.contactCount(), w.pairCount());
    CHECK(std::fabs(top.pos.y - (h + 18 * h)) < 0.03f && maxW < 0.2f, "10-box stack unstable");
}

void testStack100() {
    // Standard stability test: 100 cubes, each dropped from 1 cm above the one below.
    RigidWorld w;
    w.setDomain(AABB({-5, 0, -5}, {5, 60, 5}));
    const float h = 0.1f, gap = 0.01f;
    for (int i = 0; i < 100; ++i) w.addBox({0, h + i * (2 * h + gap), 0}, Vector3(h), Quaternion(), 500, Vector3(1));
    float maxW = 0, maxV = 0;
    for (int f = 0; f < 60 * 6; ++f) {
        for (int k = 0; k < w.params.substeps; ++k) w.step(1.0f / 60 / w.params.substeps);
        if (f >= 60 * 5) // after the impact cascade has died out: the last second
            for (auto& b : w.bodies()) {
                maxW = std::max(maxW, length(b.angVel));
                maxV = std::max(maxV, length(b.vel));
            }
    }
    const RigidBody& top = w.bodies().back();
    float expect = h + 99 * 2 * h;
    float drift = std::sqrt(top.pos.x * top.pos.x + top.pos.z * top.pos.z);
    std::printf("  100-stack: top y=%.4f (expect %.4f)  drift=%.4f  max|v|=%.4f max|w|=%.4f  contacts=%zu\n", top.pos.y,
                expect, drift, maxV, maxW, w.contactCount());
    // Stands: full height (no box slipped out), top still above the base, no body flying around.
    CHECK(std::fabs(top.pos.y - expect) < 0.05f, "stack height %f vs %f", top.pos.y, expect);
    CHECK(drift < 0.05f && maxV < 0.05f && maxW < 0.1f, "100-box stack not at rest");
}

void testStack200() {
    // Twice the classic stability test: 200 cubes, each dropped from 1 cm above the one below -
    // 40 m of stack, 199 impacts cascading down through it. It must stand at full height, stop
    // and go to sleep. The mean time of a rigid step is printed: ~200 manifolds with the block
    // solver, islands and shock propagation.
    RigidWorld w;
    w.setDomain(AABB({-5, 0, -5}, {5, 120, 5}));
    const float h = 0.1f, gap = 0.01f;
    const int n = 200;
    for (int i = 0; i < n; ++i) w.addBox({0, h + i * (2 * h + gap), 0}, Vector3(h), Quaternion(), 500, Vector3(1));
    float maxW = 0, maxV = 0;
    double stepMs = 0;
    int steps = 0;
    const int frames = 60 * 10;
    for (int f = 0; f < frames; ++f) {
        for (int k = 0; k < w.params.substeps; ++k) {
            const auto t0 = std::chrono::steady_clock::now();
            w.step(1.0f / 60 / w.params.substeps);
            stepMs += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            ++steps;
        }
        if (f >= frames - 60) // the last second
            for (auto& b : w.bodies()) {
                maxW = std::max(maxW, length(b.angVel));
                maxV = std::max(maxV, length(b.vel));
            }
    }
    const RigidBody& top = w.bodies().back();
    const float expect = h + (n - 1) * 2 * h;
    const float drift = std::sqrt(top.pos.x * top.pos.x + top.pos.z * top.pos.z);
    bool finite = true;
    float tilt = 0; // the largest lean of any box off the vertical [rad]
    for (const RigidBody& b : w.bodies()) {
        finite &= std::isfinite(b.pos.x + b.pos.y + b.pos.z + b.vel.x + b.vel.y + b.vel.z);
        tilt = std::max(tilt, std::acos(std::min(1.0f, (b.rotation() * Vector3(0, 1, 0)).y)));
    }
    std::printf("  200-stack: top y=%.4f (expect %.4f)  sideways %.4f m  max tilt %.2e rad  max|v|=%.4f max|w|=%.4f  asleep %d/%d  mean step %.2f ms\n",
                top.pos.y, expect, drift, tilt, maxV, maxW, int(w.sleepingCount()), n, stepMs / steps);
    CHECK(finite, "NaN in the 200-box stack");
    CHECK(std::fabs(top.pos.y - expect) < 0.01f, "stack height %f vs %f", top.pos.y, expect);
    // 0.11 degrees: the strict-FP build settles at 1.09e-3 rad, the FMA build at 6e-4 (same physics, other rounding).
    CHECK(tilt < 2e-3f, "boxes lean over: %e rad", tilt);
    CHECK(maxV < 0.05f && maxW < 0.1f, "200-box stack not at rest: v %f w %f", maxV, maxW);
    CHECK(int(w.sleepingCount()) == n, "only %d of %d boxes asleep", int(w.sleepingCount()), n);
    // Known limitation, measured: the boxes stay upright but slide sideways over each other during
    // the cascade of impacts (a zigzag of up to 9 cm at 40 m with 200 boxes, 1.7 cm at 20 m with
    // 100) - the friction pass cannot hold them while the impact wave unloads the contacts. The
    // sideways offset is reported above; the target for the solver work is under 1 cm.
    CHECK(drift < 0.5f, "the stack slid apart: %f m", drift);
}

void testRaycastGrab() {
    RigidWorld w;
    w.setDomain(AABB({-5, 0, -5}, {5, 10, 5}));
    int box = w.addBox({0, 0.25f, 0}, Vector3(0.25f), Quaternion(), 500, Vector3(1));
    int hull = w.addConvex(primitives::cylinder(0.2f, 0.4f, 16), {1.5f, 0.3f, 0}, Quaternion(), 500, Vector3(1));
    int body;
    float t;
    Vector3 n;
    // Straight down onto the box top face.
    CHECK(w.raycast({0.1f, 5, 0.05f}, {0, -1, 0}, 100, body, t, n) && body == box && std::fabs(t - 4.5f) < 1e-3f && n.y > 0.99f,
          "ray vs box: body %d t %f", body, t);
    // Horizontal ray at the cylinder (convex hull, Cyrus-Beck).
    CHECK(w.raycast({-2, 0.3f, 0}, {1, 0, 0}, 100, body, t, n) && body == box, "first hit must be the box, got %d", body);
    CHECK(w.raycast({5, 0.3f, 0}, {-1, 0, 0}, 100, body, t, n) && body == hull && std::fabs(t - 3.3f) < 0.02f,
          "ray vs hull: body %d t %f", body, t);
    // Grab the box at its top and pull the target 1 m up: the box must follow.
    w.grab(box, {0, 0.5f, 0});
    w.setGrabTarget({0, 1.5f, 0});
    for (int f = 0; f < 120; ++f)
        for (int k = 0; k < w.params.substeps; ++k) w.step(1.0f / 60 / w.params.substeps);
    const RigidBody& B = w.bodies()[box];
    Vector3 anchor = B.pos + B.rotation() * w.grabJoint().localAnchor;
    std::printf("  grab: anchor (%.3f %.3f %.3f) target (0 1.5 0)\n", anchor.x, anchor.y, anchor.z);
    CHECK(length(anchor - Vector3(0, 1.5f, 0)) < 0.05f, "grabbed box did not follow the target");
    w.releaseGrab();
}

void testJoints() {
    RigidWorld w;
    w.setDomain(AABB({-10, -10, -10}, {10, 10, 10}));
    w.params.sleeping = false;
    const int sub = w.params.substeps;
    auto run = [&](float seconds, const std::function<void()>& each = {}) {
        for (int f = 0; f < int(seconds * 60); ++f) {
            for (int k = 0; k < sub; ++k) w.step(1.0f / 60 / sub);
            if (each) each();
        }
    };

    // 1) Ball joint pendulum: period of a (nearly) point mass, T = 2 pi sqrt(L/g).
    int bob = w.addSphere({1.0f, 5, 0}, 0.05f, 5000, Vector3(1));
    BallJoint& ball = w.addBallJoint(bob, -1, {0, 5, 0});
    (void)ball;
    // Release from 10 degrees.
    w.bodies()[bob].pos = Vector3(std::sin(0.1745f), 5 - std::cos(0.1745f), 0);
    w.joints().back()->localAnchorA = w.bodies()[bob].rot.conjugate().rotate(Vector3(0, 5, 0) - w.bodies()[bob].pos);
    float prevVx = 0, firstCross = -1, lastCross = -1, t = 0, maxDrift = 0;
    int crossings = 0;
    for (int f = 0; f < 60 * 7; ++f) {
        for (int k = 0; k < sub; ++k) w.step(1.0f / 60 / sub);
        t += 1.0f / 60;
        float vx = w.bodies()[bob].vel.x;
        if (f > 0 && prevVx < 0 && vx >= 0) {
            if (firstCross < 0) firstCross = t;
            lastCross = t;
            ++crossings;
        }
        prevVx = vx;
        maxDrift = std::max(maxDrift, length(w.joints().back()->worldAnchorA(w.bodies()) - Vector3(0, 5, 0)));
    }
    float T = crossings > 1 ? (lastCross - firstCross) / float(crossings - 1) : 0;
    float Te = 2 * kPi * std::sqrt(1.0f / 9.81f);
    std::printf("  pendulum: T=%.4f s (theory %.4f)  anchor drift=%.2e m\n", T, Te, maxDrift);
    CHECK(std::fabs(T - Te) / Te < 0.03f, "pendulum period %f vs %f", T, Te);
    CHECK(maxDrift < 1e-3f, "ball joint drift %f", maxDrift);

    // 2) Hinge: spin a door about all axes; only the hinge axis may keep rotating.
    int door = w.addBox({3, 5, 0}, {0.3f, 0.5f, 0.03f}, Quaternion(), 500, Vector3(1));
    HingeJoint& hinge = w.addHingeJoint(door, -1, {2.7f, 5, 0}, {0, 1, 0});
    w.bodies()[door].angVel = {3, 2, 3};
    float offAxis = 0, hingeDrift = 0;
    run(2.0f, [&] {
        Vector3 wv = w.bodies()[door].angVel;
        offAxis = std::max(offAxis, std::sqrt(wv.x * wv.x + wv.z * wv.z));
        hingeDrift = std::max(hingeDrift, length(hinge.worldAnchorA(w.bodies()) - Vector3(2.7f, 5, 0)));
    });
    std::printf("  hinge: off-axis |w|=%.2e rad/s  anchor drift=%.2e m  angle=%.3f rad\n", offAxis, hingeDrift, hinge.angle(w.bodies()));
    CHECK(offAxis < 1e-2f && hingeDrift < 1e-3f, "hinge does not hold its axis");

    // 3) Hinge motor reaches its speed.
    int wheel = w.addBox({6, 5, 0}, {0.4f, 0.05f, 0.1f}, Quaternion(), 500, Vector3(1));
    HingeJoint& motor = w.addHingeJoint(wheel, -1, {6, 5, 0}, {0, 0, 1});
    motor.motorEnabled = true;
    motor.motorSpeed = 3.0f;
    motor.maxMotorTorque = 500.0f;
    run(1.0f);
    float wz = w.bodies()[wheel].angVel.z;
    std::printf("  motor: w=%.3f rad/s (target 3)\n", wz);
    CHECK(std::fabs(wz - 3.0f) < 0.15f, "motor speed %f", wz);

    // 4) Slider: gravity pulls along a tilted rail only, stops at the limit.
    Vector3 axis = normalize(Vector3(1, -1, 0));
    int cart = w.addBox({-3, 5, 0}, Vector3(0.1f), Quaternion(), 500, Vector3(1));
    SliderJoint& slider = w.addSliderJoint(cart, -1, axis);
    slider.limitEnabled = true;
    slider.lower = -0.5f;
    slider.upper = 0.5f;
    float perpErr = 0;
    run(2.0f, [&] {
        Vector3 d = w.bodies()[cart].pos - Vector3(-3, 5, 0);
        perpErr = std::max(perpErr, length(d - axis * dot(d, axis)));
    });
    float s = slider.translation(w.bodies());
    std::printf("  slider: travel=%.4f (limit %.2f)  off-rail=%.2e m  |w|=%.2e\n", s, -0.5f, perpErr, length(w.bodies()[cart].angVel));
    CHECK(perpErr < 1e-3f && std::fabs(std::fabs(s) - 0.5f) < 0.01f, "slider off rail or past limit");

    // 5) Fixed: a welded pair hitting the ground keeps its relative pose.
    RigidWorld g;
    g.setDomain(AABB({-5, 0, -5}, {5, 10, 5}));
    int a = g.addBox({0, 2, 0}, {0.3f, 0.05f, 0.05f}, Quaternion(), 500, Vector3(1));
    int b = g.addBox({-0.25f, 2.25f, 0}, {0.05f, 0.25f, 0.05f}, Quaternion(), 500, Vector3(1));
    g.addFixedJoint(a, b);
    g.bodies()[a].angVel = g.bodies()[b].angVel = Vector3(2, 1, 3);
    Vector3 rel0 = g.bodies()[b].rot.conjugate().rotate(g.bodies()[a].pos - g.bodies()[b].pos);
    float weldErr = 0;
    for (int f = 0; f < 120; ++f) {
        for (int k = 0; k < g.params.substeps; ++k) g.step(1.0f / 60 / g.params.substeps);
        Vector3 rel = g.bodies()[b].rot.conjugate().rotate(g.bodies()[a].pos - g.bodies()[b].pos);
        weldErr = std::max(weldErr, length(rel - rel0));
    }
    std::printf("  fixed: max relative offset change=%.2e m\n", weldErr);
    CHECK(weldErr < 5e-3f, "weld deformed by %f", weldErr);

    // 6) Distance rod keeps its length; rope may go slack but never stretches.
    int r1 = w.addSphere({0, 0, 5}, 0.1f, 1000, Vector3(1));
    DistanceJoint& rod = w.addDistanceJoint(r1, -1, {0, 0, 5}, {0, 1.5f, 5});
    w.bodies()[r1].vel = {3, 0, 1};
    float rodErr = 0;
    run(2.0f, [&] { rodErr = std::max(rodErr, std::fabs(length(rod.worldAnchorA(w.bodies()) - Vector3(0, 1.5f, 5)) - rod.length)); });
    std::printf("  distance rod: max length error=%.2e m\n", rodErr);
    CHECK(rodErr < 2e-3f, "rod length error %f", rodErr);
}

void testCcd() {
    auto shoot = [](bool ccd, float& finalX) {
        RigidWorld w;
        w.setDomain(AABB({-5, -5, -5}, {5, 5, 5}));
        w.params.gravity = Vector3(0.0f);
        w.params.ccd = ccd;
        w.addBox({0, 0, 0}, {0.01f, 1, 1}, Quaternion(), 0, Vector3(1)); // 2 cm static wall
        int bullet = w.addSphere({-2.23f, 0, 0}, 0.02f, 8000, Vector3(1)); // steps land at -0.23 and +0.27: across the wall
        w.bodies()[bullet].vel = {300, 0, 0}; // 0.5 m per substep
        finalX = -1e9f; // furthest x reached: > 0 means it went through the wall
        for (int f = 0; f < 30; ++f)
            for (int k = 0; k < w.params.substeps; ++k) {
                w.step(1.0f / 60 / w.params.substeps);
                finalX = std::max(finalX, w.bodies()[bullet].pos.x);
            }
        return w.ccdHits();
    };
    float xOff, xOn;
    shoot(false, xOff);
    shoot(true, xOn);
    std::printf("  bullet 300 m/s vs 2 cm wall, furthest x: without CCD %.3f, with CCD %.3f\n", xOff, xOn);
    CHECK(xOff > 0.1f, "test setup: the bullet should tunnel without CCD (x=%f)", xOff);
    CHECK(xOn < 0.0f, "CCD must stop the bullet in front of the wall (x=%f)", xOn);

    // Bullet against a free (dynamic) box: it must hit it, bounce and transfer momentum.
    {
        RigidWorld w;
        w.setDomain(AABB({-50, -50, -50}, {50, 50, 50}));
        w.params.gravity = Vector3(0.0f);
        int box = w.addBox({0, 0, 0}, {0.02f, 0.3f, 0.3f}, Quaternion(), 500, Vector3(1)); // 4 cm thick, 3.6 kg
        int bullet = w.addSphere({-2.23f, 0, 0}, 0.02f, 8000, Vector3(1));             // 0.27 kg
        w.bodies()[box].restitution = w.bodies()[bullet].restitution = 0.5f;
        w.bodies()[bullet].vel = {300, 0, 0};
        float p0 = w.bodies()[bullet].mass * 300.0f;
        float maxBulletX = -1e9f;
        for (int f = 0; f < 2; ++f) // just past the impact, before anything reaches the domain walls
            for (int k = 0; k < w.params.substeps; ++k) {
                w.step(1.0f / 60 / w.params.substeps);
                maxBulletX = std::max(maxBulletX, w.bodies()[bullet].pos.x - w.bodies()[box].pos.x);
            }
        const RigidBody &B = w.bodies()[box], &S = w.bodies()[bullet];
        float p1 = B.mass * B.vel.x + S.mass * S.vel.x;
        std::printf("  bullet vs free box: box v=%.2f m/s, bullet v=%.2f m/s, momentum %.2f -> %.2f\n", B.vel.x, S.vel.x, p0, p1);
        CHECK(maxBulletX < 0.0f, "bullet went through the moving box");
        CHECK(B.vel.x > 5.0f && S.vel.x < B.vel.x, "no impact response: box %f bullet %f", B.vel.x, S.vel.x);
        CHECK(std::fabs(p1 - p0) / p0 < 0.02f, "momentum not conserved: %f -> %f", p0, p1);
        // Restitution 0.5: relative separation speed = 0.5 * 300.
        float sep = B.vel.x - S.vel.x;
        CHECK(std::fabs(sep - 150.0f) / 150.0f < 0.1f, "restitution: separation speed %f (expected 150)", sep);
    }

    // Fast box against a thin static mesh plate (triangles through the BVH).
    RigidWorld w;
    w.setDomain(AABB({-5, -5, -5}, {5, 5, 5}));
    w.params.gravity = Vector3(0.0f);
    TriMesh plate = primitives::box({1, 0.005f, 1});
    MeshBVH bvh;
    bvh.build(plate);
    w.setStaticMesh(&bvh);
    int cube = w.addBox({0, 2, 0}, Vector3(0.05f), Quaternion::fromAxisAngle({1, 0, 1}, 0.5f), 500, Vector3(1));
    w.bodies()[cube].vel = {0, -200, 0};
    for (int f = 0; f < 30; ++f)
        for (int k = 0; k < w.params.substeps; ++k) w.step(1.0f / 60 / w.params.substeps);
    std::printf("  cube 200 m/s vs 1 cm mesh plate: final y=%.3f\n", w.bodies()[cube].pos.y);
    CHECK(w.bodies()[cube].pos.y > 0.0f, "cube tunnelled through the mesh (y=%f)", w.bodies()[cube].pos.y);

    // Time of impact accuracy: sphere r=0.5 at x=-3 moving +6 over the step towards a sphere at x=2.
    SphereShape s1(0.5f), s2(0.5f);
    SweptPose A{&s1, {-3, 0, 0}, {3, 0, 0}, Quaternion(), Vector3(0.0f)};
    SweptPose B{&s2, {2, 0, 0}, {2, 0, 0}, Quaternion(), Vector3(0.0f)};
    ToiResult r = timeOfImpact(A, B, 1e-3f);
    // contact when centre reaches x = 1  ->  s = 4/6
    std::printf("  TOI: s=%.4f (exact %.4f), %d iterations\n", r.s, 4.0f / 6.0f, r.iterations);
    CHECK(r.hit && std::fabs(r.s - 4.0f / 6.0f) < 1e-3f, "TOI %f", r.s);
}

void testCcdBodies() {
    // 1) Head-on: two thin plates at 200 m/s each (closing 400 m/s = 0.67 m per substep, 4 cm thick
    //    together): without CCD they pass through each other.
    auto headOn = [](bool ccd, float& minGap, float& worst) {
        RigidWorld w;
        w.setDomain(AABB({-50, -50, -50}, {50, 50, 50}));
        w.params.gravity = Vector3(0.0f);
        w.params.ccd = ccd;
        int a = w.addBox({-2.1f, 0, 0}, {0.01f, 0.3f, 0.3f}, Quaternion(), 1000, Vector3(1));
        int b = w.addBox({2.1f, 0.05f, 0}, {0.01f, 0.3f, 0.3f}, Quaternion(), 1000, Vector3(1));
        w.bodies()[a].vel = {200, 0, 0};
        w.bodies()[b].vel = {-200, 0, 0};
        minGap = 1e9f;
        worst = 0;
        for (int f = 0; f < 20; ++f)
            for (int k = 0; k < w.params.substeps; ++k) {
                w.step(1.0f / 60 / w.params.substeps);
                minGap = std::min(minGap, w.bodies()[b].pos.x - w.bodies()[a].pos.x); // < 0: passed through
                worst = std::max(worst, maxOverlap(w));
            }
        return w.bodies()[a].vel.x;
    };
    float gapOff, worstOff, gapOn, worstOn;
    headOn(false, gapOff, worstOff);
    float va = headOn(true, gapOn, worstOn);
    std::printf("  head-on plates 2x200 m/s: without CCD min gap %.3f; with CCD min gap %.3f, max overlap %.4f m, A bounces at %.1f m/s\n",
                gapOff, gapOn, worstOn, va);
    CHECK(gapOff < 0, "test setup: plates should pass through each other without CCD");
    CHECK(gapOn > 0 && worstOn < 0.01f, "moving bodies interpenetrated (gap %f, overlap %f)", gapOn, worstOn);
    CHECK(va < 0, "no collision response: A keeps moving forward at %f", va);

    // 2) A spray of fast bullets into a free-standing stack of boxes: never any superposition.
    RigidWorld w;
    w.setDomain(AABB({-5, 0, -5}, {5, 10, 5}));
    for (int i = 0; i < 5; ++i) w.addBox({1.0f, 0.1f + 0.2f * i, 0}, Vector3(0.1f), Quaternion(), 500, Vector3(1));
    std::vector<int> bullets;
    for (int i = 0; i < 5; ++i) {
        bullets.push_back(w.addSphere({-3.13f, 0.1f + 0.2f * i, 0.02f * i}, 0.02f, 8000, Vector3(1)));
        w.bodies()[bullets.back()].vel = {250, 0, 0};
    }
    float worst = 0;
    for (int f = 0; f < 60; ++f)
        for (int k = 0; k < w.params.substeps; ++k) {
            w.step(1.0f / 60 / w.params.substeps);
            worst = std::max(worst, maxOverlap(w));
        }
    float moved = 0;
    for (int i = 0; i < 5; ++i) moved = std::max(moved, w.bodies()[i].pos.x - 1.0f);
    std::printf("  bullets into a stack: max overlap %.4f m, stack pushed by %.2f m\n", worst, moved);
    CHECK(worst < 0.01f, "superposition between bodies: %f m", worst);
    CHECK(moved > 0.05f, "the stack must be hit (moved %f)", moved);
}

void testCcdSpinningPlate() {
    // A 1 m plate spinning at 200 rad/s around its centre, 0.3 m from a 2 cm static post: per
    // substep its edge sweeps 10 cm at the post (plate + post are 6 cm thick together).
    // The hit is measured by the angular momentum |L| the plate keeps, not by one component of w:
    // a thin plate struck off-centre tumbles (I about its long axis is 300 times smaller), and
    // the gyroscopic term then trades w between axes at constant |L|.
    auto spin = [](bool ccd, float& wFinal, float& lRatio, float& worstOverlap, int& hits) {
        RigidWorld w;
        w.setDomain(AABB({-5, -5, -5}, {5, 5, 5}));
        w.params.gravity = Vector3(0.0f);
        w.params.ccd = ccd;
        w.params.sleeping = false;
        w.addBox({0, 0, 0}, {0.01f, 1.0f, 0.01f}, Quaternion(), 0, Vector3(1));              // post
        int plate = w.addBox({0, 0, 0.3f}, {0.5f, 0.02f, 0.02f}, Quaternion(), 800, Vector3(1)); // along x
        w.bodies()[plate].angVel = {0, 200, 0};
        auto momentum = [&] { // |L| = |I_world w|
            const RigidBody& b = w.bodies()[plate];
            const Vector3 wl = b.rotation().transposed() * b.angVel;
            return length(Vector3(wl.x / b.invInertiaLocal.x, wl.y / b.invInertiaLocal.y, wl.z / b.invInertiaLocal.z));
        };
        const float L0 = momentum();
        worstOverlap = 0;
        hits = 0;
        for (int k = 0; k < 60; ++k) { // 0.1 s: ~20 revolutions without an obstacle
            w.step(1.0f / 600);
            hits += int(w.ccdHits());
            worstOverlap = std::max(worstOverlap, maxOverlap(w));
        }
        wFinal = w.bodies()[plate].angVel.y;
        lRatio = momentum() / L0;
    };
    float wOff, lOff, ovOff, wOn, lOn, ovOn;
    int hOff, hOn;
    spin(false, wOff, lOff, ovOff, hOff);
    spin(true, wOn, lOn, ovOn, hOn);
    std::printf("  spinning plate 200 rad/s vs 2 cm post: without CCD w_y=%.1f |L| %.0f%% overlap %.4f m; with CCD w_y=%.1f |L| %.0f%% overlap %.4f m, %d CCD stops\n",
                wOff, 100 * lOff, ovOff, wOn, 100 * lOn, ovOn, hOn);
    // One hit at r = 0.3 m from the centre of a free plate: the impulse that stops the point of
    // contact takes 1 - (I/m) / (I/m + r^2) ~ 52 % of the spin, more with the bounce.
    CHECK(lOn < 0.6f, "with CCD the plate must hit the post (kept %.0f%% of |L|)", 100 * lOn);
    CHECK(ovOn < 0.01f, "plate and post interpenetrated: %f m", ovOn);
}

void testGjkRandomThin() {
    // GJK distance / intersection against brute-force surface sampling for thin, long boxes in
    // random poses (the case that fooled the float tetrahedron test).
    std::mt19937 rng(11);
    std::uniform_real_distribution<float> U(-1, 1), A(0, 6.2832f);
    auto samples = [](const PosedShape& P, const Vector3& h) {
        std::vector<Vector3> v;
        const int n = 16;
        for (int i = 0; i <= n; ++i)
            for (int j = 0; j <= n; ++j)
                for (int k = 0; k <= n; ++k) {
                    if (i && j && k && i < n && j < n && k < n) continue;
                    v.push_back(P.p + P.R * Vector3(-h.x + 2 * h.x * i / n, -h.y + 2 * h.y * j / n, -h.z + 2 * h.z * k / n));
                }
        return v;
    };
    int wrong = 0, total = 300;
    for (int t = 0; t < total; ++t) {
        Vector3 ha(0.5f, 0.02f, 0.02f), hb(0.01f + 0.1f * std::fabs(U(rng)), 0.6f, 0.015f);
        BoxShape a(ha), b(hb);
        PosedShape PA{&a, Quaternion::fromAxisAngle({U(rng), U(rng), U(rng) + 1e-3f}, A(rng)).toMatrix3x3(), Vector3(U(rng), U(rng), U(rng)) * 0.4f};
        PosedShape PB{&b, Quaternion::fromAxisAngle({U(rng), U(rng), U(rng) + 1e-3f}, A(rng)).toMatrix3x3(), Vector3(U(rng), U(rng), U(rng)) * 0.4f};
        GjkResult g = gjk(PA, PB);
        ContactManifold m;
        bool sat = NarrowPhase::boxBox(PA, PB, m) && !m.points.empty() && m.points[0].depth > 0;
        float br = exactBoxDistance(PA, ha, PB, hb);
        bool bad;
        if (sat) bad = !g.intersect && g.distance > 1e-3f; // SAT says overlapping
        else bad = g.intersect ? br > 1e-3f : std::fabs(g.distance - br) > 1e-3f; // 1 mm, exact reference
        (void)samples;
        if (bad)
            std::printf("    case %d: gjk intersect=%d dist=%.4f  SAT=%d  brute=%.4f\n", t, int(g.intersect), g.distance, int(sat), br);
        wrong += bad;
    }
    std::printf("  GJK vs exact distance on %d thin-box poses (tol 1 mm): %d wrong\n", total, wrong);
    CHECK(wrong == 0, "GJK wrong in %d of %d cases", wrong, total);
}

void testConvexHullAndDecomposition() {
    // 1) Quickhull: closed 2-manifold (every directed edge has its twin), Euler F = 2V - 4, and all
    //    input points inside; also for dense nearly coplanar samples (grid points on a sphere).
    std::mt19937 rng(5);
    std::uniform_real_distribution<float> U(-1, 1);
    auto checkHull = [](const std::vector<Vector3>& pts, const TriMesh& h, const char* what) {
        std::map<std::pair<uint32_t, uint32_t>, int> edges;
        for (auto& t : h.triangles)
            for (int k = 0; k < 3; ++k) ++edges[{t[k], t[(k + 1) % 3]}];
        bool closed = true;
        for (auto& e : edges) closed &= e.second == 1 && edges.count({e.first.second, e.first.first}) == 1;
        float outside = 0;
        for (const Vector3& p : pts)
            for (size_t t = 0; t < h.triangles.size(); ++t)
                outside = std::max(outside, dot(h.faceNormal(t), p - h.positions[h.triangles[t][0]]));
        CHECK(closed && h.triangles.size() == 2 * h.positions.size() - 4, "%s: hull not a closed manifold (V=%zu F=%zu)", what,
              h.positions.size(), h.triangles.size());
        CHECK(outside < 1e-5f, "%s: input point %.2e outside the hull", what, outside);
        CHECK(h.signedVolume() > 0, "%s: hull inside out", what);
    };
    std::vector<Vector3> cloud;
    for (int i = 0; i < 5000; ++i) cloud.push_back({U(rng), U(rng), U(rng)});
    checkHull(cloud, buildConvexHull(cloud), "random cloud");
    std::vector<Vector3> grid;
    for (int i = 0; i <= 24; ++i)
        for (int j = 0; j <= 24; ++j)
            for (int k = 0; k <= 24; ++k) {
                Vector3 p(i - 12.0f, j - 12.0f, k - 12.0f);
                if (length(p) <= 12.0f) grid.push_back(p * 0.01f);
            }
    auto t0 = std::chrono::steady_clock::now();
    TriMesh gh = buildConvexHull(grid);
    double hullMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    checkHull(grid, gh, "grid sphere");
    TriMesh lim = buildConvexHull(grid, 32);
    CHECK(lim.positions.size() == 32, "vertex limit: %zu vertices", lim.positions.size());

    // 2) Teapot -> minimal set of convex parts fitted to the smooth surface.
    auto parts = primitives::teapotParts(0.3f);
    DecompositionParams prm;
    prm.resolution = 24;
    t0 = std::chrono::steady_clock::now();
    auto hulls = convexDecomposition(parts, prm);
    double decMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::vector<MeshBVH> sb(parts.size()), hb(hulls.size());
    for (size_t i = 0; i < parts.size(); ++i) sb[i].build(parts[i]);
    size_t verts = 0;
    for (size_t i = 0; i < hulls.size(); ++i) {
        hb[i].build(hulls[i]);
        verts += hulls[i].positions.size();
        CHECK(hulls[i].triangles.size() == 2 * hulls[i].positions.size() - 4, "part %zu not a valid hull", i);
        CHECK(int(hulls[i].positions.size()) <= prm.maxHullVertices, "part %zu has %zu vertices", i, hulls[i].positions.size());
    }
    AABB bb;
    for (auto& p : parts) bb.expand(p.bounds());
    int inSolid = 0, inHull = 0, both = 0;
    std::uniform_real_distribution<float> U01(0, 1);
    for (int k = 0; k < 40000; ++k) {
        Vector3 p = bb.lo + Vector3(U01(rng), U01(rng), U01(rng)) * bb.extent();
        bool s = false, h = false;
        for (auto& x : sb) if (x.bounds().contains(p) && x.isInside(p)) { s = true; break; }
        for (auto& x : hb) if (x.bounds().contains(p) && x.isInside(p)) { h = true; break; }
        inSolid += s, inHull += h, both += s && h;
    }
    float covered = float(both) / inSolid, excess = float(inHull - both) / inSolid;
    std::printf("  quickhull: grid sphere %zu pts -> %zu verts in %.1f ms\n", grid.size(), gh.positions.size(), hullMs);
    std::printf("  teapot: %zu convex parts, %zu vertices, %.0f ms; covers %.1f%% of the solid, excess %.1f%%\n", hulls.size(), verts,
                decMs, 100 * covered, 100 * excess);
    CHECK(hulls.size() >= 4 && hulls.size() <= 12, "part count %zu", hulls.size());
    CHECK(covered > 0.97f && excess < 0.08f, "decomposition does not fit the solid: covered %f excess %f", covered, excess);
    CHECK(decMs < 5000, "decomposition too slow: %.0f ms", decMs);

    // 3) Compound mass properties agree with the solid (volume, COM) within the fitting error.
    CompoundShape cs(hulls, primitives::merge(parts));
    float solidVol = 0;
    for (auto& p : parts) solidVol += p.signedVolume();
    CHECK(std::fabs(cs.volume() - solidVol) < 0.1f * solidVol, "compound volume %f vs solid %f", cs.volume(), solidVol);
}

void testTeapots() {
    // 100 non-convex bodies (compound of convex parts) dropped into a box: no deep penetration
    // between parts, everything comes to rest inside the domain.
    Simulation sim;
    loadSample(sim, Preset::RigidTeapots);
    CHECK(sim.rigid.bodies().size() == 100, "bodies %zu", sim.rigid.bodies().size());
    float worst = 0, ms = 0, peakMs = 0;
    for (int f = 0; f < 420; ++f) {
        auto t0 = std::chrono::steady_clock::now();
        sim.stepFrame();
        float dt = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
        ms += dt, peakMs = std::max(peakMs, dt);
        if (f % 10 == 9) worst = std::max(worst, maxPartOverlap(sim.rigid));
    }
    float vmax = 0;
    int sleeping = 0, outside = 0;
    AABB dom({-2.0f, 0.0f, -2.0f}, {2.0f, 5.0f, 2.0f}); // preset domain
    dom.lo -= Vector3(0.05f), dom.hi += Vector3(0.05f);
    for (const RigidBody& b : sim.rigid.bodies()) {
        vmax = std::max(vmax, length(b.vel));
        sleeping += b.sleeping;
        outside += !dom.contains(b.pos);
    }
    std::printf("  100 teapots, 7 s: %.1f ms/frame (peak %.0f), worst part overlap %.4f m, max |v| %.3f m/s, sleeping %d, outside %d\n",
                ms / 420, peakMs, worst, vmax, sleeping, outside);
    CHECK(outside == 0, "%d teapots left the box", outside);
    CHECK(worst < 0.02f, "teapots interpenetrate by %f m", worst);
    CHECK(vmax < 0.2f, "teapots do not come to rest: %f m/s", vmax);
}

void testBeamOverCubes() {
    // A 2 m beam falls flat and spinning onto a row of cubes: it must rest on / knock them, never
    // pass through them, at any substep.
    RigidWorld w;
    w.setDomain(AABB({-5, 0, -5}, {5, 10, 5}));
    w.params.sleeping = false;
    for (int i = 0; i < 6; ++i) w.addBox({-0.75f + 0.3f * i, 0.1f, 0}, Vector3(0.1f), Quaternion(), 500, Vector3(1));
    int beam = w.addBox({0, 1.5f, 0}, {1.0f, 0.03f, 0.03f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.3f), 800, Vector3(1));
    w.bodies()[beam].vel = {0, -25, 0};
    w.bodies()[beam].angVel = {0, 40, 0};
    float worst = 0;
    for (int f = 0; f < 120; ++f)
        for (int k = 0; k < w.params.substeps; ++k) {
            w.step(1.0f / 60 / w.params.substeps);
            worst = std::max(worst, maxOverlap(w));
        }
    std::printf("  beam 25 m/s + 40 rad/s onto 6 cubes: max overlap %.4f m, beam y=%.3f\n", worst, w.bodies()[beam].pos.y);
    CHECK(worst < 0.01f, "beam and cubes interpenetrated by %f m", worst);
    CHECK(w.bodies()[beam].pos.y > 0.02f, "beam fell through (y=%f)", w.bodies()[beam].pos.y);
}

void testConvexRest() {
    RigidWorld w;
    w.setDomain(AABB({-5, 0, -5}, {5, 10, 5}));
    int c = w.addConvex(primitives::cylinder(0.2f, 0.4f, 16), {0, 1.0f, 0}, Quaternion::fromAxisAngle({1, 0, 0}, kPi / 2), 700, Vector3(1));
    int b = w.addBox({1.0f, 0.3f, 0}, {0.3f, 0.3f, 0.3f}, Quaternion(), 500, Vector3(1));
    // Cone axis X -> Y with the apex up (rotation -90 deg about Z), base resting on the box.
    int t = w.addConvex(primitives::cone(0.2f, 0.4f, 12), {1.0f, 1.2f, 0}, Quaternion::fromAxisAngle({0, 0, 1}, -kPi / 2), 700, Vector3(1));
    for (int f = 0; f < 60 * 4; ++f)
        for (int k = 0; k < w.params.substeps; ++k) w.step(1.0f / 60 / w.params.substeps);
    const RigidBody &C = w.bodies()[c], &T = w.bodies()[t];
    std::printf("  convex: cylinder y=%.4f v=%.4f | cone on box y=%.4f v=%.4f\n", C.pos.y, length(C.vel), T.pos.y, length(T.vel));
    // Cylinder standing on its end (axis vertical after the rotation) rests at half its height.
    CHECK(std::fabs(C.worldBounds().lo.y) < 0.02f, "cylinder bottom %f", C.worldBounds().lo.y);
    CHECK(length(C.vel) < 0.05f && length(T.vel) < 0.05f, "convex bodies not at rest");
    CHECK(T.worldBounds().lo.y > 0.55f, "cone must rest on the box, bottom %f", T.worldBounds().lo.y);
    (void)b;
}

#include "HardContactTests.h" // industry-hard contact cases (uses CHECK and maxOverlap of TestRunner.h)

// Noether's theorem, one number each: energy is conserved because the laws do not change with
// time, momentum because they do not change with position, angular momentum because they do not
// change with orientation. A rigid solver that leaks any of them has a bias its impulses do not
// justify (Baumgarte and split impulses move positions, not energy; damping is off here).
static Vector3 angularMomentumAboutOrigin(const RigidWorld& w) {
    Vector3 L(0.0f);
    for (const RigidBody& b : w.bodies()) {
        if (b.invMass == 0) continue;
        const Matrix3x3 R = b.rotation();
        const Matrix3x3 I = R * Matrix3x3::diag({1 / b.invInertiaLocal.x, 1 / b.invInertiaLocal.y, 1 / b.invInertiaLocal.z}) * R.transposed();
        L += cross(b.pos, b.vel) * b.mass + I * b.angVel;
    }
    return L;
}

static Vector3 momentum(const RigidWorld& w) {
    Vector3 P(0.0f);
    for (const RigidBody& b : w.bodies()) if (b.invMass > 0) P += b.vel * b.mass;
    return P;
}

// Newton's cradle on the floor: six steel balls in a row 1 cm apart, e = 1, no friction; the first
// comes in at 2 m/s. Momentum and energy pass down the row: the last ball leaves at 2 m/s, the
// others stop. Found by the editor's example: the hit used to be perfectly inelastic (all six at
// 2/6 of the speed), because a contact entering the slop zone within the step was already pushing
// but was not counted as an impact (prepareContactPoints: the impact test now uses the solver's
// own target velocity).
void testNewtonCradle() {
    RigidWorld w;
    w.setDomain(AABB({-2, 0, -1}, {2, 2, 1}));
    w.params.sleeping = false;
    w.params.linearDamping = w.params.angularDamping = 0;
    w.params.rollingResistance = 0;
    const float r = 0.05f, gap = 0.01f;
    for (int i = 0; i < 6; ++i) {
        const int b = w.addSphere({-0.5f + i * (2 * r + gap), r, 0}, r, 7800, Vector3(1));
        w.bodies()[b].restitution = 1.0f;
        w.bodies()[b].friction = w.bodies()[b].staticFriction = 0;
    }
    w.bodies()[0].vel = {2, 0, 0};
    for (int f = 0; f < 30; ++f)
        for (int k = 0; k < w.params.substeps; ++k) w.step(1.0f / 60 / w.params.substeps);
    float restMax = 0;
    for (int i = 0; i < 5; ++i) restMax = std::max(restMax, std::fabs(w.bodies()[i].vel.x));
    const float last = w.bodies()[5].vel.x;
    std::printf("  Newton's cradle on the floor: last ball %.4f m/s (expected 2), others at most %.4f m/s\n", last, restMax);
    CHECK(std::fabs(last - 2.0f) < 0.02f, "the last ball leaves at %f m/s, expected 2", last);
    CHECK(restMax < 0.02f, "the other balls must stop, one moves at %f m/s", restMax);
}

void testNoetherRigid() {
    const float dt = 1.0f / 600;
    // 1) Time translation -> energy. An elastic ball (e = 1) dropped from 1 m: E = m g h + m v^2 / 2
    //    should stay what it was; a constant restitution below 1 would lose a fixed share per
    //    bounce, the position corrections must not add or take any.
    {
        RigidWorld w;
        w.setDomain(AABB({-2, 0, -2}, {2, 5, 2}));
        w.params.sleeping = false;
        w.params.linearDamping = w.params.angularDamping = 0;
        w.params.rollingResistance = 0;
        const float r = 0.1f;
        const int s = w.addSphere({0, 1.0f + r, 0}, r, 1000, Vector3(1));
        RigidBody& b = w.bodies()[s];
        b.restitution = 1.0f;
        b.friction = b.staticFriction = 0;
        const float m = b.mass, g = -w.params.gravity.y;
        auto energy = [&] { return m * g * (b.pos.y - r) + 0.5f * m * length2(b.vel); };
        const float E0 = energy();
        int bounces = 0;
        float prevVy = 0, worst = 0;
        std::printf("  elastic ball: E0 = %.4f J;", E0);
        for (int k = 0; k < 600 * 20 && bounces < 20; ++k) { // 20 bounces of a 1 m drop take ~18 s
            w.step(dt);
            if (prevVy < 0 && b.vel.y > 0) { // a bounce: energy right after it
                ++bounces;
                worst = std::max(worst, std::fabs(energy() / E0 - 1));
                if (bounces % 5 == 0) std::printf(" after %d bounces %.4f J (%+.2f%%)", bounces, energy(), 100 * (energy() / E0 - 1));
            }
            prevVy = b.vel.y;
        }
        std::printf("; worst drift %.2f%% over %d bounces\n", 100 * worst, bounces);
        // Found by this test: the ball dropped dead. An impact arrives as a speculative contact (the
        // gap is smaller than contactMargin one step before touching), and that branch applied the
        // restitution only to CCD-clamped bodies; a ball at a few m/s is never clamped, so every
        // ordinary bounce was inelastic. The bounce now lives in RigidWorld::applyRestitution().
        CHECK(bounces >= 20, "the ball stopped bouncing after %d bounces", bounces);
        CHECK(worst < 0.02f, "energy of an elastic ball drifts by %.2f%%", 100 * worst);
    }
    // 2) Space translation -> momentum; rotation -> angular momentum about a fixed point. Two boxes
    //    in zero gravity, a glancing hit that spins them both: the internal impulses (normal and
    //    friction) cancel pairwise, so P and L about the origin stay what they were.
    //    The same hit at two time steps: the contact impulses cancel exactly (same point, opposite
    //    sign), so what drift is left is the discretisation of the free rotation of the spun boxes
    //    (it grows like dt w^2: the bouncing hit spins them to 5-8 rad/s). It must shrink with dt.
    auto glancing = [&](int substeps, float& worstP, float& worstL, float& spinC) {
        RigidWorld w;
        w.setDomain(AABB({-10, -10, -10}, {10, 10, 10}));
        w.params.gravity = Vector3(0.0f);
        w.params.sleeping = false;
        w.params.linearDamping = w.params.angularDamping = 0;
        w.params.rollingResistance = 0;
        w.params.collideWithDomain = false;
        const int a = w.addBox({-1, 0, 0}, {0.15f, 0.1f, 0.2f}, Quaternion(), 800, Vector3(1));
        const int c = w.addBox({1, 0.12f, 0.05f}, {0.2f, 0.15f, 0.1f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.4f), 500, Vector3(1));
        w.bodies()[a].vel = {3, 0, 0};
        w.bodies()[a].angVel = {0, 0, 2};
        w.bodies()[c].vel = {-1, 0, 0};
        const Vector3 P0 = momentum(w), L0 = angularMomentumAboutOrigin(w);
        worstP = worstL = 0;
        for (int k = 0; k < 600 * substeps; ++k) {
            w.step(dt / float(substeps));
            worstP = std::max(worstP, length(momentum(w) - P0) / length(P0));
            worstL = std::max(worstL, length(angularMomentumAboutOrigin(w) - L0) / length(L0));
        }
        spinC = length(w.bodies()[c].angVel);
    };
    {
        float p1, l1, s1, p2, l2, s2;
        glancing(1, p1, l1, s1);
        glancing(2, p2, l2, s2);
        std::printf("  glancing collision: worst |P| drift %.1e / %.1e, worst |L| drift %.1e at dt = 1/600, %.1e at 1/1200 (spin of the hit box %.2f rad/s)\n",
                    p1, p2, l1, l2, s1);
        CHECK(p1 < 1e-4f && p2 < 1e-4f, "momentum drifts %e / %e", p1, p2);
        CHECK(l2 < l1, "angular momentum drift must shrink with the step: %e at dt, %e at dt/2", l1, l2);
        CHECK(l2 < 2e-3f, "angular momentum drifts %e at dt/2", l2);
        CHECK(s1 > 0.3f, "the glancing hit must spin the second box (%f rad/s)", s1);
    }
    // 3) Rotation symmetry of a free body: a box spun about its middle axis of inertia flips over
    //    and over (Dzhanibekov / tennis racket) while |L| and the rotational energy stay constant.
    //    The perturbation grows at sigma = w2 sqrt((I2 - I1)(I3 - I2) / (I1 I3)) (linearised Euler
    //    equations); the flips repeat with a period of the order of a few 1/sigma.
    {
        RigidWorld w;
        w.setDomain(AABB({-10, -10, -10}, {10, 10, 10}));
        w.params.gravity = Vector3(0.0f);
        w.params.sleeping = false;
        w.params.linearDamping = w.params.angularDamping = 0;
        w.params.collideWithDomain = false;
        const int i = w.addBox({0, 0, 0}, {0.05f, 0.2f, 0.1f}, Quaternion(), 1000, Vector3(1)); // I_x < I_z < I_y: middle axis z
        RigidBody& b = w.bodies()[i];
        const Vector3 I(1 / b.invInertiaLocal.x, 1 / b.invInertiaLocal.y, 1 / b.invInertiaLocal.z);
        const float w2 = 10.0f;
        b.angVel = {0.02f * w2, 0, w2}; // about the middle axis (z), nudged
        const float sigma = w2 * std::sqrt((I.z - I.x) * (I.y - I.z) / (I.x * I.y));
        const float L0 = length(angularMomentumAboutOrigin(w)), E0 = w.kineticEnergy();
        float worstL = 0, worstE = 0, prevWz = b.angVel.z, firstFlip = -1, lastFlip = -1;
        int flips = 0;
        for (int k = 0; k < 6000; ++k) {
            w.step(dt);
            worstL = std::max(worstL, std::fabs(length(angularMomentumAboutOrigin(w)) / L0 - 1));
            worstE = std::max(worstE, std::fabs(w.kineticEnergy() / E0 - 1));
            const float wz = dot(b.rotation() * Vector3(0, 0, 1), b.angVel); // spin about the body's middle axis
            if (prevWz * wz < 0) { ++flips; (firstFlip < 0 ? firstFlip : lastFlip) = k * dt; lastFlip = k * dt; }
            prevWz = wz;
        }
        const float period = flips > 1 ? 2 * (lastFlip - firstFlip) / float(flips - 1) : 0;
        std::printf("  Dzhanibekov: I = (%.4f %.4f %.4f) kg m^2, sigma %.2f 1/s, %d flips in 10 s, period %.2f s (%.1f / sigma); |L| drift %.1e, "
                    "energy drift %.1e\n", I.x, I.y, I.z, sigma, flips, period, period * sigma, worstL, worstE);
        // Found by this test: no flip. The sequential-impulse integrator (RigidWorld::step) advances
        // the angular velocity without the gyroscopic term w x (I w) of Euler's equations - only
        // the experimental XPBD path has it (XpbdSolver.cpp) - so a free body keeps w instead of L.
        // The same omission is what lets L drift in the tumbling boxes of part 2.
        CHECK(flips >= 2, "no tennis-racket flip: %d sign changes", flips);
        CHECK(worstL < 1e-3f, "|L| of a free body drifts %e", worstL);
        CHECK(worstE < 1e-2f, "rotational energy of a free body drifts %e", worstE);
    }
}

// Removing one body (the editor's meta-objects): a stack of five boxes and a ball rolling beside
// it; the ball is destroyed mid-run, then the top box of the stack. What stays keeps standing,
// the slot of a destroyed body is reused by the next add, and add/destroy cycles allocate
// nothing that grows.
static void stepFrames(RigidWorld& w, int frames) {
    for (int f = 0; f < frames; ++f)
        for (int k = 0; k < w.params.substeps; ++k) w.step(1.0f / 60 / w.params.substeps);
}

static float stackDrift(const RigidWorld& w, const std::vector<int>& boxes, const std::vector<Vector3>& start) {
    float drift = 0;
    for (size_t k = 0; k < boxes.size(); ++k)
        if (w.isAlive(boxes[k])) drift = std::max(drift, length(w.bodies()[size_t(boxes[k])].pos - start[k]));
    return drift;
}

void testDestroyBody() {
    RigidWorld w;
    w.setDomain(AABB({-3, 0, -3}, {3, 4, 3}));
    std::vector<int> boxes;
    for (int i = 0; i < 5; ++i) boxes.push_back(w.addBox({0, 0.1f + 0.2f * i, 0}, Vector3(0.1f), Quaternion(), 500, Vector3(1)));
    const int ball = w.addSphere({-1.5f, 0.1f, 0.8f}, 0.1f, 500, Vector3(1));
    w.bodies()[size_t(ball)].vel = {1.0f, 0, 0};
    stepFrames(w, 60);
    std::vector<Vector3> start;
    for (int b : boxes) start.push_back(w.bodies()[size_t(b)].pos);
    const size_t contactsBefore = w.contactCount();
    w.destroyBody(ball);
    stepFrames(w, 60);
    const float drift1 = stackDrift(w, boxes, start);
    std::printf("  destroy the ball: bodies %d, contacts %zu -> %zu, stack drift %.2e m\n", w.bodyCount(), contactsBefore,
                w.contactCount(), drift1);
    CHECK(w.bodyCount() == 5 && !w.isAlive(ball), "the ball must be gone (%d bodies)", w.bodyCount());
    CHECK(w.contactCount() < contactsBefore, "the ball's floor contact must go: %zu -> %zu", contactsBefore, w.contactCount());
    CHECK(drift1 < 1e-3f, "the stack moved %f m when the ball went", drift1);
    w.destroyBody(boxes.back()); // the top box: the rest must keep standing
    stepFrames(w, 120);
    const float drift2 = stackDrift(w, boxes, start);
    const int again = w.addSphere({1.0f, 0.5f, 0}, 0.1f, 500, Vector3(1));
    stepFrames(w, 120);
    const RigidBody& reused = w.bodies()[size_t(again)];
    std::printf("  destroy the top box: stack drift %.2e m; the next body takes slot %d (the top box was %d), rests at y %.4f\n",
                drift2, again, boxes.back(), reused.pos.y);
    CHECK(drift2 < 1e-3f, "the stack moved %f m when its top box went", drift2);
    CHECK(again == boxes.back(), "a new body must reuse the freed slot %d, got %d", boxes.back(), again);
    CHECK(std::fabs(reused.pos.y - 0.1f) < 0.01f && std::isfinite(reused.pos.x), "the new body must land normally (y %f)", reused.pos.y);
    // 100 add / destroy cycles: the number of slots does not grow, a step allocates as before.
    const size_t slots = w.bodies().size();
    const long long b0 = Probe::allocations.load();
    stepFrames(w, 10);
    const long long perFrameBefore = (Probe::allocations.load() - b0) / 10;
    for (int c = 0; c < 100; ++c) {
        const int t = w.addBox({-1.0f, 0.3f, -1.0f}, Vector3(0.05f), Quaternion(), 500, Vector3(1));
        stepFrames(w, 1);
        w.destroyBody(t);
    }
    const long long a0 = Probe::allocations.load();
    stepFrames(w, 10);
    const long long perFrame = (Probe::allocations.load() - a0) / 10;
    std::printf("  100 add/destroy cycles: %zu slots before, %zu after; %lld allocations per frame before, %lld after\n", slots,
                w.bodies().size(), perFrameBefore, perFrame);
    CHECK(w.bodies().size() == slots, "slots grew from %zu to %zu", slots, w.bodies().size());
    CHECK(perFrame <= perFrameBefore + 2, "a frame allocates %lld times after the cycles, %lld before", perFrame, perFrameBefore);
}
