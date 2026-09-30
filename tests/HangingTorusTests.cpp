// Exercise the short editor chain with the real mouse constraint: pick a rim, pull in 3D,
// abruptly move the target and release. Linking is observed, never enforced by a joint.
#include "TestRunner.h"
#include "Tests.h"
#include "samples/HangingTorusScene.h"
#include "samples/TorusChainsScene.h"
#include "core/Format.h"

#include <filesystem>
#include <fstream>

namespace {

Vector3 mouseTarget(const Vector3& start, int frame) {
    if (frame >= 150 && frame < 180) return start + Vector3(-0.22f, 0.18f, -0.08f);
    const float phase = 2 * kPi * float(frame - 60) / 72;
    return start + Vector3(0.18f * rf::sin(phase), 0.12f * (1 - rf::cos(phase)), 0.08f * rf::sin(0.5f * phase));
}

Vector3 rimCentre(const RigidBody& body) {
    const auto& shape = static_cast<const CompoundShape&>(*body.shape);
    const Matrix3x3 rotation = body.rotation() * shape.principalRotation().transposed();
    // Aim through material, not the hole at the centre of the torus.
    return body.pos + rotation * (Vector3(0, 0, TorusChainsScene::majorRadius) - shape.centerOfMass());
}

std::ofstream mouseCsv() {
    if (const char* dir = std::getenv("RF_PLOT_DIR"); dir && *dir) {
        std::filesystem::create_directories(dir);
        std::ofstream file(std::filesystem::path(dir) / "hanging-torus-mouse.csv");
        CHECK(bool(file), "cannot write mouse chain CSV");
        file << "frame,time_s,step_ms,broken_pairs,grabbed,tip_x,tip_y,tip_z,sampled_peak_overlap_m\n";
        return file;
    }
    return {};
}

} // namespace

void testHangingTorusMouse() {
    Simulation sim;
    loadSample(sim, Preset::HangingTorus);
    const auto* scene = static_cast<const HangingTorusScene*>(sim.scene());
    const auto& bodies = sim.rigid.bodies();
    CHECK(bodies.size() == 8 && sim.rigid.joints().empty() && !sim.rigid.params.sleeping, "short chain needs eight contact-only tori");
    int fixed = 0;
    for (const auto& b : bodies) fixed += b.invMass == 0;
    CHECK(fixed == 1 && bodies.front().invMass == 0, "only the upper torus may be fixed");
    const Vector3 fixedAt = bodies.front().pos;
    for (size_t k = 1; k < bodies.size(); ++k)
        CHECK(std::abs(TorusChainsScene::linkingNumber(bodies[k - 1], bodies[k])) == 1, "initial pair %zu is not linked", k);
    const int tip = int(bodies.size()) - 1;
    Vector3 grabbedAt;
    float peakOverlap = maxPartOverlap(sim.rigid), excursion = 0;
    std::vector<double> timings;
    std::ofstream csv = mouseCsv();
    for (int frame = 0; frame < 360; ++frame) {
        if (frame == 60) {
            const Vector3 direction = normalize(Vector3(1, 0, 1)), origin = rimCentre(bodies[size_t(tip)]) + 2 * direction;
            int hit = -1; float distance = 0; Vector3 normal;
            CHECK(sim.rigid.raycast(origin, -direction, 3, hit, distance, normal) && hit == tip, "the lower rim must be pickable");
            grabbedAt = origin - direction * distance;
            sim.rigid.grab(hit, grabbedAt);
            CHECK(sim.rigid.grabJoint().active, "mouse constraint did not activate");
        }
        if (frame >= 60 && frame < 240) sim.rigid.setGrabTarget(mouseTarget(grabbedAt, frame));
        if (frame == 240) sim.rigid.releaseGrab();
        const auto start = std::chrono::steady_clock::now();
        sim.stepFrame();
        timings.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count());
        CHECK(scene->brokenPairs() == 0, "mouse pull broke a link at frame %d", frame + 1);
        CHECK(length2(bodies.front().pos - fixedAt) < 1e-14f, "fixed support moved");
        for (const auto& b : bodies)
            CHECK(std::isfinite(length2(b.pos) + length2(b.vel) + length2(b.angVel) + b.rot.w * b.rot.w +
                               b.rot.x * b.rot.x + b.rot.y * b.rot.y + b.rot.z * b.rot.z), "non-finite chain state");
        if (frame % 15 == 0) peakOverlap = std::max(peakOverlap, maxPartOverlap(sim.rigid));
        const Vector3 at = bodies[size_t(tip)].pos;
        excursion = std::max(excursion, std::sqrt(at.x * at.x + at.z * at.z));
        if (csv) csv << frame + 1 << ',' << numberText(double(frame + 1) / 60) << ',' << numberText(timings.back()) << ','
                     << scene->brokenPairs() << ',' << int(sim.rigid.grabJoint().active) << ',' << numberText(at.x) << ','
                     << numberText(at.y) << ',' << numberText(at.z) << ',' << numberText(peakOverlap) << '\n';
    }
    std::sort(timings.begin(), timings.end());
    std::printf("  eight tori, mouse pull/release, 6 s: first break %d, tip excursion %.1f mm, sampled peak overlap %.3f mm; step median/p95 %.3f/%.3f ms\n",
        scene->firstBrokenFrame(), 1000 * excursion, 1000 * peakOverlap, timings[180], timings[341]);
    CHECK(scene->firstBrokenFrame() < 0 && excursion > 0.1f, "chain did not survive an effective mouse pull");
    CHECK(peakOverlap < 0.002f, "mouse-pull penetration %.3f mm exceeds 2 mm", 1000 * peakOverlap);
    CHECK(!sim.rigid.grabJoint().active, "mouse release left the constraint active");
}

// Sleep is deliberately off in the diagnostic sample. When explicitly enabled, the whole
// contact island must sleep, stay motionless, wake on a rim grab and settle after another pull.
void testHangingTorusSleep() {
    Simulation sim;
    loadSample(sim, Preset::HangingTorus);
    CHECK(!sim.rigid.params.sleeping, "diagnostic sample must expose residual motion by default");
    sim.rigid.params.sleeping = true;
    const auto* scene = static_cast<const HangingTorusScene*>(sim.scene());
    const auto& bodies = sim.rigid.bodies();
    const Vector3 fixedAt = bodies.front().pos;
    Vector3 grabbedAt, asleepTip;
    int firstSleep = -1, secondSleep = -1;
    float excursion = 0;
    double sleepingEnergy = 0;
    for (int frame = 0; frame < 3000; ++frame) {
        const int gesture = frame < 1500 ? frame : frame - 1500;
        if (gesture == 60) {
            if (frame > 1500) CHECK(sim.rigid.sleepingCount() == 7, "chain did not stay asleep before the second grab");
            const Vector3 direction = normalize(Vector3(1, 0, 1)), origin = rimCentre(bodies.back()) + 2 * direction;
            int hit = -1; float distance = 0; Vector3 normal;
            CHECK(sim.rigid.raycast(origin, -direction, 3, hit, distance, normal) && hit == 7, "sleeping rim must be pickable");
            grabbedAt = origin - direction * distance;
            sim.rigid.grab(hit, grabbedAt);
            CHECK(sim.rigid.sleepingCount() == 0 && sim.rigid.grabJoint().active, "grab did not wake the entire chain");
        }
        if (gesture >= 60 && gesture < 240) sim.rigid.setGrabTarget(mouseTarget(grabbedAt, gesture));
        if (gesture == 240) sim.rigid.releaseGrab();
        sim.stepFrame();
        CHECK(scene->brokenPairs() == 0 && length2(bodies.front().pos - fixedAt) < 1e-14f, "sleep or wake changed the chain topology/support");
        if (frame >= 240 && frame < 1500 && sim.rigid.sleepingCount() == 7 && firstSleep < 0) {
            firstSleep = frame + 1;
            asleepTip = bodies.back().pos;
        }
        if (frame >= 1560 && frame < 1740) excursion = std::max(excursion, length(bodies.back().pos - asleepTip));
        if (frame >= 1740 && sim.rigid.sleepingCount() == 7 && secondSleep < 0) secondSleep = frame + 1;
        if (frame == 1440) for (const auto& b : bodies) sleepingEnergy -= double(b.mass) * dot(sim.rigid.params.gravity, b.pos);
        if (frame > 1440 && frame < 1500) {
            double energy = sim.rigid.kineticEnergy();
            for (const auto& b : bodies) energy -= double(b.mass) * dot(sim.rigid.params.gravity, b.pos);
            CHECK(sim.rigid.sleepingCount() == 7 && sim.rigid.kineticEnergy() == 0 && energy == sleepingEnergy,
                  "sleeping island continues moving or changing energy");
        }
    }
    std::printf("  seven rings sleep %.3f / %.3f s after release; re-grab wakes all, tip moves %.1f mm\n",
        float(firstSleep - 240) / 60, float(secondSleep - 1740) / 60, 1000 * excursion);
    CHECK(firstSleep > 240 && firstSleep < 1440 && secondSleep > 1740 && sim.rigid.sleepingCount() == 7,
          "chain failed to sleep after both gestures");
    CHECK(excursion > 0.1f, "re-grab did not move the previously sleeping chain");
}
