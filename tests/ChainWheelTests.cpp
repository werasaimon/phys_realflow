// The stationary chain wheel: valid assembly, gravity settling and an actual mouse pull/release.
// Sampled GJK/EPA depths are regression checks; the separate TorusGeometryAudit tool uses SAT.
#include "TestRunner.h"
#include "Tests.h"
#include "samples/ChainWheelScene.h"
#include "samples/TorusChainsScene.h"
#include "core/Format.h"

#include <filesystem>
#include <fstream>

namespace {

double wheelEnergy(const Simulation& sim) {
    double result = sim.rigid.kineticEnergy();
    for (const auto& b : sim.rigid.bodies())
        if (b.invMass > 0) result -= double(b.mass) * dot(sim.rigid.params.gravity, b.pos);
    return result;
}

std::ofstream wheelCsv(const char* name, const char* header) {
    if (const char* dir = std::getenv("RF_PLOT_DIR"); dir && *dir) {
        std::filesystem::create_directories(dir);
        std::ofstream out(std::filesystem::path(dir) / name);
        CHECK(bool(out), "cannot open chain-wheel output %s", name);
        out << header << '\n';
        return out;
    }
    return {};
}

void wheelPoses(std::ofstream& out, int frame, const Simulation& sim) {
    if (!out) return;
    const auto& bodies = sim.rigid.bodies();
    for (size_t i = 0; i < bodies.size(); ++i) {
        const auto& b = bodies[i];
        out << frame << ',' << i;
        for (float x : {b.pos.x, b.pos.y, b.pos.z, b.rot.w, b.rot.x, b.rot.y, b.rot.z,
                        b.vel.x, b.vel.y, b.vel.z, b.angVel.x, b.angVel.y, b.angVel.z}) out << ',' << numberText(x);
        out << '\n';
    }
}

Vector3 wheelGrab(Simulation& sim) {
    const auto& body = sim.rigid.bodies().back();
    const auto& shape = static_cast<const CompoundShape&>(*body.shape);
    const Vector3 rim = body.pos + body.rotation() * shape.principalRotation().transposed()
        * (Vector3(0, 0, TorusChainsScene::majorRadius) - shape.centerOfMass());
    const Vector3 direction = normalize(Vector3(1, 0, 1)), origin = rim + 2 * direction;
    int hit = -1; float distance = 0; Vector3 normal;
    CHECK(sim.rigid.raycast(origin, -direction, 3, hit, distance, normal) && hit == ChainWheelScene::lastLink,
          "loaded ring must be selectable, hit %d", hit);
    const Vector3 at = origin - direction * distance;
    if (hit == ChainWheelScene::lastLink) sim.rigid.grab(hit, at);
    return at;
}

void wheelState(const Simulation& sim, const Vector3& wheelAt, const Vector3& anchorAt, int frame) {
    const auto* scene = static_cast<const ChainWheelScene*>(sim.scene());
    CHECK(scene->firstBrokenFrame() < 0, "chain-wheel link lost at frame %d (now %d)", scene->firstBrokenFrame(), frame);
    const auto& b = sim.rigid.bodies();
    CHECK(length2(b[0].pos - wheelAt) == 0 && length2(b[1].pos - anchorAt) == 0, "stationary support moved");
    for (const auto& body : b) {
        CHECK(std::isfinite(length2(body.pos) + length2(body.vel) + length2(body.angVel)), "non-finite body at %d", frame);
        CHECK(body.worldBounds().lo.y > 0.05f, "floor is supporting the hanging chain at %d", frame);
    }
}

} // namespace

void testChainWheelAssembly() {
    Simulation sim;
    loadSample(sim, Preset::ChainWheel);
    const auto& b = sim.rigid.bodies();
    CHECK(b.size() == 17 && sim.rigid.joints().empty(), "wheel plus sixteen rings, no hidden link joints");
    int fixed = 0;
    for (const auto& body : b) fixed += body.invMass == 0;
    CHECK(fixed == 2 && b[0].invMass == 0 && b[1].invMass == 0 && !sim.rigid.params.sleeping,
          "wheel and return end fixed; diagnostics run without sleep");
    for (int i = 1; i < ChainWheelScene::linkCount; ++i) {
        CHECK(std::abs(TorusChainsScene::linkingNumber(b[size_t(i)], b[size_t(i + 1)])) == 1, "initial link %d", i);
        CHECK(std::fabs(length(b[size_t(i)].pos - b[size_t(i + 1)].pos) - TorusChainsScene::spacing) < 2e-6f,
              "analytic chord/straight pitch changed at link %d", i);
    }
    CHECK(std::fabs(b.back().mass / b[2].mass - 2) < 1e-5f, "loaded ring must have twice the mass");
    const float depth = maxPartOverlap(sim.rigid);
    CHECK(depth < 2e-6f, "initial penetration %.6f m", depth);
    std::printf("  16 rings, wheel radius %.6f m, initial sampled depth %.3f mm, 15 links\n",
                ChainWheelScene::pitchRadius(), 1000 * depth);
    sim.scene()->setParam(0, 3);
    sim.reset();
    CHECK(std::fabs(sim.rigid.bodies().back().mass / sim.rigid.bodies()[2].mass - 3) < 1e-5f, "load reset lost parameter");
}

void testChainWheelLoaded() {
    Simulation sim;
    loadSample(sim, Preset::ChainWheel);
    auto csv = wheelCsv("chain-wheel.csv", "frame,time_s,step_ms,energy_J,kinetic_J,broken_pairs,sampled_depth_m,grabbed");
    auto poses = wheelCsv("chain-wheel-poses.csv", "frame,body,x,y,z,qw,qx,qy,qz,vx,vy,vz,wx,wy,wz");
    const Vector3 wheelAt = sim.rigid.bodies()[0].pos, anchorAt = sim.rigid.bodies()[1].pos;
    const double initialEnergy = wheelEnergy(sim);
    double passivePeak = initialEnergy, lateKinetic = 0;
    float overlap = 0, excursion = 0;
    int supportedFrames = 0;
    Vector3 grabbedAt, tipAt;
    std::vector<double> timings;
    wheelPoses(poses, 0, sim);
    for (int frame = 0; frame < 720; ++frame) {
        if (frame == 180) { grabbedAt = wheelGrab(sim); tipAt = sim.rigid.bodies().back().pos; }
        if (frame >= 180 && frame < 300) {
            const float phase = 2 * kPi * float(frame - 180) / 120;
            sim.rigid.setGrabTarget(grabbedAt + Vector3(0.08f * rf::sin(phase), 0.06f * (1 - rf::cos(phase)), 0.025f * rf::sin(phase)));
        }
        if (frame == 300) sim.rigid.releaseGrab();
        const auto start = std::chrono::steady_clock::now();
        sim.stepFrame();
        timings.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count());
        wheelState(sim, wheelAt, anchorAt, frame + 1);
        if (frame >= 60 && frame < 180) {
            bool supported = false;
            for (const auto& c : sim.rigid.debugContacts())
                supported |= (c.a == 0 || c.b == 0) && c.impulse > 0;
            supportedFrames += supported;
        }
        if (frame < 180) passivePeak = std::max(passivePeak, wheelEnergy(sim));
        if (frame >= 180 && frame < 300) excursion = std::max(excursion, length(sim.rigid.bodies().back().pos - tipAt));
        if (frame >= 660) lateKinetic += sim.rigid.kineticEnergy() / 60.0;
        if (frame % 15 == 14) overlap = std::max(overlap, maxPartOverlap(sim.rigid));
        if (frame % 60 == 59 || frame == 300) wheelPoses(poses, frame + 1, sim);
        if (csv) csv << frame + 1 << ',' << numberText(double(frame + 1) / 60) << ',' << numberText(timings.back()) << ','
            << numberText(wheelEnergy(sim)) << ',' << numberText(sim.rigid.kineticEnergy()) << ','
            << static_cast<const ChainWheelScene*>(sim.scene())->brokenPairs() << ',' << numberText(overlap) << ','
            << int(sim.rigid.grabJoint().active) << '\n';
    }
    std::sort(timings.begin(), timings.end());
    std::printf("  wheel: 12 s, excursion %.2f mm, sampled overlap %.3f mm, passive gain %.6g J, late T %.6g J, step median/p95 %.3f/%.3f ms\n",
                1000 * excursion, 1000 * overlap, passivePeak - initialEnergy, lateKinetic, timings[360], timings[684]);
    CHECK(passivePeak <= initialEnergy + 0.01, "unforced energy increase %.6f J", passivePeak - initialEnergy);
    CHECK(overlap < 0.0005f, "sampled penetration %.3f mm exceeds 0.5 mm", 1000 * overlap);
    CHECK(excursion > 0.05f && !sim.rigid.grabJoint().active, "ineffective pull or unreleased mouse");
    CHECK(lateKinetic < 0.001, "residual kinetic energy %.6g J exceeds 1 mJ", lateKinetic);
    CHECK(supportedFrames > 100, "chain must load the wheel, supported in %d/120 frames", supportedFrames);
}

void testChainDriveCycle() {
    Simulation sim;
    loadSample(sim, Preset::ChainDrive);
    CHECK(sim.rigid.bodies().size() == 17 && sim.rigid.joints().size() == 1, "drive needs 17 bodies and one axle hinge");
    for (const auto& b : sim.rigid.bodies()) CHECK(b.invMass > 0, "all drive bodies must be dynamic");
    auto& motor = static_cast<HingeJoint&>(*sim.rigid.joints().front());
    CHECK(motor.a == 0 && motor.b == -1 && motor.motorEnabled && motor.maxMotorTorque == 20,
          "wheel must have the bounded physical hinge motor");
    auto csv = wheelCsv("chain-drive.csv", "frame,time_s,angle_rad,speed_rad_s,command_rad_s,tip_y_m,energy_J,sampled_depth_m");
    auto poses = wheelCsv("chain-drive-poses.csv", "frame,body,x,y,z,qw,qx,qy,qz,vx,vy,vz,wx,wy,wz");
    const Vector3 centre = sim.rigid.bodies()[0].pos;
    float overlap = maxPartOverlap(sim.rigid), minimumLift = 1e9f, maximumReturn = 0, holdSpeed = 0;
    float baseY = 0, baseAngle = 0;
    wheelPoses(poses, 0, sim);
    CHECK(overlap < 2e-6f, "drive starts with overlap %.6f m", overlap);
    for (int frame = 0; frame < 1200; ++frame) {
        const float command = motor.motorSpeed;
        sim.stepFrame();
        const auto& bodies = sim.rigid.bodies();
        const auto* scene = static_cast<const ChainWheelScene*>(sim.scene());
        const float angle = motor.angle(bodies), y = bodies.back().pos.y, speed = bodies[0].angVel.z;
        CHECK(scene->firstBrokenFrame() < 0, "driven chain broke at %d", scene->firstBrokenFrame());
        CHECK(length(bodies[0].pos - centre) < 0.0003f, "axle drift at frame %d", frame + 1);
        for (const auto& b : bodies)
            CHECK(std::isfinite(length2(b.pos) + length2(b.vel) + length2(b.angVel)) && b.worldBounds().lo.y > 0.05f,
                  "drive must stay finite and off the floor, frame %d", frame + 1);
        if (frame % 15 == 14) overlap = std::max(overlap, maxPartOverlap(sim.rigid));
        const int cycleFrame = (frame + 1) % 600;
        if (cycleFrame == 120) { baseY = y; baseAngle = angle; }
        if (cycleFrame == 240) {
            minimumLift = std::min(minimumLift, y - baseY);
            // A tracked 0.15 rad/s command over two seconds gives 0.30 rad, allowing startup lag.
            CHECK(std::fabs(angle - baseAngle - 0.30f) < 0.01f && speed > 0.14f, "motor failed to track the two-second lift");
        }
        if (cycleFrame == 360) holdSpeed = std::max(holdSpeed, std::fabs(speed));
        if (cycleFrame == 480) {
            maximumReturn = std::max(maximumReturn, std::fabs(y - baseY));
            CHECK(std::fabs(angle - baseAngle) < 0.02f && speed < -0.14f, "motor failed to reverse");
        }
        if (frame % 60 == 59) wheelPoses(poses, frame + 1, sim);
        if (csv) csv << frame + 1 << ',' << numberText(double(frame + 1) / 60) << ',' << numberText(angle) << ','
            << numberText(speed) << ',' << numberText(command) << ',' << numberText(y) << ','
            << numberText(wheelEnergy(sim)) << ',' << numberText(overlap) << '\n';
    }
    std::printf("  driven chain: 20 s / 2 cycles, minimum lift %.2f mm, return error %.3f mm, hold speed %.6f rad/s, sampled depth %.3f mm\n",
                1000 * minimumLift, 1000 * maximumReturn, holdSpeed, 1000 * overlap);
    CHECK(minimumLift > 0.05f && maximumReturn < 0.005f, "chain did not follow the forward/reverse drive");
    CHECK(holdSpeed < 0.001f && overlap < 0.0005f, "drive hold or contact tolerance failed");
}
