// Profile unchanged short-chain physics and the snapshot publication path separately.
// Same deterministic mouse gesture as HangingTorusTests; timing excludes the overlap audit.
#include "tests/TestRunner.h"
#include "samples/HangingTorusScene.h"
#include "samples/TorusChainsScene.h"
#include "core/Probe.h"
#include "core/Format.h"
#include <iostream>

double millis(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
}

Vector3 target(const Vector3& start, int frame) {
    if (frame >= 150 && frame < 180) return start + Vector3(-0.22f, 0.18f, -0.08f);
    const float phase = 2 * kPi * float(frame - 60) / 72;
    return start + Vector3(0.18f * rf::sin(phase), 0.12f * (1 - rf::cos(phase)), 0.08f * rf::sin(0.5f * phase));
}

Vector3 rim(const RigidBody& b) {
    const auto& shape = static_cast<const CompoundShape&>(*b.shape);
    return b.pos + b.rotation() * shape.principalRotation().transposed()
        * (Vector3(0, 0, TorusChainsScene::majorRadius) - shape.centerOfMass());
}

void report(const char* name, std::vector<double> values) {
    double sum = 0;
    for (double x : values) sum += x;
    std::sort(values.begin(), values.end());
    std::cout << '"' << name << "\":{\"mean_ms\":" << numberText(sum / values.size())
        << ",\"median_ms\":" << numberText(values[values.size()/2])
        << ",\"p95_ms\":" << numberText(values[size_t(0.95 * values.size())]) << '}';
}

int main() {
    Simulation sim;
    loadSample(sim, Preset::HangingTorus);
    const auto* scene = static_cast<const HangingTorusScene*>(sim.scene());
    std::map<std::string, std::vector<double>> stages;
    std::vector<double> physics, snapshot;
    StateHash hash;
    Vector3 anchor;
    float overlap = 0;
    int broken = 0, slowFrames = 0;
    for (int frame = 0; frame < 360; ++frame) {
        if (frame == 60) {
            const Vector3 direction = normalize(Vector3(1, 0, 1)), origin = rim(sim.rigid.bodies().back()) + 2 * direction;
            int hit = -1; float distance = 0; Vector3 normal;
            if (!sim.rigid.raycast(origin, -direction, 3, hit, distance, normal) || hit != 7) return 2;
            anchor = origin - direction * distance;
            sim.rigid.grab(hit, anchor);
        }
        if (frame >= 60 && frame < 240) sim.rigid.setGrabTarget(target(anchor, frame));
        if (frame == 240) sim.rigid.releaseGrab();
        const auto start = std::chrono::steady_clock::now();
        sim.stepFrame();
        physics.push_back(millis(start));
        const auto publish = std::chrono::steady_clock::now();
        RenderSnapshot snap;
        sim.fillSnapshot(snap);
        snapshot.push_back(millis(publish));
        slowFrames += physics.back() + snapshot.back() > 1000.0 / 60;
        for (const auto& c : snap.probe.channels)
            if (c.kind == Probe::Kind::TimerMs) stages[c.name].push_back(c.value);
        for (const auto& b : sim.rigid.bodies()) { hash.add(b.pos); hash.add(b.rot); hash.add(b.vel); hash.add(b.angVel); }
        broken = std::max(broken, scene->brokenPairs());
        if (frame % 15 == 0) overlap = std::max(overlap, maxPartOverlap(sim.rigid));
    }
    std::cout << "{\"threads\":" << ThreadPool::instance().threadCount()
        << ",\"frames\":360,\"over_budget_frames\":" << slowFrames
        << ",\"broken_pairs\":" << broken << ",\"sampled_overlap_m\":" << numberText(overlap)
        << ",\"trajectory_hash\":\"" << std::hex << hash.h << std::dec << "\",";
    report("physics", physics); std::cout << ','; report("snapshot", snapshot);
    for (const auto& item : stages) { std::cout << ','; report(item.first.c_str(), item.second); }
    std::cout << "}\n";
    return broken != 0 || overlap > 0.002f;
}
