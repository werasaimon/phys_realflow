// Adapted reproduction of SimBenchmark Bouncing, with shared observers for SDK and Bullet.
// Scenario and energy metric: https://leggedrobotics.github.io/SimBenchmark/bouncing/
// This uses btDiscreteDynamicsWorld rather than the author's historical multibody wrapper.
// Build: cmake --build REFERENCE_BUILD --target rf_simbenchmark_bouncing
// Usage: rf_simbenchmark_bouncing ENGINE DT NEW_DIRECTORY [iterations=20]
// ENGINE: sdk, sdk-ccd, bullet; optional builds also jolt, box3d, box3d-native.
#include "core/Format.h"
#include "rigid_reference/BouncingSdk.h"
#include "rigid_reference/BouncingVariational.h"
#include "rigid_reference/BouncingBullet.h"
#ifdef RF_REFERENCE_JOLT
#include "rigid_reference/BouncingJolt.h"
#endif
#ifdef RF_REFERENCE_BOX3D
#include "rigid_reference/BouncingBox3d.h"
#endif

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace rf::benchmark;
struct Metrics {
    double energySq = 0, positionSq = 0, velocitySq = 0, maxRelativeEnergy = 0;
    double maxPenetration = 0, minHeight = height, meanHeight = 0, meanVelocity = 0;
    size_t samples = 0;
};

std::string number(double value) {
    // SDK serialization is float; accumulators remain double. Trace rounding is audited in Python.
    if (!std::isfinite(value) || !std::isfinite(float(value))) throw std::runtime_error("non-finite metric");
    return rf::numberText(float(value));
}

std::pair<double, double> analytic(double time) {
    // Exact free fall and e=1 reflection at centre height R; Newton's constant-g equations.
    // https://openstax.org/books/university-physics-volume-1/pages/3-5-free-fall
    constexpr double g = 9.81, r = 0.1;
    const double hit = std::sqrt(2 * (height - r) / g);
    if (time < hit) return {height - 0.5 * g * time * time, -g * time};
    const double phase = std::fmod(time - hit, 2 * hit);
    return {r + g * hit * phase - 0.5 * g * phase * phase, g * (hit - phase)};
}

double normSq(const rf::Vector3& value) {
    return double(value.x) * value.x + double(value.y) * value.y + double(value.z) * value.z;
}

void requireFinite(const State& state) {
    for (const auto v : {state.position, state.velocity, state.angularVelocity})
        if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z))
            throw std::runtime_error("non-finite body state");
}

double observe(const std::array<State, count>& states, double time, Metrics& metrics, bool accumulate) {
    const auto exact = analytic(time);
    double energy = 0;
    metrics.meanHeight = metrics.meanVelocity = 0;
    for (int i = 0; i < count; ++i) {
        const State& s = states[size_t(i)];
        requireFinite(s);
        // E=sum(m g y + m |v|^2/2 + I |w|^2/2), independent of each solver's energy helper.
        energy += mass * 9.81 * s.position.y + 0.5 * mass * normSq(s.velocity)
                  + 0.5 * sphereInertia * normSq(s.angularVelocity);
        metrics.minHeight = std::min(metrics.minHeight, double(s.position.y));
        metrics.maxPenetration = std::max(metrics.maxPenetration, std::max(0.0, 0.1 - s.position.y));
        metrics.meanHeight += s.position.y / double(count);
        metrics.meanVelocity += s.velocity.y / double(count);
        if (!accumulate) continue;
        const auto initial = startPosition(i);
        metrics.positionSq += std::pow(double(s.position.x) - initial.x, 2)
                              + std::pow(s.position.y - exact.first, 2)
                              + std::pow(double(s.position.z) - initial.z, 2);
        metrics.velocitySq += double(s.velocity.x) * s.velocity.x
                              + std::pow(s.velocity.y - exact.second, 2)
                              + double(s.velocity.z) * s.velocity.z;
    }
    metrics.maxRelativeEnergy = std::max(metrics.maxRelativeEnergy, std::fabs(energy / referenceEnergy - 1));
    if (accumulate) { metrics.energySq += std::pow(energy - referenceEnergy, 2); ++metrics.samples; }
    return energy;
}

void writeTrace(std::ofstream& trace, int frame, double time, double energy, const Metrics& metrics) {
    trace << frame << ',' << number(time) << ',' << number(energy) << ',' << number(metrics.meanHeight)
          << ',' << number(metrics.meanVelocity) << ',' << number(metrics.minHeight) << '\n';
    if (!trace) throw std::runtime_error("cannot write trace");
}

void writeMetrics(std::ofstream& out, const Metrics& m) {
    const double samples = double(m.samples);
    out << "  \"energy_mse_J2\": " << number(m.energySq / samples) << ",\n"
        << "  \"relative_energy_rmse\": " << number(std::sqrt(m.energySq / samples) / referenceEnergy) << ",\n"
        << "  \"max_relative_energy_error\": " << number(m.maxRelativeEnergy) << ",\n"
        << "  \"position_rmse_m\": " << number(std::sqrt(m.positionSq / (samples * count))) << ",\n"
        << "  \"velocity_rmse_m_s\": " << number(std::sqrt(m.velocitySq / (samples * count))) << ",\n"
        << "  \"sampled_max_penetration_m\": " << number(m.maxPenetration) << ",\n";
}

void writeSummary(const std::filesystem::path& output, const std::string& engine, float dt, int iterations,
                  int frames, const Work& accuracy, const Work& timing, const Metrics& metrics, double elapsed) {
    std::ofstream out(output / "summary.json");
    out << "{\n  \"engine\": \"" << engine << "\",\n  \"complete\": true,\n"
        << "  \"dt_s\": " << number(dt) << ",\n  \"iterations\": " << iterations << ",\n"
        << "  \"frames\": " << frames << ",\n  \"reference_energy_J\": " << number(referenceEnergy) << ",\n"
        << "  \"accuracy_accepted_time_s\": " << number(accuracy.time) << ",\n"
        << "  \"accuracy_accepted_substeps\": " << accuracy.accepted << ",\n"
        << "  \"sdk_impact_events\": " << accuracy.sdkImpactEvents << ",\n"
        << "  \"sdk_ccd_queries\": " << accuracy.sdkCcdQueries << ",\n"
        << "  \"accuracy_rejected_trials\": " << accuracy.rejected << ",\n";
    writeMetrics(out, metrics);
    out << "  \"timing_accepted_time_s\": " << number(timing.time) << ",\n"
        << "  \"timing_accepted_substeps\": " << timing.accepted << ",\n"
        << "  \"timing_rejected_trials\": " << timing.rejected << ",\n"
        << "  \"timing_pass_wall_s\": " << number(elapsed) << ",\n"
        << "  \"sdk_timing_broad_ms\": " << number(timing.sdkBroadMs) << ",\n"
        << "  \"sdk_timing_narrow_ms\": " << number(timing.sdkNarrowMs) << ",\n"
        << "  \"sdk_timing_solve_ms\": " << number(timing.sdkSolveMs) << ",\n"
        << "  \"sdk_timing_ccd_ms\": " << number(timing.sdkCcdMs) << ",\n"
        << "  \"bullet_header_version\": " << btGetVersion() << ",\n"
        << "  \"sizeof_bullet_scalar\": " << sizeof(btScalar) << "\n}\n";
    if (!out) throw std::runtime_error("cannot write summary");
}

template<class Engine>
int run(const std::string& name, float dt, int iterations, const std::filesystem::path& output) {
    const int frames = int(std::lround(duration / double(dt)));
    const bool option = name == "sdk-ccd" || name == "box3d-native";
    Metrics metrics;
    Work accuracy, timing;
    std::ofstream trace(output / "frames.csv");
    trace << "frame,time_s,energy_J,mean_height_m,mean_velocity_y_m_s,running_min_height_m\n";
    {
        Engine engine(iterations, option);
        if (std::fabs(observe(engine.states(), 0, metrics, false) - referenceEnergy) > 0.001)
            throw std::runtime_error("initial energy mismatch");
        for (int frame = 0; frame < frames; ++frame) {
            const double energy = observe(engine.states(), accuracy.time, metrics, true);
            writeTrace(trace, frame, accuracy.time, energy, metrics);
            if (!engine.step(dt, accuracy)) throw std::runtime_error("accuracy step rejected at frame " + std::to_string(frame));
        }
        observe(engine.states(), accuracy.time, metrics, false); // Geometry includes the final accepted pose.
    }
    Engine engine(iterations, option);
    const auto start = std::chrono::steady_clock::now();
    for (int frame = 0; frame < frames; ++frame)
        if (!engine.step(dt, timing)) throw std::runtime_error("timing step rejected at frame " + std::to_string(frame));
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    Metrics final;
    observe(engine.states(), timing.time, final, false);
    if (accuracy.accepted != timing.accepted || accuracy.rejected != timing.rejected || accuracy.time != timing.time)
        throw std::runtime_error("accuracy and timing workloads differ");
    writeSummary(output, name, dt, iterations, frames, accuracy, timing, metrics, elapsed);
    std::cout << name << " dt=" << number(dt) << " energy_RMSE="
              << number(std::sqrt(metrics.energySq / double(metrics.samples)) / referenceEnergy) << '\n';
    return 0;
}

void checkAnalytic() {
    for (int i = 0; i <= 40000; ++i) {
        const auto state = analytic(i * 0.0005);
        const double energy = count * mass * (9.81 * state.first + 0.5 * state.second * state.second);
        if (std::fabs(energy - referenceEnergy) > 1e-9 || state.first < 0.1 - 1e-12)
            throw std::runtime_error("analytic energy/geometry inconsistency");
    }
}
} // namespace

int main(int argc, char** argv) {
    std::filesystem::path output;
    bool created = false;
    try {
        if (argc < 4 || argc > 5) throw std::runtime_error("ENGINE DT NEW_DIRECTORY [iterations=20]");
        const std::string engine(argv[1]);
        bool supported = engine == "sdk" || engine == "sdk-ccd" || engine == "sdk-variational" || engine == "bullet";
#ifdef RF_REFERENCE_JOLT
        supported |= engine == "jolt";
#endif
#ifdef RF_REFERENCE_BOX3D
        supported |= engine == "box3d" || engine == "box3d-native";
#endif
        if (!supported) throw std::runtime_error("unknown or unavailable engine");
        float dt = 0, iterationValue = 20;
        if (!rf::parseNumber(argv[2], dt) || !std::isfinite(dt) || dt < 0.0001f || dt > 0.02f)
            throw std::runtime_error("dt must be in [0.0001,0.02]");
        if (argc == 5 && !rf::parseNumber(argv[4], iterationValue)) throw std::runtime_error("invalid iterations");
        if (!std::isfinite(iterationValue) || iterationValue < 1 || iterationValue > 1000
            || std::floor(iterationValue) != iterationValue) throw std::runtime_error("iterations must be an integer in [1,1000]");
        output = argv[3];
        if (!std::filesystem::create_directory(output)) throw std::runtime_error("output must be a new directory");
        created = true;
        checkAnalytic();
#ifdef RF_REFERENCE_JOLT
        if (engine == "jolt") return run<BouncingJolt>(engine, dt, int(iterationValue), output);
#endif
#ifdef RF_REFERENCE_BOX3D
        if (engine == "box3d" || engine == "box3d-native") return run<BouncingBox3d>(engine, dt, int(iterationValue), output);
#endif
        if (engine == "bullet") return run<BouncingBullet>(engine, dt, int(iterationValue), output);
        if (engine == "sdk-variational") return run<BouncingVariational>(engine, dt, int(iterationValue), output);
        return run<BouncingSdk>(engine, dt, int(iterationValue), output);
    } catch (const std::exception& error) {
        if (created) {
            // JSON deliberately carries no partial accuracy score; the process log has the cause.
            std::ofstream out(output / "failure.json");
            out << "{\"complete\": false}\n";
        }
        std::cerr << error.what() << '\n';
        return 2;
    }
}
