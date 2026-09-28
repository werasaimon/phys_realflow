#pragma once
// Shared by every test file: the failure counter, CHECK, the runner with the RF_TEST filter and
// the helpers several tests use. Exit code of rf_tests = number of failures.
#include "spatial/BVH.h"
#include "core/Mesh.h"
#include "gas/GasSolver.h"
#include "rigid/RigidWorld.h"
#include "scene/Simulation.h"
#include "samples/Samples.h"
#include "samples/Models.h"
#include "particles/ParticleSystem.h"
#include "rigid/BroadPhase.h"
#include "spatial/AABBTree.h"
#include "rigid/TimeOfImpact.h"
#include "rigid/ConvexDecomposition.h"
#include "rigid/NarrowPhase.h"

#include <algorithm>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <random>
#include <string>
#include <vector>


using namespace rf;

inline int g_failures = 0;
#define CHECK(cond, ...)                                                          \
    do {                                                                          \
        if (!(cond)) {                                                            \
            ++g_failures;                                                         \
            std::printf("  FAIL %s:%d: %s  ", __FILE__, __LINE__, #cond);         \
            std::printf(__VA_ARGS__);                                             \
            std::printf("\n");                                                    \
        }                                                                         \
    } while (0)


// Every test that ran: for the summary and the JUnit report (RF_JUNIT=<file>, read by CI).
struct TestResult {
    std::string name;
    double ms;
    int failures;
};
inline std::vector<TestResult> g_results;

inline void run(const char* name, const std::function<void()>& fn) {
    if (const char* only = std::getenv("RF_TEST"); only && *only && !std::strstr(name, only)) return;
    auto t0 = std::chrono::steady_clock::now();
    int before = g_failures;
    fn();
    double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::printf("[%s] %s (%.0f ms)\n", g_failures == before ? " OK " : "FAIL", name, ms);
    std::fflush(stdout); // piped output is block-buffered: a crash would lose the last tests' lines
    g_results.push_back({name, ms, g_failures - before});
}

// JUnit XML of the results: one <testcase> per test, <failure> with the number of failed CHECKs
// (their messages are in the log above).
inline void writeJUnit(const char* path) {
    FILE* f = std::fopen(path, "w");
    if (!f) return;
    double total = 0;
    int failed = 0;
    for (const TestResult& r : g_results) { total += r.ms; failed += r.failures > 0; }
    std::fprintf(f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    std::fprintf(f, "<testsuite name=\"rf_tests\" tests=\"%zu\" failures=\"%d\" time=\"%.3f\">\n", g_results.size(), failed, total / 1000);
    for (const TestResult& r : g_results) {
        std::string name;
        for (char c : r.name) name += c == '<' ? "&lt;" : c == '>' ? "&gt;" : c == '&' ? "&amp;" : c == '"' ? "&quot;" : std::string(1, c);
        std::fprintf(f, "  <testcase classname=\"rf_tests\" name=\"%s\" time=\"%.3f\"", name.c_str(), r.ms / 1000);
        if (r.failures > 0) std::fprintf(f, ">\n    <failure message=\"%d failed CHECK(s), see the log\"/>\n  </testcase>\n", r.failures);
        else std::fprintf(f, "/>\n");
    }
    std::fprintf(f, "</testsuite>\n");
    std::fclose(f);
}

// FNV-1a (Fowler, Noll, Vo): each byte is xor-ed in, then the hash is multiplied by the FNV prime.
// Components are fed one by one - never a struct's raw bytes, whose padding is not part of the state.
struct StateHash {
    uint64_t h = 14695981039346656037ull; // the FNV-1a offset basis
    void add(float v) {
        unsigned char b[4];
        std::memcpy(b, &v, 4);
        for (unsigned char c : b) {
            h ^= c;
            h *= 1099511628211ull; // the FNV prime for 64 bits
        }
    }
    void add(const Vector3& v) { add(v.x); add(v.y); add(v.z); }
    void add(const Quaternion& q) { add(q.w); add(q.x); add(q.y); add(q.z); }
    void add(const std::vector<float>& d) { for (float v : d) add(v); }
};

// Deepest overlap between any two dynamic bodies (GJK/EPA): > tolerance means superposition.
inline float maxOverlap(const RigidWorld& w) {
    float worst = 0;
    const auto& B = w.bodies();
    for (size_t i = 0; i < B.size(); ++i)
        for (size_t j = i + 1; j < B.size(); ++j) {
            if (!B[i].worldBounds().overlaps(B[j].worldBounds())) continue;
            PenetrationResult pr;
            if (penetration(B[i].posed(), B[j].posed(), pr)) worst = std::max(worst, pr.depth);
        }
    return worst;
}

// Exact distance between two separated boxes: vertex-box (both ways) and all edge-edge pairs.
inline float exactBoxDistance(const PosedShape& A, const Vector3& ha, const PosedShape& B, const Vector3& hb) {
    auto corners = [](const PosedShape& P, const Vector3& h) {
        std::vector<Vector3> c;
        for (int i = 0; i < 8; ++i) c.push_back(P.p + P.R * Vector3((i & 1) ? h.x : -h.x, (i & 2) ? h.y : -h.y, (i & 4) ? h.z : -h.z));
        return c;
    };
    auto pointBox = [](const Vector3& p, const PosedShape& P, const Vector3& h) {
        Vector3 q = P.R.transposed() * (p - P.p);
        Vector3 d = vmax(vabs(q) - h, Vector3(0.0f));
        return length(d);
    };
    auto segSeg = [](Vector3 p1, Vector3 q1, Vector3 p2, Vector3 q2) {
        Vector3 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
        float a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r), c = dot(d1, r), b = dot(d1, d2), den = a * e - b * b;
        float sP = den > 1e-12f ? clampv((b * f - c * e) / den, 0.0f, 1.0f) : 0.0f;
        float t = (b * sP + f) / e;
        if (t < 0) { t = 0; sP = clampv(-c / a, 0.0f, 1.0f); }
        else if (t > 1) { t = 1; sP = clampv((b - c) / a, 0.0f, 1.0f); }
        return length((p1 + d1 * sP) - (p2 + d2 * t));
    };
    const int edges[12][2] = {{0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 2}, {1, 3}, {4, 6}, {5, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    auto ca = corners(A, ha), cb = corners(B, hb);
    float best = 1e9f;
    for (auto& p : ca) best = std::min(best, pointBox(p, B, hb));
    for (auto& p : cb) best = std::min(best, pointBox(p, A, ha));
    for (auto& e1 : edges)
        for (auto& e2 : edges) best = std::min(best, segSeg(ca[e1[0]], ca[e1[1]], cb[e2[0]], cb[e2[1]]));
    return best;
}

// Worst penetration between bodies, compounds tested part by part (the support map of a compound
// is the hull of all its parts, so a direct GJK on it would report false overlaps).
inline float maxPartOverlap(const RigidWorld& w) {
    auto expand = [](const RigidBody& b) {
        std::vector<PosedShape> out;
        PosedShape ps = b.posed();
        if (b.type() == ShapeType::Compound) {
            for (const auto& c : static_cast<const CompoundShape*>(b.shape.get())->children())
                out.push_back({c.shape.get(), ps.R * c.R, ps.p + ps.R * c.t});
        } else {
            out.push_back(ps);
        }
        return out;
    };
    float worst = 0;
    const auto& B = w.bodies();
    for (size_t i = 0; i < B.size(); ++i)
        for (size_t j = i + 1; j < B.size(); ++j) {
            if (!B[i].worldBounds().overlaps(B[j].worldBounds())) continue;
            for (const PosedShape& a : expand(B[i]))
                for (const PosedShape& b : expand(B[j])) {
                    PenetrationResult pr;
                    if (penetration(a, b, pr)) worst = std::max(worst, pr.depth);
                }
        }
    return worst;
}

// RMS distance of a soft body's particles from the best rigid fit of their rest shape, relative
// to the body size (0 = undeformed).
inline float softShapeError(const ParticleSystem& s, const SoftBody& b, const std::vector<Vector3>& rest) {
    Vector3 c(0.0f), c0(0.0f);
    for (size_t k = 0; k < b.particles.size(); ++k) c += s.positions()[b.particles[k]], c0 += rest[k];
    c /= float(b.particles.size());
    c0 /= float(b.particles.size());
    Matrix3x3 A = Matrix3x3::zero();
    for (size_t k = 0; k < b.particles.size(); ++k) A += Matrix3x3::outer(s.positions()[b.particles[k]] - c, rest[k] - c0);
    const Matrix3x3 R = extractRotation(A, Quaternion(), 50).toMatrix3x3();
    float e = 0, size = 0;
    for (size_t k = 0; k < b.particles.size(); ++k) {
        e += length2(s.positions()[b.particles[k]] - (c + R * (rest[k] - c0)));
        size = std::max(size, length(rest[k] - c0));
    }
    return std::sqrt(e / float(b.particles.size())) / size;
}
