// Voronoi fracture (rigid/Fracture.h): the cells of a box and of a sphere must fill the body
// exactly, each cell convex and watertight, each seed in its own cell - and each cell must be a
// body the rigid solver can take (a ConvexHullShape with a volume and a finite inertia).
#include "TestRunner.h"

#include "core/Mesh.h"
#include "rigid/Fracture.h"
#include "rigid/Shapes.h"

#include <chrono>
#include <cmath>
#include <map>

using namespace rf;

namespace {

// Every edge of a closed mesh belongs to exactly two triangles.
bool watertight(const TriMesh& m) {
    std::map<std::pair<uint32_t, uint32_t>, int> edges;
    for (const auto& t : m.triangles)
        for (int k = 0; k < 3; ++k) {
            uint32_t a = t[k], b = t[(k + 1) % 3];
            edges[{std::min(a, b), std::max(a, b)}]++;
        }
    for (const auto& [e, n] : edges)
        if (n != 2) return false;
    return !edges.empty();
}

// Every vertex lies on the inner side of every face plane.
bool convex(const TriMesh& m, float tol) {
    for (size_t t = 0; t < m.triangles.size(); ++t) {
        const Vector3 n = m.faceNormal(t);
        const Vector3& o = m.positions[m.triangles[t][0]];
        for (const Vector3& p : m.positions)
            if (dot(n, p - o) > tol) return false;
    }
    return true;
}

void checkCells(const char* name, const TriMesh& hull, const std::vector<Vector3>& seeds) {
    auto t0 = std::chrono::steady_clock::now();
    const std::vector<TriMesh> cells = voronoiCells(hull, seeds);
    const float ms = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
    const float size = maxComp(hull.bounds().extent()), tol = 1e-5f * size;
    float total = 0, vmin = 1e30f, vmax = 0;
    int leaky = 0, concave = 0, badShape = 0, seedErrors = 0;
    for (size_t i = 0; i < cells.size(); ++i) {
        const float v = cells[i].signedVolume();
        total += v;
        vmin = std::min(vmin, v);
        vmax = std::max(vmax, v);
        if (!watertight(cells[i])) ++leaky;
        if (!convex(cells[i], tol)) ++concave;
        // The seed is in its own cell and in no other one.
        if (!insideConvex(cells[i], seeds[i], tol)) ++seedErrors;
        for (size_t j = 0; j < cells.size(); ++j)
            if (j != i && insideConvex(cells[j], seeds[i], -tol)) ++seedErrors;
        ConvexHullShape shape(cells[i]);
        const Vector3 I = shape.unitInertia();
        if (!(shape.volume() > 0) || !std::isfinite(I.x + I.y + I.z) || !(I.x > 0 && I.y > 0 && I.z > 0)) ++badShape;
    }
    const float hullVolume = hull.signedVolume();
    std::printf("  %s: %zu seeds -> %zu cells in %.1f ms; volume %.6f vs hull %.6f (%.1e rel), cells %.2e .. %.2e m^3, "
                "leaky %d, concave %d, seed errors %d, bad shapes %d\n",
                name, seeds.size(), cells.size(), ms, total, hullVolume, std::fabs(total / hullVolume - 1), vmin, vmax, leaky,
                concave, seedErrors, badShape);
    CHECK(cells.size() == seeds.size(), "%s: %zu cells for %zu seeds", name, cells.size(), seeds.size());
    CHECK(std::fabs(total / hullVolume - 1) < 1e-4f, "%s: cells do not fill the body: %f vs %f", name, total, hullVolume);
    CHECK(leaky == 0, "%s: %d cells are not watertight", name, leaky);
    CHECK(concave == 0, "%s: %d cells are not convex", name, concave);
    CHECK(seedErrors == 0, "%s: %d seed/cell mismatches", name, seedErrors);
    CHECK(badShape == 0, "%s: %d cells make no valid convex shape", name, badShape);
}

} // namespace

void testVoronoiFracture() {
    // A brick with seeds spread evenly, and a sphere hit at its surface: the pieces crowd
    // around the impact point as in brittle fracture.
    const TriMesh box = primitives::box({0.5f, 0.3f, 0.15f});
    checkCells("box, 40 uniform seeds", box, uniformSeeds(box, 40, 3));
    const TriMesh sphere = primitives::sphere(0.25f, 24, 12);
    const std::vector<Vector3> seeds = impactSeeds(sphere, {0.0f, 0.0f, 0.25f}, 25, 0.5f, 5);
    checkCells("sphere, 25 impact seeds", sphere, seeds);
    // Impact seeds are crowded: the median seed is closer to the impact than the mean would be.
    float near = 0;
    for (const Vector3& s : seeds) near += length(s - Vector3(0.0f, 0.0f, 0.25f)) < 0.2f ? 1.0f : 0.0f;
    std::printf("  impact seeds within 0.2 m of the hit: %.0f of %zu\n", near, seeds.size());
    CHECK(near >= 0.5f * float(seeds.size()), "impact seeds are not crowded around the hit (%.0f of %zu)", near, seeds.size());
}
