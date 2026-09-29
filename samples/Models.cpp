// Non-convex models decomposed once into convex parts (cached, shared by all bodies).
#include "samples/Models.h"

#include "core/Mesh.h"
#include "rigid/ConvexDecomposition.h"

#include <cmath>

namespace rf {

// A solid given as overlapping closed parts -> minimal set of convex hulls fitted to its smooth surface.
static std::shared_ptr<const CompoundShape> decompose(const std::vector<TriMesh>& parts) {
    DecompositionParams prm;
    prm.resolution = 24;
    prm.maxHullVertices = 64;
    return std::make_shared<CompoundShape>(convexDecomposition(parts, prm), primitives::merge(parts));
}

std::shared_ptr<const CompoundShape> teapotShape() {
    static std::shared_ptr<const CompoundShape> shape = decompose(primitives::teapotParts(0.3f));
    return shape;
}

std::shared_ptr<const CompoundShape> bunnyShape() {
    static std::shared_ptr<const CompoundShape> shape = decompose(primitives::bunnyParts(0.3f));
    return shape;
}

// With the default decomposition, finer than the teapot's: a ring's hole must stay open for the
// next ring of the chain (13 parts, reaching 2.9 mm into it).
std::shared_ptr<const CompoundShape> ringShape() {
    static std::shared_ptr<const CompoundShape> shape = [] {
        const TriMesh ring = primitives::torus(0.05f, 0.012f);
        return std::make_shared<const CompoundShape>(convexDecomposition({ring}), ring);
    }();
    return shape;
}

namespace {

// A box of half extents `half` turned by `angle` about z and put at `at`: a part of a compound.
TriMesh boxPart(const Vector3& half, float angle, const Vector3& at) {
    TriMesh b = primitives::box(half);
    b.transform(Quaternion::fromAxisAngle({0, 0, 1}, angle).toMatrix3x3(), Vector3(1.0f), at);
    return b;
}

// A solid of revolution about y from a closed profile of (height, radius) points.
TriMesh revolvedY(const std::vector<std::pair<float, float>>& profile) {
    TriMesh m = primitives::revolve(profile, 48);
    for (Vector3& p : m.positions) p = {p.y, p.x, p.z}; // the axis from x to y
    m.flipWinding();                                       // (a reflection)
    m.orientOutward();
    return m;
}

std::shared_ptr<const CompoundShape> fromParts(const std::vector<TriMesh>& parts) {
    return std::make_shared<const CompoundShape>(parts, primitives::merge(parts));
}

std::shared_ptr<const CompoundShape> decomposed(const TriMesh& m, int resolution, int maxParts) {
    DecompositionParams prm;
    prm.resolution = resolution;
    prm.maxParts = maxParts;
    return std::make_shared<const CompoundShape>(convexDecomposition({m}, prm), m);
}

} // namespace

std::shared_ptr<const CompoundShape> wheelShape() {
    static std::shared_ptr<const CompoundShape> shape = [] {
        const int n = 48;
        const float R = 0.3f, rim = 0.03f, width = 0.08f;
        std::vector<TriMesh> parts;
        for (int k = 0; k < n; ++k) {
            const float a = 2 * kPi * float(k) / n, r = R - 0.5f * rim;
            parts.push_back(boxPart({R * std::sin(kPi / n) + 0.002f, 0.5f * rim, 0.5f * width}, a + 0.5f * kPi, {r * std::cos(a), r * std::sin(a), 0}));
        }
        for (int k = 0; k < 6; ++k) {
            const float a = 2 * kPi * float(k) / 6, r = 0.5f * (R - rim + 0.05f);
            parts.push_back(boxPart({0.5f * (R - rim - 0.05f), 0.012f, 0.012f}, a, {r * std::cos(a), r * std::sin(a), 0}));
        }
        parts.push_back(boxPart({0.05f, 0.05f, 0.5f * width}, 0, Vector3(0.0f)));
        return fromParts(parts);
    }();
    return shape;
}

std::shared_ptr<const CompoundShape> tableShape() {
    static std::shared_ptr<const CompoundShape> shape = [] {
        std::vector<TriMesh> parts = {boxPart({0.3f, 0.02f, 0.2f}, 0, {0, 0.42f, 0})};
        for (int k = 0; k < 4; ++k) parts.push_back(boxPart({0.02f, 0.2f, 0.02f}, 0, {k & 1 ? 0.26f : -0.26f, 0.2f, k & 2 ? 0.16f : -0.16f}));
        return fromParts(parts);
    }();
    return shape;
}

// A shell this thin wants a finer grid than the default: 48 voxels across, 48 parts.
std::shared_ptr<const CompoundShape> bowlShape() {
    static std::shared_ptr<const CompoundShape> shape = [] {
        const float R = 0.6f, t = 0.04f;
        const int n = 24;
        std::vector<std::pair<float, float>> profile;
        for (int i = 0; i <= n; ++i) profile.push_back({-R * std::cos(0.5f * kPi * i / n), R * std::sin(0.5f * kPi * i / n)});
        for (int i = n; i >= 0; --i) profile.push_back({-(R - t) * std::cos(0.5f * kPi * i / n), (R - t) * std::sin(0.5f * kPi * i / n)});
        return decomposed(revolvedY(profile), 48, 48);
    }();
    return shape;
}

// The wall of a nesting cup is thin and sloped: a hull reaching 1 mm into the cup lifts the next cup
// by 1 mm / taper = 5 mm, so the grid is finer still (64 voxels, up to 64 parts).
std::shared_ptr<const CompoundShape> cupShape() {
    static std::shared_ptr<const CompoundShape> shape = [] {
        const float Rb = 0.04f, Rt = 0.06f, h = 0.05f, wall = 0.004f;
        const float k = (Rt - Rb) / (2 * h), across = wall * std::sqrt(1 + k * k); // the sloped wall's horizontal thickness
        return decomposed(revolvedY({{-h, 0}, {-h, Rb}, {h, Rt}, {h, Rt - across}, {-h + wall, Rb + k * wall - across}, {-h + wall, 0}}), 64, 64);
    }();
    return shape;
}

} // namespace rf
