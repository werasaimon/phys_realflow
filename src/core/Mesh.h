#pragma once
// Indexed triangle mesh, procedural primitives and OBJ/STL import.

#include "math/Math.h"

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace rf {

struct TriMesh {
    std::vector<Vector3> positions;
    std::vector<std::array<uint32_t, 3>> triangles;

    bool empty() const { return triangles.empty(); }
    AABB bounds() const;
    Vector3 faceNormal(size_t t) const;   // unit normal, (0,0,0) for degenerate
    float surfaceArea() const;
    float signedVolume() const;        // > 0 for outward oriented closed meshes

    // Angle-weighted smooth vertex normals.
    std::vector<Vector3> vertexNormals() const;

    void append(const TriMesh& other);
    void transform(const Matrix3x3& R, const Vector3& scale, const Vector3& translation); // p' = R*(scale*p)+t
    void translate(const Vector3& t);
    void flipWinding();
    // Merge vertices closer than eps (needed for STL and for pseudo-normal SDF).
    void weld(float eps);
    // Removes zero-area triangles.
    void removeDegenerate(float areaEps = 1e-14f);
    // Flips the winding if the mesh is inside-out.
    void orientOutward();
    // Uniformly scale so the largest extent equals `size` and move bounds center to `center`.
    void fitTo(const Vector3& center, float size);
};

namespace primitives {
TriMesh box(const Vector3& halfExtents);
// Surface of revolution around the X axis; profile points are (x, radius).
TriMesh revolve(const std::vector<std::pair<float, float>>& profile, int segments);
TriMesh sphere(float radius, int segments = 48, int rings = 24);
TriMesh ellipsoid(const Vector3& radii, int segments = 48, int rings = 24);
// Cylinder of given radius and length with its axis along Z (classic cross-flow case).
TriMesh cylinder(float radius, float length, int segments = 48);
// Body of revolution with a NACA 00xx thickness distribution (streamlined "teardrop").
TriMesh streamlinedBody(float length, float thicknessRatio, int segments = 48, int stations = 48);
// NACA 4-digit wing, chord along +X starting at x=0, span along Z centered at 0.
TriMesh nacaWing(const std::string& code, float chord, float span, int stations = 60);
TriMesh cone(float radius, float length, int segments = 48); // apex pointing to -X
// Torus about the Y axis: a circle of radius `minor` in the (R, y) plane, centred at R = major,
// swept around the axis.
TriMesh torus(float major, float minor, int segments = 64, int rings = 24);
// Procedural teapot as closed parts whose union is the solid (body, lid knob, spout, curved
// handle); non-convex overall. Height ~ size, centred on the origin, spout along +X.
std::vector<TriMesh> teapotParts(float size);
// Procedural sitting bunny as closed overlapping parts (haunch, chest, head, snout, ears, tail,
// paws); height ~ size, centred on the origin, facing +X.
std::vector<TriMesh> bunnyParts(float size);
TriMesh merge(const std::vector<TriMesh>& parts);
// Splits every triangle into 4 (edge midpoints shared between neighbours) until no edge is longer
// than maxEdge (at most maxLevels times). The shape is unchanged; it just gets more vertices.
TriMesh subdivided(const TriMesh& mesh, float maxEdge, int maxLevels = 5);
} // namespace primitives

// Format detected by extension (.obj, .stl). Returns false and fills `error` on failure.
bool loadMesh(const std::string& path, TriMesh& out, std::string& error);
bool saveOBJ(const std::string& path, const TriMesh& mesh);

} // namespace rf
