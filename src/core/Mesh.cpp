#include "core/Mesh.h"

#include <algorithm>
#include <map>

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace rf {

AABB TriMesh::bounds() const {
    AABB b;
    for (const Vector3& p : positions) b.expand(p);
    return b;
}

Vector3 TriMesh::faceNormal(size_t t) const {
    const auto& tri = triangles[t];
    return normalize(cross(positions[tri[1]] - positions[tri[0]], positions[tri[2]] - positions[tri[0]]));
}

float TriMesh::surfaceArea() const {
    double a = 0;
    for (const auto& t : triangles)
        a += 0.5 * length(cross(positions[t[1]] - positions[t[0]], positions[t[2]] - positions[t[0]]));
    return float(a);
}

float TriMesh::signedVolume() const {
    // Tetrahedra from a vertex of the mesh, not from the world origin: for a model far from the
    // origin (CAD files in mm) the float terms about the origin would cancel and flip the sign.
    if (positions.empty()) return 0.0f;
    const Vector3 o = positions[0];
    double v = 0;
    for (const auto& t : triangles)
        v += dot(positions[t[0]] - o, cross(positions[t[1]] - o, positions[t[2]] - o)) / 6.0;
    return float(v);
}

std::vector<Vector3> TriMesh::vertexNormals() const {
    std::vector<Vector3> n(positions.size(), Vector3(0.0f));
    for (const auto& t : triangles) {
        Vector3 fn = normalize(cross(positions[t[1]] - positions[t[0]], positions[t[2]] - positions[t[0]]));
        for (int c = 0; c < 3; ++c) {
            Vector3 e1 = normalize(positions[t[(c + 1) % 3]] - positions[t[c]]);
            Vector3 e2 = normalize(positions[t[(c + 2) % 3]] - positions[t[c]]);
            float ang = std::acos(clampv(dot(e1, e2), -1.0f, 1.0f));
            n[t[c]] += fn * ang;
        }
    }
    for (Vector3& v : n) v = normalize(v);
    return n;
}

void TriMesh::append(const TriMesh& o) {
    uint32_t base = uint32_t(positions.size());
    positions.insert(positions.end(), o.positions.begin(), o.positions.end());
    for (auto t : o.triangles) triangles.push_back({t[0] + base, t[1] + base, t[2] + base});
}

void TriMesh::transform(const Matrix3x3& R, const Vector3& scale, const Vector3& t) {
    for (Vector3& p : positions) p = R * (p * scale) + t;
    // A mirroring scale flips orientation.
    if (scale.x * scale.y * scale.z < 0) flipWinding();
}

void TriMesh::translate(const Vector3& t) {
    for (Vector3& p : positions) p += t;
}

void TriMesh::flipWinding() {
    for (auto& t : triangles) std::swap(t[1], t[2]);
}

void TriMesh::weld(float eps) {
    if (positions.empty()) return;
    const float inv = 1.0f / std::max(eps, 1e-12f);
    auto key = [&](int x, int y, int z) {
        return (uint64_t(uint32_t(x) & 0x1FFFFF) << 42) | (uint64_t(uint32_t(y) & 0x1FFFFF) << 21) |
               uint64_t(uint32_t(z) & 0x1FFFFF);
    };
    std::unordered_map<uint64_t, std::vector<uint32_t>> cells;
    std::vector<Vector3> out;
    std::vector<uint32_t> remap(positions.size());
    const float eps2 = eps * eps;
    for (size_t i = 0; i < positions.size(); ++i) {
        const Vector3& p = positions[i];
        int cx = int(std::floor(p.x * inv)), cy = int(std::floor(p.y * inv)), cz = int(std::floor(p.z * inv));
        int found = -1;
        for (int dz = -1; dz <= 1 && found < 0; ++dz)
            for (int dy = -1; dy <= 1 && found < 0; ++dy)
                for (int dx = -1; dx <= 1 && found < 0; ++dx) {
                    auto it = cells.find(key(cx + dx, cy + dy, cz + dz));
                    if (it == cells.end()) continue;
                    for (uint32_t j : it->second)
                        if (length2(out[j] - p) <= eps2) { found = int(j); break; }
                }
        if (found < 0) {
            found = int(out.size());
            out.push_back(p);
            cells[key(cx, cy, cz)].push_back(uint32_t(found));
        }
        remap[i] = uint32_t(found);
    }
    positions.swap(out);
    for (auto& t : triangles)
        for (auto& v : t) v = remap[v];
    removeDegenerate();
}

void TriMesh::removeDegenerate(float areaEps) {
    std::vector<std::array<uint32_t, 3>> kept;
    kept.reserve(triangles.size());
    for (const auto& t : triangles) {
        if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2]) continue;
        Vector3 c = cross(positions[t[1]] - positions[t[0]], positions[t[2]] - positions[t[0]]);
        if (length2(c) <= areaEps * areaEps) continue;
        kept.push_back(t);
    }
    triangles.swap(kept);
}

void TriMesh::orientOutward() {
    if (signedVolume() < 0) flipWinding();
}

void TriMesh::fitTo(const Vector3& center, float size) {
    AABB b = bounds();
    if (!b.valid()) return;
    float ext = maxComp(b.extent());
    float s = ext > 0 ? size / ext : 1.0f;
    Vector3 c = b.center();
    for (Vector3& p : positions) p = (p - c) * s + center;
}

// ---------------------------------------------------------------------------
// Primitives
// ---------------------------------------------------------------------------
namespace primitives {

TriMesh box(const Vector3& h) {
    TriMesh m;
    for (int i = 0; i < 8; ++i)
        m.positions.push_back({(i & 1) ? h.x : -h.x, (i & 2) ? h.y : -h.y, (i & 4) ? h.z : -h.z});
    // Each face as quad (a,b,c,d) counter-clockwise seen from outside.
    const int q[6][4] = {{0, 4, 6, 2}, {1, 3, 7, 5}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 2, 3, 1}, {4, 5, 7, 6}};
    for (auto& f : q) {
        m.triangles.push_back({uint32_t(f[0]), uint32_t(f[1]), uint32_t(f[2])});
        m.triangles.push_back({uint32_t(f[0]), uint32_t(f[2]), uint32_t(f[3])});
    }
    m.orientOutward();
    return m;
}

TriMesh revolve(const std::vector<std::pair<float, float>>& profile, int seg) {
    TriMesh m;
    const float rEps = 1e-7f;
    std::vector<std::vector<uint32_t>> rings;
    for (auto [x, r] : profile) {
        std::vector<uint32_t> ring(seg);
        if (r <= rEps) {
            uint32_t id = uint32_t(m.positions.size());
            m.positions.push_back({x, 0, 0});
            for (int j = 0; j < seg; ++j) ring[j] = id;
        } else {
            for (int j = 0; j < seg; ++j) {
                float th = 2.0f * kPi * j / seg;
                ring[j] = uint32_t(m.positions.size());
                m.positions.push_back({x, r * std::cos(th), r * std::sin(th)});
            }
        }
        rings.push_back(std::move(ring));
    }
    for (size_t i = 0; i + 1 < rings.size(); ++i)
        for (int j = 0; j < seg; ++j) {
            int jn = (j + 1) % seg;
            uint32_t a = rings[i][j], b = rings[i][jn], c = rings[i + 1][jn], d = rings[i + 1][j];
            m.triangles.push_back({a, b, c});
            m.triangles.push_back({a, c, d});
        }
    m.removeDegenerate();
    m.orientOutward();
    return m;
}

TriMesh ellipsoid(const Vector3& radii, int seg, int rings) {
    std::vector<std::pair<float, float>> prof;
    for (int i = 0; i <= rings; ++i) {
        float t = kPi * i / rings;
        prof.push_back({-std::cos(t), std::sin(t)});
    }
    TriMesh m = revolve(prof, seg);
    for (Vector3& p : m.positions) p = p * radii;
    return m;
}

TriMesh sphere(float r, int seg, int rings) { return ellipsoid(Vector3(r), seg, rings); }

TriMesh cylinder(float radius, float length, int seg) {
    float h = length * 0.5f;
    std::vector<std::pair<float, float>> prof = {{-h, 0}, {-h, radius}, {h, radius}, {h, 0}};
    TriMesh m = revolve(prof, seg);
    // Axis X -> Z.
    for (Vector3& p : m.positions) p = {p.z, p.y, p.x};
    m.flipWinding(); // the swap above is a reflection
    return m;
}

static float nacaThickness(float xi, float t) {
    return 5.0f * t * (0.2969f * std::sqrt(xi) - 0.1260f * xi - 0.3516f * xi * xi + 0.2843f * xi * xi * xi -
                       0.1036f * xi * xi * xi * xi);
}

TriMesh streamlinedBody(float length, float thick, int seg, int stations) {
    std::vector<std::pair<float, float>> prof;
    for (int i = 0; i <= stations; ++i) {
        float xi = 0.5f * (1.0f - std::cos(kPi * i / stations));
        float r = (i == 0 || i == stations) ? 0.0f : nacaThickness(xi, thick) * length;
        prof.push_back({xi * length - 0.5f * length, r});
    }
    return revolve(prof, seg);
}

TriMesh cone(float radius, float length, int seg) {
    std::vector<std::pair<float, float>> prof = {{-length * 0.5f, 0}, {length * 0.5f, radius}, {length * 0.5f, 0}};
    return revolve(prof, seg);
}

TriMesh nacaWing(const std::string& code, float chord, float span, int n) {
    float mC = 0.02f, pC = 0.4f, tC = 0.12f; // 2412, also for a code that is not four digits
    const bool digits = code.size() == 4 && std::all_of(code.begin(), code.end(), [](char c) { return c >= '0' && c <= '9'; });
    if (digits) {
        mC = (code[0] - '0') / 100.0f;
        pC = (code[1] - '0') / 10.0f;
        tC = ((code[2] - '0') * 10 + (code[3] - '0')) / 100.0f;
    }
    tC = std::max(tC, 0.01f);
    std::vector<Vector3> U(n + 1), L(n + 1);
    for (int i = 0; i <= n; ++i) {
        float xi = 0.5f * (1.0f - std::cos(kPi * i / n));
        float yt = (i == 0 || i == n) ? 0.0f : nacaThickness(xi, tC);
        float yc = 0, dyc = 0;
        if (mC > 0 && pC > 0) {
            if (xi < pC) {
                yc = mC / (pC * pC) * (2 * pC * xi - xi * xi);
                dyc = 2 * mC / (pC * pC) * (pC - xi);
            } else {
                yc = mC / ((1 - pC) * (1 - pC)) * ((1 - 2 * pC) + 2 * pC * xi - xi * xi);
                dyc = 2 * mC / ((1 - pC) * (1 - pC)) * (pC - xi);
            }
        }
        float th = std::atan(dyc);
        U[i] = Vector3(xi - yt * std::sin(th), yc + yt * std::cos(th), 0) * chord;
        L[i] = Vector3(xi + yt * std::sin(th), yc - yt * std::cos(th), 0) * chord;
    }
    // Closed loop: U0..Un then L(n-1)..L1.
    const int loopN = 2 * n;
    auto loopPt = [&](int k) { return k <= n ? U[k] : L[2 * n - k]; };
    auto uIdx = [&](int i) { return i; };
    auto lIdx = [&](int i) { return (i == 0 || i == n) ? i : 2 * n - i; };

    TriMesh m;
    for (int layer = 0; layer < 2; ++layer) {
        float z = layer == 0 ? -0.5f * span : 0.5f * span;
        for (int k = 0; k < loopN; ++k) {
            Vector3 p = loopPt(k);
            m.positions.push_back({p.x, p.y, z});
        }
    }
    auto id = [&](int layer, int k) { return uint32_t(layer * loopN + ((k % loopN) + loopN) % loopN); };
    for (int k = 0; k < loopN; ++k) {
        uint32_t a = id(0, k), b = id(0, k + 1), c = id(1, k + 1), d = id(1, k);
        m.triangles.push_back({a, c, b});
        m.triangles.push_back({a, d, c});
    }
    for (int layer = 0; layer < 2; ++layer)
        for (int i = 0; i < n; ++i) {
            uint32_t u0 = id(layer, uIdx(i)), u1 = id(layer, uIdx(i + 1));
            uint32_t l0 = id(layer, lIdx(i)), l1 = id(layer, lIdx(i + 1));
            if (layer == 0) {
                m.triangles.push_back({u0, u1, l1});
                m.triangles.push_back({u0, l1, l0});
            } else {
                m.triangles.push_back({u0, l1, u1});
                m.triangles.push_back({u0, l0, l1});
            }
        }
    m.removeDegenerate();
    m.orientOutward();
    return m;
}

// Tube along a circular arc (closed with caps). Every triangle is oriented against its local
// outward direction, so the part is a consistent closed surface.
static TriMesh arcTube(const Vector3& centre, float R, float r, float a0, float a1, int segs, int sides) {
    TriMesh m;
    auto ringCentre = [&](int i) {
        float a = a0 + (a1 - a0) * i / segs;
        return centre + Vector3(R * std::cos(a), R * std::sin(a), 0);
    };
    for (int i = 0; i <= segs; ++i) {
        float a = a0 + (a1 - a0) * i / segs;
        Vector3 radial(std::cos(a), std::sin(a), 0), zAxis(0, 0, 1);
        for (int k = 0; k < sides; ++k) {
            float t = 2 * kPi * k / sides;
            m.positions.push_back(ringCentre(i) + (radial * std::cos(t) + zAxis * std::sin(t)) * r);
        }
    }
    auto add = [&](uint32_t a, uint32_t b, uint32_t c, const Vector3& outward) {
        Vector3 n = cross(m.positions[b] - m.positions[a], m.positions[c] - m.positions[a]);
        if (dot(n, outward) < 0) std::swap(b, c);
        m.triangles.push_back({a, b, c});
    };
    for (int i = 0; i < segs; ++i)
        for (int k = 0; k < sides; ++k) {
            uint32_t a = i * sides + k, b = i * sides + (k + 1) % sides, c = (i + 1) * sides + (k + 1) % sides,
                     d = (i + 1) * sides + k;
            Vector3 mid = (ringCentre(i) + ringCentre(i + 1)) * 0.5f;
            Vector3 face = (m.positions[a] + m.positions[b] + m.positions[c] + m.positions[d]) * 0.25f;
            add(a, b, c, face - mid);
            add(a, c, d, face - mid);
        }
    for (int end = 0; end < 2; ++end) {
        int i = end == 0 ? 0 : segs;
        uint32_t cIdx = uint32_t(m.positions.size());
        m.positions.push_back(ringCentre(i));
        float a = a0 + (a1 - a0) * i / segs;
        Vector3 tangent(-std::sin(a), std::cos(a), 0);
        Vector3 outward = end == 0 ? -tangent * (a1 > a0 ? 1.0f : -1.0f) : tangent * (a1 > a0 ? 1.0f : -1.0f);
        for (int k = 0; k < sides; ++k) add(cIdx, i * sides + k, i * sides + (k + 1) % sides, outward);
    }
    return m;
}

std::vector<TriMesh> teapotParts(float size) {
    const float s = size / 0.75f; // the reference teapot below is 0.75 tall
    std::vector<TriMesh> parts;
    const Matrix3x3 xToY = Quaternion::fromAxisAngle({0, 0, 1}, kPi / 2).toMatrix3x3(); // revolve axis X -> Y
    // Body: surface of revolution.
    {
        std::vector<std::pair<float, float>> prof = {{0.0f, 0.0f},   {0.0f, 0.30f},  {0.03f, 0.37f}, {0.12f, 0.43f},
                                                     {0.25f, 0.46f}, {0.38f, 0.45f}, {0.48f, 0.41f}, {0.56f, 0.35f},
                                                     {0.61f, 0.28f}, {0.64f, 0.18f}, {0.64f, 0.0f}};
        TriMesh b = revolve(prof, 40);
        b.transform(xToY, Vector3(1.0f), Vector3(0.0f));
        parts.push_back(b);
    }
    // Lid knob.
    {
        TriMesh k = ellipsoid({0.07f, 0.06f, 0.07f}, 20, 10);
        k.translate({0, 0.69f, 0});
        parts.push_back(k);
    }
    // Spout: a frustum from inside the body out and up.
    {
        Vector3 a(0.30f, 0.18f, 0), b(0.74f, 0.52f, 0);
        float L = length(b - a);
        std::vector<std::pair<float, float>> prof = {{0, 0}, {0, 0.10f}, {L, 0.05f}, {L, 0}};
        TriMesh sp = revolve(prof, 24);
        Vector3 d = normalize(b - a);
        float ang = std::atan2(d.y, d.x);
        sp.transform(Quaternion::fromAxisAngle({0, 0, 1}, ang).toMatrix3x3(), Vector3(1.0f), a);
        parts.push_back(sp);
    }
    // Handle: a tube along a C-shaped arc on the -X side; both ends sink into the body.
    parts.push_back(arcTube({-0.46f, 0.34f, 0}, 0.17f, 0.045f, 1.05f, 2 * kPi - 1.05f, 20, 14));
    for (TriMesh& p : parts) {
        for (Vector3& v : p.positions) v = (v - Vector3(0, 0.36f, 0)) * s; // centre vertically
        p.removeDegenerate();
        p.orientOutward();
    }
    return parts;
}

std::vector<TriMesh> bunnyParts(float size) {
    // Ellipsoids placed in a unit-height model (feet at y = 0), then scaled and centred.
    struct Blob { Vector3 radii, centre; Quaternion rot; };
    const Quaternion earTiltL = Quaternion::fromAxisAngle({0, 0, 1}, 0.35f) * Quaternion::fromAxisAngle({1, 0, 0}, -0.2f);
    const Quaternion earTiltR = Quaternion::fromAxisAngle({0, 0, 1}, 0.35f) * Quaternion::fromAxisAngle({1, 0, 0}, 0.2f);
    const Blob blobs[] = {
        {{0.36f, 0.30f, 0.28f}, {-0.05f, 0.30f, 0.0f}, Quaternion()},   // haunch
        {{0.22f, 0.26f, 0.22f}, {0.20f, 0.36f, 0.0f}, Quaternion()},    // chest
        {{0.20f, 0.17f, 0.16f}, {0.36f, 0.64f, 0.0f}, Quaternion()},    // head
        {{0.09f, 0.08f, 0.09f}, {0.53f, 0.60f, 0.0f}, Quaternion()},    // snout
        {{0.07f, 0.24f, 0.10f}, {0.28f, 0.92f, 0.08f}, earTiltL}, // left ear
        {{0.07f, 0.24f, 0.10f}, {0.28f, 0.92f, -0.08f}, earTiltR},// right ear
        {{0.09f, 0.09f, 0.09f}, {-0.42f, 0.34f, 0.0f}, Quaternion()},   // tail
        {{0.13f, 0.06f, 0.07f}, {0.28f, 0.10f, 0.10f}, Quaternion()},   // front paws
        {{0.13f, 0.06f, 0.07f}, {0.28f, 0.10f, -0.10f}, Quaternion()},
        {{0.20f, 0.07f, 0.09f}, {0.02f, 0.07f, 0.20f}, Quaternion()},   // hind feet
        {{0.20f, 0.07f, 0.09f}, {0.02f, 0.07f, -0.20f}, Quaternion()},
    };
    std::vector<TriMesh> parts;
    AABB bb;
    for (const Blob& b : blobs) {
        TriMesh e = ellipsoid(b.radii, 24, 12);
        e.transform(b.rot.toMatrix3x3(), Vector3(1.0f), b.centre);
        bb.expand(e.bounds());
        parts.push_back(std::move(e));
    }
    const float s = size / bb.extent().y;
    const Vector3 c = bb.center();
    for (TriMesh& p : parts) {
        for (Vector3& v : p.positions) v = (v - c) * s;
        p.removeDegenerate();
        p.orientOutward();
    }
    return parts;
}

TriMesh subdivided(const TriMesh& mesh, float maxEdge, int maxLevels) {
    TriMesh m = mesh;
    for (int level = 0; level < maxLevels; ++level) {
        float longest = 0;
        for (const auto& t : m.triangles)
            for (int k = 0; k < 3; ++k) longest = std::max(longest, length(m.positions[t[k]] - m.positions[t[(k + 1) % 3]]));
        if (longest <= maxEdge) break;
        std::map<std::pair<uint32_t, uint32_t>, uint32_t> midpoint; // edge -> new vertex
        auto mid = [&](uint32_t a, uint32_t b) {
            auto key = std::make_pair(std::min(a, b), std::max(a, b));
            auto it = midpoint.find(key);
            if (it != midpoint.end()) return it->second;
            uint32_t v = uint32_t(m.positions.size());
            m.positions.push_back((m.positions[a] + m.positions[b]) * 0.5f);
            midpoint.emplace(key, v);
            return v;
        };
        std::vector<std::array<uint32_t, 3>> tris;
        tris.reserve(m.triangles.size() * 4);
        for (const auto& t : m.triangles) {
            uint32_t ab = mid(t[0], t[1]), bc = mid(t[1], t[2]), ca = mid(t[2], t[0]);
            tris.push_back({t[0], ab, ca});
            tris.push_back({ab, t[1], bc});
            tris.push_back({ca, bc, t[2]});
            tris.push_back({ab, bc, ca});
        }
        m.triangles.swap(tris);
    }
    return m;
}

TriMesh merge(const std::vector<TriMesh>& parts) {
    TriMesh m;
    for (const TriMesh& p : parts) m.append(p);
    return m;
}

} // namespace primitives

// ---------------------------------------------------------------------------
// IO
// ---------------------------------------------------------------------------
static std::string lowerExt(const std::string& path) {
    auto dot = path.find_last_of('.');
    if (dot == std::string::npos) return {};
    std::string e = path.substr(dot + 1);
    for (char& c : e) c = char(std::tolower(static_cast<unsigned char>(c)));
    return e;
}

static bool loadOBJ(const std::string& path, TriMesh& out, std::string& err) {
    std::ifstream f(std::filesystem::u8path(path));
    if (!f) { err = "cannot open file"; return false; }
    std::string line;
    std::vector<int> face;
    while (std::getline(f, line)) {
        if (line.size() < 2) continue;
        if (line[0] == 'v' && line[1] == ' ') {
            Vector3 p;
            std::istringstream ss(line.substr(2));
            ss >> p.x >> p.y >> p.z;
            out.positions.push_back(p);
        } else if (line[0] == 'f' && line[1] == ' ') {
            face.clear();
            std::istringstream ss(line.substr(2));
            std::string tok;
            while (ss >> tok) {
                int idx = std::atoi(tok.c_str()); // stops at '/'
                if (idx < 0) idx = int(out.positions.size()) + idx;
                else idx -= 1;
                if (idx < 0 || idx >= int(out.positions.size())) { err = "bad face index"; return false; }
                face.push_back(idx);
            }
            for (size_t i = 1; i + 1 < face.size(); ++i)
                out.triangles.push_back({uint32_t(face[0]), uint32_t(face[i]), uint32_t(face[i + 1])});
        }
    }
    return true;
}

static bool loadSTL(const std::string& path, TriMesh& out, std::string& err) {
    std::ifstream f(std::filesystem::u8path(path), std::ios::binary);
    if (!f) { err = "cannot open file"; return false; }
    std::vector<char> data((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    if (data.size() >= 84) {
        uint32_t n;
        std::memcpy(&n, data.data() + 80, 4);
        if (size_t(84) + size_t(n) * 50 == data.size()) {
            for (uint32_t t = 0; t < n; ++t) {
                const char* rec = data.data() + 84 + size_t(t) * 50 + 12;
                for (int v = 0; v < 3; ++v) {
                    float xyz[3];
                    std::memcpy(xyz, rec + v * 12, 12);
                    out.positions.push_back({xyz[0], xyz[1], xyz[2]});
                }
                out.triangles.push_back({3 * t, 3 * t + 1, 3 * t + 2});
            }
            return true;
        }
    }
    std::istringstream ss(std::string(data.begin(), data.end()));
    std::string tok;
    while (ss >> tok) {
        if (tok == "vertex") {
            Vector3 p;
            ss >> p.x >> p.y >> p.z;
            out.positions.push_back(p);
            if (out.positions.size() % 3 == 0) {
                uint32_t b = uint32_t(out.positions.size() - 3);
                out.triangles.push_back({b, b + 1, b + 2});
            }
        }
    }
    if (out.triangles.empty()) { err = "no triangles found"; return false; }
    return true;
}

bool loadMesh(const std::string& path, TriMesh& out, std::string& error) {
    out = TriMesh();
    std::string ext = lowerExt(path);
    bool ok = false;
    if (ext == "obj") ok = loadOBJ(path, out, error);
    else if (ext == "stl") ok = loadSTL(path, out, error);
    else error = "unsupported format (use .obj or .stl)";
    if (!ok) return false;
    AABB b = out.bounds();
    if (!b.valid() || out.triangles.empty()) { error = "empty mesh"; return false; }
    out.weld(1e-5f * length(b.extent()));
    out.orientOutward();
    return true;
}

bool saveOBJ(const std::string& path, const TriMesh& m) {
    std::ofstream f(std::filesystem::u8path(path));
    if (!f) return false;
    for (const Vector3& p : m.positions) f << "v " << p.x << ' ' << p.y << ' ' << p.z << '\n';
    for (const auto& t : m.triangles) f << "f " << t[0] + 1 << ' ' << t[1] + 1 << ' ' << t[2] + 1 << '\n';
    return bool(f);
}

} // namespace rf
