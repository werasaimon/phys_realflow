// Soft bodies: the tetrahedra cut from the particle lattice, their colours, the skin riding in them
// and the surface normals they carry. The material itself is solved in SoftBodySolver.cpp; the data
// layout is in SoftBody.h.
#include "particles/SoftBody.h"

#include <algorithm>

namespace rf {

std::vector<std::array<int, 4>> latticeTetrahedra(const std::vector<std::array<int, 3>>& cells) {
    std::vector<std::array<int, 4>> tets;
    if (cells.empty()) return tets;
    std::array<int, 3> lo = cells[0], hi = cells[0];
    for (const auto& c : cells)
        for (int a = 0; a < 3; ++a) lo[a] = std::min(lo[a], c[a]), hi[a] = std::max(hi[a], c[a]);
    const int nx = hi[0] - lo[0] + 1, ny = hi[1] - lo[1] + 1, nz = hi[2] - lo[2] + 1;
    std::vector<int> at(size_t(nx) * size_t(ny) * size_t(nz), -1); // lattice point -> its index, -1: none
    auto slot = [&](int i, int j, int k) -> int {
        i -= lo[0], j -= lo[1], k -= lo[2];
        if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz) return -1;
        return at[size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * size_t(k))];
    };
    for (size_t p = 0; p < cells.size(); ++p)
        at[size_t(cells[p][0] - lo[0]) + size_t(nx) * (size_t(cells[p][1] - lo[1]) + size_t(ny) * size_t(cells[p][2] - lo[2]))] = int(p);
    // Every cube with its lowest corner at a point: the six tetrahedra around one of its diagonals,
    // one per order of the axes: the diagonal's start, then one step along an axis, then two, then
    // its end. The diagonal alternates: a cube odd along an axis is mirrored along it. With the same
    // diagonal everywhere (Kuhn's plain lattice) the material is anisotropic - squeezed along y it
    // also shears along x = z - and a column of six soft barrels leaned over along that diagonal
    // and fell within a quarter of a second. Mirrored, neighbouring cubes still cut their shared
    // face along the same diagonal (both are mirrored alike along the face's two axes), and the
    // shear of one cube is undone by its neighbour's.
    static const int orders[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (const auto& c : cells) {
        auto corner = [&](int dx, int dy, int dz) { // local corner (0/1 per axis) of the mirrored cube
            const int d[3] = {dx, dy, dz};
            int g[3];
            for (int a = 0; a < 3; ++a) g[a] = c[a] + ((c[a] & 1) ? 1 - d[a] : d[a]);
            return slot(g[0], g[1], g[2]);
        };
        const int start = corner(0, 0, 0), end = corner(1, 1, 1);
        if (start < 0 || end < 0) continue;
        for (const auto& o : orders) {
            int one[3] = {0, 0, 0}, two[3] = {0, 0, 0};
            one[o[0]] = two[o[0]] = two[o[1]] = 1;
            const int second = corner(one[0], one[1], one[2]), third = corner(two[0], two[1], two[2]);
            if (second >= 0 && third >= 0) tets.push_back({start, second, third, end});
        }
    }
    return tets;
}

// Greedy colouring: every tetrahedron takes the lowest colour none of the tetrahedra sharing one of
// its particles has. A lattice point belongs to at most 24 tetrahedra; the colours come out 30-40.
static std::vector<int> colourTetrahedra(const std::vector<SoftTet>& tets, int firstParticle, size_t particleCount) {
    std::vector<std::vector<int>> taken(particleCount); // colours already used at every particle
    std::vector<int> colour(tets.size());
    for (size_t t = 0; t < tets.size(); ++t) {
        int c = 0;
        for (;; ++c) {
            bool used = false;
            for (int v : tets[t].v) {
                const auto& list = taken[size_t(v - firstParticle)];
                used = used || std::find(list.begin(), list.end(), c) != list.end();
            }
            if (!used) break;
        }
        colour[t] = c;
        for (int v : tets[t].v) taken[size_t(v - firstParticle)].push_back(c);
    }
    return colour;
}

void buildTetrahedra(SoftBody& body, const std::vector<std::array<int, 4>>& quads, const std::vector<Vector3>& rest) {
    std::vector<SoftTet> tets;
    for (std::array<int, 4> q : quads) {
        auto edges = [&] {
            const Vector3& x0 = rest[size_t(q[0])];
            return Matrix3x3::fromColumns(rest[size_t(q[1])] - x0, rest[size_t(q[2])] - x0, rest[size_t(q[3])] - x0);
        };
        Matrix3x3 D = edges();
        if (D.determinant() < 0) {
            std::swap(q[2], q[3]);
            D = edges();
        }
        const float det = D.determinant();
        if (det <= 0) continue; // degenerate
        SoftTet t;
        t.v = q;
        t.restVolume = det / 6.0f;
        t.restInverse = D.inverse();
        tets.push_back(t);
    }
    const int first = body.particles.empty() ? 0 : body.particles.front();
    const std::vector<int> colour = colourTetrahedra(tets, first, body.particles.size());
    const int colours = colour.empty() ? 0 : *std::max_element(colour.begin(), colour.end()) + 1;
    body.colourStart.assign(size_t(colours) + 1, 0);
    for (int c : colour) ++body.colourStart[size_t(c) + 1];
    for (int c = 0; c < colours; ++c) body.colourStart[size_t(c) + 1] += body.colourStart[size_t(c)];
    body.tets.assign(tets.size(), SoftTet());
    std::vector<int> next(body.colourStart.begin(), body.colourStart.end() - 1);
    for (size_t t = 0; t < tets.size(); ++t) body.tets[size_t(next[size_t(colour[t])]++)] = tets[t];
    buildNodes(body);
}

// Every particle's star (the tetrahedra around it) and rest volume, and the nodes' colours: greedy,
// the lowest colour no node whose star shares a particle with this one's has (27 on a lattice).
void buildNodes(SoftBody& body) {
    const size_t n = body.particles.size();
    const int first = n ? body.particles.front() : 0;
    body.nodes.assign(n, SoftNode());
    std::vector<std::vector<int>> star(n);
    for (int t = 0; t < int(body.tets.size()); ++t)
        for (int v : body.tets[size_t(t)].v) {
            star[size_t(v - first)].push_back(t);
            body.nodes[size_t(v - first)].restVolume += 0.25f * body.tets[size_t(t)].restVolume;
        }
    body.nodeStar.clear();
    body.nodeNear.clear();
    body.nodeSlot.clear();
    for (size_t k = 0; k < n; ++k) {
        SoftNode& node = body.nodes[k];
        node.starBegin = int(body.nodeStar.size());
        body.nodeStar.insert(body.nodeStar.end(), star[k].begin(), star[k].end());
        node.starEnd = int(body.nodeStar.size());
        node.nearBegin = int(body.nodeNear.size());
        for (int t : star[k])
            for (int v : body.tets[size_t(t)].v) {
                auto at = std::find(body.nodeNear.begin() + node.nearBegin, body.nodeNear.end(), v);
                if (at == body.nodeNear.end()) at = body.nodeNear.insert(body.nodeNear.end(), v);
                body.nodeSlot.push_back(uint8_t(at - (body.nodeNear.begin() + node.nearBegin)));
            }
        node.nearEnd = int(body.nodeNear.size());
    }
    body.nodeGradient.assign(body.nodeNear.size(), Vector3(0.0f));
    std::vector<std::vector<int>> taken(n); // colours of the nodes whose star holds the particle
    std::vector<int> colour(n, 0);
    int colours = 0;
    for (size_t k = 0; k < n; ++k) {
        if (star[k].empty()) continue;
        std::vector<int> particles;
        for (int t : star[k])
            for (int v : body.tets[size_t(t)].v) particles.push_back(v - first);
        int c = 0;
        for (;; ++c) {
            bool used = false;
            for (int q : particles) used = used || std::find(taken[size_t(q)].begin(), taken[size_t(q)].end(), c) != taken[size_t(q)].end();
            if (!used) break;
        }
        colour[k] = c;
        colours = std::max(colours, c + 1);
        for (int q : particles) taken[size_t(q)].push_back(c);
    }
    body.nodeColourStart.assign(size_t(colours) + 1, 0);
    for (size_t k = 0; k < n; ++k)
        if (!star[k].empty()) ++body.nodeColourStart[size_t(colour[k]) + 1];
    for (int c = 0; c < colours; ++c) body.nodeColourStart[size_t(c) + 1] += body.nodeColourStart[size_t(c)];
    body.nodeOrder.assign(size_t(body.nodeColourStart.back()), 0);
    std::vector<int> slot(body.nodeColourStart.begin(), body.nodeColourStart.end() - 1);
    for (size_t k = 0; k < n; ++k)
        if (!star[k].empty()) body.nodeOrder[size_t(slot[size_t(colour[k])]++)] = int(k);
}

Matrix3x3 deformationGradient(const SoftTet& t, const std::vector<Vector3>& x) {
    const Vector3& x0 = x[size_t(t.v[0])];
    return Matrix3x3::fromColumns(x[size_t(t.v[1])] - x0, x[size_t(t.v[2])] - x0, x[size_t(t.v[3])] - x0) * t.restInverse;
}

Matrix3x3 deformationGradient(const SoftTet& t, const std::vector<Vector3>& base, const std::vector<Vector3>& u) {
    auto edge = [&](int c) {
        const size_t a = size_t(t.v[size_t(c)]), o = size_t(t.v[0]);
        return (base[a] - base[o]) + (u[a] - u[o]);
    };
    return Matrix3x3::fromColumns(edge(1), edge(2), edge(3)) * t.restInverse;
}

// Barycentric weights of x in a tetrahedron at rest (of corners 1, 2, 3; corner 0 takes the rest)
// and how far outside it lies: the most negative of the four weights, 0 inside.
static Vector3 restWeights(const SoftTet& t, const std::vector<Vector3>& rest, const Vector3& x, float& outside) {
    const Vector3 w = t.restInverse * (x - rest[size_t(t.v[0])]);
    outside = std::max({0.0f, -w.x, -w.y, -w.z, w.x + w.y + w.z - 1.0f});
    return w;
}

void bindSurface(SoftBody& body, const TriMesh& restSurface, const std::vector<Vector3>& rest) {
    // A vertex rides in the tetrahedron around it. Outside the lattice (the particles sit half a
    // spacing inside the mesh) it takes the tetrahedron it is least outside of among those at its
    // nearest particle: the extrapolation then stays within a spacing.
    body.surface = restSurface;
    const size_t vertices = restSurface.positions.size();
    body.vertexTet.assign(vertices, 0);
    body.vertexWeights.assign(vertices, Vector3(0.0f));
    if (body.tets.empty()) return;
    const int first = body.particles.front();
    std::vector<std::vector<int>> tetsAt(body.particles.size());
    for (int t = 0; t < int(body.tets.size()); ++t)
        for (int v : body.tets[size_t(t)].v) tetsAt[size_t(v - first)].push_back(t);
    for (size_t k = 0; k < vertices; ++k) {
        const Vector3& x = restSurface.positions[k];
        int nearest = first;
        float best = kInf;
        for (int i : body.particles)
            if (!tetsAt[size_t(i - first)].empty() && length2(rest[size_t(i)] - x) < best) best = length2(rest[size_t(i)] - x), nearest = i;
        float leastOutside = kInf;
        for (int t : tetsAt[size_t(nearest - first)]) {
            float outside = 0;
            const Vector3 w = restWeights(body.tets[size_t(t)], rest, x, outside);
            if (outside < leastOutside) leastOutside = outside, body.vertexTet[k] = t, body.vertexWeights[k] = w;
        }
    }
}

void skinSurface(const SoftBody& body, const std::vector<Vector3>& x, std::vector<Vector3>& out) {
    out.resize(body.surface.positions.size());
    if (body.tets.empty()) {
        out = body.surface.positions;
        return;
    }
    for (size_t k = 0; k < out.size(); ++k) {
        const SoftTet& t = body.tets[size_t(body.vertexTet[k])];
        const Vector3& w = body.vertexWeights[k];
        const Vector3& x0 = x[size_t(t.v[0])];
        out[k] = x0 + (x[size_t(t.v[1])] - x0) * w.x + (x[size_t(t.v[2])] - x0) * w.y + (x[size_t(t.v[3])] - x0) * w.z;
    }
}

// cof(F) = det(F) F^-T, column by column: [f2 x f3, f3 x f1, f1 x f2]. It carries area normals.
static Matrix3x3 cofactor(const Matrix3x3& F) {
    const Vector3 f1 = F.col(0), f2 = F.col(1), f3 = F.col(2);
    return Matrix3x3::fromColumns(cross(f2, f3), cross(f3, f1), cross(f1, f2));
}

void turnSurfaceNormals(const SoftBody& body, const std::vector<Vector3>& x, const std::vector<Vector3>& restNormal,
                        std::vector<Vector3>& normal) {
    for (int i : body.particles) normal[size_t(i)] = Vector3(0.0f);
    for (const SoftTet& t : body.tets) {
        const Matrix3x3 C = cofactor(deformationGradient(t, x)) * t.restVolume;
        for (int i : t.v) normal[size_t(i)] += C * restNormal[size_t(i)];
    }
    for (int i : body.particles) {
        Vector3& n = normal[size_t(i)];
        n = length2(n) > 1e-30f ? normalize(n) : restNormal[size_t(i)];
    }
}

} // namespace rf
