// Soft bodies by shape matching (Müller et al. 2005) on overlapping clusters of particles: every
// cluster pulls its particles towards its best-fit rigid pose; the skin is a mesh bound to the
// clusters by smooth weights. The data layout is in SoftBody.h.
#include "particles/SoftBody.h"

namespace rf {

std::vector<SoftCluster> buildClusters(const std::vector<int>& ids, const std::vector<Vector3>& rest, float spacing,
                                       float radius) {
    AABB box;
    for (const Vector3& r : rest) box.expand(r);
    std::vector<SoftCluster> clusters;
    const Vector3 e = box.extent();
    const int nx = std::max(1, int(std::ceil(e.x / spacing))), ny = std::max(1, int(std::ceil(e.y / spacing))),
              nz = std::max(1, int(std::ceil(e.z / spacing)));
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                // Centre of this lattice cell of the body's bounding box.
                Vector3 c = box.lo + Vector3((i + 0.5f) * e.x / nx, (j + 0.5f) * e.y / ny, (k + 0.5f) * e.z / nz);
                std::vector<size_t> members; // local indices within the radius
                for (size_t m = 0; m < ids.size(); ++m)
                    if (length2(rest[m] - c) <= radius * radius) members.push_back(m);
                if (members.size() < 4) continue; // too few particles to define a rotation
                Vector3 com(0.0f);
                for (size_t m : members) com += rest[m];
                com /= float(members.size());
                SoftCluster cl;
                cl.restCentre = cl.centre = com;
                for (size_t m : members) {
                    cl.particles.push_back(ids[m]);
                    cl.restOffsets.push_back(rest[m] - com);
                }
                clusters.push_back(std::move(cl));
            }
    if (clusters.empty()) { // small body: one cluster with everything
        SoftCluster cl;
        Vector3 com(0.0f);
        for (const Vector3& r : rest) com += r;
        com /= float(rest.size());
        cl.restCentre = cl.centre = com;
        cl.particles = ids;
        for (const Vector3& r : rest) cl.restOffsets.push_back(r - com);
        clusters.push_back(std::move(cl));
    }
    return clusters;
}

void bindSurface(SoftBody& body, const TriMesh& restSurface, const std::vector<Vector3>& particleRest) {
    // Every vertex follows the clusters whose rest centre lies within 1.5 cluster radii, with the
    // weight (1 - d / R)^2 - a smooth blend, so neighbouring vertices bound to different clusters
    // do not tear the skin apart where the clusters rotate differently. A vertex with no cluster
    // in reach (a thin spike of the mesh) takes the nearest particle's clusters.
    body.surface = restSurface;
    body.vertexClusters.assign(restSurface.positions.size(), {});
    body.vertexWeights.assign(restSurface.positions.size(), {});
    const int first = body.particles.front();
    std::vector<std::vector<int>> clustersOf(body.particles.size());
    for (int k = 0; k < int(body.clusters.size()); ++k)
        for (int i : body.clusters[k].particles) clustersOf[i - first].push_back(k);
    const float R = 1.5f * std::max(body.clusterRadius, 1e-6f);
    for (size_t v = 0; v < restSurface.positions.size(); ++v) {
        const Vector3& x = restSurface.positions[v];
        for (int k = 0; k < int(body.clusters.size()); ++k) {
            const float d = length(body.clusters[k].restCentre - x);
            if (d >= R) continue;
            body.vertexClusters[v].push_back(k);
            body.vertexWeights[v].push_back(sqr(1.0f - d / R));
        }
        if (body.vertexClusters[v].empty()) {
            size_t nearest = 0;
            float best = kInf;
            for (size_t m = 0; m < particleRest.size(); ++m) {
                const float d = length2(particleRest[m] - x);
                if (d < best) { best = d; nearest = m; }
            }
            body.vertexClusters[v] = clustersOf[nearest];
            body.vertexWeights[v].assign(body.vertexClusters[v].size(), 1.0f);
        }
        float total = 0;
        for (float w : body.vertexWeights[v]) total += w;
        for (float& w : body.vertexWeights[v]) w /= std::max(total, 1e-12f);
    }
}

void skinSurface(const SoftBody& body, std::vector<Vector3>& out) {
    out.resize(body.surface.positions.size());
    for (size_t v = 0; v < out.size(); ++v) {
        const Vector3& rest = body.surface.positions[v];
        Vector3 sum(0.0f);
        for (size_t n = 0; n < body.vertexClusters[v].size(); ++n) {
            const SoftCluster& cl = body.clusters[body.vertexClusters[v][n]];
            sum += (cl.centre + cl.rotation.rotate(rest - cl.restCentre)) * body.vertexWeights[v][n];
        }
        out[v] = body.vertexClusters[v].empty() ? rest : sum;
    }
}

void solveShapeMatching(std::vector<SoftBody>& bodies, std::vector<Vector3>& p, const std::vector<float>& invMass, int passesPerStep) {
    std::vector<Vector3> goalSum;
    std::vector<int> goalCount;
    for (SoftBody& body : bodies) {
        // The pass's share of the substep's stiffness: n passes of k' leave (1 - k') ^ n = 1 - k.
        const float k = std::clamp(body.stiffness, 0.0f, 1.0f);
        const float kPass = k >= 1.0f ? 1.0f : 1.0f - std::pow(1.0f - k, 1.0f / float(std::max(1, passesPerStep)));
        // Goals of every particle, summed over the clusters it belongs to.
        goalSum.assign(body.particles.size(), Vector3(0.0f));
        goalCount.assign(body.particles.size(), 0);
        // Map global particle index -> local slot (particles of a body are contiguous).
        const int first = body.particles.front();
        for (SoftCluster& cl : body.clusters) {
            // Current centre of mass (all particles of a body have the same mass).
            Vector3 c(0.0f);
            for (int i : cl.particles) c += p[i];
            c /= float(cl.particles.size());
            // A = sum (p - c) q^T: the deformation of the cluster from rest.
            Matrix3x3 A = Matrix3x3::zero();
            for (size_t m = 0; m < cl.particles.size(); ++m) A += Matrix3x3::outer(p[cl.particles[m]] - c, cl.restOffsets[m]);
            cl.rotation = extractRotation(A, cl.rotation, 10);
            cl.centre = c;
            const Matrix3x3 R = cl.rotation.toMatrix3x3();
            for (size_t m = 0; m < cl.particles.size(); ++m) {
                int slot = cl.particles[m] - first;
                goalSum[slot] += c + R * cl.restOffsets[m];
                goalCount[slot] += 1;
            }
        }
        // Corrections towards the averaged goals. Particles belong to different numbers of clusters,
        // so these do not sum to zero on their own: remove their mean - an internal force must not
        // move the body's centre of mass (momentum conservation). Unless part of the body is held
        // (pinned / grabbed particles): that is an external support the body must follow.
        std::vector<Vector3> delta(body.particles.size(), Vector3(0.0f));
        Vector3 mean(0.0f);
        int movable = 0;
        for (size_t s = 0; s < body.particles.size(); ++s) {
            int i = body.particles[s];
            if (goalCount[s] == 0 || invMass[i] == 0) continue;
            delta[s] = (goalSum[s] / float(goalCount[s]) - p[i]) * kPass;
            mean += delta[s];
            ++movable;
        }
        const bool held = movable < int(body.particles.size());
        if (movable > 0 && !held) mean /= float(movable);
        else mean = Vector3(0.0f);
        for (size_t s = 0; s < body.particles.size(); ++s) {
            int i = body.particles[s];
            if (goalCount[s] == 0 || invMass[i] == 0) continue;
            p[i] += delta[s] - mean;
        }
    }
}

} // namespace rf
