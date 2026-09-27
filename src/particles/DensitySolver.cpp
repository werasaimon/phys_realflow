// The liquid of ParticleSystem: the neighbour grid, the density constraint of Position Based
// Fluids (Macklin & Muller 2013) with the walls as liquid at rest (Koschier & Bender 2017), the
// position correction with the artificial pressure, XSPH viscosity and vorticity confinement.
#include "particles/ParticleSystem.h"

#include "core/Parallel.h"

namespace rf {

void ParticleSystem::buildGrid(const std::vector<Vector3>& pts) {
    const int n = int(pts.size());
    const int ncell = gx_ * gy_ * gz_;
    cellOf_.resize(n);
    cellStart_.assign(ncell + 1, 0);
    const float inv = 1.0f / h_;
    parallelFor(n, [&](int i) {
        Vector3 q = (pts[i] - domain_.lo) * inv;
        int cx = clampv(int(q.x), 0, gx_ - 1), cy = clampv(int(q.y), 0, gy_ - 1), cz = clampv(int(q.z), 0, gz_ - 1);
        cellOf_[i] = cx + gx_ * (cy + gy_ * cz);
    });
    for (int i = 0; i < n; ++i) cellStart_[cellOf_[i] + 1]++;
    for (int c = 0; c < ncell; ++c) cellStart_[c + 1] += cellStart_[c];
    sorted_.resize(n);
    std::vector<int> fill(cellStart_.begin(), cellStart_.end() - 1);
    for (int i = 0; i < n; ++i) sorted_[fill[cellOf_[i]]++] = i;
}

void ParticleSystem::findNeighbors() {
    const int n = int(p_.size());
    nbrCount_.assign(n, 0);
    nbr_.resize(size_t(n) * kMaxNeighbors);
    parallelFor(n, [&](int i) {
        const Vector3 pi = p_[i];
        int c = cellOf_[i];
        int cx = c % gx_, cy = (c / gx_) % gy_, cz = c / (gx_ * gy_);
        int cnt = 0;
        int* out = &nbr_[size_t(i) * kMaxNeighbors];
        for (int z = std::max(cz - 1, 0); z <= std::min(cz + 1, gz_ - 1); ++z)
            for (int y = std::max(cy - 1, 0); y <= std::min(cy + 1, gy_ - 1); ++y)
                for (int x = std::max(cx - 1, 0); x <= std::min(cx + 1, gx_ - 1); ++x) {
                    int cell = x + gx_ * (y + gy_ * z);
                    for (int s = cellStart_[cell]; s < cellStart_[cell + 1]; ++s) {
                        int j = sorted_[s];
                        if (j == i || cnt >= kMaxNeighbors) continue;
                        if (length2(pi - p_[j]) < h2_) out[cnt++] = j;
                    }
                }
        nbrCount_[i] = cnt;
    });
}

// Domain walls in the density (Koschier & Bender 2017, "Density Maps for Improved SPH Boundary
// Handling"): the part of a particle's kernel that lies beyond a wall counts as liquid at rest
// density. For the poly6 kernel and a plane at distance d that part is closed-form,
//   Phi(d) = pi k / 4 * integral_d^h (h^2 - z^2)^4 dz,   Phi(0) = 1/2,   Phi'(d) = -pi k / 4 (h^2 - d^2)^4,
// (k the poly6 constant). Without it a particle at a wall misses neighbours, the liquid there
// clumps, and in a corner the particles line up and get pushed along the corner edge; with it the
// wall pushes along its normal. Returns sum Phi over the six walls; gradient = d(sum Phi)/dp.
float ParticleSystem::wallVolume(const Vector3& p, Vector3& gradient) const {
    const float h = h_, h2 = h2_;
    const float c = 0.25f * kPi * poly6_;
    auto F = [&](float z) { // antiderivative of (h^2 - z^2)^4
        const float z2 = z * z;
        return z * (h2 * h2 * h2 * h2 + z2 * (-4.0f / 3.0f * h2 * h2 * h2 + z2 * (1.2f * h2 * h2 + z2 * (-4.0f / 7.0f * h2 + z2 / 9.0f))));
    };
    const float Fh = F(h);
    float phi = 0;
    gradient = Vector3(0.0f);
    for (int a = 0; a < 3; ++a)
        for (int side = 0; side < 2; ++side) {
            const float d = side == 0 ? p[a] - domain_.lo[a] : domain_.hi[a] - p[a];
            if (d >= h) continue;
            const float dc = std::max(d, 0.0f);
            const float q = h2 - dc * dc;
            phi += c * (Fh - F(dc));
            gradient[a] += (side == 0 ? -1.0f : 1.0f) * c * q * q * q * q; // Phi'(d) dd/dp
        }
    return phi;
}

void ParticleSystem::computeLambda() {
    const int n = int(p_.size());
    const float invRho0 = 1.0f / params.restDensity;
    const float eps = params.relaxation / h2_;
    double errSum = parallelSum<double>(n, [&](int b, int e) {
        double err = 0;
        for (int i = b; i < e; ++i) {
            if (!isFluid(i)) {
                rho_[i] = params.restDensity;
                lambda_[i] = 0;
                continue;
            }
            const Vector3 pi = p_[i];
            float rho = mass_ * W(0);
            Vector3 gi(0.0f);
            float sum2 = 0;
            const int* nb = &nbr_[size_t(i) * kMaxNeighbors];
            for (int k = 0; k < nbrCount_[i]; ++k) {
                // Solid particles count too: they occupy the same volume as a fluid particle, so
                // the liquid is pushed out of them - and pushes back on them (buoyancy).
                Vector3 r = pi - p_[nb[k]];
                const float mj = mass_ * volume_[nb[k]];
                rho += mj * W(length2(r));
                Vector3 g = gradW(r) * (mj * invRho0);
                // Generalized masses (Macklin et al. 2014): a neighbour moves by its inverse mass
                // relative to a fluid particle's (computeDeltaP), so it weighs that much here - a
                // light cloth is not pushed 100x further than the constraint needs, a pinned one
                // not at all.
                sum2 += length2(g) * (invMass_[nb[k]] * mass_);
                gi += g;
            }
            Vector3 gWall;
            rho += params.restDensity * wallVolume(pi, gWall); // the walls as liquid at rest
            gi += gWall;
            rho_[i] = rho;
            float C = std::max(rho * invRho0 - 1.0f, 0.0f);
            err += C;
            lambda_[i] = -C / (sum2 + length2(gi) + eps);
        }
        return err;
    });
    avgDensityError_ = fluidCount_ > 0 ? float(errSum / double(fluidCount_)) : 0.0f;
}

void ParticleSystem::computeDeltaP() {
    const int n = int(p_.size());
    const float k = mass_ / params.restDensity;
    // Artificial pressure s_corr = -k (W(r) / W(dq))^4 (Macklin & Muller 2013, eq. 13) is added to
    // lambda_i + lambda_j. Their lambda is dimensionless (h = 1); ours has the units of h^2 (the
    // gradients of C go as 1/h), so k is scaled by h^2 - otherwise the same k puffs small
    // particles up far more than large ones (a resting column of 5 mm particles grew 80 %).
    const float tk = params.tensileK * h2_;
    const float invDq = deltaQW_ > 0 ? 1.0f / deltaQW_ : 0.0f;
    parallelFor(n, [&](int i) {
        const Vector3 pi = p_[i];
        Vector3 d(0.0f);
        const int* nb = &nbr_[size_t(i) * kMaxNeighbors];
        if (!isFluid(i)) {
            // A solid particle takes its share of the neighbouring fluid particles' pressure
            // corrections, scaled by its inverse mass relative to a fluid particle (two-way).
            if (invMass_[i] == 0) {
                dp_[i] = d;
                return;
            }
            for (int m = 0; m < nbrCount_[i]; ++m) {
                int j = nb[m];
                if (isFluid(j)) d += gradW(pi - p_[j]) * lambda_[j];
            }
            dp_[i] = d * (k * invMass_[i] * mass_ * volume_[i]);
            return;
        }
        for (int m = 0; m < nbrCount_[i]; ++m) {
            int j = nb[m];
            Vector3 r = pi - p_[j];
            float q = W(length2(r)) * invDq;
            float q2 = q * q;
            float scorr = -tk * q2 * q2;
            d += gradW(r) * ((lambda_[i] + lambda_[j] + scorr) * volume_[j]); // m_j = mass_ * volume_j, as in the density
        }
        Vector3 gWall;
        wallVolume(pi, gWall);
        dp_[i] = d * k + gWall * lambda_[i]; // the walls push along their normals
    });
}

void ParticleSystem::applyViscosityAndVorticity(float dt) {
    const int n = int(p_.size());
    const float c = params.viscosity;
    const float eps = params.vorticity;
    if (eps > 0) {
        parallelFor(n, [&](int i) {
            Vector3 w(0.0f);
            const int* nb = &nbr_[size_t(i) * kMaxNeighbors];
            for (int m = 0; m < nbrCount_[i]; ++m) {
                int j = nb[m];
                if (!isFluid(j)) continue;
                w += cross(gradW(p_[i] - p_[j]), v_[j] - v_[i]) * (mass_ / std::max(rho_[j], 1e-3f));
            }
            omega_[i] = w;
        });
    }
    parallelFor(n, [&](int i) {
        Vector3 xs(0.0f), eta(0.0f);
        if (!isFluid(i)) {
            vtmp_[i] = v_[i];
            return;
        }
        const float wi = eps > 0 ? length(omega_[i]) : 0.0f;
        const int* nb = &nbr_[size_t(i) * kMaxNeighbors];
        for (int m = 0; m < nbrCount_[i]; ++m) {
            int j = nb[m];
            if (!isFluid(j)) continue;
            Vector3 r = p_[i] - p_[j];
            float vol = mass_ / std::max(rho_[j], 1e-3f);
            xs += (v_[j] - v_[i]) * (W(length2(r)) * vol);
            if (eps > 0) eta += gradW(r) * ((length(omega_[j]) - wi) * vol);
        }
        Vector3 vn = v_[i] + xs * c;
        if (eps > 0 && length2(eta) > 1e-12f) vn += cross(normalize(eta), omega_[i]) * (eps * dt);
        vtmp_[i] = vn;
    });
    v_.swap(vtmp_);
}

} // namespace rf
