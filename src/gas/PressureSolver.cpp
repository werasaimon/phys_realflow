// Pressure projection of GasSolver: the Poisson equation div(w grad p) = div u* / dt on the
// MAC grid, solved by the preconditioned conjugate gradient (Jacobi), with the face weights of
// the Boris correction and the mean divergence removed per closed region (Bridson, ch. 5).
#include "gas/GasSolver.h"

#include "core/Parallel.h"

#include <chrono>

namespace rf {

// ---------------------------------------------------------------------------
// Pressure projection (PCG)
// ---------------------------------------------------------------------------
void GasSolver::project(float dt) {
    const size_t n = size_t(nx_) * ny_ * nz_;
    const int NX = nx_, NY = ny_, NZ = nz_;
    bool anyDirichlet = false;
    for (auto t : params.bc) anyDirichlet |= (t == BoundaryType::Outflow);

    const float toQ = dt / (params.fluidDensity * dx_);
    struct Acc {
        double sum = 0;
        long count = 0;
        Acc& operator+=(const Acc& o) { sum += o.sum; count += o.count; return *this; }
    };
    Acc acc = parallelSum<Acc>(NZ, [&](int k0, int k1) {
        Acc a;
        for (int k = k0; k < k1; ++k)
            for (int j = 0; j < NY; ++j)
                for (int i = 0; i < NX; ++i) {
                    size_t c = cidx(i, j, k);
                    if (solid_[c] || diag_[c] == 0) { b_[c] = 0; q_[c] = 0; continue; }
                    double div = double(u_.at(i + 1, j, k)) - u_.at(i, j, k) + v_.at(i, j + 1, k) - v_.at(i, j, k) +
                                 w_.at(i, j, k + 1) - w_.at(i, j, k);
                    // Target divergence: 0, or the expansion of burning gas (div u = e -> face sum e dx).
                    b_[c] = -div + double(expansion_[c]) * dx_;
                    q_[c] = double(p_.d[c]) * toQ; // warm start
                    a.sum += b_[c];
                    ++a.count;
                }
        return a;
    }, 1);
    (void)acc;
    (void)anyDirichlet;
    // Compatibility: zero net divergence in every closed fluid region.
    {
        const size_t nr = regionOpen_.size();
        std::vector<double> sum(nr, 0.0);
        std::vector<long> cnt(nr, 0);
        for (size_t c = 0; c < n; ++c) {
            int r = region_[c];
            if (r < 0 || regionOpen_[r]) continue;
            sum[r] += b_[c];
            ++cnt[r];
        }
        bool any = false;
        for (size_t r = 0; r < nr; ++r)
            if (cnt[r] > 0) { sum[r] /= double(cnt[r]); any = true; }
        regionMeanB_ = sum;
        if (any)
            parallelFor(int(long(n)), [&](int c) {
                int r = region_[c];
                if (r >= 0 && !regionOpen_[r]) b_[c] -= sum[r];
            }, 4096);
    }

    // Face weights of the pressure equation: 1, or 1 / (1 + v_A^2/c^2) where the Boris correction
    // makes the plasma "heavier" in a strong field - then the pressure, like the Lorentz force,
    // accelerates it by force / (rho (1 + v_A^2/c^2)): a variable-density projection,
    // div( w grad p ) = div u* (Bridson, ch. 5).
    const bool weighted = magnetic.enabled && magnetic.speedLimit > 0;
    if (weighted) magnetic.borisWeights(weightU_, weightV_, weightW_, params.fluidDensity);
    auto wu = [&](int i, int j, int k) { return weighted ? double(weightU_.at(i, j, k)) : 1.0; };
    auto wv = [&](int i, int j, int k) { return weighted ? double(weightV_.at(i, j, k)) : 1.0; };
    auto ww = [&](int i, int j, int k) { return weighted ? double(weightW_.at(i, j, k)) : 1.0; };
    // Diagonal: the weights of the cell's faces towards gas or towards an open (outflow) side.
    diagW_.assign(n, 0.0);
    parallelFor(int(NZ), [&](int k) {
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i) {
                const size_t c = cidx(i, j, k);
                if (solid_[c] || diag_[c] == 0) continue;
                if (!weighted) { diagW_[c] = diag_[c]; continue; }
                auto side = [&](bool inside, size_t nb, BoundaryType b, double w) {
                    return (inside ? !solid_[nb] : b == BoundaryType::Outflow) ? w : 0.0;
                };
                const size_t sl = size_t(NX) * NY;
                diagW_[c] = side(i > 0, c - 1, params.bc[0], wu(i, j, k)) + side(i < NX - 1, c + 1, params.bc[1], wu(i + 1, j, k)) +
                            side(j > 0, c - NX, params.bc[2], wv(i, j, k)) + side(j < NY - 1, c + NX, params.bc[3], wv(i, j + 1, k)) +
                            side(k > 0, c - sl, params.bc[4], ww(i, j, k)) + side(k < NZ - 1, c + sl, params.bc[5], ww(i, j, k + 1));
            }
    }, 1);

    auto applyA = [&](const std::vector<double>& x, std::vector<double>& out) {
        parallelFor(int(NZ), [&](int k_) {
            int k = k_;
            for (int j = 0; j < NY; ++j)
                for (int i = 0; i < NX; ++i) {
                    size_t c = cidx(i, j, k);
                    if (solid_[c] || diag_[c] == 0) { out[c] = 0; continue; }
                    double s = diagW_[c] * x[c];
                    if (i > 0 && !solid_[c - 1]) s -= wu(i, j, k) * x[c - 1];
                    if (i < NX - 1 && !solid_[c + 1]) s -= wu(i + 1, j, k) * x[c + 1];
                    if (j > 0 && !solid_[c - NX]) s -= wv(i, j, k) * x[c - NX];
                    if (j < NY - 1 && !solid_[c + NX]) s -= wv(i, j + 1, k) * x[c + NX];
                    size_t sl = size_t(NX) * NY;
                    if (k > 0 && !solid_[c - sl]) s -= ww(i, j, k) * x[c - sl];
                    if (k < NZ - 1 && !solid_[c + sl]) s -= ww(i, j, k + 1) * x[c + sl];
                    out[c] = s;
                }
        }, 1);
    };
    auto dotp = [&](const std::vector<double>& a, const std::vector<double>& b) {
        return parallelSum<double>(int(n), [&](int c0, int c1) {
            double s = 0;
            for (int c = c0; c < c1; ++c) s += a[c] * b[c];
            return s;
        }, 4096);
    };

    applyA(q_, As_);
    double rz = parallelSum<double>(int(n), [&](int c0, int c1) {
        double acc2 = 0;
        for (int c = c0; c < c1; ++c) {
            r_[c] = b_[c] - As_[c];
            z_[c] = diagW_[c] > 0 ? r_[c] / diagW_[c] : 0.0;
            s_[c] = z_[c];
            acc2 += r_[c] * z_[c];
        }
        return acc2;
    }, 4096);
    const double bnorm = std::sqrt(dotp(b_, b_));
    const double tol = std::max(1e-12, double(params.pressureTolerance) * bnorm);
    double rnorm = std::sqrt(dotp(r_, r_));
    int it = 0;
    for (; it < params.maxPressureIterations && rnorm > tol; ++it) {
        applyA(s_, As_);
        double sAs = dotp(s_, As_);
        if (std::fabs(sAs) < 1e-30) break;
        double alpha = rz / sAs;
        struct Two {
            double a = 0, b = 0;
            Two& operator+=(const Two& o) { a += o.a; b += o.b; return *this; }
        };
        Two t = parallelSum<Two>(int(n), [&](int c0, int c1) {
            Two x;
            for (int c = c0; c < c1; ++c) {
                q_[c] += alpha * s_[c];
                r_[c] -= alpha * As_[c];
                z_[c] = diagW_[c] > 0 ? r_[c] / diagW_[c] : 0.0;
                x.a += r_[c] * z_[c];
                x.b += r_[c] * r_[c];
            }
            return x;
        }, 4096);
        double rzNew = t.a, rr = t.b;
        rnorm = std::sqrt(rr);
        double beta = rzNew / rz;
        rz = rzNew;
        parallelFor(int(long(n)), [&](int c_) {
            long c = c_; s_[c] = z_[c] + beta * s_[c];
        }, 4096);
    }
    lastIters_ = it;
    lastResidual_ = bnorm > 0 ? float(rnorm / bnorm) : 0.0f;

    // Velocity update: u -= grad q (q already includes dt/(rho dx)).
    auto qAt = [&](int i, int j, int k) { return q_[cidx(i, j, k)]; };
    auto fluid = [&](int i, int j, int k) { return !solid(i, j, k) && diag_[cidx(i, j, k)] != 0; };
    const BoundaryType* bc = params.bc;
    parallelFor(int(NZ), [&](int k_) {
        int k = k_;
        for (int j = 0; j < NY; ++j) {
            for (int i = 1; i < NX; ++i)
                if (fluid(i - 1, j, k) && fluid(i, j, k)) u_.at(i, j, k) -= float(wu(i, j, k) * (qAt(i, j, k) - qAt(i - 1, j, k)));
            if (bc[0] == BoundaryType::Outflow && fluid(0, j, k)) u_.at(0, j, k) -= float(wu(0, j, k) * qAt(0, j, k));
            if (bc[1] == BoundaryType::Outflow && fluid(NX - 1, j, k)) u_.at(NX, j, k) += float(wu(NX, j, k) * qAt(NX - 1, j, k));
        }
    }, 1);
    parallelFor(int(NZ), [&](int k_) {
        int k = k_;
        for (int i = 0; i < NX; ++i) {
            for (int j = 1; j < NY; ++j)
                if (fluid(i, j - 1, k) && fluid(i, j, k)) v_.at(i, j, k) -= float(wv(i, j, k) * (qAt(i, j, k) - qAt(i, j - 1, k)));
            if (bc[2] == BoundaryType::Outflow && fluid(i, 0, k)) v_.at(i, 0, k) -= float(wv(i, 0, k) * qAt(i, 0, k));
            if (bc[3] == BoundaryType::Outflow && fluid(i, NY - 1, k)) v_.at(i, NY, k) += float(wv(i, NY, k) * qAt(i, NY - 1, k));
        }
    }, 1);
    parallelFor(int(NY), [&](int j_) {
        int j = j_;
        for (int i = 0; i < NX; ++i) {
            for (int k = 1; k < NZ; ++k)
                if (fluid(i, j, k - 1) && fluid(i, j, k)) w_.at(i, j, k) -= float(ww(i, j, k) * (qAt(i, j, k) - qAt(i, j, k - 1)));
            if (bc[4] == BoundaryType::Outflow && fluid(i, j, 0)) w_.at(i, j, 0) -= float(ww(i, j, 0) * qAt(i, j, 0));
            if (bc[5] == BoundaryType::Outflow && fluid(i, j, NZ - 1)) w_.at(i, j, NZ) += float(ww(i, j, NZ) * qAt(i, j, NZ - 1));
        }
    }, 1);

    const float toP = 1.0f / toQ;
    parallelFor(int(long(n)), [&](int c_) {
        long c = c_; p_.d[c] = float(q_[c]) * toP;
    }, 4096);
}

} // namespace rf
