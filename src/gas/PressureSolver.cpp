// Pressure projection of GasSolver: the Poisson equation div(w grad p) = div u* / dt on the
// MAC grid, solved by the preconditioned conjugate gradient - with a multigrid V-cycle as the
// preconditioner (MGPCG, Multigrid.h) or the plain Jacobi one - with the face weights of the
// Boris correction and the mean divergence removed per closed region (Bridson, ch. 5).
#include "gas/GasSolver.h"

#include "core/Parallel.h"
#include "core/Probe.h"

#include <chrono>

namespace rf {

// ---------------------------------------------------------------------------
// Pressure projection (PCG)
// ---------------------------------------------------------------------------
// Chorin's projection (Chorin 1968; Bridson 2015, ch. 5): the flow after advection, forces and
// diffusion has divergence; a pressure p is found so that u - dt/rho grad p has none, from the
// Poisson equation div(grad p) = rho/dt div u*. The steps: the right-hand side (the divergence),
// its mean removed per closed region (else the equation has no solution), the matrix diagonal
// with the face weights, the conjugate gradient solve, the velocity update, and p from q.
void GasSolver::project(float dt) {
    const size_t n = size_t(nx_) * ny_ * nz_;
    const float toQ = dt / (params.fluidDensity * dx_); // q = p dt / (rho dx): the solve's unknown
    buildPressureRightHandSide(dt);
    removeMeanDivergence();
    // Face weights of the pressure equation: 1, or 1 / (1 + v_A^2/c^2) where the Boris correction
    // makes the plasma "heavier" in a strong field - then the pressure, like the Lorentz force,
    // accelerates it by force / (rho (1 + v_A^2/c^2)): a variable-density projection,
    // div( w grad p ) = div u* (Bridson, ch. 5).
    const bool weighted = magnetic.enabled && magnetic.speedLimit > 0;
    if (weighted) magnetic.borisWeights(weightU_, weightV_, weightW_, params.fluidDensity);
    buildWeightedDiagonal(weighted);
    if (params.pressurePreconditioner != PressurePreconditioner::Multigrid) {
        solvePressurePcg(weighted);
    } else if (!solvePressureMultigridPcg(weighted)) {
        Probe::add("gas/multigrid fallbacks", 1); // counted, so a test or the panel sees it happen
        solvePressurePcg(weighted);
    }
    subtractPressureGradient(weighted);
    const float toP = 1.0f / toQ;
    parallelFor(int(long(n)), [&](int c_) {
        long c = c_; p_.d[c] = float(q_[c]) * toP;
    }, 4096);
}

// The weight of a face of the pressure equation: 1, or the Boris weight of the magnetic field.
double GasSolver::faceWeightU(int i, int j, int k, bool weighted) const { return weighted ? double(weightU_.at(i, j, k)) : 1.0; }
double GasSolver::faceWeightV(int i, int j, int k, bool weighted) const { return weighted ? double(weightV_.at(i, j, k)) : 1.0; }
double GasSolver::faceWeightW(int i, int j, int k, bool weighted) const { return weighted ? double(weightW_.at(i, j, k)) : 1.0; }

// The right-hand side b of every fluid cell: minus the divergence of the intermediate velocity
// (the sum of the six face velocities), plus the expansion of burning gas; the unknown q starts
// from last step's pressure (warm start).
void GasSolver::buildPressureRightHandSide(float dt) {
    const int NX = nx_, NY = ny_, NZ = nz_;
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
}

// Compatibility: a closed region of fluid (walls all around, no outflow) can only be solved if
// its net divergence is zero - the mean is removed from every cell of each closed region.
void GasSolver::removeMeanDivergence() {
    const size_t n = size_t(nx_) * ny_ * nz_;
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

// Diagonal of the matrix: the weights of the cell's faces towards gas or towards an open
// (outflow) side; unweighted, it is the count of such faces already in diag_.
void GasSolver::buildWeightedDiagonal(bool weighted) {
    const size_t n = size_t(nx_) * ny_ * nz_;
    const int NX = nx_, NY = ny_, NZ = nz_;
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
                diagW_[c] = side(i > 0, c - 1, params.bc[0], faceWeightU(i, j, k, weighted)) +
                            side(i < NX - 1, c + 1, params.bc[1], faceWeightU(i + 1, j, k, weighted)) +
                            side(j > 0, c - NX, params.bc[2], faceWeightV(i, j, k, weighted)) +
                            side(j < NY - 1, c + NX, params.bc[3], faceWeightV(i, j + 1, k, weighted)) +
                            side(k > 0, c - sl, params.bc[4], faceWeightW(i, j, k, weighted)) +
                            side(k < NZ - 1, c + sl, params.bc[5], faceWeightW(i, j, k + 1, weighted));
            }
    }, 1);
}

// The matrix of the pressure equation applied to a vector, face by face: the weighted 7-point
// Laplacian, with solid neighbours and non-fluid cells left out. The matrix is never stored.
void GasSolver::applyPressureMatrix(const std::vector<double>& x, std::vector<double>& out, bool weighted) const {
    const int NX = nx_, NY = ny_, NZ = nz_;
    parallelFor(int(NZ), [&](int k_) {
        int k = k_;
        for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i) {
                size_t c = cidx(i, j, k);
                if (solid_[c] || diag_[c] == 0) { out[c] = 0; continue; }
                double s = diagW_[c] * x[c];
                if (i > 0 && !solid_[c - 1]) s -= faceWeightU(i, j, k, weighted) * x[c - 1];
                if (i < NX - 1 && !solid_[c + 1]) s -= faceWeightU(i + 1, j, k, weighted) * x[c + 1];
                if (j > 0 && !solid_[c - NX]) s -= faceWeightV(i, j, k, weighted) * x[c - NX];
                if (j < NY - 1 && !solid_[c + NX]) s -= faceWeightV(i, j + 1, k, weighted) * x[c + NX];
                size_t sl = size_t(NX) * NY;
                if (k > 0 && !solid_[c - sl]) s -= faceWeightW(i, j, k, weighted) * x[c - sl];
                if (k < NZ - 1 && !solid_[c + sl]) s -= faceWeightW(i, j, k + 1, weighted) * x[c + sl];
                out[c] = s;
            }
    }, 1);
}

// The dot product of two vectors of the grid, summed in parallel blocks.
static double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    return parallelSum<double>(int(a.size()), [&](int c0, int c1) {
        double s = 0;
        for (int c = c0; c < c1; ++c) s += a[c] * b[c];
        return s;
    }, 4096);
}

// Preconditioned conjugate gradient (Jacobi preconditioner) on the 7-point Laplacian, in double:
// the matrix is never stored, applyA computes A x face by face. Stops at the relative tolerance
// of the residual or at the iteration limit; both are reported to the panel.
void GasSolver::solvePressurePcg(bool weighted) {
    const size_t n = size_t(nx_) * ny_ * nz_;
    auto applyA = [&](const std::vector<double>& x, std::vector<double>& out) { applyPressureMatrix(x, out, weighted); };
    auto dotp = [](const std::vector<double>& a, const std::vector<double>& b) { return dotProduct(a, b); };
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
}

// The finest level of the multigrid: the weight of every face of the pressure equation - w between
// two gas cells, w on an open (outflow) side of the domain, 0 at a wall or next to a solid - the
// same faces applyPressureMatrix() and buildWeightedDiagonal() use, so both see one matrix.
void GasSolver::fillMultigridFaces(bool weighted) {
    multigrid_.resize(nx_, ny_, nz_);
    const int NX = nx_, NY = ny_, NZ = nz_;
    const BoundaryType* bc = params.bc;
    auto gas = [&](int i, int j, int k) { return !solid_[cidx(i, j, k)] && diag_[cidx(i, j, k)] != 0; };
    // A face between cell a (before it) and cell b (after it); index -1 or n: outside the domain.
    auto open = [&](bool aInside, bool bInside, bool aGas, bool bGas, BoundaryType before, BoundaryType after) {
        if (!aInside) return bGas && before == BoundaryType::Outflow;
        if (!bInside) return aGas && after == BoundaryType::Outflow;
        return aGas && bGas;
    };
    parallelFor(NZ + 1, [&](int k) {
        for (int j = 0; j <= NY; ++j)
            for (int i = 0; i <= NX; ++i) {
                if (j < NY && k < NZ) {
                    const bool o = open(i > 0, i < NX, i > 0 && gas(i - 1, j, k), i < NX && gas(i, j, k), bc[0], bc[1]);
                    multigrid_.faceX(i, j, k) = o ? float(faceWeightU(i, j, k, weighted)) : 0.0f;
                }
                if (i < NX && k < NZ) {
                    const bool o = open(j > 0, j < NY, j > 0 && gas(i, j - 1, k), j < NY && gas(i, j, k), bc[2], bc[3]);
                    multigrid_.faceY(i, j, k) = o ? float(faceWeightV(i, j, k, weighted)) : 0.0f;
                }
                if (i < NX && j < NY) {
                    const bool o = open(k > 0, k < NZ, k > 0 && gas(i, j, k - 1), k < NZ && gas(i, j, k), bc[4], bc[5]);
                    multigrid_.faceZ(i, j, k) = o ? float(faceWeightW(i, j, k, weighted)) : 0.0f;
                }
            }
    }, 1);
}

// The same conjugate gradient with one V-cycle of Multigrid.h as the preconditioner, z = V(r) in
// place of z = r / diag (MGPCG: McAdams, Sifakis, Teran 2010). Same stopping rule, same warm start.
// CG needs a positive definite preconditioner: if a curvature s.As or r.z comes out not positive,
// it was not for this matrix, and false hands the rest of the solve to the Jacobi PCG.
bool GasSolver::solvePressureMultigridPcg(bool weighted) {
    const int n = int(size_t(nx_) * ny_ * nz_);
    fillMultigridFaces(weighted);
    multigrid_.build();
    applyPressureMatrix(q_, As_, weighted);
    parallelFor(n, [&](int c) { r_[size_t(c)] = b_[size_t(c)] - As_[size_t(c)]; }, 4096);
    const double bnorm = std::sqrt(dotProduct(b_, b_));
    const double tol = std::max(1e-12, double(params.pressureTolerance) * bnorm);
    double rnorm = std::sqrt(dotProduct(r_, r_)), rz = 0;
    int it = 0;
    if (rnorm > tol) {
        multigrid_.apply(r_, z_);
        parallelFor(n, [&](int c) { s_[size_t(c)] = z_[size_t(c)]; }, 4096);
        rz = dotProduct(r_, z_);
        if (!(rz > 0)) return false;
    }
    while (it < params.maxPressureIterations && rnorm > tol) {
        applyPressureMatrix(s_, As_, weighted);
        const double sAs = dotProduct(s_, As_);
        if (!(sAs > 0)) return false;
        const double alpha = rz / sAs;
        const double rr = parallelSum<double>(n, [&](int c0, int c1) {
            double sum = 0;
            for (int c = c0; c < c1; ++c) {
                q_[size_t(c)] += alpha * s_[size_t(c)];
                r_[size_t(c)] -= alpha * As_[size_t(c)];
                sum += r_[size_t(c)] * r_[size_t(c)];
            }
            return sum;
        }, 4096);
        rnorm = std::sqrt(rr);
        ++it;
        if (rnorm <= tol) break;
        multigrid_.apply(r_, z_);
        const double rzNew = dotProduct(r_, z_);
        if (!(rzNew > 0)) return false;
        const double beta = rzNew / rz;
        rz = rzNew;
        parallelFor(n, [&](int c) { s_[size_t(c)] = z_[size_t(c)] + beta * s_[size_t(c)]; }, 4096);
    }
    lastIters_ = it;
    lastResidual_ = bnorm > 0 ? float(rnorm / bnorm) : 0.0f;
    return true;
}

// Velocity update: u -= w grad q on every face between two fluid cells (q already includes
// dt/(rho dx)); on an outflow side the pressure outside is zero.
void GasSolver::subtractPressureGradient(bool weighted) {
    const int NX = nx_, NY = ny_, NZ = nz_;
    auto qAt = [&](int i, int j, int k) { return q_[cidx(i, j, k)]; };
    auto fluid = [&](int i, int j, int k) { return !solid(i, j, k) && diag_[cidx(i, j, k)] != 0; };
    auto wu = [&](int i, int j, int k) { return faceWeightU(i, j, k, weighted); };
    auto wv = [&](int i, int j, int k) { return faceWeightV(i, j, k, weighted); };
    auto ww = [&](int i, int j, int k) { return faceWeightW(i, j, k, weighted); };
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
}

} // namespace rf
