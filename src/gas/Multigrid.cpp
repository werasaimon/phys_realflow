// PressureMultigrid (Multigrid.h): the levels, the smoother, the transfers between levels and the
// V-cycle. Every loop runs over rows of cells (j, k) in parallel on the SDK's thread pool; nothing
// is allocated after resize(), so a step of the gas allocates nothing here.
#include "gas/Multigrid.h"

#include "core/Parallel.h"

#include <algorithm>

namespace rf {

namespace {

constexpr size_t kCoarsestCells = 512; // a level this small is the last one
constexpr int kCoarsestSweeps = 20;    // it is solved by this many sweeps forward, then backward
constexpr int kBandSweeps = 2;         // extra sweeps of the band around solids, before and after

// Rows of cells per parallel chunk: about 4096 cells, so a small level runs on one thread.
int rowGrain(int nx) { return std::max(1, 4096 / std::max(nx, 1)); }

// Interpolation along one axis: the coarse cells that fine cell i takes its value from, with their
// weights. Factor 1 (the axis was not coarsened): the same index, weight 1. Factor 2: the parent
// i / 2 with 3/4 and the parent's neighbour on i's side with 1/4 - linear interpolation between
// cell centres (the fine centre lies a quarter of a coarse cell from its parent's centre).
int axisInterpolation(int i, int factor, int nCoarse, int idx[2], double w[2]) {
    if (factor == 1) { idx[0] = i; w[0] = 1.0; return 1; }
    const int parent = i / 2, side = (i % 2 == 0) ? parent - 1 : parent + 1;
    idx[0] = parent;
    w[0] = 0.75;
    if (side < 0 || side >= nCoarse) return 1;
    idx[1] = side;
    w[1] = 0.25;
    return 2;
}

// The transpose of axisInterpolation: the fine cells that coarse cell I gathers from along one
// axis, with the same weights (fine cells 2I-1, 2I, 2I+1, 2I+2 with 1/4, 3/4, 3/4, 1/4).
int axisGather(int I, int factor, int nFine, int idx[4], double w[4]) {
    if (factor == 1) { idx[0] = I; w[0] = 1.0; return 1; }
    const double pattern[4] = {0.25, 0.75, 0.75, 0.25};
    int n = 0;
    for (int a = 0; a < 4; ++a) {
        const int i = 2 * I - 1 + a;
        if (i < 0 || i >= nFine) continue;
        idx[n] = i;
        w[n] = pattern[a];
        ++n;
    }
    return n;
}

} // namespace

// The hierarchy for a fine grid of nx x ny x nz: every axis with at least 4 cells is halved (the
// last cell of an odd axis has one child along it) until a level has kCoarsestCells or fewer, or
// no axis can be halved any more.
void PressureMultigrid::resize(int nx, int ny, int nz) {
    if (!levels_.empty() && levels_[0].nx == nx && levels_[0].ny == ny && levels_[0].nz == nz) return;
    levels_.clear();
    Level fine;
    fine.nx = nx;
    fine.ny = ny;
    fine.nz = nz;
    levels_.push_back(fine);
    while (levels_.back().cells() > kCoarsestCells) {
        const Level& f = levels_.back();
        Level c;
        c.fx = f.nx >= 4 ? 2 : 1;
        c.fy = f.ny >= 4 ? 2 : 1;
        c.fz = f.nz >= 4 ? 2 : 1;
        if (c.fx == 1 && c.fy == 1 && c.fz == 1) break;
        c.nx = (f.nx + c.fx - 1) / c.fx;
        c.ny = (f.ny + c.fy - 1) / c.fy;
        c.nz = (f.nz + c.fz - 1) / c.fz;
        levels_.push_back(c);
    }
    for (Level& L : levels_) allocate(L);
}

void PressureMultigrid::allocate(Level& L) {
    L.wx.assign(size_t(L.nx + 1) * L.ny * L.nz, 0.0f);
    L.wy.assign(size_t(L.nx) * (L.ny + 1) * L.nz, 0.0f);
    L.wz.assign(size_t(L.nx) * L.ny * (L.nz + 1), 0.0f);
    L.diag.assign(L.cells(), 0.0);
    L.invNorm.assign(L.cells(), 0.0);
    L.x.assign(L.cells(), 0.0);
    L.b.assign(L.cells(), 0.0);
    L.r.assign(L.cells(), 0.0);
    L.irregular.assign(L.cells(), 0);
    for (std::vector<uint32_t>& list : L.band) list.reserve(L.cells()); // filled every build, never grown
}

// The coarser levels from the finest one's face weights: their faces, their diagonals and the
// normalisation of the interpolation into the level above.
void PressureMultigrid::build() {
    computeDiagonal(levels_[0]);
    findBoundaryBand(levels_[0]);
    for (size_t l = 1; l < levels_.size(); ++l) {
        coarsenFaces(levels_[l - 1], levels_[l]);
        computeDiagonal(levels_[l]);
        findBoundaryBand(levels_[l]);
        computeInterpolationNorms(levels_[l - 1], levels_[l]);
    }
}

// The band that gets extra smoothing: every cell that takes part and is irregular - a face inside
// the grid closed (a solid or a left-out cell behind it) - or has an irregular neighbour. The walls
// of the domain are not irregular: every level represents them exactly. Split by colour.
void PressureMultigrid::findBoundaryBand(Level& L) {
    const size_t sy = size_t(L.nx), sz = size_t(L.nx) * L.ny;
    parallelFor(L.ny * L.nz, [&](int row) {
        const int j = row % L.ny, k = row / L.ny;
        for (int i = 0; i < L.nx; ++i) {
            const size_t c = L.cell(i, j, k), fx = L.faceX(i, j, k), fy = L.faceY(i, j, k), fz = L.faceZ(i, j, k);
            L.irregular[c] = L.diag[c] > 0 &&
                             ((i > 0 && L.wx[fx] == 0) || (i < L.nx - 1 && L.wx[fx + 1] == 0) || (j > 0 && L.wy[fy] == 0) ||
                              (j < L.ny - 1 && L.wy[fy + sy] == 0) || (k > 0 && L.wz[fz] == 0) || (k < L.nz - 1 && L.wz[fz + sz] == 0));
        }
    }, rowGrain(L.nx));
    L.band[0].clear();
    L.band[1].clear();
    for (int k = 0; k < L.nz; ++k)
        for (int j = 0; j < L.ny; ++j)
            for (int i = 0; i < L.nx; ++i) {
                const size_t c = L.cell(i, j, k);
                if (L.diag[c] <= 0) continue;
                const bool near = L.irregular[c] || (i > 0 && L.irregular[c - 1]) || (i < L.nx - 1 && L.irregular[c + 1]) ||
                                  (j > 0 && L.irregular[c - sy]) || (j < L.ny - 1 && L.irregular[c + sy]) ||
                                  (k > 0 && L.irregular[c - sz]) || (k < L.nz - 1 && L.irregular[c + sz]);
                if (near) L.band[(i + j + k) & 1].push_back(uint32_t(c));
            }
}

// The diagonal of every cell: the sum of its six face weights (an open side of the domain counts,
// a wall or a solid does not). A cell whose faces are all closed gets 0 and is left out.
void PressureMultigrid::computeDiagonal(Level& L) {
    parallelFor(L.ny * L.nz, [&](int row) {
        const int j = row % L.ny, k = row / L.ny;
        for (int i = 0; i < L.nx; ++i)
            L.diag[L.cell(i, j, k)] = double(L.wx[L.faceX(i, j, k)]) + L.wx[L.faceX(i + 1, j, k)] + L.wy[L.faceY(i, j, k)] +
                                      L.wy[L.faceY(i, j + 1, k)] + L.wz[L.faceZ(i, j, k)] + L.wz[L.faceZ(i, j, k + 1)];
    }, rowGrain(L.nx));
}

// A coarse face is the mean of the fine faces it covers (children outside the grid count as
// walls), divided by factor^2 along its normal: the Laplacian of a grid of spacing 2h has a
// quarter of the coupling of spacing h, in the unit form where the restriction is a mean.
void PressureMultigrid::coarsenFaces(const Level& F, Level& C) {
    auto meanX = [&](int i, int J, int K) {
        double s = 0;
        for (int k = K * C.fz; k < std::min(K * C.fz + C.fz, F.nz); ++k)
            for (int j = J * C.fy; j < std::min(J * C.fy + C.fy, F.ny); ++j) s += F.wx[F.faceX(i, j, k)];
        return s / (C.fy * C.fz) / (C.fx * C.fx);
    };
    auto meanY = [&](int I, int j, int K) {
        double s = 0;
        for (int k = K * C.fz; k < std::min(K * C.fz + C.fz, F.nz); ++k)
            for (int i = I * C.fx; i < std::min(I * C.fx + C.fx, F.nx); ++i) s += F.wy[F.faceY(i, j, k)];
        return s / (C.fx * C.fz) / (C.fy * C.fy);
    };
    auto meanZ = [&](int I, int J, int k) {
        double s = 0;
        for (int j = J * C.fy; j < std::min(J * C.fy + C.fy, F.ny); ++j)
            for (int i = I * C.fx; i < std::min(I * C.fx + C.fx, F.nx); ++i) s += F.wz[F.faceZ(i, j, k)];
        return s / (C.fx * C.fy) / (C.fz * C.fz);
    };
    // The fine face under coarse face I: I * factor, the domain's last face for the last one.
    parallelFor(C.nz + 1, [&](int K) {
        for (int J = 0; J <= C.ny; ++J)
            for (int I = 0; I <= C.nx; ++I) {
                if (J < C.ny && K < C.nz) C.wx[C.faceX(I, J, K)] = float(meanX(std::min(I * C.fx, F.nx), J, K));
                if (I < C.nx && K < C.nz) C.wy[C.faceY(I, J, K)] = float(meanY(I, std::min(J * C.fy, F.ny), K));
                if (I < C.nx && J < C.ny) C.wz[C.faceZ(I, J, K)] = float(meanZ(I, J, std::min(K * C.fz, F.nz)));
            }
    }, 1);
}

// For every fine cell, 1 / (the sum of its interpolation weights over the coarse cells that take
// part): near a wall or a solid the missing coarse cells are dropped and the rest renormalised,
// so a constant on the coarse level stays the same constant on the fine one.
void PressureMultigrid::computeInterpolationNorms(Level& F, const Level& C) {
    parallelFor(F.ny * F.nz, [&](int row) {
        const int j = row % F.ny, k = row / F.ny;
        int ci[2], cj[2], ck[2];
        double wi[2], wj[2], wk[2];
        const int nj = axisInterpolation(j, C.fy, C.ny, cj, wj), nk = axisInterpolation(k, C.fz, C.nz, ck, wk);
        for (int i = 0; i < F.nx; ++i) {
            const size_t f = F.cell(i, j, k);
            F.invNorm[f] = 0.0;
            if (F.diag[f] <= 0) continue;
            const int ni = axisInterpolation(i, C.fx, C.nx, ci, wi);
            double s = 0;
            for (int c = 0; c < nk; ++c)
                for (int b = 0; b < nj; ++b)
                    for (int a = 0; a < ni; ++a)
                        if (C.diag[C.cell(ci[a], cj[b], ck[c])] > 0) s += wi[a] * wj[b] * wk[c];
            F.invNorm[f] = s > 0 ? 1.0 / s : 0.0;
        }
    }, rowGrain(F.nx));
}

// Red-black Gauss-Seidel: x_c = (b_c + sum w_f x_neighbour) / diag_c, first on the cells with
// i + j + k even ("red"), then on the odd ones ("black") - or black first when redFirst is false.
// A red cell's neighbours are all black, so every cell of one colour is updated in parallel.
void PressureMultigrid::smooth(Level& L, int sweeps, bool redFirst) {
    for (int s = 0; s < sweeps; ++s)
        for (int pass = 0; pass < 2; ++pass) {
            const int colour = redFirst ? pass : 1 - pass;
            parallelFor(L.ny * L.nz, [&](int row) {
                const int j = row % L.ny, k = row / L.ny;
                for (int i = (j + k + colour) & 1; i < L.nx; i += 2) relaxCell(L, L.cell(i, j, k), i, j, k);
            }, rowGrain(L.nx));
        }
}

// The same sweeps on the band around solids only (its red cells, then its black ones, or reverse).
void PressureMultigrid::smoothBand(Level& L, int sweeps, bool redFirst) {
    for (int s = 0; s < sweeps; ++s)
        for (int pass = 0; pass < 2; ++pass) {
            const std::vector<uint32_t>& cells = L.band[redFirst ? pass : 1 - pass];
            parallelFor(int(cells.size()), [&](int n) {
                const size_t c = cells[size_t(n)];
                const int i = int(c % size_t(L.nx)), j = int((c / size_t(L.nx)) % size_t(L.ny)), k = int(c / (size_t(L.nx) * L.ny));
                relaxCell(L, c, i, j, k);
            }, 1024);
        }
}

// One Gauss-Seidel update of cell c = (i, j, k): the value that satisfies its own equation with
// the neighbours as they are now.
void PressureMultigrid::relaxCell(Level& L, size_t c, int i, int j, int k) {
    if (L.diag[c] <= 0) return;
    const size_t sy = size_t(L.nx), sz = size_t(L.nx) * L.ny;
    const size_t fx = L.faceX(i, j, k), fy = L.faceY(i, j, k), fz = L.faceZ(i, j, k);
    double v = L.b[c];
    if (i > 0) v += L.wx[fx] * L.x[c - 1];
    if (i < L.nx - 1) v += L.wx[fx + 1] * L.x[c + 1];
    if (j > 0) v += L.wy[fy] * L.x[c - sy];
    if (j < L.ny - 1) v += L.wy[fy + sy] * L.x[c + sy];
    if (k > 0) v += L.wz[fz] * L.x[c - sz];
    if (k < L.nz - 1) v += L.wz[fz + sz] * L.x[c + sz];
    L.x[c] = v / L.diag[c];
}

// r = b - A x on every cell that takes part (0 elsewhere).
void PressureMultigrid::computeResidual(Level& L) {
    const size_t sx = 1, sy = size_t(L.nx), sz = size_t(L.nx) * L.ny;
    parallelFor(L.ny * L.nz, [&](int row) {
        const int j = row % L.ny, k = row / L.ny;
        for (int i = 0; i < L.nx; ++i) {
            const size_t c = L.cell(i, j, k);
            if (L.diag[c] <= 0) { L.r[c] = 0; continue; }
            const size_t fx = L.faceX(i, j, k), fy = L.faceY(i, j, k), fz = L.faceZ(i, j, k);
            double v = L.b[c] - L.diag[c] * L.x[c];
            if (i > 0) v += L.wx[fx] * L.x[c - sx];
            if (i < L.nx - 1) v += L.wx[fx + 1] * L.x[c + sx];
            if (j > 0) v += L.wy[fy] * L.x[c - sy];
            if (j < L.ny - 1) v += L.wy[fy + sy] * L.x[c + sy];
            if (k > 0) v += L.wz[fz] * L.x[c - sz];
            if (k < L.nz - 1) v += L.wz[fz + sz] * L.x[c + sz];
            L.r[c] = v;
        }
    }, rowGrain(L.nx));
}

// The coarse right-hand side: the transpose of the interpolation applied to the fine residual,
// divided by 2 per halved axis (the interior column sums), so a smooth residual arrives as its mean.
void PressureMultigrid::restrictResidual(const Level& F, Level& C) {
    const double scale = 1.0 / double(C.fx * C.fy * C.fz);
    parallelFor(C.ny * C.nz, [&](int row) {
        const int J = row % C.ny, K = row / C.ny;
        int fi[4], fj[4], fk[4];
        double wi[4], wj[4], wk[4];
        const int nj = axisGather(J, C.fy, F.ny, fj, wj), nk = axisGather(K, C.fz, F.nz, fk, wk);
        for (int I = 0; I < C.nx; ++I) {
            const size_t c = C.cell(I, J, K);
            if (C.diag[c] <= 0) { C.b[c] = 0; continue; }
            const int ni = axisGather(I, C.fx, F.nx, fi, wi);
            double s = 0;
            for (int z = 0; z < nk; ++z)
                for (int y = 0; y < nj; ++y)
                    for (int x = 0; x < ni; ++x) {
                        const size_t f = F.cell(fi[x], fj[y], fk[z]);
                        s += wi[x] * wj[y] * wk[z] * F.invNorm[f] * F.r[f];
                    }
            C.b[c] = s * scale;
        }
    }, rowGrain(C.nx));
}

// The coarse correction interpolated back and added to the fine unknown (renormalised weights).
void PressureMultigrid::interpolateCorrection(const Level& C, Level& F) {
    parallelFor(F.ny * F.nz, [&](int row) {
        const int j = row % F.ny, k = row / F.ny;
        int ci[2], cj[2], ck[2];
        double wi[2], wj[2], wk[2];
        const int nj = axisInterpolation(j, C.fy, C.ny, cj, wj), nk = axisInterpolation(k, C.fz, C.nz, ck, wk);
        for (int i = 0; i < F.nx; ++i) {
            const size_t f = F.cell(i, j, k);
            if (F.invNorm[f] == 0) continue;
            const int ni = axisInterpolation(i, C.fx, C.nx, ci, wi);
            double s = 0;
            for (int c = 0; c < nk; ++c)
                for (int b = 0; b < nj; ++b)
                    for (int a = 0; a < ni; ++a) {
                        const size_t cc = C.cell(ci[a], cj[b], ck[c]);
                        if (C.diag[cc] > 0) s += wi[a] * wj[b] * wk[c] * C.x[cc];
                    }
            F.x[f] += s * F.invNorm[f];
        }
    }, rowGrain(F.nx));
}

// One V-cycle from x = 0: smooth, hand the residual down, add the coarse correction, smooth in the
// reverse colour order. The last level is solved by many sweeps forward and as many backward.
void PressureMultigrid::vcycle(int level) {
    Level& L = levels_[size_t(level)];
    std::fill(L.x.begin(), L.x.end(), 0.0);
    if (level + 1 == int(levels_.size())) {
        smooth(L, kCoarsestSweeps, true);
        smooth(L, kCoarsestSweeps, false);
        return;
    }
    Level& C = levels_[size_t(level) + 1];
    smoothBand(L, kBandSweeps, true);
    smooth(L, smoothingSweeps, true);
    computeResidual(L);
    restrictResidual(L, C);
    vcycle(level + 1);
    interpolateCorrection(C, L);
    smooth(L, smoothingSweeps, false);
    smoothBand(L, kBandSweeps, false);
}

void PressureMultigrid::apply(const std::vector<double>& r, std::vector<double>& z) {
    Level& L = levels_[0];
    parallelFor(int(L.cells()), [&](int c) { L.b[size_t(c)] = r[size_t(c)]; }, 4096);
    vcycle(0);
    parallelFor(int(L.cells()), [&](int c) { z[size_t(c)] = L.x[size_t(c)]; }, 4096);
}

} // namespace rf
