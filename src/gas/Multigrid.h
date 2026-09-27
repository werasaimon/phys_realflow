#pragma once
// A geometric multigrid V-cycle for the gas's pressure equation, used as the preconditioner of the
// conjugate gradient (MGPCG: McAdams, Sifakis, Teran 2010, "A parallel multigrid Poisson solver for
// fluids simulation on large grids"). The Jacobi preconditioner fixes the scale of each cell but not
// the long waves of the pressure: its iteration count grows with the grid (~ n per side). A V-cycle
// removes every wavelength at once - the short ones on the fine grid by smoothing, the long ones on
// coarser grids where they look short - so the count stays nearly the same on 32^3 and 128^3.
//
// The equation of one level, cell by cell (the solver's unit form, see PressureSolver.cpp):
//     diag_c x_c - sum over faces f of c:  w_f x_(neighbour across f)  =  b_c,
//     diag_c = sum of the weights of c's six faces.
// A face weight is 0 towards a solid or a wall (Neumann), w towards gas, and w on an open side of
// the domain (Dirichlet, p = 0 outside: in the diagonal only).
//
// Coarser levels have cells twice as large along every axis that still has at least 4 cells. The
// coarse face weight is the mean of the fine face weights it covers, divided by the square of the
// coarsening factor (the Laplacian of a grid of spacing 2h); a coarse cell with no open face is
// left out. Transfers: trilinear interpolation from coarse to fine (weights 3/4 and 1/4 per axis,
// renormalised over the coarse cells that exist), and its transpose, scaled, from fine to coarse.
// Smoothing: red-black Gauss-Seidel, red then black before the coarse correction, black then red
// after it, so the whole V-cycle is a symmetric operator - a valid preconditioner for CG. Next to a
// solid the coarse levels only approximate the fine one (a coarse cell half in a body is still a
// gas cell), and the residual left there would stay as spikes of divergence at the body's surface;
// a narrow band of cells around solids gets extra sweeps first, as McAdams et al. (2010, sec. 4) do.
#include <cstddef>
#include <cstdint>
#include <vector>

namespace rf {

// How the conjugate gradient of the pressure is preconditioned (GasParams::pressurePreconditioner).
enum class PressurePreconditioner { Jacobi, Multigrid };

class PressureMultigrid {
public:
    // The finest level's size; the hierarchy is allocated here, only when the size changes.
    void resize(int nx, int ny, int nz);
    // The finest level's face weights, filled by the gas solver before build():
    // faceX(i, j, k), i = 0 .. nx: the face between cells (i-1, j, k) and (i, j, k); likewise Y, Z.
    float& faceX(int i, int j, int k) { return levels_[0].wx[levels_[0].faceX(i, j, k)]; }
    float& faceY(int i, int j, int k) { return levels_[0].wy[levels_[0].faceY(i, j, k)]; }
    float& faceZ(int i, int j, int k) { return levels_[0].wz[levels_[0].faceZ(i, j, k)]; }
    // Every coarser level from the finest one's weights (after the weights are filled).
    void build();
    // z = one V-cycle applied to the residual r (both with one value per fine cell).
    void apply(const std::vector<double>& r, std::vector<double>& z);
    int levelCount() const { return int(levels_.size()); }
    // Sweeps of the smoother before and after the coarse correction (each a red and a black pass).
    int smoothingSweeps = 2;

private:
    struct Level {
        int nx = 0, ny = 0, nz = 0;
        int fx = 1, fy = 1, fz = 1;       // coarsening factor from the finer level (1 or 2)
        std::vector<float> wx, wy, wz;     // face weights
        std::vector<double> diag;          // sum of the six face weights (0: the cell is left out)
        std::vector<double> invNorm;       // 1 / sum of the interpolation weights from the coarser level
        std::vector<double> x, b, r;       // unknown, right-hand side, residual
        std::vector<uint8_t> irregular;    // 1: a cell with a closed face inside the grid (next to a solid)
        std::vector<uint32_t> band[2];     // cells at most one cell from an irregular one, red and black
        size_t cell(int i, int j, int k) const { return size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * size_t(k)); }
        size_t faceX(int i, int j, int k) const { return size_t(i) + size_t(nx + 1) * (size_t(j) + size_t(ny) * size_t(k)); }
        size_t faceY(int i, int j, int k) const { return size_t(i) + size_t(nx) * (size_t(j) + size_t(ny + 1) * size_t(k)); }
        size_t faceZ(int i, int j, int k) const { return cell(i, j, k); } // k = 0 .. nz
        size_t cells() const { return size_t(nx) * ny * nz; }
    };
    void allocate(Level& L);
    static void computeDiagonal(Level& L);
    static void coarsenFaces(const Level& fine, Level& coarse);
    static void computeInterpolationNorms(Level& fine, const Level& coarse);
    static void findBoundaryBand(Level& L);
    static void relaxCell(Level& L, size_t c, int i, int j, int k);
    static void smooth(Level& L, int sweeps, bool redFirst);
    static void smoothBand(Level& L, int sweeps, bool redFirst);
    static void computeResidual(Level& L);
    static void restrictResidual(const Level& fine, Level& coarse);
    static void interpolateCorrection(const Level& coarse, Level& fine);
    void vcycle(int level);

    std::vector<Level> levels_;
};

} // namespace rf
