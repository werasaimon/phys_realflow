// Viscosity of GasSolver: the implicit step of the viscous term for every velocity component on
// its own faces of the MAC grid,
//     (I - a L) u_new = u_old - dt/rho grad p_old,     a = nu dt / dx^2,
// L the 7-point Laplacian, then u_new + dt/rho grad p_old goes on to the pressure projection.
// Three things make its answer second order and independent of the time step (they were found by
// the verification cases gas-poiseuille and gas-couette, verification/cases/GasCases.cpp):
//  1. No-slip on the wall itself. A velocity sample whose two cells are both solid lies half a
//     cell inside the wall; taken as it is, it moves the wall half a cell out and the profile is
//     first order. It is replaced by the ghost value u_ghost = 2 u_wall - u: the straight line
//     through the gas sample and the wall's own velocity on the solid face between them (Gibou,
//     Fedkiw, Cheng & Kang 2002, J. Comput. Phys. 176, a Dirichlet boundary at a known distance -
//     here always half a cell). The matrix stays symmetric: the ghost only adds a to the diagonal.
//  2. Solved to a tolerance: the conjugate gradient (Jacobi preconditioner) runs until the
//     relative residual is params.viscosityTolerance. A fixed count of sweeps (20 Jacobi before)
//     leaves an error that depends on a, so the steady flow depended on dt.
//  3. Incremental pressure correction (Goda 1979; van Kan 1986, SIAM J. Sci. Stat. Comput. 7;
//     Guermond, Minev & Shen 2006, Comput. Methods Appl. Mech. Engrg. 195): the old pressure
//     gradient is part of the right-hand side and is given back after the solve. The projection
//     starts from the old pressure (its warm start), so it then solves exactly for the pressure's
//     change. Without it the viscous step never sees the pressure force, and after the projection
//     the wall slips by dt |grad p| / rho - first order in dt, a quarter of the mean speed on the
//     8-cell Poiseuille channel.
// The unknowns are the faces between two gas cells. Every other face - on a side of the domain or
// next to a solid - keeps the value the boundary conditions gave it. A neighbour outside the
// domain, or across an outflow side, mirrors the face (free slip on the domain's walls, as the
// Taylor-Green case needs; zero gradient through an outflow).
#include "gas/GasSolver.h"

#include "core/Parallel.h"
#include "core/Probe.h"

#include <cmath>

namespace rf {

static double sumOfProducts(const std::vector<double>& a, const std::vector<double>& b) {
    return parallelSum<double>(int(a.size()), [&](int f0, int f1) {
        double s = 0;
        for (int f = f0; f < f1; ++f) s += a[size_t(f)] * b[size_t(f)];
        return s;
    }, 4096);
}

void GasSolver::diffuse(float dt) {
    const float a = params.kinematicViscosity * dt / (dx_ * dx_);
    if (a < 1e-3f) return; // negligible next to the numerical diffusion of advection
    applyVelocityBC();     // the kept faces hold their boundary values
    Field3* components[3] = {&u_, &v_, &w_};
    int most = 0;
    for (int c = 0; c < 3; ++c) {
        Field3& F = *components[c];
        markViscousUnknowns(c, F);
        buildViscousSystem(c, F, a, dt);
        most = std::max(most, solveViscousPcg(F, a));
    }
    lastViscousIters_ = most;
    Probe::set("gas/viscous iterations", most);
}

// The two cells a face of component c lies between: (i, j, k) minus one along c, and (i, j, k).
// false for a face on a side of the domain (one of its cells is outside).
bool GasSolver::faceCells(int c, int i, int j, int k, size_t& before, size_t& after) const {
    int ib = i, jb = j, kb = k;
    (c == 0 ? ib : c == 1 ? jb : kb) -= 1;
    if (ib < 0 || jb < 0 || kb < 0 || i >= nx_ || j >= ny_ || k >= nz_) return false;
    before = cidx(ib, jb, kb);
    after = cidx(i, j, k);
    return true;
}

// 1 on every face between two gas cells: the unknowns of the viscous equation.
void GasSolver::markViscousUnknowns(int c, const Field3& F) {
    viscKind_.assign(F.d.size(), 0);
    parallelFor(F.nz, [&](int k) {
        for (int j = 0; j < F.ny; ++j)
            for (int i = 0; i < F.nx; ++i) {
                size_t before, after;
                if (faceCells(c, i, j, k, before, after) && !solid_[before] && !solid_[after]) viscKind_[F.idx(i, j, k)] = 1;
            }
    }, 1);
}

// dt/rho (p_after - p_before) / dx on a face between two gas cells: the velocity change the old
// pressure makes there, with the face's Boris weight as the projection applies it.
double GasSolver::oldPressureStep(int c, int i, int j, int k, float dt) const {
    size_t before, after;
    if (!faceCells(c, i, j, k, before, after)) return 0;
    double weight = 1;
    const Field3& W = c == 0 ? weightU_ : c == 1 ? weightV_ : weightW_;
    if (magnetic.enabled && magnetic.speedLimit > 0 && !W.d.empty()) weight = W.at(i, j, k);
    return weight * double(dt) / (double(params.fluidDensity) * dx_) * (double(p_.d[after]) - double(p_.d[before]));
}

// The diagonal and the right-hand side of every unknown face; the start of the solve is the
// velocity before it, less the old pressure's step (kept in viscG_, given back after the solve).
void GasSolver::buildViscousSystem(int c, const Field3& F, float a, float dt) {
    const size_t n = F.d.size();
    viscDiag_.assign(n, 1.0);
    viscB_.assign(n, 0.0);
    viscG_.assign(n, 0.0);
    viscX_.assign(n, 0.0);
    parallelFor(F.nz, [&](int k) {
        for (int j = 0; j < F.ny; ++j)
            for (int i = 0; i < F.nx; ++i) {
                const size_t f = F.idx(i, j, k);
                if (!viscKind_[f]) continue;
                const double g = oldPressureStep(c, i, j, k, dt);
                double diag = 1 + 6.0 * a, b = double(F.d[f]) - g;
                const int at[3] = {i, j, k};
                for (int d = 0; d < 3; ++d) {
                    addViscousNeighbour(c, F, at, d, -1, a, diag, b);
                    addViscousNeighbour(c, F, at, d, +1, a, diag, b);
                }
                viscDiag_[f] = diag;
                viscB_[f] = b;
                viscG_[f] = g;
                viscX_[f] = double(F.d[f]) - g;
            }
    }, 1);
}

// One of the six neighbours of an unknown face, a step s = -1 / +1 along axis d:
//  - another unknown: the matrix couples them (nothing here);
//  - outside the domain, or the face on an outflow side along c: a mirror, u_nb = u;
//  - a face inside a wall (both its cells solid, beside us along d != c): the ghost 2 u_wall - u,
//    u_wall the solid's velocity on the solid face half a cell back towards us;
//  - any other kept face (a domain side, a face touching one solid cell): its value, at its place.
void GasSolver::addViscousNeighbour(int c, const Field3& F, const int at[3], int d, int s, float a, double& diag,
                                    double& b) const {
    int n[3] = {at[0], at[1], at[2]};
    n[d] += s;
    const int size[3] = {F.nx, F.ny, F.nz};
    const bool outflowSide = d == c && (n[d] == 0 || n[d] == size[d] - 1) && params.bc[2 * c + (s > 0)] == BoundaryType::Outflow;
    if (n[d] < 0 || n[d] >= size[d] || outflowSide) {
        diag -= a;
        return;
    }
    const size_t fn = F.idx(n[0], n[1], n[2]);
    if (viscKind_[fn]) return;
    size_t before, after;
    if (d == c || !faceCells(c, n[0], n[1], n[2], before, after) || !(solid_[before] && solid_[after])) {
        b += a * double(F.d[fn]);
        return;
    }
    const Vector3 face{float(n[0]), float(n[1]), float(n[2])};
    Vector3 wall = origin_ + (face + F.offset) * dx_;
    wall[d] -= 0.5f * float(s) * dx_;
    float uWall = 0;
    solidFaceVelocity(long(before), long(after), wall, c, uWall);
    diag += a;
    b += 2.0 * a * double(uWall);
}

// out = A x on the unknown faces: the diagonal times x minus a times the unknown neighbours.
void GasSolver::viscousMatVec(const Field3& F, float a, const std::vector<double>& x, std::vector<double>& out) const {
    const int sx = F.nx, sy = F.ny, sz = F.nz;
    const size_t slab = size_t(sx) * size_t(sy);
    parallelFor(sz, [&](int k) {
        for (int j = 0; j < sy; ++j)
            for (int i = 0; i < sx; ++i) {
                const size_t f = F.idx(i, j, k);
                if (!viscKind_[f]) {
                    out[f] = 0;
                    continue;
                }
                double nb = 0;
                if (i > 0 && viscKind_[f - 1]) nb += x[f - 1];
                if (i < sx - 1 && viscKind_[f + 1]) nb += x[f + 1];
                if (j > 0 && viscKind_[f - size_t(sx)]) nb += x[f - size_t(sx)];
                if (j < sy - 1 && viscKind_[f + size_t(sx)]) nb += x[f + size_t(sx)];
                if (k > 0 && viscKind_[f - slab]) nb += x[f - slab];
                if (k < sz - 1 && viscKind_[f + slab]) nb += x[f + slab];
                out[f] = viscDiag_[f] * x[f] - a * nb;
            }
    }, 1);
}

// The conjugate gradient with the Jacobi preconditioner (the matrix is symmetric and strictly
// diagonally dominant, so positive definite). It stops at the relative residual
// params.viscosityTolerance - relative to the right-hand side, or to the fastest velocity of the
// grid if this component is nearly still. Then the old pressure's step is given back. Returns the
// iterations used.
int GasSolver::solveViscousPcg(Field3& F, float a) {
    const int n = int(F.d.size());
    viscR_.assign(size_t(n), 0.0);
    viscZ_.assign(size_t(n), 0.0);
    viscP_.assign(size_t(n), 0.0);
    viscAp_.assign(size_t(n), 0.0);
    viscousMatVec(F, a, viscX_, viscAp_);
    double rz = parallelSum<double>(n, [&](int f0, int f1) {
        double sum = 0;
        for (size_t f = size_t(f0); f < size_t(f1); ++f) {
            viscR_[f] = viscKind_[f] ? viscB_[f] - viscAp_[f] : 0.0;
            viscZ_[f] = viscR_[f] / viscDiag_[f];
            viscP_[f] = viscZ_[f];
            sum += viscR_[f] * viscZ_[f];
        }
        return sum;
    }, 4096);
    const double scale = std::max(std::sqrt(sumOfProducts(viscB_, viscB_)), double(maxVel_) * std::sqrt(double(n)));
    const double tolerance = double(params.viscosityTolerance) * scale;
    double rnorm = std::sqrt(sumOfProducts(viscR_, viscR_));
    int it = 0;
    for (; it < params.maxViscosityIterations && rnorm > tolerance; ++it) {
        viscousMatVec(F, a, viscP_, viscAp_);
        const double pAp = sumOfProducts(viscP_, viscAp_);
        if (!(pAp > 0)) break;
        const double alpha = rz / pAp;
        struct Two {
            double rz = 0, rr = 0;
            Two& operator+=(const Two& o) { rz += o.rz; rr += o.rr; return *this; }
        };
        const Two t = parallelSum<Two>(n, [&](int f0, int f1) {
            Two sum;
            for (size_t f = size_t(f0); f < size_t(f1); ++f) {
                viscX_[f] += alpha * viscP_[f];
                viscR_[f] -= alpha * viscAp_[f];
                viscZ_[f] = viscR_[f] / viscDiag_[f];
                sum.rz += viscR_[f] * viscZ_[f];
                sum.rr += viscR_[f] * viscR_[f];
            }
            return sum;
        }, 4096);
        rnorm = std::sqrt(t.rr);
        const double beta = t.rz / rz;
        rz = t.rz;
        parallelFor(n, [&](int f) { viscP_[size_t(f)] = viscZ_[size_t(f)] + beta * viscP_[size_t(f)]; }, 4096);
    }
    parallelFor(n, [&](int f) {
        if (viscKind_[size_t(f)]) F.d[size_t(f)] = float(viscX_[size_t(f)] + viscG_[size_t(f)]);
    }, 4096);
    return it;
}

} // namespace rf
