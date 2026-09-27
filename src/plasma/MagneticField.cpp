// Resistive magnetohydrodynamics on the gas grid: the magnetic field on the faces (constrained
// transport keeps div B = 0 to rounding), Faraday's law with the flow and the resistivity, the
// Lorentz force on the gas, the Boris correction that limits the Alfven speed, Joule heating, and
// the field's energy and divergence for the tests. The equations are in MagneticField.h.
#include "plasma/MagneticField.h"

#include "core/Parallel.h"

#include <algorithm>
#include <cmath>

namespace rf {

// Layout (cell (i, j, k) spans [i, i+1] x [j, j+1] x [k, k+1] in cell units):
//   bx(i,j,k) on the x-face at (i, j+1/2, k+1/2)    size (nx+1, ny, nz)   - as the velocity u
//   by(i,j,k) on the y-face at (i+1/2, j, k+1/2)    size (nx, ny+1, nz)
//   bz(i,j,k) on the z-face at (i+1/2, j+1/2, k)    size (nx, ny, nz+1)
//   ex, jx(i,j,k) on the x-edge at (i+1/2, j, k)    size (nx, ny+1, nz+1)
//   ey, jy(i,j,k) on the y-edge at (i, j+1/2, k)    size (nx+1, ny, nz+1)
//   ez, jz(i,j,k) on the z-edge at (i, j, k+1/2)    size (nx+1, ny+1, nz)

void MagneticField::reset(int nx, int ny, int nz, float dx, const Vector3& origin) {
    nx_ = nx; ny_ = ny; nz_ = nz;
    dx_ = dx;
    origin_ = origin;
    bx.init(nx + 1, ny, nz, {0.0f, 0.5f, 0.5f}, applied.x);
    by.init(nx, ny + 1, nz, {0.5f, 0.0f, 0.5f}, applied.y);
    bz.init(nx, ny, nz + 1, {0.5f, 0.5f, 0.0f}, applied.z);
    for (Field3* e : {&ex_, &jx_}) e->init(nx, ny + 1, nz + 1, {0.5f, 0.0f, 0.0f});
    for (Field3* e : {&ey_, &jy_}) e->init(nx + 1, ny, nz + 1, {0.0f, 0.5f, 0.0f});
    for (Field3* e : {&ez_, &jz_}) e->init(nx + 1, ny + 1, nz, {0.0f, 0.0f, 0.5f});
    b0x_ = b0y_ = b0z_ = Field3();
    etaCell_.clear();
    etaX_ = etaY_ = etaZ_ = Field3();
    etaMax_ = 0;
    conductor_.clear();
    markConductorEdges();
    updateTotal();
}

void MagneticField::setResistivityMap(const std::function<float(const Vector3&)>& eta) {
    etaCell_.assign(size_t(nx_) * ny_ * nz_, 0.0f);
    etaMax_ = 0;
    for (int k = 0; k < nz_; ++k)
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                const float e = eta(origin_ + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * dx_);
                etaCell_[size_t(i) + size_t(nx_) * (size_t(j) + size_t(ny_) * size_t(k))] = e;
                etaMax_ = std::max(etaMax_, e);
            }
    // On an edge: the mean of the 4 cells around it (clamped at the domain walls).
    auto cell = [&](int i, int j, int k) {
        i = std::clamp(i, 0, nx_ - 1); j = std::clamp(j, 0, ny_ - 1); k = std::clamp(k, 0, nz_ - 1);
        return etaCell_[size_t(i) + size_t(nx_) * (size_t(j) + size_t(ny_) * size_t(k))];
    };
    etaX_ = ex_; etaY_ = ey_; etaZ_ = ez_;
    for (int k = 0; k <= nz_; ++k)
        for (int j = 0; j <= ny_; ++j)
            for (int i = 0; i <= nx_; ++i) {
                if (i < nx_) etaX_.at(i, j, k) = 0.25f * (cell(i, j - 1, k - 1) + cell(i, j, k - 1) + cell(i, j - 1, k) + cell(i, j, k));
                if (j < ny_) etaY_.at(i, j, k) = 0.25f * (cell(i - 1, j, k - 1) + cell(i, j, k - 1) + cell(i - 1, j, k) + cell(i, j, k));
                if (k < nz_) etaZ_.at(i, j, k) = 0.25f * (cell(i - 1, j - 1, k) + cell(i, j - 1, k) + cell(i - 1, j, k) + cell(i, j, k));
            }
}

void MagneticField::updateTotal() {
    tbx_ = bx;
    tby_ = by;
    tbz_ = bz;
    if (b0x_.d.empty()) return;
    for (size_t i = 0; i < tbx_.d.size(); ++i) tbx_.d[i] += b0x_.d[i];
    for (size_t i = 0; i < tby_.d.size(); ++i) tby_.d[i] += b0y_.d[i];
    for (size_t i = 0; i < tbz_.d.size(); ++i) tbz_.d[i] += b0z_.d[i];
}

void MagneticField::setConductors(const std::vector<uint8_t>& conductor) {
    conductor_ = conductor;
    markConductorEdges();
}

void MagneticField::markConductorEdges() {
    // An edge on the surface of a conductor (one of the four cells around it is conductor, or
    // outside the domain) carries no tangential E.
    fixedEx_.assign(ex_.d.size(), 0);
    fixedEy_.assign(ey_.d.size(), 0);
    fixedEz_.assign(ez_.d.size(), 0);
    for (int k = 0; k <= nz_; ++k)
        for (int j = 0; j <= ny_; ++j)
            for (int i = 0; i <= nx_; ++i) {
                if (i < nx_)
                    fixedEx_[ex_.idx(i, j, k)] = conductor(i, j - 1, k - 1) || conductor(i, j, k - 1) || conductor(i, j - 1, k) || conductor(i, j, k);
                if (j < ny_)
                    fixedEy_[ey_.idx(i, j, k)] = conductor(i - 1, j, k - 1) || conductor(i, j, k - 1) || conductor(i - 1, j, k) || conductor(i, j, k);
                if (k < nz_)
                    fixedEz_[ez_.idx(i, j, k)] = conductor(i - 1, j - 1, k) || conductor(i, j - 1, k) || conductor(i - 1, j, k) || conductor(i, j, k);
            }
}

void MagneticField::addFromPotential(const std::function<Vector3(const Vector3&)>& A) {
    addCurl(A, bx, by, bz);
    computeCurrent(); // the diagnostics (currentAt) see it before the first step
    updateTotal();
}

void MagneticField::setBackgroundFromPotential(const std::function<Vector3(const Vector3&)>& A) {
    b0x_.init(bx.nx, bx.ny, bx.nz, bx.offset);
    b0y_.init(by.nx, by.ny, by.nz, by.offset);
    b0z_.init(bz.nx, bz.ny, bz.nz, bz.offset);
    addCurl(A, b0x_, b0y_, b0z_);
    updateTotal();
}

void MagneticField::addCurl(const std::function<Vector3(const Vector3&)>& A, Field3& bx, Field3& by, Field3& bz) const {
    // A on the edges, B += curl A on the faces - the same discrete curl as Faraday's law, so the
    // result is divergence-free to rounding.
    Field3 ax = ex_, ay = ey_, az = ez_;
    auto fill = [&](Field3& f, int comp) {
        parallelFor(f.nz, [&](int k) { // A may be costly (a sum of current loops)
            for (int j = 0; j < f.ny; ++j)
                for (int i = 0; i < f.nx; ++i)
                    f.at(i, j, k) = A(origin_ + (Vector3(float(i), float(j), float(k)) + f.offset) * dx_)[comp];
        }, 1);
    };
    fill(ax, 0);
    fill(ay, 1);
    fill(az, 2);
    const float s = 1.0f / dx_;
    for (int k = 0; k < nz_; ++k)
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i <= nx_; ++i) bx.at(i, j, k) += s * ((az.at(i, j + 1, k) - az.at(i, j, k)) - (ay.at(i, j, k + 1) - ay.at(i, j, k)));
    for (int k = 0; k < nz_; ++k)
        for (int j = 0; j <= ny_; ++j)
            for (int i = 0; i < nx_; ++i) by.at(i, j, k) += s * ((ax.at(i, j, k + 1) - ax.at(i, j, k)) - (az.at(i + 1, j, k) - az.at(i, j, k)));
    for (int k = 0; k <= nz_; ++k)
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) bz.at(i, j, k) += s * ((ay.at(i + 1, j, k) - ay.at(i, j, k)) - (ax.at(i, j + 1, k) - ax.at(i, j, k)));
}

void MagneticField::computeCurrent() {
    // J = curl B / mu0 on the edges. Across a wall the tangential B is taken as continuing
    // (index clamped): no current sheet in the wall itself.
    const float s = 1.0f / (dx_ * kMu0);
    auto ci = [&](int i) { return std::clamp(i, 0, nx_ - 1); };
    auto cj = [&](int j) { return std::clamp(j, 0, ny_ - 1); };
    auto ck = [&](int k) { return std::clamp(k, 0, nz_ - 1); };
    parallelFor(nz_ + 1, [&](int k) {
        for (int j = 0; j <= ny_; ++j)
            for (int i = 0; i <= nx_; ++i) {
                if (i < nx_)
                    jx_.at(i, j, k) = s * ((bz.at(i, cj(j), k) - bz.at(i, cj(j - 1), k)) - (by.at(i, j, ck(k)) - by.at(i, j, ck(k - 1))));
                if (j < ny_)
                    jy_.at(i, j, k) = s * ((bx.at(i, j, ck(k)) - bx.at(i, j, ck(k - 1))) - (bz.at(ci(i), j, k) - bz.at(ci(i - 1), j, k)));
                if (k < nz_)
                    jz_.at(i, j, k) = s * ((by.at(ci(i), j, k) - by.at(ci(i - 1), j, k)) - (bx.at(i, cj(j), k) - bx.at(i, cj(j - 1), k)));
            }
    }, 1);
}

void MagneticField::computeElectricField(const Field3& u, const Field3& v, const Field3& w, float substep) {
    // Ohm's law on every free edge: E = -u x B + eta_eff curl B, with u and B averaged onto the
    // edge; eta_eff = eta + the dissipation of the centred advection (numericalDissipation):
    // max(f |u| dx, |u|^2 h). The advection of B by the flow is what needs it; the Alfven waves
    // are stable without - the Lorentz force and Faraday's law are stepped one after the other,
    // symplectically. Keeping v_A out of it keeps the field frozen in near a strong magnet, where
    // v_A is large.
    updateTotal();
    const Field3 &bx = tbx_, &by = tby_, &bz = tbz_; // the whole field B0 + B1 (the current: B1 only)
    const float eta0 = resistivity(), f = numericalDissipation * dx_;
    const bool mapped = !etaCell_.empty();
    auto etaEff = [&](float eta, float a, float b) { // (a, b): the flow components on the edge
        const float u2 = a * a + b * b;
        return eta + std::max(f * std::sqrt(u2), u2 * substep);
    };
    parallelFor(nz_ + 1, [&](int k) {
        for (int j = 0; j <= ny_; ++j)
            for (int i = 0; i <= nx_; ++i) {
                if (i < nx_) { // x-edge (i+1/2, j, k)
                    const size_t e = ex_.idx(i, j, k);
                    if (fixedEx_[e]) ex_.d[e] = 0;
                    else {
                        const float vE = 0.5f * (v.at(i, j, k - 1) + v.at(i, j, k)), ByE = 0.5f * (by.at(i, j, k - 1) + by.at(i, j, k));
                        const float wE = 0.5f * (w.at(i, j - 1, k) + w.at(i, j, k)), BzE = 0.5f * (bz.at(i, j - 1, k) + bz.at(i, j, k));
                        ex_.d[e] = -(vE * BzE - wE * ByE) + etaEff(mapped ? etaX_.d[e] : eta0, vE, wE) * kMu0 * jx_.d[e];
                    }
                }
                if (j < ny_) { // y-edge (i, j+1/2, k)
                    const size_t e = ey_.idx(i, j, k);
                    if (fixedEy_[e]) ey_.d[e] = 0;
                    else {
                        const float wE = 0.5f * (w.at(i - 1, j, k) + w.at(i, j, k)), BzE = 0.5f * (bz.at(i - 1, j, k) + bz.at(i, j, k));
                        const float uE = 0.5f * (u.at(i, j, k - 1) + u.at(i, j, k)), BxE = 0.5f * (bx.at(i, j, k - 1) + bx.at(i, j, k));
                        ey_.d[e] = -(wE * BxE - uE * BzE) + etaEff(mapped ? etaY_.d[e] : eta0, wE, uE) * kMu0 * jy_.d[e];
                    }
                }
                if (k < nz_) { // z-edge (i, j, k+1/2)
                    const size_t e = ez_.idx(i, j, k);
                    if (fixedEz_[e]) ez_.d[e] = 0;
                    else {
                        const float uE = 0.5f * (u.at(i, j - 1, k) + u.at(i, j, k)), BxE = 0.5f * (bx.at(i, j - 1, k) + bx.at(i, j, k));
                        const float vE = 0.5f * (v.at(i - 1, j, k) + v.at(i, j, k)), ByE = 0.5f * (by.at(i - 1, j, k) + by.at(i, j, k));
                        ez_.d[e] = -(uE * ByE - vE * BxE) + etaEff(mapped ? etaZ_.d[e] : eta0, uE, vE) * kMu0 * jz_.d[e];
                    }
                }
            }
    }, 1);
}

void MagneticField::applyFaraday(float dt) {
    // dB/dt = -curl E: the flux through a face changes by the circulation of E around it.
    const float s = dt / dx_;
    parallelFor(nz_ + 1, [&](int k) {
        for (int j = 0; j <= ny_; ++j)
            for (int i = 0; i <= nx_; ++i) {
                if (j < ny_ && k < nz_)
                    bx.at(i, j, k) -= s * ((ez_.at(i, j + 1, k) - ez_.at(i, j, k)) - (ey_.at(i, j, k + 1) - ey_.at(i, j, k)));
                if (i < nx_ && k < nz_)
                    by.at(i, j, k) -= s * ((ex_.at(i, j, k + 1) - ex_.at(i, j, k)) - (ez_.at(i + 1, j, k) - ez_.at(i, j, k)));
                if (i < nx_ && j < ny_)
                    bz.at(i, j, k) -= s * ((ey_.at(i + 1, j, k) - ey_.at(i, j, k)) - (ex_.at(i, j + 1, k) - ex_.at(i, j, k)));
            }
    }, 1);
}

void MagneticField::induce(const Field3& u, const Field3& v, const Field3& w, float density, float dt) {
    updateTotal(); // the faces may have been set from outside
    float umax = 0;
    for (const Field3* f : {&u, &v, &w})
        for (float x : f->d) umax = std::max(umax, std::fabs(x));
    const float B = maxField();
    const float signal = umax + alfvenSpeed(B * B, density);
    const float etaMax = maxResistivity() + numericalDissipation * umax * dx_;
    // Explicit limits: signals cross at most half a cell, diffusion stays below its bound.
    const float hMax = std::min(0.5f * dx_ / std::max(signal, 1e-12f), 0.9f * dx_ * dx_ / (6.0f * std::max(etaMax, 1e-20f)));
    const int steps = std::max(1, int(std::ceil(dt / hMax)));
    for (int s = 0; s < steps; ++s) {
        computeCurrent();
        computeElectricField(u, v, w, dt / float(steps));
        applyFaraday(dt / float(steps));
    }
    computeCurrent();
    updateTotal();
}

void MagneticField::applyLorentzForce(Field3& u, Field3& v, Field3& w, const std::vector<uint8_t>& solid, float density,
                                      float dt) {
    computeCurrent();
    updateTotal();
    const Field3 &bx = tbx_, &by = tby_, &bz = tbz_; // the whole field B0 + B1 (the current: B1 only)
    const float s0 = dt / std::max(density, 1e-12f);
    const float boris = speedLimit > 0 ? 1.0f / (kMu0 * std::max(density, 1e-12f) * speedLimit * speedLimit) : 0.0f;
    // Boris correction: the field's inertia (B^2 / mu0 c^2) is added to the gas's.
    auto scale = [&](float B2) { return s0 / (1.0f + boris * B2); };
    auto isSolid = [&](int i, int j, int k) {
        return !solid.empty() && solid[size_t(i) + size_t(nx_) * (size_t(j) + size_t(ny_) * size_t(k))] != 0;
    };
    parallelFor(nz_, [&](int k) {
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                if (i > 0 && !isSolid(i - 1, j, k) && !isSolid(i, j, k)) { // x-face (i, j+1/2, k+1/2)
                    const float Jy = 0.5f * (jy_.at(i, j, k) + jy_.at(i, j, k + 1)), Jz = 0.5f * (jz_.at(i, j, k) + jz_.at(i, j + 1, k));
                    const float Bz = 0.25f * (bz.at(i - 1, j, k) + bz.at(i, j, k) + bz.at(i - 1, j, k + 1) + bz.at(i, j, k + 1));
                    const float By = 0.25f * (by.at(i - 1, j, k) + by.at(i, j, k) + by.at(i - 1, j + 1, k) + by.at(i, j + 1, k));
                    const float Bx = bx.at(i, j, k);
                    u.at(i, j, k) += scale(Bx * Bx + By * By + Bz * Bz) * (Jy * Bz - Jz * By);
                }
                if (j > 0 && !isSolid(i, j - 1, k) && !isSolid(i, j, k)) { // y-face (i+1/2, j, k+1/2)
                    const float Jz = 0.5f * (jz_.at(i, j, k) + jz_.at(i + 1, j, k)), Jx = 0.5f * (jx_.at(i, j, k) + jx_.at(i, j, k + 1));
                    const float Bx = 0.25f * (bx.at(i, j - 1, k) + bx.at(i, j, k) + bx.at(i + 1, j - 1, k) + bx.at(i + 1, j, k));
                    const float Bz = 0.25f * (bz.at(i, j - 1, k) + bz.at(i, j, k) + bz.at(i, j - 1, k + 1) + bz.at(i, j, k + 1));
                    const float By = by.at(i, j, k);
                    v.at(i, j, k) += scale(Bx * Bx + By * By + Bz * Bz) * (Jz * Bx - Jx * Bz);
                }
                if (k > 0 && !isSolid(i, j, k - 1) && !isSolid(i, j, k)) { // z-face (i+1/2, j+1/2, k)
                    const float Jx = 0.5f * (jx_.at(i, j, k) + jx_.at(i, j + 1, k)), Jy = 0.5f * (jy_.at(i, j, k) + jy_.at(i + 1, j, k));
                    const float By = 0.25f * (by.at(i, j, k - 1) + by.at(i, j, k) + by.at(i, j + 1, k - 1) + by.at(i, j + 1, k));
                    const float Bx = 0.25f * (bx.at(i, j, k - 1) + bx.at(i, j, k) + bx.at(i + 1, j, k - 1) + bx.at(i + 1, j, k));
                    const float Bz = bz.at(i, j, k);
                    w.at(i, j, k) += scale(Bx * Bx + By * By + Bz * Bz) * (Jx * By - Jy * Bx);
                }
            }
    }, 1);
}

void MagneticField::jouleHeating(std::vector<float>& heat, float dt) {
    // J^2 / sigma = mu0 eta J^2 per cell.
    heat.assign(size_t(nx_) * ny_ * nz_, 0.0f);
    parallelFor(nz_, [&](int k) {
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                const size_t c = size_t(i) + size_t(nx_) * (size_t(j) + size_t(ny_) * size_t(k));
                const float eta = etaCell_.empty() ? resistivity() : etaCell_[c];
                heat[c] = dt * kMu0 * eta * length2(cellCurrent(i, j, k));
            }
    }, 1);
}

void MagneticField::borisWeights(Field3& wx, Field3& wy, Field3& wz, float density) const {
    wx.init(nx_ + 1, ny_, nz_, bx.offset, 1.0f);
    wy.init(nx_, ny_ + 1, nz_, by.offset, 1.0f);
    wz.init(nx_, ny_, nz_ + 1, bz.offset, 1.0f);
    if (speedLimit <= 0) return;
    const float k = 1.0f / (kMu0 * std::max(density, 1e-12f) * speedLimit * speedLimit); // (v_A / c)^2 per B^2
    auto ci = [&](int i) { return std::clamp(i, 0, nx_ - 1); };
    auto cj = [&](int j) { return std::clamp(j, 0, ny_ - 1); };
    auto ck = [&](int q) { return std::clamp(q, 0, nz_ - 1); };
    // |B|^2 on a face: its own component there, the other two averaged from the 4 faces around.
    parallelFor(nz_ + 1, [&](int kk) {
        for (int j = 0; j <= ny_; ++j)
            for (int i = 0; i <= nx_; ++i) {
                if (j < ny_ && kk < nz_) { // x-face (i, j+1/2, k+1/2)
                    const float Bx = tbx_.at(i, j, kk);
                    const float By = 0.25f * (tby_.at(ci(i - 1), j, kk) + tby_.at(ci(i), j, kk) + tby_.at(ci(i - 1), j + 1, kk) + tby_.at(ci(i), j + 1, kk));
                    const float Bz = 0.25f * (tbz_.at(ci(i - 1), j, kk) + tbz_.at(ci(i), j, kk) + tbz_.at(ci(i - 1), j, kk + 1) + tbz_.at(ci(i), j, kk + 1));
                    wx.at(i, j, kk) = 1.0f / (1.0f + k * (Bx * Bx + By * By + Bz * Bz));
                }
                if (i < nx_ && kk < nz_) { // y-face (i+1/2, j, k+1/2)
                    const float By = tby_.at(i, j, kk);
                    const float Bx = 0.25f * (tbx_.at(i, cj(j - 1), kk) + tbx_.at(i, cj(j), kk) + tbx_.at(i + 1, cj(j - 1), kk) + tbx_.at(i + 1, cj(j), kk));
                    const float Bz = 0.25f * (tbz_.at(i, cj(j - 1), kk) + tbz_.at(i, cj(j), kk) + tbz_.at(i, cj(j - 1), kk + 1) + tbz_.at(i, cj(j), kk + 1));
                    wy.at(i, j, kk) = 1.0f / (1.0f + k * (Bx * Bx + By * By + Bz * Bz));
                }
                if (i < nx_ && j < ny_) { // z-face (i+1/2, j+1/2, k)
                    const float Bz = tbz_.at(i, j, kk);
                    const float Bx = 0.25f * (tbx_.at(i, j, ck(kk - 1)) + tbx_.at(i, j, ck(kk)) + tbx_.at(i + 1, j, ck(kk - 1)) + tbx_.at(i + 1, j, ck(kk)));
                    const float By = 0.25f * (tby_.at(i, j, ck(kk - 1)) + tby_.at(i, j, ck(kk)) + tby_.at(i, j + 1, ck(kk - 1)) + tby_.at(i, j + 1, ck(kk)));
                    wz.at(i, j, kk) = 1.0f / (1.0f + k * (Bx * Bx + By * By + Bz * Bz));
                }
            }
    }, 1);
}

float MagneticField::maxTimeStep(float maxFlowSpeed, float density) const {
    const float B = maxField();
    return 0.5f * dx_ / std::max(maxFlowSpeed + alfvenSpeed(B * B, density), 1e-12f);
}

// The accessors below report the whole field B0 + B1 (tbx_..., refreshed after every change).
Vector3 MagneticField::fieldAt(const Vector3& world) const {
    const Vector3 g = (world - origin_) / dx_;
    return {tbx_.sample(g), tby_.sample(g), tbz_.sample(g)};
}

Vector3 MagneticField::currentAt(const Vector3& world) const {
    const Vector3 g = (world - origin_) / dx_;
    return {jx_.sample(g), jy_.sample(g), jz_.sample(g)};
}

Vector3 MagneticField::electricFieldAt(const Vector3& world) const {
    const Vector3 g = (world - origin_) / dx_;
    return {ex_.sample(g), ey_.sample(g), ez_.sample(g)};
}

Vector3 MagneticField::cellField(int i, int j, int k) const {
    return {0.5f * (tbx_.at(i, j, k) + tbx_.at(i + 1, j, k)), 0.5f * (tby_.at(i, j, k) + tby_.at(i, j + 1, k)),
            0.5f * (tbz_.at(i, j, k) + tbz_.at(i, j, k + 1))};
}

// The four edges of a cell along each axis, averaged.
static Vector3 cellEdgeAverage(const Field3& x, const Field3& y, const Field3& z, int i, int j, int k) {
    return {0.25f * (x.at(i, j, k) + x.at(i, j + 1, k) + x.at(i, j, k + 1) + x.at(i, j + 1, k + 1)),
            0.25f * (y.at(i, j, k) + y.at(i + 1, j, k) + y.at(i, j, k + 1) + y.at(i + 1, j, k + 1)),
            0.25f * (z.at(i, j, k) + z.at(i + 1, j, k) + z.at(i, j + 1, k) + z.at(i + 1, j + 1, k))};
}

Vector3 MagneticField::cellCurrent(int i, int j, int k) const { return cellEdgeAverage(jx_, jy_, jz_, i, j, k); }
Vector3 MagneticField::cellElectricField(int i, int j, int k) const { return cellEdgeAverage(ex_, ey_, ez_, i, j, k); }

float MagneticField::maxField() const {
    return parallelMax<float>(nz_, 0.0f, [&](int k0, int k1) {
        float m = 0;
        for (int k = k0; k < k1; ++k)
            for (int j = 0; j < ny_; ++j)
                for (int i = 0; i < nx_; ++i) m = std::max(m, length2(cellField(i, j, k)));
        return std::sqrt(m);
    }, 1);
}

double MagneticField::energy() const {
    double e = 0;
    for (int k = 0; k < nz_; ++k)
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) e += length2(cellField(i, j, k));
    return e * double(dx_) * dx_ * dx_ / (2.0 * kMu0);
}

float MagneticField::maxDivergence() const {
    float m = 0;
    for (int k = 0; k < nz_; ++k)
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i)
                m = std::max(m, std::fabs(tbx_.at(i + 1, j, k) - tbx_.at(i, j, k) + tby_.at(i, j + 1, k) - tby_.at(i, j, k) +
                                          tbz_.at(i, j, k + 1) - tbz_.at(i, j, k)));
    return m / std::max(maxField(), 1e-30f);
}

} // namespace rf
