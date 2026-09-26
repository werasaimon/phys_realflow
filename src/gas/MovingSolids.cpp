// Bodies in the gas of GasSolver: rigid bodies and the liquid voxelised as moving solids, the
// no-through condition on their faces, the pressure force and torque on them, the wall function
// of the skin friction, and outside impulses put into the gas.
#include "gas/GasSolver.h"

#include "core/Parallel.h"

namespace rf {

// ---------------------------------------------------------------------------
// Moving solids (rigid bodies)
// ---------------------------------------------------------------------------
void GasSolver::voxelizeMovingSolids() {
    if (moving_.empty() && liquidPos_.empty() && !hadMoving_) return;
    const size_t nb = moving_.size();
    // 1) The cells of every body: those whose centre is inside it. A resting (sleeping) body keeps
    //    the cells of its last voxelisation - no inside tests at all. Then the liquid: the cells
    //    holding at least one particle, with the particles' mean velocity (list nb).
    std::vector<std::vector<size_t>> cells(nb + 1);
    // Cell lists are matched by index: reuse them only while the set of solids is the same size.
    const bool sameSolids = solidCells_.size() == nb + 1; // (the last list is the liquid)
    for (size_t b = 0; b < nb; ++b) {
        const MovingSolid& s = moving_[b];
        if (s.resting && sameSolids) {
            cells[b] = solidCells_[b];
            continue;
        }
        Vector3 lo = (s.bounds.lo - origin_) / dx_ - Vector3(0.5f), hi = (s.bounds.hi - origin_) / dx_ - Vector3(0.5f);
        int i0 = std::max(0, int(std::ceil(lo.x))), i1 = std::min(nx_ - 1, int(std::floor(hi.x)));
        int j0 = std::max(0, int(std::ceil(lo.y))), j1 = std::min(ny_ - 1, int(std::floor(hi.y)));
        int k0 = std::max(0, int(std::ceil(lo.z))), k1 = std::min(nz_ - 1, int(std::floor(hi.z)));
        for (int k = k0; k <= k1; ++k)
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i) {
                    size_t c = cidx(i, j, k);
                    if (solid_[c] & kStatic) continue;
                    if (s.inside(origin_ + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * dx_)) cells[b].push_back(c);
                }
    }
    if (!liquidPos_.empty()) {
        const size_t n = solid_.size();
        std::vector<int> count(n, 0);
        liquidVel_.assign(n, Vector3(0.0f));
        for (size_t p = 0; p < liquidPos_.size(); ++p) {
            Vector3 g = (liquidPos_[p] - origin_) / dx_;
            int i = int(std::floor(g.x)), j = int(std::floor(g.y)), k = int(std::floor(g.z));
            if (i < 0 || j < 0 || k < 0 || i >= nx_ || j >= ny_ || k >= nz_) continue;
            size_t c = cidx(i, j, k);
            if (solid_[c] & kStatic) continue;
            if (count[c]++ == 0) cells[nb].push_back(c);
            liquidVel_[c] += liquidVelIn_[p];
        }
        for (size_t c : cells[nb]) liquidVel_[c] /= float(count[c]);
    }

    // 2) Unmark the old cells (remembering them), mark the new ones.
    const uint32_t wasMark = ++stampId_;
    for (const auto& list : solidCells_)
        for (size_t c : list) {
            stamp_[c] = wasMark;
            solid_[c] &= uint8_t(~kMoving);
            owner_[c] = -1;
        }
    bool changed = false;
    movingCells_ = 0;
    std::vector<size_t> entered; // cells a body moved into
    for (size_t b = 0; b < cells.size(); ++b)
        for (size_t c : cells[b]) {
            if (solid_[c]) continue; // already taken by another body
            solid_[c] |= kMoving;
            owner_[c] = int(b);
            ++movingCells_;
            if (stamp_[c] != wasMark) {
                changed = true;
                entered.push_back(c);
            }
        }
    for (const auto& list : solidCells_)
        for (size_t c : list)
            if (!(solid_[c] & kMoving)) changed = true; // a cell was left

    // 3) The smoke and heat of the entered cells are pushed out to their gas neighbours instead of
    //    vanishing inside the body (conserved as long as a gas neighbour exists).
    for (size_t c : entered) {
        int i = int(c % nx_), j = int((c / nx_) % ny_), k = int(c / (size_t(nx_) * ny_));
        const int nbr[6][3] = {{i - 1, j, k}, {i + 1, j, k}, {i, j - 1, k}, {i, j + 1, k}, {i, j, k - 1}, {i, j, k + 1}};
        size_t fluid[6];
        int nf = 0;
        for (auto& q : nbr)
            if (q[0] >= 0 && q[1] >= 0 && q[2] >= 0 && q[0] < nx_ && q[1] < ny_ && q[2] < nz_ && !solid(q[0], q[1], q[2]))
                fluid[nf++] = cidx(q[0], q[1], q[2]);
        for (int f = 0; f < nf; ++f) {
            smoke_.d[fluid[f]] += smoke_.d[c] / float(nf);
            temp_.d[fluid[f]] += temp_.d[c] / float(nf);
            fuel_.d[fluid[f]] += fuel_.d[c] / float(nf);
            products_.d[fluid[f]] = std::min(1.0f, products_.d[fluid[f]] + products_.d[c] / float(nf));
        }
        smoke_.d[c] = 0;
        temp_.d[c] = 0;
        fuel_.d[c] = 0;
        products_.d[c] = 0;
    }
    solidCells_.swap(cells);
    liquidCells_ = int(solidCells_.back().size());
    if (changed) computeDiag();
    hadMoving_ = !moving_.empty() || !liquidPos_.empty();
}

bool GasSolver::solidFaceVelocity(long ca, long cb, const Vector3& fw, int comp, float& out) const {
    float sum = 0;
    int cnt = 0;
    for (long c : {ca, cb}) {
        if (c < 0 || !solid_[c]) continue;
        if (solid_[c] & kMoving)
            sum += owner_[c] < int(moving_.size()) ? moving_[owner_[c]].pointVelocity(fw)[comp] : liquidVel_[c][comp];
        ++cnt;
    }
    if (cnt == 0) return false;
    out = sum / float(cnt);
    return true;
}

void GasSolver::applyPendingImpulses() {
    // The momentum goes to the gas cells around the point - the cells the drag read the gas
    // velocity from (fluidStencil, weights summing to 1) - and in each cell to its faces that the
    // boundary conditions leave free: a face next to a solid, the liquid or the domain boundary is
    // overwritten by applyVelocityBC and would lose it. Every face stands for one cell of gas, so
    // the momentum put in is exactly the impulse.
    const float cellMass = params.fluidDensity * dx_ * dx_ * dx_;
    Field3* comps[3] = {&u_, &v_, &w_};
    const int size[3] = {nx_, ny_, nz_};
    for (const auto& [x, J] : pendingImpulses_) {
        size_t cells[8];
        float weights[8];
        const int count = fluidStencil(x, cells, weights);
        for (int m = 0; m < count; ++m) {
            const size_t c = cells[m];
            const int cell[3] = {int(c % nx_), int((c / nx_) % ny_), int(c / (size_t(nx_) * ny_))};
            const Vector3 dv = J * (weights[m] / cellMass);
            for (int a = 0; a < 3; ++a) {
                bool faceFree[2]; // the cell's lower and upper face along axis a: gas on the other side?
                for (int side = 0; side < 2; ++side) {
                    int nb[3] = {cell[0], cell[1], cell[2]};
                    nb[a] += side ? 1 : -1;
                    faceFree[side] = nb[a] >= 0 && nb[a] < size[a] && !solid(nb[0], nb[1], nb[2]);
                }
                const int nfree = int(faceFree[0]) + int(faceFree[1]);
                for (int side = 0; side < 2; ++side) {
                    if (!faceFree[side]) continue;
                    int f[3] = {cell[0], cell[1], cell[2]};
                    f[a] += side; // index of the face along axis a
                    comps[a]->at(f[0], f[1], f[2]) += dv[a] / float(nfree);
                }
            }
        }
    }
    pendingImpulses_.clear();
}

void GasSolver::computeMovingForces() {
    movForce_.assign(moving_.size(), Vector3(0.0f));
    movTorque_.assign(moving_.size(), Vector3(0.0f));
    if (moving_.empty()) return;
    if (movFricForce_.size() == moving_.size()) {
        movForce_ = movFricForce_;
        movTorque_ = movFricTorque_;
    }
    // Force of the gas on the body: -p n dA summed over the body's faces towards gas cells.
    const float A = dx_ * dx_;
    const Vector3 dir[6] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
    for (size_t o = 0; o < solidCells_.size() && o < moving_.size(); ++o)
        for (size_t c : solidCells_[o]) {
            if (owner_[c] != int(o)) continue;
            int i = int(c % nx_), j = int((c / nx_) % ny_), k = int(c / (size_t(nx_) * ny_));
            const int nb[6][3] = {{i - 1, j, k}, {i + 1, j, k}, {i, j - 1, k}, {i, j + 1, k}, {i, j, k - 1}, {i, j, k + 1}};
            for (int f = 0; f < 6; ++f) {
                int a = nb[f][0], b = nb[f][1], cc = nb[f][2];
                if (a < 0 || b < 0 || cc < 0 || a >= nx_ || b >= ny_ || cc >= nz_ || solid(a, b, cc)) continue;
                Vector3 F = dir[f] * (-p_.at(a, b, cc) * A);
                Vector3 x = origin_ + (Vector3(i + 0.5f, j + 0.5f, k + 0.5f) + dir[f] * 0.5f) * dx_;
                movForce_[o] += F;
                movTorque_[o] += cross(x - moving_[o].position, F);
            }
        }
}

void GasSolver::applyWallFriction(float dt) {
    staticFriction_ = Vector3(0.0f);
    movFricForce_.assign(moving_.size(), Vector3(0.0f));
    movFricTorque_.assign(moving_.size(), Vector3(0.0f));
    if (!params.wallFriction || params.kinematicViscosity <= 0 || (solidCount_ == 0 && movingCells_ == 0)) return;
    const float nu = params.kinematicViscosity;
    const float cellMass = params.fluidDensity * dx_ * dx_ * dx_;
    Field3* comps[3] = {&u_, &v_, &w_};
    auto fluidCell = [&](int i, int j, int k) {
        return i >= 0 && j >= 0 && k >= 0 && i < nx_ && j < ny_ && k < nz_ && !solid(i, j, k);
    };

    // Gas cells next to a wall: the static obstacle's (precomputed) and the bodies' neighbours.
    const uint32_t mark = ++stampId_;
    std::vector<size_t> wallCells;
    auto addCell = [&](size_t c) {
        if (stamp_[c] == mark || solid_[c]) return;
        stamp_[c] = mark;
        wallCells.push_back(c);
    };
    for (size_t c : staticWallCells_) addCell(c);
    for (size_t o = 0; o < solidCells_.size() && o < moving_.size(); ++o)
        for (size_t c : solidCells_[o]) {
            int i = int(c % nx_), j = int((c / nx_) % ny_), k = int(c / (size_t(nx_) * ny_));
            if (i > 0) addCell(c - 1);
            if (i < nx_ - 1) addCell(c + 1);
            if (j > 0) addCell(c - nx_);
            if (j < ny_ - 1) addCell(c + nx_);
            if (k > 0) addCell(c - size_t(nx_) * ny_);
            if (k < nz_ - 1) addCell(c + size_t(nx_) * ny_);
        }

    // Velocity changes are gathered first (faces are shared by neighbouring cells), then applied.
    struct Change { int comp; size_t face; float dv; };
    std::vector<Change> changes;
    for (size_t c : wallCells) {
        const int i = int(c % nx_), j = int((c / nx_) % ny_), k = int(c / (size_t(nx_) * ny_));
        const int cell[3] = {i, j, k};
        for (int f = 0; f < 6; ++f) {
            const int axis = f / 2, sgn = (f & 1) ? 1 : -1;
            int nb[3] = {i, j, k};
            nb[axis] += sgn;
            if (nb[axis] < 0 || nb[0] >= nx_ || nb[1] >= ny_ || nb[2] >= nz_ || !solid(nb[0], nb[1], nb[2])) continue;
            const size_t sc = cidx(nb[0], nb[1], nb[2]);
            const int owner = owner_[sc];
            if (owner >= int(moving_.size())) continue; // liquid: coupled through the particles' drag
            Vector3 d(0.0f);
            d[axis] = float(sgn); // from the gas cell towards the wall
            const Vector3 xf = origin_ + (Vector3(i + 0.5f, j + 0.5f, k + 0.5f) + d * 0.5f) * dx_;
            const Vector3 vb = owner >= 0 ? moving_[owner].pointVelocity(xf) : Vector3(0.0f);
            Vector3 rel = cellVelocity(i, j, k) - vb;
            rel[axis] = 0; // tangential part only
            const float mag = length(rel);
            if (mag < 1e-6f) continue;
            const float L = owner >= 0 ? moving_[owner].length : staticLength_;
            const float cf = skinFrictionCoefficient(mag * L / nu);
            // tau dA dt / (rho dx^3) = k u_t, taken implicitly (never reverses the flow).
            const float kk = 0.5f * cf * mag * dt / dx_;
            const Vector3 delta = rel * (-kk / (1.0f + kk));
            // Momentum lives on the faces: half of the change on each of the cell's two faces of a
            // tangential axis, unless that face is a boundary face (fixed by the boundary condition).
            Vector3 applied(0.0f);
            for (int a = 0; a < 3; ++a) {
                if (a == axis || delta[a] == 0) continue;
                for (int side = 0; side < 2; ++side) {
                    int o[3] = {cell[0], cell[1], cell[2]};
                    o[a] += side ? 1 : -1;
                    if (!fluidCell(o[0], o[1], o[2])) continue;
                    int fi[3] = {cell[0], cell[1], cell[2]};
                    fi[a] += side; // face index along axis a
                    changes.push_back({a, comps[a]->idx(fi[0], fi[1], fi[2]), 0.5f * delta[a]});
                    applied[a] += 0.5f * delta[a];
                }
            }
            const Vector3 F = applied * (-cellMass / dt); // reaction on the wall
            if (owner >= 0) {
                movFricForce_[owner] += F;
                movFricTorque_[owner] += cross(xf - moving_[owner].position, F);
            } else if (solid_[sc] & kStatic) {
                staticFriction_ += F;
            }
        }
    }
    for (const Change& ch : changes) comps[ch.comp]->d[ch.face] += ch.dv;
}

void GasSolver::computeForces() {
    if (solidCount_ == 0) { // static obstacle (moving solids are reported per body)
        force_ = Vector3(0.0f);
        cd_ = cl_ = cs_ = 0;
        return;
    }
    const float A = dx_ * dx_;
    struct F3 {
        double x = 0, y = 0, z = 0;
        F3& operator+=(const F3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    };
    F3 F = parallelSum<F3>(nz_, [&](int k0, int k1) {
        F3 acc;
        for (int k = k0; k < k1; ++k)
            for (int j = 0; j < ny_; ++j)
                for (int i = 0; i < nx_; ++i) {
                    if (!staticSolid(i, j, k)) continue;
                    const int nb[6][3] = {{i - 1, j, k}, {i + 1, j, k}, {i, j - 1, k}, {i, j + 1, k}, {i, j, k - 1}, {i, j, k + 1}};
                    const float dir[6][3] = {{-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
                    for (int f = 0; f < 6; ++f) {
                        int a = nb[f][0], b = nb[f][1], c = nb[f][2];
                        if (a < 0 || b < 0 || c < 0 || a >= nx_ || b >= ny_ || c >= nz_ || solid(a, b, c)) continue;
                        float p = p_.at(a, b, c);
                        // force on the body = -p * n_out * dA
                        acc.x -= p * dir[f][0] * A;
                        acc.y -= p * dir[f][1] * A;
                        acc.z -= p * dir[f][2] * A;
                    }
                }
        return acc;
    }, 1);
    force_ = Vector3(float(F.x), float(F.y), float(F.z)) + staticFriction_;
    refArea_ = params.usePlanformArea ? planformArea_ : frontalArea_;
    float denom = dynamicPressure() * std::max(refArea_, 1e-9f);
    cd_ = force_.x / denom;
    cl_ = force_.y / denom;
    cs_ = force_.z / denom;
}

} // namespace rf
