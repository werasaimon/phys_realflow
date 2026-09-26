#include "grid/NSGridSolver.h"

#include "core/Parallel.h"

#include <chrono>

namespace rf {

// ---------------------------------------------------------------------------
// Field sampling
// ---------------------------------------------------------------------------
static inline Vector3 sampleVel(const Field3& u, const Field3& v, const Field3& w, const Vector3& gp) {
    return {u.sample(gp), v.sample(gp), w.sample(gp)};
}

// ---------------------------------------------------------------------------
// Setup
// ---------------------------------------------------------------------------
void NSGridSolver::reset(const Vector3& origin, const MeshBVH* obstacle) {
    origin_ = origin;
    nx_ = std::max(8, params.resolutionX);
    dx_ = params.domainSize.x / float(nx_);
    ny_ = std::max(4, int(std::lround(params.domainSize.y / dx_)));
    nz_ = std::max(4, int(std::lround(params.domainSize.z / dx_)));

    float u0 = (params.bc[0] == BoundaryType::Inflow) ? params.inflowSpeed : 0.0f;
    u_.init(nx_ + 1, ny_, nz_, {0, 0.5f, 0.5f}, u0);
    v_.init(nx_, ny_ + 1, nz_, {0.5f, 0, 0.5f});
    w_.init(nx_, ny_, nz_ + 1, {0.5f, 0.5f, 0});
    u0_ = u_; v0_ = v_; w0_ = w_;
    p_.init(nx_, ny_, nz_, Vector3(0.5f));
    smoke_.init(nx_, ny_, nz_, Vector3(0.5f));
    temp_.init(nx_, ny_, nz_, Vector3(0.5f));
    fuel_.init(nx_, ny_, nz_, Vector3(0.5f));
    products_.init(nx_, ny_, nz_, Vector3(0.5f));

    const size_t n = size_t(nx_) * ny_ * nz_;
    expansion_.assign(n, 0.0f);
    heatReleaseRate_ = 0;
    pendingEmissions_.clear();
    radiators_.clear();
    solid_.assign(n, 0);
    if (obstacle && !obstacle->empty()) {
        AABB mb = obstacle->bounds();
        parallelFor(int(nz_), [&](int k_) {
            int k = k_;
            for (int j = 0; j < ny_; ++j)
                for (int i = 0; i < nx_; ++i) {
                    Vector3 c = origin_ + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * dx_;
                    if (mb.contains(c) && obstacle->isInside(c)) solid_[cidx(i, j, k)] = kStatic;
                }
        }, 1);
    }
    solidCount_ = 0;
    for (uint8_t s : solid_) solidCount_ += s != 0;
    owner_.assign(n, -1);
    movForce_.clear();
    movTorque_.clear();
    movingCells_ = 0;
    hadMoving_ = false;

    // Reference areas: projections of the voxelised body.
    int front = 0, plan = 0;
    for (int k = 0; k < nz_; ++k)
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i)
                if (solid(i, j, k)) { ++front; break; }
    for (int k = 0; k < nz_; ++k)
        for (int i = 0; i < nx_; ++i)
            for (int j = 0; j < ny_; ++j)
                if (solid(i, j, k)) { ++plan; break; }
    frontalArea_ = front * dx_ * dx_;
    planformArea_ = plan * dx_ * dx_;
    {
        AABB sb;
        for (int k = 0; k < nz_; ++k)
            for (int j = 0; j < ny_; ++j)
                for (int i = 0; i < nx_; ++i)
                    if (solid(i, j, k)) sb.expand(Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * dx_);
        staticLength_ = sb.valid() ? maxComp(sb.extent()) + dx_ : 1.0f;
    }
    staticFriction_ = Vector3(0.0f);
    movFricForce_.clear();
    movFricTorque_.clear();
    solidCells_.clear();
    stamp_.assign(n, 0);
    stampId_ = 0;
    // Gas cells touching the static obstacle (the only places where it can exert wall friction).
    staticWallCells_.clear();
    for (int k = 0; k < nz_; ++k)
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                if (solid(i, j, k)) continue;
                bool wall = (i > 0 && solid(i - 1, j, k)) || (i < nx_ - 1 && solid(i + 1, j, k)) || (j > 0 && solid(i, j - 1, k)) ||
                            (j < ny_ - 1 && solid(i, j + 1, k)) || (k > 0 && solid(i, j, k - 1)) || (k < nz_ - 1 && solid(i, j, k + 1));
                if (wall) staticWallCells_.push_back(cidx(i, j, k));
            }

    computeDiag();

    magnetic.reset(nx_, ny_, nz_, dx_, origin_);
    if (solidCount_ > 0) magnetic.setConductors(solid_); // the obstacle (a vessel wall) conducts

    q_.assign(n, 0); r_.assign(n, 0); z_.assign(n, 0); s_.assign(n, 0); As_.assign(n, 0); b_.assign(n, 0);
    force_ = Vector3(0.0f);
    cd_ = cl_ = cs_ = cdAvg_ = clAvg_ = 0;
    refArea_ = params.usePlanformArea ? planformArea_ : frontalArea_;
    lastIters_ = 0;
    lastResidual_ = 0;
    time_ = 0;
    applyVelocityBC();
}

// Pressure matrix diagonal: fluid neighbours + Dirichlet (outflow) sides.
void NSGridSolver::computeDiag() {
    diag_.assign(size_t(nx_) * ny_ * nz_, 0);
    parallelFor(nz_, [&](int k) {
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                if (solid(i, j, k)) continue;
                int d = 0;
                const int ni[6][3] = {{i - 1, j, k}, {i + 1, j, k}, {i, j - 1, k}, {i, j + 1, k}, {i, j, k - 1}, {i, j, k + 1}};
                for (int f = 0; f < 6; ++f) {
                    int a = ni[f][0], b = ni[f][1], c = ni[f][2];
                    if (a < 0 || b < 0 || c < 0 || a >= nx_ || b >= ny_ || c >= nz_) {
                        if (params.bc[f] == BoundaryType::Outflow) ++d;
                    } else if (!solid(a, b, c)) {
                        ++d;
                    }
                }
                diag_[cidx(i, j, k)] = uint8_t(d);
            }
    }, 1);

    // Fluid regions (flood fill) and whether each touches an outflow boundary.
    const size_t n = diag_.size();
    region_.assign(n, -1);
    regionOpen_.clear();
    std::vector<size_t> stack;
    for (size_t seed = 0; seed < n; ++seed) {
        if (solid_[seed] || diag_[seed] == 0 || region_[seed] >= 0) continue;
        const int r = int(regionOpen_.size());
        char open = 0;
        region_[seed] = r;
        stack.assign(1, seed);
        while (!stack.empty()) {
            size_t c = stack.back();
            stack.pop_back();
            int i = int(c % nx_), j = int((c / nx_) % ny_), k = int(c / (size_t(nx_) * ny_));
            const int ni[6][3] = {{i - 1, j, k}, {i + 1, j, k}, {i, j - 1, k}, {i, j + 1, k}, {i, j, k - 1}, {i, j, k + 1}};
            for (int f = 0; f < 6; ++f) {
                int a = ni[f][0], b = ni[f][1], cc = ni[f][2];
                if (a < 0 || b < 0 || cc < 0 || a >= nx_ || b >= ny_ || cc >= nz_) {
                    if (params.bc[f] == BoundaryType::Outflow) open = 1;
                    continue;
                }
                size_t q = cidx(a, b, cc);
                if (solid_[q] || diag_[q] == 0 || region_[q] >= 0) continue;
                region_[q] = r;
                stack.push_back(q);
            }
        }
        regionOpen_.push_back(open);
    }
}

// ---------------------------------------------------------------------------
// Moving solids (rigid bodies)
// ---------------------------------------------------------------------------
void NSGridSolver::voxelizeMovingSolids() {
    if (moving_.empty() && liquidPos_.empty() && !hadMoving_) return;
    const size_t nb = moving_.size();
    // 1) The cells of every body: those whose centre is inside it. A resting (sleeping) body keeps
    //    the cells of its last voxelisation - no inside tests at all. Then the liquid: the cells
    //    holding at least one particle, with the particles' mean velocity (list nb).
    std::vector<std::vector<size_t>> cells(nb + 1);
    for (size_t b = 0; b < nb; ++b) {
        const MovingSolid& s = moving_[b];
        if (s.resting && b + 1 < solidCells_.size()) { // (the last list is the liquid)
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

bool NSGridSolver::solidFaceVelocity(long ca, long cb, const Vector3& fw, int comp, float& out) const {
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

void NSGridSolver::applyPendingImpulses() {
    // Velocity change of the faces around the point, trilinear weights (they sum to 1, so the
    // momentum put in is exactly the impulse: every face stands for one cell of gas).
    const float cellMass = params.fluidDensity * dx_ * dx_ * dx_;
    Field3* comps[3] = {&u_, &v_, &w_};
    for (const auto& [x, J] : pendingImpulses_) {
        const Vector3 dv = J / cellMass;
        for (int a = 0; a < 3; ++a) {
            Field3& F = *comps[a];
            const Vector3 g = (x - origin_) / dx_ - F.offset;
            const int i0 = int(std::floor(g.x)), j0 = int(std::floor(g.y)), k0 = int(std::floor(g.z));
            const float tx = g.x - i0, ty = g.y - j0, tz = g.z - k0;
            for (int c = 0; c < 8; ++c) {
                int i = i0 + (c & 1), j = j0 + ((c >> 1) & 1), k = k0 + ((c >> 2) & 1);
                if (i < 0 || j < 0 || k < 0 || i >= F.nx || j >= F.ny || k >= F.nz) continue;
                float w = ((c & 1) ? tx : 1 - tx) * (((c >> 1) & 1) ? ty : 1 - ty) * (((c >> 2) & 1) ? tz : 1 - tz);
                F.at(i, j, k) += dv[a] * w;
            }
        }
    }
    pendingImpulses_.clear();
}

void NSGridSolver::computeMovingForces() {
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

void NSGridSolver::applyWallFriction(float dt) {
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
            const float Re = std::max(mag * L / nu, 1.0f);
            const float cf = Re < 5e5f ? 1.328f / std::sqrt(Re) : 0.074f * std::pow(Re, -0.2f);
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

// ---------------------------------------------------------------------------
// Advection
// ---------------------------------------------------------------------------
Vector3 NSGridSolver::sampleVelGrid(const Vector3& gp) const { return sampleVel(u0_, v0_, w0_, gp) / dx_; }

Vector3 NSGridSolver::backtrace(const Vector3& gp, float dt) const {
    Vector3 v1 = sampleVelGrid(gp);
    Vector3 mid = gp - v1 * (0.5f * dt);
    Vector3 v2 = sampleVelGrid(mid);
    return gp - v2 * dt;
}

void NSGridSolver::advect(const Field3& src, Field3& dst, float dt, bool scalar) const {
    dst.nx = src.nx; dst.ny = src.ny; dst.nz = src.nz; dst.offset = src.offset;
    dst.d.resize(src.d.size());
    const int sx = src.nx, sy = src.ny, sz = src.nz;
    const float size[3] = {float(nx_), float(ny_), float(nz_)};
    // Did the backtraced point come in from outside through an open side of the domain?
    auto fromOutside = [&](const Vector3& b) {
        for (int a = 0; a < 3; ++a) {
            if (b[a] < 0 && params.bc[2 * a] != BoundaryType::Wall) return true;
            if (b[a] > size[a] && params.bc[2 * a + 1] != BoundaryType::Wall) return true;
        }
        return false;
    };
    parallelFor(int(sz), [&](int k_) {
        int k = k_;
        for (int j = 0; j < sy; ++j)
            for (int i = 0; i < sx; ++i) {
                Vector3 gp = Vector3(float(i), float(j), float(k)) + src.offset;
                Vector3 b = backtrace(gp, dt);
                dst.at(i, j, k) = scalar && fromOutside(b) ? 0.0f : src.sample(b);
            }
    }, 1);
}

void NSGridSolver::advectField(Field3& f, Field3& t1, Field3& t2, float dt, bool scalar) {
    advect(f, t1, dt, scalar);
    if (!params.maccormack) {
        f.d.swap(t1.d);
        return;
    }
    advect(t1, t2, -dt, scalar);
    const int sx = f.nx, sy = f.ny, sz = f.nz;
    parallelFor(int(sz), [&](int k_) {
        int k = k_;
        for (int j = 0; j < sy; ++j)
            for (int i = 0; i < sx; ++i) {
                size_t id = f.idx(i, j, k);
                float val = t1.d[id] + 0.5f * (f.d[id] - t2.d[id]);
                float lo, hi;
                f.sampleMinMax(backtrace(Vector3(float(i), float(j), float(k)) + f.offset, dt), lo, hi);
                t2.d[id] = (val < lo || val > hi) ? t1.d[id] : val;
            }
    }, 1);
    f.d.swap(t2.d);
}

void NSGridSolver::advectScalars(const std::vector<Field3*>& fields, float dt) {
    // All cell-centred scalars ride the same flow: the departure points (and the arrival points
    // for MacCormack's backward step) are traced once for all of them. What arrives from outside
    // through an open side is clean ambient gas (0).
    const size_t n = size_t(nx_) * ny_ * nz_;
    departure_.resize(n);
    arrival_.resize(n);
    fromOutside_.resize(n);
    const float size[3] = {float(nx_), float(ny_), float(nz_)};
    auto outside = [&](const Vector3& b) {
        for (int a = 0; a < 3; ++a) {
            if (b[a] < 0 && params.bc[2 * a] != BoundaryType::Wall) return true;
            if (b[a] > size[a] && params.bc[2 * a + 1] != BoundaryType::Wall) return true;
        }
        return false;
    };
    parallelFor(nz_, [&](int k) {
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                const size_t c = cidx(i, j, k);
                const Vector3 gp(i + 0.5f, j + 0.5f, k + 0.5f);
                departure_[c] = backtrace(gp, dt);
                arrival_[c] = backtrace(gp, -dt);
                fromOutside_[c] = uint8_t(outside(departure_[c]) | (outside(arrival_[c]) << 1));
            }
    }, 1);
    for (Field3* field : fields) {
        Field3& f = *field;
        t1_.init(nx_, ny_, nz_, f.offset);
        t2_.init(nx_, ny_, nz_, f.offset);
        parallelFor(int(n), [&](int c) { t1_.d[c] = (fromOutside_[c] & 1) ? 0.0f : f.sample(departure_[c]); }, 4096);
        if (params.maccormack) {
            // MacCormack: trace the result back, correct by half the round-trip error, fall back to
            // the plain value where the correction leaves the range of the departure cell (limiter).
            parallelFor(int(n), [&](int c) { t2_.d[c] = (fromOutside_[c] & 2) ? 0.0f : t1_.sample(arrival_[c]); }, 4096);
            parallelFor(int(n), [&](int c) {
                if (fromOutside_[c] & 1) return;
                const float val = t1_.d[c] + 0.5f * (f.d[c] - t2_.d[c]);
                float lo, hi;
                f.sampleMinMax(departure_[c], lo, hi);
                if (val >= lo && val <= hi) t1_.d[c] = val;
            }, 4096);
        }
        f.d.swap(t1_.d);
    }
}

// ---------------------------------------------------------------------------
// Boundary conditions
// ---------------------------------------------------------------------------
void NSGridSolver::applyVelocityBC() {
    const float U = params.inflowSpeed;
    const BoundaryType* bc = params.bc;
    auto boundaryValue = [&](BoundaryType t, float interior, bool minSide) {
        switch (t) {
        case BoundaryType::Wall: return 0.0f;
        case BoundaryType::Inflow: return minSide ? U : -U;
        case BoundaryType::Outflow: return minSide ? std::min(interior, 0.0f) : std::max(interior, 0.0f);
        }
        return 0.0f;
    };
    parallelFor(int(nz_), [&](int k_) {
        int k = k_;
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i <= nx_; ++i) {
                float val;
                if (solidFaceVelocity(i > 0 ? long(cidx(i - 1, j, k)) : -1, i < nx_ ? long(cidx(i, j, k)) : -1,
                                      origin_ + Vector3(float(i), j + 0.5f, k + 0.5f) * dx_, 0, val))
                    u_.at(i, j, k) = val;
            }
            u_.at(0, j, k) = boundaryValue(bc[0], u_.at(1, j, k), true);
            u_.at(nx_, j, k) = boundaryValue(bc[1], u_.at(nx_ - 1, j, k), false);
        }
    }, 1);
    parallelFor(int(nz_), [&](int k_) {
        int k = k_;
        for (int i = 0; i < nx_; ++i) {
            for (int j = 0; j <= ny_; ++j) {
                float val;
                if (solidFaceVelocity(j > 0 ? long(cidx(i, j - 1, k)) : -1, j < ny_ ? long(cidx(i, j, k)) : -1,
                                      origin_ + Vector3(i + 0.5f, float(j), k + 0.5f) * dx_, 1, val))
                    v_.at(i, j, k) = val;
            }
            v_.at(i, 0, k) = boundaryValue(bc[2], v_.at(i, 1, k), true);
            v_.at(i, ny_, k) = boundaryValue(bc[3], v_.at(i, ny_ - 1, k), false);
        }
    }, 1);
    parallelFor(int(ny_), [&](int j_) {
        int j = j_;
        for (int i = 0; i < nx_; ++i) {
            for (int k = 0; k <= nz_; ++k) {
                float val;
                if (solidFaceVelocity(k > 0 ? long(cidx(i, j, k - 1)) : -1, k < nz_ ? long(cidx(i, j, k)) : -1,
                                      origin_ + Vector3(i + 0.5f, j + 0.5f, float(k)) * dx_, 2, val))
                    w_.at(i, j, k) = val;
            }
            w_.at(i, j, 0) = boundaryValue(bc[4], w_.at(i, j, 1), true);
            w_.at(i, j, nz_) = boundaryValue(bc[5], w_.at(i, j, nz_ - 1), false);
        }
    }, 1);
}

// ---------------------------------------------------------------------------
// Forces
// ---------------------------------------------------------------------------
Vector3 NSGridSolver::cellVelocity(int i, int j, int k) const {
    return {0.5f * (u_.at(i, j, k) + u_.at(i + 1, j, k)), 0.5f * (v_.at(i, j, k) + v_.at(i, j + 1, k)),
            0.5f * (w_.at(i, j, k) + w_.at(i, j, k + 1))};
}

static Vector3 curlAt(const NSGridSolver& s, int i, int j, int k, Vector3 (NSGridSolver::*cv)(int, int, int) const,
                   int nx, int ny, int nz, float dx) {
    int im = std::max(i - 1, 0), ip = std::min(i + 1, nx - 1);
    int jm = std::max(j - 1, 0), jp = std::min(j + 1, ny - 1);
    int km = std::max(k - 1, 0), kp = std::min(k + 1, nz - 1);
    Vector3 dX = ((s.*cv)(ip, j, k) - (s.*cv)(im, j, k)) / (float(ip - im) * dx + 1e-12f);
    Vector3 dY = ((s.*cv)(i, jp, k) - (s.*cv)(i, jm, k)) / (float(jp - jm) * dx + 1e-12f);
    Vector3 dZ = ((s.*cv)(i, j, kp) - (s.*cv)(i, j, km)) / (float(kp - km) * dx + 1e-12f);
    return {dY.z - dZ.y, dZ.x - dX.z, dX.y - dY.x};
}

void NSGridSolver::addForces(float dt) {
    const float hb = params.heatBuoyancy, sb = params.smokeBuoyancy;
    const bool fire = combustion.enabled;
    if (hb != 0 || sb != 0 || fire) {
        parallelFor(int(nz_), [&](int k_) {
            int k = k_;
            for (int j = 1; j < ny_; ++j)
                for (int i = 0; i < nx_; ++i) {
                    float T = 0.5f * (temp_.at(i, j - 1, k) + temp_.at(i, j, k));
                    float S = 0.5f * (smoke_.at(i, j - 1, k) + smoke_.at(i, j, k));
                    float lift = fire ? combustion.buoyancy(T) : hb * T; // fire: ideal gas, T in kelvin
                    v_.at(i, j, k) += dt * (lift - sb * S);
                }
        }, 1);
    }
    const float eps = params.vorticityConfinement;
    if (eps > 0) {
        const size_t n = size_t(nx_) * ny_ * nz_;
        std::vector<Vector3> omega(n), f(n, Vector3(0.0f));
        std::vector<float> mag(n);
        parallelFor(int(nz_), [&](int k_) {
            int k = k_;
            for (int j = 0; j < ny_; ++j)
                for (int i = 0; i < nx_; ++i) {
                    Vector3 w = curlAt(*this, i, j, k, &NSGridSolver::cellVelocity, nx_, ny_, nz_, dx_);
                    omega[cidx(i, j, k)] = w;
                    mag[cidx(i, j, k)] = length(w);
                }
        }, 1);
        parallelFor(int(nz_ - 1) - (1), [&](int k_) {
            int k = k_ + (1);
            for (int j = 1; j < ny_ - 1; ++j)
                for (int i = 1; i < nx_ - 1; ++i) {
                    if (solid(i, j, k)) continue;
                    Vector3 g(mag[cidx(i + 1, j, k)] - mag[cidx(i - 1, j, k)], mag[cidx(i, j + 1, k)] - mag[cidx(i, j - 1, k)],
                           mag[cidx(i, j, k + 1)] - mag[cidx(i, j, k - 1)]);
                    float l = length(g);
                    if (l < 1e-9f) continue;
                    f[cidx(i, j, k)] = cross(g / l, omega[cidx(i, j, k)]) * (eps * dx_);
                }
        }, 1);
        parallelFor(int(nz_), [&](int k_) {
            int k = k_;
            for (int j = 0; j < ny_; ++j)
                for (int i = 0; i < nx_; ++i) {
                    const Vector3& fc = f[cidx(i, j, k)];
                    if (i > 0) u_.at(i, j, k) += 0.5f * dt * (fc.x + f[cidx(i - 1, j, k)].x);
                    if (j > 0) v_.at(i, j, k) += 0.5f * dt * (fc.y + f[cidx(i, j - 1, k)].y);
                    if (k > 0) w_.at(i, j, k) += 0.5f * dt * (fc.z + f[cidx(i, j, k - 1)].z);
                }
        }, 1);
    }
}

void NSGridSolver::diffuse(float dt) {
    const float a = params.kinematicViscosity * dt / (dx_ * dx_);
    if (a < 1e-3f) return; // negligible next to the numerical diffusion of advection
    const int iters = 20;
    Field3* comps[3] = {&u_, &v_, &w_};
    for (int c = 0; c < 3; ++c) {
        Field3& F = *comps[c];
        t1_ = F; // right-hand side (value before diffusion)
        for (int it = 0; it < iters; ++it) {
            t2_ = F;
            const int sx = F.nx, sy = F.ny, sz = F.nz;
            parallelFor(int(sz), [&](int k_) {
                int k = k_;
                for (int j = 0; j < sy; ++j)
                    for (int i = 0; i < sx; ++i) {
                        // boundary-normal faces are fixed by the boundary conditions
                        if ((c == 0 && (i == 0 || i == sx - 1)) || (c == 1 && (j == 0 || j == sy - 1)) ||
                            (c == 2 && (k == 0 || k == sz - 1)))
                            continue;
                        size_t id = F.idx(i, j, k);
                        float cur = t2_.d[id];
                        auto nb = [&](int a2, int b2, int c2) {
                            if (a2 < 0 || b2 < 0 || c2 < 0 || a2 >= sx || b2 >= sy || c2 >= sz) return cur;
                            return t2_.at(a2, b2, c2);
                        };
                        float sum = nb(i - 1, j, k) + nb(i + 1, j, k) + nb(i, j - 1, k) + nb(i, j + 1, k) +
                                    nb(i, j, k - 1) + nb(i, j, k + 1);
                        F.d[id] = (t1_.d[id] + a * sum) / (1.0f + 6.0f * a);
                    }
            }, 1);
            applyVelocityBC(); // keeps solid faces at zero (no-slip)
        }
    }
}

// ---------------------------------------------------------------------------
// Pressure projection (PCG)
// ---------------------------------------------------------------------------
void NSGridSolver::project(float dt) {
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

    auto applyA = [&](const std::vector<double>& x, std::vector<double>& out) {
        parallelFor(int(NZ), [&](int k_) {
            int k = k_;
            for (int j = 0; j < NY; ++j)
                for (int i = 0; i < NX; ++i) {
                    size_t c = cidx(i, j, k);
                    if (solid_[c] || diag_[c] == 0) { out[c] = 0; continue; }
                    double s = diag_[c] * x[c];
                    if (i > 0 && !solid_[c - 1]) s -= x[c - 1];
                    if (i < NX - 1 && !solid_[c + 1]) s -= x[c + 1];
                    if (j > 0 && !solid_[c - NX]) s -= x[c - NX];
                    if (j < NY - 1 && !solid_[c + NX]) s -= x[c + NX];
                    size_t sl = size_t(NX) * NY;
                    if (k > 0 && !solid_[c - sl]) s -= x[c - sl];
                    if (k < NZ - 1 && !solid_[c + sl]) s -= x[c + sl];
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
            z_[c] = diag_[c] ? r_[c] / diag_[c] : 0.0;
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
                z_[c] = diag_[c] ? r_[c] / diag_[c] : 0.0;
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
                if (fluid(i - 1, j, k) && fluid(i, j, k)) u_.at(i, j, k) -= float(qAt(i, j, k) - qAt(i - 1, j, k));
            if (bc[0] == BoundaryType::Outflow && fluid(0, j, k)) u_.at(0, j, k) -= float(qAt(0, j, k));
            if (bc[1] == BoundaryType::Outflow && fluid(NX - 1, j, k)) u_.at(NX, j, k) += float(qAt(NX - 1, j, k));
        }
    }, 1);
    parallelFor(int(NZ), [&](int k_) {
        int k = k_;
        for (int i = 0; i < NX; ++i) {
            for (int j = 1; j < NY; ++j)
                if (fluid(i, j - 1, k) && fluid(i, j, k)) v_.at(i, j, k) -= float(qAt(i, j, k) - qAt(i, j - 1, k));
            if (bc[2] == BoundaryType::Outflow && fluid(i, 0, k)) v_.at(i, 0, k) -= float(qAt(i, 0, k));
            if (bc[3] == BoundaryType::Outflow && fluid(i, NY - 1, k)) v_.at(i, NY, k) += float(qAt(i, NY - 1, k));
        }
    }, 1);
    parallelFor(int(NY), [&](int j_) {
        int j = j_;
        for (int i = 0; i < NX; ++i) {
            for (int k = 1; k < NZ; ++k)
                if (fluid(i, j, k - 1) && fluid(i, j, k)) w_.at(i, j, k) -= float(qAt(i, j, k) - qAt(i, j, k - 1));
            if (bc[4] == BoundaryType::Outflow && fluid(i, j, 0)) w_.at(i, j, 0) -= float(qAt(i, j, 0));
            if (bc[5] == BoundaryType::Outflow && fluid(i, j, NZ - 1)) w_.at(i, j, NZ) += float(qAt(i, j, NZ - 1));
        }
    }, 1);

    const float toP = 1.0f / toQ;
    parallelFor(int(long(n)), [&](int c_) {
        long c = c_; p_.d[c] = float(q_[c]) * toP;
    }, 4096);
}

void NSGridSolver::computeForces() {
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

void NSGridSolver::injectSources() {
    if (params.smokeRake && params.bc[0] == BoundaryType::Inflow) {
        int k = nz_ / 2;
        for (int j = 1; j < ny_ - 1; ++j)
            if (j % 3 == 1)
                for (int i = 0; i < std::min(2, nx_); ++i) smoke_.at(i, j, k) = 1.0f;
    }
    if (source.enabled) {
        Vector3 c = (source.center - origin_) / dx_;
        float r = source.radius / dx_;
        int i0 = std::max(0, int(c.x - r)), i1 = std::min(nx_ - 1, int(c.x + r) + 1);
        int j0 = std::max(0, int(c.y - r)), j1 = std::min(ny_ - 1, int(c.y + r) + 1);
        int k0 = std::max(0, int(c.z - r)), k1 = std::min(nz_ - 1, int(c.z + r) + 1);
        for (int k = k0; k <= k1; ++k)
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i) {
                    if (solid(i, j, k)) continue;
                    if (length2(Vector3(i + 0.5f, j + 0.5f, k + 0.5f) - c) > r * r) continue;
                    temp_.at(i, j, k) = std::max(temp_.at(i, j, k), source.temperature);
                    smoke_.at(i, j, k) = std::max(smoke_.at(i, j, k), source.smoke);
                    fuel_.at(i, j, k) = std::max(fuel_.at(i, j, k), source.fuel);
                }
        if (length2(source.velocity) > 0) { // a jet: the gas in the source moves with it
            Disturbance jet;
            jet.center = source.center;
            jet.radius = source.radius;
            jet.velocity = source.velocity;
            applyDisturbance(jet);
        }
    }
    applyPendingEmissions();
}

void NSGridSolver::conductHeat(float dt) {
    // Fourier's law in the gas: dT/dt = div(alpha grad T). The diffusivity of air grows with
    // temperature (~T^1.75); on a face it is the mean of its two cells. No flux into solids or
    // through the walls. Explicit, in as many sub-steps as stability needs (6 alpha h / dx^2 <= 1).
    const float Tair = combustion.ambientTemperature;
    const size_t n = temp_.d.size();
    std::vector<float> alpha(n);
    const float alphaMax = parallelMax<float>(int(n), 0.0f, [&](int b, int e) {
        float m = 0;
        for (int c = b; c < e; ++c) {
            alpha[c] = combustion.thermalDiffusivity * std::pow((temp_.d[c] + Tair) / Tair, 1.75f);
            m = std::max(m, alpha[c]);
        }
        return m;
    }, 4096);
    const int steps = std::max(1, int(std::ceil(6.0f * alphaMax * dt / (dx_ * dx_) / 0.9f)));
    const float h = dt / float(steps) / (dx_ * dx_);
    const long sx = 1, sy = nx_, sz = long(nx_) * ny_;
    std::vector<float> T;
    for (int s = 0; s < steps; ++s) {
        T = temp_.d;
        parallelFor(nz_, [&](int k) {
            for (int j = 0; j < ny_; ++j)
                for (int i = 0; i < nx_; ++i) {
                    const size_t c = cidx(i, j, k);
                    if (solid_[c]) continue;
                    float flux = 0;
                    auto face = [&](bool inside, long step) {
                        if (!inside || solid_[c + step]) return;
                        flux += 0.5f * (alpha[c] + alpha[c + step]) * (T[c + step] - T[c]);
                    };
                    face(i > 0, -sx); face(i < nx_ - 1, sx);
                    face(j > 0, -sy); face(j < ny_ - 1, sy);
                    face(k > 0, -sz); face(k < nz_ - 1, sz);
                    temp_.d[c] = T[c] + h * flux;
                }
        }, 1);
    }
}

void NSGridSolver::collectRadiators() {
    const float sigma = 5.670e-8f, Tair = combustion.ambientTemperature;
    const float kappaV = combustion.absorptionCoefficient * dx_ * dx_ * dx_;
    radiators_.clear();
    for (int k = 0; k < nz_; ++k)
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                const float T = temp_.at(i, j, k) + Tair;
                if (T < 700.0f) continue; // below that the flame's radiation is negligible
                const float power = 4.0f * kappaV * sigma * (T * T * T * T - Tair * Tair * Tair * Tair);
                radiators_.push_back({origin_ + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * dx_, power});
            }
}

float NSGridSolver::irradianceAt(const Vector3& world) const {
    const float minR2 = dx_ * dx_; // a point source is only fair from about a cell away
    float q = 0;
    for (const Radiator& r : radiators_) q += r.power / (4.0f * kPi * std::max(length2(r.position - world), minR2));
    return q;
}

void NSGridSolver::applyPendingEmissions() {
    const float cellVolume = dx_ * dx_ * dx_;
    const float heatCapacity = params.fluidDensity * combustion.specificHeat * cellVolume; // [J/K] of a cell
    for (const GasEmission& e : pendingEmissions_) {
        const Vector3 g = (e.position - origin_) / dx_;
        const int i = std::clamp(int(std::floor(g.x)), 0, nx_ - 1);
        const int j = std::clamp(int(std::floor(g.y)), 0, ny_ - 1);
        const int k = std::clamp(int(std::floor(g.z)), 0, nz_ - 1);
        if (solid(i, j, k)) continue;
        fuel_.at(i, j, k) += e.fuel / cellVolume;
        smoke_.at(i, j, k) += e.smoke / cellVolume;
        temp_.at(i, j, k) = std::max(0.0f, temp_.at(i, j, k) + e.heat / heatCapacity);
    }
    pendingEmissions_.clear();
}

// ---------------------------------------------------------------------------
// Time step
// ---------------------------------------------------------------------------
float NSGridSolver::step(float maxDt) {
    auto absMax = [](const std::vector<float>& d) {
        return parallelMax<float>(int(d.size()), 0.0f, [&](int b, int e) {
            float m = 0;
            for (int i = b; i < e; ++i) m = std::max(m, std::fabs(d[i]));
            return m;
        }, 4096);
    };
    float umax = std::max({absMax(u_.d), absMax(v_.d), absMax(w_.d)});
    maxVel_ = umax;
    float vref = std::max({umax, params.inflowSpeed, 0.05f});
    float dt = std::min(maxDt, params.cfl * dx_ / vref);
    if (magnetic.enabled) dt = std::min(dt, magnetic.maxTimeStep(umax, params.fluidDensity)); // Alfven waves

    using Clock = std::chrono::steady_clock;
    auto ms = [](Clock::time_point a, Clock::time_point b) { return std::chrono::duration<float, std::milli>(b - a).count(); };
    auto t0 = Clock::now();
    voxelizeMovingSolids();
    solidMs_ = ms(t0, Clock::now());
    applyVelocityBC();
    injectSources();

    u0_ = u_; v0_ = v_; w0_ = w_;
    advectField(u_, t1_, t2_, dt);
    advectField(v_, t1_, t2_, dt);
    advectField(w_, t1_, t2_, dt);
    std::vector<Field3*> scalars = {&smoke_};
    if (source.enabled || params.heatBuoyancy != 0 || combustion.enabled || magnetic.enabled) scalars.push_back(&temp_);
    if (combustion.enabled) {
        scalars.push_back(&fuel_);
        scalars.push_back(&products_); // fresh air (0) comes in through the open sides
    }
    advectScalars(scalars, dt);
    if (combustion.enabled) {
        const double burnt = combustion.react(fuel_.d, products_.d, temp_.d, smoke_.d, expansion_, solid_, dt);
        // Power of the flame: every unit of burnt fuel heated its cell by heatRelease.
        const float cellHeatCapacity = params.fluidDensity * combustion.specificHeat * dx_ * dx_ * dx_;
        heatReleaseRate_ = float(burnt * combustion.heatRelease * cellHeatCapacity / dt);
        conductHeat(dt);
        collectRadiators();
    }
    applyVelocityBC();

    addForces(dt);
    if (magnetic.enabled) magnetic.applyLorentzForce(u_, v_, w_, solid_, params.fluidDensity, dt);
    applyPendingImpulses();
    diffuse(dt);
    applyVelocityBC();
    auto t1 = Clock::now();
    applyWallFriction(dt);
    auto t2 = Clock::now();
    project(dt);
    if (magnetic.enabled) {
        // Faraday with the new, divergence-free flow; the current's Joule heat warms the gas.
        magnetic.induce(u_, v_, w_, params.fluidDensity, dt);
        std::vector<float> joule;
        magnetic.jouleHeating(joule, dt);
        const float rhoCp = params.fluidDensity * combustion.specificHeat;
        for (size_t c = 0; c < joule.size(); ++c)
            if (!solid_[c]) temp_.d[c] += joule[c] / rhoCp;
    }
    auto t3 = Clock::now();
    computeForces();
    computeMovingForces();
    solidMs_ += ms(t1, t2) + ms(t3, Clock::now());
    pressureMs_ = ms(t2, t3);
    computeDiagnostics();

    if (params.smokeDissipation > 0) {
        float f = std::exp(-params.smokeDissipation * dt);
        for (float& s : smoke_.d) s *= f;
    }
    if (params.temperatureDissipation > 0) {
        float f = std::exp(-params.temperatureDissipation * dt);
        for (float& t : temp_.d) t *= f;
    }
    lastDt_ = dt;

    // Exponential moving average of the coefficients over ~1 flow-through time.
    float tau = std::max(params.domainSize.x / std::max(params.inflowSpeed, 0.1f), 0.05f);
    float a = std::min(1.0f, dt / tau);
    if (time_ == 0) { cdAvg_ = cd_; clAvg_ = cl_; }
    else { cdAvg_ += a * (cd_ - cdAvg_); clAvg_ += a * (cl_ - clAvg_); }
    time_ += dt;
    return dt;
}

void NSGridSolver::applyDisturbance(const Disturbance& d) {
    const float R = std::max(d.radius, 0.5f * dx_);
    const Vector3 cg = (d.center - origin_) / dx_; // centre in cell units
    const float rg = R / dx_;
    auto weight = [&](const Vector3& gp) {
        float q = length2(gp - cg) / (rg * rg);
        return q < 1.0f ? (1.0f - q) * (1.0f - q) : 0.0f;
    };
    // Velocity components live on faces: blend each within its own staggered layout.
    Field3* comps[3] = {&u_, &v_, &w_};
    for (int c = 0; c < 3; ++c) {
        Field3& F = *comps[c];
        int i0 = std::max(0, int(cg.x - rg) - 1), i1 = std::min(F.nx - 1, int(cg.x + rg) + 1);
        int j0 = std::max(0, int(cg.y - rg) - 1), j1 = std::min(F.ny - 1, int(cg.y + rg) + 1);
        int k0 = std::max(0, int(cg.z - rg) - 1), k1 = std::min(F.nz - 1, int(cg.z + rg) + 1);
        for (int k = k0; k <= k1; ++k)
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i) {
                    float w = weight(Vector3(float(i), float(j), float(k)) + F.offset);
                    if (w > 0) F.at(i, j, k) += (d.velocity[c] - F.at(i, j, k)) * (w * d.velocityBlend);
                }
    }
    if (d.smoke > 0 || d.heat > 0) {
        int i0 = std::max(0, int(cg.x - rg)), i1 = std::min(nx_ - 1, int(cg.x + rg) + 1);
        int j0 = std::max(0, int(cg.y - rg)), j1 = std::min(ny_ - 1, int(cg.y + rg) + 1);
        int k0 = std::max(0, int(cg.z - rg)), k1 = std::min(nz_ - 1, int(cg.z + rg) + 1);
        for (int k = k0; k <= k1; ++k)
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i) {
                    if (solid(i, j, k)) continue;
                    float w = weight(Vector3(i + 0.5f, j + 0.5f, k + 0.5f));
                    if (w <= 0) continue;
                    smoke_.at(i, j, k) = std::min(1.0f, smoke_.at(i, j, k) + d.smoke * w);
                    temp_.at(i, j, k) += d.heat * w;
                }
    }
    applyVelocityBC(); // walls and solids stay impermeable
}

void NSGridSolver::computeDiagnostics() {
    // Max |div u| over fluid cells [1/s] and total smoke volume [m^3 * density].
    struct D {
        float maxDiv = 0;
        double smoke = 0;
        D& operator+=(const D& o) { maxDiv = std::max(maxDiv, o.maxDiv); smoke += o.smoke; return *this; }
    };
    D d = parallelSum<D>(nz_, [&](int k0, int k1) {
        D acc;
        for (int k = k0; k < k1; ++k)
            for (int j = 0; j < ny_; ++j)
                for (int i = 0; i < nx_; ++i) {
                    if (solid(i, j, k)) continue;
                    acc.smoke += smoke_.at(i, j, k);
                    if (diag_[cidx(i, j, k)] == 0) continue;
                    // The uniform part of a closed region (gas squeezed by a moving body) cannot be
                    // projected away; measure the part the projection is responsible for.
                    int r = region_[cidx(i, j, k)];
                    double m = (r >= 0 && r < int(regionMeanB_.size())) ? regionMeanB_[r] : 0.0;
                    float div = float(u_.at(i + 1, j, k) - u_.at(i, j, k) + v_.at(i, j + 1, k) - v_.at(i, j, k) +
                                      w_.at(i, j, k + 1) - w_.at(i, j, k) + m) / dx_;
                    acc.maxDiv = std::max(acc.maxDiv, std::fabs(div));
                }
        return acc;
    }, 1);
    maxDivergence_ = d.maxDiv;
    totalSmoke_ = float(d.smoke * dx_ * dx_ * dx_);
}

// ---------------------------------------------------------------------------
// Queries
// ---------------------------------------------------------------------------
bool NSGridSolver::touchesGas(const Vector3& x) const {
    Vector3 g = (x - origin_) / dx_;
    int i = int(std::floor(g.x)), j = int(std::floor(g.y)), k = int(std::floor(g.z));
    const int nb[7][3] = {{0, 0, 0}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    for (auto& d : nb) {
        int a = i + d[0], b = j + d[1], c = k + d[2];
        if (a < 0 || b < 0 || c < 0 || a >= nx_ || b >= ny_ || c >= nz_) continue;
        if (!solid(a, b, c)) return true;
    }
    return false;
}

int NSGridSolver::fluidStencil(const Vector3& world, size_t cells[8], float weights[8]) const {
    const Vector3 g = (world - origin_) / dx_ - Vector3(0.5f); // in cell-centre coordinates
    const int i0 = int(std::floor(g.x)), j0 = int(std::floor(g.y)), k0 = int(std::floor(g.z));
    const float tx = g.x - i0, ty = g.y - j0, tz = g.z - k0;
    int n = 0;
    float sum = 0;
    for (int c = 0; c < 8; ++c) {
        int i = i0 + (c & 1), j = j0 + ((c >> 1) & 1), k = k0 + ((c >> 2) & 1);
        if (i < 0 || j < 0 || k < 0 || i >= nx_ || j >= ny_ || k >= nz_ || solid(i, j, k)) continue;
        float w = ((c & 1) ? tx : 1 - tx) * (((c >> 1) & 1) ? ty : 1 - ty) * (((c >> 2) & 1) ? tz : 1 - tz);
        if (w <= 0) continue;
        cells[n] = cidx(i, j, k);
        weights[n] = w;
        sum += w;
        ++n;
    }
    if (n > 0) {
        for (int m = 0; m < n; ++m) weights[m] /= sum;
        return n;
    }
    // Deep inside a solid corner: nearest gas cell.
    const int ic = int(std::lround(g.x)), jc = int(std::lround(g.y)), kc = int(std::lround(g.z));
    float best = kInf;
    for (int k = kc - 2; k <= kc + 2; ++k)
        for (int j = jc - 2; j <= jc + 2; ++j)
            for (int i = ic - 2; i <= ic + 2; ++i) {
                if (i < 0 || j < 0 || k < 0 || i >= nx_ || j >= ny_ || k >= nz_ || solid(i, j, k)) continue;
                float d = length2(Vector3(float(i), float(j), float(k)) - g);
                if (d < best) { best = d; cells[0] = cidx(i, j, k); }
            }
    if (best == kInf) return 0;
    weights[0] = 1;
    return 1;
}

float NSGridSolver::fluidPressureAt(const Vector3& world) const {
    size_t c[8];
    float w[8];
    int n = fluidStencil(world, c, w);
    float p = 0;
    for (int m = 0; m < n; ++m) p += w[m] * p_.d[c[m]];
    return p;
}

Vector3 NSGridSolver::fluidVelocityAt(const Vector3& world) const {
    size_t c[8];
    float w[8];
    int n = fluidStencil(world, c, w);
    Vector3 u(0.0f);
    for (int m = 0; m < n; ++m) {
        int i = int(c[m] % nx_), j = int((c[m] / nx_) % ny_), k = int(c[m] / (size_t(nx_) * ny_));
        u += cellVelocity(i, j, k) * w[m];
    }
    return u;
}

Vector3 NSGridSolver::velocityAt(const Vector3& world) const { return sampleVel(u_, v_, w_, (world - origin_) / dx_); }

float NSGridSolver::cellValue(GridField f, int i, int j, int k) const {
    switch (f) {
    case GridField::Speed: return length(cellVelocity(i, j, k));
    case GridField::Pressure: return p_.at(i, j, k);
    case GridField::PressureCoeff: return p_.at(i, j, k) / dynamicPressure();
    case GridField::Vorticity: return length(curlAt(*this, i, j, k, &NSGridSolver::cellVelocity, nx_, ny_, nz_, dx_));
    case GridField::Smoke: return smoke_.at(i, j, k);
    case GridField::Temperature: return temp_.at(i, j, k);
    case GridField::VelocityX: return cellVelocity(i, j, k).x;
    case GridField::VelocityY: return cellVelocity(i, j, k).y;
    case GridField::MagneticFlux: return magnetic.enabled ? length(magnetic.cellField(i, j, k)) : 0.0f;
    case GridField::CurrentDensity: return magnetic.enabled ? length(magnetic.cellCurrent(i, j, k)) : 0.0f;
    }
    return 0;
}

} // namespace rf
