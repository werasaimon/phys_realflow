// GasSolver: the grid, its setup, the time step and the forces on the gas. The other parts live
// in Advection.cpp, PressureSolver.cpp, MovingSolids.cpp and Heat.cpp.
#include "gas/GasSolver.h"
#include "gas/GasSolverInternal.h"

#include "core/Parallel.h"
#include "core/Probe.h"

#include <chrono>

namespace rf {

// ---------------------------------------------------------------------------
// Setup
// ---------------------------------------------------------------------------
// A new grid for the current parameters: its size, its fields, the static obstacle voxelised into
// it, the obstacle's reference areas, the wall cells for the skin friction, then the pressure
// matrix, the magnetic field and the solver's work vectors.
void GasSolver::reset(const Vector3& origin, const MeshBVH* obstacle) {
    origin_ = origin;
    chooseGrid();
    allocateFields();
    voxelizeStaticObstacle(obstacle);
    measureObstacle();
    findStaticWallCells();
    computeDiag();
    magnetic.reset(nx_, ny_, nz_, dx_, origin_);
    if (solidCount_ > 0) magnetic.setConductors(solid_); // the obstacle (a vessel wall) conducts
    const size_t n = size_t(nx_) * ny_ * nz_;
    q_.assign(n, 0); r_.assign(n, 0); z_.assign(n, 0); s_.assign(n, 0); As_.assign(n, 0); b_.assign(n, 0);
    force_ = Vector3(0.0f);
    cd_ = cl_ = cs_ = cdAvg_ = clAvg_ = 0;
    refArea_ = params.usePlanformArea ? planformArea_ : frontalArea_;
    lastIters_ = 0;
    lastResidual_ = 0;
    time_ = 0;
    applyVelocityBC();
}

// The number of cells: resolutionX cells along x fix the spacing dx, the other sides follow the
// domain size; a cap keeps a thin, tall domain from asking for 10^8 cells.
void GasSolver::chooseGrid() {
    nx_ = std::max(8, params.resolutionX);
    dx_ = params.domainSize.x / float(nx_);
    ny_ = std::max(4, int(std::lround(params.domainSize.y / dx_)));
    nz_ = std::max(4, int(std::lround(params.domainSize.z / dx_)));
    // At most 16 million cells (~3 GB with all fields, as 400 x 200 x 200): a thin, tall domain would
    // otherwise ask for 10^8 cells and fail to allocate. The spacing then grows to fit.
    const double kMaxCells = 16e6, cells = double(nx_) * ny_ * nz_;
    if (cells > kMaxCells) {
        nx_ = std::max(8, int(float(nx_) / float(std::cbrt(cells / kMaxCells))));
        dx_ = params.domainSize.x / float(nx_);
        ny_ = std::max(4, int(std::lround(params.domainSize.y / dx_)));
        nz_ = std::max(4, int(std::lround(params.domainSize.z / dx_)));
    }
}

// The MAC grid: velocities on the faces (u on x-faces, v on y-faces, w on z-faces), pressure,
// smoke, temperature, fuel and products in the cell centres; an inflow starts the u faces at the
// inflow speed. Everything handed in for the old scene goes: the bodies, the liquid, pending
// impulses and emissions (a later scene without bodies would otherwise keep the old water as
// solid cells).
void GasSolver::allocateFields() {
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
    pendingImpulses_.clear();
    radiators_.clear();
    moving_.clear();
    liquidPos_.clear();
    liquidVelIn_.clear();
    liquidVel_.clear();
    liquidCells_ = 0;
}

// The static obstacle: every cell whose centre is inside the mesh is a wall (kStatic); with a
// vessel, every cell outside it. The moving-solid bookkeeping starts empty.
void GasSolver::voxelizeStaticObstacle(const MeshBVH* obstacle) {
    const size_t n = size_t(nx_) * ny_ * nz_;
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
    if (vessel) { // a closed vessel: everything outside it is wall
        parallelFor(int(nz_), [&](int k) {
            for (int j = 0; j < ny_; ++j)
                for (int i = 0; i < nx_; ++i) {
                    const Vector3 c = origin_ + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * dx_;
                    if (!vessel(c)) solid_[cidx(i, j, k)] = kStatic;
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
}

// Reference areas of the voxelised body, the projections it presents to the flow (frontal) and
// from above (planform, for a wing), and the length of its boundary layer along the flow.
void GasSolver::measureObstacle() {
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
        // Running length of the boundary layer: the extent along the flow (x), not the largest
        // extent (that is the span of a cylinder or a wing).
        staticLength_ = sb.valid() ? sb.extent().x + dx_ : 1.0f;
    }
}

// Gas cells touching the static obstacle: the only places where it can exert wall friction. The
// friction bookkeeping of the moving bodies starts empty too.
void GasSolver::findStaticWallCells() {
    const size_t n = size_t(nx_) * ny_ * nz_;
    staticFriction_ = Vector3(0.0f);
    movFricForce_.clear();
    movFricTorque_.clear();
    solidCells_.clear();
    stamp_.assign(n, 0);
    stampId_ = 0;
    staticWallCells_.clear();
    for (int k = 0; k < nz_; ++k)
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                if (solid(i, j, k)) continue;
                bool wall = (i > 0 && solid(i - 1, j, k)) || (i < nx_ - 1 && solid(i + 1, j, k)) || (j > 0 && solid(i, j - 1, k)) ||
                            (j < ny_ - 1 && solid(i, j + 1, k)) || (k > 0 && solid(i, j, k - 1)) || (k < nz_ - 1 && solid(i, j, k + 1));
                if (wall) staticWallCells_.push_back(cidx(i, j, k));
            }
}

// Pressure matrix diagonal: fluid neighbours + Dirichlet (outflow) sides.
void GasSolver::computeDiag() {
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
// Boundary conditions
// ---------------------------------------------------------------------------
void GasSolver::applyVelocityBC() {
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
Vector3 GasSolver::cellVelocity(int i, int j, int k) const {
    return {0.5f * (u_.at(i, j, k) + u_.at(i + 1, j, k)), 0.5f * (v_.at(i, j, k) + v_.at(i, j + 1, k)),
            0.5f * (w_.at(i, j, k) + w_.at(i, j, k + 1))};
}

static Vector3 curlAt(const GasSolver& s, int i, int j, int k, Vector3 (GasSolver::*cv)(int, int, int) const,
                   int nx, int ny, int nz, float dx) {
    int im = std::max(i - 1, 0), ip = std::min(i + 1, nx - 1);
    int jm = std::max(j - 1, 0), jp = std::min(j + 1, ny - 1);
    int km = std::max(k - 1, 0), kp = std::min(k + 1, nz - 1);
    Vector3 dX = ((s.*cv)(ip, j, k) - (s.*cv)(im, j, k)) / (float(ip - im) * dx + 1e-12f);
    Vector3 dY = ((s.*cv)(i, jp, k) - (s.*cv)(i, jm, k)) / (float(jp - jm) * dx + 1e-12f);
    Vector3 dZ = ((s.*cv)(i, j, kp) - (s.*cv)(i, j, km)) / (float(kp - km) * dx + 1e-12f);
    return {dY.z - dZ.y, dZ.x - dX.z, dX.y - dY.x};
}

void GasSolver::addForces(float dt) {
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
                    Vector3 w = curlAt(*this, i, j, k, &GasSolver::cellVelocity, nx_, ny_, nz_, dx_);
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

void GasSolver::diffuse(float dt) {
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

void GasSolver::injectSources() {
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

// ---------------------------------------------------------------------------
// Time step
// ---------------------------------------------------------------------------
// One time step of the gas, the operator splitting of Stam 1999 / Bridson 2015 in order: the
// bodies become walls, sources add gas, everything is carried by the flow (advection), the fire
// burns, the forces push, viscosity diffuses, the walls rub, the pressure makes the flow
// divergence-free, the magnetic field follows the flow. Every stage is a named step below; the
// probe measures each.
float GasSolver::step(float maxDt) {
    const float dt = chooseTimeStep(maxDt);
    using Clock = std::chrono::steady_clock;
    auto ms = [](Clock::time_point a, Clock::time_point b) { return std::chrono::duration<float, std::milli>(b - a).count(); };
    auto t0 = Clock::now();
    {
        Probe::Timer timer("gas/solids ms");
        voxelizeMovingSolids();
    }
    solidMs_ = ms(t0, Clock::now());
    applyVelocityBC();
    injectSources();
    {
        Probe::Timer timer("gas/advect ms");
        advectAll(dt);
    }
    burn(dt);
    applyVelocityBC();
    {
        Probe::Timer timer("gas/forces ms"); // buoyancy, confinement, Lorentz, the bodies' impulses
        applyForces(dt);
    }
    {
        Probe::Timer timer("gas/diffuse ms");
        diffuse(dt);
    }
    applyVelocityBC();
    auto t1 = Clock::now();
    applyWallFriction(dt);
    auto t2 = Clock::now();
    {
        Probe::Timer timer("gas/pressure ms");
        project(dt);
    }
    if (magnetic.enabled) {
        Probe::Timer timer("gas/mhd ms");
        induceMagneticField(dt);
    }
    auto t3 = Clock::now();
    computeForces();
    computeMovingForces();
    solidMs_ += ms(t1, t2) + ms(t3, Clock::now());
    pressureMs_ = ms(t2, t3);
    computeDiagnostics();
    reportStep(dt);
    dissipateScalars(dt);
    lastDt_ = dt;
    averageCoefficients(dt);
    time_ += dt;
    return dt;
}

// The CFL condition: no fluid may cross more than `cfl` cells in one step, measured by the fastest
// velocity on the grid (or the inflow); with a magnetic field the Alfven waves set a limit too.
float GasSolver::chooseTimeStep(float maxDt) {
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
    return dt;
}

// Semi-Lagrangian advection (MacCormack) of the velocity itself and of every scalar the scene
// carries: smoke always; temperature when something heats the gas; fuel and products in a fire.
void GasSolver::advectAll(float dt) {
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
}

// The fire: the reaction of fuel with air (Combustion), the power of the flame, heat conduction
// and the radiating cells. A fire that was switched off must not leave its last expansion as a
// source of the pressure projection.
void GasSolver::burn(float dt) {
    if (combustion.enabled) {
        Probe::Timer timer("gas/heat ms"); // the reaction, conduction and the radiators
        const double burnt = combustion.react(fuel_.d, products_.d, temp_.d, smoke_.d, expansion_, solid_, dt);
        // Power of the flame: every unit of burnt fuel heated its cell by heatRelease.
        const float cellHeatCapacity = params.fluidDensity * combustion.specificHeat * dx_ * dx_ * dx_;
        heatReleaseRate_ = float(burnt * combustion.heatRelease * cellHeatCapacity / dt);
        conductHeat(dt);
        collectRadiators();
    } else if (heatReleaseRate_ != 0 || !radiators_.empty()) {
        std::fill(expansion_.begin(), expansion_.end(), 0.0f);
        heatReleaseRate_ = 0;
        radiators_.clear();
    }
}

// Body forces on the gas: buoyancy and vorticity confinement (addForces), the Lorentz force of
// the magnetic field, and the impulses the bodies, cloth and liquid put into it this frame.
void GasSolver::applyForces(float dt) {
    addForces(dt);
    if (magnetic.enabled) magnetic.applyLorentzForce(u_, v_, w_, solid_, params.fluidDensity, dt);
    applyPendingImpulses();
}

// Faraday's law with the new, divergence-free flow (constrained transport in MagneticField); the
// current's Joule heat warms the gas.
void GasSolver::induceMagneticField(float dt) {
    magnetic.induce(u_, v_, w_, params.fluidDensity, dt);
    std::vector<float> joule;
    magnetic.jouleHeating(joule, dt);
    const float rhoCp = params.fluidDensity * combustion.specificHeat;
    for (size_t c = 0; c < joule.size(); ++c)
        if (!solid_[c]) temp_.d[c] += joule[c] / rhoCp;
}

// The step's numbers for the probe.
void GasSolver::reportStep(float dt) const {
    Probe::set("gas/dt", dt);
    Probe::set("gas/pressure iterations", lastIters_);
    Probe::set("gas/pressure residual", lastResidual_);
    Probe::set("gas/max speed", maxVel_);
    Probe::set("gas/cells", double(nx_) * ny_ * nz_);
    if (magnetic.enabled) {
        Probe::set("mhd/max B", magnetic.maxField());
        Probe::set("mhd/energy", magnetic.energy());
        Probe::set("mhd/div B", magnetic.maxDivergence());
    }
}

// Smoke and heat fade with the rates the scene asked for (exponential decay per step).
void GasSolver::dissipateScalars(float dt) {
    if (params.smokeDissipation > 0) {
        float f = std::exp(-params.smokeDissipation * dt);
        for (float& s : smoke_.d) s *= f;
    }
    if (params.temperatureDissipation > 0) {
        float f = std::exp(-params.temperatureDissipation * dt);
        for (float& t : temp_.d) t *= f;
    }
}

// Exponential moving average of the drag and lift coefficients over ~1 flow-through time.
void GasSolver::averageCoefficients(float dt) {
    float tau = std::max(params.domainSize.x / std::max(params.inflowSpeed, 0.1f), 0.05f);
    float a = std::min(1.0f, dt / tau);
    if (time_ == 0) { cdAvg_ = cd_; clAvg_ = cl_; }
    else { cdAvg_ += a * (cd_ - cdAvg_); clAvg_ += a * (cl_ - clAvg_); }
}

void GasSolver::setTracer(const std::function<float(const Vector3&)>& density) {
    parallelFor(int(nz_), [&](int k) {
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                const size_t c = cidx(i, j, k);
                smoke_.d[c] = solid_[c] ? 0.0f : density(origin_ + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * dx_);
            }
    }, 1);
}

void GasSolver::applyDisturbance(const Disturbance& d) {
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

void GasSolver::computeDiagnostics() {
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
                    // projected away, and burning gas expands on purpose (its target divergence):
                    // measure the part the projection is responsible for.
                    const size_t c = cidx(i, j, k);
                    int r = region_[c];
                    double m = (r >= 0 && r < int(regionMeanB_.size())) ? regionMeanB_[r] : 0.0;
                    float div = float(u_.at(i + 1, j, k) - u_.at(i, j, k) + v_.at(i, j + 1, k) - v_.at(i, j, k) +
                                      w_.at(i, j, k + 1) - w_.at(i, j, k) + m - double(expansion_[c]) * dx_) / dx_;
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
bool GasSolver::touchesGas(const Vector3& x) const {
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

int GasSolver::fluidStencil(const Vector3& world, size_t cells[8], float weights[8]) const {
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

float GasSolver::fluidPressureAt(const Vector3& world) const {
    size_t c[8];
    float w[8];
    int n = fluidStencil(world, c, w);
    float p = 0;
    for (int m = 0; m < n; ++m) p += w[m] * p_.d[c[m]];
    return p;
}

Vector3 GasSolver::fluidVelocityAt(const Vector3& world) const {
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

Vector3 GasSolver::velocityAt(const Vector3& world) const { return sampleVel(u_, v_, w_, (world - origin_) / dx_); }

float GasSolver::cellValue(GridField f, int i, int j, int k) const {
    switch (f) {
    case GridField::Speed: return length(cellVelocity(i, j, k));
    case GridField::Pressure: return p_.at(i, j, k);
    case GridField::PressureCoeff: return p_.at(i, j, k) / dynamicPressure();
    case GridField::Vorticity: return length(curlAt(*this, i, j, k, &GasSolver::cellVelocity, nx_, ny_, nz_, dx_));
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
