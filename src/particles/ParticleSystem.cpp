#include "particles/ParticleSystem.h"

#include "core/Parallel.h"

namespace rf {

void ParticleSystem::reset(const AABB& domain) {
    domain_ = domain;
    const float r = params.particleRadius;
    const float s = 2.0f * r;
    h_ = 4.0f * r;
    h2_ = h_ * h_;
    poly6_ = 315.0f / (64.0f * kPi * std::pow(h_, 9.0f));
    spikyGrad_ = -45.0f / (kPi * std::pow(h_, 6.0f));

    // Calibrate the particle mass so that a particle in a full lattice has exactly the rest density.
    double sum = 0;
    for (int i = -3; i <= 3; ++i)
        for (int j = -3; j <= 3; ++j)
            for (int k = -3; k <= 3; ++k) sum += W(length2(Vector3(float(i), float(j), float(k)) * s));
    mass_ = float(params.restDensity / sum);
    deltaQW_ = W(0.04f * h2_); // |dq| = 0.2 h

    x_.clear(); v_.clear(); p_.clear(); dp_.clear(); omega_.clear(); vtmp_.clear();
    phase_.clear(); object_.clear(); invMass_.clear(); volume_.clear(); rest_.clear();
    contacts_.clear();
    softBodies_.clear(); cloths_.clear();
    grab_ = ParticleGrab();
    nextObject_ = 0;
    fluidCount_ = 0;
    rho_.clear(); lambda_.clear(); nbrCount_.clear(); nbr_.clear();
    contactBody_.clear(); contactNormal_.clear(); contactPoint_.clear(); contactDepth_.clear();
    bodyShift_.clear(); bodyTurn_.clear();
    emitter.accumulated = 0;
    avgDensityError_ = maxSpeed_ = 0;

    Vector3 e = domain_.extent();
    gx_ = std::max(1, int(std::ceil(e.x / h_)));
    gy_ = std::max(1, int(std::ceil(e.y / h_)));
    gz_ = std::max(1, int(std::ceil(e.z / h_)));
}

static bool blockedBySolid(const Vector3& p, float r, const MeshBVH* mesh, const RigidWorld* rigid) {
    if (mesh && !mesh->empty()) {
        AABB mb = mesh->bounds();
        mb.lo -= Vector3(r);
        mb.hi += Vector3(r);
        // Unbounded closest-point query: deep inside the obstacle is blocked too (a bounded
        // signedDistance reports "far outside" there).
        ClosestHit hit;
        if (mb.contains(p) && mesh->closestPoint(p, kInf, hit) && hit.signedDistance < r) return true;
    }
    if (rigid)
        for (const RigidBody& b : rigid->bodies()) {
            Vector3 n;
            if (b.signedDistance(p, n) < r) return true;
        }
    return false;
}

void ParticleSystem::addBlock(const AABB& box, const Vector3& vel) {
    const float r = params.particleRadius, s = 2 * r;
    AABB b(vmax(box.lo, domain_.lo + Vector3(r)), vmin(box.hi, domain_.hi - Vector3(r)));
    if (!b.valid()) return;
    int nx = int(b.extent().x / s) + 1, ny = int(b.extent().y / s) + 1, nz = int(b.extent().z / s) + 1;
    uint32_t seed = 12345;
    auto rnd = [&]() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) * (1.0f / 16777216.0f) - 0.5f; };
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                if (int(x_.size()) >= params.maxParticles) return;
                Vector3 p = b.lo + Vector3(float(i), float(j), float(k)) * s + Vector3(rnd(), rnd(), rnd()) * (0.02f * s);
                if (blockedBySolid(p, r, mesh_, rigid_)) continue;
                addParticle(p, vel, ParticlePhase::Fluid, -1, 1.0f / mass_);
            }
}

void ParticleSystem::addParticle(const Vector3& x, const Vector3& v, ParticlePhase phase, int object, float invMass, float volume) {
    x_.push_back(x);
    v_.push_back(v);
    phase_.push_back(uint8_t(phase));
    object_.push_back(object);
    invMass_.push_back(invMass);
    volume_.push_back(volume);
    rest_.push_back(x);
    if (phase == ParticlePhase::Fluid) ++fluidCount_;
}

int ParticleSystem::addSoftBody(const TriMesh& shape, float density, float stiffness, const Vector3& color, const Vector3& velocity) {
    const float s = spacing();
    MeshBVH bvh;
    bvh.build(shape);
    const AABB b = shape.bounds();
    SoftBody body;
    body.object = nextObject_++;
    body.stiffness = stiffness;
    body.color = color;
    const float invMass = 1.0f / (density * s * s * s);
    std::vector<Vector3> rest;
    for (float z = b.lo.z + 0.5f * s; z < b.hi.z; z += s)
        for (float y = b.lo.y + 0.5f * s; y < b.hi.y; y += s)
            for (float x = b.lo.x + 0.5f * s; x < b.hi.x; x += s) {
                Vector3 p(x, y, z);
                if (int(x_.size()) >= params.maxParticles || !domain_.contains(p) || !bvh.isInside(p)) continue;
                body.particles.push_back(int(x_.size()));
                rest.push_back(p);
                addParticle(p, velocity, ParticlePhase::Soft, body.object, invMass);
            }
    if (body.particles.empty()) return -1;
    // Clusters every 3 particle spacings, overlapping (radius 4 spacings): the body can bend.
    body.clusters = buildClusters(body.particles, rest, 3.0f * s, 4.0f * s);
    bindSurface(body, primitives::subdivided(shape, 1.5f * s), rest);
    softBodies_.push_back(std::move(body));
    return int(softBodies_.size()) - 1;
}

int ParticleSystem::addCloth(const Vector3& origin, const Vector3& u, const Vector3& v, const ClothMaterial& material, int pinMask,
                        const Vector3& color) {
    const float s = params.clothSpacing * params.particleRadius;
    const float sheetVolume = sqr(0.5f * params.clothSpacing); // s x s x 2r relative to (2r)^3
    Cloth c;
    c.object = nextObject_++;
    c.color = color;
    c.material = material;
    c.width = std::max(2, int(std::lround(length(u) / s)) + 1);
    c.height = std::max(2, int(std::lround(length(v) / s)) + 1);
    c.spacing = length(u) / float(c.width - 1);
    c.firstParticle = int(x_.size());
    // Every particle stands for the same share of the sheet: its area (for the fire's heat and fuel
    // and the gas drag) and its mass come from the same number.
    c.particleArea = length(cross(u, v)) / float(c.width * c.height);
    const float invMass = 1.0f / (material.areaDensity * c.particleArea);
    std::vector<Vector3> rest;
    for (int y = 0; y < c.height; ++y)
        for (int x = 0; x < c.width; ++x) {
            Vector3 p = origin + u * (float(x) / (c.width - 1)) + v * (float(y) / (c.height - 1));
            bool corner00 = x == 0 && y == 0, corner10 = x == c.width - 1 && y == 0;
            bool corner01 = x == 0 && y == c.height - 1, corner11 = x == c.width - 1 && y == c.height - 1;
            bool pinned = (corner00 && (pinMask & 1)) || (corner10 && (pinMask & 2)) || (corner01 && (pinMask & 4)) ||
                          (corner11 && (pinMask & 8)) || (y == 0 && (pinMask & 16)) || (y == c.height - 1 && (pinMask & 32)) ||
                          (x == 0 && (pinMask & 64)) || (x == c.width - 1 && (pinMask & 128));
            rest.push_back(p);
            addParticle(p, Vector3(0.0f), ParticlePhase::Cloth, c.object, pinned ? 0.0f : invMass, sheetVolume);
        }
    buildClothConstraints(c, rest);
    buildTethers(c, invMass_);
    cloths_.push_back(std::move(c));
    return int(cloths_.size()) - 1;
}

bool ParticleSystem::grab(const Vector3& point) {
    releaseGrab();
    const float s = spacing();
    int picked = -1;
    float best = sqr(3.0f * s);
    for (int i = 0; i < int(x_.size()); ++i) {
        float d = length2(x_[i] - point);
        if (d < best) { best = d; picked = i; }
    }
    if (picked < 0) return false;
    // Cloth: a pinch of fabric (~3 cm, like fingers) - after the first threads give, the load keeps
    // flowing through the pinch and the crack grows; soft body / liquid: a small ball.
    const float radius = phase_[picked] == uint8_t(ParticlePhase::Cloth) ? 2.0f * s
                         : phase_[picked] == uint8_t(ParticlePhase::Soft) ? 1.5f * s
                                                                         : 2.5f * s;
    grab_.particles.push_back(picked);
    for (int i = 0; i < int(x_.size()); ++i)
        if (i != picked && object_[i] == object_[picked] && length2(x_[i] - x_[picked]) <= radius * radius)
            grab_.particles.push_back(i);
    for (int i : grab_.particles) {
        grab_.offsets.push_back(x_[i] - x_[picked]);
        grab_.savedInvMass.push_back(invMass_[i]);
        invMass_[i] = 0; // kinematic while held
    }
    grab_.target = x_[picked];
    return true;
}

void ParticleSystem::releaseGrab() {
    for (size_t k = 0; k < grab_.particles.size(); ++k)
        if (grab_.particles[k] < int(invMass_.size())) invMass_[grab_.particles[k]] = grab_.savedInvMass[k];
    grab_ = ParticleGrab();
}

void ParticleSystem::burnCloths(float dt, const std::function<GasHeat(const Vector3&)>& gas, float ambientTemperature,
                                std::vector<FireOutput>& out) {
    for (Cloth& c : cloths_) {
        const ClothMaterial& m = c.material;
        if (!m.flammable) continue;
        const int count = c.width * c.height;
        std::vector<float> gasT(count), irradiance(count), heat(count, 0.0f), fuel(count, 0.0f), unburntBefore = c.unburnt;
        parallelFor(count, [&](int k) {
            const GasHeat g = gas(x_[c.firstParticle + k]);
            gasT[k] = g.temperature;
            irradiance[k] = g.irradiance;
        }, 64);
        burnCloth(c, gasT, irradiance, ambientTemperature, dt, heat, fuel);
        for (int k = 0; k < count; ++k) {
            const int i = c.firstParticle + k;
            // Burnt fabric is lighter: mass = fresh mass * (char + (1 - char) * unburnt).
            if (c.unburnt[k] != unburntBefore[k]) {
                const float before = m.charMassFraction + (1.0f - m.charMassFraction) * unburntBefore[k];
                const float now = m.charMassFraction + (1.0f - m.charMassFraction) * c.unburnt[k];
                invMass_[i] *= before / now;
                // A particle held by the mouse has inverse mass 0 until release: its real one waits
                // in the grab record and must get lighter too.
                for (size_t g = 0; g < grab_.particles.size(); ++g)
                    if (grab_.particles[g] == i) grab_.savedInvMass[g] *= before / now;
            }
            if (heat[k] != 0 || fuel[k] != 0) out.push_back({x_[i], heat[k], fuel[k]});
        }
    }
}

void ParticleSystem::stepClothsInSmallSteps(float dt) {
    const int m = std::max(1, params.clothSubsteps);
    const float h = dt / float(m);
    const Vector3 g = params.gravity;
    for (Cloth& c : cloths_) {
        const int first = c.firstParticle, count = c.width * c.height;
        // Start of the substep: positions x, velocities before this substep's gravity; the
        // pinned / grabbed particles move linearly to where the substep puts them (p_ now).
        std::vector<Vector3> q(count), u(count), end(count);
        for (int k = 0; k < count; ++k) {
            const int i = first + k;
            q[k] = x_[i];
            u[k] = invMass_[i] > 0 ? v_[i] - g * dt : Vector3(0.0f);
            end[k] = p_[i];
        }
        for (int s = 0; s < m; ++s) {
            if (s == m / 2) updateTethers(c); // torn during the first half: re-measure mid-step
            for (int k = 0; k < count; ++k) {
                const int i = first + k;
                if (invMass_[i] == 0) {
                    p_[i] = q[k] + (end[k] - q[k]) * (1.0f / float(m - s)); // linear to the end point
                    continue;
                }
                u[k] += g * h;
                p_[i] = q[k] + u[k] * h;
            }
            for (DistanceConstraint& dc : c.constraints) dc.lambda = 0;
            solveCloth(c, p_, invMass_, h);
            if (tearCloth(c, p_, h) > 0 && s + 1 == m) updateTethers(c);
            for (int k = 0; k < count; ++k) {
                const int i = first + k;
                if (invMass_[i] > 0) {
                    collide(i, p_[i], q[k], false, h); // walls, obstacle, rigid bodies (friction on this small step)
                    u[k] = (p_[i] - q[k]) / h;
                }
                q[k] = p_[i];
            }
        }
    }
}

void ParticleSystem::findParticleContacts() {
    const int n = int(p_.size());
    const float d0 = spacing(), reach2 = sqr(1.5f * d0);
    // Per particle in parallel, then concatenated in particle order (deterministic).
    std::vector<std::vector<ParticleContact>> local(n);
    parallelFor(n, [&](int i) {
        const int* nb = &nbr_[size_t(i) * kMaxNeighbors];
        for (int m = 0; m < nbrCount_[i]; ++m) {
            const int j = nb[m];
            if (j <= i) continue;                   // each pair once
            if (isFluid(i) && isFluid(j)) continue; // fluid-fluid: the density constraint
            if (invMass_[i] + invMass_[j] == 0) continue;
            if (object_[i] >= 0 && object_[i] == object_[j]) {
                if (phase_[i] == uint8_t(ParticlePhase::Soft)) continue;         // held by shape matching
                if (length2(rest_[i] - rest_[j]) < 2.25f * d0 * d0) continue;     // cloth neighbours at rest
            }
            if (length2(p_[i] - p_[j]) >= reach2) continue; // cannot touch within this substep
            // The normal comes from the positions at the start of the step and stays fixed: the
            // pair may not swap sides within a step, even if other constraints (a cloth snapping
            // back, a heavy body pushing) move one particle past the other.
            Vector3 nrm = x_[i] - x_[j];
            if (length2(nrm) < 1e-12f) nrm = p_[i] - p_[j];
            if (length2(nrm) < 1e-12f) continue;
            local[i].push_back({i, j, normalize(nrm)});
        }
    });
    contacts_.clear();
    for (auto& l : local) contacts_.insert(contacts_.end(), l.begin(), l.end());
}

void ParticleSystem::solveParticleContacts() {
    // Gauss-Seidel: each contact is resolved at once, so a stack of contacts (a body resting on
    // cloth resting on ...) passes its correction through within one sweep. The pair's inverse
    // masses split the correction (momentum conserving).
    // Limitation of all position-based solvers (FleX included): with a large mass ratio between
    // touching particles (a heavy body on a very light cloth, beyond ~1:10) the light side takes
    // almost the whole correction and the support converges too slowly - use realistic materials
    // (canvas ~1-2 kg/m^2 under foam-like bodies) or more solid iterations.
    const float d0 = spacing(), d02 = d0 * d0;
    for (const ParticleContact& c : contacts_) {
        const Vector3 r = p_[c.i] - p_[c.j];
        const float along = dot(r, c.normal);               // separation along the normal
        if (along >= d0) continue;
        if (length2(r - c.normal * along) >= d02) continue;  // side by side, not touching
        const float wi = invMass_[c.i], wj = invMass_[c.j];
        const Vector3 corr = c.normal * ((d0 - along) / (wi + wj));
        p_[c.i] += corr * wi;
        p_[c.j] -= corr * wj;
    }
}

void ParticleSystem::emitParticles(float dt) {
    if (!emitter.enabled) return;
    const float s = spacing();
    emitter.accumulated += emitter.speed * dt;
    Vector3 d = normalize(emitter.direction);
    Vector3 u = normalize(std::fabs(d.x) > 0.9f ? cross(d, Vector3(0, 1, 0)) : cross(d, Vector3(1, 0, 0)));
    Vector3 w = cross(d, u);
    int n = int(emitter.radius / s);
    while (emitter.accumulated >= s) {
        emitter.accumulated -= s;
        Vector3 c = emitter.position + d * emitter.accumulated;
        for (int a = -n; a <= n; ++a)
            for (int b = -n; b <= n; ++b) {
                if (float(a * a + b * b) * s * s > emitter.radius * emitter.radius) continue;
                if (int(x_.size()) >= params.maxParticles) return;
                Vector3 p = c + u * (a * s) + w * (b * s);
                if (!domain_.contains(p)) continue;
                addParticle(p, d * emitter.speed, ParticlePhase::Fluid, -1, 1.0f / mass_);
            }
    }
}

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
    const float tk = params.tensileK;
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

void ParticleSystem::solveBodyContacts(float dt) {
    if (!rigid_ || bodyShift_.empty()) return;
    const auto& bodies = rigid_->bodies();
    // Contacts of each body, solved one after another (Gauss-Seidel): every contact sees the body
    // already moved by the ones before it.
    std::vector<std::vector<int>> contacts(bodies.size());
    for (int i = 0; i < int(contactBody_.size()); ++i)
        if (contactBody_[i] >= 0) {
            contacts[contactBody_[i]].push_back(i);
            contactBody_[i] = -1;
        }
    for (size_t b = 0; b < bodies.size(); ++b) {
        if (contacts[b].empty()) continue;
        const RigidBody& body = bodies[b];
        const Vector3 shift0 = bodyShift_[b], turn0 = bodyTurn_[b]; // the pose the contacts were found at
        for (int i : contacts[b]) {
            const Vector3 n = contactNormal_[i], c = contactPoint_[i];
            const Vector3 rc = c - (body.pos + bodyShift_[b]);
            // Penetration now: what it was, minus how far the body has moved away from it since.
            const Vector3 moved = (bodyShift_[b] - shift0) + cross(bodyTurn_[b] - turn0, rc);
            const float depth = contactDepth_[i] + dot(moved, n);
            if (depth <= 0) continue;
            // Generalized inverse masses (Mueller et al. 2020, eq. 2-3): the particle, the body at c.
            const float wp = invMass_[i];
            const Vector3 rn = cross(rc, n);
            const float wb = body.invMass + dot(rn, body.applyInvInertiaWorld(rn));
            const float lambda = depth / (wp + wb); // [kg m]
            p_[i] += n * (lambda * wp);
            bodyShift_[b] -= n * (lambda * body.invMass);
            bodyTurn_[b] -= body.applyInvInertiaWorld(rn) * lambda;
            // Friction: the tangential slip of the particle against the body, shared the same way.
            if (params.wallFriction > 0 && wp > 0) {
                const Vector3 bodyMove = body.velocityAt(c) * dt + (bodyShift_[b] + cross(bodyTurn_[b], rc));
                const Vector3 slip = (p_[i] - x_[i]) - bodyMove;
                const Vector3 t = (slip - n * dot(slip, n)) * params.wallFriction;
                const Vector3 lt = t / (wp + wb);
                p_[i] -= lt * wp;
                bodyShift_[b] += lt * body.invMass;
                bodyTurn_[b] += body.applyInvInertiaWorld(cross(rc, lt));
            }
        }
    }
}

void ParticleSystem::collide(int i, Vector3& p, const Vector3& start, bool record, float dt) {
    const float r = params.particleRadius;
    // Domain walls.
    p = vmax(domain_.lo + Vector3(r), vmin(p, domain_.hi - Vector3(r)));

    // Static obstacle mesh.
    if (mesh_ && !mesh_->empty()) {
        AABB mb = mesh_->bounds();
        mb.lo -= Vector3(r);
        mb.hi += Vector3(r);
        if (mb.contains(p)) {
            ClosestHit hit;
            if (mesh_->closestPoint(p, 3 * r, hit) && hit.signedDistance < r) {
                Vector3 n = hit.normal;
                p += n * (r - hit.signedDistance);
                Vector3 dx = p - start;
                p -= (dx - n * dot(dx, n)) * params.wallFriction;
            }
        }
    }

    // Rigid bodies. Fixed ones (and every body in the cloth's small steps, record = false): the
    // particle is pushed out. Movable ones in the main passes: the contact is only recorded here
    // and solved together with the body in solveBodyContacts. A body is met where this substep
    // has moved it: its start pose shifted by bodyShift_ and turned by bodyTurn_ (small angle).
    if (rigid_) {
        auto& bodies = rigid_->bodies();
        const bool coupled = record && bodyShift_.size() == bodies.size();
        for (int b = 0; b < int(bodies.size()); ++b) {
            const RigidBody& body = bodies[b];
            const Vector3 shift = coupled ? bodyShift_[b] : Vector3(0.0f), turn = coupled ? bodyTurn_[b] : Vector3(0.0f);
            const Vector3 centre = body.pos + shift;
            if (length2(p - centre) > sqr(body.boundingRadius() + r)) continue;
            Vector3 n;
            const float sd = body.signedDistance(body.pos + (p - centre) - cross(turn, p - centre), n);
            if (sd >= r) continue;
            n = normalize(n + cross(turn, n));
            if (coupled && body.invMass > 0) {
                contactBody_[i] = b;
                contactNormal_[i] = n;
                contactDepth_[i] = r - sd;
                contactPoint_[i] = p - n * r;
                continue;
            }
            p += n * (r - sd);
            const Vector3 rel = (p - start) - body.velocityAt(p) * dt;
            p -= (rel - n * dot(rel, n)) * params.wallFriction;
        }
    }
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

void ParticleSystem::step(float dt) {
    emitParticles(dt);
    const int n = int(x_.size());
    if (n == 0) return;
    p_.resize(n); dp_.resize(n); rho_.resize(n); lambda_.resize(n); omega_.resize(n); vtmp_.resize(n);
    contactBody_.assign(n, -1);
    contactNormal_.resize(n);
    contactPoint_.resize(n);
    contactDepth_.resize(n);
    const size_t nb = rigid_ ? rigid_->bodies().size() : 0;
    bodyShift_.assign(nb, Vector3(0.0f));
    bodyTurn_.assign(nb, Vector3(0.0f));

    const Vector3 g = params.gravity;
    parallelFor(n, [&](int i) {
        if (invMass_[i] == 0) { // pinned
            v_[i] = Vector3(0.0f);
            p_[i] = x_[i];
            return;
        }
        v_[i] += g * dt;
        Vector3 p = x_[i] + v_[i] * dt;
        collide(i, p, x_[i], true, dt);
        p_[i] = p;
    });
    solveBodyContacts(dt);
    for (size_t k = 0; k < grab_.particles.size(); ++k) p_[grab_.particles[k]] = grab_.target + grab_.offsets[k];
    stepClothsInSmallSteps(dt);
    const bool solids = fluidCount_ < size_t(n);
    for (Cloth& c : cloths_)
        for (DistanceConstraint& dc : c.constraints) dc.lambda = 0;

    buildGrid(p_);
    findNeighbors();

    for (int it = 0; it < params.solverIterations; ++it) {
        computeLambda();
        computeDeltaP();
        parallelFor(n, [&](int i) {
            if (invMass_[i] == 0) return;
            Vector3 p = p_[i] + dp_[i];
            collide(i, p, x_[i], true, dt);
            p_[i] = p;
        });
        solveBodyContacts(dt);
        if (!solids) continue;
        if (it == 0) findParticleContacts();
        for (int pass = 0; pass < std::max(1, params.solidIterations); ++pass) {
            solveParticleContacts();
            // Contacts push cloth particles around (a body resting on a sheet): the cloth is
            // re-satisfied after every contact pass so the two converge together.
            for (Cloth& c : cloths_) solveCloth(c, p_, invMass_, dt);
            solveShapeMatching(softBodies_, p_, invMass_);
            parallelFor(n, [&](int i) {
                if (invMass_[i] == 0 || isFluid(i)) return;
                collide(i, p_[i], x_[i], true, dt);
            });
            solveBodyContacts(dt);
        }
    }

    for (Cloth& c : cloths_) {
        tearCloth(c, p_, dt);
        updateTethers(c);
    }

    const float vmaxAllowed = 0.5f * h_ / dt;
    float vmaxSeen = parallelMax<float>(n, 0.0f, [&](int b, int e) {
        float m = 0;
        for (int i = b; i < e; ++i) {
            Vector3 v = (p_[i] - x_[i]) / dt;
            float l = length(v);
            if (l > vmaxAllowed) v *= vmaxAllowed / l;
            v_[i] = v;
            m = std::max(m, std::min(l, vmaxAllowed));
        }
        return m;
    });
    maxSpeed_ = vmaxSeen;

    applyViscosityAndVorticity(dt);
    x_.swap(p_);

    // The bodies' motion from the particle contacts becomes their velocity change (the rigid
    // solver moves them in its next step).
    for (size_t b = 0; b < bodyShift_.size(); ++b)
        if (length2(bodyShift_[b]) + length2(bodyTurn_[b]) > 0)
            rigid_->applyVelocityChange(int(b), bodyShift_[b] / dt, bodyTurn_[b] / dt);
}

} // namespace rf
