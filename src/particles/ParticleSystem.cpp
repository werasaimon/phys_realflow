// ParticleSystem: the particles, what is made of them (blocks of liquid, soft bodies, cloth),
// the emitter and the step. The liquid's density solver is in DensitySolver.cpp, the contacts in
// ParticleContacts.cpp, cloth and soft bodies in Cloth.cpp and SoftBody.cpp.
#include "particles/ParticleSystem.h"

#include "core/Parallel.h"
#include "core/Probe.h"

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
    // Clusters every 1.5 particle spacings, each 2 spacings in radius (as FleX): a cluster spans
    // ~4 particles, so a body a few particles thick bends and squashes between its clusters. Bigger
    // clusters (3 / 4 spacings) covered a small body whole and made it rigid.
    // A body must be at least 3 particles across: a cluster of a flat sheet of particles has no
    // definite rotation, and its skin flies apart.
    body.clusterRadius = 2.0f * s;
    body.clusters = buildClusters(body.particles, rest, 1.5f * s, body.clusterRadius);
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

int ParticleSystem::pinParticles(const std::function<bool(const Vector3&)>& region) {
    int pinned = 0;
    for (size_t i = 0; i < x_.size(); ++i) {
        if (isFluid(int(i)) || invMass_[i] == 0 || !region(x_[i])) continue;
        invMass_[i] = 0;
        ++pinned;
    }
    return pinned;
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
    prepareBodyQuery(false); // the small steps meet the bodies where the substep started
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

// One step of position based fluids (Macklin & Müller 2013) with the soft bodies and cloth in the
// same loop: predict where every particle goes, then a few iterations that push the particles
// apart until the density is right and the solids keep their shape, then the velocities are
// read off the corrected positions.
void ParticleSystem::step(float dt) {
    emitParticles(dt);
    const int n = int(x_.size());
    if (n == 0) return;
    beginStep(n);
    predictPositions(dt);
    {
        Probe::Timer timer("particles/cloth ms");
        stepClothsInSmallSteps(dt);
    }
    const bool solids = fluidCount_ < size_t(n);
    for (Cloth& c : cloths_)
        for (DistanceConstraint& dc : c.constraints) dc.lambda = 0;
    {
        Probe::Timer timer("particles/neighbors ms");
        buildGrid(p_);
        findNeighbors();
    }
    for (int it = 0; it < params.solverIterations; ++it) solveIteration(it, solids, dt);
    finishStep(dt);
}

// Work arrays for this step's particle count, and the contact records cleared.
void ParticleSystem::beginStep(int n) {
    p_.resize(n); dp_.resize(n); rho_.resize(n); lambda_.resize(n); omega_.resize(n); vtmp_.resize(n);
    contactBody_.assign(n, -1);
    contactNormal_.resize(n);
    contactPoint_.resize(n);
    contactDepth_.resize(n);
    const size_t nb = rigid_ ? rigid_->bodies().size() : 0;
    bodyShift_.assign(nb, Vector3(0.0f));
    bodyTurn_.assign(nb, Vector3(0.0f));
}

// Explicit prediction: gravity into the velocity, the position one step ahead, pushed out of the
// walls and bodies it would enter (a pinned particle stays). The bodies met on the way get their
// contacts solved, and the grabbed particles follow the mouse.
void ParticleSystem::predictPositions(float dt) {
    const int n = int(x_.size());
    const Vector3 g = params.gravity;
    prepareBodyQuery(true);
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
    Probe::add("rigid/tree queries", n);
    {
        int touching = 0; // particles that met a movable body in the prediction (for the probe)
        for (int b : contactBody_) touching += b >= 0;
        Probe::set("particles/body contacts", touching);
    }
    solveBodyContacts(dt);
    for (size_t k = 0; k < grab_.particles.size(); ++k) p_[grab_.particles[k]] = grab_.target + grab_.offsets[k];
}

// One solver iteration: the density constraint of every fluid particle (lambda, then the
// position correction), the walls and bodies again, and - with solids present - the particle
// contacts, the cloth constraints and the shape matching, each followed by the bodies.
void ParticleSystem::solveIteration(int it, bool solids, float dt) {
    const int n = int(x_.size());
    {
        Probe::Timer timer("particles/density ms"); // summed over the iterations
        computeLambda();
        computeDeltaP();
    }
    prepareBodyQuery(true);
    parallelFor(n, [&](int i) {
        if (invMass_[i] == 0) return;
        Vector3 p = p_[i] + dp_[i];
        collide(i, p, x_[i], true, dt);
        p_[i] = p;
    });
    Probe::add("rigid/tree queries", n);
    solveBodyContacts(dt);
    if (!solids) return;
    Probe::Timer contactTimer("particles/contacts ms"); // the solid passes of this iteration
    if (it == 0) findParticleContacts();
    for (int pass = 0; pass < std::max(1, params.solidIterations); ++pass) {
        solveParticleContacts();
        // Contacts push cloth particles around (a body resting on a sheet): the cloth is
        // re-satisfied after every contact pass so the two converge together.
        for (Cloth& c : cloths_) solveCloth(c, p_, invMass_, dt);
        solveShapeMatching(softBodies_, p_, invMass_, params.solverIterations * std::max(1, params.solidIterations));
        prepareBodyQuery(true);
        parallelFor(n, [&](int i) {
            if (invMass_[i] == 0 || isFluid(i)) return;
            collide(i, p_[i], x_[i], true, dt);
        });
        Probe::add("rigid/tree queries", n - int(fluidCount_));
        solveBodyContacts(dt);
    }
}

// The step's end: cloth tears where it was over-stretched, the velocity is the position change
// over the step (capped at half a kernel radius per step, the limit of the neighbour search),
// viscosity and vorticity confinement, and the bodies get the velocity change the particles
// gave them (the rigid solver moves them in its next step).
void ParticleSystem::finishStep(float dt) {
    const int n = int(x_.size());
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
    Probe::set("particles/count", n);
    Probe::set("particles/fluid", double(fluidCount_));
    Probe::set("particles/max speed", maxSpeed_);
    applyViscosityAndVorticity(dt);
    x_.swap(p_);
    for (size_t b = 0; b < bodyShift_.size(); ++b)
        if (length2(bodyShift_[b]) + length2(bodyTurn_[b]) > 0)
            rigid_->applyVelocityChange(int(b), bodyShift_[b] / dt, bodyTurn_[b] / dt);
}

} // namespace rf
