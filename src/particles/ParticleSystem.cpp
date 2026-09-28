// ParticleSystem: the particles, what is made of them (blocks of liquid, cloth), the emitter, the
// step, and taking one group out again. The liquid's density solver is in DensitySolver.cpp, the
// contacts in ParticleContacts.cpp, the soft bodies in ParticleSoftBodies.cpp (building, small
// steps), SoftTets.cpp (the Neo-Hookean material) and SoftBody.cpp (the legacy clusters), the
// cloth in Cloth.cpp.
#include "particles/ParticleSystem.h"

#include "core/Parallel.h"
#include "core/Probe.h"

#include <algorithm>

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
    phase_.clear(); object_.clear(); group_.clear(); invMass_.clear(); volume_.clear(); rest_.clear();
    surfaceDepth_.clear(); restSurfaceNormal_.clear(); surfaceNormal_.clear(); friction_.clear();
    contacts_.clear(); intersections_.clear(); unresolved_.clear();
    softBodies_.clear(); cloths_.clear();
    walls_ = domain;
    tetBatchesStale_ = true;
    softSmallSteps_ = 0;
    grab_ = ParticleGrab();
    nextObject_ = 0;
    nextGroup_ = 0;
    emitterGroup_ = -1;
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

int ParticleSystem::addBlock(const AABB& box, const Vector3& vel) {
    const float r = params.particleRadius, s = 2 * r;
    const int group = nextGroup_++;
    AABB b(vmax(box.lo, domain_.lo + Vector3(r)), vmin(box.hi, domain_.hi - Vector3(r)));
    if (!b.valid()) return group;
    int nx = int(b.extent().x / s) + 1, ny = int(b.extent().y / s) + 1, nz = int(b.extent().z / s) + 1;
    uint32_t seed = 12345;
    auto rnd = [&]() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) * (1.0f / 16777216.0f) - 0.5f; };
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                if (int(x_.size()) >= params.maxParticles) return group;
                Vector3 p = b.lo + Vector3(float(i), float(j), float(k)) * s + Vector3(rnd(), rnd(), rnd()) * (0.02f * s);
                if (blockedBySolid(p, r, mesh_, rigid_)) continue;
                addParticle(p, vel, ParticlePhase::Fluid, -1, group, 1.0f / mass_);
            }
    return group;
}

void ParticleSystem::addParticle(const Vector3& x, const Vector3& v, ParticlePhase phase, int object, int group, float invMass,
                                 float volume) {
    x_.push_back(x);
    group_.push_back(group);
    v_.push_back(v);
    phase_.push_back(uint8_t(phase));
    object_.push_back(object);
    invMass_.push_back(invMass);
    volume_.push_back(volume);
    rest_.push_back(x);
    surfaceDepth_.push_back(-1.0f); // no surface; addSoftBody measures its own particles'
    restSurfaceNormal_.push_back(Vector3(0.0f));
    surfaceNormal_.push_back(Vector3(0.0f));
    friction_.push_back(0.0f); // addSoftBody gives its particles the material's
    if (phase == ParticlePhase::Fluid) ++fluidCount_;
}

int ParticleSystem::addCloth(const Vector3& origin, const Vector3& u, const Vector3& v, const ClothMaterial& material, int pinMask,
                        const Vector3& color) {
    const float s = params.clothSpacing * params.particleRadius;
    const float sheetVolume = sqr(0.5f * params.clothSpacing); // s x s x 2r relative to (2r)^3
    Cloth c;
    c.object = nextObject_++;
    c.group = nextGroup_++;
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
            addParticle(p, Vector3(0.0f), ParticlePhase::Cloth, c.object, c.group, pinned ? 0.0f : invMass, sheetVolume);
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

// How many small steps a cloth's step dt is cut into: params.clothSubsteps as a rule ("small
// steps", Macklin et al. 2019), more while its threads are pulled so hard that such a step is
// longer than their tension allows (the string limit dt <= sqrt(m l / T), stableClothStep: a hand
// yanking the cloth, the free edge of a swinging sheet snapping) - at most kMaxClothSmallSteps,
// by then the threads are near their strength and about to break anyway.
static constexpr int kMaxClothSmallSteps = 64;
int ParticleSystem::clothSmallSteps(const Cloth& c, float dt) const {
    const int usual = std::max(1, params.clothSubsteps);
    const float longestStableStep = stableClothStep(c, invMass_); // infinity while nothing pulls
    const int needed = int(std::ceil(dt / longestStableStep));    // 0 then
    return std::min(std::max(needed, usual), std::max(usual, kMaxClothSmallSteps));
}

// The cloth's own motion, in m small steps of h = dt / m each ("small steps", Macklin et al. 2019).
// Every small step is a recipe of four:
//   1. move: a free particle flies on under gravity, a pinned or grabbed one slides in a straight
//      line towards where this substep puts it;
//   2. solve the cloth once: every thread exactly, shear and bending by one pass, the tethers;
//   3. tear the threads loaded beyond their strength;
//   4. collide with the walls, the obstacle and the bodies, and read the velocity off the move.
// The start points, velocities and end points live in the cloth's work space (Cloth::step*), kept
// from step to step, so the cloth's step allocates nothing.
void ParticleSystem::stepClothsInSmallSteps(float dt) {
    const Vector3 g = params.gravity;
    prepareBodyQuery(false); // the small steps meet the bodies where the substep started
    for (Cloth& c : cloths_) {
        const int m = c.smallSteps = clothSmallSteps(c, dt);
        const float h = dt / float(m);
        const int first = c.firstParticle, count = c.width * c.height;
        std::vector<Vector3>& start = c.stepStart;
        std::vector<Vector3>& velocity = c.stepVelocity;
        std::vector<Vector3>& end = c.stepEnd;
        start.resize(size_t(count));
        velocity.resize(size_t(count));
        end.resize(size_t(count));
        for (int k = 0; k < count; ++k) { // where the substep starts, the velocity before its gravity
            const int i = first + k;
            start[k] = x_[i];
            velocity[k] = invMass_[i] > 0 ? v_[i] - g * dt : Vector3(0.0f);
            end[k] = p_[i];
        }
        for (int s = 0; s < m; ++s) {
            if (s == m / 2) updateTethers(c); // torn during the first half: re-measure mid-step
            for (int k = 0; k < count; ++k) { // 1. move
                const int i = first + k;
                if (invMass_[i] == 0) {
                    p_[i] = start[k] + (end[k] - start[k]) * (1.0f / float(m - s));
                    continue;
                }
                velocity[k] += g * h;
                p_[i] = start[k] + velocity[k] * h;
            }
            for (DistanceConstraint& dc : c.constraints) dc.lambda = 0;
            solveCloth(c, p_, invMass_, h, ThreadSolve::ExactLines); // 2. solve
            if (tearCloth(c, p_, h) > 0 && s + 1 == m) updateTethers(c); // 3. tear
            for (int k = 0; k < count; ++k) { // 4. collide, then the velocity of this small step
                const int i = first + k;
                if (invMass_[i] > 0) {
                    collide(i, p_[i], start[k], false, h);
                    velocity[k] = (p_[i] - start[k]) / h;
                }
                start[k] = p_[i];
            }
        }
    }
}

void ParticleSystem::emitParticles(float dt) {
    if (!emitter.enabled) return;
    Probe::Timer timer("particles/emit ms");
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
                if (emitterGroup_ < 0) emitterGroup_ = nextGroup_++; // everything the nozzle releases
                addParticle(p, d * emitter.speed, ParticlePhase::Fluid, -1, emitterGroup_, 1.0f / mass_);
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
    // The contact passes below solve the cloth again over the whole step dt. They start from the
    // force the cloth carries now, not from zero: in XPBD the force is lambda / dt^2, so the last
    // small step's lambda_h becomes lambda_h (dt / h)^2 = lambda_h m^2. Started from zero, the
    // pass took the elastic stretch that holds the cloth's weight for an error and pulled it out
    // (~90 % of it per pass): the threads then read a third of their real load at the end of the
    // step and jerked between the two states every step.
    for (Cloth& c : cloths_) {
        const float smallSteps = float(c.smallSteps); // dt / h
        for (DistanceConstraint& dc : c.constraints) dc.lambda *= smallSteps * smallSteps;
    }
    {
        Probe::Timer timer("particles/neighbors ms");
        buildGrid(p_);
        findNeighbors();
    }
    if (solids) { // the solid contacts of the step; bodies inside each other pulled apart first (no speed)
        Probe::Timer timer("particles/contacts ms");
        findParticleContacts();
        preStabilizeContacts();
        setMainSolveTargets();
    }
    {
        // After the contacts are known: the small steps solve the ones between elastic bodies too.
        Probe::Timer timer("particles/soft tetrahedra ms");
        stepSoftTetsInSmallSteps(dt);
    }
    for (int it = 0; it < params.solverIterations; ++it) solveIteration(solids, dt);
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
    Probe::Timer timer("particles/predict ms");
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
void ParticleSystem::solveIteration(bool solids, float dt) {
    const int n = int(x_.size());
    {
        Probe::Timer timer("particles/density ms"); // summed over the iterations
        computeLambda();
        computeDeltaP();
    }
    {
        Probe::Timer timer("particles/bodies ms"); // the density's push, kept out of the walls and bodies
        prepareBodyQuery(true);
        parallelFor(n, [&](int i) {
            if (invMass_[i] == 0) return;
            Vector3 p = p_[i] + dp_[i];
            collide(i, p, x_[i], true, dt);
            p_[i] = p;
        });
        Probe::add("rigid/tree queries", n);
        solveBodyContacts(dt);
    }
    if (!solids) return;
    Probe::Timer contactTimer("particles/contacts ms"); // the solid passes of this iteration
    for (int pass = 0; pass < std::max(1, params.solidIterations); ++pass) {
        solveParticleContacts();
        // Contacts push cloth particles around (a body resting on a sheet): the cloth is
        // re-satisfied after every contact pass so the two converge together.
        for (Cloth& c : cloths_) solveCloth(c, p_, invMass_, dt, ThreadSolve::Projection);
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
    if (!cloths_.empty()) {
        Probe::Timer timer("particles/cloth ms");
        for (Cloth& c : cloths_) {
            tearCloth(c, p_, dt);
            updateTethers(c);
        }
    }
    {
        Probe::Timer timer("particles/velocity ms"); // v = dx / dt, viscosity, vorticity, the bodies' share
        const float vmaxAllowed = 0.5f * h_ / dt;
        const bool smallSteps = softSmallSteps_ > 0 && softStart_.size() == size_t(n) && tetParticle_.size() == size_t(n);
        float vmaxSeen = parallelMax<float>(n, 0.0f, [&](int b, int e) {
            float m = 0;
            for (int i = b; i < e; ++i) {
                // A free particle of a tetrahedral body: its last small step's velocity, plus what
                // the main passes moved it since (a contact's push). Its mean velocity over the
                // substep lags the true one by a dt / 2 - restarted from it, a swinging beam felt
                // an extra push of half its spring and swung 30 % slow.
                Vector3 v = smallSteps && tetParticle_[i] && invMass_[i] > 0 ? softVelocity_[i] + (p_[i] - softStart_[i]) / dt
                                                                              : (p_[i] - x_[i]) / dt;
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
    if (Probe::drawEnabled()) {
        Probe::Timer timer("particles/debug draw ms");
        drawDebug(dt); // the research layers (ParticleDebugDraw.cpp)
    }
}

// ---------------------------------------------------------------------------
// Groups: removing one object without touching the rest
// ---------------------------------------------------------------------------
size_t ParticleSystem::groupSize(int group) const {
    return size_t(std::count(group_.begin(), group_.end(), group));
}

void ParticleSystem::removeGroup(int group) {
    if (group < 0 || std::find(group_.begin(), group_.end(), group) == group_.end()) return;
    releaseGrab(); // the grabbed particles may be among the removed ones
    const std::vector<int> newIndex = renumberWithout(group);
    size_t kept = 0;
    for (int k : newIndex) kept += k >= 0 ? 1 : 0;
    renumberSolids(group, newIndex); // uses the old numbering of the cloths' first particles
    compactParticles(newIndex, kept);
    if (group == emitterGroup_) emitterGroup_ = -1; // the nozzle starts a new group next time
}

// Where every particle goes: -1 for the removed group, the others close up in their order, so
// what remains is the same system with the gaps taken out.
std::vector<int> ParticleSystem::renumberWithout(int group) const {
    std::vector<int> newIndex(x_.size(), -1);
    int next = 0;
    for (size_t i = 0; i < x_.size(); ++i)
        if (group_[i] != group) newIndex[i] = next++;
    return newIndex;
}

// Every per-particle array keeps the entries of the particles that stay, in order. The arrays a
// step rebuilds from scratch (neighbours, the grid, contacts) are just emptied.
void ParticleSystem::compactParticles(const std::vector<int>& newIndex, size_t kept) {
    const size_t n = newIndex.size();
    auto compact = [&](auto& values) {
        if (values.size() != n) return;
        for (size_t i = 0; i < n; ++i)
            if (newIndex[i] >= 0) values[size_t(newIndex[i])] = values[i];
        values.resize(kept);
    };
    compact(x_); compact(v_); compact(p_); compact(dp_); compact(omega_); compact(vtmp_);
    compact(rho_); compact(lambda_); compact(phase_); compact(object_); compact(group_);
    compact(invMass_); compact(volume_); compact(rest_);
    compact(surfaceDepth_); compact(restSurfaceNormal_); compact(surfaceNormal_); compact(friction_);
    softStart_.clear(); softVelocity_.clear(); softEnd_.clear(); // rebuilt by the next small steps
    compact(contactBody_); compact(contactNormal_); compact(contactPoint_); compact(contactDepth_);
    nbrCount_.clear(); nbr_.clear(); cellOf_.clear(); sorted_.clear(); contacts_.clear();
    intersections_.clear(); unresolved_.clear(); // they name soft bodies by their index
    fluidCount_ = size_t(std::count(phase_.begin(), phase_.end(), uint8_t(ParticlePhase::Fluid)));
}

// Soft bodies and cloths of the group go; the others follow the new numbering. A cloth's
// particles are one block, so it moves down as a whole by the particles removed before it.
void ParticleSystem::renumberSolids(int group, const std::vector<int>& newIndex) {
    auto ofGroup = [group](const auto& solid) { return solid.group == group; };
    softBodies_.erase(std::remove_if(softBodies_.begin(), softBodies_.end(), ofGroup), softBodies_.end());
    for (SoftBody& b : softBodies_) {
        for (int& i : b.particles) i = newIndex[size_t(i)];
        for (SoftCluster& c : b.clusters)
            for (int& i : c.particles) i = newIndex[size_t(i)];
        for (SoftTet& t : b.tets)
            for (int& i : t.v) i = newIndex[size_t(i)];
    }
    tetBatchesStale_ = true;
    cloths_.erase(std::remove_if(cloths_.begin(), cloths_.end(), ofGroup), cloths_.end());
    for (Cloth& c : cloths_) shiftCloth(c, c.firstParticle - newIndex[size_t(c.firstParticle)]);
}

} // namespace rf
