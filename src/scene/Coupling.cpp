// How the solvers of Simulation act on each other within a frame: bodies and particles
// interleaved in time, the gas pushing bodies, cloth and liquid and they in turn being its
// moving walls, and the loads on the obstacle.
#include "scene/Simulation.h"

#include "core/Parallel.h"
#include "core/Probe.h"

namespace rf {

void Simulation::updateSurfaceLoads() {
    if (!grid.hasObstacle() || !obstacleMesh_ || obstacleMesh_->empty()) {
        surfaceLoads_ = SurfaceLoads();
        return;
    }
    // Moment about the body's reference point (for the wing: ~ the aerodynamic centre), per
    // unit of the body size (diameter / chord).
    surfaceLoads_ = computeSurfaceLoads(grid, *obstacleMesh_, obstacle.position, grid.referenceArea(),
                                        std::max(obstacle.size, 1e-3f));
}

std::vector<MovingSolid> Simulation::movingSolids() const {
    std::vector<MovingSolid> out;
    out.reserve(rigid.bodies().size());
    for (const RigidBody& b : rigid.bodies()) {
        MovingSolid m;
        // A destroyed body keeps its slot (the gas forces are indexed like the bodies): an empty
        // solid that marks no cell.
        if (!b.alive) {
            m.inside = [](const Vector3&) { return false; };
            m.resting = true;
            out.push_back(std::move(m));
            continue;
        }
        m.bounds = b.worldBounds();
        m.position = b.pos;
        m.velocity = b.vel;
        m.angularVelocity = b.angVel;
        m.length = 2.0f * b.boundingRadius();
        std::shared_ptr<const ConvexShape> shape = b.shape;
        const Matrix3x3 Rt = b.rotation().transposed();
        const Vector3 p = b.pos;
        m.inside = [shape, Rt, p](const Vector3& x) { return shape->contains(Rt * (x - p)); };
        m.resting = b.sleeping; // has not moved: the gas reuses its cells
        out.push_back(std::move(m));
    }
    // Soft bodies (after the rigid ones): cells within one particle spacing of a particle are
    // inside; the body moves with the mean velocity of its particles.
    const float s = particles.spacing();
    for (const SoftBody& sb : particles.softBodies()) {
        MovingSolid m;
        std::vector<Vector3> pts;
        Vector3 c(0.0f), v(0.0f);
        for (int i : sb.particles) {
            pts.push_back(particles.positions()[i]);
            m.bounds.expand(particles.positions()[i]);
            c += particles.positions()[i];
            v += particles.velocities()[i];
        }
        m.bounds.lo -= Vector3(s);
        m.bounds.hi += Vector3(s);
        m.position = c / float(pts.size());
        m.velocity = v / float(pts.size());
        m.length = maxComp(m.bounds.extent());
        const float r2 = s * s;
        m.inside = [pts = std::move(pts), r2](const Vector3& x) {
            for (const Vector3& p : pts)
                if (length2(p - x) < r2) return true;
            return false;
        };
        out.push_back(std::move(m));
    }
    return out;
}

void Simulation::applyGasDragOnCloth(float dt) {
    // Aerodynamic force of the gas on every cloth particle (its patch of the sheet): pressure drag
    // of a flat plate on the part of the relative wind normal to the sheet (Cd 1.2), skin friction
    // along it (Cf 0.02). Taken implicitly (the particle cannot overtake the wind in one step); the
    // reaction goes into the gas at the same point, so momentum is conserved.
    const float rho = grid.params.fluidDensity;
    const auto& x = particles.positions();
    const auto& v = particles.velocities();
    const auto& w = particles.invMasses();
    for (const Cloth& c : particles.cloths()) {
        const float area = c.particleArea; // the share of the sheet one particle stands for
        for (int y = 0; y < c.height; ++y)
            for (int xg = 0; xg < c.width; ++xg) {
                const int i = c.particle(xg, y);
                if (w[i] == 0) continue;
                const Vector3 du = x[c.particle(std::min(xg + 1, c.width - 1), y)] - x[c.particle(std::max(xg - 1, 0), y)];
                const Vector3 dvv = x[c.particle(xg, std::min(y + 1, c.height - 1))] - x[c.particle(xg, std::max(y - 1, 0))];
                const Vector3 n = normalize(cross(du, dvv));
                const Vector3 rel = grid.velocityAt(x[i]) - v[i];
                const float wn = dot(rel, n);
                const Vector3 relT = rel - n * wn;
                const float kn = 0.5f * rho * 1.2f * area * std::fabs(wn) * dt * w[i];
                const float kt = 0.5f * rho * 0.02f * area * length(relT) * dt * w[i];
                const Vector3 dv = n * (wn * kn / (1.0f + kn)) + relT * (kt / (1.0f + kt));
                particles.addVelocity(i, dv);
                grid.addImpulse(x[i], -dv / w[i]);
            }
    }
}

void Simulation::applyGasDragOnLiquid(float dt) {
    // Drag of the air on the liquid particles that touch air (the surface, spray, drops): a sphere
    // of the particle's radius (Cd 0.47) in the relative wind, implicit; the reaction goes into the
    // air. The wind blows spray away and drags the surface; flying drops are slowed down.
    if (particles.fluidCount() == 0) return;
    const float rho = grid.params.fluidDensity, r = particles.params.particleRadius;
    const float area = kPi * r * r;
    const auto& x = particles.positions();
    const auto& v = particles.velocities();
    const auto& w = particles.invMasses();
    const auto& phase = particles.phases();
    for (size_t i = 0; i < x.size(); ++i) {
        if (phase[i] != uint8_t(ParticlePhase::Fluid) || w[i] == 0 || !grid.touchesGas(x[i])) continue;
        const Vector3 rel = grid.fluidVelocityAt(x[i]) - v[i];
        const float k = 0.5f * rho * 0.47f * area * length(rel) * dt * w[i];
        const Vector3 dv = rel * (k / (1.0f + k));
        particles.addVelocity(int(i), dv);
        grid.addImpulse(x[i], -dv / w[i]);
    }
}

void Simulation::passLiquidToGas() {
    if (particles.fluidCount() == 0) {
        grid.setLiquid({}, {});
        return;
    }
    std::vector<Vector3> pos, vel;
    pos.reserve(particles.fluidCount());
    vel.reserve(particles.fluidCount());
    for (size_t i = 0; i < particles.size(); ++i)
        if (particles.phases()[i] == uint8_t(ParticlePhase::Fluid)) {
            pos.push_back(particles.positions()[i]);
            vel.push_back(particles.velocities()[i]);
        }
    grid.setLiquid(std::move(pos), std::move(vel));
}

void Simulation::applyGasOnSoftBodies() {
    // Soft bodies: the gas pressure impulses of the last gas steps (plus the buoyancy of the
    // displaced gas) as a velocity change of the whole body.
    const auto& bodies = particles.softBodies();
    for (size_t b = 0; b < bodies.size() && b < softGasImpulse_.size(); ++b) {
        float mass = 0;
        for (int i : bodies[b].particles) mass += particles.invMasses()[i] > 0 ? 1.0f / particles.invMasses()[i] : 0.0f;
        if (mass <= 0) continue;
        const float volume = float(bodies[b].particles.size()) * std::pow(particles.spacing(), 3.0f);
        const Vector3 J = softGasImpulse_[b] - gravity() * (grid.params.fluidDensity * volume * frameDt);
        for (int i : bodies[b].particles)
            if (particles.invMasses()[i] > 0) particles.addVelocity(i, J / mass);
    }
}

// Step 1 of stepGasWithBodies: the pressure impulses the gas gave every body over the last
// frame's gas steps, plus the buoyancy of the displaced gas (the Boussinesq grid carries no
// hydrostatic pressure), as a wrench; the soft bodies take theirs as a velocity change.
void Simulation::pushBodiesByGas() {
    Probe::Timer timer("scene/coupling ms"); // what the solvers pass each other
    const float rho = grid.params.fluidDensity;
    for (int i = 0; i < int(rigid.bodies().size()); ++i) {
        RigidBody& b = rigid.bodies()[i];
        if (b.invMass == 0) continue;
        Vector3 J = i < int(gasImpulse_.size()) ? gasImpulse_[i] : Vector3(0.0f);
        Vector3 L = i < int(gasAngularImpulse_.size()) ? gasAngularImpulse_[i] : Vector3(0.0f);
        J -= gravity() * (rho * b.shape->volume() * frameDt); // Archimedes in the gas
        // A resting body is not woken by pressure noise far below its sleep threshold.
        if (b.sleeping && length(J) * b.invMass < 0.5f * rigid.params.sleepLinear) continue;
        rigid.applyExternalWrench(i, J, L);
    }
    applyGasOnSoftBodies();
}

void Simulation::stepBodiesAndParticles(bool gasDrag) {
    // Every rigid step has the length the rigid solver asks for (frameDt / its substeps), with or
    // without particles; before each particle step come the rigid steps that end within it
    // (10 rigid and 3 particle steps: 3, 3, 4). With no particles the particle steps cost nothing
    // (only the emitter runs).
    const int nr = std::max(1, rigid.params.substeps), np = std::max(1, particles.params.substeps);
    const float hr = frameDt / float(nr), hp = frameDt / float(np);
    int r = 0;
    for (int p = 0; p < np; ++p) {
        const int rEnd = (p + 1) * nr / np; // rigid steps done by the end of this particle step
        for (; r < rEnd; ++r) rigid.step(hr);
        if (gasDrag) {
            Probe::Timer timer("scene/coupling ms");
            applyGasDragOnCloth(hp);
            applyGasDragOnLiquid(hp);
        }
        particles.step(hp);
    }
}

void Simulation::stepGasWithBodies() {
    // Weak (staggered) two-way coupling, one exchange per frame:
    //  1) gas -> bodies: pressure impulses integrated over the previous frame's gas steps, plus the
    //     buoyancy of the displaced gas (the Boussinesq grid carries no hydrostatic pressure);
    //  2) rigid substeps;
    //  3) bodies -> gas: the bodies at their new poses and velocities are the moving boundaries of
    //     the gas steps that cover the same frame time.
    if (rigid.anyHeld() && time_ >= releaseTime) rigid.releaseHeld();
    if (gasPushesBodies) pushBodiesByGas();
    stepBodiesAndParticles(gasPushesBodies); // rigid bodies interleaved with the particles

    // Fire: the cloths take heat from the gas (or burn and give it heat and fuel gas).
    if (grid.combustion.enabled) {
        Probe::Timer timer("particles/heat ms");
        std::vector<FireOutput> fire;
        auto gasHeat = [this](const Vector3& x) { return GasHeat{grid.temperatureAt(x), grid.irradianceAt(x)}; };
        particles.burnCloths(frameDt, gasHeat, grid.combustion.ambientTemperature, fire);
        for (const FireOutput& f : fire) grid.addEmission({f.position, f.fuel, f.heat, 0.0f});
    }

    gasImpulse_.assign(rigid.bodies().size(), Vector3(0.0f));
    gasAngularImpulse_.assign(rigid.bodies().size(), Vector3(0.0f));
    softGasImpulse_.assign(particles.softBodies().size(), Vector3(0.0f));
    float remaining = frameDt;
    for (int it = 0; it < 8 && remaining > 1e-6f; ++it) {
        {
            Probe::Timer timer("scene/coupling ms");
            grid.setMovingSolids(movingSolids());
            passLiquidToGas();
        }
        float dt = grid.step(remaining);
        const auto& F = grid.movingForces();
        const auto& T = grid.movingTorques();
        for (size_t i = 0; i < F.size(); ++i) {
            if (i < gasImpulse_.size()) {
                gasImpulse_[i] += F[i] * dt;
                gasAngularImpulse_[i] += T[i] * dt;
            } else if (i - gasImpulse_.size() < softGasImpulse_.size()) {
                softGasImpulse_[i - gasImpulse_.size()] += F[i] * dt;
            }
        }
        remaining -= dt;
        lastGridDt_ = dt;
    }
    gasForceMax_ = 0;
    for (const Vector3& J : gasImpulse_) gasForceMax_ = std::max(gasForceMax_, length(J) / frameDt);
}

} // namespace rf
