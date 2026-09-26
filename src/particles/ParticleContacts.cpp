// Contacts of ParticleSystem: solid particles against each other, particles against the walls
// and the static mesh, and the two-way XPBD contacts with rigid bodies (Muller et al. 2020).
#include "particles/ParticleSystem.h"

#include "core/Parallel.h"

namespace rf {

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

} // namespace rf
