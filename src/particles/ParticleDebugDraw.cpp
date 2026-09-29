// The research view of the particle solver, drawn by layer into the Probe after the step (see
// core/Probe.h): the neighbour list of the particle nearest to the probe point (what the density
// constraint of Macklin & Müller 2013 sums over), the density error of the liquid (rho / rho0 - 1:
// what PBF drives to zero), the tetrahedra of the soft bodies coloured by their volume change (what
// the hydrostatic constraint holds) and the tension of every cloth thread against its strength
// (what tears it).
#include "particles/ParticleSystem.h"

#include "core/Format.h"
#include "core/Probe.h"

#include <algorithm>
#include <cmath>

namespace rf {

// Every particle layer, at the end of a step; the frame's last substep replaces the ones before.
void ParticleSystem::drawDebug(float dt) const {
    Probe::clearLayers(Probe::bits(DrawLayer::ParticleNeighbours, DrawLayer::ClothTension));
    drawNeighbours();
    drawDensityError();
    drawSoftTetrahedra();
    drawClothTension(dt);
}

// The particle nearest to Probe::probePoint(), ringed, and a line to each of its neighbours.
void ParticleSystem::drawNeighbours() const {
    if (!Probe::layerOn(DrawLayer::ParticleNeighbours) || x_.empty() || nbrCount_.size() != x_.size()) return;
    const Vector3 at = Probe::probePoint();
    int best = 0;
    for (int i = 1; i < int(x_.size()); ++i)
        if (length2(x_[size_t(i)] - at) < length2(x_[size_t(best)] - at)) best = i;
    const Vector3 c = x_[size_t(best)];
    Probe::point(DrawLayer::ParticleNeighbours, c, Vector3(1.0f), 2.5f * params.particleRadius);
    const int count = nbrCount_[size_t(best)];
    const int* nb = &nbr_[size_t(best) * kMaxNeighbors];
    for (int k = 0; k < count; ++k) Probe::line(DrawLayer::ParticleNeighbours, c, x_[size_t(nb[k])], Vector3(0.3f, 1.0f, 0.6f));
    Probe::label(DrawLayer::ParticleNeighbours, c, format("частица %d: %d соседей, rho %.0f", best, count,
                                                          double(best < int(rho_.size()) ? rho_[size_t(best)] : 0.0f)));
}

// Every liquid particle coloured by its density error: blue 10 % too light, green exact, red 10 %
// too dense.
void ParticleSystem::drawDensityError() const {
    if (!Probe::layerOn(DrawLayer::DensityError) || rho_.size() != x_.size()) return;
    const float rho0 = std::max(params.restDensity, 1e-6f);
    for (size_t i = 0; i < x_.size(); ++i) {
        if (phase_[i] != uint8_t(ParticlePhase::Fluid)) continue;
        const float error = rho_[i] / rho0 - 1.0f;
        Probe::point(DrawLayer::DensityError, x_[i], heatColor(0.5f + 5.0f * error), params.particleRadius);
    }
}

// Every soft body's tetrahedra: their edges coloured by the volume ratio det F - blue squeezed to
// half, green at rest, red blown up to one and a half; a tetrahedron turned inside out is drawn white.
void ParticleSystem::drawSoftTetrahedra() const {
    if (!Probe::layerOn(DrawLayer::SoftTetrahedra)) return;
    static const int edges[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    for (const SoftBody& b : softBodies_)
        for (const SoftTet& t : b.tets) {
            const float J = deformationGradient(t, x_).determinant();
            const Vector3 col = J <= 0 ? Vector3(1.0f) : heatColor(std::clamp(J - 0.5f, 0.0f, 1.0f));
            for (const auto& e : edges) Probe::line(DrawLayer::SoftTetrahedra, x_[size_t(t.v[e[0]])], x_[size_t(t.v[e[1]])], col);
        }
}

// Every cloth thread (warp, weft) coloured by its tension over its strength: blue slack, red about
// to tear. A thread that never tears (strength 0) is scaled by the largest tension of its cloth.
void ParticleSystem::drawClothTension(float dt) const {
    if (!Probe::layerOn(DrawLayer::ClothTension)) return;
    for (const Cloth& c : cloths_) {
        float largest = 1e-6f;
        for (const DistanceConstraint& d : c.constraints)
            if (!d.broken && d.strength <= 0) largest = std::max(largest, threadTension(d, x_, dt));
        for (const DistanceConstraint& d : c.constraints) {
            if (d.broken || (d.kind != DistanceConstraint::Warp && d.kind != DistanceConstraint::Weft)) continue;
            const float tension = std::max(threadTension(d, x_, dt), 0.0f);
            const float t = d.strength > 0 ? tension / d.strength : tension / largest;
            Probe::line(DrawLayer::ClothTension, x_[size_t(d.a)], x_[size_t(d.b)], heatColor(t));
        }
    }
}

} // namespace rf
