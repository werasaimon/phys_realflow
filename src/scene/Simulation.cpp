// Simulation: the solvers as one scene - construction, reset and the frame step. The scenes
// themselves are in Presets.cpp, how the solvers act on each other in Coupling.cpp, and what the
// viewer gets to draw in Snapshot.cpp.
#include "scene/Simulation.h"

#include <chrono>

namespace rf {

Simulation::Simulation() {
    obstacleMesh_ = std::make_shared<TriMesh>();
    fluidDomain_ = AABB({-1.0f, 0.0f, -0.4f}, {1.0f, 1.2f, 0.4f});
    rigidDomain_ = AABB({-2.0f, 0.0f, -2.0f}, {2.0f, 5.0f, 2.0f});
    loadPreset(Preset::DamBreak);
}

void Simulation::setGravity(const Vector3& g) {
    rigid.params.gravity = g;
    particles.params.gravity = g;
    grid.combustion.gravity = length(g); // the hot gas rises against it (along +y)
}

void Simulation::reset() {
    time_ = 0;
    frame_ = 0;
    surfaceLoads_ = SurfaceLoads(); // the loads of the old geometry (filled again by the tunnel steps)
    switch (mode_) {
    case SimMode::Fluid: setupFluidScene(); break;
    case SimMode::WindTunnel: setupTunnelScene(); break;
    case SimMode::Rigid: setupRigidScene(); break;
    }
}

// ---------------------------------------------------------------------------
// Stepping
// ---------------------------------------------------------------------------
void Simulation::stepFrame() {
    auto t0 = std::chrono::steady_clock::now();
    switch (mode_) {
    case SimMode::Fluid:
    case SimMode::Rigid: // (bodies, and whatever particles were added: the same stepping)
        stepBodiesAndParticles(false);
        time_ += frameDt;
        break;
    case SimMode::WindTunnel:
        if (rigid.bodies().empty() && particles.size() == 0) { // gas alone
            grid.setMovingSolids({});
            lastGridDt_ = grid.step(0.05f);
            // Plasma: the Alfven waves keep the steps short - take a few per frame (up to 1/60 s).
            for (int extra = 0; extra < 3 && grid.magnetic.enabled && grid.time() - time_ < frameDt; ++extra)
                lastGridDt_ = grid.step(frameDt - (grid.time() - time_));
            time_ = grid.time();
            if (preset_ == Preset::Tokamak) controlTokamakPosition();
        } else {
            stepGasWithBodies();
            time_ += frameDt;
        }
        updateSurfaceLoads();
        break;
    }
    ++frame_;
    lastStepMs_ = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

int Simulation::sliceLayer() const {
    int n[3] = {grid.nx(), grid.ny(), grid.nz()};
    int axis = clampv(vis.sliceAxis, 0, 2);
    return clampv(int(vis.slicePosition * n[axis]), 0, n[axis] - 1);
}

} // namespace rf
