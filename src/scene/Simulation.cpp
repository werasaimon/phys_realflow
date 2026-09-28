// Simulation: the solvers as one scene - loading a Scene, the boxes it builds in, reset and the
// frame step. How the solvers act on each other is in Coupling.cpp, the obstacle geometry in
// Obstacle.cpp, and what the viewer gets to draw in Snapshot.cpp.
#include "scene/Simulation.h"

#include <chrono>

namespace rf {

const AABB Simulation::kDefaultTank({-1.0f, 0.0f, -0.4f}, {1.0f, 1.2f, 0.4f});
const AABB Simulation::kDefaultArena({-2.0f, 0.0f, -2.0f}, {2.0f, 5.0f, 2.0f});

Simulation::Simulation() {
    obstacleMesh_ = std::make_shared<TriMesh>();
    useLiquidTank(kDefaultTank); // an empty tank until a scene is loaded
}

void Simulation::setGravity(const Vector3& g) {
    rigid.params.gravity = g;
    particles.params.gravity = g;
    grid.combustion.gravity = length(g); // the hot gas rises against it (along +y)
}

// ---------------------------------------------------------------------------
// The scene
// ---------------------------------------------------------------------------
void Simulation::load(std::unique_ptr<Scene> scene) {
    scene_ = std::move(scene);
    obstacle = ObstacleSettings();
    vis = VisSettings();
    particles.params = ParticleParams();
    particles.emitter = ParticleEmitter();
    grid.params = GasParams();
    grid.source = HeatSource();
    grid.combustion = Combustion();
    grid.magnetic = MagneticField();
    grid.vessel = nullptr;
    rigid.params = RigidParams();
    releaseTime = 0;
    if (scene_) scene_->configure(*this);
    buildObstacleGeometry();
    reset();
    ++paramsVersion_;
}

void Simulation::reset() {
    time_ = 0;
    frame_ = 0;
    surfaceLoads_ = SurfaceLoads(); // the loads of the old geometry (filled again by the tunnel steps)
    if (scene_) scene_->build(*this);
    else useLiquidTank(kDefaultTank);
}

void Simulation::setSceneParam(int index, float value) {
    if (!scene_) return;
    scene_->setParam(index, value);
    reset();
    ++paramsVersion_;
}

// ---------------------------------------------------------------------------
// The boxes a scene builds in
// ---------------------------------------------------------------------------
void Simulation::useLiquidTank(const AABB& box) {
    mode_ = SimMode::Fluid;
    rigid.clear();
    rigid.setDomain(box);
    rigid.setStaticMesh(&obstacleBVH_);
    particles.setStaticMesh(&obstacleBVH_);
    particles.setRigidWorld(&rigid);
    particles.reset(box);
}

void Simulation::useGasBox(const Vector3& originFraction) {
    mode_ = SimMode::WindTunnel;
    const Vector3 size = grid.params.domainSize;
    const Vector3 origin(-originFraction.x * size.x, -originFraction.y * size.y, -originFraction.z * size.z);
    // A vessel (grid.vessel set by the scene) is the wall of the gas itself; the obstacle mesh is
    // then only drawn, not voxelised as a body in the flow.
    grid.reset(origin, grid.vessel ? nullptr : &obstacleBVH_);
    // Rigid bodies live in the gas box (walls = domain, the obstacle is static geometry); the
    // particle system (soft bodies, cloth, liquid) shares it.
    rigid.clear();
    rigid.setDomain(grid.domain());
    rigid.setStaticMesh(&obstacleBVH_);
    particles.setStaticMesh(&obstacleBVH_);
    particles.setRigidWorld(&rigid);
    particles.reset(grid.domain());
    gasImpulse_.clear();
    gasAngularImpulse_.clear();
    softGasImpulse_.clear();
    gasForceMax_ = 0;
}

void Simulation::useRigidArena(const AABB& box) {
    mode_ = SimMode::Rigid;
    rigid.clear();
    rigid.setDomain(box);
    rigid.setStaticMesh(&obstacleBVH_);
    particles.setStaticMesh(&obstacleBVH_); // soft bodies / cloth can be added in this mode too
    particles.setRigidWorld(&rigid);
    particles.reset(box);
}

// ---------------------------------------------------------------------------
// Stepping
// ---------------------------------------------------------------------------
void Simulation::stepFrame() {
    auto t0 = std::chrono::steady_clock::now();
    Probe::beginFrame();
    const long long allocations0 = Probe::allocations.load();
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
                lastGridDt_ = grid.step(float(frameDt - (grid.time() - time_)));
            time_ = grid.time();
        } else {
            stepGasWithBodies();
            time_ += frameDt;
        }
        {
            Probe::Timer timer("gas/surface loads ms"); // the pressure and friction on the obstacle's triangles
            updateSurfaceLoads();
        }
        break;
    }
    if (scene_) {
        Probe::Timer timer("scene/controllers ms"); // the scene's own work: magnets, emitters, timed events
        scene_->afterStep(*this);
    }
    ++frame_;
    lastStepMs_ = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
    Probe::set("frame/step ms", lastStepMs_);
    Probe::set("frame/time s", time_);
    Probe::set("memory/allocations per frame", double(Probe::allocations.load() - allocations0));
    publishChannels(); // after the counters: measuring is not part of the frame's physics
}

// The physical channels (scene/Channels.h): measured only when asked, then each single value is
// also reported to the Probe under its id, so the charts and the CSV find it by name. A probe
// line's profile stays in the measurements (a chart of a profile is not a time series).
void Simulation::publishChannels() {
    if (channels.empty()) {
        measurements_.clear(); // nothing asked: nothing measured (clear keeps the memory)
        return;
    }
    measurements_ = measureAll(*this, channels);
    for (const Measurement& m : measurements_)
        if (m.profile.empty()) Probe::set(m.info.id.c_str(), m.value);
}

int Simulation::sliceLayer() const {
    int n[3] = {grid.nx(), grid.ny(), grid.nz()};
    int axis = clampv(vis.sliceAxis, 0, 2);
    return clampv(int(vis.slicePosition * n[axis]), 0, n[axis] - 1);
}

} // namespace rf
