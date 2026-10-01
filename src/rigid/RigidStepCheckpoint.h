#pragma once
// Persistent mechanical state needed to undo a rigid impulse-solver trial. Geometry and joint
// parameters are immutable during a step. Broad-phase caches are refitted after restoration.
#include "rigid/RigidWorld.h"

namespace rf {

class RigidStepCheckpoint {
public:
    explicit RigidStepCheckpoint(const RigidWorld& world) { save(world); }
    void save(const RigidWorld& world);
    void restore(RigidWorld& world) const;
    void restoreLoads(RigidWorld& world) const;

private:
    std::vector<RigidBody> bodies_;
    std::vector<RigidWorld::Manifold> manifolds_;
    // A contiguous snapshot avoids one allocation per newly appearing contact-cache node.
    std::vector<std::pair<uint64_t, RigidWorld::CachedPair>> cache_;
    std::vector<RigidWorld::Frozen> frozen_, held_;
    std::vector<std::vector<JacobianRow>> rows_;
    std::vector<std::vector<float>> warm_;
    std::vector<AABB> boxes_;
    std::vector<std::pair<int, int>> pairs_;
    std::vector<int> levels_, islandParent_, treeProxies_;
    std::vector<char> clamped_;
    AABBTree tree_;
    RigidWorld::GrabJoint grab_;
    RigidWorld::WarmStartStats warmStats_;
    RigidWorld::Timings timings_;
    CcdDiagnostics ccd_;
    uint32_t cacheStamp_ = 0;
    int nextIsland_ = 0;
    size_t contacts_ = 0, hits_ = 0;
    float lastDt_ = 0;
};

} // namespace rf
