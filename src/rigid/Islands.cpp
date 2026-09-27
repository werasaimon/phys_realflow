// Islands of RigidWorld: the connected groups of touching bodies (union-find over contacts and
// joints), their sleep and waking, and bodies frozen for a step (sleepers, held bodies).
#include "rigid/RigidWorld.h"

#include <algorithm>

namespace rf {

bool RigidWorld::wakeIsland(int i) {
    const int island = bodies_[i].sleepIsland;
    bool any = false;
    for (RigidBody& b : bodies_) {
        bool member = &b == &bodies_[i] || (island >= 0 && b.sleepIsland == island);
        if (!member) continue;
        any |= b.sleeping;
        b.sleeping = false;
        b.sleepIsland = -1;
    }
    return any;
}

void RigidWorld::wake(int i) {
    wakeIsland(i);
    bodies_[i].sleepTimer = 0;
}

size_t RigidWorld::sleepingCount() const {
    size_t n = 0;
    for (const RigidBody& b : bodies_) n += b.sleeping;
    return n;
}

RigidWorld::Frozen RigidWorld::makeStatic(int i) {
    RigidBody& b = bodies_[i];
    const Frozen f{i, b.invMass, b.invInertiaLocal};
    b.invMass = 0;
    b.invInertiaLocal = Vector3(0.0f);
    b.updateInertia();
    return f;
}

void RigidWorld::restore(const Frozen& f) {
    RigidBody& b = bodies_[f.body];
    b.invMass = f.invMass;
    b.invInertiaLocal = f.invInertiaLocal;
    b.updateInertia();
}

void RigidWorld::hold(int i) {
    if (i < 0 || i >= int(bodies_.size()) || bodies_[i].invMass == 0) return; // static already
    held_.push_back(makeStatic(i));
}

void RigidWorld::releaseHeld() {
    for (const Frozen& f : held_) {
        if (f.body >= int(bodies_.size())) continue;
        restore(f);
        wake(f.body);
    }
    held_.clear();
}

// Sleeping bodies take part in the step as static bodies: their inverse mass is set to zero for
// the duration of the step (restored in unfreezeAll), so every solver path treats them as fixed.
void RigidWorld::freezeSleepers() {
    frozen_.clear();
    for (int i = 0; i < int(bodies_.size()); ++i) {
        RigidBody& b = bodies_[i];
        if (!b.sleeping || b.invMass == 0) continue;
        b.vel = b.angVel = Vector3(0.0f);
        frozen_.push_back(makeStatic(i));
    }
}

void RigidWorld::unfreezeAll(bool onlyAwake) {
    std::vector<Frozen>& keep = frozenKeep_; // kept between steps (swapped with frozen_ below)
    keep.clear();
    for (const Frozen& f : frozen_) {
        if (onlyAwake && bodies_[f.body].sleeping) { keep.push_back(f); continue; }
        restore(f);
    }
    frozen_.swap(keep);
}

// Islands = connected components of the contact graph among dynamic bodies (union-find).
// Before the solve: an island containing any moving body wakes up entirely (and its frozen
// members get their mass back). After the solve: an island whose bodies have all been slow for
// sleepTime falls asleep.
bool RigidWorld::updateIslands(bool decideSleep, float dt) {
    const int n = int(bodies_.size());
    islandParent_.resize(n);
    for (int i = 0; i < n; ++i) islandParent_[i] = i;
    auto find = [&](int x) {
        while (islandParent_[x] != x) x = islandParent_[x] = islandParent_[islandParent_[x]];
        return x;
    };
    for (const Manifold& m : manifolds_) {
        if (m.b < 0 || bodies_[m.a].mass <= 0 || bodies_[m.b].mass <= 0) continue;
        bool touching = false;
        for (const SolverPoint& p : m.points) touching |= p.depth > -params.slop;
        if (touching) islandParent_[find(m.a)] = find(m.b);
    }
    for (const auto& j : joints_)
        if (j->b >= 0 && bodies_[j->a].mass > 0 && bodies_[j->b].mass > 0) islandParent_[find(j->a)] = find(j->b);
    std::vector<float> minTimer(n, kInf);
    for (int i = 0; i < n; ++i)
        if (bodies_[i].mass > 0) minTimer[find(i)] = std::min(minTimer[find(i)], bodies_[i].sleepTimer);
    bool woke = false;
    std::vector<int> newIsland(n, -1);
    for (int i = 0; i < n; ++i) {
        RigidBody& b = bodies_[i];
        if (b.mass <= 0) continue;
        bool islandRests = minTimer[find(i)] >= params.sleepTime;
        if (!decideSleep) {
            // Touched by a moving body: the whole sleeping island wakes (Box2D).
            if (b.sleeping && !islandRests) woke |= wakeIsland(i);
        } else if (islandRests && !b.sleeping) {
            int r = find(i);
            if (newIsland[r] < 0) newIsland[r] = nextIsland_++;
            b.sleeping = true;
            b.sleepIsland = newIsland[r];
            b.vel = b.angVel = Vector3(0.0f);
        }
    }
    if (woke) unfreezeAll(true); // woken bodies get their mass back for this very step
    return woke;
    (void)dt;
}

} // namespace rf
