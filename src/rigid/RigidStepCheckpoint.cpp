// Save/restore for rejected impulse-solver trials. The same external force and torque act
// throughout every accepted subdivision of the caller's interval; rejected trials consume none.
#include "rigid/RigidStepCheckpoint.h"
#include "core/Probe.h"

namespace rf {

void RigidStepCheckpoint::save(const RigidWorld& w) {
    Probe::Timer timer("rigid/checkpoint ms");
    bodies_ = w.bodies_; manifolds_ = w.manifolds_;
    cache_.assign(w.cache_.begin(), w.cache_.end());
    frozen_ = w.frozen_; held_ = w.held_; grab_ = w.grab_;
    boxes_ = w.boxes_; pairs_ = w.pairs_; levels_ = w.levels_;
    islandParent_ = w.islandParent_; treeProxies_ = w.treeProxies_; tree_ = w.worldTree_;
    clamped_ = w.ccdClamped_; ccd_ = w.ccdDiagnostics_; hits_ = w.ccdHits_;
    cacheStamp_ = w.cacheStamp_; nextIsland_ = w.nextIsland_; lastDt_ = w.lastDt_;
    contacts_ = w.contactCount_; warmStats_ = w.warmStats_; timings_ = w.timings_;
    rows_.resize(w.joints_.size()); warm_.resize(w.joints_.size());
    for (size_t i = 0; i < w.joints_.size(); ++i) {
        rows_[i] = w.joints_[i]->rows_;
        warm_[i] = w.joints_[i]->warm_;
    }
}

void RigidStepCheckpoint::restore(RigidWorld& w) const {
    Probe::Timer timer("rigid/checkpoint ms");
    w.bodies_ = bodies_; w.manifolds_ = manifolds_;
    w.cache_.clear();
    w.cache_.insert(cache_.begin(), cache_.end());
    w.frozen_ = frozen_; w.held_ = held_; w.grab_ = grab_;
    w.boxes_ = boxes_; w.pairs_ = pairs_; w.levels_ = levels_;
    w.islandParent_ = islandParent_; w.treeProxies_ = treeProxies_; w.worldTree_ = tree_;
    w.ccdClamped_ = clamped_; w.ccdDiagnostics_ = ccd_; w.ccdHits_ = hits_;
    w.cacheStamp_ = cacheStamp_; w.nextIsland_ = nextIsland_; w.lastDt_ = lastDt_;
    w.contactCount_ = contacts_; w.warmStats_ = warmStats_; w.timings_ = timings_;
    for (size_t i = 0; i < w.joints_.size(); ++i) {
        w.joints_[i]->rows_ = rows_[i];
        w.joints_[i]->warm_ = warm_[i];
    }
    w.broadphase_->update(boxes_);
}

void RigidStepCheckpoint::restoreLoads(RigidWorld& w) const {
    for (size_t i = 0; i < bodies_.size(); ++i) {
        w.bodies_[i].force = bodies_[i].force;
        w.bodies_[i].torque = bodies_[i].torque;
    }
}

} // namespace rf
