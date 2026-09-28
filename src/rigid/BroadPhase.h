#pragma once
// Broad phase: finds candidate body pairs whose (fattened) AABBs overlap. The collision pipeline is
//     BroadPhase (pairs)  ->  NarrowPhase (manifolds)  ->  ContactCache (warm start)  ->  solver
// Implementations are interchangeable (Strategy): brute force for reference/tests, static BVH
// (rebuilt every step), dynamic AABB tree, and incremental sweep and prune (default).

#include "spatial/AABBTree.h"
#include "spatial/BVH.h"
#include "math/Math.h"

#include <cstdint>
#include <unordered_set>
#include <utility>
#include <vector>

namespace rf {

class BroadPhase {
public:
    virtual ~BroadPhase() = default;
    virtual const char* name() const = 0;
    // Rebuild / refit from the current body bounds.
    virtual void update(const std::vector<AABB>& bounds) = 0;
    // Overlapping pairs (i < j), each reported once.
    virtual void findPairs(std::vector<std::pair<int, int>>& pairs) const = 0;
};

// O(n^2) reference implementation.
class BruteForceBroadPhase final : public BroadPhase {
public:
    const char* name() const override { return "brute force"; }
    void update(const std::vector<AABB>& bounds) override { bounds_ = bounds; }
    void findPairs(std::vector<std::pair<int, int>>& pairs) const override;

private:
    std::vector<AABB> bounds_;
};

// Binned-SAH BVH over the body AABBs; each body queries the tree.
class BvhBroadPhase final : public BroadPhase {
public:
    const char* name() const override { return "BVH"; }
    void update(const std::vector<AABB>& bounds) override;
    void findPairs(std::vector<std::pair<int, int>>& pairs) const override;

private:
    std::vector<AABB> bounds_;
    BVH bvh_;
};

// Dynamic AABB tree (as Box2D's b2BroadPhase): one persistent leaf per body with a fat box and a
// persistent set of pairs. A step re-inserts only the bodies that left their fat boxes ("moved"),
// drops the pairs whose fat boxes no longer overlap, and queries the tree for the moved bodies
// only - bodies at rest cost nothing.
class AABBTreeBroadPhase final : public BroadPhase {
public:
    explicit AABBTreeBroadPhase(float fatten = 0.0f) : tree_(fatten) {}
    const char* name() const override { return "AABB tree"; }
    void update(const std::vector<AABB>& bounds) override;
    void findPairs(std::vector<std::pair<int, int>>& pairs) const override { pairs = pairs_; }
    const AABBTree& tree() const { return tree_; }

private:
    AABBTree tree_;
    std::vector<int> proxies_;                 // leaf of every body
    std::vector<std::pair<int, int>> pairs_;   // overlapping fat boxes, sorted, a < b
};

// Incremental sweep and prune on three axes (Baraff 1992; Bullet btAxisSweep3, PhysX SAP).
// Each axis keeps the sorted list of interval endpoints (min/max of every box). Between frames the
// boxes move little, so re-sorting with insertion sort is nearly linear. Every swap of a min
// endpoint past a max endpoint is an event: "B.min moves below A.max" may start an overlap - the
// full 3D overlap of the new boxes is tested and the pair is added; "B.max moves below A.min" ends
// one - the pair is removed. Each pair of endpoints swaps at most once per update, so after the
// three axes are sorted the pair set is exactly the set of overlapping boxes.
//
// Fat boxes (as in Bullet's dbvt, Box2D, PhysX): with fatten > 0 every box is stored enlarged by
// fatten * its size, and its endpoints move only when the real box leaves the fat one. Bodies that
// jiggle in place (piles, stacks) then cause no endpoint swaps at all; the reported pairs are the
// overlapping fat boxes (a superset - the narrow phase rejects the extra ones).
class SweepAndPruneBroadPhase final : public BroadPhase {
public:
    explicit SweepAndPruneBroadPhase(float fatten = 0.0f) : fatten_(fatten) {}
    const char* name() const override { return "sweep and prune"; }
    void update(const std::vector<AABB>& bounds) override;
    void findPairs(std::vector<std::pair<int, int>>& pairs) const override;
    // Diagnostics: endpoint swaps of the last update (small when the motion is coherent).
    size_t lastSwaps() const { return swaps_; }

private:
    struct Endpoint {
        float value;
        uint32_t data; // body << 1 | isMax
    };
    void rebuild();
    void sortAxis(int axis);
    static uint64_t key(uint32_t a, uint32_t b) { return a < b ? (uint64_t(a) << 32) | b : (uint64_t(b) << 32) | a; }

    float fatten_;
    std::vector<AABB> bounds_; // stored (fat) boxes
    std::vector<Endpoint> axes_[3];
    std::unordered_set<uint64_t> pairs_;
    mutable std::vector<uint64_t> keys_; // scratch of findPairs (the sorted output), kept between steps
    mutable std::vector<std::pair<int, int>> sortedPairs_; // the last sorted output
    mutable bool pairsChanged_ = true;   // the set changed since sortedPairs_ was made
    size_t swaps_ = 0;
};

} // namespace rf
