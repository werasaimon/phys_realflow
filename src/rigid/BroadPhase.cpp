#include "rigid/BroadPhase.h"

#include "core/Parallel.h"

#include <algorithm>

namespace rf {

void BruteForceBroadPhase::findPairs(std::vector<std::pair<int, int>>& pairs) const {
    pairs.clear();
    for (int i = 0; i < int(bounds_.size()); ++i)
        for (int j = i + 1; j < int(bounds_.size()); ++j)
            if (bounds_[i].overlaps(bounds_[j])) pairs.emplace_back(i, j);
}

void BvhBroadPhase::update(const std::vector<AABB>& bounds) {
    bounds_ = bounds;
    bvh_.build(bounds_, 2);
}

void BvhBroadPhase::findPairs(std::vector<std::pair<int, int>>& pairs) const {
    // Queries are independent: run them in parallel, then concatenate in body order (deterministic).
    const int n = int(bounds_.size());
    std::vector<std::vector<std::pair<int, int>>> local(n);
    parallelFor(n, [&](int i) {
        bvh_.queryAABB(bounds_[i], [&](uint32_t j) {
            // Leaves report all their primitives: test the primitive's own box too.
            if (int(j) > i && bounds_[i].overlaps(bounds_[j])) local[i].emplace_back(i, int(j));
        });
    }, 64);
    pairs.clear();
    for (auto& l : local) pairs.insert(pairs.end(), l.begin(), l.end());
}

// ---------------------------------------------------------------------------
// Dynamic AABB tree
// ---------------------------------------------------------------------------
void AABBTreeBroadPhase::update(const std::vector<AABB>& bounds) {
    if (bounds.size() != proxies_.size()) { // bodies added or removed: start over
        tree_.clear();
        proxies_.resize(bounds.size());
        for (size_t i = 0; i < bounds.size(); ++i) proxies_[i] = tree_.insert(bounds[i], int(i));
        tree_.findPairs(pairs_);
        return;
    }
    std::vector<int> moved;
    for (size_t i = 0; i < bounds.size(); ++i)
        if (tree_.update(proxies_[i], bounds[i])) moved.push_back(int(i));
    if (moved.empty()) return;
    // Pairs whose fat boxes separated (only a moved body's box can have changed).
    std::vector<char> isMoved(bounds.size(), 0);
    for (int i : moved) isMoved[i] = 1;
    std::vector<std::pair<int, int>> kept;
    kept.reserve(pairs_.size());
    for (const auto& pr : pairs_)
        if ((!isMoved[pr.first] && !isMoved[pr.second]) ||
            tree_.fatBox(proxies_[pr.first]).overlaps(tree_.fatBox(proxies_[pr.second])))
            kept.push_back(pr);
    // New pairs: only the moved bodies ask the tree.
    for (int i : moved)
        tree_.query(tree_.fatBox(proxies_[i]), [&](int j) {
            if (j != i) kept.emplace_back(std::min(i, j), std::max(i, j));
        });
    std::sort(kept.begin(), kept.end());
    kept.erase(std::unique(kept.begin(), kept.end()), kept.end());
    pairs_.swap(kept);
}

// ---------------------------------------------------------------------------
// Sweep and prune
// ---------------------------------------------------------------------------
void SweepAndPruneBroadPhase::update(const std::vector<AABB>& bounds) {
    const bool sameBodies = bounds.size() == bounds_.size() && !axes_[0].empty();
    auto fat = [&](const AABB& b) {
        Vector3 m = b.extent() * fatten_;
        return AABB(b.lo - m, b.hi + m);
    };
    if (!sameBodies) {
        bounds_.resize(bounds.size());
        for (size_t i = 0; i < bounds.size(); ++i) bounds_[i] = fat(bounds[i]);
        rebuild();
        return;
    }
    bool moved = false;
    for (size_t i = 0; i < bounds.size(); ++i) {
        const AABB& s = bounds_[i];
        const AABB& b = bounds[i];
        if (b.lo.x >= s.lo.x && b.lo.y >= s.lo.y && b.lo.z >= s.lo.z && b.hi.x <= s.hi.x && b.hi.y <= s.hi.y &&
            b.hi.z <= s.hi.z && (fatten_ > 0 || (b.lo.x == s.lo.x && b.lo.y == s.lo.y && b.lo.z == s.lo.z && b.hi.x == s.hi.x &&
                             b.hi.y == s.hi.y && b.hi.z == s.hi.z)))
            continue; // still inside its fat box: endpoints stay
        bounds_[i] = fat(b);
        moved = true;
    }
    swaps_ = 0;
    if (!moved) return;
    for (int a = 0; a < 3; ++a) {
        for (Endpoint& e : axes_[a]) {
            const AABB& b = bounds_[e.data >> 1];
            e.value = (e.data & 1) ? b.hi[a] : b.lo[a];
        }
        sortAxis(a);
    }
}

void SweepAndPruneBroadPhase::sortAxis(int axis) {
    std::vector<Endpoint>& E = axes_[axis];
    for (size_t i = 1; i < E.size(); ++i) {
        const Endpoint e = E[i];
        size_t j = i;
        // Ties keep their order; a min equal to a max still counts as touching (AABB::overlaps
        // is inclusive), so min endpoints sort before max endpoints of the same value.
        auto before = [&](const Endpoint& x, const Endpoint& y) {
            return x.value < y.value || (x.value == y.value && !(x.data & 1) && (y.data & 1));
        };
        while (j > 0 && before(e, E[j - 1])) {
            const Endpoint& f = E[j - 1];
            const uint32_t be = e.data >> 1, bf = f.data >> 1;
            if (be != bf) {
                const bool eMax = e.data & 1, fMax = f.data & 1;
                if (!eMax && fMax) {
                    // e's min passes below f's max: the intervals start to overlap on this axis.
                    if (bounds_[be].overlaps(bounds_[bf])) pairs_.insert(key(be, bf));
                } else if (eMax && !fMax) {
                    // e's max passes below f's min: the intervals separate.
                    pairs_.erase(key(be, bf));
                }
            }
            E[j] = f;
            --j;
            ++swaps_;
        }
        E[j] = e;
    }
}

void SweepAndPruneBroadPhase::rebuild() {
    const uint32_t n = uint32_t(bounds_.size());
    pairs_.clear();
    swaps_ = 0;
    for (int a = 0; a < 3; ++a) {
        std::vector<Endpoint>& E = axes_[a];
        E.clear();
        E.reserve(2 * n);
        for (uint32_t i = 0; i < n; ++i) {
            E.push_back({bounds_[i].lo[a], i << 1});
            E.push_back({bounds_[i].hi[a], (i << 1) | 1});
        }
        std::stable_sort(E.begin(), E.end(), [](const Endpoint& x, const Endpoint& y) {
            return x.value < y.value || (x.value == y.value && !(x.data & 1) && (y.data & 1));
        });
    }
    if (n == 0) return;
    // Initial pairs: one sweep along the axis of largest spread, full 3D test for the candidates.
    int best = 0;
    float spread[3] = {0, 0, 0};
    {
        Vector3 mean(0.0f), mean2(0.0f);
        for (const AABB& b : bounds_) {
            Vector3 c = b.center();
            mean += c;
            mean2 += Vector3(c.x * c.x, c.y * c.y, c.z * c.z);
        }
        for (int a = 0; a < 3; ++a) spread[a] = mean2[a] / n - sqr(mean[a] / n);
        best = spread[1] > spread[best] ? 1 : best;
        best = spread[2] > spread[best] ? 2 : best;
    }
    std::vector<uint32_t> active;
    for (const Endpoint& e : axes_[best]) {
        const uint32_t b = e.data >> 1;
        if (e.data & 1) {
            active.erase(std::find(active.begin(), active.end(), b));
            continue;
        }
        for (uint32_t o : active)
            if (bounds_[o].overlaps(bounds_[b])) pairs_.insert(key(o, b));
        active.push_back(b);
    }
}

void SweepAndPruneBroadPhase::findPairs(std::vector<std::pair<int, int>>& pairs) const {
    // Sorted output: the solver order (and so the result) does not depend on hash-set order.
    std::vector<uint64_t> keys(pairs_.begin(), pairs_.end());
    std::sort(keys.begin(), keys.end());
    pairs.clear();
    pairs.reserve(keys.size());
    for (uint64_t k : keys) pairs.emplace_back(int(k >> 32), int(k & 0xffffffffu));
}

} // namespace rf
