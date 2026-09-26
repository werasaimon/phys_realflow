#pragma once
// Dynamic AABB tree: a bounding volume hierarchy that is updated incrementally while objects move
// (Box2D b2DynamicTree, Bullet btDbvt; E. Catto, "Dynamic Bounding Volume Hierarchies", GDC 2019).
//
//  * every object is a leaf holding a "fat" box (the object's box enlarged by a margin), so an
//    object that jiggles in place does not touch the tree at all;
//  * insert: the leaf goes next to the sibling that increases the total surface area the least
//    (surface area heuristic, SAH), found by descending the tree;
//  * remove: the leaf's parent is replaced by the leaf's sibling;
//  * after every change the path to the root is refitted and rebalanced by tree rotations (as in
//    AVL trees), so the height stays O(log n).
//
// Static triangle meshes use the BVH (built once, top-down) in BVH.h instead.

#include "math/Math.h"

#include <utility>
#include <vector>

namespace rf {

class AABBTree {
public:
    struct Node {
        AABB box;          // leaf: fat box of its object; interior: union of the two children
        int parent = -1;
        int left = -1;     // children, -1 for a leaf
        int right = -1;
        int object = -1;   // leaf: user id of the object
        int height = 0;    // leaf = 0, interior = 1 + max(children)
        bool isLeaf() const { return left < 0; }
    };

    // fatten: the stored box is the object's box enlarged by fatten * its size on every side
    // (0 = exact boxes).
    explicit AABBTree(float fatten = 0.1f) : fatten_(fatten) {}

    // Adds an object; returns its proxy (leaf node index) for update() / remove().
    int insert(const AABB& box, int object);
    void remove(int proxy);
    // New box of the object; the tree changes only if it left its fat box. Returns true if it did.
    bool update(int proxy, const AABB& box);
    void clear();

    // fn(object) for every object whose (fat) box overlaps the query box.
    template <class F> void query(const AABB& box, F&& fn) const;
    // fn(object) for every object whose (fat) box the ray o + t d, 0 <= t <= maxT, passes through.
    template <class F> void raycast(const Vector3& origin, const Vector3& dir, float maxT, F&& fn) const;
    // All pairs of objects (a < b) with overlapping fat boxes.
    void findPairs(std::vector<std::pair<int, int>>& pairs) const;

    int root() const { return root_; }
    int height() const { return root_ < 0 ? 0 : nodes_[root_].height; }
    int objectCount() const { return leafCount_; }
    const Node& node(int i) const { return nodes_[i]; }
    const AABB& fatBox(int proxy) const { return nodes_[proxy].box; }
    // Checks every invariant (links, heights, boxes contain children, balance); for tests.
    bool validate() const;

private:
    int allocateNode();
    void freeNode(int i);
    void insertLeaf(int leaf);
    void removeLeaf(int leaf);
    int balance(int a);              // one rotation at node a if unbalanced; returns the new subtree root
    void refitUpwards(int i);        // boxes and heights from i to the root, rebalancing on the way
    AABB fat(const AABB& b) const;
    bool validateNode(int i) const;

    float fatten_;
    std::vector<Node> nodes_;
    int root_ = -1;
    int freeList_ = -1; // free nodes are chained through Node::parent
    int leafCount_ = 0;
};

// ---------------------------------------------------------------------------
template <class F> void AABBTree::query(const AABB& box, F&& fn) const {
    if (root_ < 0) return;
    std::vector<int> stack{root_};
    while (!stack.empty()) {
        const Node& n = nodes_[stack.back()];
        stack.pop_back();
        if (!n.box.overlaps(box)) continue;
        if (n.isLeaf()) {
            fn(n.object);
        } else {
            stack.push_back(n.left);
            stack.push_back(n.right);
        }
    }
}

template <class F> void AABBTree::raycast(const Vector3& origin, const Vector3& dir, float maxT, F&& fn) const {
    if (root_ < 0) return;
    const Vector3 inv(1.0f / dir.x, 1.0f / dir.y, 1.0f / dir.z); // +-inf for axis-parallel rays is fine
    std::vector<int> stack{root_};
    while (!stack.empty()) {
        const Node& n = nodes_[stack.back()];
        stack.pop_back();
        if (n.box.rayHit(origin, inv, maxT) == kInf) continue;
        if (n.isLeaf()) {
            fn(n.object);
        } else {
            stack.push_back(n.left);
            stack.push_back(n.right);
        }
    }
}

} // namespace rf
