#include "spatial/AABBTree.h"

#include <algorithm>

namespace rf {

static AABB unite(const AABB& a, const AABB& b) {
    AABB r = a;
    r.expand(b);
    return r;
}

static bool containsBox(const AABB& outer, const AABB& inner) {
    return outer.lo.x <= inner.lo.x && outer.lo.y <= inner.lo.y && outer.lo.z <= inner.lo.z && inner.hi.x <= outer.hi.x &&
           inner.hi.y <= outer.hi.y && inner.hi.z <= outer.hi.z;
}

AABB AABBTree::fat(const AABB& b) const {
    Vector3 m = b.extent() * fatten_;
    return AABB(b.lo - m, b.hi + m);
}

// ---------------------------------------------------------------------------
// Node pool
// ---------------------------------------------------------------------------
int AABBTree::allocateNode() {
    int i;
    if (freeList_ >= 0) {
        i = freeList_;
        freeList_ = nodes_[i].parent;
    } else {
        i = int(nodes_.size());
        nodes_.emplace_back();
    }
    nodes_[i] = Node();
    return i;
}

void AABBTree::freeNode(int i) {
    nodes_[i].parent = freeList_;
    nodes_[i].height = -1; // marks a free node
    freeList_ = i;
}

void AABBTree::clear() {
    nodes_.clear();
    root_ = freeList_ = -1;
    leafCount_ = 0;
}

// ---------------------------------------------------------------------------
// Public operations
// ---------------------------------------------------------------------------
int AABBTree::insert(const AABB& box, int object) {
    int leaf = allocateNode();
    nodes_[leaf].box = fat(box);
    nodes_[leaf].object = object;
    insertLeaf(leaf);
    ++leafCount_;
    return leaf;
}

void AABBTree::remove(int proxy) {
    removeLeaf(proxy);
    freeNode(proxy);
    --leafCount_;
}

bool AABBTree::update(int proxy, const AABB& box) {
    const AABB& stored = nodes_[proxy].box;
    if (fatten_ > 0 && containsBox(stored, box)) return false; // still inside its fat box
    if (fatten_ == 0 && containsBox(stored, box) && containsBox(box, stored)) return false; // exact: unchanged
    removeLeaf(proxy);
    nodes_[proxy].box = fat(box);
    insertLeaf(proxy);
    return true;
}

void AABBTree::findPairs(std::vector<std::pair<int, int>>& pairs) const {
    pairs.clear();
    for (const Node& n : nodes_) {
        if (n.height != 0) continue; // leaves only (free nodes have height -1)
        query(n.box, [&](int other) {
            if (n.object < other) pairs.emplace_back(n.object, other); // each pair once
        });
    }
    std::sort(pairs.begin(), pairs.end());
}

// ---------------------------------------------------------------------------
// Insertion: descend to the sibling with the lowest SAH cost
// ---------------------------------------------------------------------------
void AABBTree::insertLeaf(int leaf) {
    if (root_ < 0) {
        root_ = leaf;
        nodes_[leaf].parent = -1;
        return;
    }
    const AABB box = nodes_[leaf].box;
    int index = root_;
    while (!nodes_[index].isLeaf()) {
        const Node& n = nodes_[index];
        const float area = n.box.surfaceArea();
        const float combined = unite(n.box, box).surfaceArea();
        // Making a new parent for this node and the leaf here costs:
        const float costHere = 2.0f * combined;
        // Going further down, every node on the way grows by this much:
        const float inherited = 2.0f * (combined - area);
        auto costDown = [&](int child) {
            const Node& c = nodes_[child];
            float grown = unite(box, c.box).surfaceArea();
            return (c.isLeaf() ? grown : grown - c.box.surfaceArea()) + inherited;
        };
        const float costLeft = costDown(n.left), costRight = costDown(n.right);
        if (costHere < costLeft && costHere < costRight) break;
        index = costLeft < costRight ? n.left : n.right;
    }
    const int sibling = index;

    // A new parent takes the sibling's place and gets the sibling and the leaf as children.
    const int oldParent = nodes_[sibling].parent;
    const int parent = allocateNode();
    nodes_[parent].parent = oldParent;
    nodes_[parent].box = unite(box, nodes_[sibling].box);
    nodes_[parent].height = nodes_[sibling].height + 1;
    nodes_[parent].left = sibling;
    nodes_[parent].right = leaf;
    nodes_[sibling].parent = parent;
    nodes_[leaf].parent = parent;
    if (oldParent < 0) root_ = parent;
    else if (nodes_[oldParent].left == sibling) nodes_[oldParent].left = parent;
    else nodes_[oldParent].right = parent;

    refitUpwards(nodes_[leaf].parent);
}

// ---------------------------------------------------------------------------
// Removal: the leaf's sibling takes the parent's place
// ---------------------------------------------------------------------------
void AABBTree::removeLeaf(int leaf) {
    if (leaf == root_) {
        root_ = -1;
        return;
    }
    const int parent = nodes_[leaf].parent;
    const int grandParent = nodes_[parent].parent;
    const int sibling = nodes_[parent].left == leaf ? nodes_[parent].right : nodes_[parent].left;
    nodes_[sibling].parent = grandParent;
    if (grandParent < 0) {
        root_ = sibling;
    } else {
        if (nodes_[grandParent].left == parent) nodes_[grandParent].left = sibling;
        else nodes_[grandParent].right = sibling;
    }
    freeNode(parent);
    nodes_[leaf].parent = -1;
    if (grandParent >= 0) refitUpwards(grandParent);
}

void AABBTree::refitUpwards(int i) {
    while (i >= 0) {
        i = balance(i);
        Node& n = nodes_[i];
        n.height = 1 + std::max(nodes_[n.left].height, nodes_[n.right].height);
        n.box = unite(nodes_[n.left].box, nodes_[n.right].box);
        i = n.parent;
    }
}

// ---------------------------------------------------------------------------
// Balancing: if one child of A is more than one level taller, rotate it up (Box2D).
//
//          A                    C
//        /   \                /   \
//       B     C      ->      A     F     (G taller than F: the other way round)
//            / \            / \
//           F   G          B   G
// ---------------------------------------------------------------------------
int AABBTree::balance(int iA) {
    Node& A = nodes_[iA];
    if (A.isLeaf() || A.height < 2) return iA;
    const int iB = A.left, iC = A.right;
    const int diff = nodes_[iC].height - nodes_[iB].height;
    if (diff >= -1 && diff <= 1) return iA;

    // Rotate the taller child U (C if diff > 0, else B) up into A's place; A keeps the other
    // child S and takes the shorter grandchild; U keeps its taller grandchild.
    const bool rightTaller = diff > 0;
    const int iU = rightTaller ? iC : iB, iS = rightTaller ? iB : iC;
    Node& U = nodes_[iU];
    const int iX = U.left, iY = U.right;
    const int iTall = nodes_[iX].height > nodes_[iY].height ? iX : iY;
    const int iShort = iTall == iX ? iY : iX;

    // U replaces A under A's parent.
    U.parent = A.parent;
    if (U.parent < 0) root_ = iU;
    else if (nodes_[U.parent].left == iA) nodes_[U.parent].left = iU;
    else nodes_[U.parent].right = iU;

    // A becomes a child of U, keeping S and adopting U's shorter child.
    U.left = iA;
    U.right = iTall;
    A.parent = iU;
    if (rightTaller) A.right = iShort; // A = (S, short)
    else A.left = iShort;              // A = (short, S)
    nodes_[iShort].parent = iA;

    A.box = unite(nodes_[iS].box, nodes_[iShort].box);
    A.height = 1 + std::max(nodes_[iS].height, nodes_[iShort].height);
    U.box = unite(A.box, nodes_[iTall].box);
    U.height = 1 + std::max(A.height, nodes_[iTall].height);
    return iU;
}

// ---------------------------------------------------------------------------
// Validation (tests)
// ---------------------------------------------------------------------------
bool AABBTree::validateNode(int i) const {
    const Node& n = nodes_[i];
    if (n.isLeaf()) return n.height == 0 && n.right < 0 && n.object >= 0;
    const Node &l = nodes_[n.left], &r = nodes_[n.right];
    if (l.parent != i || r.parent != i) return false;
    if (n.height != 1 + std::max(l.height, r.height)) return false;
    if (!containsBox(n.box, l.box) || !containsBox(n.box, r.box)) return false;
    return validateNode(n.left) && validateNode(n.right);
}

bool AABBTree::validate() const {
    if (root_ < 0) return leafCount_ == 0;
    if (nodes_[root_].parent != -1 || !validateNode(root_)) return false;
    int leaves = 0;
    for (const Node& n : nodes_) leaves += n.height == 0;
    return leaves == leafCount_;
}

} // namespace rf
