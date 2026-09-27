// The research view of the rigid solver: what it computed in the last step, drawn by layer into
// the Probe (see core/Probe.h for the layers and their analogues in PhysX / Bullet / Box2D) - the
// contacts with their normals, forces and depths, the bodies' boxes, mass frames and velocities,
// the contact islands, the world's AABB tree and the static mesh's BVH, the joints, and GJK / EPA
// of one watched body pair. Every group tests its layers first: with the layers off nothing here
// runs.
#include "rigid/RigidWorld.h"

#include "core/Format.h"
#include "core/Probe.h"
#include "rigid/GjkEpa.h"

#include <algorithm>
#include <cmath>

namespace rf {

namespace {

// A colour per index, well spread (golden-ratio hue steps), for islands.
Vector3 indexColor(int i) {
    const float h = std::fmod(0.618034f * float(i), 1.0f) * 6.0f;
    const float x = 1.0f - std::fabs(std::fmod(h, 2.0f) - 1.0f);
    Vector3 c = h < 1 ? Vector3(1, x, 0) : h < 2 ? Vector3(x, 1, 0) : h < 3 ? Vector3(0, 1, x)
              : h < 4 ? Vector3(0, x, 1) : h < 5 ? Vector3(x, 0, 1) : Vector3(1, 0, x);
    return c * 0.8f + Vector3(0.2f);
}

// Union-find root with path halving.
int findRoot(std::vector<int>& parent, int i) {
    while (parent[size_t(i)] != i) {
        parent[size_t(i)] = parent[size_t(parent[size_t(i)])];
        i = parent[size_t(i)];
    }
    return i;
}

} // namespace

// Every rigid layer, after the contacts are solved. A frame has several steps (substeps): each
// replaces the drawing of the one before, so the picture is that of the frame's last step.
void RigidWorld::drawDebug() const {
    Probe::clearLayers(Probe::bits(DrawLayer::ContactPoints, DrawLayer::WitnessPoints));
    drawContacts();
    drawBodyFrames();
    drawIslands();
    drawTrees();
    drawJoints();
    drawWatchedPair();
}

// Contact points, their normals (0.1 m), the normal force (impulse / dt, kForceScale metres per
// newton) and the penetration depth as a red segment along the normal.
void RigidWorld::drawContacts() const {
    const bool points = Probe::layerOn(DrawLayer::ContactPoints), normals = Probe::layerOn(DrawLayer::ContactNormals);
    const bool forces = Probe::layerOn(DrawLayer::ContactImpulses), depths = Probe::layerOn(DrawLayer::PenetrationDepth);
    if (!points && !normals && !forces && !depths) return;
    const float dt = std::max(lastDt_, 1e-6f);
    for (const Manifold& m : manifolds_)
        for (const SolverPoint& p : m.points) {
            if (points) Probe::point(DrawLayer::ContactPoints, p.position, Vector3(1.0f, 0.3f, 0.2f), 0.012f);
            if (normals) Probe::arrow(DrawLayer::ContactNormals, p.position, p.normal * 0.1f, Vector3(1.0f, 0.6f, 0.2f));
            if (forces && p.jn > 0)
                Probe::arrow(DrawLayer::ContactImpulses, p.position, p.normal * (p.jn / dt * Probe::kForceScale), Vector3(1.0f, 1.0f, 0.3f));
            if (depths && p.depth > 0)
                Probe::line(DrawLayer::PenetrationDepth, p.position, p.position - p.normal * p.depth, Vector3(1.0f, 0.1f, 0.1f));
        }
}

// Per body: its world box (awake: blue, sleeping: dim grey), the centre of mass as a small cross,
// the principal axes scaled by sqrt(I / m) (the radius of gyration), and the velocities: the
// linear one as the distance covered in 0.1 s, the angular one as 0.05 m per rad/s.
void RigidWorld::drawBodyFrames() const {
    const bool boxes = Probe::layerOn(DrawLayer::BodyAabbs), sleeping = Probe::layerOn(DrawLayer::Sleeping);
    const bool com = Probe::layerOn(DrawLayer::CentreOfMass), inertia = Probe::layerOn(DrawLayer::InertiaAxes);
    const bool vel = Probe::layerOn(DrawLayer::Velocities);
    if (!boxes && !sleeping && !com && !inertia && !vel) return;
    for (const RigidBody& b : bodies_) {
        if (!b.alive || b.invMass == 0) continue;
        if (boxes && !b.sleeping) Probe::box(DrawLayer::BodyAabbs, b.worldBounds(), Vector3(0.3f, 0.8f, 1.0f));
        if (sleeping && b.sleeping) Probe::box(DrawLayer::Sleeping, b.worldBounds(), Vector3(0.35f));
        const Matrix3x3 R = b.rotation();
        const float s = 0.15f * b.boundingRadius();
        for (int axis = 0; axis < 3 && com; ++axis) {
            Vector3 e(0.0f);
            e[axis] = s;
            Probe::line(DrawLayer::CentreOfMass, b.pos - R * e, b.pos + R * e, Vector3(1.0f));
        }
        const float invI[3] = {b.invInertiaLocal.x, b.invInertiaLocal.y, b.invInertiaLocal.z};
        for (int axis = 0; axis < 3 && inertia; ++axis) {
            Vector3 e(0.0f);
            e[axis] = invI[axis] > 0 ? std::sqrt(1.0f / (invI[axis] * b.mass)) : 0.0f; // radius of gyration
            Vector3 col(0.2f);
            col[axis] = 1.0f;
            Probe::line(DrawLayer::InertiaAxes, b.pos - R * e, b.pos + R * e, col);
        }
        if (vel) {
            Probe::arrow(DrawLayer::Velocities, b.pos, b.vel * 0.1f, Vector3(0.2f, 1.0f, 1.0f));
            Probe::arrow(DrawLayer::Velocities, b.pos, b.angVel * 0.05f, Vector3(1.0f, 0.3f, 1.0f));
        }
    }
}

// Contact islands: bodies joined by contacts (or joints) between two dynamic bodies share a
// colour. Found here with a union-find over the manifolds - only while the layer is on.
void RigidWorld::drawIslands() const {
    if (!Probe::layerOn(DrawLayer::Islands)) return;
    std::vector<int> parent(bodies_.size());
    for (size_t i = 0; i < parent.size(); ++i) parent[i] = int(i);
    auto join = [&](int a, int b) {
        if (a < 0 || b < 0 || bodies_[size_t(a)].invMass == 0 || bodies_[size_t(b)].invMass == 0) return;
        parent[size_t(findRoot(parent, a))] = findRoot(parent, b);
    };
    for (const Manifold& m : manifolds_) join(m.a, m.b);
    for (const auto& j : joints_) join(j->a, j->b);
    for (size_t i = 0; i < bodies_.size(); ++i) {
        const RigidBody& b = bodies_[i];
        if (!b.alive || b.invMass == 0) continue;
        Probe::box(DrawLayer::Islands, b.worldBounds(), indexColor(findRoot(parent, int(i))));
    }
}

// The world's dynamic AABB tree (every node, coloured by depth: the root red, the leaves blue) and
// the static mesh's BVH down to Probe::kBvhDepth levels.
void RigidWorld::drawTrees() const {
    if (Probe::layerOn(DrawLayer::WorldTree) && worldTree_.root() >= 0) {
        const float height = float(std::max(worldTree_.height(), 1));
        std::vector<std::pair<int, int>> stack{{worldTree_.root(), 0}};
        while (!stack.empty()) {
            auto [i, depth] = stack.back();
            stack.pop_back();
            const AABBTree::Node& n = worldTree_.node(i);
            Probe::box(DrawLayer::WorldTree, n.box, heatColor(1.0f - float(depth) / height));
            if (!n.isLeaf()) {
                stack.push_back({n.left, depth + 1});
                stack.push_back({n.right, depth + 1});
            }
        }
    }
    if (Probe::layerOn(DrawLayer::MeshBvh) && mesh_ && !mesh_->empty()) {
        const std::vector<BVHNode>& nodes = mesh_->bvh().nodes();
        std::vector<std::pair<int, int>> stack{{0, 0}};
        while (!stack.empty()) {
            auto [i, depth] = stack.back();
            stack.pop_back();
            const BVHNode& n = nodes[size_t(i)];
            Probe::box(DrawLayer::MeshBvh, n.box, heatColor(1.0f - float(depth) / float(Probe::kBvhDepth)));
            if (!n.isLeaf() && depth < Probe::kBvhDepth) {
                stack.push_back({n.first, depth + 1});
                stack.push_back({n.first + 1, depth + 1});
            }
        }
    }
}

// Joints: both anchors, the link between them (it should have zero length for a ball joint) and
// the joint's axis (hinge / slider) as a 0.2 m arrow.
void RigidWorld::drawJoints() const {
    if (!Probe::layerOn(DrawLayer::JointFrames)) return;
    for (const auto& j : joints_) {
        const Vector3 a = j->worldAnchorA(bodies_), b = j->worldAnchorB(bodies_);
        Probe::point(DrawLayer::JointFrames, a, Vector3(0.2f, 1.0f, 0.4f), 0.015f);
        Probe::point(DrawLayer::JointFrames, b, Vector3(0.2f, 0.6f, 1.0f), 0.015f);
        Probe::line(DrawLayer::JointFrames, a, b, Vector3(1.0f, 1.0f, 0.3f));
        const Vector3 axis = j->worldAxis(bodies_);
        if (length2(axis) > 0) Probe::arrow(DrawLayer::JointFrames, a, normalize(axis) * 0.2f, Vector3(0.9f));
    }
}

} // namespace rf
