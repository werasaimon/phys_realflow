// Non-convex models decomposed once into convex parts (cached, shared by all bodies).
#include "samples/Models.h"

#include "core/Mesh.h"
#include "rigid/ConvexDecomposition.h"

namespace rf {

// A solid given as overlapping closed parts -> minimal set of convex hulls fitted to its smooth surface.
static std::shared_ptr<const CompoundShape> decompose(const std::vector<TriMesh>& parts) {
    DecompositionParams prm;
    prm.resolution = 24;
    prm.maxHullVertices = 64;
    return std::make_shared<CompoundShape>(convexDecomposition(parts, prm), primitives::merge(parts));
}

std::shared_ptr<const CompoundShape> teapotShape() {
    static std::shared_ptr<const CompoundShape> shape = decompose(primitives::teapotParts(0.3f));
    return shape;
}

std::shared_ptr<const CompoundShape> bunnyShape() {
    static std::shared_ptr<const CompoundShape> shape = decompose(primitives::bunnyParts(0.3f));
    return shape;
}

// With the default decomposition, finer than the teapot's: a ring's hole must stay open for the
// next ring of the chain (13 parts, reaching 2.9 mm into it).
std::shared_ptr<const CompoundShape> ringShape() {
    static std::shared_ptr<const CompoundShape> shape = [] {
        const TriMesh ring = primitives::torus(0.05f, 0.012f);
        return std::make_shared<const CompoundShape>(convexDecomposition({ring}), ring);
    }();
    return shape;
}

} // namespace rf
