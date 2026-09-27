// The shape of a scene-graph entity as a triangle mesh: the primitives (box, sphere, cylinder,
// cone, plane) built from the entity's size, and models read from OBJ/STL files. GraphScene builds
// the bodies from these meshes and an editor draws the same meshes, so what is seen is what is
// simulated. Model files and their convex decompositions are cached: a scene is rebuilt at every
// edit and on every Play, and reading a model or decomposing it (seconds) must happen only once.
#include "scene/SceneGraph.h"

#include "rigid/ConvexDecomposition.h"
#include "rigid/Shapes.h"

#include <cstdio>
#include <map>
#include <mutex>

namespace rf {

namespace {

// A relative model path is found from the scene file's folder; an absolute one is used as it is.
std::string resolvePath(const std::string& file, const std::string& baseDirectory) {
    const bool absolute = !file.empty() && (file[0] == '/' || file[0] == '\\' || (file.size() > 1 && file[1] == ':'));
    if (absolute || baseDirectory.empty()) return file;
    const char last = baseDirectory.back();
    return baseDirectory + (last == '/' || last == '\\' ? "" : "/") + file;
}

// The editor (drawing) and the simulation thread (building) read the caches at the same time.
std::mutex& cacheMutex() {
    static std::mutex m;
    return m;
}

// A model file as read from disk, once per path. An unreadable file is remembered as its error.
struct ModelFile {
    std::shared_ptr<const TriMesh> mesh;
    std::string error;
};

const ModelFile& modelFile(const std::string& path) {
    static std::map<std::string, ModelFile> cache;
    std::lock_guard<std::mutex> lock(cacheMutex());
    auto it = cache.find(path);
    if (it != cache.end()) return it->second;
    ModelFile f;
    auto mesh = std::make_shared<TriMesh>();
    if (loadMesh(path, *mesh, f.error) && !mesh->empty()) f.mesh = mesh;
    else if (f.error.empty()) f.error = "the file holds no triangles";
    return cache.emplace(path, f).first->second;
}

// The primitive shapes, centred on the origin: size is the full size, round shapes take the
// diameter from size.x and the height (along y) from size.y.
TriMesh primitiveMesh(ShapeKind shape, const Vector3& s) {
    TriMesh m;
    switch (shape) {
    case ShapeKind::Box: m = primitives::box(s * 0.5f); break;
    case ShapeKind::Plane: m = primitives::box(Vector3(0.5f * s.x, 0.01f, 0.5f * s.z)); break;
    case ShapeKind::Sphere: m = primitives::sphere(0.5f * s.x, 24, 12); break;
    case ShapeKind::Cylinder: // the primitive's axis is z: turn it to y
        m = primitives::cylinder(0.5f * s.x, s.y, 24);
        m.transform(Quaternion::fromAxisAngle({1, 0, 0}, -0.5f * kPi).toMatrix3x3(), Vector3(1.0f), Vector3(0.0f));
        break;
    case ShapeKind::Cone: // the primitive's apex points to -x: turn it up (+y)
        m = primitives::cone(0.5f * s.x, s.y, 24);
        m.transform(Quaternion::fromAxisAngle({0, 0, 1}, -0.5f * kPi).toMatrix3x3(), Vector3(1.0f), Vector3(0.0f));
        break;
    case ShapeKind::Mesh: break; // read from a file: see entityLocalMesh
    }
    m.translate(-m.bounds().center()); // centred on the entity's position
    return m;
}

} // namespace

Quaternion entityRotation(const Entity& e) {
    const Vector3& d = e.rotationDeg;
    return Quaternion::fromAxisAngle({0, 1, 0}, degToRad(d.y)) * Quaternion::fromAxisAngle({1, 0, 0}, degToRad(d.x)) *
           Quaternion::fromAxisAngle({0, 0, 1}, degToRad(d.z));
}

TriMesh entityLocalMesh(const Entity& e, const std::string& baseDirectory, std::string* error) {
    if (e.shape != ShapeKind::Mesh) return primitiveMesh(e.shape, e.size);
    const ModelFile& file = modelFile(resolvePath(e.meshFile, baseDirectory));
    if (!file.mesh) {
        if (error) *error = e.meshFile + ": " + file.error;
        return TriMesh();
    }
    TriMesh m = *file.mesh;
    m.fitTo(Vector3(0.0f), maxComp(e.size)); // uniform scale: the model keeps its proportions
    return m;
}

TriMesh entityMesh(const Entity& e, const std::string& baseDirectory) {
    TriMesh m = entityLocalMesh(e, baseDirectory);
    m.transform(entityRotation(e).toMatrix3x3(), Vector3(1.0f), e.position);
    return m;
}

bool entityIsGeometryOnly(const Entity& e) {
    return !e.rigid.enabled && !e.soft.enabled && !e.liquid.enabled && !e.cloth.enabled && !e.magnet.enabled &&
           !e.emitter.enabled && !e.heat.enabled && !e.flammable.enabled;
}

std::shared_ptr<const CompoundShape> entityCompound(const Entity& e, const std::string& baseDirectory) {
    const TriMesh local = entityLocalMesh(e, baseDirectory);
    if (local.empty()) return nullptr;
    // One decomposition per file and size: the voxel resolution is relative to the model, but the
    // hulls are in metres, so a resized model needs its own.
    char key[64];
    std::snprintf(key, sizeof(key), "|%.9g", double(maxComp(e.size)));
    const std::string id = resolvePath(e.meshFile, baseDirectory) + key;
    static std::map<std::string, std::shared_ptr<const CompoundShape>> cache;
    {
        std::lock_guard<std::mutex> lock(cacheMutex());
        auto it = cache.find(id);
        if (it != cache.end()) return it->second;
    }
    DecompositionParams params; // as for the teapot and the bunny of the samples
    params.resolution = 24;
    params.maxHullVertices = 64;
    auto shape = std::make_shared<const CompoundShape>(convexDecomposition({local}, params), local);
    std::lock_guard<std::mutex> lock(cacheMutex());
    return cache.emplace(id, shape).first->second;
}

} // namespace rf
