// The static obstacle of a Simulation: one triangle mesh built from ObstacleSettings (a primitive,
// a wing, or a loaded model - a whole terrain will do), with a BVH over it. The rigid bodies, the
// particles and the gas all collide with the same mesh.
#include "scene/Simulation.h"

namespace rf {

void Simulation::buildObstacleGeometry() {
    TriMesh m;
    const float s = std::max(obstacle.size, 1e-3f);
    switch (obstacle.shape) {
    case ObstacleShape::None: break;
    case ObstacleShape::Sphere: m = primitives::sphere(0.5f * s); break;
    case ObstacleShape::Cube: m = primitives::box(Vector3(0.5f * s)); break;
    case ObstacleShape::Cylinder: m = primitives::cylinder(0.5f * s, obstacle.span); break;
    case ObstacleShape::Wing:
        m = primitives::nacaWing(obstacle.nacaCode, s, obstacle.span);
        m.translate({-0.35f * s, 0, 0}); // rotate about ~ the aerodynamic centre
        break;
    case ObstacleShape::Streamlined: m = primitives::streamlinedBody(s, obstacle.thickness); break;
    case ObstacleShape::Ellipsoid: m = primitives::ellipsoid({0.5f * s, 0.25f * s, 0.25f * s}); break;
    case ObstacleShape::Cone: m = primitives::cone(0.5f * s, s); break;
    case ObstacleShape::Custom:
        if (obstacle.customMesh) {
            m = *obstacle.customMesh;
            m.fitTo(Vector3(0.0f), s);
        }
        break;
    }
    if (!m.empty()) {
        Quaternion q = Quaternion::fromEuler(degToRad(obstacle.yawDeg), -degToRad(obstacle.angleOfAttackDeg), degToRad(obstacle.rollDeg));
        m.transform(q.toMatrix3x3(), Vector3(1.0f), obstacle.position);
    }
    obstacleMesh_ = std::make_shared<TriMesh>(std::move(m));
    obstacleBVH_.build(*obstacleMesh_);
    ++obstacleVersion_;
}

void Simulation::rebuildObstacle() {
    buildObstacleGeometry();
    reset();
    ++paramsVersion_;
}

bool Simulation::loadCustomMesh(const std::string& path, std::string& error) {
    auto m = std::make_shared<TriMesh>();
    if (!loadMesh(path, *m, error)) return false;
    obstacle.customMesh = m;
    auto slash = path.find_last_of("/\\");
    obstacle.customName = slash == std::string::npos ? path : path.substr(slash + 1);
    obstacle.shape = ObstacleShape::Custom;
    rebuildObstacle();
    return true;
}

} // namespace rf
