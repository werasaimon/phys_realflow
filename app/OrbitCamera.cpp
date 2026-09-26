#include "OrbitCamera.h"

#include <QtMath>

#include <algorithm>
#include <cmath>

QVector3D OrbitCamera::eye() const {
    float cy = std::cos(qDegreesToRadians(yaw_)), sy = std::sin(qDegreesToRadians(yaw_));
    float cp = std::cos(qDegreesToRadians(pitch_)), sp = std::sin(qDegreesToRadians(pitch_));
    return target_ + distance_ * QVector3D(cp * cy, sp, cp * sy);
}

QMatrix4x4 OrbitCamera::view() const {
    QMatrix4x4 v;
    v.lookAt(eye(), target_, QVector3D(0, 1, 0));
    return v;
}

QMatrix4x4 OrbitCamera::projection(float aspect) const {
    QMatrix4x4 p;
    p.perspective(fov_, aspect, std::max(0.005f, distance_ * 0.01f), distance_ * 50.0f + 50.0f);
    return p;
}

void OrbitCamera::rotate(float dx, float dy) {
    yaw_ += dx * 0.4f;
    pitch_ = std::clamp(pitch_ + dy * 0.4f, -89.0f, 89.0f);
}

void OrbitCamera::pan(float dx, float dy) {
    QMatrix4x4 v = view();
    QVector3D right(v(0, 0), v(0, 1), v(0, 2)), up(v(1, 0), v(1, 1), v(1, 2));
    float s = distance_ * 0.0015f;
    target_ += (-right * dx + up * dy) * s;
}

void OrbitCamera::zoom(float steps) {
    distance_ = std::clamp(distance_ * std::pow(0.88f, steps), 0.05f, 500.0f);
}

void OrbitCamera::frame(const rf::AABB& box, float yawDeg, float pitchDeg) {
    if (!box.valid()) return;
    rf::Vector3 c = box.center();
    target_ = QVector3D(c.x, c.y, c.z);
    float radius = 0.5f * rf::length(box.extent());
    distance_ = radius / std::sin(qDegreesToRadians(fov_ * 0.5f)) * 0.95f;
    yaw_ = yawDeg;
    pitch_ = pitchDeg;
}

Ray OrbitCamera::screenRay(const QPointF& px, const QSize& vp) const {
    const float w = std::max(1, vp.width()), h = std::max(1, vp.height());
    const float x = 2.0f * float(px.x()) / w - 1.0f;
    const float y = 1.0f - 2.0f * float(px.y()) / h;
    QMatrix4x4 inv = (projection(w / h) * view()).inverted();
    QVector3D nearP = inv.map(QVector3D(x, y, -1.0f)); // map() performs the perspective divide
    QVector3D farP = inv.map(QVector3D(x, y, 1.0f));
    QVector3D d = (farP - nearP).normalized();
    return Ray{rf::Vector3(nearP.x(), nearP.y(), nearP.z()), rf::Vector3(d.x(), d.y(), d.z())};
}

float intersectAxisPlane(const Ray& r, int axis, float coord) {
    if (std::fabs(r.dir[axis]) < 1e-8f) return -1.0f;
    float t = (coord - r.origin[axis]) / r.dir[axis];
    return t >= 0 ? t : -1.0f;
}

bool intersectBox(const Ray& r, const rf::AABB& b, float& tNear, float& tFar) {
    tNear = 0.0f;
    tFar = rf::kInf;
    for (int a = 0; a < 3; ++a) {
        if (std::fabs(r.dir[a]) < 1e-12f) {
            if (r.origin[a] < b.lo[a] || r.origin[a] > b.hi[a]) return false;
            continue;
        }
        float t0 = (b.lo[a] - r.origin[a]) / r.dir[a], t1 = (b.hi[a] - r.origin[a]) / r.dir[a];
        if (t0 > t1) std::swap(t0, t1);
        tNear = std::max(tNear, t0);
        tFar = std::min(tFar, t1);
        if (tNear > tFar) return false;
    }
    return true;
}
