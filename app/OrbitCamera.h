#pragma once
// Orbit camera around a target point + screen-to-world ray unprojection.

#include "math/Math.h"

#include <QMatrix4x4>
#include <QPointF>
#include <QSize>
#include <QVector3D>

struct Ray {
    rf::Vector3 origin;
    rf::Vector3 dir; // unit
    rf::Vector3 at(float t) const { return origin + dir * t; }
};

class OrbitCamera {
public:
    QMatrix4x4 view() const;
    QMatrix4x4 projection(float aspect) const;
    QVector3D eye() const;
    QVector3D forward() const { return (target_ - eye()).normalized(); }
    float fovDeg() const { return fov_; }
    float distance() const { return distance_; }

    void rotate(float dxPixels, float dyPixels);
    void pan(float dxPixels, float dyPixels);
    void zoom(float wheelSteps);
    // Fit the whole box in view; yaw/pitch in degrees.
    void frame(const rf::AABB& box, float yawDeg, float pitchDeg);

    // Ray through a pixel (logical coordinates, origin top-left): the inverse view-projection
    // maps the pixel at NDC depth -1 (near plane) and +1 (far plane) back to world space.
    Ray screenRay(const QPointF& pixel, const QSize& viewport) const;

private:
    QVector3D target_{0, 0.5f, 0};
    float distance_ = 3.5f, yaw_ = 55.0f, pitch_ = 22.0f, fov_ = 40.0f;
};

// Ray / axis-aligned plane intersection: plane {p : p[axis] == coord}. Returns t or -1.
float intersectAxisPlane(const Ray& r, int axis, float coord);
// Ray / box: entry and exit distances (tNear clamped to 0). False if missed.
bool intersectBox(const Ray& r, const rf::AABB& b, float& tNear, float& tFar);
