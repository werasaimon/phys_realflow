#pragma once
// 3D vector (float): positions, velocities, forces.

#include "math/Scalar.h"

namespace rf {

struct Vector3 {
    float x = 0, y = 0, z = 0;

    constexpr Vector3() = default;
    constexpr Vector3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
    explicit constexpr Vector3(float s) : x(s), y(s), z(s) {}

    float& operator[](int i) { return (&x)[i]; }
    float operator[](int i) const { return (&x)[i]; }

    constexpr Vector3 operator+(const Vector3& b) const { return {x + b.x, y + b.y, z + b.z}; }
    constexpr Vector3 operator-(const Vector3& b) const { return {x - b.x, y - b.y, z - b.z}; }
    constexpr Vector3 operator*(const Vector3& b) const { return {x * b.x, y * b.y, z * b.z}; }
    constexpr Vector3 operator/(const Vector3& b) const { return {x / b.x, y / b.y, z / b.z}; }
    constexpr Vector3 operator*(float s) const { return {x * s, y * s, z * s}; }
    constexpr Vector3 operator/(float s) const { return {x / s, y / s, z / s}; }
    constexpr Vector3 operator-() const { return {-x, -y, -z}; }
    Vector3& operator+=(const Vector3& b) { x += b.x; y += b.y; z += b.z; return *this; }
    Vector3& operator-=(const Vector3& b) { x -= b.x; y -= b.y; z -= b.z; return *this; }
    Vector3& operator*=(float s) { x *= s; y *= s; z *= s; return *this; }
    Vector3& operator/=(float s) { x /= s; y /= s; z /= s; return *this; }
};

constexpr Vector3 operator*(float s, const Vector3& v) { return v * s; }
constexpr float dot(const Vector3& a, const Vector3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
constexpr Vector3 cross(const Vector3& a, const Vector3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline float length2(const Vector3& a) { return dot(a, a); }
inline float length(const Vector3& a) { return std::sqrt(dot(a, a)); }
inline Vector3 normalize(const Vector3& a) {
    float l = length(a);
    return l > 1e-20f ? a / l : Vector3(0, 0, 0);
}
inline Vector3 vmin(const Vector3& a, const Vector3& b) { return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)}; }
inline Vector3 vmax(const Vector3& a, const Vector3& b) { return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)}; }
inline Vector3 vabs(const Vector3& a) { return {std::fabs(a.x), std::fabs(a.y), std::fabs(a.z)}; }
inline Vector3 lerp(const Vector3& a, const Vector3& b, float t) { return a + (b - a) * t; }
inline float maxComp(const Vector3& a) { return std::max(a.x, std::max(a.y, a.z)); }
inline float minComp(const Vector3& a) { return std::min(a.x, std::min(a.y, a.z)); }
// Some unit vector perpendicular to a (a need not be normalised).
inline Vector3 anyPerpendicular(const Vector3& a) {
    return normalize(std::fabs(a.x) > 0.57f ? Vector3(a.y, -a.x, 0) : Vector3(0, a.z, -a.y));
}

} // namespace rf
