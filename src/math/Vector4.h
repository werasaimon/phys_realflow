#pragma once
// 4D vector (float): homogeneous coordinates, RGBA.

#include "math/Vector3.h"

namespace rf {

struct Vector4 {
    float x = 0, y = 0, z = 0, w = 0;

    constexpr Vector4() = default;
    constexpr Vector4(float x_, float y_, float z_, float w_) : x(x_), y(y_), z(z_), w(w_) {}
    constexpr Vector4(const Vector3& v, float w_) : x(v.x), y(v.y), z(v.z), w(w_) {}
    explicit constexpr Vector4(float s) : x(s), y(s), z(s), w(s) {}

    float& operator[](int i) { return (&x)[i]; }
    float operator[](int i) const { return (&x)[i]; }
    constexpr Vector3 xyz() const { return {x, y, z}; }

    constexpr Vector4 operator+(const Vector4& b) const { return {x + b.x, y + b.y, z + b.z, w + b.w}; }
    constexpr Vector4 operator-(const Vector4& b) const { return {x - b.x, y - b.y, z - b.z, w - b.w}; }
    constexpr Vector4 operator*(const Vector4& b) const { return {x * b.x, y * b.y, z * b.z, w * b.w}; }
    constexpr Vector4 operator*(float s) const { return {x * s, y * s, z * s, w * s}; }
    constexpr Vector4 operator/(float s) const { return {x / s, y / s, z / s, w / s}; }
    constexpr Vector4 operator-() const { return {-x, -y, -z, -w}; }
    Vector4& operator+=(const Vector4& b) { x += b.x; y += b.y; z += b.z; w += b.w; return *this; }
    Vector4& operator-=(const Vector4& b) { x -= b.x; y -= b.y; z -= b.z; w -= b.w; return *this; }
    Vector4& operator*=(float s) { x *= s; y *= s; z *= s; w *= s; return *this; }
};

constexpr Vector4 operator*(float s, const Vector4& v) { return v * s; }
constexpr float dot(const Vector4& a, const Vector4& b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }
inline float length2(const Vector4& a) { return dot(a, a); }
inline float length(const Vector4& a) { return std::sqrt(dot(a, a)); }
inline Vector4 normalize(const Vector4& a) {
    float l = length(a);
    return l > 1e-20f ? a / l : Vector4(0.0f);
}
inline Vector4 lerp(const Vector4& a, const Vector4& b, float t) { return a + (b - a) * t; }

} // namespace rf
