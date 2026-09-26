#pragma once
// 2D vector (float).

#include "math/Scalar.h"

namespace rf {

struct Vector2 {
    float x = 0, y = 0;

    constexpr Vector2() = default;
    constexpr Vector2(float x_, float y_) : x(x_), y(y_) {}
    explicit constexpr Vector2(float s) : x(s), y(s) {}

    float& operator[](int i) { return (&x)[i]; }
    float operator[](int i) const { return (&x)[i]; }

    constexpr Vector2 operator+(const Vector2& b) const { return {x + b.x, y + b.y}; }
    constexpr Vector2 operator-(const Vector2& b) const { return {x - b.x, y - b.y}; }
    constexpr Vector2 operator*(const Vector2& b) const { return {x * b.x, y * b.y}; }
    constexpr Vector2 operator/(const Vector2& b) const { return {x / b.x, y / b.y}; }
    constexpr Vector2 operator*(float s) const { return {x * s, y * s}; }
    constexpr Vector2 operator/(float s) const { return {x / s, y / s}; }
    constexpr Vector2 operator-() const { return {-x, -y}; }
    Vector2& operator+=(const Vector2& b) { x += b.x; y += b.y; return *this; }
    Vector2& operator-=(const Vector2& b) { x -= b.x; y -= b.y; return *this; }
    Vector2& operator*=(float s) { x *= s; y *= s; return *this; }
    Vector2& operator/=(float s) { x /= s; y /= s; return *this; }
};

constexpr Vector2 operator*(float s, const Vector2& v) { return v * s; }
constexpr float dot(const Vector2& a, const Vector2& b) { return a.x * b.x + a.y * b.y; }
// z component of the 3D cross product (signed parallelogram area).
constexpr float cross(const Vector2& a, const Vector2& b) { return a.x * b.y - a.y * b.x; }
// Counter-clockwise perpendicular.
constexpr Vector2 perp(const Vector2& a) { return {-a.y, a.x}; }
inline float length2(const Vector2& a) { return dot(a, a); }
inline float length(const Vector2& a) { return std::sqrt(dot(a, a)); }
inline Vector2 normalize(const Vector2& a) {
    float l = length(a);
    return l > 1e-20f ? a / l : Vector2(0, 0);
}
inline Vector2 vmin(const Vector2& a, const Vector2& b) { return {std::min(a.x, b.x), std::min(a.y, b.y)}; }
inline Vector2 vmax(const Vector2& a, const Vector2& b) { return {std::max(a.x, b.x), std::max(a.y, b.y)}; }
inline Vector2 lerp(const Vector2& a, const Vector2& b, float t) { return a + (b - a) * t; }

} // namespace rf
