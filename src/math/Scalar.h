#pragma once
// Scalar constants and helpers shared by the math types.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace rf {

constexpr float kPi = 3.14159265358979323846f;
constexpr float kInf = std::numeric_limits<float>::infinity();

template <class T> inline T clampv(T v, T lo, T hi) { return v < lo ? lo : (v > hi ? hi : v); }
inline float sqr(float v) { return v * v; }
inline float lerp(float a, float b, float t) { return a + (b - a) * t; }
inline float degToRad(float d) { return d * (kPi / 180.0f); }
inline float radToDeg(float r) { return r * (180.0f / kPi); }

} // namespace rf
