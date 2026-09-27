#pragma once
// The elementary functions of the SDK: rf::sin, cos, asin, acos, atan, atan2, exp, pow.
//
// Step 1 of 2. For now each one simply calls the compiler's std:: function, so every result is
// bit for bit what it was before. What this step gives: the rigid-body simulation path (the
// quaternions, the narrow phase, EPA, the joints, fracture, the meshes and BVH that shapes are
// made of) calls these, not std::, so there is ONE place to change.
//
// Step 2 (not done yet): our own series here - the argument reduced to a small interval (Cody &
// Waite, "Software Manual for the Elementary Functions", 1980), then a Taylor or minimax
// polynomial built only from + - * / and sqrt, which IEEE 754 rounds the same way on every
// machine. Box2D v3 does this for the same reason (b2ComputeCosSin, b2Atan2; E. Catto,
// "Determinism", box2d.org, 2024): the C++ standard does not say how exact std::sin must be, so
// MinGW, glibc and MSVC round its last bit differently, and a pile of boxes amplifies one bit into
// a different pile after a few seconds. Until step 2 the results are the same for any number of
// threads and from run to run, but bit-equal across compilers only with the same C runtime.
// Gas, particles and relativity still call std:: directly (docs/10-verification.md lists them).

#include <cmath>

namespace rf {

inline float sin(float x) { return std::sin(x); }
inline double sin(double x) { return std::sin(x); }
inline float cos(float x) { return std::cos(x); }
inline double cos(double x) { return std::cos(x); }
inline float asin(float x) { return std::asin(x); }
inline double asin(double x) { return std::asin(x); }
inline float acos(float x) { return std::acos(x); }
inline double acos(double x) { return std::acos(x); }
inline float atan(float x) { return std::atan(x); }
inline double atan(double x) { return std::atan(x); }
inline float atan2(float y, float x) { return std::atan2(y, x); }
inline double atan2(double y, double x) { return std::atan2(y, x); }
inline float exp(float x) { return std::exp(x); }
inline double exp(double x) { return std::exp(x); }
inline float pow(float x, float y) { return std::pow(x, y); }
inline double pow(double x, double y) { return std::pow(x, y); }

} // namespace rf
