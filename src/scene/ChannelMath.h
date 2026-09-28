#pragma once
// The small tools the channel measurements share (Channels.cpp, ObjectChannels.cpp, FieldProbes.cpp).
//
// Sums of energy and momentum are kept in double: the energy of a scene is a sum of thousands of
// terms, and a change of one part in 10^5 (a numerical loss a test looks for) would drown in float
// rounding - the MHD energy lesson of docs/06 (a wave 4 million times weaker than the field).
#include "scene/Channels.h"

#include <cmath>
#include <string>

namespace rf {

// A vector in double, for sums.
struct Double3 {
    double x = 0, y = 0, z = 0;
    Double3() = default;
    Double3(double a, double b, double c) : x(a), y(b), z(c) {}
    explicit Double3(const Vector3& v) : x(v.x), y(v.y), z(v.z) {}
    Double3& operator+=(const Double3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Double3 operator+(const Double3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Double3 operator-(const Double3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Double3 operator*(double s) const { return {x * s, y * s, z * s}; }
    double length() const { return std::sqrt(x * x + y * y + z * z); }
};

inline double dot(const Double3& a, const Double3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Double3 cross(const Double3& a, const Double3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

// A measurement with its name card, in one line.
inline Measurement measured(const std::string& id, const char* unit, const std::string& name, ChannelGroup group, double value) {
    Measurement m;
    m.info.id = id;
    m.info.unit = unit;
    m.info.name = name;
    m.info.group = group;
    m.value = value;
    return m;
}

// The height reference of potential energy: the bottom corner of the scene box. The potential
// energy of a mass m at x in gravity g is  U = -m g . (x - floor)  -  m |g| h for g pointing down,
// h the height above the floor of the box.
Vector3 potentialFloor(const Simulation& sim);

} // namespace rf
