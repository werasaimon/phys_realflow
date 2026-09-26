#pragma once
// Row-major 4x4 matrix (float): affine transforms (translation, rotation, scale) and projections.
// Points are column vectors: p' = M * (p, 1).

#include "math/Quaternion.h"
#include "math/Vector4.h"

namespace rf {

struct Matrix4x4 {
    float m[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

    static Matrix4x4 identity() { return Matrix4x4(); }
    static Matrix4x4 zero() { Matrix4x4 r; for (auto& row : r.m) for (float& v : row) v = 0; return r; }
    static Matrix4x4 translation(const Vector3& t) {
        Matrix4x4 r;
        r.m[0][3] = t.x; r.m[1][3] = t.y; r.m[2][3] = t.z;
        return r;
    }
    static Matrix4x4 scale(const Vector3& s) {
        Matrix4x4 r;
        r.m[0][0] = s.x; r.m[1][1] = s.y; r.m[2][2] = s.z;
        return r;
    }
    static Matrix4x4 rotation(const Matrix3x3& R) {
        Matrix4x4 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = R.m[i][j];
        return r;
    }
    static Matrix4x4 rotation(const Quaternion& q) { return rotation(q.toMatrix3x3()); }
    // Translation * Rotation * Scale: scale first, then rotate, then move.
    static Matrix4x4 trs(const Vector3& t, const Quaternion& q, const Vector3& s) {
        Matrix3x3 R = q.toMatrix3x3();
        Matrix4x4 r;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) r.m[i][j] = R.m[i][j] * s[j];
            r.m[i][3] = t[i];
        }
        return r;
    }
    // OpenGL-style perspective projection (right-handed, camera looks along -Z, depth to [-1, 1]).
    static Matrix4x4 perspective(float fovY, float aspect, float zNear, float zFar) {
        float f = 1.0f / std::tan(0.5f * fovY);
        Matrix4x4 r = zero();
        r.m[0][0] = f / aspect;
        r.m[1][1] = f;
        r.m[2][2] = (zFar + zNear) / (zNear - zFar);
        r.m[2][3] = 2.0f * zFar * zNear / (zNear - zFar);
        r.m[3][2] = -1.0f;
        return r;
    }
    // View matrix of a camera at `eye` looking at `target`.
    static Matrix4x4 lookAt(const Vector3& eye, const Vector3& target, const Vector3& up) {
        Vector3 f = normalize(target - eye), s = normalize(cross(f, up)), u = cross(s, f);
        Matrix4x4 r;
        for (int j = 0; j < 3; ++j) { r.m[0][j] = s[j]; r.m[1][j] = u[j]; r.m[2][j] = -f[j]; }
        r.m[0][3] = -dot(s, eye);
        r.m[1][3] = -dot(u, eye);
        r.m[2][3] = dot(f, eye);
        return r;
    }

    Matrix4x4 operator*(const Matrix4x4& b) const {
        Matrix4x4 r = zero();
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                for (int k = 0; k < 4; ++k) r.m[i][j] += m[i][k] * b.m[k][j];
        return r;
    }
    Vector4 operator*(const Vector4& v) const {
        Vector4 r;
        for (int i = 0; i < 4; ++i) r[i] = m[i][0] * v.x + m[i][1] * v.y + m[i][2] * v.z + m[i][3] * v.w;
        return r;
    }
    // Point (w = 1), with the perspective divide when the last row is not (0, 0, 0, 1).
    Vector3 transformPoint(const Vector3& p) const {
        Vector4 r = *this * Vector4(p, 1.0f);
        return r.w != 0 && r.w != 1 ? r.xyz() / r.w : r.xyz();
    }
    // Direction (w = 0): rotation and scale only.
    Vector3 transformDirection(const Vector3& d) const { return (*this * Vector4(d, 0.0f)).xyz(); }
    Matrix3x3 upper3x3() const {
        Matrix3x3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = m[i][j];
        return r;
    }
    Matrix4x4 transposed() const {
        Matrix4x4 r;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) r.m[i][j] = m[j][i];
        return r;
    }
    // General inverse by cofactors (Laplace expansion over 2x2 minors); zero matrix if singular.
    Matrix4x4 inverse() const {
        const auto& a = m;
        float s0 = a[0][0] * a[1][1] - a[1][0] * a[0][1], s1 = a[0][0] * a[1][2] - a[1][0] * a[0][2];
        float s2 = a[0][0] * a[1][3] - a[1][0] * a[0][3], s3 = a[0][1] * a[1][2] - a[1][1] * a[0][2];
        float s4 = a[0][1] * a[1][3] - a[1][1] * a[0][3], s5 = a[0][2] * a[1][3] - a[1][2] * a[0][3];
        float c5 = a[2][2] * a[3][3] - a[3][2] * a[2][3], c4 = a[2][1] * a[3][3] - a[3][1] * a[2][3];
        float c3 = a[2][1] * a[3][2] - a[3][1] * a[2][2], c2 = a[2][0] * a[3][3] - a[3][0] * a[2][3];
        float c1 = a[2][0] * a[3][2] - a[3][0] * a[2][2], c0 = a[2][0] * a[3][1] - a[3][0] * a[2][1];
        float det = s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0;
        if (std::fabs(det) < 1e-30f) return zero();
        float id = 1.0f / det;
        Matrix4x4 r;
        r.m[0][0] = (a[1][1] * c5 - a[1][2] * c4 + a[1][3] * c3) * id;
        r.m[0][1] = (-a[0][1] * c5 + a[0][2] * c4 - a[0][3] * c3) * id;
        r.m[0][2] = (a[3][1] * s5 - a[3][2] * s4 + a[3][3] * s3) * id;
        r.m[0][3] = (-a[2][1] * s5 + a[2][2] * s4 - a[2][3] * s3) * id;
        r.m[1][0] = (-a[1][0] * c5 + a[1][2] * c2 - a[1][3] * c1) * id;
        r.m[1][1] = (a[0][0] * c5 - a[0][2] * c2 + a[0][3] * c1) * id;
        r.m[1][2] = (-a[3][0] * s5 + a[3][2] * s2 - a[3][3] * s1) * id;
        r.m[1][3] = (a[2][0] * s5 - a[2][2] * s2 + a[2][3] * s1) * id;
        r.m[2][0] = (a[1][0] * c4 - a[1][1] * c2 + a[1][3] * c0) * id;
        r.m[2][1] = (-a[0][0] * c4 + a[0][1] * c2 - a[0][3] * c0) * id;
        r.m[2][2] = (a[3][0] * s4 - a[3][1] * s2 + a[3][3] * s0) * id;
        r.m[2][3] = (-a[2][0] * s4 + a[2][1] * s2 - a[2][3] * s0) * id;
        r.m[3][0] = (-a[1][0] * c3 + a[1][1] * c1 - a[1][2] * c0) * id;
        r.m[3][1] = (a[0][0] * c3 - a[0][1] * c1 + a[0][2] * c0) * id;
        r.m[3][2] = (-a[3][0] * s3 + a[3][1] * s1 - a[3][2] * s0) * id;
        r.m[3][3] = (a[2][0] * s3 - a[2][1] * s1 + a[2][2] * s0) * id;
        return r;
    }
};

} // namespace rf
