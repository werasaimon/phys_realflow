#pragma once
// Row-major 3x3 matrix (float): rotations, inertia tensors, effective-mass blocks.

#include "math/Vector3.h"

namespace rf {

struct Matrix3x3 {
    float m[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    static Matrix3x3 identity() { return Matrix3x3(); }
    static Matrix3x3 zero() { Matrix3x3 r; for (auto& row : r.m) for (float& v : row) v = 0; return r; }
    static Matrix3x3 diag(const Vector3& d) { Matrix3x3 r = zero(); r.m[0][0] = d.x; r.m[1][1] = d.y; r.m[2][2] = d.z; return r; }
    static Matrix3x3 fromColumns(const Vector3& a, const Vector3& b, const Vector3& c) {
        Matrix3x3 r;
        for (int i = 0; i < 3; ++i) { r.m[i][0] = a[i]; r.m[i][1] = b[i]; r.m[i][2] = c[i]; }
        return r;
    }
    static Matrix3x3 fromRows(const Vector3& a, const Vector3& b, const Vector3& c) {
        Matrix3x3 r;
        for (int j = 0; j < 3; ++j) { r.m[0][j] = a[j]; r.m[1][j] = b[j]; r.m[2][j] = c[j]; }
        return r;
    }
    // a b^T
    static Matrix3x3 outer(const Vector3& a, const Vector3& b) {
        Matrix3x3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = a[i] * b[j];
        return r;
    }
    // Cross-product matrix: skew(a) * v = cross(a, v).
    static Matrix3x3 skew(const Vector3& a) {
        Matrix3x3 r = zero();
        r.m[0][1] = -a.z; r.m[0][2] = a.y;
        r.m[1][0] = a.z;  r.m[1][2] = -a.x;
        r.m[2][0] = -a.y; r.m[2][1] = a.x;
        return r;
    }

    Vector3 operator*(const Vector3& v) const {
        return {m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
                m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
                m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
    }
    Matrix3x3 operator*(const Matrix3x3& b) const {
        Matrix3x3 r = zero();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k) r.m[i][j] += m[i][k] * b.m[k][j];
        return r;
    }
    Matrix3x3 operator*(float s) const {
        Matrix3x3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = m[i][j] * s;
        return r;
    }
    Matrix3x3 operator+(const Matrix3x3& b) const {
        Matrix3x3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = m[i][j] + b.m[i][j];
        return r;
    }
    Matrix3x3 operator-(const Matrix3x3& b) const {
        Matrix3x3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = m[i][j] - b.m[i][j];
        return r;
    }
    Matrix3x3& operator+=(const Matrix3x3& b) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) m[i][j] += b.m[i][j];
        return *this;
    }
    Matrix3x3 transposed() const {
        Matrix3x3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = m[j][i];
        return r;
    }
    Vector3 col(int j) const { return {m[0][j], m[1][j], m[2][j]}; }
    Vector3 row(int i) const { return {m[i][0], m[i][1], m[i][2]}; }
    float trace() const { return m[0][0] + m[1][1] + m[2][2]; }
    float determinant() const {
        const auto& a = m;
        return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
               a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    }
    // Inverse by the adjugate; the zero matrix when |det| <= singularDet (e.g. two static bodies).
    Matrix3x3 inverse(float singularDet = 1e-20f) const {
        const auto& a = m;
        float det = determinant();
        Matrix3x3 r = zero();
        if (std::fabs(det) <= singularDet) return r;
        float id = 1.0f / det;
        r.m[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * id;
        r.m[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * id;
        r.m[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * id;
        r.m[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * id;
        r.m[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * id;
        r.m[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * id;
        r.m[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * id;
        r.m[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * id;
        r.m[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * id;
        return r;
    }
};

// Eigen-decomposition of a symmetric 3x3 matrix (cyclic Jacobi, double precision internally):
// A = V * diag(eigenvalues) * V^T, columns of V are the eigenvectors.
inline void symmetricEigen(const Matrix3x3& A, Vector3& eig, Matrix3x3& V) {
    double a[3][3], v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    double scale = 0; // |A|^2: the stop tests are relative (unit-density inertia of a mm part is ~1e-13)
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            a[i][j] = A.m[i][j];
            scale += a[i][j] * a[i][j];
        }
    for (int sweep = 0; sweep < 50; ++sweep) {
        double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if (off <= 1e-24 * scale) break;
        const int pairs[3][2] = {{0, 1}, {0, 2}, {1, 2}};
        for (auto& pq : pairs) {
            int p = pq[0], q = pq[1];
            if (a[p][q] * a[p][q] <= 1e-30 * scale) continue;
            double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
            double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
            // a <- J^T a J,  V <- V J  with J = rotation in the (p, q) plane
            for (int k = 0; k < 3; ++k) {
                double akp = a[k][p], akq = a[k][q];
                a[k][p] = c * akp - s * akq;
                a[k][q] = s * akp + c * akq;
            }
            for (int k = 0; k < 3; ++k) {
                double apk = a[p][k], aqk = a[q][k];
                a[p][k] = c * apk - s * aqk;
                a[q][k] = s * apk + c * aqk;
            }
            for (int k = 0; k < 3; ++k) {
                double vkp = v[k][p], vkq = v[k][q];
                v[k][p] = c * vkp - s * vkq;
                v[k][q] = s * vkp + c * vkq;
            }
        }
    }
    for (int i = 0; i < 3; ++i) {
        eig[i] = float(a[i][i]);
        for (int j = 0; j < 3; ++j) V.m[i][j] = float(v[i][j]);
    }
}

} // namespace rf
