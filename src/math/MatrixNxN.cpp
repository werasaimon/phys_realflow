#include "math/MatrixNxN.h"

#include <cmath>
#include <utility>

namespace rf {

MatrixNxN::MatrixNxN(int rows, int cols, double value) : rows_(rows), cols_(cols), a_(size_t(rows) * cols, value) {}

MatrixNxN MatrixNxN::identity(int n) {
    MatrixNxN m(n, n);
    for (int i = 0; i < n; ++i) m(i, i) = 1.0;
    return m;
}

MatrixNxN MatrixNxN::operator+(const MatrixNxN& b) const {
    MatrixNxN r = *this;
    for (size_t k = 0; k < a_.size(); ++k) r.a_[k] += b.a_[k];
    return r;
}

MatrixNxN MatrixNxN::operator-(const MatrixNxN& b) const {
    MatrixNxN r = *this;
    for (size_t k = 0; k < a_.size(); ++k) r.a_[k] -= b.a_[k];
    return r;
}

MatrixNxN MatrixNxN::operator*(const MatrixNxN& b) const {
    MatrixNxN r(rows_, b.cols_);
    for (int i = 0; i < rows_; ++i)
        for (int k = 0; k < cols_; ++k) {
            double aik = (*this)(i, k);
            for (int j = 0; j < b.cols_; ++j) r(i, j) += aik * b(k, j);
        }
    return r;
}

MatrixNxN MatrixNxN::operator*(double s) const {
    MatrixNxN r = *this;
    for (double& v : r.a_) v *= s;
    return r;
}

std::vector<double> MatrixNxN::operator*(const std::vector<double>& x) const {
    std::vector<double> y(rows_, 0.0);
    for (int i = 0; i < rows_; ++i)
        for (int j = 0; j < cols_; ++j) y[i] += (*this)(i, j) * x[j];
    return y;
}

MatrixNxN MatrixNxN::transposed() const {
    MatrixNxN r(cols_, rows_);
    for (int i = 0; i < rows_; ++i)
        for (int j = 0; j < cols_; ++j) r(j, i) = (*this)(i, j);
    return r;
}

// ---------------------------------------------------------------------------
// LU decomposition with partial pivoting: P A = L U (L unit lower, U upper, stored together)
// ---------------------------------------------------------------------------
bool MatrixNxN::decomposeLU(MatrixNxN& lu, std::vector<int>& perm, int& sign) const {
    if (rows_ != cols_) return false;
    const int n = rows_;
    lu = *this;
    perm.resize(n);
    for (int i = 0; i < n; ++i) perm[i] = i;
    sign = 1;
    double scale = 0;
    for (double v : a_) scale = std::max(scale, std::fabs(v));
    const double tiny = 1e-14 * std::max(scale, 1e-300);
    for (int c = 0; c < n; ++c) {
        int piv = c; // largest entry in the column below the diagonal
        for (int i = c + 1; i < n; ++i)
            if (std::fabs(lu(i, c)) > std::fabs(lu(piv, c))) piv = i;
        if (std::fabs(lu(piv, c)) <= tiny) return false;
        if (piv != c) {
            for (int j = 0; j < n; ++j) std::swap(lu(c, j), lu(piv, j));
            std::swap(perm[c], perm[piv]);
            sign = -sign;
        }
        for (int i = c + 1; i < n; ++i) {
            lu(i, c) /= lu(c, c); // multiplier, stored in L
            for (int j = c + 1; j < n; ++j) lu(i, j) -= lu(i, c) * lu(c, j);
        }
    }
    return true;
}

bool MatrixNxN::solveLU(const std::vector<double>& b, std::vector<double>& x) const {
    MatrixNxN lu;
    std::vector<int> perm;
    int sign;
    if (int(b.size()) != rows_ || !decomposeLU(lu, perm, sign)) return false;
    const int n = rows_;
    x.assign(n, 0.0);
    for (int i = 0; i < n; ++i) { // forward: L y = P b
        double s = b[perm[i]];
        for (int j = 0; j < i; ++j) s -= lu(i, j) * x[j];
        x[i] = s;
    }
    for (int i = n - 1; i >= 0; --i) { // backward: U x = y
        double s = x[i];
        for (int j = i + 1; j < n; ++j) s -= lu(i, j) * x[j];
        x[i] = s / lu(i, i);
    }
    return true;
}

double MatrixNxN::determinant() const {
    MatrixNxN lu;
    std::vector<int> perm;
    int sign;
    if (!decomposeLU(lu, perm, sign)) return 0.0;
    double d = sign;
    for (int i = 0; i < rows_; ++i) d *= lu(i, i);
    return d;
}

bool MatrixNxN::inverse(MatrixNxN& out) const {
    if (rows_ != cols_) return false;
    const int n = rows_;
    out = MatrixNxN(n, n);
    std::vector<double> e(n), col;
    for (int j = 0; j < n; ++j) { // column j of the inverse solves A x = e_j
        std::fill(e.begin(), e.end(), 0.0);
        e[j] = 1.0;
        if (!solveLU(e, col)) return false;
        for (int i = 0; i < n; ++i) out(i, j) = col[i];
    }
    return true;
}

// ---------------------------------------------------------------------------
// Cholesky: A = L L^T
// ---------------------------------------------------------------------------
bool MatrixNxN::solveCholesky(const std::vector<double>& b, std::vector<double>& x) const {
    if (rows_ != cols_ || int(b.size()) != rows_) return false;
    const int n = rows_;
    MatrixNxN L(n, n);
    for (int j = 0; j < n; ++j) {
        double d = (*this)(j, j);
        for (int k = 0; k < j; ++k) d -= L(j, k) * L(j, k);
        if (d <= 0) return false; // not positive definite
        L(j, j) = std::sqrt(d);
        for (int i = j + 1; i < n; ++i) {
            double s = (*this)(i, j);
            for (int k = 0; k < j; ++k) s -= L(i, k) * L(j, k);
            L(i, j) = s / L(j, j);
        }
    }
    x = b;
    for (int i = 0; i < n; ++i) { // L y = b
        for (int k = 0; k < i; ++k) x[i] -= L(i, k) * x[k];
        x[i] /= L(i, i);
    }
    for (int i = n - 1; i >= 0; --i) { // L^T x = y
        for (int k = i + 1; k < n; ++k) x[i] -= L(k, i) * x[k];
        x[i] /= L(i, i);
    }
    return true;
}

// ---------------------------------------------------------------------------
// Cyclic Jacobi eigen-decomposition of a symmetric matrix
// ---------------------------------------------------------------------------
void MatrixNxN::symmetricEigen(std::vector<double>& values, MatrixNxN& V) const {
    const int n = rows_;
    MatrixNxN a = *this;
    V = identity(n);
    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0, total = 0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                total += a(i, j) * a(i, j);
                if (i != j) off += a(i, j) * a(i, j);
            }
        if (off <= 1e-30 * std::max(total, 1e-300)) break;
        for (int p = 0; p < n - 1; ++p)
            for (int q = p + 1; q < n; ++q) {
                if (std::fabs(a(p, q)) < 1e-300) continue;
                // Rotation in the (p, q) plane that zeroes a(p, q).
                double theta = (a(q, q) - a(p, p)) / (2.0 * a(p, q));
                double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < n; ++k) { // a <- a J
                    double akp = a(k, p), akq = a(k, q);
                    a(k, p) = c * akp - s * akq;
                    a(k, q) = s * akp + c * akq;
                }
                for (int k = 0; k < n; ++k) { // a <- J^T a
                    double apk = a(p, k), aqk = a(q, k);
                    a(p, k) = c * apk - s * aqk;
                    a(q, k) = s * apk + c * aqk;
                }
                for (int k = 0; k < n; ++k) { // V <- V J
                    double vkp = V(k, p), vkq = V(k, q);
                    V(k, p) = c * vkp - s * vkq;
                    V(k, q) = s * vkp + c * vkq;
                }
            }
    }
    values.resize(n);
    for (int i = 0; i < n; ++i) values[i] = a(i, i);
}

// ---------------------------------------------------------------------------
// Small dense solve (n <= 4, float, no allocation)
// ---------------------------------------------------------------------------
bool solveSmall(int n, const float M[4][4], const float r[4], float x[4]) {
    float a[4][5];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) a[i][j] = M[i][j];
        a[i][n] = r[i];
    }
    for (int c = 0; c < n; ++c) {
        int piv = c;
        for (int i = c + 1; i < n; ++i)
            if (std::fabs(a[i][c]) > std::fabs(a[piv][c])) piv = i;
        if (std::fabs(a[piv][c]) < 1e-12f) return false;
        if (piv != c)
            for (int j = 0; j <= n; ++j) std::swap(a[c][j], a[piv][j]);
        for (int i = c + 1; i < n; ++i) {
            float f = a[i][c] / a[c][c];
            for (int j = c; j <= n; ++j) a[i][j] -= f * a[c][j];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        float v = a[i][n];
        for (int j = i + 1; j < n; ++j) v -= a[i][j] * x[j];
        x[i] = v / a[i][i];
    }
    return true;
}

} // namespace rf
