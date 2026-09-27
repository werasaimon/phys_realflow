// The linear-algebra toolbox of MatrixNxN beyond solving systems: outer and Kronecker products,
// trace and norms, the QR decomposition (Householder), the singular value decomposition
// (one-sided Jacobi) and the matrix exponential (scaling and squaring, Pade [13/13]).
//
// What each one is for, in one line:
//   outer(u, v)   - u v^T, the matrix that sends v's direction to u: the rank-one building block.
//   kronecker     - A (x) B, the matrix of "A on one factor, B on the other" (coupled systems).
//   qr            - an orthonormal basis Q of the columns and the coordinates R in it: least
//                   squares without squaring the condition number.
//   svd           - every matrix is a rotation, a stretch along axes (sigma), a rotation: rank,
//                   condition number, the closest matrix of lower rank, pseudo-inverse.
//   expm          - exp(A t) solves dx/dt = A x exactly; exp of a skew matrix is a rotation (the
//                   exponential map so(3) -> SO(3) of rigid rotation).
// Declared in MatrixNxN.h. Double precision, allocation allowed: not for the solvers' hot loops.
#include "math/MatrixNxN.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace rf {

MatrixNxN MatrixNxN::outer(const std::vector<double>& u, const std::vector<double>& v) {
    MatrixNxN r(int(u.size()), int(v.size()));
    for (size_t i = 0; i < u.size(); ++i)
        for (size_t j = 0; j < v.size(); ++j) r(int(i), int(j)) = u[i] * v[j];
    return r;
}

MatrixNxN MatrixNxN::kronecker(const MatrixNxN& A, const MatrixNxN& B) {
    MatrixNxN r(A.rows() * B.rows(), A.cols() * B.cols());
    for (int i = 0; i < A.rows(); ++i)
        for (int j = 0; j < A.cols(); ++j)
            for (int k = 0; k < B.rows(); ++k)
                for (int l = 0; l < B.cols(); ++l) r(i * B.rows() + k, j * B.cols() + l) = A(i, j) * B(k, l);
    return r;
}

double MatrixNxN::trace() const {
    double t = 0;
    for (int i = 0; i < std::min(rows_, cols_); ++i) t += (*this)(i, i);
    return t;
}

double MatrixNxN::frobeniusNorm() const {
    double s = 0;
    for (double v : a_) s += v * v;
    return std::sqrt(s);
}

double MatrixNxN::norm1() const {
    double m = 0;
    for (int j = 0; j < cols_; ++j) {
        double s = 0;
        for (int i = 0; i < rows_; ++i) s += std::fabs((*this)(i, j));
        m = std::max(m, s);
    }
    return m;
}

// ---------------------------------------------------------------------------------------------
// QR by Householder reflections: column k is mirrored onto the axis e_k by H = I - 2 v v^T / v^T v
// (the mirror whose normal is v = x - alpha e_k, alpha = -sign(x_k) |x| to avoid cancellation).
// R = H_n ... H_1 A, Q = H_1 ... H_n.
// ---------------------------------------------------------------------------------------------
void MatrixNxN::qr(MatrixNxN& Q, MatrixNxN& R) const {
    const int m = rows_, n = cols_;
    R = *this;
    Q = identity(m);
    std::vector<double> v(static_cast<size_t>(m));
    for (int k = 0; k < std::min(m - 1, n); ++k) {
        double norm = 0;
        for (int i = k; i < m; ++i) norm += R(i, k) * R(i, k);
        norm = std::sqrt(norm);
        if (norm == 0) continue;
        const double alpha = R(k, k) > 0 ? -norm : norm;
        double vv = 0;
        for (int i = k; i < m; ++i) {
            v[size_t(i)] = R(i, k) - (i == k ? alpha : 0.0);
            vv += v[size_t(i)] * v[size_t(i)];
        }
        if (vv == 0) continue;
        for (int j = 0; j < n; ++j) { // R <- H R
            double s = 0;
            for (int i = k; i < m; ++i) s += v[size_t(i)] * R(i, j);
            s *= 2.0 / vv;
            for (int i = k; i < m; ++i) R(i, j) -= s * v[size_t(i)];
        }
        for (int i = 0; i < m; ++i) { // Q <- Q H
            double s = 0;
            for (int l = k; l < m; ++l) s += Q(i, l) * v[size_t(l)];
            s *= 2.0 / vv;
            for (int l = k; l < m; ++l) Q(i, l) -= s * v[size_t(l)];
        }
        for (int i = k + 1; i < m; ++i) R(i, k) = 0.0; // exactly zero below the diagonal
    }
}

// ---------------------------------------------------------------------------------------------
// SVD by one-sided Jacobi: rotate pairs of columns of U = A (and of V = I alongside) until every
// two columns are orthogonal. Then A V = U, the column lengths are the singular values.
// ---------------------------------------------------------------------------------------------
static void orthogonalizeColumns(MatrixNxN& U, MatrixNxN& V) {
    const int m = U.rows(), n = U.cols();
    for (int sweep = 0; sweep < 60; ++sweep) {
        bool rotated = false;
        for (int p = 0; p < n - 1; ++p)
            for (int q = p + 1; q < n; ++q) {
                double alpha = 0, beta = 0, gamma = 0;
                for (int i = 0; i < m; ++i) {
                    alpha += U(i, p) * U(i, p);
                    beta += U(i, q) * U(i, q);
                    gamma += U(i, p) * U(i, q);
                }
                if (std::fabs(gamma) <= 1e-15 * std::sqrt(alpha * beta) || gamma == 0) continue;
                rotated = true;
                // The rotation that makes columns p and q orthogonal (as in the Jacobi eigen method).
                const double zeta = (beta - alpha) / (2.0 * gamma);
                const double t = (zeta >= 0 ? 1.0 : -1.0) / (std::fabs(zeta) + std::sqrt(1.0 + zeta * zeta));
                const double c = 1.0 / std::sqrt(1.0 + t * t), s = c * t;
                for (MatrixNxN* M : {&U, &V})
                    for (int i = 0; i < M->rows(); ++i) {
                        const double a = (*M)(i, p), b = (*M)(i, q);
                        (*M)(i, p) = c * a - s * b;
                        (*M)(i, q) = s * a + c * b;
                    }
            }
        if (!rotated) return;
    }
}

void MatrixNxN::svd(MatrixNxN& U, std::vector<double>& sigma, MatrixNxN& V) const {
    if (rows_ < cols_) { // A^T = V S U^T
        transposed().svd(V, sigma, U);
        return;
    }
    const int m = rows_, n = cols_;
    MatrixNxN W = *this, R = identity(n);
    orthogonalizeColumns(W, R);
    std::vector<int> order(static_cast<size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    std::vector<double> len(size_t(n), 0.0);
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) len[size_t(j)] += W(i, j) * W(i, j);
        len[size_t(j)] = std::sqrt(len[size_t(j)]);
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) { return len[size_t(a)] > len[size_t(b)]; });
    U = MatrixNxN(m, n);
    V = MatrixNxN(n, n);
    sigma.assign(size_t(n), 0.0);
    for (int k = 0; k < n; ++k) {
        const int j = order[size_t(k)];
        sigma[size_t(k)] = len[size_t(j)];
        for (int i = 0; i < m; ++i) U(i, k) = len[size_t(j)] > 0 ? W(i, j) / len[size_t(j)] : 0.0;
        for (int i = 0; i < n; ++i) V(i, k) = R(i, j);
    }
}

// ---------------------------------------------------------------------------------------------
// The matrix exponential (Higham 2005). exp(A) = (exp(A / 2^s))^(2^s): scale A down until its
// 1-norm is below theta_13 = 5.37, where the [13/13] Pade approximant r(A) = q(A)^-1 p(A) equals
// exp to double precision, then square s times. p(A) = U + V, q(A) = V - U with U odd, V even.
// ---------------------------------------------------------------------------------------------
MatrixNxN MatrixNxN::expm() const {
    static const double b[14] = {64764752532480000.0, 32382376266240000.0, 7771770303897600.0, 1187353796428800.0,
                                 129060195264000.0,   10559470521600.0,    670442572800.0,     33522128640.0,
                                 1323241920.0,        40840800.0,          960960.0,           16380.0,
                                 182.0,               1.0};
    const double theta13 = 5.371920351148152;
    const int n = rows_;
    const double norm = norm1();
    const int s = norm > theta13 ? int(std::ceil(std::log2(norm / theta13))) : 0;
    const MatrixNxN A = *this * std::ldexp(1.0, -s);
    const MatrixNxN I = identity(n), A2 = A * A, A4 = A2 * A2, A6 = A4 * A2;
    const MatrixNxN U = A * (A6 * (A6 * b[13] + A4 * b[11] + A2 * b[9]) + A6 * b[7] + A4 * b[5] + A2 * b[3] + I * b[1]);
    const MatrixNxN V = A6 * (A6 * b[12] + A4 * b[10] + A2 * b[8]) + A6 * b[6] + A4 * b[4] + A2 * b[2] + I * b[0];
    MatrixNxN qInverse;
    (V - U).inverse(qInverse); // q(A) is well conditioned for |A| <= theta_13 (Higham 2005, sec. 2)
    MatrixNxN E = qInverse * (V + U);
    for (int k = 0; k < s; ++k) E = E * E;
    return E;
}

} // namespace rf
