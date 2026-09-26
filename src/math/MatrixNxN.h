#pragma once
// Dense rows x cols matrix (double) with the standard direct solvers:
//   solveLU        - Gaussian elimination with partial pivoting (any non-singular matrix)
//   solveCholesky  - L L^T (symmetric positive definite, e.g. mass / stiffness matrices)
//   symmetricEigen - cyclic Jacobi rotations (all eigenvalues and eigenvectors)
// plus solveSmall(): the same elimination for n <= 4 in float without heap allocation, for hot
// loops such as the contact block solver.

#include <cstddef>
#include <vector>

namespace rf {

class MatrixNxN {
public:
    MatrixNxN() = default;
    MatrixNxN(int rows, int cols, double value = 0.0);
    static MatrixNxN identity(int n);

    int rows() const { return rows_; }
    int cols() const { return cols_; }
    double& operator()(int i, int j) { return a_[size_t(i) * cols_ + j]; }
    double operator()(int i, int j) const { return a_[size_t(i) * cols_ + j]; }

    MatrixNxN operator+(const MatrixNxN& b) const;
    MatrixNxN operator-(const MatrixNxN& b) const;
    MatrixNxN operator*(const MatrixNxN& b) const;
    MatrixNxN operator*(double s) const;
    std::vector<double> operator*(const std::vector<double>& x) const;
    MatrixNxN transposed() const;

    // A x = b. False if the matrix is singular (or not square / wrong sizes).
    bool solveLU(const std::vector<double>& b, std::vector<double>& x) const;
    double determinant() const;
    bool inverse(MatrixNxN& out) const;
    // A x = b for symmetric positive definite A. False if A is not positive definite.
    bool solveCholesky(const std::vector<double>& b, std::vector<double>& x) const;
    // Symmetric A = V diag(values) V^T; the columns of V are the eigenvectors.
    void symmetricEigen(std::vector<double>& values, MatrixNxN& vectors) const;

private:
    // In-place LU with partial pivoting: returns false if singular; perm = row order, sign = det sign.
    bool decomposeLU(MatrixNxN& lu, std::vector<int>& perm, int& sign) const;

    int rows_ = 0, cols_ = 0;
    std::vector<double> a_;
};

// M x = r for n <= 4 (float, row-major, no allocation), Gaussian elimination with partial pivoting.
bool solveSmall(int n, const float M[4][4], const float r[4], float x[4]);

} // namespace rf
