// The math toolbox against its definitions: the Einstein convention against explicit loops (and
// its refusals of formulas that break the rule), raising / lowering / permuting / tracing, and the
// matrix functions - outer and Kronecker products, QR, SVD, the matrix exponential - against the
// identities that define them (Q^T Q = I, U S V^T = A, exp of a skew matrix = the Rodrigues
// rotation = our quaternion's).
//
// Registration (tests/Tests.h and tests/main.cpp):
//   run("math: tensors - Einstein summation vs loops, raise / lower / trace, refusals", testTensorEinstein);
//   run("math: outer, Kronecker, QR, SVD, matrix exponential", testMatrixToolbox);
#include "TestRunner.h"

#include "math/MatrixNxN.h"
#include "math/Quaternion.h"
#include "math/Tensor.h"

#include <cmath>
#include <random>

using namespace rf;

namespace {

Tensor randomTensor(int dim, const std::string& variance, std::mt19937& rng) {
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    Tensor t(dim, variance);
    for (double& v : t.data()) v = u(rng);
    return t;
}

MatrixNxN randomMatrix(int rows, int cols, std::mt19937& rng) {
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    MatrixNxN m(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j) m(i, j) = u(rng);
    return m;
}

// The message of the invalid_argument einstein() throws, or "" if it did not throw.
std::string refusal(const std::function<void()>& f) {
    try {
        f();
    } catch (const std::invalid_argument& e) {
        return e.what();
    }
    return "";
}

double maxAbsDiff(const MatrixNxN& a, const MatrixNxN& b) {
    double m = 0;
    for (int i = 0; i < a.rows(); ++i)
        for (int j = 0; j < a.cols(); ++j) m = std::max(m, std::fabs(a(i, j) - b(i, j)));
    return m;
}

} // namespace

// C^a = A^a_bc B^bc and T^a_a, by einstein() and by loops; raising then lowering gives the tensor
// back; the symmetric and antisymmetric parts add up to it; the rule refuses what it must.
void testTensorEinstein() {
    std::mt19937 rng(7);
    const int n = 3;
    const Tensor A = randomTensor(n, "^__", rng), B = randomTensor(n, "^^", rng), T = randomTensor(n, "^_", rng);
    const Tensor C = einstein("^a_bc ^bc -> ^a", A, B);
    double worst = 0;
    for (int a = 0; a < n; ++a) {
        double s = 0;
        for (int b = 0; b < n; ++b)
            for (int c = 0; c < n; ++c) s += A(a, b, c) * B(b, c);
        worst = std::max(worst, std::fabs(C(a) - s));
    }
    double tr = 0;
    for (int a = 0; a < n; ++a) tr += T(a, a);
    const double trE = einstein("^a_a ->", T).scalar(), trT = trace(T, 0, 1).scalar();
    // A metric (symmetric, invertible) to raise and lower with: g = M M^T + I (a plain matrix
    // product: sum over c of M_ac M_bc is not an Einstein contraction, both c are down).
    const MatrixNxN Mm = randomMatrix(n, n, rng), gm = Mm * Mm.transposed() + MatrixNxN::identity(n);
    Tensor g(n, "__"), ginv(n, "^^");
    MatrixNxN gi;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) g(i, j) = gm(i, j);
    gm.inverse(gi);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) ginv(i, j) = gi(i, j);
    const double roundTrip = (raise(lower(A, 0, g), 0, ginv) - A).maxAbs();
    const Tensor S = randomTensor(n, "___", rng);
    const double parts = (symmetrize(S, 1, 2) + antisymmetrize(S, 1, 2) - S).maxAbs();
    const double symmetric = (symmetrize(S, 1, 2) - permute(symmetrize(S, 1, 2), "abc", "acb")).maxAbs();
    std::printf("  A^a_bc B^bc: einstein vs loops %.1e; trace T^a_a %.6f (einstein %.6f, trace() %.6f); raise(lower(A)) - A %.1e; "
                "sym + antisym - S %.1e, sym part symmetric to %.1e; outer rank %d\n",
                worst, tr, trE, trT, roundTrip, parts, symmetric, outer(A, B).rank());
    CHECK(worst < 1e-14 && std::fabs(trE - tr) < 1e-14 && std::fabs(trT - tr) < 1e-14, "Einstein summation differs from the loops");
    CHECK(roundTrip < 1e-12 && parts < 1e-15 && symmetric < 1e-15 && outer(A, B).rank() == 5, "index operations broken");
    // The rule: an index twice up, a free index missing, a variance that does not match, an index
    // that changes position - each refused with a message.
    const Tensor v = randomTensor(n, "^", rng);
    const std::string twiceUp = refusal([&] { einstein("^a ^a ->", v, v); });
    const std::string missing = refusal([&] { einstein("^a_bc ^b -> ^a", A, v); });
    const std::string variance = refusal([&] { einstein("_a ^a ->", v, v); });
    const std::string moved = refusal([&] { einstein("^a -> _a", v); });
    std::printf("  refused: \"%s\"\n", twiceUp.c_str());
    CHECK(twiceUp.find("дважды сверху") != std::string::npos && twiceUp.find("up twice") != std::string::npos, "two uppers not refused");
    CHECK(!missing.empty() && !variance.empty() && !moved.empty(), "a bad formula was accepted: '%s' '%s' '%s'", missing.c_str(),
          variance.c_str(), moved.c_str());
}

// QR, SVD and exp against their definitions (tolerances: a few hundred times the double epsilon
// of these 3 x 3 .. 6 x 6 problems).
void testMatrixToolbox() {
    std::mt19937 rng(11);
    const MatrixNxN A = randomMatrix(5, 3, rng), B = randomMatrix(3, 3, rng), C = randomMatrix(3, 3, rng), D = randomMatrix(3, 3, rng);
    const std::vector<double> u = {1, 2, 3}, w = {4, 5};
    const MatrixNxN uv = MatrixNxN::outer(u, w);
    const double kronTrace = std::fabs(MatrixNxN::kronecker(B, C).trace() - B.trace() * C.trace());
    const double mixed = maxAbsDiff(MatrixNxN::kronecker(B, C) * MatrixNxN::kronecker(D, B), MatrixNxN::kronecker(B * D, C * B));
    MatrixNxN Q, R;
    A.qr(Q, R);
    double below = 0;
    for (int i = 0; i < R.rows(); ++i)
        for (int j = 0; j < std::min(i, R.cols()); ++j) below = std::max(below, std::fabs(R(i, j)));
    const double qrErr = maxAbsDiff(Q * R, A), qOrtho = maxAbsDiff(Q.transposed() * Q, MatrixNxN::identity(5));
    double svdErr = 0, svdOrtho = 0;
    bool sorted = true;
    for (const MatrixNxN& X : {A, A.transposed()}) {
        MatrixNxN U, V;
        std::vector<double> s;
        X.svd(U, s, V);
        MatrixNxN S(int(s.size()), int(s.size()));
        for (size_t k = 0; k < s.size(); ++k) S(int(k), int(k)) = s[k];
        svdErr = std::max(svdErr, maxAbsDiff(U * S * V.transposed(), X));
        svdOrtho = std::max({svdOrtho, maxAbsDiff(U.transposed() * U, MatrixNxN::identity(int(s.size()))),
                             maxAbsDiff(V.transposed() * V, MatrixNxN::identity(int(s.size())))});
        for (size_t k = 1; k < s.size(); ++k) sorted &= s[k] <= s[k - 1];
    }
    std::printf("  outer (2,1) %.0f; tr(B x C) - trB trC %.1e; (BxC)(DxB) - (BD)x(CB) %.1e; QR: |QR - A| %.1e, |Q^TQ - I| %.1e, below R %.1e; "
                "SVD: |USV^T - A| %.1e, orthogonality %.1e, sorted %d\n",
                uv(2, 1), kronTrace, mixed, qrErr, qOrtho, below, svdErr, svdOrtho, int(sorted));
    CHECK(uv(2, 1) == 15 && kronTrace < 1e-13 && mixed < 1e-13, "outer / Kronecker products wrong");
    CHECK(qrErr < 1e-13 && qOrtho < 1e-13 && below == 0, "QR wrong");
    CHECK(svdErr < 1e-12 && svdOrtho < 1e-12 && sorted, "SVD wrong");
    // exp of the skew matrix [w]x is the rotation by |w| about w (Rodrigues); our quaternion is
    // float, so it agrees to float precision only.
    const double wx = 0.3, wy = -0.7, wz = 1.1, angle = std::sqrt(wx * wx + wy * wy + wz * wz);
    MatrixNxN K(3, 3);
    K(0, 1) = -wz; K(0, 2) = wy; K(1, 0) = wz; K(1, 2) = -wx; K(2, 0) = -wy; K(2, 1) = wx;
    const MatrixNxN Kn = K * (1.0 / angle);
    const MatrixNxN rodrigues = MatrixNxN::identity(3) + Kn * std::sin(angle) + Kn * Kn * (1.0 - std::cos(angle));
    const MatrixNxN E = K.expm();
    const Quaternion q = Quaternion::fromAxisAngle(Vector3(float(wx), float(wy), float(wz)), float(angle));
    double quatErr = 0;
    for (int j = 0; j < 3; ++j) {
        Vector3 e(0.0f);
        e[j] = 1.0f;
        const Vector3 r = q.rotate(e);
        for (int i = 0; i < 3; ++i) quatErr = std::max(quatErr, std::fabs(double(r[i]) - E(i, j)));
    }
    MatrixNxN big(3, 3), nil(2, 2);
    big(0, 0) = 10; big(1, 1) = -4; big(2, 2) = 2; // |A| > theta_13: the scaling-and-squaring path
    nil(0, 1) = 1;                                 // nilpotent: exp = I + N exactly
    const MatrixNxN Eb = big.expm(), En = nil.expm();
    const double bigErr = std::max({std::fabs(Eb(0, 0) / std::exp(10.0) - 1), std::fabs(Eb(1, 1) / std::exp(-4.0) - 1), std::fabs(Eb(2, 2) / std::exp(2.0) - 1)});
    const double rotErr = maxAbsDiff(E, rodrigues), nilErr = std::fabs(En(0, 1) - 1) + std::fabs(En(0, 0) - 1) + std::fabs(En(1, 0));
    std::printf("  expm: skew -> Rodrigues %.1e, -> our float quaternion %.1e; diag(10, -4, 2) relative %.1e; nilpotent %.1e\n",
                rotErr, quatErr, bigErr, nilErr);
    CHECK(rotErr < 1e-12 && bigErr < 1e-13 && nilErr < 1e-15, "matrix exponential wrong");
    CHECK(quatErr < 1e-6, "expm and the quaternion disagree: %e", quatErr);
}
