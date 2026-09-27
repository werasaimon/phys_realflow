#pragma once
// Curvature from the metric, by the textbook definitions, for any Spacetime. Conventions of Misner,
// Thorne & Wheeler "Gravitation" (MTW), signature -+++, G = c = 1:
//
//   Gamma^l_mn = 1/2 g^ls (d_m g_sn + d_n g_sm - d_s g_mn)                (MTW 8.24b)
//   R^r_smn   = d_m Gamma^r_ns - d_n Gamma^r_ms
//               + Gamma^r_ml Gamma^l_ns - Gamma^r_nl Gamma^l_ms          (MTW 11.12)
//   R_sn      = R^r_srn,   R = g^sn R_sn,   G_mn = R_mn - 1/2 R g_mn      (MTW 14.5, 14.7, 14.9)
//   K         = R_abcd R^abcd  (Kretschmann: curvature no coordinate change can hide)
//
// What each one means. Gamma says how the coordinate grid turns and stretches from point to point
// (it is not zero even in flat space in spherical coordinates). Riemann says whether a small loop
// brings a vector back turned: that, and only that, is curvature - gravity's tides. Ricci is the
// part that changes volumes; Einstein's equations G_mn + Lambda g_mn = 8 pi T_mn say that matter
// (T) is what makes it (docs/08-relativity.md, "Кривизна из метрики"; the action they come from:
// docs/11-action.md). Each tensor is computed twice, by explicit loops and by einstein(); the
// tests compare the two.
//
// Error budget of the finite differences (4th-order central, step h = 1e-3 of the coordinate's
// size, Spacetime::differenceStep), measured on Schwarzschild at r = 7M (the test "curvature:
// Schwarzschild ..." prints the whole curve, docs/08 plots it). Truncation falls as h^4: at
// h = 3e-2 / 1e-2 / 3e-3 the Kretschmann error is 1.2e-5 / 1.5e-7 / 1.2e-9. Rounding grows as
// ~1e-16 / h^2 for Riemann, which differentiates Gamma once more: 2e-9 at h = 3e-4, 1.5e-7 at 3e-5.
// They meet at h = 1e-3 with 2e-10 for the curvature and 1.6e-11 for Gamma. The tests check at
// 1e-8 .. 1e-7: a margin of 50 or more against another point, metric or machine.
#include "math/Tensor.h"
#include "relativity/Spacetime.h"

#include <vector>

namespace rf {

Tensor metricAt(const Spacetime& st, const double* x);            // g_mn, "__"
Tensor inverseMetricAt(const Spacetime& st, const double* x);     // g^mn, "^^"
Tensor metricDerivativesAt(const Spacetime& st, const double* x); // d_k g_mn, "___" (k first)

Tensor christoffel(const Spacetime& st, const double* x);         // Gamma^l_mn, "^__", loops
Tensor christoffelEinstein(const Spacetime& st, const double* x); // the same by einstein()
Tensor riemann(const Spacetime& st, const double* x);             // R^r_smn, "^___", loops
Tensor riemannEinstein(const Spacetime& st, const double* x);     // the same by einstein()
Tensor ricci(const Tensor& riemann);                              // R_sn, "__", loops
Tensor ricciEinstein(const Tensor& riemann);                      // the same by einstein()
double ricciScalar(const Tensor& ricci, const Tensor& ginv);      // R = g^sn R_sn
Tensor einsteinTensor(const Tensor& ricci, double R, const Tensor& g); // G_mn = R_mn - R g_mn / 2
double kretschmann(const Tensor& riemann, const Tensor& g, const Tensor& ginv);        // loops
double kretschmannEinstein(const Tensor& riemann, const Tensor& g, const Tensor& ginv); // einstein()

// Everything at one point.
struct CurvatureAt {
    Tensor g, ginv, gamma, riemann, ricci, einstein;
    double R = 0, K = 0;
};
CurvatureAt curvatureAt(const Spacetime& st, const double* x);

// The geodesic equation in Lagrangian form, d^2 x^l / dtau^2 = -Gamma^l_mn u^m u^n (MTW 13.37):
// works for ANY Spacetime, also a time-dependent one (the Hamiltonian Geodesic.h needs a
// stationary, axisymmetric metric). x and u have dimension() components.
void geodesicAcceleration(const Spacetime& st, const double* x, const double* u, double* a);
void geodesicStepRK4(const Spacetime& st, std::vector<double>& x, std::vector<double>& u, double h);
// Geodesic deviation (MTW 11.10): the relative acceleration of two neighbouring free fallers
// separated by xi, D^2 xi^a / dtau^2 = -R^a_bcd u^b xi^c u^d - tides.
std::vector<double> geodesicDeviation(const Tensor& riemann, const std::vector<double>& u, const std::vector<double>& xi);

} // namespace rf
