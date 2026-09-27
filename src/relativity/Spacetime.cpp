// The built-in spacetimes: each metric() writes g_mn exactly as the formula in Spacetime.h reads.
// Coordinates x = (x^0, x^1, ...) in the order of coordinates(); geometric units G = c = 1.
#include "relativity/Spacetime.h"

namespace rf {

namespace {

// g = diag(d0, d1, ...) in dimension n.
void diagonal(double* g, int n, const double* d) {
    for (int i = 0; i < n * n; ++i) g[i] = 0.0;
    for (int i = 0; i < n; ++i) g[i * n + i] = d[i];
}

// ds^2 = -f dt^2 + dr^2 / f + r^2 dOmega^2: every static spherical metric of this family
// (Schwarzschild, Reissner-Nordstrom, de Sitter) differs only in f(r).
void staticSpherical(double f, const double* x, double* g) {
    const double r = x[1], s = std::sin(x[2]);
    const double d[4] = {-f, 1.0 / f, r * r, r * r * s * s};
    diagonal(g, 4, d);
}

} // namespace

void MinkowskiCartesian::metric(const double*, double* g) const {
    const double d[4] = {-1.0, 1.0, 1.0, 1.0};
    diagonal(g, 4, d);
}

void MinkowskiSpherical::metric(const double* x, double* g) const { staticSpherical(1.0, x, g); }

void TwoSphere::metric(const double* x, double* g) const {
    const double s = std::sin(x[0]);
    const double d[2] = {R * R, R * R * s * s};
    diagonal(g, 2, d);
}

void SchwarzschildSpacetime::metric(const double* x, double* g) const { staticSpherical(1.0 - 2.0 * M / x[1], x, g); }

void KerrSpacetime::metric(const double* x, double* g) const {
    double g4[4][4];
    kerr.covariant(x, g4);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) g[i * 4 + j] = g4[i][j];
}

void ReissnerNordstrom::metric(const double* x, double* g) const {
    const double r = x[1];
    staticSpherical(1.0 - 2.0 * M / r + Q * Q / (r * r), x, g);
}

void DeSitterStatic::metric(const double* x, double* g) const {
    const double r = x[1];
    staticSpherical(1.0 - Lambda * r * r / 3.0, x, g);
}

void FlatFLRW::metric(const double* x, double* g) const {
    const double s = a(x[0]);
    const double d[4] = {-1.0, s * s, s * s, s * s};
    diagonal(g, 4, d);
}

void EllisWormhole::metric(const double* x, double* g) const {
    const double r2 = x[1] * x[1] + b0 * b0, s = std::sin(x[2]);
    const double d[4] = {-1.0, 1.0, r2, r2 * s * s};
    diagonal(g, 4, d);
}

} // namespace rf
