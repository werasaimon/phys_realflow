#include "relativity/Metric.h"

#include <algorithm>

namespace rf {

void Metric::covariant(const double x[4], double g[4][4]) const {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) g[i][j] = 0.0;
    const double r = x[1], th = x[2];
    const double s = std::sin(th), c = std::cos(th), s2 = s * s;
    if (coordinates == Coordinates::EddingtonFinkelstein) {
        // Schwarzschild in ingoing Eddington-Finkelstein coordinates (v, r, theta, phi).
        g[0][0] = -(1.0 - 2.0 * M / r);
        g[0][1] = g[1][0] = 1.0;
        g[2][2] = r * r;
        g[3][3] = r * r * s2;
        return;
    }
    const double Sigma = r * r + a * a * c * c;
    const double Delta = r * r - 2.0 * M * r + a * a;
    g[0][0] = -(1.0 - 2.0 * M * r / Sigma);
    g[0][3] = g[3][0] = -2.0 * M * a * r * s2 / Sigma;
    g[1][1] = Sigma / Delta;
    g[2][2] = Sigma;
    g[3][3] = (r * r + a * a + 2.0 * M * a * a * r * s2 / Sigma) * s2;
}

void Metric::contravariant(const double x[4], double ginv[4][4]) const {
    double g[4][4];
    covariant(x, g);
    invert4(g, ginv);
}

double Metric::photonSphereRadius(bool prograde) const {
    const double sgn = prograde ? -1.0 : 1.0;
    return 2.0 * M * (1.0 + std::cos(2.0 / 3.0 * std::acos(sgn * a / M)));
}

double Metric::iscoRadius(bool prograde) const {
    const double q = a / M;
    const double z1 = 1.0 + std::cbrt(1.0 - q * q) * (std::cbrt(1.0 + q) + std::cbrt(1.0 - q));
    const double z2 = std::sqrt(3.0 * q * q + z1 * z1);
    const double root = std::sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2));
    return M * (3.0 + z2 + (prograde ? -root : root));
}

void Metric::circularOrbit(double r, double& E, double& L, bool prograde) const {
    // Bardeen, Press & Teukolsky 1972, eq. 2.12-2.13 (the upper sign is prograde); a retrograde
    // orbit is the prograde one about a hole of spin -a with the sign of L flipped.
    const double sgn = prograde ? 1.0 : -1.0, aa = sgn * a;
    const double sM = std::sqrt(M), sr = std::sqrt(r);
    const double denom = std::pow(r, 0.75) * std::sqrt(r * sr - 3.0 * M * sr + 2.0 * aa * sM);
    E = (r * sr - 2.0 * M * sr + aa * sM) / denom;
    L = sgn * sM * (r * r - 2.0 * aa * sM * sr + aa * aa) / denom;
}

double Metric::keplerianOmega(double r, bool prograde) const {
    const double sgn = prograde ? 1.0 : -1.0;
    return sgn * std::sqrt(M) / (std::pow(r, 1.5) + sgn * a * std::sqrt(M));
}

void Metric::constants(const double x[4], const double p[4], double mu, double& E, double& L, double& Q) const {
    const double th = x[2], s = std::sin(th), c = std::cos(th);
    E = -p[0];
    L = p[3];
    Q = p[2] * p[2] + c * c * (a * a * (mu * mu - E * E) + L * L / (s * s));
}

double Metric::hamiltonian(const double x[4], const double p[4]) const {
    double ginv[4][4];
    contravariant(x, ginv);
    double h = 0.0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) h += ginv[i][j] * p[i] * p[j];
    return 0.5 * h;
}

double Metric::schwarzschildRadiusSI(double massKg) {
    const double G = 6.67430e-11, c = 299792458.0;
    return 2.0 * G * massKg / (c * c);
}

double Metric::timeUnitSI(double massKg) {
    const double G = 6.67430e-11, c = 299792458.0;
    return G * massKg / (c * c * c);
}

bool invert4(const double m[4][4], double out[4][4]) {
    double a[4][8];
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            a[i][j] = m[i][j];
            a[i][4 + j] = i == j ? 1.0 : 0.0;
        }
    }
    for (int col = 0; col < 4; ++col) {
        int pivot = col;
        for (int i = col + 1; i < 4; ++i)
            if (std::fabs(a[i][col]) > std::fabs(a[pivot][col])) pivot = i;
        if (std::fabs(a[pivot][col]) < 1e-300) return false;
        if (pivot != col)
            for (int j = 0; j < 8; ++j) std::swap(a[pivot][j], a[col][j]);
        const double inv = 1.0 / a[col][col];
        for (int j = 0; j < 8; ++j) a[col][j] *= inv;
        for (int i = 0; i < 4; ++i) {
            if (i == col) continue;
            const double f = a[i][col];
            if (f == 0.0) continue;
            for (int j = 0; j < 8; ++j) a[i][j] -= f * a[col][j];
        }
    }
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) out[i][j] = a[i][4 + j];
    return true;
}

} // namespace rf
