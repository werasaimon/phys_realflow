// Geodesics of the Kerr metric in Hamiltonian form: the equations of motion, the RK4 and the
// adaptive Dormand-Prince RK45 integrators, the integration until capture or escape, and the
// 4-momentum of a photon arriving from a direction in an observer's frame. See Geodesic.h.
#include "relativity/Geodesic.h"

#include <algorithm>
#include <cmath>

namespace rf {

// g^{mu nu} at x and its derivatives d g^{mu nu} / d x^k for k = 1 (r) and k = 2 (theta) by the
// 5-point central difference (truncation ~ h^4 / 30 g'''''); k = 0 and 3 are zero (stationary,
// axisymmetric).
void Geodesic::contravariantDerivatives(const double x[4], double ginv[4][4], double dg[4][4][4]) const {
    metric_.contravariant(x, ginv);
    for (int k = 0; k < 4; ++k)
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) dg[k][i][j] = 0.0;
    for (int k = 1; k <= 2; ++k) {
        double h = k == 1 ? 1e-4 * std::max(std::fabs(x[1]), 1.0) : 1e-4;
        // Boyer-Lindquist coordinates are singular on the horizon: the stencil must stay outside it.
        if (k == 1 && metric_.coordinates == Metric::Coordinates::BoyerLindquist) {
            const double gap = x[1] - metric_.horizonRadius();
            if (gap > 0) h = std::min(h, 0.2 * gap);
        }
        double xp[4] = {x[0], x[1], x[2], x[3]};
        double g1[4][4], g2[4][4], g3[4][4], g4[4][4];
        xp[k] = x[k] - 2.0 * h;
        metric_.contravariant(xp, g1);
        xp[k] = x[k] - h;
        metric_.contravariant(xp, g2);
        xp[k] = x[k] + h;
        metric_.contravariant(xp, g3);
        xp[k] = x[k] + 2.0 * h;
        metric_.contravariant(xp, g4);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j) dg[k][i][j] = (g1[i][j] - 8.0 * g2[i][j] + 8.0 * g3[i][j] - g4[i][j]) / (12.0 * h);
    }
}

void Geodesic::derivative(const GeodesicState& s, double dx[4], double dp[4]) const {
    double ginv[4][4], dg[4][4][4];
    contravariantDerivatives(s.x, ginv, dg);
    for (int mu = 0; mu < 4; ++mu) {
        dx[mu] = 0.0;
        for (int nu = 0; nu < 4; ++nu) dx[mu] += ginv[mu][nu] * s.p[nu];
    }
    for (int k = 0; k < 4; ++k) {
        double sum = 0.0;
        for (int a = 0; a < 4; ++a)
            for (int b = 0; b < 4; ++b) sum += dg[k][a][b] * s.p[a] * s.p[b];
        dp[k] = -0.5 * sum;
    }
}

namespace {

struct Deriv {
    double dx[4], dp[4];
};

void evaluate(const Geodesic& g, const GeodesicState& s, Deriv& d) { g.derivative(s, d.dx, d.dp); }

GeodesicState advanced(const GeodesicState& s, double h, const Deriv& d) {
    GeodesicState out = s;
    for (int i = 0; i < 4; ++i) {
        out.x[i] += h * d.dx[i];
        out.p[i] += h * d.dp[i];
    }
    out.lambda += h;
    return out;
}

} // namespace

void Geodesic::stepRK4(GeodesicState& s, double h) const {
    Deriv k1, k2, k3, k4;
    evaluate(*this, s, k1);
    GeodesicState s2 = advanced(s, 0.5 * h, k1);
    evaluate(*this, s2, k2);
    GeodesicState s3 = advanced(s, 0.5 * h, k2);
    evaluate(*this, s3, k3);
    GeodesicState s4 = advanced(s, h, k3);
    evaluate(*this, s4, k4);
    for (int i = 0; i < 4; ++i) {
        s.x[i] += h / 6.0 * (k1.dx[i] + 2.0 * k2.dx[i] + 2.0 * k3.dx[i] + k4.dx[i]);
        s.p[i] += h / 6.0 * (k1.dp[i] + 2.0 * k2.dp[i] + 2.0 * k3.dp[i] + k4.dp[i]);
    }
    s.lambda += h;
}

// Dormand & Prince 1980 (the RK5(4)7M pair): fifth-order solution, fourth-order error estimate.
int Geodesic::stepRK45(GeodesicState& s, double& h, double tol, double hMax) const {
    static const double c2 = 1.0 / 5, c3 = 3.0 / 10, c4 = 4.0 / 5, c5 = 8.0 / 9;
    static const double a21 = 1.0 / 5;
    static const double a31 = 3.0 / 40, a32 = 9.0 / 40;
    static const double a41 = 44.0 / 45, a42 = -56.0 / 15, a43 = 32.0 / 9;
    static const double a51 = 19372.0 / 6561, a52 = -25360.0 / 2187, a53 = 64448.0 / 6561, a54 = -212.0 / 729;
    static const double a61 = 9017.0 / 3168, a62 = -355.0 / 33, a63 = 46732.0 / 5247, a64 = 49.0 / 176, a65 = -5103.0 / 18656;
    static const double b1 = 35.0 / 384, b3 = 500.0 / 1113, b4 = 125.0 / 192, b5 = -2187.0 / 6784, b6 = 11.0 / 84;
    static const double e1 = 71.0 / 57600, e3 = -71.0 / 16695, e4 = 71.0 / 1920, e5 = -17253.0 / 339200, e6 = 22.0 / 525, e7 = -1.0 / 40;
    (void)c2; (void)c3; (void)c4; (void)c5; // the stages' abscissae: unused, the system is autonomous
    int attempts = 0;
    Deriv k1;
    evaluate(*this, s, k1);
    for (;;) {
        ++attempts;
        const double sign = h < 0 ? -1.0 : 1.0;
        if (std::fabs(h) > hMax) h = sign * hMax;
        GeodesicState t = s;
        Deriv k2, k3, k4, k5, k6, k7;
        auto stage = [&](Deriv& out, const double* w, const Deriv* const* ks, int n) {
            GeodesicState q = s;
            for (int i = 0; i < 4; ++i) {
                for (int m = 0; m < n; ++m) {
                    q.x[i] += h * w[m] * ks[m]->dx[i];
                    q.p[i] += h * w[m] * ks[m]->dp[i];
                }
            }
            evaluate(*this, q, out);
        };
        { const double w[] = {a21}; const Deriv* ks[] = {&k1}; stage(k2, w, ks, 1); }
        { const double w[] = {a31, a32}; const Deriv* ks[] = {&k1, &k2}; stage(k3, w, ks, 2); }
        { const double w[] = {a41, a42, a43}; const Deriv* ks[] = {&k1, &k2, &k3}; stage(k4, w, ks, 3); }
        { const double w[] = {a51, a52, a53, a54}; const Deriv* ks[] = {&k1, &k2, &k3, &k4}; stage(k5, w, ks, 4); }
        { const double w[] = {a61, a62, a63, a64, a65}; const Deriv* ks[] = {&k1, &k2, &k3, &k4, &k5}; stage(k6, w, ks, 5); }
        // Fifth-order solution (the seventh stage is its derivative: FSAL).
        for (int i = 0; i < 4; ++i) {
            t.x[i] = s.x[i] + h * (b1 * k1.dx[i] + b3 * k3.dx[i] + b4 * k4.dx[i] + b5 * k5.dx[i] + b6 * k6.dx[i]);
            t.p[i] = s.p[i] + h * (b1 * k1.dp[i] + b3 * k3.dp[i] + b4 * k4.dp[i] + b5 * k5.dp[i] + b6 * k6.dp[i]);
        }
        t.lambda = s.lambda + h;
        evaluate(*this, t, k7);
        double err = 0.0;
        for (int i = 0; i < 4; ++i) {
            const double ex = h * (e1 * k1.dx[i] + e3 * k3.dx[i] + e4 * k4.dx[i] + e5 * k5.dx[i] + e6 * k6.dx[i] + e7 * k7.dx[i]);
            const double ep = h * (e1 * k1.dp[i] + e3 * k3.dp[i] + e4 * k4.dp[i] + e5 * k5.dp[i] + e6 * k6.dp[i] + e7 * k7.dp[i]);
            const double sx = tol * (1.0 + std::max(std::fabs(s.x[i]), std::fabs(t.x[i])));
            const double sp = tol * (1.0 + std::max(std::fabs(s.p[i]), std::fabs(t.p[i])));
            err = std::max(err, std::max(std::fabs(ex) / sx, std::fabs(ep) / sp));
        }
        if (err <= 1.0 || std::fabs(h) < 1e-14) {
            s = t;
            const double grow = err > 1e-300 ? std::min(5.0, 0.9 * std::pow(err, -0.2)) : 5.0;
            h *= std::max(grow, 0.2);
            if (std::fabs(h) > hMax) h = sign * hMax;
            return attempts;
        }
        h *= std::max(0.2, 0.9 * std::pow(err, -0.25));
    }
}

GeodesicResult Geodesic::integrate(GeodesicState& s, const Options& opt,
                                   const std::function<bool(const GeodesicState&, const GeodesicState&)>& stop) const {
    GeodesicResult r;
    metric_.constants(s.x, s.p, opt.mu, r.E0, r.L0, r.Q0);
    r.H0 = metric_.hamiltonian(s.x, s.p);
    const double rMin = opt.rMin > 0 ? opt.rMin : metric_.horizonRadius() * (1.0 + 1e-3);
    double h = (opt.backward ? -1.0 : 1.0) * (opt.adaptive ? opt.hInit : opt.hFixed);
    const double lambdaEnd = s.lambda + (opt.backward ? -opt.lambdaMax : opt.lambdaMax);
    while (r.end == GeodesicEnd::Running) {
        const GeodesicState before = s;
        if (opt.adaptive) stepRK45(s, h, opt.tol, opt.hMax);
        else stepRK4(s, h);
        ++r.steps;
        if (s.x[1] < rMin) r.end = GeodesicEnd::Captured;
        else if (s.x[1] > opt.rMax) r.end = GeodesicEnd::Escaped;
        else if ((opt.backward ? s.lambda <= lambdaEnd : s.lambda >= lambdaEnd)) r.end = GeodesicEnd::LambdaLimit;
        else if (r.steps >= opt.maxSteps) r.end = GeodesicEnd::StepLimit;
        else if (stop && stop(before, s)) r.end = GeodesicEnd::Stopped;
    }
    metric_.constants(s.x, s.p, opt.mu, r.E1, r.L1, r.Q1);
    r.H1 = metric_.hamiltonian(s.x, s.p);
    return r;
}

void Geodesic::momentumFromLocal(const Metric& m, const double x[4], const double dir[3], double speed,
                                 double energyLocal, double p[4]) {
    double g[4][4], ginv[4][4];
    m.covariant(x, g);
    m.contravariant(x, ginv);
    // The ZAMO tetrad: alpha = 1 / sqrt(-g^{tt}), omega = -g_{t phi} / g_{phi phi}.
    const double alpha = 1.0 / std::sqrt(-ginv[0][0]);
    const double omega = -g[0][3] / g[3][3];
    const double gamma = speed < 1.0 ? 1.0 / std::sqrt(1.0 - speed * speed) : 1.0; // light: energyLocal is the energy
    const double scale = speed < 1.0 ? gamma : energyLocal;
    // Contravariant p^mu = scale (e_t + v^i e_i), v = speed * dir.
    double u[4] = {0, 0, 0, 0};
    u[0] += 1.0 / alpha;
    u[3] += omega / alpha;
    u[1] += speed * dir[0] / std::sqrt(g[1][1]);
    u[2] += speed * dir[1] / std::sqrt(g[2][2]);
    u[3] += speed * dir[2] / std::sqrt(g[3][3]);
    for (int i = 0; i < 4; ++i) {
        p[i] = 0.0;
        for (int j = 0; j < 4; ++j) p[i] += g[i][j] * u[j];
        p[i] *= scale;
    }
}

} // namespace rf
