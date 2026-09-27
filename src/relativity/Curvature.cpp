// Curvature from the metric: Christoffel symbols, the Riemann, Ricci and Einstein tensors, the
// Kretschmann scalar, the geodesic equation with Gamma and geodesic deviation. Every tensor twice:
// explicit loops (the formula term by term) and einstein() (the formula as text). See Curvature.h.
#include "relativity/Curvature.h"

#include "math/MatrixNxN.h"

#include <functional>

namespace rf {

namespace {

using Field = std::function<std::vector<double>(const double*)>;

// d_k F_c at x for every coordinate k and component c, by the 4th-order central difference
// (-F(x + 2h) + 8 F(x + h) - 8 F(x - h) + F(x - 2h)) / 12h; out[k * count + c].
std::vector<double> differentiate(const Spacetime& st, const double* x, const Field& F, size_t count) {
    const int n = st.dimension();
    std::vector<double> out(size_t(n) * count), xs(x, x + n);
    for (int k = 0; k < n; ++k) {
        const double h = st.differenceStep(k, x);
        auto at = [&](double shift) {
            xs[size_t(k)] = x[k] + shift;
            std::vector<double> f = F(xs.data());
            xs[size_t(k)] = x[k];
            return f;
        };
        const std::vector<double> p2 = at(2 * h), p1 = at(h), m1 = at(-h), m2 = at(-2 * h);
        for (size_t c = 0; c < count; ++c) out[size_t(k) * count + c] = (-p2[c] + 8.0 * p1[c] - 8.0 * m1[c] + m2[c]) / (12.0 * h);
    }
    return out;
}

std::vector<double> metricComponents(const Spacetime& st, const double* x) {
    std::vector<double> g(size_t(st.dimension() * st.dimension()));
    st.metric(x, g.data());
    return g;
}

// A tensor of the given variance whose components are `values` (row-major, as Tensor stores them).
Tensor fromValues(int dim, const std::string& variance, const std::vector<double>& values) {
    Tensor t(dim, variance);
    t.data() = values;
    return t;
}

} // namespace

Tensor metricAt(const Spacetime& st, const double* x) {
    return fromValues(st.dimension(), "__", metricComponents(st, x));
}

Tensor inverseMetricAt(const Spacetime& st, const double* x) {
    const int n = st.dimension();
    const Tensor g = metricAt(st, x);
    MatrixNxN m(n, n), inv;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) m(i, j) = g(i, j);
    m.inverse(inv);
    Tensor gi(n, "^^");
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) gi(i, j) = inv(i, j);
    return gi;
}

Tensor metricDerivativesAt(const Spacetime& st, const double* x) {
    const int n = st.dimension();
    std::vector<double> dg(size_t(n * n * n));
    if (!st.metricDerivatives(x, dg.data()))
        dg = differentiate(st, x, [&](const double* y) { return metricComponents(st, y); }, size_t(n * n));
    return fromValues(n, "___", dg);
}

// ---------------------------------------------------------------------------------------------
// Christoffel symbols: Gamma^l_mn = 1/2 g^ls (d_m g_sn + d_n g_sm - d_s g_mn)
// ---------------------------------------------------------------------------------------------
Tensor christoffel(const Spacetime& st, const double* x) {
    const int n = st.dimension();
    const Tensor ginv = inverseMetricAt(st, x), dg = metricDerivativesAt(st, x); // dg(k, i, j) = d_k g_ij
    Tensor gamma(n, "^__");
    for (int l = 0; l < n; ++l)
        for (int m = 0; m < n; ++m)
            for (int nu = 0; nu < n; ++nu) {
                double sum = 0;
                for (int s = 0; s < n; ++s) sum += ginv(l, s) * (dg(m, s, nu) + dg(nu, s, m) - dg(s, m, nu));
                gamma(l, m, nu) = 0.5 * sum;
            }
    return gamma;
}

Tensor christoffelEinstein(const Spacetime& st, const double* x) {
    const Tensor ginv = inverseMetricAt(st, x), dg = metricDerivativesAt(st, x);
    // The bracket B_smn = d_m g_sn + d_n g_sm - d_s g_mn: the slots of dg are (derivative, i, j).
    const Tensor B = permute(dg, "msn", "smn") + permute(dg, "nsm", "smn") - dg;
    return einstein("^ls _smn -> ^l_mn", ginv, B) * 0.5;
}

// ---------------------------------------------------------------------------------------------
// Riemann: R^r_smn = d_m Gamma^r_ns - d_n Gamma^r_ms + Gamma^r_ml Gamma^l_ns - Gamma^r_nl Gamma^l_ms
// ---------------------------------------------------------------------------------------------
// d_k Gamma^r_ab as a tensor with slots (k, r, a, b), variance "_^__".
static Tensor christoffelDerivatives(const Spacetime& st, const double* x, bool byEinstein) {
    const int n = st.dimension();
    Field gammaAt = [&](const double* y) { return (byEinstein ? christoffelEinstein(st, y) : christoffel(st, y)).data(); };
    return fromValues(n, "_^__", differentiate(st, x, gammaAt, size_t(n * n * n)));
}

Tensor riemann(const Spacetime& st, const double* x) {
    const int n = st.dimension();
    const Tensor G = christoffel(st, x), dG = christoffelDerivatives(st, x, false);
    Tensor R(n, "^___");
    for (int r = 0; r < n; ++r)
        for (int s = 0; s < n; ++s)
            for (int m = 0; m < n; ++m)
                for (int nu = 0; nu < n; ++nu) {
                    double v = dG(m, r, nu, s) - dG(nu, r, m, s); // d_m Gamma^r_ns - d_n Gamma^r_ms
                    for (int l = 0; l < n; ++l) v += G(r, m, l) * G(l, nu, s) - G(r, nu, l) * G(l, m, s);
                    R(r, s, m, nu) = v;
                }
    return R;
}

Tensor riemannEinstein(const Spacetime& st, const double* x) {
    const Tensor G = christoffelEinstein(st, x), dG = christoffelDerivatives(st, x, true);
    return permute(dG, "mrns", "rsmn") - permute(dG, "nrms", "rsmn") // d_m Gamma^r_ns - d_n Gamma^r_ms
           + einstein("^r_ml ^l_ns -> ^r_smn", G, G) - einstein("^r_nl ^l_ms -> ^r_smn", G, G);
}

// ---------------------------------------------------------------------------------------------
// Ricci, scalar curvature, Einstein, Kretschmann
// ---------------------------------------------------------------------------------------------
Tensor ricci(const Tensor& R) {
    const int n = R.dim();
    Tensor ric(n, "__");
    for (int s = 0; s < n; ++s)
        for (int nu = 0; nu < n; ++nu)
            for (int r = 0; r < n; ++r) ric(s, nu) += R(r, s, r, nu);
    return ric;
}

Tensor ricciEinstein(const Tensor& R) { return einstein("^r_srn -> _sn", R); }

double ricciScalar(const Tensor& ric, const Tensor& ginv) { return einstein("^sn _sn ->", ginv, ric).scalar(); }

Tensor einsteinTensor(const Tensor& ric, double R, const Tensor& g) { return ric - g * (0.5 * R); }

// R_abcd = g_ra R^a_bcd: the first index lowered (loops).
static Tensor loweredRiemann(const Tensor& R, const Tensor& g) {
    const int n = R.dim();
    Tensor down(n, "____");
    for (int r = 0; r < n; ++r)
        for (int s = 0; s < n; ++s)
            for (int m = 0; m < n; ++m)
                for (int nu = 0; nu < n; ++nu)
                    for (int a = 0; a < n; ++a) down(r, s, m, nu) += g(r, a) * R(a, s, m, nu);
    return down;
}

double kretschmann(const Tensor& R, const Tensor& g, const Tensor& ginv) {
    const int n = R.dim();
    const Tensor down = loweredRiemann(R, g);
    // K = sum over r s m n of R_rsmn R^rsmn, R^rsmn = g^ra g^sb g^mc g^nd R_abcd.
    double K = 0;
    for (int r = 0; r < n; ++r)
        for (int s = 0; s < n; ++s)
            for (int m = 0; m < n; ++m)
                for (int nu = 0; nu < n; ++nu) {
                    double up = 0;
                    for (int a = 0; a < n; ++a)
                        for (int b = 0; b < n; ++b)
                            for (int c = 0; c < n; ++c)
                                for (int d = 0; d < n; ++d) up += ginv(r, a) * ginv(s, b) * ginv(m, c) * ginv(nu, d) * down(a, b, c, d);
                    K += down(r, s, m, nu) * up;
                }
    return K;
}

double kretschmannEinstein(const Tensor& R, const Tensor& g, const Tensor& ginv) {
    const Tensor down = einstein("_ra ^a_smn -> _rsmn", g, R);
    // K = R_rsmn R^r_bcd g^sb g^mc g^nd: the other three indices raised by the inverse metric.
    return einstein("_rsmn ^r_bcd ^sb ^mc ^nd ->", std::vector<const Tensor*>{&down, &R, &ginv, &ginv, &ginv}).scalar();
}

CurvatureAt curvatureAt(const Spacetime& st, const double* x) {
    CurvatureAt c;
    c.g = metricAt(st, x);
    c.ginv = inverseMetricAt(st, x);
    c.gamma = christoffel(st, x);
    c.riemann = riemann(st, x);
    c.ricci = ricci(c.riemann);
    c.R = ricciScalar(c.ricci, c.ginv);
    c.einstein = einsteinTensor(c.ricci, c.R, c.g);
    c.K = kretschmann(c.riemann, c.g, c.ginv);
    return c;
}

// ---------------------------------------------------------------------------------------------
// Geodesics with Gamma, geodesic deviation
// ---------------------------------------------------------------------------------------------
void geodesicAcceleration(const Spacetime& st, const double* x, const double* u, double* a) {
    const int n = st.dimension();
    const Tensor G = christoffel(st, x);
    for (int l = 0; l < n; ++l) {
        double sum = 0;
        for (int m = 0; m < n; ++m)
            for (int nu = 0; nu < n; ++nu) sum += G(l, m, nu) * u[m] * u[nu];
        a[l] = -sum;
    }
}

// Classical RK4 on the first-order system dx/dtau = u, du/dtau = -Gamma u u.
void geodesicStepRK4(const Spacetime& st, std::vector<double>& x, std::vector<double>& u, double h) {
    const size_t n = x.size();
    std::vector<double> kx[4], ku[4], xs(n), us(n);
    const double w[4] = {0.0, 0.5, 0.5, 1.0};
    for (int stage = 0; stage < 4; ++stage) {
        for (size_t i = 0; i < n; ++i) {
            xs[i] = x[i] + (stage ? w[stage] * h * kx[stage - 1][i] : 0.0);
            us[i] = u[i] + (stage ? w[stage] * h * ku[stage - 1][i] : 0.0);
        }
        kx[stage] = us;
        ku[stage].assign(n, 0.0);
        geodesicAcceleration(st, xs.data(), us.data(), ku[stage].data());
    }
    for (size_t i = 0; i < n; ++i) {
        x[i] += h / 6.0 * (kx[0][i] + 2.0 * kx[1][i] + 2.0 * kx[2][i] + kx[3][i]);
        u[i] += h / 6.0 * (ku[0][i] + 2.0 * ku[1][i] + 2.0 * ku[2][i] + ku[3][i]);
    }
}

std::vector<double> geodesicDeviation(const Tensor& R, const std::vector<double>& u, const std::vector<double>& xi) {
    const int n = R.dim();
    Tensor U(n, "^"), X(n, "^");
    U.data() = u;
    X.data() = xi;
    return (einstein("^a_bcd ^b ^c ^d -> ^a", R, U, X, U) * -1.0).data();
}

} // namespace rf
