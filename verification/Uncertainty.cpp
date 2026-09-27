// Monte Carlo propagation of input uncertainty with Latin hypercube draws (see Uncertainty.h).
#include "Uncertainty.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>

namespace rf::verify {

double normalQuantile(double p) {
    // P. J. Acklam, "An algorithm for computing the inverse normal cumulative distribution
    // function" (2003): a rational function in the centre, another one in the tails.
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                               1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                               6.680131188771972e+01, -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                               -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
    p = std::clamp(p, 1e-12, 1.0 - 1e-12);
    const double pLow = 0.02425;
    if (p < pLow || p > 1 - pLow) { // the tails
        const double q = std::sqrt(-2 * std::log(p < pLow ? p : 1 - p));
        const double x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
                         ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
        return p < pLow ? x : -x;
    }
    const double q = p - 0.5, r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}

// One value of the input from a probability p in (0, 1).
static double drawInput(const UncertainInput& in, double p) {
    if (in.distribution == UncertainInput::Distribution::Uniform) return in.nominal + in.spread * (2 * p - 1);
    return in.nominal + in.spread * normalQuantile(p);
}

// The Latin hypercube: draws[s][i] = the value of input i in sample s.
static std::vector<std::vector<double>> latinHypercube(const std::vector<UncertainInput>& inputs, int samples, uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::vector<std::vector<double>> draws(size_t(samples), std::vector<double>(inputs.size()));
    std::vector<int> strata(static_cast<size_t>(samples));
    for (size_t i = 0; i < inputs.size(); ++i) {
        std::iota(strata.begin(), strata.end(), 0);
        std::shuffle(strata.begin(), strata.end(), rng);
        for (int s = 0; s < samples; ++s) {
            const double p = (strata[size_t(s)] + unit(rng)) / double(samples); // one point in its stratum
            draws[size_t(s)][i] = drawInput(inputs[i], p);
        }
    }
    return draws;
}

// The q-quantile of sorted values, linear between neighbours.
static double quantile(const std::vector<double>& sorted, double q) {
    if (sorted.empty()) return 0;
    const double x = q * double(sorted.size() - 1);
    const size_t i = size_t(std::floor(x));
    if (i + 1 >= sorted.size()) return sorted.back();
    return sorted[i] + (x - double(i)) * (sorted[i + 1] - sorted[i]);
}

UqSummary propagate(const std::vector<UncertainInput>& inputs, const std::function<double(const std::vector<double>&)>& model,
                    int samples, uint32_t seed) {
    UqSummary out;
    if (samples < 2) return out;
    for (const std::vector<double>& x : latinHypercube(inputs, samples, seed)) out.outputs.push_back(model(x));
    out.samples = samples;
    double sum = 0;
    for (double y : out.outputs) sum += y;
    out.mean = sum / samples;
    double ss = 0;
    for (double y : out.outputs) ss += (y - out.mean) * (y - out.mean);
    out.stddev = std::sqrt(ss / (samples - 1));
    std::vector<double> sorted = out.outputs;
    std::sort(sorted.begin(), sorted.end());
    out.lo95 = quantile(sorted, 0.025);
    out.hi95 = quantile(sorted, 0.975);
    return out;
}

} // namespace rf::verify
