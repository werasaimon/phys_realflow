// Pulls, their classes, and the chi-squared test of a set of pulls (see Significance.h).
#include "Significance.h"

#include <algorithm>
#include <cmath>

namespace rf::verify {

SigmaClass sigmaClass(double z) {
    if (!std::isfinite(z)) return SigmaClass::Undefined;
    const double a = std::fabs(z);
    if (a < 2) return SigmaClass::Agrees;
    if (a < 3) return SigmaClass::Tension;
    if (a < 5) return SigmaClass::Evidence;
    return SigmaClass::Established;
}

const char* sigmaWords(SigmaClass s, bool validation) {
    switch (s) {
    case SigmaClass::Agrees: return "согласие (|z| < 2)";
    case SigmaClass::Tension: return "напряжение (2 ≤ |z| < 3)";
    case SigmaClass::Evidence: return validation ? "свидетельство ошибки модели (3 ≤ |z| < 5)" : "свидетельство численной ошибки сверх оценки (3 ≤ |z| < 5)";
    case SigmaClass::Established: return validation ? "ошибка модели установлена (|z| ≥ 5)" : "численная ошибка сверх оценки установлена (|z| ≥ 5)";
    case SigmaClass::Undefined: return "z не определён (u_val = 0)";
    }
    return "";
}

double twoSidedP(double z) { return std::erfc(std::fabs(z) / std::sqrt(2.0)); }

// The series gamma(a, x) = e^-x x^a sum x^n / (a (a+1) ... (a+n)), for x < a + 1.
static double gammaSeries(double a, double x) {
    double term = 1.0 / a, sum = term;
    for (int n = 1; n < 1000; ++n) {
        term *= x / (a + n);
        sum += term;
        if (std::fabs(term) < std::fabs(sum) * 1e-15) break;
    }
    return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
}

// Q(a, x) by its continued fraction (modified Lentz), for x >= a + 1.
static double gammaContinuedFraction(double a, double x) {
    const double tiny = 1e-300;
    double b = x + 1 - a, c = 1 / tiny, d = 1 / b, h = d;
    for (int i = 1; i < 1000; ++i) {
        const double an = -i * (i - a);
        b += 2;
        d = an * d + b;
        if (std::fabs(d) < tiny) d = tiny;
        c = b + an / c;
        if (std::fabs(c) < tiny) c = tiny;
        d = 1 / d;
        const double delta = d * c;
        h *= delta;
        if (std::fabs(delta - 1) < 1e-15) break;
    }
    return h * std::exp(-x + a * std::log(x) - std::lgamma(a));
}

double regularizedGammaP(double a, double x) {
    if (x <= 0) return 0;
    return x < a + 1 ? gammaSeries(a, x) : 1 - gammaContinuedFraction(a, x);
}

double chiSquaredPValue(double chi2, int ndf) {
    if (ndf <= 0) return 1;
    const double a = 0.5 * ndf, x = 0.5 * chi2;
    if (x <= 0) return 1;
    return x < a + 1 ? 1 - gammaSeries(a, x) : gammaContinuedFraction(a, x);
}

SetSummary summarize(const std::vector<std::pair<std::string, double>>& pulls) {
    SetSummary s;
    for (const auto& [id, z] : pulls) {
        if (!std::isfinite(z)) continue;
        ++s.n;
        s.chi2 += z * z;
        s.observedOver2 += std::fabs(z) > 2;
        if (std::fabs(z) > std::fabs(s.largestPull)) { s.largestPull = z; s.largestId = id; }
    }
    s.pValue = chiSquaredPValue(s.chi2, s.n);
    s.bonferroniP = std::min(1.0, s.n * twoSidedP(s.largestPull));
    s.expectedOver2 = s.n * twoSidedP(2.0);
    return s;
}

} // namespace rf::verify
