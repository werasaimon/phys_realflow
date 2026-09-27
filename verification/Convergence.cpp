// Observed order, Richardson extrapolation, the grid convergence index and the least-squares
// order of a refinement study (formulas and sources in Convergence.h).
#include "Convergence.h"

#include <algorithm>
#include <cmath>

namespace rf::verify {

RichardsonEstimate richardson(double f1, double f2, double f3, double r) {
    RichardsonEstimate out;
    const double e21 = f2 - f1, e32 = f3 - f2;
    const double Fs = 1.25;
    if (e21 == 0) { // the two finest answers agree to the last bit: nothing left to estimate
        out.order = kNaN;
        out.extrapolated = f1;
        out.gciFine = out.gciCoarse = 0;
        return out;
    }
    if (e32 / e21 < 0) {
        // Oscillating convergence: no order. Roache's safety factor for an unproven order (3)
        // times the largest change is the error band.
        out.oscillating = true;
        out.extrapolated = f1;
        out.gciFine = 3.0 * std::max(std::fabs(e21), std::fabs(e32));
        return out;
    }
    const double p = std::log(std::fabs(e32 / e21)) / std::log(r);
    out.order = p;
    if (!(p > 0.05)) { // diverging or stalled: the refinement does not shrink the change
        out.extrapolated = f1;
        out.gciFine = 3.0 * std::max(std::fabs(e21), std::fabs(e32));
        return out;
    }
    const double rp = std::pow(r, p);
    out.extrapolated = f1 - e21 / (rp - 1.0);
    out.gciFine = Fs * std::fabs(e21) / (rp - 1.0);
    out.gciCoarse = Fs * std::fabs(e32) / (rp - 1.0);
    out.asymptoticRatio = out.gciFine > 0 ? out.gciCoarse / (rp * out.gciFine) : kNaN;
    return out;
}

double orderFromErrors(double eCoarse, double eFine, double r) {
    if (!(eCoarse > 0) || !(eFine > 0)) return kNaN;
    return std::log(eCoarse / eFine) / std::log(r);
}

OrderFit fitOrder(const std::vector<ConvergencePoint>& points) {
    std::vector<double> x, y;
    for (const ConvergencePoint& p : points)
        if (std::isfinite(p.error) && std::fabs(p.error) > 0 && p.h > 0) {
            x.push_back(std::log(p.h));
            y.push_back(std::log(std::fabs(p.error)));
        }
    OrderFit fit;
    const size_t n = x.size();
    if (n < 2) return fit;
    double mx = 0, my = 0;
    for (size_t i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
    mx /= double(n);
    my /= double(n);
    double sxx = 0, sxy = 0;
    for (size_t i = 0; i < n; ++i) { sxx += (x[i] - mx) * (x[i] - mx); sxy += (x[i] - mx) * (y[i] - my); }
    if (sxx <= 0) return fit;
    fit.order = sxy / sxx;
    fit.lnC = my - fit.order * mx;
    if (n > 2) { // standard error of the slope from the residuals
        double ss = 0;
        for (size_t i = 0; i < n; ++i) { const double r = y[i] - (fit.lnC + fit.order * x[i]); ss += r * r; }
        fit.stdError = std::sqrt(ss / double(n - 2) / sxx);
    }
    return fit;
}

} // namespace rf::verify
