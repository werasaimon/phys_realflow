#pragma once
// Refinement studies: how fast the error falls as the grid (or the time step) is refined, and how
// large the numerical error of the finest answer is.
//
// With a constant refinement ratio r = h_coarse / h_fine and three answers f1 (fine), f2, f3
// (coarse), the observed order is
//     p = ln |(f3 - f2) / (f2 - f1)| / ln r,
// the Richardson extrapolation (the answer of an infinitely fine grid) is
//     f_exact ~ f1 + (f1 - f2) / (r^p - 1),
// and the grid convergence index of the fine answer (Roache 1998, "Verification and Validation in
// Computational Science and Engineering"; Celik et al. 2008, J. Fluids Eng. 130, 078001) is
//     GCI_fine = Fs |f1 - f2| / (r^p - 1),   Fs = 1.25 for three or more grids,
// an error band that holds the exact answer with ~95 % confidence. The answers are in the
// asymptotic range - where the order means something - when GCI_coarse / (r^p GCI_fine) ~ 1.
// ASME V&V 20-2009 turns the GCI into a standard uncertainty: u_num = GCI / 1.15 (for Fs = 1.25).
//
// When the exact answer is known, the errors e(h) give the order directly: p = ln(e_coarse /
// e_fine) / ln r for two levels, and for more levels the slope of the least-squares line
// through (ln h, ln e), with its standard error.

#include "Benchmark.h"

#include <vector>

namespace rf::verify {

struct RichardsonEstimate {
    double order = kNaN;        // observed order p
    double extrapolated = kNaN; // Richardson extrapolation of the answer
    double gciFine = kNaN;      // absolute, in the units of the answer
    double gciCoarse = kNaN;    // of the medium answer against the coarse one
    double asymptoticRatio = kNaN; // GCI_coarse / (r^p GCI_fine); 1 = asymptotic range
    bool oscillating = false;   // (f3 - f2) and (f2 - f1) of opposite sign: no order, GCI from the spread
};

// f1 fine, f2 medium, f3 coarse, r = h2 / h1 = h3 / h2.
RichardsonEstimate richardson(double f1, double f2, double f3, double r);

// The standard numerical uncertainty from a GCI (ASME V&V 20: u_num = GCI / 1.15 for Fs = 1.25).
inline double numericalUncertaintyFromGci(double gci) { return gci / 1.15; }

// Order from two errors: ln(e_coarse / e_fine) / ln r.
double orderFromErrors(double eCoarse, double eFine, double r);

struct OrderFit {
    double order = kNaN;    // slope of ln e against ln h
    double stdError = kNaN; // its standard error (kNaN for two points)
    double lnC = kNaN;      // intercept: e ~ C h^p
};

// Least squares through (ln h, ln |error|) of the points that have an error.
OrderFit fitOrder(const std::vector<ConvergencePoint>& points);

} // namespace rf::verify
