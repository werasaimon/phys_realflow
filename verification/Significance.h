#pragma once
// Significance in standard deviations, as particle physics reports it (PDG, "Statistics" review;
// Cowan, "Statistical Data Analysis", 1998):
//   the pull     z = (S - D) / u_val   - how many standard uncertainties apart ours and the reference are;
//   its reading  |z| < 2  agrees,  2..3  tension,  3..5  evidence of a model error,  >= 5  established
//                (the 5-sigma convention of discovery: a two-sided p of 5.7e-7 under "no error");
//   the set      chi^2 = sum z_i^2 over N independent cases, chi^2 / ndf with ndf = N, and the p-value
//                P(chi^2_N >= observed) from the regularized upper incomplete gamma Q(N/2, chi^2/2);
//   look elsewhere: among N cases some pull is large by chance. Bonferroni: the global p of the
//                largest pull is at most N p_local; and the expected count of |z| > 2 among N
//                agreeing cases is N * 0.0455 - two of forty is normal, not a finding.

#include <string>
#include <vector>

namespace rf::verify {

enum class SigmaClass { Agrees, Tension, Evidence, Established, Undefined };

SigmaClass sigmaClass(double z);
const char* sigmaWords(SigmaClass s, bool validation); // Russian, for the board

// Two-sided p-value of a standard normal pull: P(|Z| >= |z|) = erfc(|z| / sqrt 2).
double twoSidedP(double z);

// Regularized lower incomplete gamma P(a, x) = gamma(a, x) / Gamma(a): the series for x < a + 1,
// the continued fraction (modified Lentz) for Q = 1 - P beyond (Press et al., Numerical Recipes 6.2).
double regularizedGammaP(double a, double x);

// P(chi^2 with ndf degrees of freedom >= chi2) = Q(ndf / 2, chi2 / 2).
double chiSquaredPValue(double chi2, int ndf);

struct SetSummary {
    int n = 0;                  // cases with a finite pull
    double chi2 = 0;
    double pValue = 1;
    double largestPull = 0;     // signed
    std::string largestId;
    double bonferroniP = 1;     // min(1, n * two-sided p of the largest pull)
    double expectedOver2 = 0;   // n * 0.0455
    int observedOver2 = 0;
};

// Over pairs (case id, pull) with finite pulls.
SetSummary summarize(const std::vector<std::pair<std::string, double>>& pulls);

} // namespace rf::verify
