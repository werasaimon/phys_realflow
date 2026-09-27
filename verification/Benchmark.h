#pragma once
// The benchmark registry: every verification or validation case of the SDK as data - what is
// computed, what it is compared with (a reference with its source and its uncertainty), how the
// case runs and when it counts as passed. The structure follows ASME V&V 10-2019 / V&V 20-2009 and
// Oberkampf & Roy, "Verification and Validation in Scientific Computing" (2010):
//   code verification     - is the equation solved right? An exact answer (analytic or
//                           manufactured) and the observed order of accuracy against the
//                           theoretical one;
//   solution verification - how large is the numerical error of this answer? Refinement in space
//                           or time, Richardson extrapolation, the grid convergence index;
//   validation            - is it the right equation? The answer against an experiment, with the
//                           experiment's uncertainty, the input uncertainty and the numerical one.
// rf_verify runs the registry, writes the results and generates the board (docs/09).

#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace rf::verify {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Where the reference number comes from and how sure it is (a standard uncertainty, 1 sigma).
// kind: "analytic" (u = 0: an exact solution), "experiment", "code" (another solver's result).
struct Reference {
    std::string source;     // author, year, what
    std::string doi_or_url; // how to find it
    double value = 0;
    double uncertainty = 0;
    std::string kind = "analytic";
    std::string note;       // how the value or its uncertainty was obtained
};

// One level of a refinement study: the step h (grid spacing, time step, particle spacing), the
// computed value and, when an exact answer is known, the error of the level.
struct ConvergencePoint {
    double h = 0;
    double value = 0;
    double error = kNaN;
};

// How a case is run: quick (the fast cases on coarse levels, minutes) or full (all cases on
// finer levels); the number of refinement levels; Monte Carlo samples for the input uncertainty
// (0 = none); the seed of every random choice.
struct RunOptions {
    bool full = false;
    int levels = 3;
    int uqSamples = 0;
    uint32_t seed = 1;
};

// What a case returns. value is the quantity compared with the reference; numericalUncertainty
// (u_num) comes from the refinement study (GCI / 1.15, ASME V&V 20), inputUncertainty (u_input)
// from the Monte Carlo propagation of uncertain inputs. Both are standard uncertainties.
struct Result {
    double value = kNaN;
    double numericalUncertainty = 0;
    double inputUncertainty = 0;
    std::vector<ConvergencePoint> convergence;
    std::string unit;
    double observedOrder = kNaN;    // of the refinement study, if there is one
    double theoreticalOrder = kNaN; // what the scheme should give
    std::string hLabel = "h";       // what h is in the convergence chart ("dx", "dt", ...)
    std::string detail;             // one line of the numbers behind the value
    double seconds = 0;             // wall time of the run
    int runs = 1;                   // seeds run (a stochastic case: the value is their mean)
    double spread = 0;              // standard deviation of the value over the seeds
};

// When the case passes. |E| = |value - reference| (relative to |reference| if relative) at most
// tolerance: pass; at most warnTolerance: warning; more: failed. The comparison error E is also
// judged against the validation uncertainty u_val (ValidationMetric.h).
struct Acceptance {
    double tolerance = 0;
    double warnTolerance = 0;
    bool relative = false;
};

struct Case {
    std::string id;       // short, stable: the key of the result history
    std::string title;    // for the board (Russian, as the docs)
    std::string category; // "code-verification" | "solution-verification" | "validation"
    std::string test;     // the rf_tests test of the same physics, if any
    Reference ref;
    std::function<Result(const RunOptions&)> run;
    Acceptance acceptance;
    bool slow = false;    // only in --full
    // The split of the double blind (like training and hold-out sets): "tuning" - developers may
    // look at the comparison and tune on it; "holdout" - blinded until a release (Blinding.h).
    std::string split = "tuning";
    // SHA-256 of the sealed reference (verification/sealed/<id>.ref); empty = an open case whose
    // reference is in `ref`. A blinded case keeps ref.value and ref.uncertainty at NaN.
    std::string commitment;
    // Sensitive to its initial state (a splash, a wake): run with several seeds that perturb the
    // start slightly; the spread of the answers is the run-to-run scatter (Runner.cpp).
    bool stochastic = false;
};

// Every case, in the order of the board. Registered by the case files (cases/).
const std::vector<Case>& registry();

// The case files add their cases here (one function per file, called by registry()).
void addRigidCases(std::vector<Case>& cases);
void addGasCases(std::vector<Case>& cases);
void addParticleCases(std::vector<Case>& cases);

} // namespace rf::verify
