#pragma once
// The validation metric of ASME V&V 20-2009: the comparison error
//     E = S - D                                  (simulation minus reference data)
// and the validation uncertainty
//     u_val = sqrt(u_num^2 + u_input^2 + u_D^2)
// (numerical, input and reference-data standard uncertainties; u_D = 0 for an exact solution).
// The model error delta_model = E - (errors of the simulation and of the data) then lies in
// [E - u_val, E + u_val] (one standard uncertainty). Its significance is the pull z = E / u_val,
// read in the language of particle physics (Significance.h): |z| < 2 agrees, 2..3 tension,
// 3..5 evidence of a model error, >= 5 established.
// The acceptance of the case (Benchmark.h) then says whether |E| is good enough for us: pass,
// warning or fail. A blinded case gets neither E nor a verdict (Blinding.h).

#include "Benchmark.h"
#include "Significance.h"

#include <string>

namespace rf::verify {

enum class Verdict { Pass, Warn, Fail, Error, Blinded };

struct Assessment {
    double E = kNaN;      // S - D
    double relativeE = kNaN;
    double uVal = kNaN;   // validation uncertainty
    double z = kNaN;      // pull E / u_val (NaN when u_val = 0)
    SigmaClass sigma = SigmaClass::Undefined;
    double deltaLo = kNaN, deltaHi = kNaN; // the model error interval E -+ u_val
    Verdict verdict = Verdict::Error;
    std::string wording;   // the reading in words (Russian, for the board)
};

// ref: the reference in force - the registry's for an open case, the sealed one once unblinded.
Assessment assess(const Case& c, const Reference& ref, const Result& r, bool blinded);

const char* verdictSymbol(Verdict v); // board: ✅ 🟡 ❌ ⚠️ 🔒
const char* verdictWord(Verdict v);   // console and JSON: PASS WARN FAIL ERROR BLIND

} // namespace rf::verify
