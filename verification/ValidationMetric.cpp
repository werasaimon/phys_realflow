// The ASME V&V 20 comparison of a result with its reference, its pull, and the verdict of the case.
// The verdict looks at |E| only (against the case's tolerance), never at |E| - u_val: a large
// uncertainty must not turn a poor answer into a pass. u_val and z say how far E can be trusted.
#include "ValidationMetric.h"

#include "core/Format.h"

#include <cmath>

namespace rf::verify {

static std::string number(double x) {
    if (!std::isfinite(x)) return "—";
    const double a = std::fabs(x);
    if (a != 0 && (a < 1e-3 || a >= 1e5)) return rf::format("%.2e", x);
    return rf::format(a < 1 ? "%.4f" : "%.3f", x);
}

static Verdict acceptanceVerdict(const Acceptance& acc, const Assessment& a) {
    const double measure = acc.relative ? std::fabs(a.relativeE) : std::fabs(a.E);
    if (!std::isfinite(measure)) return Verdict::Error;
    if (measure <= acc.tolerance) return Verdict::Pass;
    if (measure <= acc.warnTolerance) return Verdict::Warn;
    return Verdict::Fail;
}

Assessment assess(const Case& c, const Reference& ref, const Result& r, bool blinded) {
    Assessment a;
    if (blinded) {
        a.verdict = Verdict::Blinded;
        a.wording = "слепой случай: эталон запечатан, E и вердикт скрыты до --unblind";
        return a;
    }
    if (!std::isfinite(r.value)) {
        a.verdict = Verdict::Error;
        a.wording = "расчёт не дал числа";
        return a;
    }
    const double D = ref.value;
    a.E = r.value - D;
    a.relativeE = D != 0 ? a.E / std::fabs(D) : kNaN;
    a.uVal = std::sqrt(r.numericalUncertainty * r.numericalUncertainty + r.inputUncertainty * r.inputUncertainty +
                       ref.uncertainty * ref.uncertainty);
    a.deltaLo = a.E - a.uVal;
    a.deltaHi = a.E + a.uVal;
    a.z = a.uVal > 0 ? a.E / a.uVal : kNaN;
    a.sigma = sigmaClass(a.z);
    a.verdict = acceptanceVerdict(c.acceptance, a);
    // Only a validation case compares a model with nature; in verification E is a numerical
    // discrepancy against an exact answer, and u_val is the uncertainty we estimated for it.
    const bool validation = c.category == "validation";
    a.wording = std::isfinite(a.z) ? "z = " + rf::format("%+.1f", a.z) + " σ: " + sigmaWords(a.sigma, validation) : sigmaWords(a.sigma, validation);
    if (validation && std::isfinite(a.z) && a.sigma != SigmaClass::Agrees)
        a.wording += ", δ_model ∈ [" + number(a.deltaLo) + ", " + number(a.deltaHi) + "]";
    return a;
}

const char* verdictSymbol(Verdict v) {
    switch (v) {
    case Verdict::Pass: return "✅";
    case Verdict::Warn: return "🟡";
    case Verdict::Fail: return "❌";
    case Verdict::Error: return "⚠️";
    case Verdict::Blinded: return "🔒";
    }
    return "?";
}

const char* verdictWord(Verdict v) {
    switch (v) {
    case Verdict::Pass: return "PASS";
    case Verdict::Warn: return "WARN";
    case Verdict::Fail: return "FAIL";
    case Verdict::Error: return "ERROR";
    case Verdict::Blinded: return "BLIND";
    }
    return "?";
}

} // namespace rf::verify
