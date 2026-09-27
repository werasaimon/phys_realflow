// Smoke test of the verification registry (verification/, rf_verify): the refinement formulas give
// the known answers on made-up sequences, the Monte Carlo spread of a linear model is its known
// standard deviation, SHA-256 matches the FIPS 180-4 test vectors, a sealed reference goes through
// seal -> freeze -> unblind (and a changed sealed file is refused), chi^2 p-values match the
// tables, the pull classes switch at 2, 3 and 5 sigma, two fast cases run and pass, and the
// board block goes between its markers of a scratch board without touching the hand-kept rest.
#include "TestRunner.h"
#include "Tests.h"

#include "Convergence.h"
#include "Procedure.h"
#include "Runner.h"
#include "Sha256.h"
#include "Significance.h"
#include "Uncertainty.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

using namespace rf::verify;

static void checkFormulas() {
    // f(h) = 1 + 0.5 h^2 on h = 0.4, 0.2, 0.1: order 2, extrapolation 1, GCI = 1.25 |f1 - f2| / 3.
    auto f = [](double h) { return 1 + 0.5 * h * h; };
    const RichardsonEstimate re = richardson(f(0.1), f(0.2), f(0.4), 2.0);
    CHECK(std::fabs(re.order - 2) < 1e-9, "Richardson order %f, expected 2", re.order);
    CHECK(std::fabs(re.extrapolated - 1) < 1e-12, "Richardson extrapolation %f, expected 1", re.extrapolated);
    CHECK(std::fabs(re.gciFine - 1.25 * std::fabs(f(0.1) - f(0.2)) / 3) < 1e-12, "GCI %e", re.gciFine);
    CHECK(std::fabs(re.asymptoticRatio - 1) < 1e-9, "asymptotic ratio %f, expected 1", re.asymptoticRatio);
    std::vector<ConvergencePoint> pts;
    for (double h : {0.4, 0.2, 0.1, 0.05}) pts.push_back({h, 0, 3 * h * h * h});
    const OrderFit fit = fitOrder(pts);
    CHECK(std::fabs(fit.order - 3) < 1e-9, "least-squares order %f, expected 3", fit.order);
    // y = 2 x1 + x2, x1 ~ N(0, 0.1), x2 ~ U(-0.3, 0.3): sigma = sqrt(0.04 + 0.03) = 0.2646.
    const UqSummary uq = propagate({{"x1", 0, 0.1, UncertainInput::Distribution::Normal}, {"x2", 0, 0.3, UncertainInput::Distribution::Uniform}},
                                   [](const std::vector<double>& x) { return 2 * x[0] + x[1]; }, 400, 7);
    CHECK(std::fabs(uq.stddev / std::sqrt(0.07) - 1) < 0.05, "Monte Carlo sigma %f, expected %f", uq.stddev, std::sqrt(0.07));
    CHECK(std::fabs(normalQuantile(0.975) - 1.959964) < 1e-5, "normal quantile %f", normalQuantile(0.975));
}

// FIPS 180-4 examples (NIST CSRC "Examples with intermediate values"): one block, two blocks,
// the empty message and a million times 'a'.
static void checkSha256() {
    const std::pair<std::string, const char*> vectors[] = {
        {"abc", "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"},
        {"", "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"},
        {"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq", "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"},
        {std::string(1000000, 'a'), "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0"}};
    for (const auto& [message, expected] : vectors)
        CHECK(sha256Hex(message) == expected, "SHA-256 of a %zu-byte message: %s", message.size(), sha256Hex(message).c_str());
}

// Chi-squared tail probabilities from the tables (Abramowitz & Stegun 26.4) and closed forms.
static void checkStatistics() {
    // {chi2, ndf, p, tolerance}; the last row is the normal approximation's neighbourhood (~0.481).
    const double table[][4] = {{3.841459, 1, 0.05, 2e-6}, {18.307038, 10, 0.05, 2e-6}, {11.344867, 3, 0.01, 2e-6}, {2.0, 2, 0.36787944, 1e-8},
                               {100, 100, 0.4812, 5e-4}};
    for (const auto& row : table) {
        const double p = chiSquaredPValue(row[0], int(row[1]));
        CHECK(std::fabs(p - row[2]) < row[3], "chi2 = %g, ndf = %g: p = %.7f, table %.7f", row[0], row[1], p, row[2]);
    }
    CHECK(std::fabs(regularizedGammaP(1, 2.5) - (1 - std::exp(-2.5))) < 1e-12, "P(1, x) = 1 - e^-x");
    CHECK(std::fabs(twoSidedP(1.959964) - 0.05) < 1e-6, "two-sided p of 1.96 sigma");
    CHECK(sigmaClass(1.99) == SigmaClass::Agrees && sigmaClass(-2.0) == SigmaClass::Tension && sigmaClass(2.99) == SigmaClass::Tension,
          "pull classes around 2 sigma");
    CHECK(sigmaClass(3.0) == SigmaClass::Evidence && sigmaClass(-4.99) == SigmaClass::Evidence && sigmaClass(5.0) == SigmaClass::Established,
          "pull classes around 3 and 5 sigma");
    CHECK(regressionMark({{1.0, 0.01}}, 1.05, 0.01).rfind("↑ 3.5σ", 0) == 0, "a 3.5 sigma rise: %s", regressionMark({{1.0, 0.01}}, 1.05, 0.01).c_str());
    CHECK(regressionMark({{1.0, 0.01}}, 1.02, 0.01).rfind("≈", 0) == 0, "a 1.4 sigma change is not flagged");
}

// Seal -> (unblind refused before a freeze) -> freeze -> unblind -> status unblinded; a changed
// acceptance goes into the log; a sealed file changed afterwards no longer matches.
static void checkSealedRoundTrip() {
    namespace fs = std::filesystem;
    const fs::path dir = fs::temp_directory_path() / "rf_verify_blind_test";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const ProcedurePaths paths{dir.string(), (dir / "log.md").string(), "test123"};
    SealedReference sealed{2.5, 0.1, newSalt()};
    CHECK(writeSealed(paths.sealedDir, "toy", sealed), "cannot write a sealed file into %s", dir.string().c_str());
    Case toy{"toy", "toy", "validation", "", {"toy", "", kNaN, kNaN, "experiment", ""},
             [](const RunOptions&) { Result r; r.value = 2.6; r.numericalUncertainty = 0.05; return r; }, {0.1, 0.2, true}, false, "holdout",
             commitmentOf(sealed)};
    CHECK(statusOf(toy, readLog(paths.logPath)) == BlindStatus::Blinded, "a sealed case starts blinded");
    CHECK(unblindCase(toy, paths, RunOptions(), 1) != 0, "unblinding without a freeze must be refused");
    CHECK(freezeCase(toy, paths) == 0 && unblindCase(toy, paths, RunOptions(), 1) == 0, "freeze then unblind");
    const std::vector<LogEntry> log = readLog(paths.logPath);
    CHECK(log.size() == 2 && log.back().event == "unblind" && statusOf(toy, log) == BlindStatus::Unblinded, "the log holds freeze and unblind");
    toy.acceptance.tolerance = 0.2;
    CHECK(logAcceptanceChanges({&toy}, paths) == 1 && readLog(paths.logPath).back().event == "acceptance-change", "a changed rule is logged");
    sealed.value = 2.55; // tampering: same salt, another value
    writeSealed(paths.sealedDir, "toy", sealed);
    Reference ref;
    std::string why;
    CHECK(!referenceInForce(toy, BlindStatus::Unblinded, paths.sealedDir, ref, why), "a changed sealed file must not match its commitment");
    fs::remove_all(dir);
}

// A pass with a 7.6 sigma pull (the incline's stick angle: E = 0.032 deg within a 0.25 deg
// tolerance, u_num = 0.0042 deg) must read "resolved bias" with a footnote, not a bare tick.
static void checkResolvedBias() {
    Case c{"toy-bias", "toy", "code-verification", "", {"toy", "", 21.8014, 0, "analytic", ""},
           [](const RunOptions&) { return Result(); }, {0.25, 1.0, false}};
    Outcome o;
    o.c = &c;
    o.ref = c.ref;
    o.r.value = 21.8335;
    o.r.numericalUncertainty = 0.0042;
    o.a = assess(c, o.ref, o.r, false);
    const std::string block = boardBlock(RunInfo{"2026-01-01T00:00:00Z", "test", "quick", 1}, {o}, History(), "img");
    CHECK(o.a.verdict == Verdict::Pass && o.a.z > 7, "the toy case must pass with z > 7 (z = %f)", o.a.z);
    CHECK(block.find("✅ смещение разрешено¹") != std::string::npos && block.find("- ¹ `toy-bias`: E = 0.0321") != std::string::npos,
          "a pass with |z| >= 2 must show the resolved-bias status and its footnote");
}

void testVerificationSmoke() {
    checkFormulas();
    checkSha256();
    checkStatistics();
    checkSealedRoundTrip();
    checkResolvedBias();
    std::vector<const Case*> chosen;
    for (const Case& c : registry())
        if (c.id == "rigid-projectile" || c.id == "rigid-incline-slide") chosen.push_back(&c);
    CHECK(chosen.size() == 2, "the registry lost a fast case (%zu found)", chosen.size());
    const std::vector<Outcome> outcomes = runCases(chosen, RunOptions(), Environment(), false);
    for (const Outcome& o : outcomes)
        std::printf("  %s: %s (reference %s), %s\n", o.c->id.c_str(), num(o.r.value).c_str(), num(o.ref.value).c_str(), verdictWord(o.a.verdict));
    for (const Outcome& o : outcomes) CHECK(o.a.verdict == Verdict::Pass, "%s: %s", o.c->id.c_str(), o.r.detail.c_str());

    // The board: a scratch file with a hand-kept section; written twice, the block is replaced,
    // not added again, and the hand-kept text stays.
    const std::string path = (std::filesystem::temp_directory_path() / "rf_verify_smoke_board.md").string();
    { std::ofstream f(path, std::ios::binary); f << "# Board\n\nIntro.\n\n## Hand kept\n\n| row | 1 |\n"; }
    RunInfo info{"2026-01-01T00:00:00Z", "test", "quick", 1};
    const std::string block = boardBlock(info, outcomes, History(), "img/verification");
    const bool written = updateBoard(path, block) && updateBoard(path, block);
    std::stringstream text;
    text << std::ifstream(path, std::ios::binary).rdbuf();
    const std::string s = text.str();
    const size_t first = s.find("rf_verify: начало"), second = s.find("rf_verify: начало", first + 1);
    CHECK(written && first != std::string::npos && second == std::string::npos, "the block must be in the board exactly once");
    CHECK(s.find("| row | 1 |") != std::string::npos && s.find("`rigid-projectile`") != std::string::npos, "board lost a row");
    std::filesystem::remove(path);
}
