// Running the cases, the reference in force, the UTC time and the commit hash (see Runner.h).
#include "Runner.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>

namespace rf::verify {

static Result runOnce(const Case& c, const RunOptions& options) {
    try {
        return c.run(options);
    } catch (const std::exception& e) {
        Result r;
        r.detail = std::string("исключение: ") + e.what();
        return r;
    }
}

// A stochastic case K times with seeds seed, seed + 1, ...: the value is the mean, the spread the
// standard deviation of the K values; its standard error spread / sqrt(K) joins u_num in
// quadrature. The input uncertainty (Monte Carlo) is computed on the first seed only.
static Result runSeeds(const Case& c, const RunOptions& options, int seeds) {
    if (!c.stochastic || seeds <= 1) return runOnce(c, options);
    std::vector<Result> runs;
    for (int k = 0; k < seeds; ++k) {
        RunOptions o = options;
        o.seed = options.seed + uint32_t(k);
        if (k > 0) o.uqSamples = 0;
        runs.push_back(runOnce(c, o));
    }
    Result r = runs.front();
    double mean = 0, uNum = 0;
    for (const Result& x : runs) { mean += x.value; uNum += x.numericalUncertainty; }
    mean /= seeds;
    uNum /= seeds;
    double ss = 0;
    for (const Result& x : runs) ss += (x.value - mean) * (x.value - mean);
    r.value = mean;
    r.runs = seeds;
    r.spread = std::sqrt(ss / (seeds - 1));
    r.numericalUncertainty = std::sqrt(uNum * uNum + r.spread * r.spread / seeds);
    for (size_t i = 0; i < r.convergence.size(); ++i) { // the refinement levels averaged over the seeds
        double v = 0, e = 0;
        for (const Result& x : runs) { v += x.convergence[i].value; e += x.convergence[i].error; }
        r.convergence[i].value = v / seeds;
        r.convergence[i].error = e / seeds;
    }
    char buf[96];
    std::snprintf(buf, sizeof(buf), "; %d зёрен: среднее %.4g, разброс σ = %.3g", seeds, mean, r.spread);
    r.detail += buf;
    return r;
}

bool referenceInForce(const Case& c, BlindStatus status, const std::string& sealedDir, Reference& out, std::string& why) {
    out = c.ref;
    if (status == BlindStatus::Open) return true;
    if (status == BlindStatus::Blinded) { out.value = out.uncertainty = kNaN; return true; }
    SealedReference sealed;
    if (!readSealed(sealedDir, c.id, sealed)) { why = "нет запечатанного файла " + c.id + ".ref"; return false; }
    if (commitmentOf(sealed) != c.commitment) { why = "запечатанный файл не совпадает с обязательством (изменён?)"; return false; }
    out.value = sealed.value;
    out.uncertainty = sealed.uncertainty;
    return true;
}

std::vector<Outcome> runCases(const std::vector<const Case*>& cases, const RunOptions& options, const Environment& env, bool verbose) {
    std::vector<Outcome> out;
    for (const Case* c : cases) {
        Outcome o;
        o.c = c;
        o.status = statusOf(*c, env.log);
        std::string why;
        const bool refOk = referenceInForce(*c, o.status, env.sealedDir, o.ref, why);
        const auto t0 = std::chrono::steady_clock::now();
        o.r = runSeeds(*c, options, env.seeds);
        o.r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        o.a = assess(*c, o.ref, o.r, o.status == BlindStatus::Blinded);
        if (!refOk) { o.a = Assessment(); o.a.wording = why; }
        if (verbose) {
            std::printf("[%-5s] %-24s %s (reference %s)  %.1f s\n", verdictWord(o.a.verdict), c->id.c_str(), num(o.r.value).c_str(),
                        o.status == BlindStatus::Blinded ? "sealed" : num(o.ref.value).c_str(), o.r.seconds);
            std::fflush(stdout);
        }
        out.push_back(std::move(o));
    }
    return out;
}

std::string utcNow(bool forFileName) {
    const std::time_t t = std::time(nullptr);
    std::tm utc{};
#ifdef _WIN32
    gmtime_s(&utc, &t);
#else
    gmtime_r(&t, &utc);
#endif
    char buf[32];
    std::strftime(buf, sizeof(buf), forFileName ? "%Y-%m-%dT%H-%M-%SZ" : "%Y-%m-%dT%H:%M:%SZ", &utc);
    return buf;
}

// Runs a shell command and gives its standard output; false if it failed. Plain std: std::system
// with the output sent to a temporary file (standard C++ has no pipes to a child process).
static bool commandOutput(const std::string& command, std::string& out) {
    namespace fs = std::filesystem;
    std::random_device device;
    std::error_code ec;
    const fs::path file = fs::temp_directory_path(ec) / ("rf_verify_" + std::to_string(device()) + ".txt");
#ifdef _WIN32
    const char* nowhere = "NUL";
#else
    const char* nowhere = "/dev/null";
#endif
    const int status = std::system((command + " > \"" + file.string() + "\" 2> " + nowhere).c_str());
    std::stringstream text;
    {
        std::ifstream f(file, std::ios::binary);
        text << f.rdbuf();
    }
    fs::remove(file, ec);
    out = text.str();
    return status == 0;
}

// "<short hash>" or "<short hash>+dirty" through one way of calling git; "" if it gives no hash.
// Dirty: `git status --porcelain` lists anything - a changed, added or untracked file. With
// --no-optional-locks the status does not rewrite the index, so it never collides with a commit
// made at the same moment.
static std::string gitCommit(const std::string& git) {
    std::string hash, status;
    if (!commandOutput(git + " rev-parse --short HEAD", hash)) return "";
    while (!hash.empty() && (hash.back() == '\n' || hash.back() == '\r' || hash.back() == ' ')) hash.pop_back();
    if (hash.empty() || hash.find_first_of(" \n") != std::string::npos) return "";
    if (!commandOutput(git + " --no-optional-locks status --porcelain", status)) return hash + "+unknown";
    return status.empty() ? hash : hash + "+dirty";
}

std::string commitHash(const std::string& sourceDir) {
    if (const char* commit = std::getenv("RF_COMMIT"); commit && *commit) return commit;
    // git in the source tree; git itself honours GIT_DIR (and GIT_WORK_TREE) if they are set.
    std::string commit = gitCommit("git -C \"" + sourceDir + "\"");
    if (commit.empty())
        if (const char* dir = std::getenv("RF_GIT_DIR"); dir && *dir) // a git directory kept apart from the sources
            commit = gitCommit("git --git-dir=\"" + std::string(dir) + "\" --work-tree=\"" + sourceDir + "\"");
    return commit.empty() ? "nogit" : commit;
}

} // namespace rf::verify
