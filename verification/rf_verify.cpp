// rf_verify: runs the verification and validation registry of the SDK and writes the results.
//
//   rf_verify [--quick | --full] [--case <id>] [--board] [--uq <n>] [--seeds <k>]
//   rf_verify --list | --replot <results.json> [--replot <later.json> ...] | --seal <id> <value> <uncertainty> | --freeze <id> | --unblind <id> | --note <id> <text>
//   common:   [--root <dir>] [--out <dir>] [--docs <dir>]
//
//   --quick    the fast cases on coarse levels (default; a few minutes)
//   --full     every case, finer levels, Monte Carlo input uncertainty; writes the board
//   --case id  one case only (never writes the board)
//   --board    write the board and the charts after a quick run too
//   --uq n     Monte Carlo samples of the input uncertainty (default: 0 quick, 8 full)
//   --seeds k  runs of a stochastic case (default: 2 quick, 5 full)
//   --list     print the registry with the blind status of every case and stop
//   --replot f the board and the charts again from saved result files, without running; given
//              more than once, a later file's cases replace the same cases of the earlier ones (a case
//              re-run alone after a fix) and the board names where they come from (a run
//              saved without a commit, "nogit", is shown with the current one)
//   --seal     write verification/sealed/<id>.ref and print its commitment (blind analysis)
//   --freeze   write the case's acceptance rule and the commit into the unblinding log
//   --unblind  check the sealed reference, compute E, z and the verdict once, log them
//   --note     write a note into the unblinding log: a change of how a case measures (why, before, after)
//   --root dir the source tree the defaults below live in (default: this build's sources)
//   --out dir  where the result file goes (default: <root>/verification/results)
//   --docs dir where the board and the charts go (default: <root>/docs)
//
// Writes results/<UTC time>-<commit>.json always; with the board, docs/09-benchmarks.md (the
// generated block only) and the charts in docs/img/verification/. The exit code is the number
// of failed cases (errors included), so a script or CI can stop on it.
//
// Environment: RF_COMMIT (the commit, verbatim), GIT_DIR / GIT_WORK_TREE (read by git itself),
// RF_GIT_DIR (a git directory kept apart from the sources) - the order is in Runner.h.
#include "Procedure.h"
#include "Report.h"
#include "Runner.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>

#ifndef RF_SOURCE_DIR
#define RF_SOURCE_DIR "."
#endif

using namespace rf::verify;

struct Options {
    RunOptions run;
    std::vector<std::string> replot;
    std::string only, root = RF_SOURCE_DIR, outDir, docsDir, freeze, unblind, sealId, noteId, noteText;
    double sealValue = kNaN, sealUncertainty = kNaN;
    bool board = false, list = false;
    int uq = -1, seeds = -1;
    std::string sealedDir() const { return root + "/verification/sealed"; }
    std::string logPath() const { return root + "/verification/unblinding-log.md"; }
};

static bool parseOne(int argc, char** argv, int& i, Options& o) {
    const std::string a = argv[i];
    const bool v = i + 1 < argc;
    if (a == "--quick") o.run.full = false;
    else if (a == "--full") o.run.full = true;
    else if (a == "--board") o.board = true;
    else if (a == "--list") o.list = true;
    else if (a == "--case" && v) o.only = argv[++i];
    else if (a == "--out" && v) o.outDir = argv[++i];
    else if (a == "--docs" && v) o.docsDir = argv[++i];
    else if (a == "--root" && v) o.root = argv[++i];
    else if (a == "--replot" && v) o.replot.push_back(argv[++i]);
    else if (a == "--uq" && v) o.uq = std::atoi(argv[++i]);
    else if (a == "--seeds" && v) o.seeds = std::atoi(argv[++i]);
    else if (a == "--freeze" && v) o.freeze = argv[++i];
    else if (a == "--unblind" && v) o.unblind = argv[++i];
    else if (a == "--note" && i + 2 < argc) { o.noteId = argv[++i]; o.noteText = argv[++i]; }
    else if (a == "--seal" && i + 3 < argc) { o.sealId = argv[++i]; o.sealValue = std::atof(argv[++i]); o.sealUncertainty = std::atof(argv[++i]); }
    else return false;
    return true;
}

static bool parse(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; ++i)
        if (!parseOne(argc, argv, i, o)) {
            std::printf("unknown argument %s (see the head of verification/rf_verify.cpp for the usage)\n", argv[i]);
            return false;
        }
    if (o.outDir.empty()) o.outDir = o.root + "/verification/results";
    if (o.docsDir.empty()) o.docsDir = o.root + "/docs";
    o.run.uqSamples = o.uq >= 0 ? o.uq : (o.run.full ? 8 : 0);
    o.run.levels = o.run.full ? 4 : 3;
    if (o.seeds < 0) o.seeds = o.run.full ? 5 : 2;
    if (o.run.full) o.board = true;
    if (!o.only.empty()) o.board = false;
    return true;
}

static const Case* findCase(const std::string& id) {
    for (const Case& c : registry())
        if (c.id == id) return &c;
    std::printf("no case %s in the registry\n", id.c_str());
    return nullptr;
}

static std::vector<const Case*> chooseCases(const Options& o) {
    std::vector<const Case*> chosen;
    for (const Case& c : registry())
        if (!o.only.empty() ? c.id == o.only : (o.run.full || !c.slow)) chosen.push_back(&c);
    return chosen;
}

static void writeBoard(const RunInfo& info, const std::vector<Outcome>& outcomes, const Options& o, const std::vector<std::string>& skipInHistory) {
    const std::string images = o.docsDir + "/img/verification";
    std::error_code ec;
    std::filesystem::create_directories(images, ec); // an existing folder or a refusal: the writes below report it
    const std::vector<std::string> charts = writeCharts(images, outcomes);
    const History history = readHistory(o.outDir, skipInHistory);
    const bool ok = updateBoard(o.docsDir + "/09-benchmarks.md", boardBlock(info, outcomes, history, "img/verification"));
    std::printf("board %s: %s/09-benchmarks.md; %zu charts in %s\n", ok ? "written" : "NOT WRITTEN", o.docsDir.c_str(), charts.size(), images.c_str());
}

static int failedCount(const std::vector<Outcome>& outcomes) {
    int failed = 0;
    for (const Outcome& x : outcomes) failed += x.a.verdict == Verdict::Fail || x.a.verdict == Verdict::Error;
    return failed;
}

// One result file into the outcomes: a case already there is replaced (a later file wins) and
// the replacement is named in `from`.
static bool mergeResults(const std::string& file, RunInfo& info, std::vector<Outcome>& outcomes, std::string& from) {
    RunInfo later;
    std::vector<Outcome> read;
    if (!readResults(file, later, read)) return false;
    if (outcomes.empty()) { info = later; outcomes = read; return true; }
    for (Outcome& o : read) {
        bool replaced = false;
        for (Outcome& old : outcomes)
            if (old.c == o.c) { old = o; replaced = true; }
        if (!replaced) outcomes.push_back(o);
        from += (from.empty() ? "" : ", ") + std::string("`") + o.c->id + "`";
    }
    from += " — из прогона " + later.utc + " (`--" + later.mode + "`, коммит `" + later.commit + "`)";
    return true;
}

// --replot a.json [--replot b.json ...]: the board and the charts from saved runs, nothing is
// simulated; later files replace the same cases of earlier ones.
static int replot(const Options& o) {
    RunInfo info;
    std::vector<Outcome> outcomes;
    std::vector<std::string> names;
    std::string from;
    for (const std::string& file : o.replot) {
        std::string one;
        if (!mergeResults(file, info, outcomes, one)) {
            std::printf("cannot read results from %s\n", file.c_str());
            return 2;
        }
        if (!names.empty()) from += (from.empty() ? "" : "; ") + one;
        names.push_back(std::filesystem::path(file).filename().string());
    }
    if (info.commit.empty() || info.commit == "nogit") {
        info.commit = commitHash(o.root);
        std::printf("the run was saved without a commit; shown with the current one: %s\n", info.commit.c_str());
    }
    if (!from.empty()) info.note = "Заменены более поздним прогоном: " + from + ".";
    printTable(outcomes);
    writeBoard(info, outcomes, o, names);
    return failedCount(outcomes);
}

static int list(const Options& o) {
    const std::vector<LogEntry> log = readLog(o.logPath());
    for (const Case& c : registry())
        std::printf("%-24s %-22s %-9s %-7s %s%s\n", c.id.c_str(), c.category.c_str(), statusWord(statusOf(c, log)), c.split.c_str(),
                    c.slow ? "[full] " : "", c.ref.source.c_str());
    return 0;
}

static int runRegistry(const Options& o, const ProcedurePaths& paths) {
    const std::vector<const Case*> cases = chooseCases(o);
    if (cases.empty()) {
        std::printf("no case matches\n");
        return 2;
    }
    logAcceptanceChanges(cases, paths);
    RunInfo info{utcNow(false), paths.commit, o.run.full ? "full" : "quick", 0};
    std::printf("rf_verify --%s: %zu cases, %d seeds for stochastic ones, commit %s, %s\n", info.mode.c_str(), cases.size(), o.seeds,
                info.commit.c_str(), info.utc.c_str());
    const auto t0 = std::chrono::steady_clock::now();
    const std::vector<Outcome> outcomes = runCases(cases, o.run, Environment{paths.sealedDir, readLog(paths.logPath), o.seeds}, true);
    info.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    printTable(outcomes);
    std::error_code ec;
    std::filesystem::create_directories(o.outDir, ec);
    const std::string resultName = utcNow(true) + "-" + info.commit + (o.only.empty() ? "" : "-" + o.only) + ".json";
    const bool saved = writeJson(o.outDir + "/" + resultName, info, outcomes);
    std::printf("\nresults: %s/%s %s(%.0f s)\n", o.outDir.c_str(), resultName.c_str(), saved ? "" : "NOT WRITTEN - the folder is not writable ",
                info.seconds);
    if (o.board) writeBoard(info, outcomes, o, {resultName});
    return failedCount(outcomes);
}

int main(int argc, char** argv) {
    Options o;
    if (!parse(argc, argv, o)) return 2;
    const ProcedurePaths paths{o.sealedDir(), o.logPath(), commitHash(o.root)};
    if (o.list) return list(o);
    if (!o.replot.empty()) return replot(o);
    if (!o.sealId.empty()) return sealCase(o.sealId, o.sealValue, o.sealUncertainty, paths);
    if (!o.freeze.empty()) {
        const Case* c = findCase(o.freeze);
        return c ? freezeCase(*c, paths) : 2;
    }
    if (!o.noteId.empty()) {
        const Case* c = findCase(o.noteId);
        return c ? noteCase(*c, o.noteText, paths) : 2;
    }
    if (!o.unblind.empty()) {
        const Case* c = findCase(o.unblind);
        return c ? unblindCase(*c, paths, o.run, o.seeds) : 2;
    }
    return runRegistry(o, paths);
}
