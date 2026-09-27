#pragma once
// What rf_verify writes after a run:
//   results/<UTC time>-<commit>.json - every case's numbers, one case per line (the history of
//                                      every later run reads them back); a blinded case carries
//                                      no reference, no E and no verdict;
//   the console table;
//   the board, docs/09-benchmarks.md - only the block between the rf_verify markers is replaced,
//                                      everything else in the file (the rows still kept by hand)
//                                      is read from disk right before writing and kept as it is
//                                      (Board.cpp);
//   SVG charts (Charts.cpp)         - the refinement studies on log-log axes with the slope of
//                                      the theoretical order, E +- u_val of the validation cases,
//                                      and the pulls of every case against the 1/2/3 sigma bands.

#include "Benchmark.h"
#include "Blinding.h"
#include "ValidationMetric.h"

#include <map>
#include <string>
#include <vector>

namespace rf::verify {

struct Outcome {
    const Case* c = nullptr;
    Reference ref;                         // in force: the registry's, or the sealed one once unblinded
    BlindStatus status = BlindStatus::Open;
    Result r;
    Assessment a;
};

struct RunInfo {
    std::string utc;    // 2026-09-27T21:05:00Z
    std::string commit; // short hash or "nogit"
    std::string mode;   // quick | full
    double seconds = 0;
    std::string note;   // said after the run line of the board (e.g. cases taken from a later run)
};

// One earlier result of a case: the value and its run-to-run standard error (spread / sqrt(runs)).
struct HistoryPoint {
    double value = kNaN;
    double sigma = 0;
};
// id -> the results of earlier runs, oldest first (from the result files in a directory).
using History = std::map<std::string, std::vector<HistoryPoint>>;

bool writeJson(const std::string& path, const RunInfo& info, const std::vector<Outcome>& outcomes); // false: not written
History readHistory(const std::string& resultsDir, const std::vector<std::string>& skipFiles);
// A result file read back (for --replot: the board and the charts again without running): the
// cases found in the registry, their results, the assessment recomputed. False if unreadable.
bool readResults(const std::string& path, RunInfo& info, std::vector<Outcome>& outcomes);
void printTable(const std::vector<Outcome>& outcomes);

// Board.cpp: the generated block of the board (markdown), with charts under imagesDir (relative
// to the board), and putting it into the board file between the markers (else before the first
// section). updateBoard returns false if the file cannot be read or written.
std::string boardBlock(const RunInfo& info, const std::vector<Outcome>& outcomes, const History& history, const std::string& imagesDir);
bool updateBoard(const std::string& boardPath, const std::string& block);
// The regression mark of a value against the last earlier result: ↑ / ↓ with the change in units
// of sqrt(sigma_a^2 + sigma_b^2) when it exceeds 3 of them, ≈ otherwise.
std::string regressionMark(const std::vector<HistoryPoint>& earlier, double value, double sigma);

// Charts.cpp: one SVG per case with a refinement study; one for all validation cases; the pull
// plot. Returns the file names written (relative to dir).
std::vector<std::string> writeCharts(const std::string& dir, const std::vector<Outcome>& outcomes);

// Numbers for people: 4 significant digits, "—" for not-a-number.
std::string num(double x);

} // namespace rf::verify
