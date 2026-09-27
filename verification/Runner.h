#pragma once
// Running cases: each one timed, a crash of one case (an exception) recorded as an error of that
// case instead of stopping the run, a stochastic case run with several seeds, the reference in
// force chosen by the blind status (open: the registry's; unblinded: the sealed file, checked
// against its commitment; blinded: none), every result judged by the validation metric. Also the
// two facts a result file is named by: the UTC time and the commit of the sources.

#include "Report.h"

#include <string>
#include <vector>

namespace rf::verify {

struct Environment {
    std::string sealedDir;      // verification/sealed
    std::vector<LogEntry> log;  // verification/unblinding-log.md, read
    int seeds = 1;              // runs of a stochastic case
};

// Runs the cases in order; prints one line per case as it finishes when verbose.
std::vector<Outcome> runCases(const std::vector<const Case*>& cases, const RunOptions& options, const Environment& env, bool verbose);

// The reference in force for a case; false (with the reason) if a sealed file is missing or does
// not match its commitment.
bool referenceInForce(const Case& c, BlindStatus status, const std::string& sealedDir, Reference& out, std::string& why);

// "2026-09-27T21:05:00Z" (display) and "2026-09-27T21-05-00Z" (file names: no colons on Windows).
std::string utcNow(bool forFileName);

// The commit of the sources, in this order:
//   $RF_COMMIT, verbatim, if set (a CI job that knows its commit);
//   `git -C <sourceDir> rev-parse --short HEAD` - git honours $GIT_DIR / $GIT_WORK_TREE if set;
//   `git --git-dir=$RF_GIT_DIR --work-tree=<sourceDir> ...` - a git directory kept apart from the
//   work tree (on the author's machine the SDK's is C:/Users/PC/rf-git);
//   else "nogit".
// From git, "+dirty" is appended when `git status --porcelain` is not empty: "1be1f22+dirty".
std::string commitHash(const std::string& sourceDir);

} // namespace rf::verify
