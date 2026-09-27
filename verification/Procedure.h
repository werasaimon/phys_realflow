#pragma once
// The steps of the blind analysis as commands of rf_verify (the rules in Blinding.h and
// docs/10-verification.md):
//   seal     - write verification/sealed/<id>.ref with a fresh salt and print the commitment to
//              put into the registry (Case::commitment); the value then leaves the sources;
//   freeze   - write the case's acceptance rule and the commit into the unblinding log;
//   unblind  - only after a freeze with the same acceptance: check the sealed file against the
//              commitment, run the case, compute E, z and the verdict, write them into the log;
//   acceptance changes of cases the log already knows are written into the log on every run.
// Every command returns the process exit code (0 = done) and says what it did or why not.

#include "Blinding.h"
#include "Runner.h"

#include <string>
#include <vector>

namespace rf::verify {

struct ProcedurePaths {
    std::string sealedDir; // verification/sealed
    std::string logPath;   // verification/unblinding-log.md
    std::string commit;    // of the sources
};

int sealCase(const std::string& id, double value, double uncertainty, const ProcedurePaths& p);
int freezeCase(const Case& c, const ProcedurePaths& p);
int unblindCase(const Case& c, const ProcedurePaths& p, const RunOptions& options, int seeds);
// Appends an acceptance-change entry for every given case whose current acceptance differs from
// the last one in the log; returns how many were written.
int logAcceptanceChanges(const std::vector<const Case*>& cases, const ProcedurePaths& p);
// A "note" entry: a change to how a case measures or what it compares (why, before, after),
// so that a change of method is on record next to the acceptance history, never silent.
int noteCase(const Case& c, const std::string& text, const ProcedurePaths& p);

} // namespace rf::verify
