#pragma once
// Blind analysis (Klein & Roodman 2005, "Blind analysis in nuclear and particle physics",
// Annu. Rev. Nucl. Part. Sci. 55:141): whoever tunes the code must not see how far the answer is
// from the reference, or the tuning drifts towards the reference and the comparison proves
// nothing. So the reference of a blinded case is not in the sources:
//   verification/sealed/<id>.ref   holds value, uncertainty and a random salt;
//   the registry holds only the commitment SHA-256("value|uncertainty|salt") (%.17g numbers);
//   while blinded, rf_verify shows S, u_num and u_input, never D, E or the verdict;
//   --freeze <id> writes the acceptance rule and the commit into the unblinding log;
//   --unblind <id> checks the sealed file against the commitment (a changed file is refused),
//   computes E and the verdict once, and writes them into the log. From then on the case is
//   unblinded; any later change of its acceptance rule is written into the log as a new entry.
// Status of a case: open (no commitment), blinded, unblinded (the log has its unblind entry).

#include "Benchmark.h"

#include <string>
#include <vector>

namespace rf::verify {

enum class BlindStatus { Open, Blinded, Unblinded };
const char* statusWord(BlindStatus s); // open | blinded | unblinded

struct SealedReference {
    double value = kNaN;
    double uncertainty = kNaN;
    std::string salt;
};

std::string commitmentOf(const SealedReference& s);
std::string newSalt(); // 32 hex digits from std::random_device
bool writeSealed(const std::string& dir, const std::string& id, const SealedReference& s);
bool readSealed(const std::string& dir, const std::string& id, SealedReference& out);

// One line of verification/unblinding-log.md:
//   - <UTC> · <event> · <id> · commit <hash> · commitment <sha> · <acceptance> · <details>
// events: freeze, unblind, acceptance-change, note (a change of how a case measures).
struct LogEntry {
    std::string date, event, id, commit, commitment, acceptance, details;
};

std::vector<LogEntry> readLog(const std::string& path);
bool appendLog(const std::string& path, const LogEntry& e);

// The acceptance rule in words, as the log keeps it: "|E|/|D| ≤ 0.08; предупреждение ≤ 0.15".
std::string acceptanceText(const Acceptance& a);

BlindStatus statusOf(const Case& c, const std::vector<LogEntry>& log);
// The acceptance last written into the log for this case ("" if the log never mentions it).
std::string lastLoggedAcceptance(const std::string& id, const std::vector<LogEntry>& log);

} // namespace rf::verify
