// Seal, freeze, unblind and the logging of acceptance changes (see Procedure.h).
#include "Procedure.h"

#include <cmath>
#include <cstdio>

namespace rf::verify {

int sealCase(const std::string& id, double value, double uncertainty, const ProcedurePaths& p) {
    SealedReference s{value, uncertainty, newSalt()};
    if (!writeSealed(p.sealedDir, id, s)) {
        std::printf("cannot write %s/%s.ref\n", p.sealedDir.c_str(), id.c_str());
        return 2;
    }
    std::printf("sealed %s into %s/%s.ref\ncommitment (put it into the case's registry entry):\n%s\n", id.c_str(), p.sealedDir.c_str(), id.c_str(),
                commitmentOf(s).c_str());
    return 0;
}

int freezeCase(const Case& c, const ProcedurePaths& p) {
    if (c.commitment.empty()) {
        std::printf("%s is an open case: nothing to freeze\n", c.id.c_str());
        return 2;
    }
    const LogEntry e{utcNow(false), "freeze", c.id, p.commit, c.commitment, acceptanceText(c.acceptance), "правило приёмки заморожено до снятия слепоты"};
    if (!appendLog(p.logPath, e)) {
        std::printf("cannot write the log %s\n", p.logPath.c_str());
        return 2;
    }
    std::printf("frozen: %s, %s, commit %s\n", c.id.c_str(), e.acceptance.c_str(), p.commit.c_str());
    return 0;
}

// The last freeze of this case under its current commitment ("" acceptance if none).
static const LogEntry* lastFreeze(const Case& c, const std::vector<LogEntry>& log) {
    const LogEntry* found = nullptr;
    for (const LogEntry& e : log)
        if (e.id == c.id && e.event == "freeze" && e.commitment == c.commitment) found = &e;
    return found;
}

// Why the case may not be unblinded now ("" = it may).
static std::string unblindRefusal(const Case& c, const std::vector<LogEntry>& log) {
    const BlindStatus status = statusOf(c, log);
    if (status == BlindStatus::Open) return "открытый случай: снимать нечего";
    if (status == BlindStatus::Unblinded) return "слепота уже снята (запись в журнале)";
    const LogEntry* freeze = lastFreeze(c, log);
    if (!freeze) return "нет записи freeze: сначала rf_verify --freeze " + c.id;
    if (freeze->acceptance != acceptanceText(c.acceptance))
        return "правило приёмки изменилось после freeze (" + freeze->acceptance + " -> " + acceptanceText(c.acceptance) + "): заморозьте снова";
    return "";
}

int unblindCase(const Case& c, const ProcedurePaths& p, const RunOptions& options, int seeds) {
    const std::vector<LogEntry> log = readLog(p.logPath);
    if (const std::string why = unblindRefusal(c, log); !why.empty()) {
        std::printf("unblind %s refused: %s\n", c.id.c_str(), why.c_str());
        return 2;
    }
    Reference ref;
    std::string why;
    if (!referenceInForce(c, BlindStatus::Unblinded, p.sealedDir, ref, why)) {
        std::printf("unblind %s refused: %s\n", c.id.c_str(), why.c_str());
        return 2;
    }
    const std::vector<Outcome> run = runCases({&c}, options, Environment{p.sealedDir, log, seeds}, true);
    const Result& r = run.front().r;
    const Assessment a = assess(c, ref, r, false);
    char details[400];
    std::snprintf(details, sizeof(details), "S = %.6g ± %.3g ± %.3g (UQ); D = %.6g ± %.3g; E = %.4g; u_val = %.3g; z = %+.2f; вердикт %s",
                  r.value, r.numericalUncertainty, r.inputUncertainty, ref.value, ref.uncertainty, a.E, a.uVal, a.z, verdictWord(a.verdict));
    const LogEntry e{utcNow(false), "unblind", c.id, p.commit, c.commitment, acceptanceText(c.acceptance), details};
    if (!appendLog(p.logPath, e)) {
        std::printf("cannot write the log %s\n", p.logPath.c_str());
        return 2;
    }
    std::printf("unblinded %s: %s\n%s\n", c.id.c_str(), details, a.wording.c_str());
    return 0;
}

int noteCase(const Case& c, const std::string& text, const ProcedurePaths& p) {
    const std::string commitment = c.commitment.empty() ? "— (открытый случай)" : c.commitment;
    const LogEntry e{utcNow(false), "note", c.id, p.commit, commitment, acceptanceText(c.acceptance), text};
    if (!appendLog(p.logPath, e)) {
        std::printf("cannot write the log %s\n", p.logPath.c_str());
        return 2;
    }
    std::printf("noted for %s in %s\n", c.id.c_str(), p.logPath.c_str());
    return 0;
}

int logAcceptanceChanges(const std::vector<const Case*>& cases, const ProcedurePaths& p) {
    const std::vector<LogEntry> log = readLog(p.logPath);
    int written = 0;
    for (const Case* c : cases) {
        const std::string last = lastLoggedAcceptance(c->id, log), now = acceptanceText(c->acceptance);
        if (last.empty() || last == now) continue;
        const LogEntry e{utcNow(false), "acceptance-change", c->id, p.commit, c->commitment, now, "было: " + last};
        if (appendLog(p.logPath, e)) ++written;
        std::printf("WARNING: the acceptance of %s changed (%s -> %s); written into %s\n", c->id.c_str(), last.c_str(), now.c_str(), p.logPath.c_str());
    }
    return written;
}

} // namespace rf::verify
