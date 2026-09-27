// Sealed references, commitments and the unblinding log (see Blinding.h).
#include "Blinding.h"

#include "Sha256.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>

namespace rf::verify {

const char* statusWord(BlindStatus s) {
    switch (s) {
    case BlindStatus::Open: return "open";
    case BlindStatus::Blinded: return "blinded";
    case BlindStatus::Unblinded: return "unblinded";
    }
    return "?";
}

std::string commitmentOf(const SealedReference& s) {
    char buf[128];
    std::snprintf(buf, sizeof(buf), "%.17g|%.17g|", s.value, s.uncertainty);
    return sha256Hex(buf + s.salt);
}

std::string newSalt() {
    std::random_device device;
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%08x%08x%08x%08x", device(), device(), device(), device());
    return buf;
}

bool writeSealed(const std::string& dir, const std::string& id, const SealedReference& s) {
    std::ofstream f(dir + "/" + id + ".ref", std::ios::binary);
    char buf[160];
    std::snprintf(buf, sizeof(buf), "value=%.17g\nuncertainty=%.17g\nsalt=", s.value, s.uncertainty);
    f << "# Sealed reference of the blinded case " << id << ". Do not read while tuning (docs/10-verification.md).\n"
      << buf << s.salt << "\n";
    return bool(f);
}

bool readSealed(const std::string& dir, const std::string& id, SealedReference& out) {
    std::ifstream f(dir + "/" + id + ".ref", std::ios::binary);
    if (!f) return false;
    for (std::string line; std::getline(f, line);) {
        while (!line.empty() && (line.back() == '\r' || line.back() == ' ')) line.pop_back();
        if (line.rfind("value=", 0) == 0) out.value = std::strtod(line.c_str() + 6, nullptr);
        else if (line.rfind("uncertainty=", 0) == 0) out.uncertainty = std::strtod(line.c_str() + 12, nullptr);
        else if (line.rfind("salt=", 0) == 0) out.salt = line.substr(5);
    }
    return std::isfinite(out.value) && !out.salt.empty();
}

// Splits "a · b · c" at the middle dots.
static std::vector<std::string> fields(const std::string& line) {
    const std::string sep = " · ";
    std::vector<std::string> out;
    size_t from = 0;
    for (size_t at; (at = line.find(sep, from)) != std::string::npos; from = at + sep.size()) out.push_back(line.substr(from, at - from));
    out.push_back(line.substr(from));
    return out;
}

std::vector<LogEntry> readLog(const std::string& path) {
    std::vector<LogEntry> log;
    std::ifstream f(path, std::ios::binary);
    for (std::string line; std::getline(f, line);) {
        while (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.rfind("- ", 0) != 0) continue;
        const std::vector<std::string> v = fields(line.substr(2));
        if (v.size() < 6) continue;
        LogEntry e{v[0], v[1], v[2], v[3], v[4], v[5], v.size() > 6 ? v[6] : ""};
        if (e.commit.rfind("commit ", 0) == 0) e.commit = e.commit.substr(7);
        if (e.commitment.rfind("commitment ", 0) == 0) e.commitment = e.commitment.substr(11);
        log.push_back(e);
    }
    return log;
}

bool appendLog(const std::string& path, const LogEntry& e) {
    std::ofstream f(path, std::ios::binary | std::ios::app);
    f << "- " << e.date << " · " << e.event << " · " << e.id << " · commit " << e.commit << " · commitment " << e.commitment << " · "
      << e.acceptance << (e.details.empty() ? "" : " · " + e.details) << "\n";
    return bool(f);
}

std::string acceptanceText(const Acceptance& a) {
    char buf[128];
    std::snprintf(buf, sizeof(buf), "%s ≤ %.6g; предупреждение ≤ %.6g", a.relative ? "|E|/|D|" : "|E|", a.tolerance, a.warnTolerance);
    return buf;
}

BlindStatus statusOf(const Case& c, const std::vector<LogEntry>& log) {
    if (c.commitment.empty()) return BlindStatus::Open;
    for (const LogEntry& e : log)
        if (e.id == c.id && e.event == "unblind" && e.commitment == c.commitment) return BlindStatus::Unblinded;
    return BlindStatus::Blinded;
}

std::string lastLoggedAcceptance(const std::string& id, const std::vector<LogEntry>& log) {
    std::string last;
    for (const LogEntry& e : log)
        if (e.id == id) last = e.acceptance;
    return last;
}

} // namespace rf::verify
