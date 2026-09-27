// The result file, the history of earlier runs and the console table (see Report.h). The board
// is in Board.cpp, the charts in Charts.cpp.
#include "Report.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace rf::verify {

std::string num(double x) {
    if (!std::isfinite(x)) return "—";
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.4g", x);
    return buf;
}

static std::string jsonNumber(double x) {
    if (!std::isfinite(x)) return "null";
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.10g", x);
    return buf;
}

static std::string jsonString(const std::string& s) {
    std::string out = "\"";
    for (char ch : s) {
        if (ch == '"' || ch == '\\') { out += '\\'; out += ch; }
        else if (ch == '\n') out += "\\n";
        else if (static_cast<unsigned char>(ch) < 0x20) out += ' ';
        else out += ch;
    }
    return out + "\"";
}

// One case on one line: the history reader needs no JSON parser. A blinded case has null
// reference, E, u_val and z (its Outcome carries NaN there) and the verdict BLIND.
static std::string jsonCase(const Outcome& o) {
    const Result& r = o.r;
    std::string s = "    {\"id\": " + jsonString(o.c->id) + ", \"category\": " + jsonString(o.c->category);
    s += ", \"status\": " + jsonString(statusWord(o.status)) + ", \"split\": " + jsonString(o.c->split);
    s += ", \"value\": " + jsonNumber(r.value) + ", \"u_num\": " + jsonNumber(r.numericalUncertainty);
    s += ", \"u_input\": " + jsonNumber(r.inputUncertainty) + ", \"runs\": " + std::to_string(r.runs) + ", \"spread\": " + jsonNumber(r.spread);
    s += ", \"reference\": " + jsonNumber(o.ref.value) + ", \"u_ref\": " + jsonNumber(o.ref.uncertainty);
    s += ", \"E\": " + jsonNumber(o.a.E) + ", \"u_val\": " + jsonNumber(o.a.uVal) + ", \"z\": " + jsonNumber(o.a.z);
    s += ", \"order\": " + jsonNumber(r.observedOrder) + ", \"order_theory\": " + jsonNumber(r.theoreticalOrder);
    s += ", \"verdict\": " + jsonString(verdictWord(o.a.verdict)) + ", \"seconds\": " + jsonNumber(r.seconds);
    s += ", \"unit\": " + jsonString(r.unit) + ", \"h_label\": " + jsonString(r.hLabel) + ", \"convergence\": [";
    for (size_t i = 0; i < r.convergence.size(); ++i) {
        const ConvergencePoint& p = r.convergence[i];
        s += (i ? ", [" : "[") + jsonNumber(p.h) + ", " + jsonNumber(p.value) + ", " + jsonNumber(p.error) + "]";
    }
    return s + "], \"detail\": " + jsonString(r.detail) + "}";
}

bool writeJson(const std::string& path, const RunInfo& info, const std::vector<Outcome>& outcomes) {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    f << "{\n  \"date\": " << jsonString(info.utc) << ",\n  \"commit\": " << jsonString(info.commit) << ",\n  \"mode\": "
      << jsonString(info.mode) << ",\n  \"seconds\": " << jsonNumber(info.seconds) << ",\n  \"cases\": [\n";
    for (size_t i = 0; i < outcomes.size(); ++i) f << jsonCase(outcomes[i]) << (i + 1 < outcomes.size() ? ",\n" : "\n");
    f << "  ]\n}\n";
    return bool(f);
}

// The value after `"key": ` on a line (NaN for null or missing).
static double fieldValue(const std::string& line, const std::string& key) {
    const size_t at = line.find("\"" + key + "\": ");
    if (at == std::string::npos) return kNaN;
    const char* p = line.c_str() + at + key.size() + 4;
    char* end = nullptr;
    const double v = std::strtod(p, &end);
    return end == p ? kNaN : v;
}

// The string after `"key": "` on a line, with \" \\ \n undone.
static std::string fieldString(const std::string& line, const std::string& key) {
    const size_t at = line.find("\"" + key + "\": \"");
    std::string out;
    if (at == std::string::npos) return out;
    for (size_t i = at + key.size() + 5; i < line.size() && line[i] != '"'; ++i) {
        if (line[i] == '\\' && i + 1 < line.size()) { ++i; out += line[i] == 'n' ? '\n' : line[i]; }
        else out += line[i];
    }
    return out;
}

History readHistory(const std::string& resultsDir, const std::vector<std::string>& skipFiles) {
    namespace fs = std::filesystem;
    History h;
    std::error_code ec;
    std::vector<fs::path> files;
    for (const auto& e : fs::directory_iterator(resultsDir, ec))
        if (e.path().extension() == ".json" &&
            std::find(skipFiles.begin(), skipFiles.end(), e.path().filename().string()) == skipFiles.end())
            files.push_back(e.path());
    std::sort(files.begin(), files.end()); // the names start with the UTC time: oldest first
    for (const fs::path& file : files) {
        std::ifstream f(file);
        for (std::string line; std::getline(f, line);) {
            const std::string id = fieldString(line, "id");
            if (id.empty()) continue;
            const double spread = fieldValue(line, "spread"), runs = fieldValue(line, "runs");
            const double sigma = std::isfinite(spread) && runs >= 1 ? spread / std::sqrt(runs) : 0;
            h[id].push_back({fieldValue(line, "value"), sigma});
        }
    }
    return h;
}

// "convergence": [[h, value, error], ...] of one case line.
static std::vector<ConvergencePoint> fieldConvergence(const std::string& line) {
    std::vector<ConvergencePoint> out;
    size_t at = line.find("\"convergence\": [");
    if (at == std::string::npos) return out;
    at += 16;
    while (at < line.size() && line[at] == '[') {
        const size_t end = line.find(']', at);
        std::stringstream ss(line.substr(at + 1, end - at - 1));
        double v[3] = {kNaN, kNaN, kNaN};
        std::string item;
        for (int k = 0; k < 3 && std::getline(ss, item, ','); ++k) v[k] = item.find("null") != std::string::npos ? kNaN : std::stod(item);
        out.push_back({v[0], v[1], v[2]});
        at = end + 1;
        while (at < line.size() && (line[at] == ',' || line[at] == ' ')) ++at;
    }
    return out;
}

// One case line back into an Outcome (the assessment recomputed from the reference in force).
static Outcome outcomeFromLine(const Case& c, const std::string& line) {
    Outcome o;
    o.c = &c;
    const std::string status = fieldString(line, "status");
    o.status = status == "blinded" ? BlindStatus::Blinded : status == "unblinded" ? BlindStatus::Unblinded : BlindStatus::Open;
    o.ref = c.ref;
    if (o.status == BlindStatus::Unblinded) { o.ref.value = fieldValue(line, "reference"); o.ref.uncertainty = fieldValue(line, "u_ref"); }
    Result& r = o.r;
    r.value = fieldValue(line, "value");
    r.numericalUncertainty = std::isfinite(fieldValue(line, "u_num")) ? fieldValue(line, "u_num") : 0;
    r.inputUncertainty = std::isfinite(fieldValue(line, "u_input")) ? fieldValue(line, "u_input") : 0;
    r.runs = std::isfinite(fieldValue(line, "runs")) ? int(fieldValue(line, "runs")) : 1;
    r.spread = std::isfinite(fieldValue(line, "spread")) ? fieldValue(line, "spread") : 0;
    r.observedOrder = fieldValue(line, "order");
    r.theoreticalOrder = fieldValue(line, "order_theory");
    r.seconds = fieldValue(line, "seconds");
    r.unit = fieldString(line, "unit");
    r.detail = fieldString(line, "detail");
    r.convergence = fieldConvergence(line);
    if (const std::string h = fieldString(line, "h_label"); !h.empty()) r.hLabel = h;
    o.a = assess(c, o.ref, r, o.status == BlindStatus::Blinded);
    return o;
}

bool readResults(const std::string& path, RunInfo& info, std::vector<Outcome>& outcomes) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    for (std::string line; std::getline(f, line);) {
        const std::string id = fieldString(line, "id");
        if (id.empty()) {
            if (line.find("\"date\": ") != std::string::npos) info.utc = fieldString(line, "date");
            if (line.find("\"commit\": ") != std::string::npos) info.commit = fieldString(line, "commit");
            if (line.find("\"mode\": ") != std::string::npos) info.mode = fieldString(line, "mode");
            if (line.find("\"seconds\": ") != std::string::npos) info.seconds = fieldValue(line, "seconds");
            continue;
        }
        for (const Case& c : registry())
            if (c.id == id) outcomes.push_back(outcomeFromLine(c, line));
    }
    return !outcomes.empty();
}

void printTable(const std::vector<Outcome>& outcomes) {
    std::printf("\n%-24s %12s %12s %10s %10s %7s %8s %6s %8s\n", "case", "ours", "reference", "E", "u_val", "z", "order", "", "seconds");
    for (const Outcome& o : outcomes) // a blinded reference is NaN here: printed as a dash
        std::printf("%-24s %12s %12s %10s %10s %7s %8s %6s %8.1f\n", o.c->id.c_str(), num(o.r.value).c_str(), num(o.ref.value).c_str(),
                    num(o.a.E).c_str(), num(o.a.uVal).c_str(), num(o.a.z).c_str(), num(o.r.observedOrder).c_str(), verdictWord(o.a.verdict),
                    o.r.seconds);
}

} // namespace rf::verify
