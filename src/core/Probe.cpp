#include "core/Probe.h"

#include <algorithm>
#include <cstring>
#include <mutex>
#include <unordered_map>

namespace rf {

std::atomic<long long> Probe::allocations{0};
std::atomic<bool> Probe::drawEnabled_{false};

namespace {

struct Store {
    std::mutex mutex;
    std::unordered_map<std::string, Probe::Channel> channels;
    std::vector<Probe::Line> lines;
    std::vector<Probe::Point> points;
};

// One store for the process, made on first use (safe from static initialisation order).
Store& store() {
    static Store s;
    return s;
}

Probe::Channel& channel(Store& s, const char* name, Probe::Kind kind) {
    auto it = s.channels.find(name);
    if (it == s.channels.end()) it = s.channels.emplace(name, Probe::Channel{name, 0.0, kind}).first;
    it->second.kind = kind;
    return it->second;
}

} // namespace

void Probe::set(const char* name, double value) {
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    channel(s, name, Kind::Value).value = value;
}

void Probe::add(const char* name, double delta) {
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    channel(s, name, Kind::Counter).value += delta;
}

Probe::Timer::~Timer() {
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0_).count();
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    channel(s, name_, Kind::TimerMs).value += ms;
}

void Probe::line(const Vector3& a, const Vector3& b, const Vector3& color) {
    if (!drawEnabled()) return;
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.lines.push_back({a, b, color});
}

void Probe::arrow(const Vector3& p, const Vector3& v, const Vector3& color) {
    if (!drawEnabled()) return;
    const float len = length(v);
    if (len < 1e-9f) return;
    const Vector3 tip = p + v, d = v / len;
    // A head of two short strokes in a plane through the arrow.
    const Vector3 side = normalize(std::fabs(d.y) < 0.9f ? cross(d, Vector3(0, 1, 0)) : cross(d, Vector3(1, 0, 0)));
    const float h = 0.15f * len;
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.lines.push_back({p, tip, color});
    s.lines.push_back({tip, tip - d * h + side * (0.5f * h), color});
    s.lines.push_back({tip, tip - d * h - side * (0.5f * h), color});
}

void Probe::point(const Vector3& p, const Vector3& color, float size) {
    if (!drawEnabled()) return;
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.points.push_back({p, color, size});
}

void Probe::box(const AABB& b, const Vector3& color) {
    if (!drawEnabled()) return;
    const Vector3& lo = b.lo;
    const Vector3& hi = b.hi;
    const Vector3 c[8] = {{lo.x, lo.y, lo.z}, {hi.x, lo.y, lo.z}, {hi.x, hi.y, lo.z}, {lo.x, hi.y, lo.z},
                          {lo.x, lo.y, hi.z}, {hi.x, lo.y, hi.z}, {hi.x, hi.y, hi.z}, {lo.x, hi.y, hi.z}};
    static const int e[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    for (const auto& k : e) s.lines.push_back({c[k[0]], c[k[1]], color});
}

void Probe::beginFrame() {
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    for (auto& [name, c] : s.channels)
        if (c.kind != Kind::Value) c.value = 0;
    s.lines.clear();
    s.points.clear();
}

Probe::Snapshot Probe::snapshot() {
    Store& s = store();
    Snapshot out;
    std::lock_guard<std::mutex> lock(s.mutex);
    out.channels.reserve(s.channels.size());
    for (const auto& [name, c] : s.channels) out.channels.push_back(c);
    std::sort(out.channels.begin(), out.channels.end(), [](const Channel& a, const Channel& b) { return a.name < b.name; });
    out.lines = s.lines;
    out.points = s.points;
    return out;
}

std::vector<std::string> Probe::channels() {
    Store& s = store();
    std::vector<std::string> names;
    std::lock_guard<std::mutex> lock(s.mutex);
    for (const auto& [name, c] : s.channels) names.push_back(name);
    std::sort(names.begin(), names.end());
    return names;
}

void Probe::clear() {
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.channels.clear();
    s.lines.clear();
    s.points.clear();
}

double Probe::Snapshot::value(const std::string& name, double fallback) const {
    auto it = std::lower_bound(channels.begin(), channels.end(), name,
                               [](const Channel& c, const std::string& n) { return c.name < n; });
    return it != channels.end() && it->name == name ? it->value : fallback;
}

} // namespace rf
