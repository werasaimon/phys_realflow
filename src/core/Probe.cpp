// The probe's store: one map of channels behind a mutex, the debug drawings of the frame by layer
// (each capped at kLayerCap), the research pointers (watched pair, probe point) and the allocation
// counter. What the probe is for and how to use it is explained in Probe.h.
#include "core/Probe.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <map>
#include <mutex>
#include <string_view>

namespace rf {

std::atomic<long long> Probe::allocations{0};
std::atomic<uint32_t> Probe::layers_{0};
std::atomic<uint64_t> Probe::frames_{1};

namespace {

constexpr int kLayers = int(DrawLayer::Count);

struct Store {
    std::mutex mutex;
    // A map with a transparent comparator: find(const char*) compares without building a
    // std::string, so reporting a channel allocates only the first time its name is seen.
    std::map<std::string, Probe::Channel, std::less<>> channels;
    std::vector<Probe::Line> lines;
    std::vector<Probe::Point> points;
    std::vector<Probe::Label> labels;
    std::array<int, kLayers> used{};    // primitives of each layer this frame
    std::array<bool, kLayers> capped{}; // the "(обрезано)" label is out for this layer
    int watchA = -1, watchB = -1;
    Vector3 probePoint{0.0f};
};

// One store for the process, made on first use (safe from static initialisation order).
Store& store() {
    static Store s;
    return s;
}

Probe::Channel& channel(Store& s, const char* name, Probe::Kind kind) {
    auto it = s.channels.find(std::string_view(name));
    if (it == s.channels.end()) it = s.channels.emplace(std::string(name), Probe::Channel{name, 0.0, kind}).first;
    it->second.kind = kind;
    return it->second;
}

// Room for n more primitives of layer l this frame (the store is locked). The first refusal leaves
// a label at `at` so the viewer knows the picture is incomplete.
bool reserve(Store& s, DrawLayer l, int n, const Vector3& at) {
    int& used = s.used[size_t(l)];
    if (used + n <= Probe::kLayerCap) {
        used += n;
        return true;
    }
    if (!s.capped[size_t(l)]) {
        s.capped[size_t(l)] = true;
        s.labels.push_back({at, std::string(Probe::layerName(l)) + " (обрезано)", l});
    }
    return false;
}

} // namespace

const char* Probe::layerName(DrawLayer l) {
    static const char* names[] = {"Misc", "ContactPoints", "ContactNormals", "ContactImpulses", "PenetrationDepth",
                                  "BodyAabbs", "WorldTree", "MeshBvh", "CentreOfMass", "InertiaAxes", "Velocities",
                                  "Islands", "Sleeping", "JointFrames", "GjkSimplex", "EpaPolytope", "WitnessPoints",
                                  "ParticleNeighbours", "DensityError", "SoftClusters", "ClothTension", "GasGrid",
                                  "GasVelocity", "PressureGradient", "Divergence", "Vorticity", "FieldLinesB",
                                  "CurrentDensity"};
    static_assert(sizeof(names) / sizeof(names[0]) == size_t(DrawLayer::Count), "a name per layer");
    return uint32_t(l) < uint32_t(DrawLayer::Count) ? names[uint32_t(l)] : "?";
}

void Probe::enableDraw(bool on) {
    const uint32_t basic = bit(DrawLayer::Misc) | bit(DrawLayer::ContactPoints) | bit(DrawLayer::ContactNormals) |
                           bit(DrawLayer::BodyAabbs);
    setLayers(on ? basic : 0u);
}

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

void Probe::line(DrawLayer l, const Vector3& a, const Vector3& b, const Vector3& color) {
    if (!layerOn(l)) return;
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (reserve(s, l, 1, a)) s.lines.push_back({a, b, color, l});
}

void Probe::arrow(DrawLayer l, const Vector3& p, const Vector3& v, const Vector3& color) {
    if (!layerOn(l)) return;
    const float len = length(v);
    if (len < 1e-9f) return;
    const Vector3 tip = p + v, d = v / len;
    // A head of two short strokes in a plane through the arrow.
    const Vector3 side = normalize(std::fabs(d.y) < 0.9f ? cross(d, Vector3(0, 1, 0)) : cross(d, Vector3(1, 0, 0)));
    const float h = 0.15f * len;
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (!reserve(s, l, 3, p)) return;
    s.lines.push_back({p, tip, color, l});
    s.lines.push_back({tip, tip - d * h + side * (0.5f * h), color, l});
    s.lines.push_back({tip, tip - d * h - side * (0.5f * h), color, l});
}

void Probe::point(DrawLayer l, const Vector3& p, const Vector3& color, float size) {
    if (!layerOn(l)) return;
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (reserve(s, l, 1, p)) s.points.push_back({p, color, size, l});
}

void Probe::box(DrawLayer l, const AABB& b, const Vector3& color) {
    if (!layerOn(l)) return;
    const Vector3& lo = b.lo;
    const Vector3& hi = b.hi;
    const Vector3 c[8] = {{lo.x, lo.y, lo.z}, {hi.x, lo.y, lo.z}, {hi.x, hi.y, lo.z}, {lo.x, hi.y, lo.z},
                          {lo.x, lo.y, hi.z}, {hi.x, lo.y, hi.z}, {hi.x, hi.y, hi.z}, {lo.x, hi.y, hi.z}};
    static const int e[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (!reserve(s, l, 12, lo)) return;
    for (const auto& k : e) s.lines.push_back({c[k[0]], c[k[1]], color, l});
}

void Probe::label(DrawLayer l, const Vector3& at, const std::string& text) {
    if (!layerOn(l)) return;
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (reserve(s, l, 1, at)) s.labels.push_back({at, text, l});
}

void Probe::watchPair(int bodyA, int bodyB) {
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.watchA = bodyA;
    s.watchB = bodyB;
}

void Probe::watchedPair(int& bodyA, int& bodyB) {
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    bodyA = s.watchA;
    bodyB = s.watchB;
}

void Probe::setProbePoint(const Vector3& p) {
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.probePoint = p;
}

Vector3 Probe::probePoint() {
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    return s.probePoint;
}

void Probe::beginFrame() {
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    for (auto& [name, c] : s.channels)
        if (c.kind != Kind::Value) c.value = 0;
    s.lines.clear();
    s.points.clear();
    s.labels.clear();
    s.used.fill(0);
    s.capped.fill(false);
    ++frames_;
}

void Probe::clearLayers(uint32_t mask) {
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    auto inMask = [mask](const auto& primitive) { return (mask & bit(primitive.layer)) != 0; };
    s.lines.erase(std::remove_if(s.lines.begin(), s.lines.end(), inMask), s.lines.end());
    s.points.erase(std::remove_if(s.points.begin(), s.points.end(), inMask), s.points.end());
    s.labels.erase(std::remove_if(s.labels.begin(), s.labels.end(), inMask), s.labels.end());
    for (int l = 0; l < kLayers; ++l)
        if (mask & bit(DrawLayer(l))) {
            s.used[size_t(l)] = 0;
            s.capped[size_t(l)] = false;
        }
}

Probe::Snapshot Probe::snapshot() {
    Store& s = store();
    Snapshot out;
    std::lock_guard<std::mutex> lock(s.mutex);
    out.channels.reserve(s.channels.size());
    for (const auto& [name, c] : s.channels) out.channels.push_back(c); // a std::map: already sorted by name
    out.lines = s.lines;
    out.points = s.points;
    out.labels = s.labels;
    return out;
}

std::vector<std::string> Probe::channels() {
    Store& s = store();
    std::vector<std::string> names;
    std::lock_guard<std::mutex> lock(s.mutex);
    for (const auto& [name, c] : s.channels) names.push_back(name); // sorted: it is a std::map
    return names;
}

void Probe::clear() {
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.channels.clear();
    s.lines.clear();
    s.points.clear();
    s.labels.clear();
    s.used.fill(0);
    s.capped.fill(false);
    ++frames_;
}

size_t Probe::primitiveCount() {
    Store& s = store();
    std::lock_guard<std::mutex> lock(s.mutex);
    return s.lines.size() + s.points.size() + s.labels.size();
}

double Probe::Snapshot::value(const std::string& name, double fallback) const {
    auto it = std::lower_bound(channels.begin(), channels.end(), name,
                               [](const Channel& c, const std::string& n) { return c.name < n; });
    return it != channels.end() && it->name == name ? it->value : fallback;
}

Vector3 heatColor(float t) {
    t = std::min(std::max(t, 0.0f), 1.0f);
    if (t < 0.5f) return Vector3(0.0f, 2.0f * t, 1.0f - 2.0f * t) + Vector3(0.1f);     // blue -> green
    return Vector3(2.0f * (t - 0.5f), 1.0f - 2.0f * (t - 0.5f), 0.0f) + Vector3(0.1f); // green -> red
}

} // namespace rf
