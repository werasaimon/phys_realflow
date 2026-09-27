// The research view of the gas solver, drawn by layer into the Probe on one slice of cells (the
// viewer's slice plane): the grid itself, the velocity u, the acceleration the pressure gives
// (-grad p / rho: the projection's correction), the divergence div u that the projection removes
// (Chorin 1968; Bridson, "Fluid Simulation for Computer Graphics", ch. 5), the vorticity curl u,
// and in a plasma the current density J = curl B / mu0. Arrows are scaled so the longest one spans
// the spacing of the drawn cells; a label gives the largest value. Cells are skipped with a stride
// that keeps a slice to about 2500 cells. Only the public state of the solver is read.
#include "gas/GasSolver.h"

#include "core/Format.h"
#include "core/Probe.h"

#include <algorithm>
#include <cmath>
#include <functional>

namespace rf {

namespace {

using CellFn = std::function<void(int, int, int)>;
using VectorFn = std::function<Vector3(int, int, int)>;

// The cells of one slice, visited with a stride: axis = the slice normal, layer = its index.
struct Slice {
    int axis = 2, layer = 0, n[3] = {1, 1, 1}, stride = 1;
    int u() const { return axis == 0 ? 1 : 0; } // the two in-plane grid axes
    int v() const { return axis == 2 ? 1 : 2; }
    void forEach(const CellFn& fn) const {
        for (int b = 0; b < n[v()]; b += stride)
            for (int a = 0; a < n[u()]; a += stride) {
                int c[3];
                c[axis] = layer;
                c[u()] = a;
                c[v()] = b;
                fn(c[0], c[1], c[2]);
            }
    }
};

Slice makeSlice(const GasSolver& g, int axis, int layer) {
    Slice s;
    s.axis = std::min(std::max(axis, 0), 2);
    s.n[0] = g.nx();
    s.n[1] = g.ny();
    s.n[2] = g.nz();
    s.layer = std::min(std::max(layer, 0), s.n[s.axis] - 1);
    s.stride = std::max(1, int(std::ceil(std::sqrt(double(s.n[s.u()]) * s.n[s.v()] / 2500.0))));
    return s;
}

Vector3 unit(int axis) { return Vector3(axis == 0 ? 1.0f : 0.0f, axis == 1 ? 1.0f : 0.0f, axis == 2 ? 1.0f : 0.0f); }

Vector3 cellCentre(const GasSolver& g, int i, int j, int k) {
    return g.origin() + Vector3(float(i) + 0.5f, float(j) + 0.5f, float(k) + 0.5f) * g.dx();
}

// The central difference of f along `axis` at cell c [f per metre]; one-sided at the grid border
// and next to a solid cell (whose value is not the gas's).
float derivative(const GasSolver& g, const std::function<float(int, int, int)>& f, int axis, int i, int j, int k) {
    int lo[3] = {i, j, k}, hi[3] = {i, j, k};
    const int n[3] = {g.nx(), g.ny(), g.nz()};
    if (lo[axis] > 0) --lo[axis];
    if (hi[axis] < n[axis] - 1) ++hi[axis];
    if (g.solid(lo[0], lo[1], lo[2])) lo[axis] = (axis == 0 ? i : axis == 1 ? j : k);
    if (g.solid(hi[0], hi[1], hi[2])) hi[axis] = (axis == 0 ? i : axis == 1 ? j : k);
    const int cells = hi[axis] - lo[axis];
    if (cells == 0) return 0.0f;
    return (f(hi[0], hi[1], hi[2]) - f(lo[0], lo[1], lo[2])) / (float(cells) * g.dx());
}

// -grad p / rho of a gas cell [m/s^2]: the acceleration the pressure gives the gas there.
Vector3 pressureAcceleration(const GasSolver& g, int i, int j, int k) {
    if (g.solid(i, j, k)) return Vector3(0.0f);
    auto p = [&g](int a, int b, int c) { return g.cellValue(GridField::Pressure, a, b, c); };
    const Vector3 grad(derivative(g, p, 0, i, j, k), derivative(g, p, 1, i, j, k), derivative(g, p, 2, i, j, k));
    return grad * (-1.0f / std::max(g.params.fluidDensity, 1e-6f));
}

// curl u of a gas cell [1/s] from the cell-centred velocities around it.
Vector3 vorticity(const GasSolver& g, int i, int j, int k) {
    if (g.solid(i, j, k)) return Vector3(0.0f);
    auto u = [&g](int c) { return [&g, c](int a, int b, int d) { return g.cellVelocity(a, b, d)[c]; }; };
    const float dwdy = derivative(g, u(2), 1, i, j, k), dvdz = derivative(g, u(1), 2, i, j, k);
    const float dudz = derivative(g, u(0), 2, i, j, k), dwdx = derivative(g, u(2), 0, i, j, k);
    const float dvdx = derivative(g, u(1), 0, i, j, k), dudy = derivative(g, u(0), 1, i, j, k);
    return Vector3(dwdy - dvdz, dudz - dwdx, dvdx - dudy);
}

// The lines of the slice plane through the cell centres of its layer, every stride-th.
void drawGrid(const GasSolver& g, const Slice& s) {
    const Vector3 du = unit(s.u()) * g.dx(), dv = unit(s.v()) * g.dx(), color(0.45f, 0.5f, 0.6f);
    const Vector3 base = g.origin() + unit(s.axis) * ((float(s.layer) + 0.5f) * g.dx());
    for (int a = 0; a <= s.n[s.u()]; a += s.stride)
        Probe::line(DrawLayer::GasGrid, base + du * float(a), base + du * float(a) + dv * float(s.n[s.v()]), color);
    for (int b = 0; b <= s.n[s.v()]; b += s.stride)
        Probe::line(DrawLayer::GasGrid, base + dv * float(b), base + dv * float(b) + du * float(s.n[s.u()]), color);
}

// Arrows of a vector field on the slice, scaled so the longest spans the spacing of drawn cells.
void drawArrows(DrawLayer l, const GasSolver& g, const Slice& s, const VectorFn& value, const Vector3& color, const char* unitName) {
    float largest = 0;
    s.forEach([&](int i, int j, int k) { largest = std::max(largest, length(value(i, j, k))); });
    if (largest <= 0) return;
    const float scale = 0.9f * g.dx() * float(s.stride) / largest;
    Vector3 top = cellCentre(g, 0, 0, 0);
    s.forEach([&](int i, int j, int k) {
        const Vector3 v = value(i, j, k);
        if (length2(v) > 0) Probe::arrow(l, cellCentre(g, i, j, k), v * scale, color);
        top = cellCentre(g, i, j, k);
    });
    Probe::label(l, top, format("%s: макс. %.3g %s", Probe::layerName(l), double(largest), unitName));
}

// Every gas cell of the slice as a point coloured by div u: green 0, red outflow, blue inflow;
// the scale is the largest |div u| of the slice.
void drawDivergence(const GasSolver& g, const Slice& s) {
    float largest = 1e-9f;
    s.forEach([&](int i, int j, int k) {
        if (!g.solid(i, j, k)) largest = std::max(largest, std::fabs(g.cellDivergence(i, j, k)));
    });
    Vector3 top = cellCentre(g, 0, 0, 0);
    s.forEach([&](int i, int j, int k) {
        if (g.solid(i, j, k)) return;
        const float t = 0.5f + 0.5f * g.cellDivergence(i, j, k) / largest;
        Probe::point(DrawLayer::Divergence, cellCentre(g, i, j, k), heatColor(t), 0.3f * g.dx() * float(s.stride));
        top = cellCentre(g, i, j, k);
    });
    Probe::label(DrawLayer::Divergence, top, format("div u: макс. |%.3g| 1/с", double(largest)));
}

} // namespace

float GasSolver::cellDivergence(int i, int j, int k) const {
    const float du = u_.at(i + 1, j, k) - u_.at(i, j, k);
    const float dv = v_.at(i, j + 1, k) - v_.at(i, j, k);
    const float dw = w_.at(i, j, k + 1) - w_.at(i, j, k);
    return (du + dv + dw) / dx_;
}

void GasSolver::drawDebug(int axis, int layer) const {
    const uint32_t mask = Probe::bit(DrawLayer::GasGrid) | Probe::bit(DrawLayer::GasVelocity) |
                          Probe::bit(DrawLayer::PressureGradient) | Probe::bit(DrawLayer::Divergence) |
                          Probe::bit(DrawLayer::Vorticity) | Probe::bit(DrawLayer::CurrentDensity);
    if ((Probe::layers() & mask) == 0 || nx_ * ny_ * nz_ == 0) return;
    const GasSolver& g = *this;
    const Slice s = makeSlice(g, axis, layer);
    if (Probe::layerOn(DrawLayer::GasGrid)) drawGrid(g, s);
    if (Probe::layerOn(DrawLayer::GasVelocity))
        drawArrows(DrawLayer::GasVelocity, g, s,
                   [&g](int i, int j, int k) { return g.solid(i, j, k) ? Vector3(0.0f) : g.cellVelocity(i, j, k); },
                   Vector3(0.3f, 0.9f, 1.0f), "м/с");
    if (Probe::layerOn(DrawLayer::PressureGradient))
        drawArrows(DrawLayer::PressureGradient, g, s, [&g](int i, int j, int k) { return pressureAcceleration(g, i, j, k); },
                   Vector3(1.0f, 0.55f, 0.2f), "м/с²");
    if (Probe::layerOn(DrawLayer::Divergence)) drawDivergence(g, s);
    if (Probe::layerOn(DrawLayer::Vorticity))
        drawArrows(DrawLayer::Vorticity, g, s, [&g](int i, int j, int k) { return vorticity(g, i, j, k); },
                   Vector3(1.0f, 0.5f, 1.0f), "1/с");
    if (Probe::layerOn(DrawLayer::CurrentDensity) && magnetic.enabled)
        drawArrows(DrawLayer::CurrentDensity, g, s, [&g](int i, int j, int k) { return g.magnetic.currentAt(cellCentre(g, i, j, k)); },
                   Vector3(1.0f, 0.8f, 0.2f), "А/м²");
}

} // namespace rf
