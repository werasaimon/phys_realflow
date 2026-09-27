// The SVG charts of the board, plain text written by hand (no library):
//   <case>.svg     - a refinement study. With an exact answer: the error against the step h on
//                    log-log axes, the least-squares line through the points (its slope is the
//                    observed order) and a triangle with the slope of the theoretical order.
//                    Without one: the value against h (log h) and the reference as a dashed line.
//   validation.svg - every validation case as E / D in %, with the bar E -+ u_val and the
//                    case's tolerance as a grey band around zero (zero = the reference); blinded
//                    cases are left out;
//   pulls.svg      - every case with a finite pull z = E / u_val against the 1/2/3 sigma bands.
// Colours: one data colour (blue), text and guides in greys; a light surface like the other
// figures of the docs. Every point carries a <title>: its numbers show on hover.
#include "Convergence.h"
#include "Report.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>

namespace rf::verify {

namespace {

const char* kSurface = "#fcfcfb";
const char* kInk = "#0b0b0b";
const char* kInk2 = "#52514e";
const char* kGrid = "#e5e4df";
const char* kBlue = "#2a78d6";
const char* kBand = "#f0efec";

// One axis: data range [lo, hi] onto pixels [p0, p1], linear or logarithmic.
struct Axis {
    double lo = 0, hi = 1, p0 = 0, p1 = 1;
    bool log = false;
    double map(double v) const {
        const double a = log ? std::log10(v) : v, l = log ? std::log10(lo) : lo, h = log ? std::log10(hi) : hi;
        return p0 + (a - l) / (h - l) * (p1 - p0);
    }
    // Tick values: 1, 2, 5 x 10^k inside the range (log), or ~5 round steps (linear).
    std::vector<double> ticks() const {
        std::vector<double> t;
        if (log) { // 1, 2, 5 per decade; every integer multiple when the range is too narrow for 3
            for (const std::vector<double>& ms : {std::vector<double>{1, 2, 5}, std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9}}) {
                t.clear();
                for (int k = int(std::floor(std::log10(lo))); k <= int(std::ceil(std::log10(hi))); ++k)
                    for (double m : ms) {
                        const double v = m * std::pow(10.0, k);
                        if (v >= lo && v <= hi) t.push_back(v);
                    }
                if (t.size() >= 3) break;
            }
            return t;
        }
        const double raw = (hi - lo) / 5, mag = std::pow(10.0, std::floor(std::log10(raw)));
        const double step = raw / mag < 2 ? 2 * mag : raw / mag < 5 ? 5 * mag : 10 * mag;
        for (double v = std::ceil(lo / step) * step; v <= hi + 1e-12; v += step) t.push_back(std::fabs(v) < step * 1e-9 ? 0 : v);
        return t;
    }
};

// A range padded on both sides (by a share of its width; in decades for a log axis).
Axis paddedAxis(double lo, double hi, bool log, double pad, double p0, double p1) {
    Axis a;
    a.log = log;
    a.p0 = p0;
    a.p1 = p1;
    if (log) {
        const double l = std::log10(lo), h = std::log10(hi), w = std::max(h - l, 0.3);
        a.lo = std::pow(10.0, l - pad * w);
        a.hi = std::pow(10.0, h + pad * w);
    } else {
        const double w = std::max(hi - lo, 1e-12 + 1e-3 * std::fabs(hi));
        a.lo = lo - pad * w;
        a.hi = hi + pad * w;
    }
    return a;
}

std::string f1(double v) {
    char b[32];
    std::snprintf(b, sizeof(b), "%.1f", v);
    return b;
}

std::string svgOpen(int w, int h, const std::string& title) {
    return "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 " + std::to_string(w) + " " + std::to_string(h) + "\" width=\"" +
           std::to_string(w) + "\" height=\"" + std::to_string(h) +
           "\" font-family=\"Segoe UI, Helvetica, Arial, sans-serif\" font-size=\"13\">\n<rect width=\"" + std::to_string(w) + "\" height=\"" +
           std::to_string(h) + "\" fill=\"" + kSurface + "\"/>\n<text x=\"" + std::to_string(w / 2) +
           "\" y=\"24\" text-anchor=\"middle\" font-size=\"15\" font-weight=\"600\" fill=\"" + kInk + "\">" + title + "</text>\n";
}

// Grid lines and tick labels of both axes, and the axis titles.
std::string axes(const Axis& x, const Axis& y, const std::string& xTitle, const std::string& yTitle) {
    std::string s;
    for (double v : x.ticks())
        s += "<line x1=\"" + f1(x.map(v)) + "\" y1=\"" + f1(y.p1) + "\" x2=\"" + f1(x.map(v)) + "\" y2=\"" + f1(y.p0) + "\" stroke=\"" + kGrid +
             "\"/><text x=\"" + f1(x.map(v)) + "\" y=\"" + f1(y.p0 + 18) + "\" text-anchor=\"middle\" fill=\"" + kInk2 + "\">" + num(v) + "</text>\n";
    for (double v : y.ticks())
        s += "<line x1=\"" + f1(x.p0) + "\" y1=\"" + f1(y.map(v)) + "\" x2=\"" + f1(x.p1) + "\" y2=\"" + f1(y.map(v)) + "\" stroke=\"" + kGrid +
             "\"/><text x=\"" + f1(x.p0 - 8) + "\" y=\"" + f1(y.map(v) + 4) + "\" text-anchor=\"end\" fill=\"" + kInk2 + "\">" + num(v) + "</text>\n";
    s += "<text x=\"" + f1(0.5 * (x.p0 + x.p1)) + "\" y=\"" + f1(y.p0 + 40) + "\" text-anchor=\"middle\" fill=\"" + kInk2 + "\">" + xTitle + "</text>\n";
    s += "<text transform=\"translate(18 " + f1(0.5 * (y.p0 + y.p1)) + ") rotate(-90)\" text-anchor=\"middle\" fill=\"" + kInk2 + "\">" + yTitle +
         "</text>\n";
    return s;
}

std::string dot(double x, double y, const std::string& tip) {
    return "<circle cx=\"" + f1(x) + "\" cy=\"" + f1(y) + "\" r=\"5\" fill=\"" + kBlue + "\" stroke=\"" + kSurface + "\" stroke-width=\"2\"><title>" +
           tip + "</title></circle>\n";
}

std::string line(double x1, double y1, double x2, double y2, const char* colour, const char* extra = "") {
    return "<line x1=\"" + f1(x1) + "\" y1=\"" + f1(y1) + "\" x2=\"" + f1(x2) + "\" y2=\"" + f1(y2) + "\" stroke=\"" + colour + "\" " + extra + "/>\n";
}

// The slope triangle of the theoretical order p under the finest point: a leg along h of one
// refinement ratio, the other leg the error change r^p.
std::string slopeTriangle(const Axis& x, const Axis& y, const ConvergencePoint& fine, double r, double p) {
    const double x1 = fine.h, x2 = fine.h * r, y1 = std::fabs(fine.error) * 0.35, y2 = y1 * std::pow(r, p);
    if (x2 > x.hi || y2 > y.hi || y1 < y.lo) return "";
    const std::string poly = f1(x.map(x1)) + "," + f1(y.map(y1)) + " " + f1(x.map(x2)) + "," + f1(y.map(y1)) + " " + f1(x.map(x2)) + "," + f1(y.map(y2));
    return "<polygon points=\"" + poly + "\" fill=\"none\" stroke=\"" + kInk2 + "\" stroke-dasharray=\"4 3\"/>\n<text x=\"" + f1(x.map(x2) + 6) +
           "\" y=\"" + f1(0.5 * (y.map(y1) + y.map(y2)) + 4) + "\" fill=\"" + kInk2 + "\">теория p = " + num(p) + "</text>\n";
}

// Error against h (log-log), the fitted line, the theoretical slope.
std::string errorChart(const Outcome& o) {
    const std::vector<ConvergencePoint>& pts = o.r.convergence;
    double hLo = 1e300, hHi = 0, eLo = 1e300, eHi = 0;
    for (const ConvergencePoint& p : pts) {
        hLo = std::min(hLo, p.h); hHi = std::max(hHi, p.h);
        eLo = std::min(eLo, std::fabs(p.error)); eHi = std::max(eHi, std::fabs(p.error));
    }
    const OrderFit fit = fitOrder(pts);
    const Axis x = paddedAxis(hLo, hHi, true, 0.15, 80, 600), y = paddedAxis(eLo * 0.3, eHi, true, 0.12, 330, 50);
    char order[64];
    std::snprintf(order, sizeof(order), " (p = %.2f)", fit.order);
    std::string s = svgOpen(640, 400, o.c->title + order);
    s += axes(x, y, o.r.hLabel + " (лог.)", "ошибка (лог.)");
    if (std::isfinite(fit.order)) {
        auto yFit = [&](double h) { return std::exp(fit.lnC + fit.order * std::log(h)); };
        s += line(x.map(hLo), y.map(yFit(hLo)), x.map(hHi), y.map(yFit(hHi)), kBlue, "stroke-width=\"2\"");
    }
    for (const ConvergencePoint& p : pts) s += dot(x.map(p.h), y.map(std::fabs(p.error)), "h = " + num(p.h) + ", ошибка " + num(p.error));
    const ConvergencePoint& fine = *std::min_element(pts.begin(), pts.end(), [](auto& a, auto& b) { return a.h < b.h; });
    if (std::isfinite(o.r.theoreticalOrder)) s += slopeTriangle(x, y, fine, 2.0, o.r.theoreticalOrder);
    return s + "</svg>\n";
}

// Value against h (log h), the reference as a dashed line - none for a blinded case.
std::string valueChart(const Outcome& o) {
    const std::vector<ConvergencePoint>& pts = o.r.convergence;
    const double D = o.ref.value;
    double hLo = 1e300, hHi = 0, vLo = std::isfinite(D) ? D : 1e300, vHi = std::isfinite(D) ? D : -1e300;
    for (const ConvergencePoint& p : pts) {
        hLo = std::min(hLo, p.h); hHi = std::max(hHi, p.h);
        if (std::isfinite(p.value)) { vLo = std::min(vLo, p.value); vHi = std::max(vHi, p.value); }
    }
    const Axis x = paddedAxis(hLo, hHi, true, 0.15, 80, 600), y = paddedAxis(vLo, vHi, false, 0.2, 330, 50);
    std::string s = svgOpen(640, 400, o.c->title + ": сходимость" + (std::isfinite(D) ? "" : " (эталон запечатан)"));
    s += axes(x, y, o.r.hLabel + " (лог.)", "значение" + (o.r.unit.empty() ? "" : ", " + o.r.unit));
    if (std::isfinite(D)) {
        const double yr = y.map(D);
        s += line(x.p0, yr, x.p1, yr, kInk2, "stroke-dasharray=\"6 4\"");
        s += "<text x=\"" + f1(x.p1 - 4) + "\" y=\"" + f1(yr - 6) + "\" text-anchor=\"end\" fill=\"" + kInk2 + "\">эталон " + num(D) + "</text>\n";
    }
    for (size_t i = 1; i < pts.size(); ++i)
        if (std::isfinite(pts[i].value) && std::isfinite(pts[i - 1].value))
            s += line(x.map(pts[i - 1].h), y.map(pts[i - 1].value), x.map(pts[i].h), y.map(pts[i].value), kBlue, "stroke-width=\"2\"");
    for (const ConvergencePoint& p : pts)
        if (std::isfinite(p.value)) s += dot(x.map(p.h), y.map(p.value), "h = " + num(p.h) + ", значение " + num(p.value));
    return s + "</svg>\n";
}

// Every validation case as a row: E / D in %, the bar E -+ u_val, the tolerance band.
std::string validationChart(const std::vector<const Outcome*>& rows) {
    const int h = 90 + 56 * int(rows.size());
    double lim = 5;
    for (const Outcome* o : rows) {
        const double D = std::fabs(o->ref.value);
        if (D > 0 && std::isfinite(o->a.E)) lim = std::max(lim, 100 * (std::fabs(o->a.E) + o->a.uVal) / D);
        if (o->c->acceptance.relative) lim = std::max(lim, 100 * o->c->acceptance.warnTolerance);
    }
    Axis x = paddedAxis(-lim, lim, false, 0.08, 330, 620); // the row titles need ~300 px on the left
    Axis y;
    y.p0 = h - 50;
    y.p1 = 44;
    std::string s = svgOpen(640, h, "Валидация: E = S − D в % от эталона, полоса ± u_val");
    for (double v : x.ticks())
        s += line(x.map(v), y.p1, x.map(v), y.p0, kGrid) + "<text x=\"" + f1(x.map(v)) + "\" y=\"" + f1(y.p0 + 18) +
             "\" text-anchor=\"middle\" fill=\"" + kInk2 + "\">" + num(v) + " %</text>\n";
    s += line(x.map(0), y.p1, x.map(0), y.p0, kInk2) + "<text x=\"" + f1(x.map(0)) + "\" y=\"" + f1(y.p1 - 4) +
         "\" text-anchor=\"middle\" fill=\"" + kInk2 + "\">эталон</text>\n";
    for (size_t i = 0; i < rows.size(); ++i) {
        const Outcome& o = *rows[i];
        const double yc = y.p1 + 30 + 56 * double(i), D = std::fabs(o.ref.value);
        if (o.c->acceptance.relative) {
            const double t = 100 * o.c->acceptance.tolerance;
            s += "<rect x=\"" + f1(x.map(-t)) + "\" y=\"" + f1(yc - 14) + "\" width=\"" + f1(x.map(t) - x.map(-t)) + "\" height=\"28\" fill=\"" + kBand +
                 "\"><title>допуск ± " + num(t) + " %</title></rect>\n";
        }
        s += "<text x=\"318\" y=\"" + f1(yc + 4) + "\" text-anchor=\"end\" fill=\"" + kInk + "\">" + o.c->title + "</text>\n";
        if (D == 0 || !std::isfinite(o.a.E)) continue;
        const double e = 100 * o.a.E / D, u = 100 * o.a.uVal / D;
        s += line(x.map(e - u), yc, x.map(e + u), yc, kBlue, "stroke-width=\"2\"");
        s += line(x.map(e - u), yc - 6, x.map(e - u), yc + 6, kBlue, "stroke-width=\"2\"") + line(x.map(e + u), yc - 6, x.map(e + u), yc + 6, kBlue, "stroke-width=\"2\"");
        s += dot(x.map(e), yc, "E = " + num(o.a.E) + " (" + num(e) + " %), u_val = " + num(o.a.uVal));
    }
    return s + "</svg>\n";
}

// The pull plot: one row per case with a finite pull, z on a common axis from -6 to 6 sigma, the
// 1, 2 and 3 sigma bands in greys (darkest innermost); a pull beyond the axis sits at its edge with
// its number printed.
std::string pullChart(const std::vector<const Outcome*>& rows) {
    const int h = 90 + 26 * int(rows.size());
    Axis x;
    x.lo = -6; x.hi = 6; x.p0 = 330; x.p1 = 620;
    const double top = 44, bottom = h - 50;
    std::string s = svgOpen(640, h, "Пулы z = E / u_val, полосы 1, 2, 3 σ");
    const char* bands[3] = {"#f0efec", "#e5e4df", "#d9d8d2"};
    for (int k = 3; k >= 1; --k)
        s += "<rect x=\"" + f1(x.map(-k)) + "\" y=\"" + f1(top) + "\" width=\"" + f1(x.map(k) - x.map(-k)) + "\" height=\"" + f1(bottom - top) +
             "\" fill=\"" + bands[3 - k] + "\"><title>± " + std::to_string(k) + " σ</title></rect>\n";
    for (int v = -6; v <= 6; v += 2)
        s += "<text x=\"" + f1(x.map(v)) + "\" y=\"" + f1(bottom + 18) + "\" text-anchor=\"middle\" fill=\"" + kInk2 + "\">" + std::to_string(v) + " σ</text>\n";
    s += line(x.map(0), top, x.map(0), bottom, kInk2);
    for (size_t i = 0; i < rows.size(); ++i) {
        const Outcome& o = *rows[i];
        const double yc = top + 20 + 26 * double(i), z = std::clamp(o.a.z, -6.0, 6.0);
        s += "<text x=\"318\" y=\"" + f1(yc + 4) + "\" text-anchor=\"end\" fill=\"" + kInk + "\" font-size=\"12\">" + o.c->title + "</text>\n";
        s += dot(x.map(z), yc, o.c->id + ": z = " + num(o.a.z));
        if (z != o.a.z) s += "<text x=\"" + f1(x.map(z) + (z > 0 ? -10 : 10)) + "\" y=\"" + f1(yc - 8) + "\" text-anchor=\"" + (z > 0 ? "end" : "start") +
                             "\" fill=\"" + kInk2 + "\" font-size=\"11\">" + num(o.a.z) + " σ</text>\n";
    }
    return s + "</svg>\n";
}

bool writeText(const std::string& path, const std::string& text) {
    std::ofstream f(path, std::ios::binary);
    f << text;
    return bool(f);
}

} // namespace

std::vector<std::string> writeCharts(const std::string& dir, const std::vector<Outcome>& outcomes) {
    std::vector<std::string> written;
    std::vector<const Outcome*> validation, pulls;
    for (const Outcome& o : outcomes) {
        if (o.c->category == "validation" && o.status != BlindStatus::Blinded) validation.push_back(&o);
        if (std::isfinite(o.a.z)) pulls.push_back(&o);
        if (o.r.convergence.size() < 2) continue;
        bool errors = true;
        for (const ConvergencePoint& p : o.r.convergence) errors &= std::isfinite(p.error) && std::fabs(p.error) > 0;
        if (writeText(dir + "/" + o.c->id + ".svg", errors ? errorChart(o) : valueChart(o))) written.push_back(o.c->id + ".svg");
    }
    if (!validation.empty() && writeText(dir + "/validation.svg", validationChart(validation))) written.push_back("validation.svg");
    if (writeText(dir + "/pulls.svg", pullChart(pulls))) written.push_back("pulls.svg");
    return written;
}

} // namespace rf::verify
