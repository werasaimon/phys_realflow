#!/usr/bin/env python3
"""The plots of docs/11-action.md, drawn from real runs of tests/ActionTests.cpp.

Run the tests with a folder for their curves, then this script turns every CSV into an SVG:

    RF_PLOT_DIR=plots RF_TEST="action:" ./rf_tests            (the tests write plots/*.csv)
    python tools/plot_action.py plots docs/images/action      (this script writes the SVGs)

Python standard library only (no matplotlib), so the plots rebuild anywhere. Every SVG carries
its own light and dark palette (prefers-color-scheme), axes with units, a legend, and a caption
naming the test that produced the numbers. Nothing here is drawn by hand: every point is a row
of a CSV written by a test.
"""
import csv
import math
import os
import sys

# The validated categorical palette (light, dark) of the plots, by series slot, and the inks.
SERIES = [("#2a78d6", "#3987e5"), ("#eb6834", "#d95926"), ("#1baf7a", "#199e70")]
INK = {"surface": ("#fcfcfb", "#1a1a19"), "text": ("#0b0b0b", "#ffffff"), "muted": ("#52514e", "#c3c2b7"),
       "grid": ("#e4e3df", "#33332f"), "axis": ("#8a8983", "#6b6a64")}
DASHES = ["", "7 4", "2 3"]  # the second encoding of a series, beside its colour

W, H = 760, 430
LEFT, RIGHT, TOP, BOTTOM = 84, 24, 64, 86


def read_csv(folder, name):
    with open(os.path.join(folder, name + ".csv"), newline="") as f:
        rows = list(csv.reader(f))
    header, data = rows[0], [[float(x) for x in r] for r in rows[1:] if r]
    return {h: [r[i] for r in data] for i, h in enumerate(header)}


def nice_step(span, count):
    """A round tick step (1, 2 or 5 times a power of ten) giving about `count` ticks."""
    raw = span / max(count, 1)
    power = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 5, 10):
        if m * power >= raw:
            return m * power
    return 10 * power


def linear_ticks(lo, hi, count=5):
    step = nice_step(hi - lo, count)
    first = math.ceil(lo / step - 1e-9) * step
    ticks, t = [], first
    while t <= hi + 1e-9 * step:
        ticks.append(round(t, 12))
        t += step
    return ticks


def fmt(v):
    if v == 0:
        return "0"
    if abs(v) >= 1e4 or abs(v) < 1e-3:
        m, e = f"{v:.0e}".split("e")
        return f"{m}·10<tspan dy='-6' font-size='10'>{int(e)}</tspan>"
    text = f"{v:.6g}"
    return text.replace(".", ",").replace("-", "−")


def style():
    light = "".join(f"--{k}:{v[0]};" for k, v in INK.items()) + "".join(f"--s{i}:{c[0]};" for i, c in enumerate(SERIES))
    dark = "".join(f"--{k}:{v[1]};" for k, v in INK.items()) + "".join(f"--s{i}:{c[1]};" for i, c in enumerate(SERIES))
    return (f"<style>svg{{{light}}}@media (prefers-color-scheme: dark){{svg{{{dark}}}}}"
            "text{font-family:'Segoe UI',system-ui,-apple-system,Roboto,sans-serif;fill:var(--text)}"
            ".muted{fill:var(--muted)}.tick{font-size:12px;fill:var(--muted);font-variant-numeric:tabular-nums}"
            ".title{font-size:16px;font-weight:600}.cap{font-size:11.5px;fill:var(--muted)}</style>")


class Axis:
    """A linear or logarithmic map from data to pixels."""

    def __init__(self, lo, hi, a, b, log=False):
        self.log = log
        self.lo, self.hi = (math.log10(lo), math.log10(hi)) if log else (lo, hi)
        self.a, self.b = a, b

    def __call__(self, v):
        v = math.log10(v) if self.log else v
        return self.a + (v - self.lo) / (self.hi - self.lo) * (self.b - self.a)


def frame(title, xlabel, ylabel, caption):
    out = [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {W} {H}' width='{W}' height='{H}' role='img'>",
           f"<title>{title}</title>", style(),
           f"<rect width='{W}' height='{H}' rx='10' fill='var(--surface)'/>",
           f"<text x='{LEFT}' y='28' class='title'>{title}</text>",
           f"<text x='{(LEFT + W - RIGHT) / 2}' y='{H - BOTTOM + 38}' text-anchor='middle' class='muted' font-size='13'>{xlabel}</text>",
           f"<text transform='translate(20 {(TOP + H - BOTTOM) / 2}) rotate(-90)' text-anchor='middle' class='muted' font-size='13'>{ylabel}</text>",
           ]
    # The caption on one or two lines (about 110 characters fit the width).
    words, lines = caption.split(" "), [""]
    for word in words:
        if len(lines[-1]) + len(word) + 1 > 110 and lines[-1]:
            lines.append("")
        lines[-1] = (lines[-1] + " " + word).strip()
    for i, line in enumerate(lines[-2:]):
        out.append(f"<text x='{LEFT}' y='{H - 26 + 15 * i if len(lines) > 1 else H - 14}' class='cap'>{line}</text>")
    return out


def grid(out, xa, ya, xticks, yticks):
    x0, x1, y0, y1 = LEFT, W - RIGHT, H - BOTTOM, TOP
    for t in yticks:
        y = ya(t)
        out.append(f"<line x1='{x0}' x2='{x1}' y1='{y:.1f}' y2='{y:.1f}' stroke='var(--grid)' stroke-width='1'/>")
        out.append(f"<text x='{x0 - 8}' y='{y + 4:.1f}' text-anchor='end' class='tick'>{fmt(t)}</text>")
    for t in xticks:
        x = xa(t)
        out.append(f"<line x1='{x:.1f}' x2='{x:.1f}' y1='{y0}' y2='{y0 + 5}' stroke='var(--axis)' stroke-width='1'/>")
        out.append(f"<text x='{x:.1f}' y='{y0 + 20}' text-anchor='middle' class='tick'>{fmt(t)}</text>")
    out.append(f"<line x1='{x0}' x2='{x1}' y1='{y0}' y2='{y0}' stroke='var(--axis)' stroke-width='1'/>")
    out.append(f"<line x1='{x0}' x2='{x0}' y1='{y0}' y2='{y1}' stroke='var(--axis)' stroke-width='1'/>")


def legend(out, names):
    x = LEFT
    for i, name in enumerate(names):
        dash = f" stroke-dasharray='{DASHES[i]}'" if DASHES[i] else ""
        out.append(f"<line x1='{x:.1f}' x2='{x + 26:.1f}' y1='47' y2='47' stroke='var(--s{i})' stroke-width='2.5'{dash}/>")
        out.append(f"<text x='{x + 32:.1f}' y='51' font-size='12.5'>{name}</text>")
        x += 44 + 7.2 * len(name)


def polyline(out, xa, ya, xs, ys, slot):
    step = max(1, math.ceil(len(xs) / 300))  # at most ~300 points a line: finer than the pixels
    keep = list(range(0, len(xs), step)) + ([len(xs) - 1] if (len(xs) - 1) % step else [])
    pts = " ".join(f"{xa(xs[i]):.1f},{ya(ys[i]):.1f}" for i in keep)
    dash = f" stroke-dasharray='{DASHES[slot]}'" if DASHES[slot] else ""
    out.append(f"<polyline points='{pts}' fill='none' stroke='var(--s{slot})' stroke-width='2' "
               f"stroke-linejoin='round' stroke-linecap='round'{dash}/>")


def line_chart(path, title, xlabel, ylabel, caption, x, series, scale=1.0):
    """Lines over a shared x; `series` = [(name, values)], all in the unit of ylabel."""
    ys_all = [v * scale for _, s in series for v in s]
    lo, hi = min(0.0, min(ys_all)), max(ys_all)
    yticks = linear_ticks(lo, hi)
    lo, hi = min(lo, yticks[0]), max(hi, yticks[-1])
    xticks = linear_ticks(min(x), max(x))
    xa = Axis(min(x), max(xticks[-1], max(x)), LEFT, W - RIGHT)
    ya = Axis(lo, hi, H - BOTTOM, TOP)
    out = frame(title, xlabel, ylabel, caption)
    grid(out, xa, ya, xticks, yticks)
    for i, (_, s) in enumerate(series):
        polyline(out, xa, ya, x, [v * scale for v in s], i)
    legend(out, [n for n, _ in series])
    out.append("</svg>")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")


def log_range(lo, hi):
    """The axis range around the data (a third of an octave of margin) and its 1-2-5 ticks inside."""
    lo, hi = lo / 1.25, hi * 1.25
    ticks = [m * 10 ** e for e in range(math.floor(math.log10(lo)), math.ceil(math.log10(hi)) + 1) for m in (1, 2, 5)]
    return lo, hi, [t for t in ticks if lo <= t <= hi]


def convergence_chart(path, title, caption, dt, series, order):
    """Log-log drift against the step, with the triangle of the expected slope `order`."""
    vals = [v for _, s in series for v in s if v > 0]
    xlo, xhi, xt = log_range(min(dt), max(dt))
    ylo, yhi, yt = log_range(min(vals), max(vals))
    xa = Axis(xlo, xhi, LEFT, W - RIGHT, log=True)
    ya = Axis(ylo, yhi, H - BOTTOM, TOP, log=True)
    out = frame(title, "шаг по времени Δt, с", "наибольший относительный дрейф за 5 с", caption)
    grid(out, xa, ya, xt, yt)
    for i, (_, s) in enumerate(series):
        polyline(out, xa, ya, dt, s, i)
        for x, y in zip(dt, s):
            out.append(f"<circle cx='{xa(x):.1f}' cy='{ya(y):.1f}' r='4.5' fill='var(--s{i})' stroke='var(--surface)' stroke-width='2'/>")
    # The slope triangle under the first series: one decade in dt, `order` decades in the drift.
    s0 = series[0][1]
    xr, yr = dt[1], s0[1] * 0.45
    x2 = xr * 2
    y2 = yr * 2 ** order
    pts = f"{xa(xr):.1f},{ya(yr):.1f} {xa(x2):.1f},{ya(yr):.1f} {xa(x2):.1f},{ya(y2):.1f}"
    out.append(f"<polygon points='{pts}' fill='none' stroke='var(--muted)' stroke-width='1.2'/>")
    out.append(f"<text x='{xa(x2) + 6:.1f}' y='{(ya(yr) + ya(y2)) / 2 + 4:.1f}' class='cap'>наклон {order}</text>")
    legend(out, [n for n, _ in series])
    out.append("</svg>")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")


def main():
    src, dst = sys.argv[1], sys.argv[2]
    os.makedirs(dst, exist_ok=True)
    run = "RF_PLOT_DIR=… rf_tests, tests/ActionTests.cpp"
    d = read_csv(src, "friction")
    line_chart(os.path.join(dst, "friction.svg"), "Трение: потерянная энергия = работа трения", "время, с", "энергия, Дж",
               f"Брусок 32 кг, v₀ = 3 м/с, μ = 0,6. Тест «action: friction…» ({run})", d["t_s"],
               [("потеряно энергии", d["energy_lost_J"]), ("μ m g × путь", d["friction_work_J"])])
    d = read_csv(src, "damped_pendulum")
    line_chart(os.path.join(dst, "damped-pendulum.svg"), "Маятник с затуханием: сустав работы не совершает", "время, с",
               "энергия, Дж", f"Маятник 1 м на шаровом суставе, c = 0,3 1/с. Тест «action: a damped pendulum…» ({run})",
               d["t_s"], [("потеряно энергии", d["energy_lost_J"]), ("∫ 2R dt (затухание)", d["damping_work_J"])])
    d = read_csv(src, "viscous")
    line_chart(os.path.join(dst, "viscous.svg"), "Вихрь Тейлора–Грина: вязкость и численная диссипация", "время, с",
               "энергия, мДж", f"32 ячейки, ν = 0,01 м²/с. Тест «action: viscous dissipation…» ({run})", d["t_s"],
               [("потеряно решателем", d["energy_lost_J"]), ("∫ 2νρ|S|² dV dt", d["viscous_work_J"]),
                ("точно: K₀(1 − e^(−4νπ²t))", d["exact_loss_J"])], scale=1e3)
    d = read_csv(src, "joule")
    line_chart(os.path.join(dst, "joule.svg"), "Магнитная энергия уходит в джоулево тепло", "время, мс", "энергия, мДж",
               f"Bz = B₁cos(πx), σ = 10⁵ См/м, стенки — идеальные проводники. Тест «action: magnetic energy…» ({run})",
               d["t_ms"], [("потеряно магнитной энергии", d["magnetic_energy_lost_uJ"]),
                           ("∫ J²/σ по свободным рёбрам", d["joule_heat_uJ"]), ("jouleHeating() по ячейкам", d["cell_diagnostic_uJ"])], scale=1e-3)
    d = read_csv(src, "free_rotation")
    line_chart(os.path.join(dst, "free-rotation-drift.svg"), "Свободное вращение: дрейф инвариантов Нётер", "время, с",
               "относительный дрейф, ×10⁻³", f"Брусок вращается вокруг средней оси, Δt = 1/600 с. Тест «action: a free spinning box…» ({run})",
               d["t_s"], [("энергия E/E₀ − 1", d["energy_drift"]), ("|L − L₀| / |L₀|", d["momentum_drift"])], scale=1e3)
    line_chart(os.path.join(dst, "free-rotation-omega.svg"), "Эффект Джанибекова: угловая скорость в осях тела", "время, с",
               "ω, рад/с", f"Та же прогонка: ω_y меняет знак — тело переворачивается. Тест «action: a free spinning box…» ({run})",
               d["t_s"], [("ω_x", d["wx"]), ("ω_y (средняя ось)", d["wy"]), ("ω_z", d["wz"])])
    d = read_csv(src, "free_rotation_convergence")
    convergence_chart(os.path.join(dst, "free-rotation-convergence.svg"), "Дрейф L сходится с первым порядком по Δt",
                      f"5 с свободного вращения при четырёх шагах. Тест «action: a free spinning box…» ({run})", d["dt_s"],
                      [("|L − L₀| / |L₀|", d["momentum_drift"]), ("|E/E₀ − 1|", d["energy_drift"])], order=1)
    print("plots written to", dst)


if __name__ == "__main__":
    main()
