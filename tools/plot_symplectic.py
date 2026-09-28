#!/usr/bin/env python3
"""The plot of docs/08-relativity.md, "Симплектические геодесические", drawn from a real run of
tests/SymplecticTests.cpp: the largest energy error so far against the number of orbits, for RK4
and for two symplectic methods, on an eccentric inclined orbit around a Kerr hole (a = 0.9 M).

    RF_LONG_ORBITS=1000000 RF_PLOT_DIR=plots RF_TEST="symplectic geodesics: energy" ./rf_tests
    python tools/plot_symplectic.py plots docs/img          (writes docs/img/gr-symplectic-drift.svg)

On log-log axes a drift that grows linearly in time is a line of slope 1; a bounded wobble is flat.
RK4 at the step of the symplectic methods is drawn thin and dashed; RK4 at a third of that step
costs the same as Tao-4 (12 gradients per unit of time) and is the fair comparison. Python standard
library only; the axes, grid and both-theme palette are those of tools/plot_action.py.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_action import H, LEFT, RIGHT, BOTTOM, TOP, W, Axis, frame, grid, legend, log_range, polyline, read_csv  # noqa: E402


def decade_ticks(lo, hi):
    """Powers of ten inside [lo, hi]: many decades need nothing finer."""
    return [10.0 ** e for e in range(math.ceil(math.log10(lo)), math.floor(math.log10(hi)) + 1)]


def crossing(xs, rising, flat):
    """The first x where the rising curve passes the flat one (linear in log-log between samples)."""
    for i in range(1, len(xs)):
        if rising[i - 1] < flat[i - 1] and rising[i] >= flat[i]:
            a, b = math.log(rising[i - 1] / flat[i - 1]), math.log(rising[i] / flat[i])
            t = a / (a - b)
            return math.exp(math.log(xs[i - 1]) + t * (math.log(xs[i]) - math.log(xs[i - 1])))
    return None


def main():
    src, dst = sys.argv[1], sys.argv[2]
    os.makedirs(dst, exist_ok=True)
    d = read_csv(src, "symplectic_drift")
    orbits = d["orbits"]
    series = [("RK4 той же цены (h = M/3)", d["rk4_same_cost"]), ("неявная средняя точка (h = M)", d["implicit_midpoint"]),
              ("Тао-4 (h = M, ω = 0,1)", d["tao4"])]
    values = [v for _, s in series for v in s if v > 0] + [v for v in d["rk4_h1"] if v > 0]
    xlo, xhi, _ = log_range(min(orbits), max(orbits))
    ylo, yhi, _ = log_range(min(values), max(values))
    xa = Axis(xlo, xhi, LEFT, W - RIGHT, log=True)
    ya = Axis(ylo, yhi, H - BOTTOM, TOP, log=True)
    out = frame(f"Дрейф энергии за {max(orbits):,.0f} витков: RK4 против симплектических".replace(",", " "),
                "витки вокруг дыры", "наибольшая |H − H₀| к этому витку",
                "Керр a = 0,9 M, орбита r = 8…15 M, Q = 4 M². Тест «symplectic geodesics: energy…» "
                "(RF_LONG_ORBITS=… RF_PLOT_DIR=… rf_tests, tests/SymplecticTests.cpp)")
    grid(out, xa, ya, decade_ticks(xlo, xhi), decade_ticks(ylo, yhi))
    # RK4 at the symplectic step: thin, dashed, muted - it leaves the frame of the others quickly.
    pts = " ".join(f"{xa(x):.1f},{ya(y):.1f}" for x, y in zip(orbits, d["rk4_h1"]) if y > 0)
    out.append(f"<polyline points='{pts}' fill='none' stroke='var(--muted)' stroke-width='1.3' stroke-dasharray='4 4'/>")
    # Its label in the middle of the line, just left of it, where nothing else is drawn.
    mid = min(range(len(orbits)), key=lambda i: abs(math.log(orbits[i] / 300.0)))
    out.append(f"<text x='{xa(orbits[mid]) - 8:.1f}' y='{ya(d['rk4_h1'][mid]) - 8:.1f}' text-anchor='end' class='cap'>"
               "RK4, h = M: уходит линейно</text>")
    for i, (_, s) in enumerate(series):
        polyline(out, xa, ya, orbits, s, i)
    # Where RK4 of the same cost passes Tao-4: from there on the symplectic method is ahead, for good.
    # Inside the run: a dotted vertical line. Beyond it: RK4's error grows linearly, so the crossing
    # is extrapolated from the last point and written at the right edge.
    x = crossing(orbits, d["rk4_same_cost"], d["tao4"])
    if x:
        out.append(f"<line x1='{xa(x):.1f}' y1='{TOP}' x2='{xa(x):.1f}' y2='{H - BOTTOM}' stroke='var(--muted)' stroke-dasharray='2 4'/>")
        out.append(f"<text x='{xa(x) - 6:.1f}' y='{TOP + 14}' text-anchor='end' class='cap'>≈ {x:,.0f} витков</text>".replace(",", " "))
    else:  # in the empty lower right corner, under the rising line
        ahead = orbits[-1] * d["tao4"][-1] / d["rk4_same_cost"][-1]
        text = f"RK4 той же цены догонит Тао-4 ≈ на {ahead / 1e6:.2f} млн витков".replace(".", ",")
        out.append(f"<text x='{W - RIGHT - 4}' y='{H - BOTTOM - 30}' text-anchor='end' class='cap'>{text}</text>")
        out.append(f"<text x='{W - RIGHT - 4}' y='{H - BOTTOM - 14}' text-anchor='end' class='cap'>"
                   "(его ошибка растёт линейно, у Тао-4 — нет)</text>")
    legend(out, [n for n, _ in series])
    out.append("</svg>")
    path = os.path.join(dst, "gr-symplectic-drift.svg")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    print("plot written to", path)


if __name__ == "__main__":
    main()
