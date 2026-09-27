#!/usr/bin/env python3
"""The plot of docs/08-relativity.md, "Кривизна из метрики", drawn from a real run of
tests/CurvatureTests.cpp: the error of the Christoffel symbols and of the Kretschmann scalar
against the finite-difference step, on Schwarzschild at r = 7M.

    RF_PLOT_DIR=plots RF_TEST="curvature: Schwarzschild" ./rf_tests   (writes plots/curvature_step.csv)
    python tools/plot_curvature.py plots docs/img                      (writes docs/img/curvature-step.svg)

The error first falls as h^4 (the truncation of the 4th-order differences), then rises as 1/h^2
(the rounding of the nested differences the Riemann tensor needs): the default step 1e-3 sits at
the bottom of the V. Python standard library only; the drawing helpers (axes, grid, palette in
both themes) are those of tools/plot_action.py.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_action import H, LEFT, RIGHT, BOTTOM, TOP, W, Axis, frame, grid, legend, log_range, polyline, read_csv  # noqa: E402


def main():
    src, dst = sys.argv[1], sys.argv[2]
    os.makedirs(dst, exist_ok=True)
    d = read_csv(src, "curvature_step")
    h = d["step_fraction"]
    series = [("символы Кристоффеля Γ", d["gamma_rel_err"]), ("инвариант Кречмана K = R_abcd R^abcd", d["kretschmann_rel_err"])]
    vals = [v for _, s in series for v in s if v > 0]
    xlo, xhi, xt = log_range(min(h), max(h))
    ylo, yhi, yt = log_range(min(vals), max(vals))
    yt = [t for t in yt if abs(math.log10(t) - round(math.log10(t))) < 1e-9]  # powers of ten only: 7 decades
    xa = Axis(xlo, xhi, LEFT, W - RIGHT, log=True)
    ya = Axis(ylo, yhi, H - BOTTOM, TOP, log=True)
    out = frame("Кривизна из метрики: ошибка против шага разности", "шаг h, доля размера координаты",
                "относительная ошибка",
                "Шварцшильд, r = 7M. Тест «curvature: Schwarzschild…» (RF_PLOT_DIR=… rf_tests, tests/CurvatureTests.cpp)")
    grid(out, xa, ya, xt, yt)
    for i, (_, s) in enumerate(series):
        polyline(out, xa, ya, h, s, i)
        for x, y in zip(h, s):
            out.append(f"<circle cx='{xa(x):.1f}' cy='{ya(y):.1f}' r='4.5' fill='var(--s{i})' stroke='var(--surface)' stroke-width='2'/>")
    # The slope-4 triangle of the truncation over the two largest steps, and the default step.
    x1, x2 = h[1], h[0]
    y1 = d["kretschmann_rel_err"][1] * 0.3
    y2 = y1 * (x2 / x1) ** 4
    pts = f"{xa(x1):.1f},{ya(y1):.1f} {xa(x2):.1f},{ya(y1):.1f} {xa(x2):.1f},{ya(y2):.1f}"
    out.append(f"<polygon points='{pts}' fill='none' stroke='var(--muted)' stroke-width='1.2'/>")
    out.append(f"<text x='{xa(x2) - 6:.1f}' y='{(ya(y1) + ya(y2)) / 2 + 4:.1f}' text-anchor='end' class='cap'>наклон 4</text>")
    xd = xa(1e-3)
    out.append(f"<line x1='{xd:.1f}' y1='{TOP}' x2='{xd:.1f}' y2='{H - BOTTOM}' stroke='var(--muted)' stroke-dasharray='3 4'/>")
    out.append(f"<text x='{xd + 6:.1f}' y='{TOP + 14}' class='cap'>шаг по умолчанию 10⁻³</text>")
    legend(out, [n for n, _ in series])
    out.append("</svg>")
    path = os.path.join(dst, "curvature-step.svg")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    print("plot written to", path)


if __name__ == "__main__":
    main()
