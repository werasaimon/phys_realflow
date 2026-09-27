#!/usr/bin/env python3
"""The pictures of the book "Язык природы" (docs/math/), drawn from real runs of the engine.

Every picture of the book is a CSV written by one test of tests/MathBookTests.cpp. Run the tests
with a folder for their curves, then this script turns every CSV into an SVG:

    RF_PLOT_DIR=plots RF_TEST=mathbook ./rf_tests          (the tests write plots/*.csv)
    python tools/plot_mathbook.py plots docs/img/math        (this script writes the SVGs)

The drawing tools - axes, ticks, legend, the light and dark palette - are those of
tools/plot_action.py, so every plot of the documentation looks the same. Python standard library
only. Nothing here is drawn by hand: every point is a row of a CSV written by a test.
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_action import (W, H, LEFT, RIGHT, TOP, BOTTOM, Axis, frame, grid, legend, linear_ticks,  # noqa: E402
                         log_range, polyline, read_csv)

RUN = "RF_PLOT_DIR=… rf_tests, tests/MathBookTests.cpp"


def axis_range(values, log):
    """The axis range around the data and its ticks.

    Logarithmic: 1-2-5 in every decade while the data spans a few decades; over a wide span only
    the powers of ten, and every second or third of them, so the labels never crowd.
    """
    if not log:
        ticks = linear_ticks(min(values), max(values))
        return min(min(values), ticks[0]), max(max(values), ticks[-1]), ticks
    lo, hi, ticks = log_range(min(v for v in values if v > 0), max(values))
    decades = math.log10(hi / lo)
    if decades > 4:
        every = 1 if decades <= 8 else 2 if decades <= 14 else 3
        first, last = math.ceil(math.log10(lo)), math.floor(math.log10(hi))
        ticks = [10.0 ** e for e in range(first, last + 1) if (e - first) % every == 0]
    return lo, hi, ticks


def chart(name, dst, title, xlabel, ylabel, caption, series, xlog=False, ylog=False, dots=False, square=False):
    """One chart. `series` = [(legend name, xs, ys)], at most three.

    1. the ranges of both axes cover every series (logarithmic where asked);
    2. with `square` one metre is as long across as it is up, so a circle looks round;
    3. the grid, the lines (and dots), the legend.
    """
    xs = [x for _, sx, _ in series for x in sx]
    ys = [y for _, _, sy in series for y in sy]
    xlo, xhi, xticks = axis_range(xs, xlog)
    ylo, yhi, yticks = axis_range(ys, ylog)
    if square:  # widen the x range until a metre is as long across as it is up
        metres_across = (yhi - ylo) * (W - RIGHT - LEFT) / (H - BOTTOM - TOP)
        middle = (xlo + xhi) / 2
        xlo, xhi = middle - metres_across / 2, middle + metres_across / 2
        xticks = [t for t in linear_ticks(xlo, xhi) if xlo <= t <= xhi]
    xa = Axis(xlo, xhi, LEFT, W - RIGHT, log=xlog)
    ya = Axis(ylo, yhi, H - BOTTOM, TOP, log=ylog)
    out = frame(title, xlabel, ylabel, caption)
    grid(out, xa, ya, xticks, yticks)
    for slot, (_, sx, sy) in enumerate(series):
        keep = [(x, y) for x, y in zip(sx, sy) if (not xlog or x > 0) and (not ylog or y > 0)]
        polyline(out, xa, ya, [x for x, _ in keep], [y for _, y in keep], slot)
        if dots:
            for x, y in keep:
                out.append(f"<circle cx='{xa(x):.1f}' cy='{ya(y):.1f}' r='4' fill='var(--s{slot})' "
                           f"stroke='var(--surface)' stroke-width='1.5'/>")
    legend(out, [n for n, _, _ in series])
    out.append("</svg>")
    with open(os.path.join(dst, name + ".svg"), "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")


def circle(radius, centre_y, count=200):
    """The exact circle of the Lorentz chapter: x = r sin(a), y = centre + r cos(a)."""
    angles = [2 * math.pi * k / count for k in range(count + 1)]
    return [radius * math.sin(a) for a in angles], [centre_y + radius * math.cos(a) for a in angles]


def main():
    src, dst = sys.argv[1], sys.argv[2]
    os.makedirs(dst, exist_ok=True)

    d = read_csv(src, "rounding")
    chart("rounding", dst, "0,1 + 0,1 + … десять миллионов раз", "сколько раз сложили", "относительная ошибка суммы",
          f"float хранит ~7 значащих цифр, double ~16. Тест «mathbook: float and double…» ({RUN})",
          [("float", d["n"], d["float_relative_error"]), ("double", d["n"], d["double_relative_error"])],
          xlog=True, ylog=True, dots=True)

    d = read_csv(src, "lorentz")
    cx, cy = circle(1.0, -1.0)
    chart("lorentz", dst, "Заряд в магнитном поле идёт по кругу", "x, м", "y, м",
          f"F = q v × B, радиус m v / (q B) = 1 м, шаг Бориса h = 0,05 с, 200 шагов. Тест «mathbook: vectors…» ({RUN})",
          [("путь заряда (движок)", d["x_m"], d["y_m"]), ("точный круг", cx, cy)], square=True)

    d = read_csv(src, "derivative")
    chart("derivative", dst, "Наклон высоты — это скорость", "время, с", "скорость, м/с",
          f"Мяч падает 2 с, шаг 1/60 с. Две линии лежат одна на другой. Тест «mathbook: the slope…» ({RUN})",
          [("скорость v из движка", d["t_s"], d["velocity_mps"]), ("наклон высоты Δy/Δt", d["t_s"], d["slope_mps"])])

    d = read_csv(src, "step_order")
    chart("step-order", dst, "Вдвое меньше шаг — вдвое меньше ошибка", "шаг по времени Δt, с", "ошибка высоты через 2 с, м",
          f"Точная высота y₀ − g t²/2. Первый порядок: прямая с наклоном 1. Тест «mathbook: the slope…» ({RUN})",
          [("ошибка движка", d["dt_s"], d["height_error_m"])], xlog=True, ylog=True, dots=True)

    d = read_csv(src, "schemes")
    chart("schemes", dst, "Маятник тремя способами: энергия", "время, с", "энергия E / E₀",
          f"Маятник 1 м, 60°, шаг 0,05 с. Схема Эйлера разгоняет его. Тест «mathbook: a pendulum…» ({RUN})",
          [("Эйлер", d["t_s"], d["euler"]), ("симплектический Эйлер", d["t_s"], d["symplectic_euler"]),
           ("Рунге–Кутта 4", d["t_s"], d["runge_kutta_4"])])
    chart("schemes-zoom", dst, "Те же схемы крупно: без Эйлера", "время, с", "энергия E / E₀",
          f"Симплектическая схема дышит около 1, но не уходит; РК-4 держит 1 до 6·10⁻⁵. Тест «mathbook: a pendulum…» ({RUN})",
          [("симплектический Эйлер", d["t_s"], d["symplectic_euler"]), ("Рунге–Кутта 4", d["t_s"], d["runge_kutta_4"])])

    d = read_csv(src, "linear_systems")
    solved = next(k for k, r in enumerate(d["conjugate_gradients"]) if r < 1e-15)  # CG is exact from here on
    chart("linear-systems", dst, "1024 уравнения давления: Якоби и сопряжённые градиенты", "итерация",
          "невязка |f − A p| / |f|",
          f"Сетка 32 × 32, уравнение Пуассона. Тест «mathbook: a big linear system…» ({RUN})",
          [("Якоби", d["iteration"], d["jacobi"]),
           ("сопряжённые градиенты", d["iteration"][:solved], d["conjugate_gradients"][:solved])],
          ylog=True)

    d = read_csv(src, "flip")
    chart("flip", dst, "Эффект Джанибекова: вращение вокруг средней оси", "время, с", "ω вокруг средней оси, рад/с",
          f"Брусок 0,2 × 0,6 × 1 м, ω₀ = 4 рад/с, шаг 1/60 с. Точно — уравнения Эйлера. Тест «mathbook: the tennis-racket…» ({RUN})",
          [("движок", d["t_s"], d["engine_omega_y"]), ("точно (уравнения Эйлера)", d["t_s"], d["exact_omega_y"])])

    d = read_csv(src, "divergence")
    chart("divergence", dst, "Дивергенция скорости до и после проекции давлением", "x, м", "div u, 1/с",
          f"u = (sin πx · sin πy, 0, 0), средний ряд сетки 32. Тест «mathbook: div u…» ({RUN})",
          [("до проекции", d["x_m"], d["divergence_before"]), ("после проекции", d["x_m"], d["divergence_after"])],
          dots=True)

    d = read_csv(src, "sigma")
    chart("sigma", dst, "Разброс среднего падает как 1/√N", "сколько кубиков усредняли, N", "разброс среднего σ",
          f"Среднее N случайных чисел из [0, 1), 400 повторов. Теория σ = 1/√(12N). Тест «mathbook: the spread…» ({RUN})",
          [("измерено", d["n"], d["measured_spread"]), ("теория 1/√(12N)", d["n"], d["theory_spread"])],
          xlog=True, ylog=True, dots=True)
    print("pictures written to", dst)


if __name__ == "__main__":
    main()
