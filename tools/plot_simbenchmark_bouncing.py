#!/usr/bin/env python3
"""Plot the adapted benchmark using recorded metrics and raw traces; Python stdlib only.

python3 tools/plot_simbenchmark_bouncing.py /tmp/simbenchmark-final-bouncing-v2-20261001 \
    docs/book/evidence/simbenchmark-bouncing-20261001
These are specific contact configurations, not an overall engine ranking.
"""
import argparse
import csv
import json
from pathlib import Path

from plot_action import (Axis, BOTTOM, H, LEFT, RIGHT, TOP, W, frame, grid,
                         legend, line_chart, log_range, polyline)


def metric_plot(path, runs, engines, metric, title, ylabel, factor=1):
    groups = []
    for engine, label in engines:
        rows = sorted((r for r in runs if r["engine"] == engine), key=lambda r: r["requested_dt_s"])
        if len(rows) != 5 or any("failure" in row for row in rows):
            raise ValueError("incomplete five-step sweep: " + engine)
        groups.append((label, [r["requested_dt_s"] for r in rows],
                       [r["summary"][metric] * factor for r in rows]))
    xs = [x for _, x, _ in groups][0]
    if any(x != xs for _, x, _ in groups):
        raise ValueError("timestep mismatch")
    values = [v for _, _, y in groups for v in y]
    xlo, xhi, xticks = log_range(min(xs), max(xs))
    ylo, yhi, yticks = log_range(min(values), max(values))
    xa = Axis(xlo, xhi, LEFT, W - RIGHT, log=True)
    ya = Axis(ylo, yhi, H - BOTTOM, TOP, log=True)
    out = frame(title, "Шаг времени, с (логарифмическая шкала)", ylabel,
                "49 шаров, e=1, 20 с. Из report.json. Различия контактных профилей описаны в главе 13.")
    grid(out, xa, ya, xticks, yticks)
    for slot, (_, x, y) in enumerate(groups):
        polyline(out, xa, ya, x, y, slot)
        for px, py in zip(x, y):
            out.append(f"<circle cx='{xa(px):.1f}' cy='{ya(py):.1f}' r='3' fill='var(--s{slot})'/>")
    legend(out, [name for name, _, _ in groups])
    out.append("</svg>")
    path.write_text("\n".join(out) + "\n")


def energy_trace(raw, output):
    series = []
    times = None
    for engine, label in (("sdk", "SDK"), ("sdk-ccd", "SDK CCD"), ("jolt", "Jolt")):
        with (raw / (engine + "-0.01") / "frames.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        t = [float(r["time_s"]) for r in rows]
        if times is not None and t != times:
            raise ValueError("trace sampling mismatch")
        times = t
        series.append((label, [(float(r["energy_J"]) / 24034.5 - 1) * 100 for r in rows]))
    line_chart(str(output / "bouncing-energy-time.svg"), "Упругий отскок: энергия на шаге 10 мс",
               "Модельное время, с", "Отклонение энергии, %",
               "Из frames.csv; SDK и Jolt почти совпадают. CCD-профиль показывает потерю энергии при дроблении.",
               times, series)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    report = json.loads((args.raw / "report.json").read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    common = (("sdk", "SDK"), ("bullet", "Bullet"), ("jolt", "Jolt"))
    metric_plot(args.output / "bouncing-energy-step.svg", report["runs"], common,
                "relative_energy_rmse", "SimBenchmark Bouncing: ошибка энергии", "Относительная RMSE, %", 100)
    metric_plot(args.output / "bouncing-position-step.svg", report["runs"], common,
                "position_rmse_m", "SimBenchmark Bouncing: ошибка траектории", "RMSE положения, м")
    alternatives = (("sdk-ccd", "SDK CCD"), ("box3d", "Box3D 0 Hz"), ("box3d-native", "Box3D native"))
    metric_plot(args.output / "bouncing-contact-profiles.svg", report["runs"], alternatives,
                "relative_energy_rmse", "Bouncing: дополнительные контактные профили", "Относительная RMSE, %", 100)
    energy_trace(args.raw, args.output)
    print("plots written to", args.output)


if __name__ == "__main__":
    main()
