#!/usr/bin/env python3
"""Plot measured nut motion, pitch error and energy, using the standard library.

python3 tools/plot_nut_bolt.py INPUT/nut-bolt.csv OUTPUT_DIRECTORY
"""
import argparse
import csv
import math
from pathlib import Path

from plot_action import line_chart


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with args.csv.open(newline="", encoding="utf-8") as stream:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]
    if not rows or any(not math.isfinite(value) for row in rows for value in row.values()):
        raise ValueError("empty or non-finite observations")
    times = [row["time_s"] for row in rows]
    if any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError("simulation time must increase")
    args.output.mkdir(parents=True, exist_ok=True)
    line_chart(str(args.output / "nut-lift.svg"), "Гайка: подъём и обратный ход", "Время, с", "Высота, мм",
               "Свободная гайка; привод прикладывает только ограниченный угловой импульс.", times,
               [("Начало координат гайки", [1000 * row["height_m"] for row in rows])])
    line_chart(str(args.output / "nut-pitch-error.svg"), "Отклонение от идеального шага", "Время, с", "Ошибка, мм",
               "Начало отсчёта — после 3 с оседания; зазор и фасетки остаются частью геометрии.", times,
               [("y − y₀ − p(θ − θ₀)/(2π)", [1000 * row["screw_error_m"] for row in rows])])
    line_chart(str(args.output / "nut-residual.svg"), "Энергия после вычитания работы привода", "Время, с", "Остаток, Дж",
               "Отрицательный остаток включает физические и численные потери; работа мыши не измерена.", times,
               [("E − E₀ − W", [row["residual_J"] for row in rows])])


if __name__ == "__main__":
    main()
