#!/usr/bin/env python3
"""Plot the measured stationary chain-wheel test; standard library only.

python3 tools/plot_chain_wheel.py INPUT/chain-wheel.csv OUTPUT_DIRECTORY
Input comes from RF_TEST="chain wheel" RF_PLOT_DIR=INPUT build-core/rf_tests.
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
        rows = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
    if not rows or any(not math.isfinite(v) for row in rows for v in row.values()):
        raise ValueError("empty or non-finite measurements")
    times = [row["time_s"] for row in rows]
    if any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError("simulation times must increase")
    args.output.mkdir(parents=True, exist_ok=True)
    caption = "16 колец, сон выключен. Мышь действует от 3 до 5 с и вносит работу; это не замкнутая система."
    line_chart(str(args.output / "energy.svg"), "Цепь на неподвижной звёздочке", "Время симуляции, с", "E, Дж", caption,
               times, [("Полная энергия", [r["energy_J"] for r in rows])])
    line_chart(str(args.output / "kinetic.svg"), "Остаточное движение после отпускания", "Время симуляции, с", "log10(T / Дж)",
               "В 5 с мышь отпущена. Малое движение не доказывает глобальный минимум энергии при трении.",
               times, [("Кинетическая энергия", [math.log10(max(r["kinetic_J"], 1e-12)) for r in rows])])


if __name__ == "__main__":
    main()
