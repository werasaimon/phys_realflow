#!/usr/bin/env python3
"""Plot measured chain motor work and energy residual, using the standard library.

python3 tools/plot_motor_work.py INPUT/chain-energy.csv OUTPUT_DIRECTORY
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
    initial = rows[0]["energy_J"] - rows[0]["motor_work_J"] - rows[0]["residual_J"]
    line_chart(str(args.output / "chain-work.svg"), "Цепь: энергия и работа двигателя", "Время, с", "Энергия, Дж",
               "Итоговый моторный импульс каждого принятого подшага; промежуточные итерации не считаются временем.", times,
               [("E − E₀", [row["energy_J"] - initial for row in rows]),
                ("Работа мотора", [row["motor_work_J"] for row in rows])])
    line_chart(str(args.output / "chain-residual.svg"), "Цепь: остаток энергетического баланса", "Время, с", "Остаток, Дж",
               "E − E₀ − W; отрицательный остаток включает трение, демпфирование и численные потери.", times,
               [("Остаток", [row["residual_J"] for row in rows])])


if __name__ == "__main__":
    main()
