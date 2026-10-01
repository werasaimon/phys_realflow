#!/usr/bin/env python3
"""Plot the measured drive motion; standard library only.

python3 tools/plot_chain_drive.py INPUT/chain-drive.csv OUTPUT_DIRECTORY
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
    if not rows or any(not math.isfinite(v) for r in rows for v in r.values()):
        raise ValueError("empty or non-finite observations")
    times = [r["time_s"] for r in rows]
    if any(b <= a for a, b in zip(times, times[1:])):
        raise ValueError("simulation times must increase")
    args.output.mkdir(parents=True, exist_ok=True)
    line_chart(str(args.output / "angle.svg"), "Вращение колеса: два цикла", "Время симуляции, с", "Угол, градусы",
               "Измеренный угол шарнира. Привод ограничен 20 Н·м; удержание и обратный ход каждые 10 с.",
               times, [("Угол колеса", [math.degrees(r["angle_rad"]) for r in rows])])
    line_chart(str(args.output / "lift.svg"), "Подъём и опускание груза", "Время симуляции, с", "Высота центра кольца, мм",
               "Цепь движется благодаря контактам; позы колец не задаются контроллером.",
               times, [("Жёлтое звено", [1000 * r["tip_y_m"] for r in rows])])


if __name__ == "__main__":
    main()
