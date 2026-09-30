#!/usr/bin/env python3
"""Summarize measured chain settling and plot energy/height (Python standard library only).

Usage: python3 tools/plot_chain_relaxation.py /tmp/phys-chain-relaxation-20260930 OUTPUT_DIR
Each completed subdirectory must contain ChainRelaxation's energy.csv and final.csv.
Time zero is release for shaken, initial placement for straight, and restart for branches.
The free-damped endpoint is a comparison state, never a certified global minimum.
"""
import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

from plot_action import line_chart


def read(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(stream)]
    if not rows or any(not math.isfinite(v) for r in rows for v in r.values()):
        raise ValueError(f"invalid observations: {path}")
    return rows


def summarize(rows):
    final = rows[-1]
    late = [r for r in rows if r["time_s"] >= final["time_s"] - 10]
    nonnegative = [r for r in rows if r["time_s"] >= 0]
    energy = [r["kinetic_j"] + r["potential_j"] for r in nonnegative]
    speed = [r["max_speed_mps"] for r in late]
    tmean = statistics.mean(r["time_s"] for r in late)
    ymean = statistics.mean(r["com_y"] for r in late)
    slope = sum((r["time_s"] - tmean) * (r["com_y"] - ymean) for r in late)
    slope /= sum((r["time_s"] - tmean) ** 2 for r in late)
    return {
        "frames": len(rows) - 1, "final": final,
        "late10_mean_kinetic_j": statistics.mean(r["kinetic_j"] for r in late),
        "late10_max_kinetic_j": max(r["kinetic_j"] for r in late),
        "late10_max_speed_mps": max(speed),
        "late10_mean_max_speed_mps": statistics.mean(speed),
        "late10_com_height_range_mm": 1000 * (max(r["com_y"] for r in late) - min(r["com_y"] for r in late)),
        "late10_com_height_slope_mm_per_s": 1000 * slope,
        "tip_radial_mm": 1000 * math.hypot(final["tip_x"], final["tip_z"]),
        "max_contact_depth_mm": 1000 * max(r["contact_depth_m"] for r in rows),
        "max_energy_rise_per_frame_j": max(b - a for a, b in zip(energy, energy[1:])),
        "sampled_ccd_hits": sum(r["ccd_hits_last_substep"] for r in rows),
        "max_broken": max(r["broken"] for r in rows),
        "max_sleepers": max(r["sleepers"] for r in rows),
    }


def plot_group(output, data, names, title, filename, key, ylabel, transform=lambda x: x):
    if not all(name in data for name, _ in names):
        return
    curves = [[r for r in data[name] if r["time_s"] >= 0] for name, _ in names]
    # Compare common frame times, without interpolating different time grids.
    times = sorted(set.intersection(*(set(r["time_s"] for r in c) for c in curves)))
    series = []
    for (_, label), curve in zip(names, curves):
        values = {r["time_s"]: transform(r[key]) for r in curve}
        series.append((label, [values[t] for t in times]))
    line_chart(str(output / filename), title, "Время после запуска / отпускания, с", ylabel,
               "8 торов; сон выключен. Данные ChainRelaxation.cpp. Сравнение состояний не доказывает глобальный минимум.",
               times, series)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data, report = {}, {}
    for path in sorted(args.runs.glob("*/energy.csv")):
        final = path.with_name("final.csv")
        if not final.exists() or len(final.read_text().splitlines()) != 9:
            continue
        name = path.parent.name
        measured = read(path)
        if len(measured) < 2:
            continue  # A zero-frame snapshot inspection is geometry, not a settling trajectory.
        data[name] = measured
        report[name] = summarize(data[name])
        report[name]["csv_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        report[name]["poses"] = read(final)
        sat = path.parent / "sat" / "body-checkpoints.csv"
        if sat.exists():
            samples = read(sat)
            report[name]["sat_frames_checked"] = len(samples)
            report[name]["sampled_max_child_sat_depth_mm"] = 1000 * max(r["max_child_sat_depth_m"] for r in samples)
            report[name]["final_child_sat_depth_mm"] = 1000 * samples[-1]["max_child_sat_depth_m"]
    if "free-damped" in report:
        reference = report["free-damped"]["final"]["potential_j"]
        for r in report.values():
            r["energy_above_free_damped_j"] = r["final"]["potential_j"] - reference
            r["height_above_free_damped_mm"] = 1000 * r["energy_above_free_damped_j"] / (9.81 * r["final"]["mass_kg"])
    (args.output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    plot_group(args.output, data, [("shaken", "После рывка"), ("straight", "Без рывка")],
               "Остаточная кинетическая энергия", "kinetic.svg", "kinetic_j", "log10(K / Дж)",
               lambda k: math.log10(max(k, 1e-12)))
    groups = [([("control", "Контроль"), ("no-roll", "Без качения"), ("no-coulomb", "Без μ")], "friction"),
              ([("control", "Контроль"), ("free-damped", "Без трения"), ("free-fine", "Без трения, шаг / 2")], "reference")]
    base = report.get("shaken", {}).get("final", {}).get("com_y", 0)
    for names, name in groups:
        plot_group(args.output, data, names, "Изменение высоты центра масс после перезапуска",
                   name + ".svg", "com_y", "Δh центра масс, мм", lambda y: 1000 * (y - base))
    for name, r in report.items():
        f = r["final"]
        print(f"{name:15} U={f['potential_j']:.6f} J, K={f['kinetic_j']:.3g} J, "
              f"COM={1000*f['com_y']:.3f} mm, tip={r['tip_radial_mm']:.3f} mm, "
              f"late dh/dt={r['late10_com_height_slope_mm_per_s']:.5f} mm/s")


if __name__ == "__main__":
    main()
