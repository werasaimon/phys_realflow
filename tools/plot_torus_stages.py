#!/usr/bin/env python3
"""Plot measured contacts and one-step motion checks of the isolated torus stage observer.

Usage: python3 tools/plot_torus_stages.py TRACE_DIR OUTPUT_DIR
TRACE_DIR contains stages.csv, contacts.csv and motion.csv. Each motion mode starts
from the same current state; these checks are not independent simulation trajectories.
Python standard library only. Source data and observer hooks are retained with the report.
"""
import argparse
import csv
import math
from pathlib import Path

from plot_action import line_chart


def read(path):
    with path.open(encoding='utf-8', newline='') as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f'no measurements in {path}')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trace', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    stages = {(int(r['substep']), int(r['phase'])): r for r in read(args.trace / 'stages.csv')}
    motion = {(int(r['substep']), int(r['mode'])): r for r in read(args.trace / 'motion.csv')}
    contacts = {}
    for r in read(args.trace / 'contacts.csv'):
        contacts.setdefault((int(r['substep']), int(r['phase'])), []).append(r)
    steps = sorted(s for s, p in stages if p == 2 and 470 <= s <= 495)
    if len(steps) < 2 or any((s, mode) not in motion for s in steps for mode in (0, 3)):
        parser.error('trace must include matching stages and motion modes 0/3 for substeps 470..495')
    args.output.mkdir(parents=True, exist_ok=True)
    raw = [float(stages[s, 2]['raw_depth_m']) * 1000 for s in steps]
    kept = [max((float(r['depth_m']) * 1000 for r in contacts.get((s, 2), [])), default=0) for s in steps]
    line_chart(str(args.output / 'torus-contact-depth.svg'),
        'Перед обновлением поз: какие проникновения видит решатель', 'Подшаг', 'Глубина, мм',
        'Пара 67/68. Геометрия и ограничения наблюдаются на одном этапе; ноль означает отсутствие точек.',
        steps, [('Все части', raw), ('В решателе', kept)])
    physical, full = [], []
    for s in steps:
        start = float(stages[s, 2]['distance_m'])
        physical.append(1000 * (float(motion[s, 0]['distance_m']) - start))
        full.append(1000 * (float(motion[s, 3]['distance_m']) - start))
    if not all(math.isfinite(v) for v in raw + kept + physical + full):
        parser.error('nonfinite measurements')
    line_chart(str(args.output / 'torus-split-motion.svg'),
        'Один подшаг: вклад коррекции проникновения в расхождение', 'Подшаг', 'Прирост расстояния, мм',
        'Каждый расчёт начинает со штатного состояния. На шаге 490 зацепление исчезает только с поступательной коррекцией.',
        steps, [('Физические скорости', physical), ('С коррекцией', full)])
    print('plots written to', args.output)


if __name__ == '__main__':
    main()
