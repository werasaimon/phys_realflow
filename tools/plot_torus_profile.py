#!/usr/bin/env python3
"""Plot recorded torus profiling CSVs using the repository's standard-library SVG helpers.

Usage: python3 tools/plot_torus_profile.py FRAME_DIR PAIR_DIR OUTPUT_DIR
Raw pair contacts are observed after integration; frame timings cover the actual Simulation.
"""
import argparse
import csv
from pathlib import Path

from plot_action import line_chart


def read(path):
    with path.open(newline='', encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f'no measurements in {path}')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('frames', type=Path)
    parser.add_argument('pair', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    frames = read(args.frames / 'frames.csv')
    recorded_pair = read(args.pair / 'pair.csv')
    pair = [r for r in recorded_pair if 450 <= int(r['substep']) <= 510]
    if len(pair) < 2:
        parser.error('pair.csv must contain at least two measurements in substeps 450..510')
    args.output.mkdir(parents=True, exist_ok=True)
    first_break = next((r for r in recorded_pair if r['linking'] != recorded_pair[0]['linking']), None)
    break_note = (f"Изменение зацепления: шаг {first_break['substep']}, t={float(first_break['time_s']):.4f} с."
                  if first_break else 'Изменение зацепления в записанных отсчётах не найдено.')
    x = [float(r['time_s']) for r in pair]
    series = [(label, [float(r[key]) for r in pair]) for label, key in (
        ('Исходные', 'raw_pairs'), ('Отброшены', 'rejected_pairs'), ('Оставлены', 'retained_pairs'))]
    line_chart(str(args.output / 'torus-seam-pairs.svg'),
        'Наблюдаемая пара: отбор контактных участков', 'Время симуляции, с', 'Пар выпуклых частей',
        break_note + ' Контакты пересчитаны после интеграции.', x, series)
    series = [(label, [float(r[key]) for r in pair]) for label, key in (
        ('Прошли проверку', 'safe_points'), ('Теряются с участком', 'lost_safe_points'))]
    line_chart(str(args.output / 'torus-lost-points.svg'),
        'Наблюдаемая пара: первая точка определяет судьбу остальных', 'Время симуляции, с', 'Точек',
        'Индивидуальный предикат внутренних швов. Это проверка отбора, а не доказательство правильности нормалей.', x, series)
    x = [float(r['time_s']) for r in frames]
    series = [(label, [float(r[key]) for r in frames]) for label, key in (
        ('Поиск контактов', 'collide_ms'), ('Решатель', 'solve_ms'), ('CCD', 'ccd_ms'))]
    line_chart(str(args.output / 'torus-stage-times.svg'),
        'Сцена 39: время этапов на кадр', 'Время симуляции, с', 'Время, мс',
        'Из frames.csv профилировщика. Это времена записанного запуска, не рейтинг движков.', x, series)
    print('plots written to', args.output)


if __name__ == '__main__':
    main()
