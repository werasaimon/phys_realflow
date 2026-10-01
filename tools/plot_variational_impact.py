#!/usr/bin/env python3
"""Plot measured variational impact results, without synthesizing data.

python3 tools/plot_variational_impact.py /tmp/phys-impact-research-20261001/bouncing \
    docs/book/evidence/variational-impact-20261001
"""
import argparse
import csv
import json
from pathlib import Path
from plot_action import line_chart
from plot_simbenchmark_bouncing import metric_plot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('raw', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    report = json.loads((args.raw / 'report.json').read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    engines = (('sdk', 'SDK SI'), ('sdk-ccd', 'SDK CCD'),
               ('sdk-variational', 'Вариационный опыт'))
    metric_plot(args.output / 'energy-step.svg', report['runs'], engines,
                'relative_energy_rmse', 'Упругие сферы: ошибка энергии', 'RMSE энергии, %', 100)
    metric_plot(args.output / 'position-step.svg', report['runs'], engines,
                'position_rmse_m', 'Упругие сферы: ошибка траектории', 'RMSE положения, м')
    series = []
    times = None
    for engine, label in engines:
        with (args.raw / (engine + '-0.01') / 'frames.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        t = [float(row['time_s']) for row in rows]
        if times is not None and t != times:
            raise ValueError('incompatible observation times')
        times = t
        series.append((label, [(float(row['energy_J']) / 24034.5 - 1) * 100 for row in rows]))
    line_chart(str(args.output / 'energy-time.svg'), 'Энергия: один опыт, разные шаговые схемы',
               'Модельное время, с', 'Отклонение энергии, %',
               '49 сфер, e=1, без трения; шаг 10 мс, CCD tolerance 2 мм. Режим сфер не применён к цепям.',
               times, series)


if __name__ == '__main__':
    main()
