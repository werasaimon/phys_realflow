#!/usr/bin/env python3
"""Run the standard rigid suite and preserve its log, JUnit, frame CSVs and Markdown report.

Usage: python3 tools/run_rigid_standard.py --build-dir build-core --threads 4
This runs PhysRealFlow only; comparisons with other engines require an external shared harness.
"""
import argparse
import csv
import datetime
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError('must be positive')
    return number


def percentile(values, fraction):
    return sorted(values)[max(0, math.ceil(len(values) * fraction) - 1)]


def frame_summary(path):
    with path.open(encoding='utf-8', newline='') as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f'empty frame report: {path}')
    times = [float(row['step_ms']) for row in rows]
    active = [float(row['step_ms']) for row in rows if int(row['asleep']) < int(row['dynamic'])]
    active_median = f'{percentile(active, 0.5):.3f}' if active else '—'
    last = rows[-1]
    return (f"| {path.stem} | {len(rows)} | {percentile(times, 0.5):.3f} | "
            f"{active_median} | {percentile(times, 0.95):.3f} | {last['asleep']}/{last['dynamic']} | "
            f"{1000 * float(last['max_floor_penetration_m']):.3f} | "
            f"{1000 * float(last['late_max_path_m']):.3f} |")


def torus_summary(path):
    with path.open(encoding='utf-8', newline='') as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f'empty torus report: {path}')
    times = [float(row['step_ms']) for row in rows]
    broken = max(int(row['broken_pairs']) for row in rows)
    first = next((row['frame'] for row in rows if int(row['broken_pairs'])), '—')
    last = rows[-1]
    return (f"| {path.stem} | {len(rows)} | {percentile(times, 0.5):.3f} | "
            f"{percentile(times, 0.95):.3f} | {broken} | {first} | "
            f"{1000 * float(last['max_distance_m']):.3f} | "
            f"{float(last['energy_gain_j']):.6f} |")


def write_report(output, metadata, returncode):
    suite = ET.parse(output / 'junit.xml').getroot()
    cases = suite.findall('testcase')
    failures = [case for case in cases if case.find('failure') is not None]
    report = [
        '# Стандартные тесты твёрдых тел', '',
        f"UTC: {metadata['utc']}; commit: `{metadata['commit']}`; потоков: {metadata['threads']}.", '',
        f"Условия: {(metadata.get('note') or 'дополнительная нагрузка не записана').rstrip('.')}.", '',
        f"Результат: {len(cases) - len(failures)}/{len(cases)} тестов прошли; код выхода {returncode}.", '',
        'Времена — измерения этого запуска PhysRealFlow, без сравнения с другими движками. '
        'CSV содержит отсчёты раз в кадр; узкие пики между подшагами могут быть пропущены.', '',
    ]
    paths = sorted((output / 'frames').glob('*.csv'))
    regular_paths = [path for path in paths if not path.stem.startswith('torus-chains-')]
    if regular_paths:
        report += ['| Сцена | Кадры | Медиана, мс | Медиана неспящих кадров, мс | p95, мс | Спят | Макс. проникновение пола, мм | Путь за последние 2 с, мм |',
                   '|---|---:|---:|---:|---:|---:|---:|---:|']
        report += [frame_summary(path) for path in regular_paths]
    torus_paths = [path for path in paths if path.stem.startswith('torus-chains-')]
    if torus_paths:
        report += ['', '| Цепи торов | Кадры | Медиана, мс | p95, мс | Макс. разрывов | Первый разрыв, кадр | Макс. расстояние, мм | Прирост энергии, Дж |',
                   '|---|---:|---:|---:|---:|---:|---:|---:|']
        report += [torus_summary(path) for path in torus_paths]
    report += ['', '| Проверка | Результат | Время, с |', '|---|---|---:|']
    for case in cases:
        state = 'FAIL' if case.find('failure') is not None else 'PASS'
        name = case.attrib['name'].replace('|', '\\|')
        report.append(f"| {name} | {state} | {case.attrib['time']} |")
    report += ['', 'Ошибки и остальные физические измерения: [run.log](run.log). '
               'Состояние сборки и машины: [metadata.json](metadata.json).', '']
    (output / 'report.md').write_text('\n'.join(report), encoding='utf-8')
    return len(cases), len(failures)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', type=Path, default=Path('build-core'))
    parser.add_argument('--threads', type=positive_int, default=4)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--note', default='', help='record timing conditions, e.g. concurrent builds or other tests')
    parser.add_argument('--filter', default='standard rigid:', help='RF_TEST substring; e.g. "torus chains:"')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    build = args.build_dir.resolve()
    binary = next((path for path in (build / 'rf_tests', build / 'rf_tests.exe', build / 'Release/rf_tests.exe')
                   if path.is_file()), None)
    if binary is None:
        parser.error(f'build rf_tests first: no executable in {build}')
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    output = (args.output or build / f'rigid-standard-{stamp}').resolve()
    if output.exists():
        parser.error(f'output already exists; choose a new run directory: {output}')
    output.mkdir(parents=True)
    (output / 'frames').mkdir()
    metadata = {
        'utc': stamp, 'platform': platform.platform(), 'processor': platform.processor(),
        'logical_cpus': os.cpu_count(), 'threads': args.threads, 'binary': str(binary),
        'note': args.note,
        'filter': args.filter,
        'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
        'worktree': subprocess.check_output(['git', 'status', '--short'], cwd=root, text=True),
    }
    cache = build / 'CMakeCache.txt'
    metadata['cmake'] = [line for line in cache.read_text().splitlines()
                         if line.startswith(('CMAKE_BUILD_TYPE:', 'CMAKE_CXX_COMPILER:', 'RF_STRICT_FP:', 'RF_WERROR:', 'RF_SANITIZE:'))]
    (output / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    env = dict(os.environ, RF_TEST=args.filter, RF_THREADS=str(args.threads),
               RF_JUNIT=str(output / 'junit.xml'), RF_PLOT_DIR=str(output / 'frames'))
    with (output / 'run.log').open('w', encoding='utf-8') as log:
        process = subprocess.Popen([str(binary)], cwd=root, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            log.write(line)
            log.flush()
            print(line, end='', flush=True)
        returncode = process.wait()
    try:
        count, failed = write_report(output, metadata, returncode)
    except (OSError, ValueError, ET.ParseError) as error:
        print(f'Report incomplete: {error}; log: {output / "run.log"}', file=sys.stderr)
        return returncode or 1
    print(f'Report: {output / "report.md"}')
    return returncode or (0 if count and not failed else 1)


if __name__ == '__main__':
    sys.exit(main())
