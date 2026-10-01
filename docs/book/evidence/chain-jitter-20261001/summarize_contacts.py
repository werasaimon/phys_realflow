"""Summarize a single recorded frame; no trajectory or solver changes.

Usage: python3 summarize_contacts.py PATH_PREFIX
Reads PREFIX-contacts.csv, PREFIX-energy.csv and PREFIX.json.
"""
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
import sys


def summarize(prefix):
    with Path(str(prefix) + '-contacts.csv').open() as file:
        rows = [row for row in csv.DictReader(file) if row['stage'] == 'warm']
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[key] for key in ('substep', 'a', 'b', 'patch'))].append(row)
    reused = duplicates = loaded_aliases = 0
    for group in groups.values():
        counts = Counter(row['chosen_old'] for row in group if int(row['chosen_old']) >= 0)
        reused += any(count > 1 for count in counts.values())
        duplicates += sum(count - 1 for count in counts.values())
        loaded_aliases += sum(float(row['normal_impulse_Ns']) > 0 and
                              counts[row['chosen_old']] > 1 for row in group)
    stages = defaultdict(lambda: dict(delta_T_J=0.0, delta_U_J=0.0, delta_E_J=0.0))
    with Path(str(prefix) + '-energy.csv').open() as file:
        for row in csv.DictReader(file):
            for key in stages[row['stage']]:
                stages[row['stage']][key] += float(row[key])
    actual = json.loads(Path(str(prefix) + '.json').read_text())
    return dict(frame=actual['frame'], time_s=actual['frame']/60, actual=actual,
                warm_rows=len(rows), patch_observations=len(groups),
                patches_with_reused_old=reused, excess_assignments=duplicates,
                loaded_rows_in_reused_groups=loaded_aliases,
                hash_multi_hits=sum(int(row['hash_hits']) > 1 for row in rows),
                chosen_not_nearest=sum(row['chosen_old'] != row['nearest_old'] for row in rows),
                max_chosen_distance_m=max(float(row['chosen_distance_m']) for row in rows),
                max_normal_angle_deg=max(float(row['normal_angle_rad']) for row in rows)*180/math.pi,
                stage_totals=dict(stages),
                scope='two chosen frames; reused old points observed, causal effect not isolated')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    print(json.dumps(summarize(Path(sys.argv[1])), indent=2))
