"""Bounded replay of two late frames of the isolated RP3D-style prototype.

Run from the repo root after building build-core libraries. Linux/GCC only.
Writes temporary artifacts; no SDK source changes. Requires CPUs 1-5.
"""
import json
import os
from pathlib import Path
import resource
import subprocess
import tempfile
import time

from summarize_contacts import summarize

evidence = Path(__file__).resolve().parent
repo = evidence.parents[3]
output = Path(tempfile.mkdtemp(prefix='phys-chain-contacts-'))
env = dict(os.environ, RF_THREADS='5')
for key in ('RF_TEST', 'RF_CHAIN_UNLOCK', 'RF_CHAIN_NO_SHOCK', 'RF_PLOT_DIR'):
    env.pop(key, None)
source = (repo/'src/rigid/ContactSolver.cpp').read_text()
old = '        if (p.positionBias <= 0) continue;'
new = '        if (p.depth < 0) continue; // Keep zero-target rows for touching contacts (RP3D-style).'
if source.count(old) != 1:
    raise SystemExit('This replay requires the recorded baseline split guard')
(output/'ContactSolver-rp3.cpp').write_text(source.replace(old, new))
subprocess.run(['nice', '-n', '10', 'taskset', '-c', '1-5', 'c++', '-std=gnu++17',
                '-O3', '-DNDEBUG', '-march=native', '-Isrc', '-c',
                str(output/'ContactSolver-rp3.cpp'), '-o', str(output/'ContactSolver-rp3.o')],
               cwd=repo, check=True)
subprocess.run(['nice', '-n', '10', 'taskset', '-c', '1-5', 'c++', '-std=c++17',
                '-O3', '-fno-access-control', '-Isrc', '-I.', str(evidence/'late_contact_probe.cpp'),
                str(output/'ContactSolver-rp3.o'), 'build-core/librfsamples.a',
                'build-core/librfcore.a', '-pthread', '-o', str(output/'late-contact-rp3')],
               cwd=repo, check=True)


def limit():
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))


checks = []
for frame, pull in ((1399, 0), (1728, 1)):
    prefix = output/('late-contacts-' + str(frame))
    command = ['timeout', '--kill-after=5s', '180s', 'nice', '-n', '10',
               'taskset', '-c', '1-5', '/usr/bin/time', '-v', '-o', str(prefix)+'.resources',
               str(output/'late-contact-rp3'), str(frame), str(pull),
               str(prefix)+'-energy.csv', str(prefix)+'-contacts.csv']
    started = time.monotonic()
    with Path(str(prefix)+'.json').open('w') as stdout:
        with Path(str(prefix)+'.stderr.txt').open('w') as stderr:
            result = subprocess.run(command, cwd=repo, env=env, stdout=stdout,
                                    stderr=stderr, preexec_fn=limit)
    checks.append(dict(frame=frame, command=command, exit_code=result.returncode,
                       seconds=time.monotonic()-started))
    (output/'checks.json').write_text(json.dumps(checks, indent=2)+'\n')
    if result.returncode:
        raise SystemExit(result.returncode)
    summary = summarize(prefix)
    Path(str(prefix)+'-summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(frame, summary['actual'], flush=True)
print('Artifacts:', output)
