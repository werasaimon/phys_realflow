import os, subprocess, time, json, hashlib
from pathlib import Path

root = Path('/home/wera_n/GIT/PHYSICS/phys_realflow')
out = Path('/tmp/phys-impact-research-20261001')
records = []

def run(name, command, extra=None, timeout=600):
    env = dict(os.environ, RF_THREADS='4')
    env.pop('RF_TEST', None)
    if extra: env.update(extra)
    full = ['nice', '-n', '15', 'taskset', '-c', '2-5', *map(str, command)]
    start = time.monotonic()
    print('START', name, flush=True)
    with (out/(name+'.log')).open('w') as log:
        try:
            code = subprocess.run(full, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired: code = 124
    records.append(dict(name=name, command=full, environment=extra or {}, exit_code=code, seconds=time.monotonic()-start))
    (out/'angular-trial-checks.json').write_text(json.dumps(records, indent=2)+'\n')
    print('END', name, code, flush=True)

for label, build in [('baseline', 'build-strict'), ('candidate', 'build-angular-trial')]:
    for kind, source, args in [('chain', 'angular_chain_probe.cpp', []), ('sandbox', 'sandbox_probe_control.cpp', ['substeps20-slop0.1mm'])]:
        binary = out/('angular-'+kind+'-'+label)
        cmd = ['c++', '-std=c++17', '-O3', '-ffp-contract=off', '-fno-fast-math', '-Wall', '-Wextra', '-Werror', '-Isrc', '-I.', str(out/source), str(out/build/'librfsamples.a'), str(out/build/'librfcore.a'), '-pthread', '-o', str(binary)]
        subprocess.run(cmd, cwd=root, check=True)
        run('angular-'+kind+'-'+label, [binary, *args])

patterns = ['Noether:', 'box stack', 'tall stack', 'mixed non-convex sandbox', 'six cups', 'pile of 50', 'short torus chain', 'continuous', 'CCD']
for i, pattern in enumerate(patterns):
    name = 'angular-regression-'+str(i).zfill(2)
    run(name, [out/'build-angular-trial/rf_tests'], {'RF_TEST':pattern, 'RF_JUNIT':str(out/(name+'.xml'))})
print('FINISHED ANGULAR TRIAL', flush=True)
