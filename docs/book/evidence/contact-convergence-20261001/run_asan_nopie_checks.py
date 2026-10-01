import os, subprocess, time, json
from pathlib import Path

root=Path('/home/wera_n/GIT/PHYSICS/phys_realflow')
out=Path('/tmp/phys-impact-research-20261001')
env=dict(os.environ, RF_THREADS='10')
receipts=[]

def run(name, command, extra=None, timeout=1800):
    start=time.monotonic()
    print('START',name,flush=True)
    values=dict(env)
    values.pop('RF_TEST',None)
    if extra: values.update(extra)
    full=['nice','-n','10','taskset','-c','1-5,7-11',*command]
    with (out/(name+'.log')).open('w') as log:
        try:
            code=subprocess.run(full,cwd=root,env=values,stdout=log,stderr=subprocess.STDOUT,timeout=timeout).returncode
        except subprocess.TimeoutExpired: code=124
    receipts.append(dict(name=name,command=full,exit_code=code,seconds=time.monotonic()-start))
    (out/'asan-nopie-checks.json').write_text(json.dumps(receipts,indent=2)+'\n')
    print('END',name,code,flush=True)
    return code == 0

asan=out/'build-sanitize'
patterns=['variational:', 'CCD', 'math', 'primitives', 'bvh', 'mass', 'gjk', 'rigid bodies', 'box stack',
          'tall stack', 'joints', 'continuous', 'Noether:', 'Newton', 'short torus chain survives', 'shock pass']
for i,pattern in enumerate(patterns):
    name='asan-nopie-'+str(i).zfill(2)
    run(name,[str(asan/'rf_tests')],{'RF_TEST':pattern,'RF_JUNIT':str(out/(name+'.xml')),
        'ASAN_OPTIONS':'detect_leaks=0','UBSAN_OPTIONS':'print_stacktrace=1:halt_on_error=1'},600)
print('FINISHED SANITIZERS',flush=True)
