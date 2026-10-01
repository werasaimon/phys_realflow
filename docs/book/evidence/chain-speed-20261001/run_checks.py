import json, os, subprocess, time
from pathlib import Path

root=Path('/home/wera_n/GIT/PHYSICS/phys_realflow')
out=Path('/tmp/phys-chain-speed-20261001')
receipts=[]

def run(name, command, pattern=None, timeout=1200):
    env=dict(os.environ, RF_THREADS='5', RF_JUNIT=str(out/(name+'.xml')))
    env.pop('RF_TEST',None)
    if pattern: env['RF_TEST']=pattern
    cmd=['nice','-n','10','taskset','-c','1-5',*map(str,command)]
    print('START',name,flush=True)
    start=time.monotonic()
    with (out/(name+'.log')).open('w') as log:
        try: code=subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=timeout).returncode
        except subprocess.TimeoutExpired: code=124
    receipts.append(dict(name=name,command=cmd,filter=pattern,threads=5,exit_code=code,seconds=time.monotonic()-start))
    (out/'checks.json').write_text(json.dumps(receipts,indent=2)+'\n')
    print('END',name,code,flush=True)

run('short-chain',['build-core/rf_tests'],'short torus chain')
run('sandbox',['build-core/rf_tests'],'mixed non-convex sandbox')
run('rigid-determinism',['build-core/rf_tests'],'new samples give')
run('chain-determinism',['build-core/rf_tests'],'torus chains: determinism')
run('default',['build-core/rf_tests'])
run('readability',['python3','tools/readability.py','.'])
run('diff-check',['git','diff','--check'])
print('FINISHED',flush=True)
