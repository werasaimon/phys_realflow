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
    (out/'sandbox-final-checks.json').write_text(json.dumps(receipts,indent=2)+'\n')
    print('END',name,code,flush=True)
    return code == 0

for build,label in [('build-core','native'),(str(out/'build-strict'),'strict')]:
    for name,pattern in [('sandbox','mixed non-convex sandbox'),('determinism','new samples give'),('presets','presets'),('coherence','coherence:')]:
        run(label+'-final-'+name,[build+'/rf_tests'],{'RF_TEST':pattern,'RF_JUNIT':str(out/(label+'-final-'+name+'.xml'))},600)
run('final-readability',['python3','tools/readability.py','.'])
run('final-diff-check',['git','diff','--check'])
