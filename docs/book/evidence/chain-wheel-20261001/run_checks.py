import json, os, resource, subprocess, time
from pathlib import Path

root=Path('/home/wera_n/GIT/PHYSICS/phys_realflow')
out=Path('/tmp/phys-chain-wheel-20261001')
receipts=[]

def limits():
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))

def run(name,command,threads=5,pattern=None,seconds=180):
    env=dict(os.environ,RF_THREADS=str(threads),RF_JUNIT=str(out/(name+'.xml')))
    env.pop('RF_TEST',None)
    env.pop('RF_PLOT_DIR',None)
    if pattern: env['RF_TEST']=pattern
    cmd=['timeout','--kill-after=5s',str(seconds)+'s','nice','-n','10','taskset','-c','1-5',
         '/usr/bin/time','-v','-o',str(out/(name+'.resources')),*map(str,command)]
    print('START',name,flush=True)
    start=time.monotonic()
    with (out/(name+'.log')).open('w') as log:
        result=subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limits)
    receipts.append(dict(name=name,command=cmd,threads=threads,filter=pattern,exit_code=result.returncode,seconds=time.monotonic()-start,memory_limit_bytes=2*1024**3))
    (out/'checks.json').write_text(json.dumps(receipts,indent=2)+'\n')
    print('END',name,result.returncode,flush=True)
    return result.returncode

run('sat-audit',[out/'geometry-audit',out/'sat-audit',out/'native/chain-wheel-poses.csv','-2','41'])
for threads in [1,2,5]:
    run('threads-'+str(threads),[out/'study','40','120','2'],threads)
for substeps in [40,80,160]:
    run('refine-'+str(substeps),[out/'study',str(substeps),'240','2'])
for load in [1,3]:
    run('load-'+str(load),[out/'study','40','360',str(load)])
run('default',['build-core/rf_tests'],seconds=1200)
run('short-chain',['build-core/rf_tests'],pattern='short torus chain',seconds=180)
run('readability',['python3','tools/readability.py','.'])
run('diff-check',['git','diff','--check'])
print('FINISHED',flush=True)
