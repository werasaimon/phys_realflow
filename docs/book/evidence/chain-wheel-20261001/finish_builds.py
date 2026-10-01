import json, os, subprocess, time
from pathlib import Path
out=Path('/tmp/phys-chain-wheel-20261001')
root=Path('/home/wera_n/GIT/PHYSICS/phys_realflow')
strict=Path('/tmp/phys-impact-research-20261001/build-strict')
receipts=[]
print('Waiting for the bounded regression queue; no concurrent build load',flush=True)
deadline=time.monotonic()+1500
while time.monotonic()<deadline:
    try:
        previous=json.loads((out/'checks.json').read_text())
        if any(r['name']=='diff-check' for r in previous):
            if any(r['exit_code'] for r in previous): raise SystemExit('Regression failure: inspect checks.json before continuing')
            break
    except (FileNotFoundError,json.JSONDecodeError):pass
    time.sleep(2)
else: raise SystemExit('Regression queue did not finish in time')

def run(name,cmd,extra=None,seconds=240):
    env=dict(os.environ,RF_THREADS='5')
    env.pop('RF_TEST',None)
    if extra:env.update(extra)
    full=['timeout','--kill-after=5s',str(seconds)+'s','nice','-n','10','taskset','-c','1-5',*map(str,cmd)]
    print('START',name,flush=True)
    start=time.monotonic()
    with (out/(name+'.log')).open('w') as log:
        code=subprocess.run(full,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
    receipts.append(dict(name=name,command=full,exit_code=code,seconds=time.monotonic()-start))
    (out/'build-checks.json').write_text(json.dumps(receipts,indent=2)+'\n')
    print('END',name,code,flush=True)
    if code:raise SystemExit(code)

run('strict-samples-build',['cmake','--build',strict,'--target','rfsamples','--parallel','2'])
run('strict-check-build',['c++','-std=c++17','-O3','-ffp-contract=off','-DRF_STRICT_FP=1','-Wall','-Wextra','-Werror','-Isrc','-I.',out/'strict_main.cpp','tests/ChainWheelTests.cpp',strict/'librfsamples.a',strict/'librfcore.a','-pthread','-o',out/'strict-checks'])
run('strict-chain-wheel',[out/'strict-checks'],{'RF_JUNIT':str(out/'strict-chain-wheel.xml'),'RF_PLOT_DIR':str(out/'strict')})
run('editor-build',['cmake','--build','/tmp/phys-realflow-editor-standard-build','--parallel','2'],seconds=360)
print('BUILDS FINISHED',flush=True)
