import json, os, subprocess, resource, time
from pathlib import Path
out=Path('/tmp/phys-chain-drive-20261001')
strict=Path('/tmp/phys-impact-research-20261001/build-strict')
rows=[]
def limit():resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
def run(name,cmd,pattern=None,extra=None,timeout=180,limited=True):
    env=dict(os.environ,RF_THREADS='5',RF_JUNIT=str(out/(name+'.xml')))
    env.pop('RF_TEST',None);env.pop('RF_PLOT_DIR',None)
    if pattern:env['RF_TEST']=pattern
    if extra:env.update(extra)
    command=['timeout','--kill-after=5s',str(timeout)+'s','nice','-n','10','taskset','-c','1-5','/usr/bin/time','-v','-o',str(out/(name+'.resources')),*map(str,cmd)]
    print('START',name,flush=True);start=time.monotonic()
    with (out/(name+'.log')).open('w') as log:
        code=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limit if limited else None).returncode
    rows.append(dict(name=name,command=command,exit_code=code,seconds=time.monotonic()-start))
    (out/'checks.json').write_text(json.dumps(rows,indent=2)+'\n')
    print('END',name,code,flush=True)
    if code:raise SystemExit(code)
run('drive',['build-core/rf_tests'],'chain drive',{'RF_PLOT_DIR':str(out/'native')})
run('stationary',['build-core/rf_tests'],'chain wheel')
run('presets',['build-core/rf_tests'],'presets')
run('joints',['build-core/rf_tests'],'standard rigid: joints:')
run('audit-build',['c++','-std=c++17','-O3','-Isrc','-I.','tools/TorusGeometryAudit.cpp','build-core/librfsamples.a','build-core/librfcore.a','-pthread','-o',out/'audit'])
run('audit',[out/'audit',out/'sat',out/'native/chain-drive-poses.csv','-2','42'])
run('strict-build',['cmake','--build',strict,'--target','rfsamples','--parallel','2'],limited=False)
run('strict-test-build',['c++','-std=c++17','-O3','-ffp-contract=off','-DRF_STRICT_FP=1','-Wall','-Wextra','-Werror','-Isrc','-I.',out/'strict_main.cpp','tests/ChainWheelTests.cpp',strict/'librfsamples.a',strict/'librfcore.a','-pthread','-o',out/'strict-check'])
run('strict-drive',[out/'strict-check'],extra={'RF_PLOT_DIR':str(out/'strict')})
run('editor-build',['cmake','--build','/tmp/phys-realflow-editor-standard-build','--parallel','2'],timeout=360,limited=False)
run('readability',['python3','tools/readability.py','.'])
run('diff-check',['git','diff','--check'])
print('FINISHED',flush=True)
