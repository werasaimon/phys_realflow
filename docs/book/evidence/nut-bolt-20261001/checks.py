import json, os, resource, subprocess, time
from pathlib import Path
out=Path('/tmp/phys-nut-bolt-20261001')
strict=Path('/tmp/phys-impact-research-20261001/build-strict')
rows=[]
def limited(): resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
def run(name,args,pattern=None,extra=None,timeout=180,limit=True):
    env=dict(os.environ,RF_THREADS='5',RF_JUNIT=str(out/(name+'.xml')))
    for key in ['RF_TEST','RF_PLOT_DIR','RF_CHAIN_UNLOCK','RF_CHAIN_NO_SHOCK']: env.pop(key,None)
    if pattern: env['RF_TEST']=pattern
    if extra: env.update(extra)
    cmd=['timeout','--kill-after=5s',str(timeout)+'s','nice','-n','10','taskset','-c','1-5','/usr/bin/time','-v','-o',str(out/(name+'.resources')),*map(str,args)]
    print('START',name,flush=True); start=time.monotonic()
    with (out/(name+'.output.txt')).open('w') as log:
        code=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limited if limit else None).returncode
    rows.append(dict(name=name,command=cmd,pattern=pattern,extra=extra,exit_code=code,seconds=time.monotonic()-start))
    (out/'final-checks.json').write_text(json.dumps(rows,indent=2)+'\n')
    print('END',name,code,flush=True)
    if code: raise SystemExit(code)
run('motor-final',['build-core/rf_tests'],'motor work')
run('nut-final',['build-core/rf_tests'],'nut bolt',{'RF_PLOT_DIR':str(out/'native')},240)
run('energy-build',['c++','-std=c++17','-O3','-Wall','-Wextra','-Werror','-Isrc','-I.','tools/ChainDriveStudy.cpp','build-core/librfsamples.a','build-core/librfcore.a','-pthread','-o',out/'chain-energy-final'])
run('energy-30s',[out/'chain-energy-final','1800','2','40','1',out/'chain-energy-30s.csv'],timeout=240)
run('energy-disabled',[out/'chain-energy-final','600','2','40','0',out/'chain-energy-disabled.csv'])
run('energy-heavy',[out/'chain-energy-final','600','3','40','1',out/'chain-energy-heavy.csv'])
run('energy-refined',[out/'chain-energy-final','600','2','80','1',out/'chain-energy-refined.csv'])
run('chain-final',['build-core/rf_tests'],'chain drive',{'RF_PLOT_DIR':str(out/'native')})
run('audit-build',['c++','-std=c++17','-O3','-Isrc','-I.','tools/TorusGeometryAudit.cpp','build-core/librfsamples.a','build-core/librfcore.a','-pthread','-o',out/'audit'])
run('nut-sat',[out/'audit',out/'nut-sat',out/'native/nut-bolt-poses.csv','-2','43'])
run('chain-sat',[out/'audit',out/'chain-sat',out/'native/chain-drive-poses.csv','-2','42'])
run('default',['build-core/rf_tests'],timeout=1000)
run('strict-build',['cmake','--build',strict,'--target','rfsamples','--parallel','2'],timeout=360,limit=False)
run('strict-test-build',['c++','-std=c++17','-O3','-ffp-contract=off','-DRF_STRICT_FP=1','-Wall','-Wextra','-Werror','-Isrc','-I.',out/'strict_checks.cpp','tests/MotorWorkTests.cpp','tests/NutBoltTests.cpp','tests/ChainWheelTests.cpp',strict/'librfsamples.a',strict/'librfcore.a','-pthread','-o',out/'strict-check'])
run('strict-tests',[out/'strict-check'],extra={'RF_PLOT_DIR':str(out/'strict')},timeout=360)
run('readability',['python3','tools/readability.py','.'])
run('diff-check',['git','diff','--check'])
print('FINISHED',flush=True)
