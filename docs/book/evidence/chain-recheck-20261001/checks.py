import os, json, subprocess, resource, time, xml.etree.ElementTree as ET
from pathlib import Path
out=Path('/tmp/phys-chain-recheck-20261001')
rows=[]
def limit(): resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
def run(name,args,pattern=None,seconds=240):
    env=dict(os.environ,RF_THREADS='5',RF_JUNIT=str(out/(name+'.xml')),RF_PLOT_DIR=str(out))
    for key in ['RF_TEST','RF_CHAIN_UNLOCK','RF_CHAIN_NO_SHOCK']:env.pop(key,None)
    if pattern:env['RF_TEST']=pattern
    cmd=['timeout','--kill-after=5s',str(seconds)+'s','nice','-n','10','taskset','-c','1-5','/usr/bin/time','-v','-o',str(out/(name+'.resources')),*args]
    print('START',name,flush=True);start=time.monotonic()
    with (out/(name+'.output.txt')).open('w') as log:
        code=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,preexec_fn=limit).returncode
    row=dict(name=name,command=cmd,pattern=pattern,exit_code=code,seconds=time.monotonic()-start)
    if (out/(name+'.xml')).exists():row['junit']=ET.parse(out/(name+'.xml')).getroot().attrib
    rows.append(row);(out/'checks.json').write_text(json.dumps(rows,indent=2)+'\n')
    print('END',name,code,flush=True)
    if code:raise SystemExit(code)
run('chain-wheel',['build-core/rf_tests'],'chain wheel')
run('chain-drive',['build-core/rf_tests'],'chain drive')
run('hanging-chain',['build-core/rf_tests'],'short torus chain')
run('audit-build',['c++','-std=c++17','-O3','-Isrc','-I.','tools/TorusGeometryAudit.cpp','build-core/librfsamples.a','build-core/librfcore.a','-pthread','-o',str(out/'audit')])
run('wheel-sat',[str(out/'audit'),str(out/'wheel-sat'),str(out/'chain-wheel-poses.csv'),'-2','41'])
run('drive-sat',[str(out/'audit'),str(out/'drive-sat'),str(out/'chain-drive-poses.csv'),'-2','42'])
print('FINISHED',flush=True)
