import os,json,subprocess,time,shutil,re
from pathlib import Path
p=Path(__file__).resolve().parent
env=dict(os.environ,RF_THREADS="5")
for k in ("RF_TEST","RF_PLOT_DIR","RF_CHAIN_UNLOCK","RF_CHAIN_NO_SHOCK"):env.pop(k,None)
log=(p/"editor-rp3.log").open("w")
cmd=["nice","-n","5","taskset","-c","1-5",str(p/"realflow-rp3"),"--preset","40","--size","1400x900"]
app=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
window=None
for i in range(40):
    if app.poll() is not None:raise SystemExit("Editor exited: "+str(app.returncode)+"; see "+str(p/"editor-rp3.log"))
    result=subprocess.run(["xprop","-root","_NET_CLIENT_LIST"],text=True,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
    for candidate in re.findall(r"0x[0-9a-fA-F]+",result.stdout):
        props=subprocess.run(["xprop","-id",candidate,"_NET_WM_PID","WM_NAME"],text=True,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL).stdout
        if f'= {app.pid}\n' in props:
            geometry=subprocess.run(["xwininfo","-id",candidate],text=True,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL).stdout
            if "IsViewable" in geometry:window=candidate;break
    if window:break
    time.sleep(0.25)
if window and shutil.which("wmctrl"):
    subprocess.run(["wmctrl","-ir",window,"-T","PhysRealFlow — цепь RP3D, эксперимент"],check=False)
    subprocess.run(["wmctrl","-ia",window],check=False)
time.sleep(3)
shot=p/"editor-rp3.png"
if window and shutil.which("import"):
    subprocess.run(["import","-window",window,str(shot)],timeout=10,check=False)
info=dict(pid=app.pid,window=window,preset=40,command=cmd,running=app.poll() is None,log=str(p/"editor-rp3.log"),screenshot=str(shot) if shot.exists() else None,solver="RP3D zero-target split rows; 12 iterations; sleeping off")
(p/"editor-rp3.json").write_text(json.dumps(info,ensure_ascii=False,indent=2)+"\n")
print(json.dumps(info,ensure_ascii=False,indent=2))
if not info["running"]:raise SystemExit(1)

if not window:raise SystemExit("Editor process runs but no mapped window was found")
print("Editor remains open until the user closes it.",flush=True)
code=app.wait();info["exit_code"]=code;info["running"]=False
(p/"editor-rp3.json").write_text(json.dumps(info,ensure_ascii=False,indent=2)+"\n")
print("Editor exit:",code,flush=True)
