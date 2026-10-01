import json
import os
from pathlib import Path
import re
import subprocess
import time

binary = Path('/tmp/phys-realflow-editor-standard-build/realflow')
log_path = Path('/tmp/phys-chain-drive-20261001/editor-chain.log')
env = dict(os.environ, RF_THREADS='5')
with log_path.open('w') as log:
    process = subprocess.Popen(
        ['nice', '-n', '5', 'taskset', '-c', '1-5', str(binary), '--preset', '42', '--size', '1400x900'],
        env=env, stdout=log, stderr=subprocess.STDOUT,
        stdin=subprocess.DEVNULL, start_new_session=True)

deadline = time.monotonic() + 20
window = None
while time.monotonic() < deadline:
    if process.poll() is not None:
        raise SystemExit(f'Editor exited {process.returncode}: {log_path.read_text()}')
    clients = subprocess.run(['xprop', '-root', '_NET_CLIENT_LIST'],
                             text=True, capture_output=True)
    for candidate in re.findall(r'0x[0-9a-fA-F]+', clients.stdout):
        props = subprocess.run(['xprop', '-id', candidate, '_NET_WM_PID', 'WM_NAME'],
                               text=True, capture_output=True).stdout
        if f'= {process.pid}\n' in props:
            info = subprocess.run(['xwininfo', '-id', candidate],
                                  text=True, capture_output=True).stdout
            if 'IsViewable' in info:
                window = candidate
                break
    if window:
        break
    time.sleep(0.3)
if not window:
    raise SystemExit(f'No mapped editor window found for PID {process.pid}; log: {log_path}')

shots = []
resources = []
for index in range(2):
    time.sleep(3)
    shot = Path(f'/tmp/phys-chain-drive-20261001/editor-chain-{index + 1}.png')
    subprocess.run(['import', '-window', window, str(shot)], check=True, timeout=10)
    shots.append(str(shot))
    status_path = Path(f"/proc/{process.pid}/status")
    if status_path.exists():
        fields = status_path.read_text().splitlines()
        resources.append({line.split(":", 1)[0]: line.split(":", 1)[1].strip() for line in fields if line.startswith(("VmRSS:", "VmHWM:", "Threads:"))})
result = {'pid': process.pid, 'window': window, 'preset': 42,
          'binary': str(binary), 'log': str(log_path), 'screenshots': shots,
          'running': process.poll() is None, 'resources': resources}
Path('/tmp/phys-chain-drive-20261001/editor-chain.json').write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
