"""From the repository root: python3 docs/book/evidence/chain-jitter-20261001/reproduce.py.

Linux/GCC diagnostic only. Writes to /tmp; runs sequentially on CPU 1-5 with
five SDK workers, nice 10, a 2 GiB address-space limit and per-run timeout.
Requires an existing Release build-core with the recorded numerical settings.
"""
import hashlib
import csv
from collections import defaultdict
import json
import os
from pathlib import Path
import resource
import shlex
import subprocess
import sys
import tempfile
import time

evidence = Path(__file__).resolve().parent
repo = evidence.parents[3]
output = Path(tempfile.mkdtemp(prefix="phys-chain-jitter-"))
env = dict(os.environ, RF_THREADS="5")
for key in ("RF_TEST", "RF_CHAIN_UNLOCK", "RF_CHAIN_NO_SHOCK", "RF_PLOT_DIR"):
    env.pop(key, None)


def limit():
    resource.setrlimit(resource.RLIMIT_AS, (2 * 1024**3, 2 * 1024**3))


for source, binary in (("jitter_probe.cpp", "jitter-probe"),
                       ("jitter_ab.cpp", "jitter-ab"),
                       ("stage_probe.cpp", "stage-probe")):
    subprocess.run([
        "nice", "-n", "10", "taskset", "-c", "1-5", "c++", "-std=c++17",
        "-O3", "-Wall", "-Wextra", "-Werror", "-fno-access-control", "-Isrc", "-I.",
        str(evidence / source), "build-core/librfsamples.a", "build-core/librfcore.a",
        "-pthread", "-o", str(output / binary),
    ], cwd=repo, check=True)

solver_source = (repo / "src/rigid/ContactSolver.cpp").read_text()
old_guard = "        if (p.positionBias <= 0) continue;"
new_guard = "        if (p.depth < 0) continue; // Keep zero-target rows for touching contacts (RP3D-style)."
if solver_source.count(old_guard) != 1:
    raise SystemExit("RP3D prototype requires the recorded baseline ContactSolver")
for variant in ("control", "rp3"):
    source = output / ("ContactSolver-" + variant + ".cpp")
    source.write_text(solver_source if variant == "control" else solver_source.replace(old_guard, new_guard))
    obj = output / ("ContactSolver-" + variant + ".o")
    subprocess.run(["nice", "-n", "10", "taskset", "-c", "1-5", "c++",
                    "-std=gnu++17", "-O3", "-DNDEBUG", "-march=native", "-Isrc",
                    "-c", str(source), "-o", str(obj)], cwd=repo, check=True)
    for probe, binary in (("jitter_probe.cpp", "jitter-" + variant),
                          ("jitter_ab.cpp", "jitter-ab-" + variant),
                          ("split_pair_probe.cpp", "split-pair-" + variant),
                          ("stage_probe.cpp", "stage-" + variant)):
        subprocess.run(["nice", "-n", "10", "taskset", "-c", "1-5", "c++",
                        "-std=c++17", "-O3", "-fno-access-control", "-Isrc", "-I.",
                        str(evidence / probe), str(obj), "build-core/librfsamples.a",
                        "build-core/librfcore.a", "-pthread", "-o", str(output / binary)],
                       cwd=repo, check=True)
    result = subprocess.check_output([str(output / ("split-pair-" + variant))], cwd=repo)
    (output / ("split-pair-" + variant + ".json")).write_bytes(result)

jobs = [("rest", "jitter-probe", ["1800", "0", "0"]),
        ("pull", "jitter-probe", ["1800", "1", "0"]),
        ("counterfactual-785", "stage-probe", ["785", "0"])]
for label, variant in (("beta005", "5"), ("iter24", "6")):
    for name, pull in (("rest", "0"), ("pull", "1")):
        jobs.append((label + "-" + name, "jitter-ab", ["1800", pull, variant]))
jobs.extend([("control-rest", "jitter-control", ["1800", "0", "0"]),
             ("rp3-rest", "jitter-rp3", ["1800", "0", "0"]),
             ("rp3-pull", "jitter-rp3", ["1800", "1", "0"]),
             ("rp3-stage-1011", "stage-rp3", ["1011", "1"]),
             ("rp3-iter24-rest", "jitter-ab-rp3", ["1800", "0", "6"]),
             ("rp3-iter24-pull", "jitter-ab-rp3", ["1800", "1", "6"])])

checks = []
for name, binary, args in jobs:
    cmd = ["timeout", "--kill-after=5s", "180s", "nice", "-n", "10",
           "taskset", "-c", "1-5", "/usr/bin/time", "-v", "-o",
           str(output / (name + ".resources")), str(output / binary),
           *args, str(output / (name + ".csv"))]
    start = time.monotonic()
    with (output / (name + ".json")).open("w") as out:
        with (output / (name + ".stderr.txt")).open("w") as err:
            result = subprocess.run(cmd, cwd=repo, env=env, stdout=out, stderr=err,
                                    preexec_fn=limit)
    checks.append(dict(name=name, command=cmd, exit_code=result.returncode,
                       seconds=time.monotonic() - start))
    (output / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(name, result.returncode, (output / (name + ".json")).read_text(), flush=True)
    if result.returncode:
        raise SystemExit(result.returncode)
    if name == "control-rest":
        reference = json.loads((output / "rest.json").read_text())
        control = json.loads((output / "control-rest.json").read_text())
        if control != reference:
            raise SystemExit("Control override differs from the baseline")

for name in ("counterfactual-785", "rp3-stage-1011"):
    with (output / (name + ".csv")).open() as file:
        rows = list(csv.DictReader(file))
    sums = defaultdict(lambda: dict(delta_T_J=0.0, delta_U_J=0.0, delta_E_J=0.0))
    for row in rows:
        for key in sums[row["stage"]]:
            sums[row["stage"]][key] += float(row[key])
    actual = json.loads((output / (name + ".json")).read_text())
    split = sums["split-counterfactual"]
    summary = dict(actual=actual, stage_totals=dict(sums), counterfactual_split_delta=split,
                   local_no_split_delta_E_J=actual["delta_E_J"] - split["delta_E_J"],
                   predicted_split_U_J=sum(float(row["predicted_split_U_J"]) for row in rows
                                           if row["stage"] == "rest-damping"),
                   max_corrected_points=max(int(row["corrected_points"]) for row in rows))
    (output / (name + "-summary.json")).write_text(json.dumps(summary, indent=2) + "\n")

files = [*(evidence / source for source in
           ("jitter_probe.cpp", "jitter_ab.cpp", "stage_probe.cpp", "split_pair_probe.cpp")),
         repo / "build-core/librfcore.a", repo / "build-core/librfsamples.a",
         *(output / binary for binary in ("jitter-probe", "jitter-ab", "stage-probe",
                                          "jitter-control", "jitter-rp3",
                                          "jitter-ab-control", "jitter-ab-rp3",
                                          "split-pair-control", "split-pair-rp3",
                                          "stage-control", "stage-rp3")),
         output / "ContactSolver-rp3.cpp", output / "ContactSolver-rp3.o"]
hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in files}
(output / "hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")

if "--regress" in sys.argv:
    commands = subprocess.check_output(["ninja", "-C", "build-core", "-t", "commands", "rf_tests"],
                                       cwd=repo, text=True).splitlines()
    line = next(line for line in reversed(commands) if " -o rf_tests " in line)
    link = shlex.split(line.split("&&")[1])
    link[link.index("-o") + 1] = str(output / "rf-tests-rp3")
    link.insert(link.index("librfcore.a"), str(output / "ContactSolver-rp3.o"))
    subprocess.run(link, cwd=repo / "build-core", check=True)
    for name, test_filter in (("default", None), ("hanging", "short torus chain"),
                              ("wheel", "chain wheel loaded"),
                              ("drive", "chain drive lifts"), ("nut", "nut bolt lifts")):
        test_env = dict(env, RF_JUNIT=str(output / ("rp3-" + name + ".xml")))
        if test_filter:
            test_env["RF_TEST"] = test_filter
        cmd = ["timeout", "--kill-after=5s", "180s", "nice", "-n", "10",
               "taskset", "-c", "1-5", "/usr/bin/time", "-v", "-o",
               str(output / ("rp3-" + name + ".resources")), str(output / "rf-tests-rp3")]
        with (output / ("rp3-" + name + ".log")).open("w") as out:
            result = subprocess.run(cmd, cwd=repo, env=test_env, stdout=out,
                                    stderr=subprocess.STDOUT, preexec_fn=limit)
        checks.append(dict(name="rp3-" + name, filter=test_filter,
                           command=cmd, exit_code=result.returncode))
        (output / "checks.json").write_text(json.dumps(checks, indent=2) + "\n")
        print(name, result.returncode, flush=True)
    hashes[str(output / "rf-tests-rp3")] = hashlib.sha256((output / "rf-tests-rp3").read_bytes()).hexdigest()
    (output / "hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")

print("Artifacts:", output)
