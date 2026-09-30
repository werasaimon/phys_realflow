#!/usr/bin/env python3
"""Run native torus references sequentially; separate stepping and observer timings.

python3 tools/profile_rigid_references.py --output /tmp/rigid-profile-new \
    --references /tmp/phys-engine-reference-build --sdk /tmp/rf_torus_profile_current
Physics failures (exit 1) are data. Import/runtime errors stop the experiment.
These runs compare specific configurations, not engine quality or equal workloads.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import time


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def distribution(values):
    values = sorted(values)
    return {"mean": sum(values) / len(values), "p50": values[len(values) // 2],
            "p95": values[math.ceil(.95 * len(values)) - 1], "max": values[-1]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--references", required=True, type=Path)
    parser.add_argument("--sdk", required=True, type=Path)
    parser.add_argument("--frames", type=int, default=120)
    args = parser.parse_args()
    if args.frames < 1:
        parser.error("frames must be positive")
    args.output.mkdir(exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    configs = [("sdk", [str(args.sdk), "--output", "OUT", "--frames", str(args.frames)]),
               ("bullet", [str(args.references / "rf_bullet_reference"), "OUT", "6", str(args.frames), "1"]),
               ("jolt", [str(args.references / "rf_jolt_reference"), "OUT", "6", str(args.frames), ".004", "1"]),
               ("box3d_soft", [str(args.references / "rf_box3d_reference"), "OUT", "10", str(args.frames), "0"]),
               ("box3d_full", [str(args.references / "rf_box3d_reference"), "OUT", "10", str(args.frames), "1"])]
    sources = [Path(__file__), root / "tools/TorusProfile.cpp", *root.glob("tools/*TorusReference.cpp"),
               root / "tools/rigid_reference/ReferenceTiming.h", root / "build-core/librfcore.a"]
    report = {"platform": platform.platform(), "frames": args.frames, "rf_threads": 4,
              "source_sha256": {str(p): digest(p) for p in sources}, "runs": []}
    for name, command in configs:
        command = [str(args.output / name) if word == "OUT" else word for word in command]
        print("Starting", name, flush=True)
        started = time.monotonic()
        with (args.output / (name + ".log")).open("w") as log:
            result = subprocess.run(command, cwd=root, env={**os.environ, "RF_THREADS": "4"},
                                    stdout=log, stderr=subprocess.STDOUT, check=False)
        elapsed = time.monotonic() - started
        if result.returncode not in (0, 1):
            raise RuntimeError(f"{name} failed to run: {result.returncode}; inspect log")
        with (args.output / name / "frames.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        if len(rows) != args.frames:
            raise RuntimeError(f"{name}: incomplete CSV")
        timings = {k: distribution([float(row[k]) for row in rows]) for k in rows[0]
                   if k.endswith("_ms")}
        entry = {"name": name, "command": command, "binary_sha256": digest(command[0]),
                 "exit_code": result.returncode, "process_wall_s": elapsed, "timings_ms": timings,
                 "summary": (args.output / name / "summary.txt").read_text()}
        report["runs"].append(entry)
        (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
        print(name, "exit", result.returncode, "wall_s", round(elapsed, 2), flush=True)


if __name__ == "__main__":
    main()
