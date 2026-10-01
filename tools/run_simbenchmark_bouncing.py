#!/usr/bin/env python3
"""Run an adapted SimBenchmark Bouncing comparison, without changing solver settings.

RF_THREADS=10 nice -n 10 taskset -c 1-5,7-11 python3 tools/run_simbenchmark_bouncing.py \
    --binary /tmp/phys-engine-reference-build/rf_simbenchmark_bouncing \
    --upstream /tmp/phys-ccd-retry-20261001 --output /tmp/simbenchmark-bouncing-new
Raw traces stay in the new output directory. Every failure is retained in report.json.
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


AUTHOR_COMMIT = "3fd17646352f4c1c2c3108f93dae00f2a2d60097"
STEPS = ("0.01", "0.005", "0.002", "0.001", "0.0005")
ENGINES = ("sdk", "bullet", "jolt", "box3d", "box3d-native", "sdk-ccd")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_trace(directory, summary):
    """Recompute the author's MSE from rounded CSV; check completion and finiteness."""
    with (directory / "frames.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != summary["frames"] or not summary["complete"]:
        raise ValueError("incomplete simulation")
    squared = 0.0
    worst = 0.0
    e0 = summary["reference_energy_J"]
    for index, row in enumerate(rows):
        if int(row["frame"]) != index or not all(math.isfinite(float(v)) for v in row.values()):
            raise ValueError("invalid trace")
        if abs(float(row["time_s"]) - index * summary["dt_s"]) > 3e-6:
            raise ValueError("trace consumed the wrong time interval")
        error = float(row["energy_J"]) - e0
        squared += error * error
        worst = max(worst, abs(error))
    mse = squared / len(rows)
    # CSV energies are binary32 round trips: <= half ULP, plus decimal spelling of that value.
    max_energy = max(abs(float(row["energy_J"])) for row in rows)
    quantization = 2 ** (math.frexp(max_energy)[1] - 25) + 0.00005
    mse_bound = 2 * worst * quantization + quantization ** 2 + abs(mse) * 1e-6
    if abs(mse - summary["energy_mse_J2"]) > mse_bound:
        raise ValueError("MSE disagrees with trace beyond float serialization error")
    if abs(summary["accuracy_accepted_time_s"] - 20) > 3e-6:
        raise ValueError("accuracy pass did not consume 20 seconds")
    if summary["accuracy_accepted_time_s"] != summary["timing_accepted_time_s"]:
        raise ValueError("timing pass consumed a different interval")
    return {"trace_sha256": digest(directory / "frames.csv"), "trace_rows": len(rows),
            "recomputed_energy_mse_J2": mse, "mse_rounding_bound_J2": mse_bound}


def run_one(args, engine, dt):
    name = engine + "-" + dt
    directory = args.output / name
    command = [str(args.binary), engine, dt, str(directory), "20"]
    entry = {"engine": engine, "requested_dt_s": float(dt), "command": command}
    entry["binary_sha256"] = digest(args.binary)
    if entry["binary_sha256"] != args.binary_sha256:
        entry["failure"] = "binary changed before the run"
        return entry
    started = time.monotonic()
    with (args.output / (name + ".log")).open("w") as log:
        try:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                    timeout=args.timeout, check=False)
            entry["exit_code"] = result.returncode
        except subprocess.TimeoutExpired:
            entry.update(exit_code=None, failure="timeout")
    entry["process_wall_s"] = time.monotonic() - started
    if digest(args.binary) != args.binary_sha256:
        entry["failure"] = "binary changed during the run"
    entry["log_sha256"] = digest(args.output / (name + ".log"))
    if entry["exit_code"] == 0:
        try:
            summary = json.loads((directory / "summary.json").read_text())
            entry.update(summary=summary, audit=audit_trace(directory, summary))
        except (OSError, ValueError, KeyError) as error:
            entry["failure"] = str(error)
    else:
        entry.setdefault("failure", "process failed; inspect log")
    return entry


def upstream_manifest(directory):
    tree = json.loads((directory / "simbenchmark-tree.json").read_text())
    if tree["sha"] != AUTHOR_COMMIT:
        raise ValueError("unexpected upstream tree commit")
    names = ("bouncing.yaml", "BouncingBenchmark.hpp", "BouncingBullet.cpp",
             "bouncing-ball.urdf", "bouncing-plane.urdf")
    entries = []
    for name in names:
        if name.startswith("bouncing-"):
            suffix = "/ball.urdf" if "ball" in name else "/plane.urdf"
            candidates = [p for p in tree["tree"] if "bouncing" in p["path"].lower() and p["path"].endswith(suffix)]
        elif name == "BouncingBullet.cpp":
            candidates = [p for p in tree["tree"] if p["path"] == "benchmark/bulletSim/BouncingBenchmark.cpp"]
        else:
            candidates = [p for p in tree["tree"] if p["path"].endswith("/" + name)]
        if len(candidates) != 1:
            raise ValueError("ambiguous upstream source: " + name)
        path = directory / name
        data = path.read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        if blob != candidates[0]["sha"]:
            raise ValueError("upstream blob hash mismatch: " + name)
        entries.append({"author_path": candidates[0]["path"], "git_blob_sha1": blob,
                        "sha256": digest(path), "bytes": len(data)})
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--upstream", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--engines", nargs="+", choices=ENGINES + ("sdk-variational",), default=ENGINES)
    parser.add_argument("--provenance", type=Path, help="JSON with reference library versions/build receipts")
    args = parser.parse_args()
    args.binary = args.binary.resolve(strict=True)
    args.binary_sha256 = digest(args.binary)
    if args.timeout <= 0:
        parser.error("timeout must be positive")
    upstream = upstream_manifest(args.upstream)
    args.output.mkdir(exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    sources = [Path(__file__), root / "tools/SimBenchmarkBouncing.cpp",
               *sorted(root.glob("tools/rigid_reference/Bouncing*.h")),
               root / "tools/rigid_reference/JoltRuntime.h",
               root / "tools/rigid_reference/CMakeLists.txt", root / "build-core/librfcore.a",
               *sorted(root.glob("src/rigid/Variational*")),
               root / "src/rigid/ConservativeAdvancement.cpp", root / "src/rigid/TimeOfImpact.h"]
    report = {"benchmark": "adapted SimBenchmark Bouncing", "author_commit": AUTHOR_COMMIT,
              "engines": args.engines,
              "upstream_sources": upstream, "binary_sha256": args.binary_sha256,
              "source_sha256": {str(p.relative_to(root)): digest(p) for p in sources},
              "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
              "platform": platform.platform(), "rf_threads": os.getenv("RF_THREADS", "default"),
              "timestep_sweep_s": [float(s) for s in STEPS], "iterations": 20,
              "metric": "mean((E_before_step - 24034.5 J)^2)",
              "additional_metrics": ["relative energy RMSE", "analytic position/velocity RMSE",
                                     "sampled penetration", "separate timing pass wall time"],
              "runs": []}
    if args.provenance:
        report["reference_provenance"] = json.loads(args.provenance.read_text())
    if hasattr(os, "sched_getaffinity"):
        report["allowed_logical_cpus"] = sorted(os.sched_getaffinity(0))
    for dt in STEPS:
        for engine in args.engines:
            print("Starting", engine, dt, flush=True)
            entry = run_one(args, engine, dt)
            report["runs"].append(entry)
            (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
            print(engine, dt, "FAILED: " + entry["failure"] if "failure" in entry
                  else "energy RMSE=" + str(entry["summary"]["relative_energy_rmse"]), flush=True)
    if any("failure" in entry for entry in report["runs"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
