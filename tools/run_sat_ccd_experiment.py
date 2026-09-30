#!/usr/bin/env python3
"""Build and measure GJK -> SAT + rotational CCD in an isolated copy of the dirty SDK.

python3 tools/run_sat_ccd_experiment.py --output /tmp/new-sat-study --frames 120
Optional --prepare-only prepares/builds/checks without simulation; --resume runs that prepared copy.
No repository solver changes, external engines, hidden geometry changes or golden-hash updates.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = {"baseline": {}, "sat": {"RF_EXPERIMENT_SAT": "1"},
            "sat-ccd": {"RF_EXPERIMENT_SAT": "1", "RF_EXPERIMENT_SAT_CCD": "1"}}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_once(path, old, new):
    text = path.read_text()
    if text.count(old) != 1:
        raise RuntimeError(f"patch anchor mismatch in {path}: {old[:80]}")
    path.write_text(text.replace(old, new))


def prepare(output):
    output.mkdir()
    source = output / "source"
    source.mkdir()
    for directory in ("src", "samples", "tests", "tools", "verification"):
        shutil.copytree(ROOT / directory, source / directory)
    for file in ("CMakeLists.txt", "AGENTS.md", "CLAUDE.md"):
        shutil.copy2(ROOT / file, source / file)
    manifest = {"root": str(ROOT), "core_sha256": digest(ROOT / "build-core/librfcore.a"),
                "baseline": {str(p.relative_to(source)): digest(p)
                             for p in (source / "src").rglob("*") if p.is_file()}}
    (output / "baseline.json").write_text(json.dumps(manifest, indent=2) + "\n")
    for file in ("SatCcd.h", "SatCcd.cpp"):
        shutil.copy2(source / "tools/rigid_reference" / file, source / "src/rigid" / file)
    patch(source)
    return source


def patch(source):
    narrow = source / "src/rigid/NarrowPhase.cpp"
    replace_once(narrow, '#include "rigid/NarrowPhase.h"',
                 '#include "rigid/NarrowPhase.h"\n#include "rigid/SatCcd.h"')
    anchor = "bool NarrowPhase::convexConvex(const PosedShape& A, const PosedShape& B, ContactManifold& m) {"
    replace_once(narrow, anchor, anchor + "\n    if (experimental::narrowEnabled() && experimental::isPolytope(A.shape) && experimental::isPolytope(B.shape))\n        return experimental::collideSat(A, B, m);")
    toi = source / "src/rigid/TimeOfImpact.cpp"
    replace_once(toi, '#include "rigid/TimeOfImpact.h"',
                 '#include "rigid/TimeOfImpact.h"\n#include "rigid/SatCcd.h"')
    anchor = "ToiResult timeOfImpact(const SweptPose& A, const SweptPose& B, float tol, int maxIt) {"
    replace_once(toi, anchor, anchor + "\n    if (experimental::ccdEnabled() && experimental::isPolytope(A.shape) && experimental::isPolytope(B.shape))\n        return experimental::conservativeToi(A, B, tol, maxIt);")
    anchor = "ToiResult partsTimeOfImpact(const SweptPose& A, const SweptPose& B, float tol) {"
    replace_once(toi, anchor, anchor + "\n    if (experimental::ccdEnabled()) return experimental::compoundToi(A, B, tol);")
    profile = source / "tools/TorusProfile.cpp"
    replace_once(profile, '#include "rigid/NarrowPhase.h"',
                 '#include "rigid/NarrowPhase.h"\n#include "rigid/SatCcd.h"')
    replace_once(profile, "    return 0; // this is a measurement tool, not a passing physics acceptance test",
                 "    experimental::writeStats(summary);\n    return 0; // this is a measurement tool, not a passing physics acceptance test")
    cmake = source / "CMakeLists.txt"
    replace_once(cmake, "add_library(rfcore STATIC", "add_library(rfcore STATIC\n  src/rigid/SatCcd.h src/rigid/SatCcd.cpp")
    with cmake.open("a") as stream:
        stream.write("\nadd_executable(rf_sat_checks tools/rigid_reference/SatCcdChecks.cpp)\ntarget_include_directories(rf_sat_checks PRIVATE tools/rigid_reference)\ntarget_link_libraries(rf_sat_checks PRIVATE rfsamples)\nadd_executable(rf_sat_profile tools/TorusProfile.cpp)\ntarget_link_libraries(rf_sat_profile PRIVATE rfsamples)\n")


def environment(variant):
    env = dict(os.environ, RF_THREADS="4")
    for key in ("RF_EXPERIMENT_SAT", "RF_EXPERIMENT_SAT_CCD", "RF_TEST"):
        env.pop(key, None)
    env.update(VARIANTS[variant])
    return env


def run(command, log, env, required=True, timeout=600):
    start = time.monotonic()
    with log.open("w") as stream:
        try:
            result = subprocess.run([str(x) for x in command], stdout=stream, stderr=subprocess.STDOUT, env=env, timeout=timeout)
            code, expired = result.returncode, False
        except subprocess.TimeoutExpired:
            code, expired = None, True
    record = {"command": [str(x) for x in command], "returncode": code, "timed_out": expired,
              "seconds": time.monotonic() - start, "log": str(log)}
    print(json.dumps(record), flush=True)
    if required and (expired or code):
        raise RuntimeError(f"command failed; see {log}")
    return record


def summarize(output):
    with (output / "frames.csv").open() as stream:
        frames = list(csv.DictReader(stream))
    with (output / "pair.csv").open() as stream:
        pairs = list(csv.DictReader(stream))
    wall = sorted(float(r["wall_ms"]) for r in frames)
    breaks = [int(r["substep"]) for r in pairs if abs(float(r["linking"])) != 1]
    fields = ("collide_ms", "solve_ms", "ccd_ms")
    return {"first_observed_broken_frame": next((int(r["frame"]) for r in frames
            if any(int(r[k]) for k in ("broken_hanging", "broken_released", "broken_bridge"))), -1),
            "worst_broken": max(sum(int(r[k]) for k in ("broken_hanging", "broken_released", "broken_bridge")) for r in frames),
            "pair_first_broken_substep": min(breaks, default=-1),
            "pair_max_raw_depth_m": max((float(r["raw_max_depth_m"]) for r in pairs), default=None),
            "frame_p50_ms": wall[len(wall)//2], "frame_p95_ms": wall[int((len(wall)-1)*.95)],
            "stage_mean_ms": {k: sum(float(r[k]) for r in frames)/len(frames) for k in fields},
            "summary": (output / "summary.txt").read_text()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--frames", type=int, default=120)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--variants", nargs="+", choices=list(VARIANTS), default=list(VARIANTS))
    args = parser.parse_args()
    if args.frames < 1:
        parser.error("frames must be positive")
    output = args.output.resolve()
    source = output / "source" if args.resume else prepare(output)
    build = output / "build"
    env = environment("baseline")
    run(["cmake", "-S", source, "-B", build, "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release",
         "-DRF_WERROR=ON", "-DRF_BUILD_VERIFY=OFF"], output / "configure.log", env)
    run(["cmake", "--build", build, "--parallel", "4", "--target", "rf_sat_checks", "rf_sat_profile", "rf_tests"], output / "build.log", env)
    run([build / "rf_sat_checks"], output / "analytic.log", env)
    if args.prepare_only:
        return
    report = {"frames": args.frames, "binary_sha256": digest(build / "rf_sat_profile"),
              "sources_sha256": {str(p.relative_to(source)): digest(p) for p in
                                  [source / "src/rigid/SatCcd.cpp", source / "src/rigid/SatCcd.h",
                                   source / "src/rigid/NarrowPhase.cpp", source / "src/rigid/TimeOfImpact.cpp",
                                   source / "tools/rigid_reference/SatCcdChecks.cpp"]}, "variants": {}}
    for variant in args.variants:
        env = environment(variant)
        folder = output / f"{variant}-{args.frames}frames"
        result = run([build / "rf_sat_profile", "--output", folder, "--frames", args.frames,
                      "--mode", "substep"], output / f"{variant}-{args.frames}frames.log", env)
        cold = run([build / "rf_tests"], output / f"{variant}-cold.log",
                   dict(env, RF_TEST="intersecting stationary"), required=False, timeout=120)
        regression = run([build / "rf_tests"], output / f"{variant}-cups.log",
                         dict(env, RF_TEST="six cups"), required=False, timeout=120)
        report["variants"][variant] = {"regression": regression, "cold": cold, "run": result, **summarize(folder)}
        (output / f"report-{args.frames}frames.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
