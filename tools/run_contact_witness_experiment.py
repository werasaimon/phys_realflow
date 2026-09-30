#!/usr/bin/env python3
"""Measure true surface witnesses + per-point seam filtering in an isolated SDK copy.

python3 tools/run_contact_witness_experiment.py --output /tmp/new-witness-study --frames 120
Baseline and candidate use the same binary, scene and settings. No repository solver modifications.
"""
import argparse
import json
from pathlib import Path
import shutil
from run_sat_ccd_experiment import ROOT, digest, replace_once, run, summarize, environment


def patch(source, face_normal=False):
    header = source / "src/rigid/NarrowPhase.h"
    replace_once(header, "    float depth = 0; // > 0 penetration, < 0 speculative gap (within the contact margin)",
                 """    float depth = 0; // > 0 penetration, < 0 speculative gap (within the contact margin)
    Vector3 pointA{0},pointB{0};
    bool witnesses=false;
    static ContactPoint surface(Vector3 a,Vector3 b,Vector3 n) {
        ContactPoint p{(a+b)*0.5f,n,dot(b-a,n)};
        p.pointA=a; p.pointB=b; p.witnesses=true; return p;
    }""")
    replace_once(header, "class NarrowPhase {", "bool witnessExperiment();\n\nclass NarrowPhase {")
    narrow = source / "src/rigid/NarrowPhase.cpp"
    replace_once(narrow, '#include "rigid/NarrowPhase.h"', '#include "rigid/NarrowPhase.h"\n#include <cstdlib>')
    replace_once(narrow, "namespace rf {", """namespace rf {
bool witnessExperiment() {
    static const bool on=[] { const char* v=std::getenv("RF_EXPERIMENT_WITNESS"); return v && v[0]=='1'; }();
    return on;
}""")
    replace_once(narrow, "    for (const ContactPoint& c : src.points) dst.add(c.position, -c.normal, c.depth);", """    for (const ContactPoint& c : src.points) {
        if(!c.witnesses) { dst.add(c.position,-c.normal,c.depth); continue; }
        auto p=c; p.normal=-p.normal; std::swap(p.pointA,p.pointB); dst.points.push_back(p);
    }""")
    replace_once(narrow, "        pts.push_back({(c1 + c2) * 0.5f, n, depth});", "        pts.push_back(witnessExperiment() ? ContactPoint::surface(c1,c2,n) : ContactPoint{(c1+c2)*0.5f,n,depth});")
    replace_once(narrow, "    clipIncidentFace(ref, nr, inc, clipped);", """    // Side planes parallel to n preserve polygon membership under projection onto the reference plane.
    clipIncidentFace(ref,witnessExperiment() ? n : nr,inc,clipped);""")
    replace_once(narrow, "        const float d = depths[i] + shift;", "        const float d = depths[i] + (witnessExperiment() ? 0 : shift);")
    replace_once(narrow, "        pts.push_back({inc[i] + n * (refIsB ? 0.5f * d : -0.5f * d), n, d}); // midway between the surfaces", """        if(witnessExperiment()) {
            const Vector3 a=refIsB ? inc[i] : inc[i]-n*d;
            const Vector3 b=refIsB ? inc[i]+n*d : inc[i];
            pts.push_back(ContactPoint::surface(a,b,n));
        } else pts.push_back({inc[i]+n*(refIsB ? 0.5f*d : -0.5f*d),n,d});""")
    # Set exact GJK/EPA witnesses before any face/perturbation fallback.
    replace_once(narrow, "    // 1) Boundary simplices: clip the supporting faces (exact multi-point manifold).", """    if(witnessExperiment()) {
        if(!g.intersect) single=ContactPoint::surface(g.pointA,g.pointB,n);
        else {
            const auto pr=epa(A,B,g);
            if(pr.valid) single=ContactPoint::surface(pr.pointA,pr.pointB,n);
        }
    }
    // 1) Boundary simplices: clip the supporting faces (exact multi-point manifold).""")
    replace_once(narrow, "        if (!dup) pts.push_back({pos, n, depth});", "        if(!dup) pts.push_back(witnessExperiment() ? ContactPoint::surface(pa,pb,n) : ContactPoint{pos,n,depth});")
    if face_normal:
        replace_once(narrow, "bool NarrowPhase::faceManifold(const PosedShape& A, const PosedShape& B, const Vector3& n, float depthLow, float depthHigh,\n                               std::vector<ContactPoint>& pts) {",
                     "bool NarrowPhase::faceManifold(const PosedShape& A, const PosedShape& B, const Vector3& axis, float depthLow, float depthHigh,\n                               std::vector<ContactPoint>& pts) {\n    Vector3 n=axis;")
        replace_once(narrow, "    // Side planes parallel to n preserve polygon membership under projection onto the reference plane.",
                     "    // For a point/plane gap, its translational gradient is the actual plane normal.\n    if(witnessExperiment()) n=refIsB ? nr : -nr;\n    // Side planes parallel to n preserve polygon membership under projection onto the reference plane.")
    contact = source / "src/rigid/ContactSolver.cpp"
    old = "            if (inOtherPart(S.partsA, S.boundsA, u, onA - n * kSeamProbe) || inOtherPart(S.partsB, S.boundsB, v, onB + n * kSeamProbe)) continue;"
    replace_once(contact, old, """            if(witnessExperiment()) {
                part.points.erase(std::remove_if(part.points.begin(),part.points.end(),[&](const ContactPoint& point) {
                    const Vector3 a=point.witnesses ? point.pointA : point.position-point.normal*(0.5f*point.depth);
                    const Vector3 b=point.witnesses ? point.pointB : point.position+point.normal*(0.5f*point.depth);
                    return inOtherPart(S.partsA,S.boundsA,u,a-point.normal*kSeamProbe) ||
                           inOtherPart(S.partsB,S.boundsB,v,b+point.normal*kSeamProbe);
                }),part.points.end());
                if(part.points.empty()) continue;
            } else if(inOtherPart(S.partsA,S.boundsA,u,onA-n*kSeamProbe) || inOtherPart(S.partsB,S.boundsB,v,onB+n*kSeamProbe)) continue;""")
    observer = source / "tools/CompoundWitnessAudit.cpp"
    replace_once(observer, "            const auto onA=c.position-na*(.5f*c.depth),onB=c.position+na*(.5f*c.depth);",
                 "            const auto onA=c.witnesses ? c.pointA : c.position-na*(.5f*c.depth),onB=c.witnesses ? c.pointB : c.position+na*(.5f*c.depth);")
    with (source / "CMakeLists.txt").open("a") as f:
        f.write("\nadd_executable(rf_witness_profile tools/TorusProfile.cpp)\ntarget_link_libraries(rf_witness_profile PRIVATE rfsamples)\nadd_executable(rf_witness_audit tools/CompoundWitnessAudit.cpp)\ntarget_link_libraries(rf_witness_audit PRIVATE rfsamples)\nadd_executable(rf_face_checks tools/FaceWitnessChecks.cpp)\ntarget_link_libraries(rf_face_checks PRIVATE rfcore)\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--frames", default=120, type=int)
    parser.add_argument("--face-normal", action="store_true", help="use actual reference-plane normal in face contacts")
    args = parser.parse_args()
    if args.frames < 1:
        parser.error("positive frames required")
    output = args.output.resolve(); output.mkdir()
    source = output / "source"; source.mkdir()
    for directory in ("src", "samples", "tests", "tools", "verification"):
        shutil.copytree(ROOT / directory, source / directory)
    for file in ("CMakeLists.txt", "AGENTS.md", "CLAUDE.md"):
        shutil.copy2(ROOT / file, source / file)
    baseline = {str(p.relative_to(source)): digest(p) for p in (source / "src").rglob("*") if p.is_file()}
    (output / "baseline.json").write_text(json.dumps(baseline, indent=2) + "\n")
    patch(source, args.face_normal); build = output / "build"; env = environment("baseline")
    env.pop("RF_EXPERIMENT_WITNESS", None)
    run(["cmake", "-S", source, "-B", build, "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release", "-DRF_WERROR=ON", "-DRF_BUILD_VERIFY=OFF"], output / "configure.log", env)
    run(["cmake", "--build", build, "--parallel", "4", "--target", "rf_tests", "rf_witness_profile", "rf_witness_audit", "rf_face_checks"], output / "build.log", env)
    report = {"hypothesis": "True clipped witnesses and individual seam filtering remove synthetic anchors without losing whole patches.",
              "limits": "Retains finite seam probe, patch reduction, midpoint impulses and existing CCD; not exact CSG boundary or island rollback.",
              "frames": args.frames, "reference_face_normal": args.face_normal, "binary_sha256": digest(build / "rf_witness_profile"),
              "sources_sha256": {str(p.relative_to(source)): digest(p) for p in (source / "src").rglob("*") if p.is_file()}, "variants": {}}
    for name, flag in (("baseline", "0"), ("witness", "1")):
        active = dict(env, RF_EXPERIMENT_WITNESS=flag)
        geometry = run([build / "rf_face_checks"] + (["--expect-consistent"] if args.face_normal and flag == "1" else []), output / (name + "-geometry.log"), active)
        run([build / "rf_witness_audit", output / (name + "-audit")], output / (name + "-audit.log"), active)
        cold = run([build / "rf_tests"], output / (name + "-cold.log"), dict(active, RF_TEST="intersecting stationary"), required=False, timeout=120)
        cups = run([build / "rf_tests"], output / (name + "-cups.log"), dict(active, RF_TEST="six cups"), required=False, timeout=120)
        folder = output / (name + "-frames")
        profile = run([build / "rf_witness_profile", "--output", folder, "--frames", args.frames, "--mode", "substep"], output / (name + "-frames.log"), active)
        report["variants"][name] = {"geometry": geometry, "cold": cold, "cups": cups, "run": profile, **summarize(folder)}
        (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print("Witness experiment complete", flush=True)


if __name__ == "__main__":
    main()
