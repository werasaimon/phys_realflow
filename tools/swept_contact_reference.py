#!/usr/bin/env python3
"""Analytic research reference for bounded penetration along a rigid trajectory.

python3 tools/swept_contact_reference.py --output /tmp/swept-reference-new.json
Spheres may be offset from a rotating body's origin. This is a trajectory guard,
not an integrator, general mesh CCD, IPC solver, or a floating-point certificate.
The bound below follows the speed bound and the reverse triangle inequality.
"""
import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path


@dataclass(frozen=True)
class Motion:
    origin: tuple
    translation: tuple
    offset: tuple
    angle: float
    radius: float

    def point(self, t):
        c, s = math.cos(self.angle*t), math.sin(self.angle*t)
        x, y = self.offset
        return (self.origin[0]+t*self.translation[0]+c*x-s*y,
                self.origin[1]+t*self.translation[1]+s*x+c*y)

    def speed_bound(self):
        return math.hypot(*self.translation)+abs(self.angle)*math.hypot(*self.offset)


def gap(a, b, t):
    x, y = a.point(t), b.point(t)
    return math.dist(x, y)-a.radius-b.radius


def guard(a, b, epsilon, budget=4096, roundoff=1e-12):
    # |g(t)-g(mid)| <= L * |t-mid|. Hence this bound covers each whole interval,
    # including rotation; endpoint/chord-only CCD cannot establish the same claim.
    speed = a.speed_bound()+b.speed_bound()
    pending = [(0., 1.)]
    calls, minimum = 0, math.inf
    while pending and calls < budget:
        left, right = pending.pop()
        mid = (left+right)/2
        value = gap(a, b, mid)
        calls += 1
        minimum = min(minimum, value)
        if value < -epsilon-roundoff:
            return {"status": "violation", "witness_t": mid, "witness_gap_m": value, "queries": calls}
        if value-speed*(right-left)/2-roundoff >= -epsilon:
            continue
        if mid == left or mid == right:
            return {"status": "unresolved", "queries": calls}
        pending.extend([(mid, right), (left, mid)])
    return {"status": "unresolved" if pending else "admissible", "queries": calls,
            "sampled_min_gap_m": minimum}


def validate():
    fixed = Motion((0., 0.), (0., 0.), (0., 0.), 0., .1)
    crossing = Motion((-1., 0.), (2., 0.), (0., 0.), 0., .1)
    rotating = Motion((0., 0.), (0., 0.), (1., 0.), 2*math.pi, .05)
    obstacle = Motion((0., 1.), (0., 0.), (0., 0.), 0., .05)
    near = Motion((.1997, 0.), (0., 0.), (0., 0.), 0., .1)
    cases = [
        ("separated_motion", Motion((-1., 1.), (2., 0.), (0., 0.), 0., .1), fixed, .0005, 4096, "admissible"),
        ("translational_tunneling", crossing, fixed, .0005, 4096, "violation"),
        ("full_rotation_same_endpoints", rotating, obstacle, .0005, 4096, "violation"),
        ("0.3mm_overlap_allowed_at_0.5mm", near, fixed, .0005, 4096, "admissible"),
        ("0.3mm_overlap_rejected_at_0.1mm", near, fixed, .0001, 4096, "violation"),
        ("budget_exhaustion_is_not_success", rotating, obstacle, .0005, 1, "unresolved"),
    ]
    report = []
    for name, a, b, epsilon, budget, expected in cases:
        result = guard(a, b, epsilon, budget)
        reverse = guard(b, a, epsilon, budget)
        if result["status"] != expected or result != reverse:
            raise AssertionError((name, result, reverse, expected))
        if "witness_t" in result and gap(a, b, result["witness_t"]) >= -epsilon:
            raise AssertionError("invalid rejection witness")
        if name in ("translational_tunneling", "full_rotation_same_endpoints"):
            if min(gap(a,b,0),gap(a,b,1)) <= 0:
                raise AssertionError("fixture must fool endpoint-only detection")
        report.append({"case": name, "epsilon_m": epsilon, "expected": expected, **result})
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = {"scope": "analytic offset-sphere trajectories; no dynamics or mesh CCD",
              "roundoff_note": "1e-12 guard is heuristic, not directed interval arithmetic",
              "cases": validate()}
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print("6 analytic trajectory checks PASS (including allowed overlap and unresolved budget)")


if __name__ == "__main__":
    main()
