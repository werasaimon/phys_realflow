#!/usr/bin/env python3
"""Reproduce the small analytic examples in docs/book, using only the standard library.

Run: python3 tools/check_book_math.py
This checks the book's calculations, not the production C++ implementation or CCD safety.
The angular-lock case deliberately demonstrates energy creation by the current formula.
"""
from fractions import Fraction as F
import json
import math


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def mul(matrix, vector):
    return [dot(row, vector) for row in matrix]


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def angular_lock():
    off = F(2, 5)
    inverse = [[F(1), off, F(0)], [off, F(1), F(0)], [F(0), F(0), F(1)]]
    det = 1 - off * off
    inertia = [[1/det, -off/det, F(0)], [-off/det, 1/det, F(0)], [F(0), F(0), F(1)]]
    moments = [F(5, 7), F(1), F(5, 3)]
    axes = [[1, 1, 0], [0, 0, 1], [1, -1, 0]]
    for moment, axis in zip(moments, axes):
        require(mul(inertia, axis) == [moment*x for x in axis], 'principal inertia')
    require(min(moments) > 0 and 2*max(moments) <= sum(moments), 'physical inertia')
    spin, normal = [F(1, 10), F(0), F(0)], [F(0), F(1), F(0)]
    old_impulse = [-x * dot(normal, mul(inertia, spin)) for x in normal]
    scalar = -dot(normal, spin) / dot(normal, mul(inverse, normal))
    gain = dot(spin, old_impulse) + F(1, 2)*dot(old_impulse, mul(inverse, old_impulse))
    require(old_impulse == [0, F(1, 21), 0] and scalar == 0, 'angular impulse')
    require(gain == F(1, 882), 'angular energy gain')
    return {'old_gain_j': float(gain), 'scalar_impulse': float(scalar)}


def sat_containment():
    a, b = (-2, 2), (-1, 1)
    overlap_length = min(a[1], b[1]) - max(a[0], b[0])
    separating_shift = min(a[1]-b[0], b[1]-a[0])
    require(overlap_length == 2 and separating_shift == 3, 'SAT containment')
    return {'overlap_length': overlap_length, 'separating_shift': separating_shift}


def rotation_arc():
    radius = 0.1
    checked = 0
    for theta in (0.01, 0.2, math.pi, 2*math.pi):
        bound = radius * theta * theta / 8
        for k in range(101):
            s = k / 100
            point = (radius*math.cos(s*theta), radius*math.sin(s*theta))
            chord = (radius*(1-s+s*math.cos(theta)), radius*s*math.sin(theta))
            error = math.hypot(point[0]-chord[0], point[1]-chord[1])
            require(error <= bound + 1e-15, 'sampled rotation arc bound')
            checked += 1
    return {'samples': checked, 'bound_for_0_2_rad_m': radius*0.2**2/8}


def friction_square():
    limit = F(3, 10)*F(2)
    squared_norm = limit**2 + limit**2
    require(squared_norm == 2*limit**2, 'square corner exceeds Coulomb circle')
    return {'corner_to_circle_ratio': math.sqrt(2)}


def inelastic_chain():
    mass, before = [F(1), F(2), F(3)], [F(3), F(0), F(0)]
    impulses = [F(5, 2), F(3, 2)]
    after = [before[0]-impulses[0]/mass[0],
             before[1]+(impulses[0]-impulses[1])/mass[1], before[2]+impulses[1]/mass[2]]
    energy = lambda velocity: sum(m*v*v/2 for m, v in zip(mass, velocity))
    require(after == [F(1, 2)]*3 and dot(mass, before) == dot(mass, after), 'chain momentum')
    require(energy(before) == F(9, 2) and energy(after) == F(3, 4), 'chain energy')
    return {'energy_before_j': float(energy(before)), 'energy_after_j': float(energy(after))}


def pair_momentum():
    point, impulse = [1, 2, 3], [2, -3, 5]
    centres = [[-1, 0, 2], [3, 1, -2]]
    total = [0, 0, 0]
    for sign, centre in zip((1, -1), centres):
        applied = [sign*x for x in impulse]
        orbital = cross(centre, applied)
        spin = cross([p-x for p, x in zip(point, centre)], applied)
        total = [a+b+c for a, b, c in zip(total, orbital, spin)]
    require(total == [0, 0, 0], 'internal angular momentum balance')
    return {'angular_momentum_change': total}


def main():
    checks = [angular_lock, sat_containment, rotation_arc, friction_square, inelastic_chain, pair_momentum]
    results = {check.__name__: check() for check in checks}
    print(json.dumps({'passed': len(results), 'scope': 'book examples only', 'results': results}, indent=2))


if __name__ == '__main__':
    main()
