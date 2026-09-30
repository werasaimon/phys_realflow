#!/usr/bin/env python3
"""Independent normal-impulse reference at fixed rigid-body poses, in SI units.

Usage: python3 tools/contact_projection_reference.py --output NEW_FILE.json
Enumerates active sets for <=12 touching contacts; double precision, no CFM.
This is an instantaneous frictionless, perfectly inelastic reference, not a
geometry detector, integrator, CCD implementation, or replacement for rfcore.
Background: https://www.cse.lehigh.edu/~trink/Papers/STicra00.pdf
The mass-metric projection and its energy balance are derived in
docs/18-rigid-contact-model.md.
"""
import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def linear_solve(matrix, rhs):
    n = len(rhs)
    if n == 0:
        return []
    work = [list(row)+[value] for row, value in zip(matrix, rhs)]
    scale = max(abs(x) for row in matrix for x in row)
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(work[i][k]))
        if abs(work[pivot][k]) <= 1e-12 * max(scale, 1e-30):
            return None
        work[k], work[pivot] = work[pivot], work[k]
        divisor = work[k][k]
        work[k] = [x/divisor for x in work[k]]
        for i in range(n):
            if i != k:
                multiplier = work[i][k]
                work[i] = [x-multiplier*y for x, y in zip(work[i], work[k])]
    return [row[-1] for row in work]


@dataclass
class Body:
    mass: float
    position: tuple
    velocity: tuple
    spin: tuple = (0, 0, 0)
    # A scalar isotropic inertia in this small reference (spherical inertia tensor).
    inertia: float = 1


def jacobian(bodies, contacts):
    rows = []
    for a, b, point, normal in contacts:
        if abs(dot(normal, normal)-1) > 1e-10:
            raise ValueError('contact normal must be a unit vector')
        row = [0.0]*(6*len(bodies))
        for index, sign in ((a, 1), (b, -1)):
            if index < 0:
                continue  # prescribed static support, outside the dynamic degrees of freedom
            torque = cross(sub(point, bodies[index].position), normal)
            row[6*index:6*index+6] = [sign*x for x in (*normal, *torque)]
        rows.append(row)
    return rows


def project(bodies, contacts):
    if not contacts or len(contacts) > 12 or any(b.mass <= 0 or b.inertia <= 0 for b in bodies):
        raise ValueError('reference requires positive masses/inertias and 1..12 contacts')
    rows = jacobian(bodies, contacts)
    inverse_mass = [x for b in bodies for x in (1/b.mass,)*3+(1/b.inertia,)*3]
    free = [x for b in bodies for x in (*b.velocity, *b.spin)]
    matrix = [[sum(x*y*m for x, y, m in zip(a, b, inverse_mass)) for b in rows] for a in rows]
    rhs = [dot(row, free) for row in rows]
    tolerance = 1e-10*max(1, *(abs(x) for x in rhs))
    for mask in range(1 << len(rows)):
        active = [i for i in range(len(rows)) if mask & (1 << i)]
        solution = linear_solve([[matrix[i][j] for j in active] for i in active], [-rhs[i] for i in active])
        if solution is None or any(x < -tolerance for x in solution):
            continue
        impulses = [0.0]*len(rows)
        for i, value in zip(active, solution):
            impulses[i] = max(0.0, value)
        residual = [b+dot(row, impulses) for b, row in zip(rhs, matrix)]
        if min(residual) < -tolerance or any(abs(residual[i]) > tolerance for i in active):
            continue
        final = [v+m*sum(row[k]*impulse for row, impulse in zip(rows, impulses))
                 for k, (v, m) in enumerate(zip(free, inverse_mass))]
        energy_before = 0.5*sum(v*v/m for v, m in zip(free, inverse_mass))
        energy_after = 0.5*sum(v*v/m for v, m in zip(final, inverse_mass))
        dissipated = 0.5*sum(x*dot(row, impulses) for x, row in zip(impulses, matrix))
        assert abs(energy_after-energy_before+dissipated) < 1e-9*max(1, energy_before)
        return final, impulses, {'energy_before_j': energy_before, 'energy_after_j': energy_after,
                                 'dissipated_j': dissipated, 'min_normal_velocity_m_s': min(residual)}
    raise RuntimeError('no complementary active set found; no result accepted')


def momenta(bodies, velocities):
    linear, angular = [0.0]*3, [0.0]*3
    for i, body in enumerate(bodies):
        v, w = velocities[6*i:6*i+3], velocities[6*i+3:6*i+6]
        p = [body.mass*x for x in v]
        orbital = cross(body.position, p)
        for k in range(3):
            linear[k] += p[k]
            angular[k] += orbital[k]+body.inertia*w[k]
    return linear, angular


def check_closed(bodies, final):
    before = [x for b in bodies for x in (*b.velocity, *b.spin)]
    for initial, result in zip(momenta(bodies, before), momenta(bodies, final)):
        assert max(abs(a-b) for a, b in zip(initial, result)) < 1e-10


def chain_example(duplicate=False):
    bodies = [Body(m, (i, 0, 0), (3 if i == 0 else 0, 0, 0)) for i, m in enumerate((1, 2, 3))]
    contacts = [(0, 1, (0.5, 0, 0), (-1, 0, 0)), (1, 2, (1.5, 0, 0), (-1, 0, 0))]
    if duplicate:
        contacts.append(contacts[0])
    final, impulses, report = project(bodies, contacts)
    assert max(abs(final[6*i]-0.5) for i in range(3)) < 1e-10
    assert abs(report['energy_after_j']-0.75) < 1e-10
    check_closed(bodies, final)
    if not duplicate:
        assert max(abs(x-y) for x, y in zip(impulses, (2.5, 1.5))) < 1e-10
    return report


def glancing_example(angle=0, shift=(0, 0, 0)):
    def rotate(v):
        c, s = math.cos(angle), math.sin(angle)
        return (c*v[0]-s*v[1], s*v[0]+c*v[1], v[2])
    def position(v):
        return tuple(x+y for x, y in zip(rotate(v), shift))
    bodies = [Body(1, position((-0.5, 0, 0)), rotate((1, 0, 0))),
              Body(1, position((0.5, 0, 0)), (0, 0, 0))]
    final, impulses, report = project(bodies, [(0, 1, position((0, 0.25, 0)), rotate((-1, 0, 0)))])
    assert abs(impulses[0]-1/2.125) < 1e-10
    expected = rotate((1-1/2.125, 0, 0))
    assert max(abs(x-y) for x, y in zip(final[:3], expected)) < 1e-10
    check_closed(bodies, final)
    return report


def gravity_example():
    step, gravity = 1/600, 9.81
    bodies = [Body(m, (0, i+0.5, 0), (0, -gravity*step, 0)) for i, m in enumerate((1, 2, 3))]
    contacts = [(0, -1, (0, 0, 0), (0, 1, 0)),
                (1, 0, (0, 1, 0), (0, 1, 0)), (2, 1, (0, 2, 0), (0, 1, 0))]
    final, impulses, report = project(bodies, contacts)
    assert max(abs(x) for x in final) < 1e-10
    assert max(abs(x-y) for x, y in zip(impulses, [m*gravity*step for m in (6, 5, 3)])) < 1e-10
    report['support_impulse_ns'] = impulses[0]
    report['gravity_impulse_ns'] = -sum(b.mass for b in bodies)*gravity*step
    assert abs(report['support_impulse_ns']+report['gravity_impulse_ns']) < 1e-10
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    results = {'scope': 'fixed poses, touching contacts, frictionless plastic impulses, isotropic inertias',
               'closed_chain': chain_example(), 'duplicate_constraint': chain_example(True),
               'glancing': glancing_example(), 'transformed_glancing': glancing_example(0.73, (2, -3, 1)),
               'gravity_supported_stack': gravity_example()}
    data = json.dumps(results, indent=2, allow_nan=False)+'\n'
    if args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            stream.write(data)
    print(data, end='')


if __name__ == '__main__':
    main()
