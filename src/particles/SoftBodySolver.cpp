// The material of the soft bodies: corotated linear elasticity - the energy of Unreal's Chaos Flesh
// and PhysX's deformable volumes - solved as XPBD constraints, one pass per small step. The step
// itself, with the contacts, is in SoftBodyStep.cpp; the tetrahedra and nodes are built in SoftBody.cpp.
//
// The energy of a body, with Lamé's mu and lambda (from Young's E and Poisson's nu):
//   U = sum_t V_t mu |F_t - R_t|^2 + sum_i V_i lambda/2 (J_i - 1)^2,
// F_t the deformation gradient of tetrahedron t, R_t its rotation (the polar decomposition, found
// by Mueller et al. 2016 from the last pass's), J_i the volume ratio at particle i: the mean of det F
// over the tetrahedra around it (Bonet & Burton 1998). For small strain eps it is
// mu eps:eps + lambda/2 (tr eps)^2 - linear elasticity, so a beam bends and a bar stretches by E
// and nu. Both terms are zero at rest and under any rotation.
//
// As XPBD constraints (Macklin, Muller, Chentanez 2016), each of energy C^2 / (2 alpha):
//   shape    C_t = F_t - R_t (nine components),   alpha = 1 / (2 mu V_t)
//   volume   C_i = J_i - 1,                       alpha = 1 / (lambda V_i)
// Three things make one pass per small step settle exactly in the energy's equilibrium - the
// static answer the same at every step (0.02 % apart between Courant numbers 1 and 1/4, the
// principle of virtual work met to 0.2 % of the weight's work), a bar at nu = 0.49 as long and as
// thin as theory says:
//  1. the shape constraint as a matrix, solved exactly per tetrahedron (solveTet). As the scalar
//     |F - R| its gradient - the direction of F - R - turned round within a pass when a loaded body
//     sat near its rest shape, and the passes settled where the forces they had applied balanced the
//     weight, not the forces of the shape they left: a clamped beam at rest had its material's
//     forces 46 times its weight out of balance;
//  2. the volume per particle, not per tetrahedron (solveNode). With one constraint per
//     tetrahedron there are more of them than the particles have degrees of freedom: hundreds of
//     patterns of pressure push on no particle, nothing damps them, and carried from step to step
//     they grew until a nearly incompressible bar blew up. The tetrahedra also locked, stiffer the
//     nearer nu came to 1/2;
//  3. the multipliers carried from small step to small step and applied first (warm start). From
//     zero one pass leaves every constraint short, and the body came out 25-40 % too soft.
//
// Positions come in two parts, base + u: the particles where the step began and how far this step
// has moved them. A small step of a stiff body moves a particle by less than the float step of its
// position (6e-8 m at 1 m from the origin): added to the position, the moves were rounded away, and
// the finer the step, the softer a rubber beam came out (+61 % at 53 small steps). The edges of a
// tetrahedron are summed from the edges of the base and of u, both exact.
//
// Tried and dropped: the stable Neo-Hookean split of Macklin & Muller 2021 (its two constraints
// carry a prestress each; stiff against the step their 2 x 2 system blew a resting cube to 900 %
// strain); Vertex Block Descent (Chen et al. 2024: sluggish short of convergence - a rubber
// cantilever swung at half its frequency - and 2.4 times the cost of a pass).
#include "particles/SoftBody.h"

#include "core/Parallel.h"

#include <algorithm>
#include <cmath>

namespace rf {

namespace {

using PassConstants = SoftPassConstants;

PassConstants passConstants(const SoftMaterial& m, float dt) {
    PassConstants k;
    k.shapeCompliance = 1.0f / (2.0f * shearModulus(m));
    const float lambda = lameLambda(m);
    k.volumeCompliance = lambda > 0 ? 1.0f / lambda : 0.0f;
    k.invDt2 = 1.0f / (dt * dt);
    return k;
}

// The corners' shares of F: F = sum_c x_c b_c^T, b_1..3 the rows of D_m^-1, b_0 = -(b_1 + b_2 + b_3).
// A constraint with dC/dF = Q has the gradient Q b_c at corner c.
void shares(const SoftTet& t, Vector3 b[4]) {
    const Matrix3x3 T = t.restInverse.transposed();
    b[1] = T.col(0), b[2] = T.col(1), b[3] = T.col(2), b[0] = -(b[1] + b[2] + b[3]);
}

// A tetrahedron's shape constraint C = F - R with its 3 x 3 multiplier Lambda (XPBD, eq. 18 of
// Macklin et al. 2016, for a vector constraint). With R held, C is linear in the corners
// (dC/dx_c = b_c), and one update solves the tetrahedron exactly:
//   dLambda = -(C + a Lambda) (A + a I)^-1,   A = sum_c w_c b_c b_c^T,   dx_c = w_c dLambda b_c,
// a = alpha / dt^2. Lambda stays in world axes: kept in the tetrahedron's turned frame, it swung
// with the rotation of a tetrahedron turned inside out - which jumps - and the warm start pushed the
// wrong way.
void solveTet(SoftTet& t, const PassConstants& k, const std::vector<Vector3>& base, std::vector<Vector3>& u, const std::vector<float>& invMass) {
    const float w[4] = {invMass[size_t(t.v[0])], invMass[size_t(t.v[1])], invMass[size_t(t.v[2])], invMass[size_t(t.v[3])]};
    if (w[0] + w[1] + w[2] + w[3] == 0) return;
    Vector3 b[4];
    shares(t, b);
    const Matrix3x3 F = deformationGradient(t, base, u);
    t.rotation = extractRotation(F, t.rotation, 3);
    const Matrix3x3 R = t.rotation.toMatrix3x3();
    const float a = k.shapeCompliance / t.restVolume * k.invDt2;
    Matrix3x3 A = Matrix3x3::identity() * a;
    for (int c = 0; c < 4; ++c) A += Matrix3x3::outer(b[c], b[c]) * w[c];
    const Matrix3x3 dL = (F - R + t.lambdaShape * a) * A.inverse(1e-30f) * -1.0f;
    for (int c = 0; c < 4; ++c) u[size_t(t.v[c])] += dL * b[c] * w[c];
    t.lambdaShape += dL;
}

// A node's volume constraint C = J_i - 1 and its gradient at the particles around it (at most 27 on
// the lattice; a particle in several of the tetrahedra takes the sum of their shares), grad[q] at
// nodeNear[nearBegin + q]. With s_t = V_t / (4 V_i), J_i = sum_t s_t det F_t = sum_t v_t / (4 V_i):
// V_t det F_t is simply the tetrahedron's volume now, v_t = e1 . (e2 x e3) / 6 of its edges from
// corner 0, whose gradient at corner 1 is e2 x e3 / 6 (and so on round; corner 0 takes minus the
// sum) - three cross products per tetrahedron where cof(F) b_c took three times the work.
constexpr int kNodeParticles = 27;

float nodeConstraint(const SoftBody& body, const SoftNode& n, const std::vector<Vector3>& base, const std::vector<Vector3>& u,
                     Vector3 grad[kNodeParticles]) {
    for (int q = 0; q < n.nearEnd - n.nearBegin; ++q) grad[q] = Vector3(0.0f);
    const float per = 1.0f / (24.0f * n.restVolume); // (1 / 6) / (4 V_i)
    double sixV = 0;
    for (int s = n.starBegin; s < n.starEnd; ++s) {
        const SoftTet& t = body.tets[size_t(body.nodeStar[size_t(s)])];
        auto edge = [&](int c) {
            const size_t a = size_t(t.v[size_t(c)]), o = size_t(t.v[0]);
            return (base[a] - base[o]) + (u[a] - u[o]);
        };
        const Vector3 e1 = edge(1), e2 = edge(2), e3 = edge(3);
        const Vector3 g1 = cross(e2, e3), g2 = cross(e3, e1), g3 = cross(e1, e2);
        sixV += double(dot(e1, g1));
        const uint8_t* slot = &body.nodeSlot[size_t(4 * s)];
        grad[slot[0]] -= (g1 + g2 + g3) * per;
        grad[slot[1]] += g1 * per;
        grad[slot[2]] += g2 * per;
        grad[slot[3]] += g3 * per;
    }
    return float(sixV * double(per) - 1.0);
}

// The node's scalar XPBD update: dl = (-C - a l) / (sum w |grad C|^2 + a). The gradient is kept for
// the next warm start (warmStartNode).
void solveNode(SoftBody& body, SoftNode& n, const PassConstants& k, const std::vector<Vector3>& base, std::vector<Vector3>& u,
               const std::vector<float>& invMass) {
    Vector3 grad[kNodeParticles];
    const float C = nodeConstraint(body, n, base, u, grad);
    const int count = n.nearEnd - n.nearBegin;
    const int* near = &body.nodeNear[size_t(n.nearBegin)];
    float wg = 0;
    for (int q = 0; q < count; ++q) wg += invMass[size_t(near[q])] * length2(grad[q]);
    const float a = k.volumeCompliance / n.restVolume * k.invDt2;
    if (wg + a < 1e-30f) return;
    const float dl = (-C - a * n.lambda) / (wg + a);
    n.lambda += dl;
    Vector3* kept = &body.nodeGradient[size_t(n.nearBegin)];
    for (int q = 0; q < count; ++q) {
        u[size_t(near[q])] += grad[q] * (invMass[size_t(near[q])] * dl);
        kept[q] = grad[q];
    }
}

// The warm start's moves - the forces of the last small step, x = x~ + M^-1 grad C^T lambda, the
// start a Gauss-Seidel pass needs to go on from those multipliers - into `move` (slot q = particle
// first + q), not into the positions: see warmStart. A node pushes along the gradient its last
// update found (at rest it is the gradient now, and the static answer stays exact).
void warmStartTet(const SoftTet& t, const std::vector<float>& invMass, std::vector<Vector3>& move, int first) {
    Vector3 b[4];
    shares(t, b);
    for (int c = 0; c < 4; ++c) move[size_t(t.v[c] - first)] += t.lambdaShape * b[c] * invMass[size_t(t.v[c])];
}

void warmStartNode(const SoftBody& body, const SoftNode& n, const std::vector<float>& invMass, std::vector<Vector3>& move, int first) {
    if (n.lambda == 0) return;
    const int* near = &body.nodeNear[size_t(n.nearBegin)];
    const Vector3* kept = &body.nodeGradient[size_t(n.nearBegin)];
    for (int q = 0; q < n.nearEnd - n.nearBegin; ++q) move[size_t(near[q] - first)] += kept[q] * (invMass[size_t(near[q])] * n.lambda);
}

// Runs f(body, i) for every index i of every colour of every body, colour after colour: one colour
// of all the bodies in one parallel loop - bodies share no particle, and neither do two indices of
// one colour of a body - so the result is the same on any number of threads, and a big body and
// many small ones fill the threads alike. colourStart(body) gives a body's colour table.
template <class ColourStart, class Visit>
void byColour(const std::vector<SoftBody*>& bodies, std::vector<SoftColourRun>& runs, ColourStart colourStart, Visit f, int grain = 64) {
    size_t colours = 0;
    for (const SoftBody* b : bodies) colours = std::max(colours, colourStart(*b).size());
    for (size_t c = 0; c + 1 < colours; ++c) {
        runs.clear();
        int total = 0;
        for (int k = 0; k < int(bodies.size()); ++k) {
            const std::vector<int>& start = colourStart(*bodies[size_t(k)]);
            if (c + 1 >= start.size() || start[c + 1] == start[c]) continue;
            runs.push_back({k, start[c], total});
            total += start[c + 1] - start[c];
        }
        parallelRanges(total, [&](int from, int to) {
            size_t r = size_t(std::upper_bound(runs.begin(), runs.end(), from, [](int i, const SoftColourRun& run) { return i < run.offset; }) -
                              runs.begin()) - 1;
            for (int i = from; i < to; ++i) {
                while (r + 1 < runs.size() && i >= runs[r + 1].offset) ++r;
                f(*bodies[size_t(runs[r].body)], runs[r].begin + (i - runs[r].offset));
            }
        }, grain);
    }
}

// The warm start of the bodies: every tetrahedron's and node's move found from the same
// configuration (Jacobi), then added. Applied one after another, each would meet its neighbours
// already moved by forces far larger than the strain they hold (the pressure wave of a nearly
// incompressible body crosses several cells a small step), and would push along the gradients of a
// shape that is not the one the pass starts from.
void warmStart(const std::vector<SoftBody*>& bodies, std::vector<Vector3>& u, const std::vector<float>& invMass,
               std::vector<SoftColourRun>& runs) {
    byColour(bodies, runs, [](const SoftBody& b) -> const std::vector<int>& { return b.colourStart; },
             [&](SoftBody& b, int t) { warmStartTet(b.tets[size_t(t)], invMass, b.warmMove, b.particles.front()); });
    byColour(bodies, runs, [](const SoftBody& b) -> const std::vector<int>& { return b.nodeColourStart; }, [&](SoftBody& b, int q) {
        if (b.pass.volumeCompliance > 0) warmStartNode(b, b.nodes[size_t(b.nodeOrder[size_t(q)])], invMass, b.warmMove, b.particles.front());
    });
    parallelFor(int(bodies.size()), [&](int k) {
        const SoftBody& b = *bodies[size_t(k)];
        for (size_t q = 0; q < b.particles.size(); ++q) u[size_t(b.particles.front()) + q] += b.warmMove[q];
    }, 1);
}

// The material's forces are internal: they must not turn the body. Each projection keeps the
// angular momentum to first order in its correction, but a corrections-squared remainder stays,
// and in a loaded body it has the same sign pass after pass: a jelly resting on the floor turned
// about the vertical by 10 degrees in 3 s. So after a pass the body is turned back, about its
// centre of mass, by the turn the pass gave it - w = I^-1 sum m r x dx - applied as an exact
// rotation (the shape is untouched). Not while a particle is held: the pin's torque is an outside one.
void undoTurn(const SoftBody& body, const std::vector<Vector3>& base, std::vector<Vector3>& u, const std::vector<float>& invMass,
              const std::vector<Vector3>& before) {
    Vector3 centreBase(0.0f), centreMove(0.0f); // the centre of mass, in the same two parts
    float mass = 0;
    for (int i : body.particles) {
        const float w = invMass[size_t(i)];
        if (w == 0) return;
        centreBase += base[size_t(i)] / w;
        centreMove += u[size_t(i)] / w;
        mass += 1.0f / w;
    }
    centreBase /= mass;
    centreMove /= mass;
    auto arm = [&](int i) { return (base[size_t(i)] - centreBase) + (u[size_t(i)] - centreMove); };
    Matrix3x3 inertia = Matrix3x3::zero();
    Vector3 turned(0.0f);
    const int first = body.particles.front();
    for (int i : body.particles) {
        const float m = 1.0f / invMass[size_t(i)];
        const Vector3 r = arm(i);
        inertia += (Matrix3x3::identity() * length2(r) - Matrix3x3::outer(r, r)) * m;
        turned += cross(r, u[size_t(i)] - before[size_t(i - first)]) * m;
    }
    const Vector3 w = inertia.inverse(1e-30f) * turned;
    const float angle = length(w);
    if (angle < 1e-12f) return;
    const Matrix3x3 turnBack = Quaternion::fromAxisAngle(w / angle, -angle).toMatrix3x3() - Matrix3x3::identity();
    for (int i : body.particles) u[size_t(i)] += turnBack * arm(i);
}

} // namespace

void solveSoftBodies(const std::vector<SoftBody*>& bodies, const std::vector<Vector3>& base, std::vector<Vector3>& u,
                     const std::vector<float>& invMass, float dt, bool withWarmStart, std::vector<SoftColourRun>& runs) {
    parallelFor(int(bodies.size()), [&](int k) {
        SoftBody& b = *bodies[size_t(k)];
        b.pass = passConstants(b.material, dt);
        b.passStart.assign(u.begin() + b.particles.front(), u.begin() + b.particles.front() + int(b.particles.size()));
        if (withWarmStart) b.warmMove.assign(b.particles.size(), Vector3(0.0f));
    }, 1);
    if (withWarmStart) warmStart(bodies, u, invMass, runs);
    byColour(bodies, runs, [](const SoftBody& b) -> const std::vector<int>& { return b.colourStart; },
             [&](SoftBody& b, int t) { solveTet(b.tets[size_t(t)], b.pass, base, u, invMass); });
    byColour(bodies, runs, [](const SoftBody& b) -> const std::vector<int>& { return b.nodeColourStart; }, [&](SoftBody& b, int q) {
        if (b.pass.volumeCompliance > 0) solveNode(b, b.nodes[size_t(b.nodeOrder[size_t(q)])], b.pass, base, u, invMass);
    }, 8); // a node is heavy: small runs keep every thread busy
    parallelFor(int(bodies.size()), [&](int k) { undoTurn(*bodies[size_t(k)], base, u, invMass, bodies[size_t(k)]->passStart); }, 1);
}

void scaleSoftMultipliers(SoftBody& body, float f) {
    for (SoftTet& t : body.tets) t.lambdaShape = t.lambdaShape * f;
    for (SoftNode& n : body.nodes) n.lambda *= f;
}

// Internal friction (Müller, Heidelberger, Hennix, Ratcliff 2007, "Position Based Dynamics", sec.
// 3.5): at the end of every substep the part of each particle's velocity that is not the motion of
// the body as a whole - its centre's velocity plus its turning, w x r - decays by exp(-rate dt).
// The body's momentum and angular momentum stay exactly (the removed parts sum to zero in both);
// a thrown jelly flies on and a spinning one keeps spinning, only its wobble dies. It acts on the
// velocities after the step, so the material is as stiff as before: tried inside the position
// solve, as a Rayleigh term on each constraint's rate, it undid part of the elastic correction of
// the same step, and a rubber beam swung wider and slower than without it.
void dampSoftBody(const SoftBody& b, std::vector<Vector3>& v, const std::vector<Vector3>& x, const std::vector<float>& invMass, float dt) {
    if (b.material.damping <= 0) return;
    const float keep = std::exp(-b.material.damping * dt);
    // A body held by pins or by the mouse cannot move as a whole: its wobble is its motion against
    // what holds it (a clamped beam swings as if it turned about the clamp - taken for a turn of
    // the whole body, its swing was never damped).
    Vector3 held(0.0f);
    int holders = 0;
    for (int i : b.particles)
        if (invMass[size_t(i)] == 0) held += v[size_t(i)], ++holders;
    if (holders > 0) {
        held /= float(holders);
        for (int i : b.particles)
            if (invMass[size_t(i)] > 0) v[size_t(i)] = held + (v[size_t(i)] - held) * keep;
        return;
    }
    Vector3 centre(0.0f), momentum(0.0f);
    float mass = 0;
    for (int i : b.particles) {
        centre += x[size_t(i)] / invMass[size_t(i)];
        momentum += v[size_t(i)] / invMass[size_t(i)];
        mass += 1.0f / invMass[size_t(i)];
    }
    centre /= mass;
    const Vector3 velocity = momentum / mass;
    Matrix3x3 inertia = Matrix3x3::zero();
    Vector3 L(0.0f);
    for (int i : b.particles) {
        const float w = invMass[size_t(i)];
        const Vector3 r = x[size_t(i)] - centre;
        inertia += (Matrix3x3::identity() * length2(r) - Matrix3x3::outer(r, r)) * (1.0f / w);
        L += cross(r, v[size_t(i)] - velocity) / w;
    }
    const Vector3 omega = inertia.inverse(1e-30f) * L;
    for (int i : b.particles) {
        const Vector3 rigid = velocity + cross(omega, x[size_t(i)] - centre);
        v[size_t(i)] = rigid + (v[size_t(i)] - rigid) * keep;
    }
}

} // namespace rf
