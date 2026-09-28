// The Neo-Hookean soft body: the particle lattice cut into tetrahedra, each one a small piece of
// an elastic solid, solved by XPBD. The chain from the physics to the code:
//
//  1. The material stores the energy density of the STABLE Neo-Hookean model (Smith, de Goes,
//     Kim 2018, "Stable Neo-Hookean Flesh Simulation", ACM TOG 37(2)) in the form of Macklin &
//     Mueller 2021 ("A Constraint-based Formulation of Stable Neo-Hookean Materials", MIG '21):
//         Psi(F) = mu/2 (tr(F^T F) - 3) + lambda/2 (det F - gamma)^2,   gamma = 1 + mu / lambda,
//     F the deformation gradient, mu and lambda the Lame parameters from E and nu (lambda shifted
//     by mu so that small strains follow linear elasticity: lameParameters, SoftBody.h). The shift
//     gamma makes the rest shape F = I free of stress; unlike the classic Neo-Hookean energy it
//     has no log(det F), so an inverted tetrahedron (det F < 0) is pushed back, not blown up.
//  2. Per tetrahedron of rest volume V, that energy is two constraints with compliances:
//         C_D = sqrt(tr(F^T F)),  alpha_D = 1 / (mu V)       (deviatoric: shape)
//         C_H = det F - gamma,    alpha_H = 1 / (lambda V)   (hydrostatic: volume)
//     since U = C_D^2 / (2 alpha_D) + C_H^2 / (2 alpha_H) = V Psi (up to a constant).
//  3. XPBD (Macklin, Mueller, Chentanez 2016, "XPBD: Position-Based Simulation of Compliant
//     Constrained Dynamics", eq. 18 and 26) turns every constraint into a position update with
//     alpha~ = alpha / h^2. The two constraints of a tetrahedron are solved TOGETHER (a 2 x 2
//     system): at rest their gradients are parallel and their pulls cancel exactly - solved one
//     after the other, the first would crush the tetrahedron and the second blow it up again.
//  4. Small steps (Macklin, Storey, Lu, Terdiman, Chentanez, Jeschke, Mueller 2019, "Small Steps
//     in Physics Simulation"): many short steps with one pass each converge far better than
//     long steps with many passes (ParticleSystem::stepSoftTetsInSmallSteps).
#include "particles/SoftBody.h"

#include <algorithm>

namespace rf {

Lame lameParameters(const SoftMaterial& m) {
    const float E = std::max(m.youngModulus, 1.0f);
    const float nu = std::clamp(m.poissonRatio, 0.05f, 0.495f);
    Lame l;
    l.mu = E / (2.0f * (1.0f + nu));
    l.lambda = E * nu / ((1.0f + nu) * (1.0f - 2.0f * nu)) + l.mu; // lambda_LE + mu (see SoftBody.h)
    return l;
}

float youngFromStiffness(float stiffness) {
    const float k = std::clamp(stiffness, 0.0f, 0.7f);
    return std::pow(10.0f, 3.5f + 5.0f * k); // [Pa]: 0.05 -> 5.6 kPa, 0.3 -> 100 kPa, 0.5 -> 1 MPa, 0.7 and up -> 10 MPa
}

float highestFrequency(const SoftMaterial& material, float spacing) {
    const Lame l = lameParameters(material); // lambda_LE + 2 mu = l.lambda + l.mu
    const float soundSpeed = std::sqrt((l.lambda + l.mu) / std::max(material.density, 1e-3f));
    return 2.0f * soundSpeed / std::max(spacing, 1e-6f);
}

// ---------------------------------------------------------------------------------------------
// Building: the lattice's cubes into tetrahedra
// ---------------------------------------------------------------------------------------------

// The five tetrahedra of one cube (the "5-tetrahedra" cut): the four corners of one parity -
// no two of them joined by an edge of the cube - make the big central tetrahedron (a third of
// the cube), and each of the four other corners with its three neighbours cuts off one corner (a
// sixth each). Corners are numbered q = a + 2b + 4c for the offset (a, b, c). The parity is the
// lattice's, (i + j + k + a + b + c) even: every face diagonal then joins two even points, and a
// face shared by two cubes is cut the same way from both sides - the tetrahedra fit face to face.
static int cubeTets(int i, int j, int k, int corners[5][4]) {
    const bool evenCube = (i + j + k) % 2 == 0;
    int centre[4], centres = 0, count = 0;
    for (int q = 0; q < 8; ++q) {
        const int bits = (q & 1) + ((q >> 1) & 1) + ((q >> 2) & 1);
        const bool evenPoint = (bits % 2 == 0) == evenCube;
        if (evenPoint) {
            centre[centres++] = q;
            continue;
        }
        corners[count][0] = q; // an odd corner and its three neighbours along the cube's edges
        corners[count][1] = q ^ 1;
        corners[count][2] = q ^ 2;
        corners[count][3] = q ^ 4;
        ++count;
    }
    for (int m = 0; m < 4; ++m) corners[count][m] = centre[m];
    return count + 1;
}

// The tetrahedra of the cube at lattice point (i, j, k) whose corners are all inside the body,
// appended to body.tets (corners as lattice node numbers). Returns how many.
static int addCubeTets(const SoftLattice& lattice, int i, int j, int k, SoftBody& body) {
    int node[8];
    for (int q = 0; q < 8; ++q) node[q] = lattice.nodeAt(i + (q & 1), j + ((q >> 1) & 1), k + ((q >> 2) & 1));
    int corners[5][4];
    const int tets = cubeTets(i, j, k, corners);
    int added = 0;
    for (int t = 0; t < tets; ++t) {
        SoftTet tet;
        bool inside = true;
        for (int m = 0; m < 4; ++m) {
            tet.v[m] = node[corners[t][m]];
            inside = inside && tet.v[m] >= 0;
        }
        if (!inside) continue;
        body.tets.push_back(tet);
        ++added;
    }
    return added;
}

// The cubes in eight colours by the parity of their lattice point, (i mod 2) + 2 (j mod 2) +
// 4 (k mod 2): two cubes of one colour are two or more points apart along some axis and share no
// corner, so the cubes of a colour can be solved at the same time without touching each other's
// particles (a graph colouring, as for cloth in Cloth.cpp).
void buildLatticeTets(SoftLattice& lattice, SoftBody& body) {
    body.tets.clear();
    body.cubeStart.clear();
    lattice.cubeIndex.assign(size_t(lattice.nx) * size_t(lattice.ny) * size_t(lattice.nz), -1);
    for (int colour = 0; colour < 8; ++colour) {
        body.colourStart[colour] = int(body.cubeStart.size());
        for (int k = (colour >> 2) & 1; k + 1 < lattice.nz; k += 2)
            for (int j = (colour >> 1) & 1; j + 1 < lattice.ny; j += 2)
                for (int i = colour & 1; i + 1 < lattice.nx; i += 2) {
                    const int first = int(body.tets.size());
                    if (addCubeTets(lattice, i, j, k, body) == 0) continue;
                    lattice.cubeIndex[size_t(lattice.point(i, j, k))] = int(body.cubeStart.size());
                    body.cubeStart.push_back(first);
                }
    }
    body.colourStart[8] = int(body.cubeStart.size());
    body.cubeStart.push_back(int(body.tets.size())); // the end of the last cube
}

void setRestShape(SoftBody& body, const std::vector<Vector3>& rest) {
    for (SoftTet& t : body.tets) {
        auto edges = [&]() {
            const Vector3 x0 = rest[size_t(t.v[0])];
            return Matrix3x3::fromColumns(rest[size_t(t.v[1])] - x0, rest[size_t(t.v[2])] - x0, rest[size_t(t.v[3])] - x0);
        };
        if (edges().determinant() < 0) std::swap(t.v[2], t.v[3]); // positive orientation: det F = +1 at rest
        const Matrix3x3 Dm = edges();
        t.restVolume = Dm.determinant() / 6.0f;
        t.restInverse = Dm.inverse();
        t.lambdaD = t.lambdaH = 0;
    }
}

// Barycentric coordinates of x in tetrahedron t at rest: (b1, b2, b3) = D_m^-1 (x - x0), b0 = 1 -
// b1 - b2 - b3. All four >= 0 inside; the smallest of them is how far inside (negative: outside).
static void restBarycentric(const SoftTet& t, const std::vector<Vector3>& rest, const Vector3& x, float w[4]) {
    const Vector3 b = t.restInverse * (x - rest[size_t(t.v[0])]);
    w[0] = 1.0f - b.x - b.y - b.z;
    w[1] = b.x;
    w[2] = b.y;
    w[3] = b.z;
}

// The tetrahedron that holds x best among those of the cubes within `reach` lattice points of it:
// the one whose smallest barycentric coordinate is largest (inside it, or the least outside).
static bool bestTetNear(const SoftBody& body, const SoftLattice& lattice, const std::vector<Vector3>& rest, const Vector3& x,
                        int reach, SkinBinding& out) {
    const Vector3 cell = (x - lattice.origin) * (1.0f / lattice.spacing);
    const int ci = int(std::floor(cell.x)), cj = int(std::floor(cell.y)), ck = int(std::floor(cell.z));
    float best = -kInf;
    for (int k = ck - reach; k <= ck + reach; ++k)
        for (int j = cj - reach; j <= cj + reach; ++j)
            for (int i = ci - reach; i <= ci + reach; ++i) {
                if (i < 0 || j < 0 || k < 0 || i >= lattice.nx || j >= lattice.ny || k >= lattice.nz) continue;
                const int cube = lattice.cubeIndex[size_t(lattice.point(i, j, k))];
                if (cube < 0) continue;
                for (int t = body.cubeStart[size_t(cube)]; t < body.cubeStart[size_t(cube) + 1]; ++t) {
                    float w[4];
                    restBarycentric(body.tets[size_t(t)], rest, x, w);
                    const float inside = std::min(std::min(w[0], w[1]), std::min(w[2], w[3]));
                    if (inside <= best) continue;
                    best = inside;
                    out.tet = t;
                    std::copy(w, w + 4, out.weight);
                }
            }
    return best > -kInf;
}

void bindSurfaceToTets(SoftBody& body, const SoftLattice& lattice, const TriMesh& restSurface, const std::vector<Vector3>& rest) {
    // A mesh vertex lies at most about a spacing from the particles' hull, so the cubes one point
    // around it hold its tetrahedron; a far spike of the mesh looks further out.
    body.surface = restSurface;
    body.skin.assign(restSurface.positions.size(), SkinBinding());
    for (size_t v = 0; v < body.skin.size(); ++v) {
        const Vector3& x = restSurface.positions[v];
        for (int reach = 1; reach < 64; reach *= 2)
            if (bestTetNear(body, lattice, rest, x, reach, body.skin[v])) break;
    }
}

// ---------------------------------------------------------------------------------------------
// Solving: two XPBD constraints per tetrahedron
// ---------------------------------------------------------------------------------------------

// The deformation gradient F = D_s D_m^-1 of a tetrahedron now (D_s its edge matrix now).
static Matrix3x3 deformation(const SoftTet& t, const std::vector<Vector3>& p) {
    const Vector3 x0 = p[size_t(t.v[0])];
    const Matrix3x3 Ds = Matrix3x3::fromColumns(p[size_t(t.v[1])] - x0, p[size_t(t.v[2])] - x0, p[size_t(t.v[3])] - x0);
    return Ds * t.restInverse;
}

// cof(F) = det(F) F^-T = [f1 x f2, f2 x f0, f0 x f1] (columns f of F): the gradient of det F.
static Matrix3x3 cofactor(const Matrix3x3& F) {
    const Vector3 f0 = F.col(0), f1 = F.col(1), f2 = F.col(2);
    return Matrix3x3::fromColumns(cross(f1, f2), cross(f2, f0), cross(f0, f1));
}

// The two constraints of a tetrahedron and their gradients by its four corners.
struct TetConstraints {
    double deviatoric = 0, hydrostatic = 0; // C_D = sqrt(tr F^T F), C_H = det F - gamma
    Vector3 gradD[4], gradH[4];
};

// C_D, C_H and their gradients. By the chain rule through F = D_s D_m^-1: the gradient by the
// corners x1, x2, x3 are the columns of (dC/dF) D_m^-T and by x0 minus their sum, with dC_D/dF =
// F / C_D and dC_H/dF = cof(F). False for a tetrahedron crushed to a point (no way to push).
static bool evaluateTet(const SoftTet& t, const std::vector<Vector3>& p, float gamma, TetConstraints& c) {
    const Matrix3x3 F = deformation(t, p);
    const float ic = length2(F.col(0)) + length2(F.col(1)) + length2(F.col(2)); // tr(F^T F)
    if (ic < 1e-12f) return false;
    c.deviatoric = std::sqrt(double(ic));
    c.hydrostatic = double(F.determinant()) - gamma;
    const Matrix3x3 DmT = t.restInverse.transposed();
    const Matrix3x3 GD = F * float(1.0 / c.deviatoric) * DmT, GH = cofactor(F) * DmT;
    for (int k = 0; k < 3; ++k) {
        c.gradD[k + 1] = GD.col(k);
        c.gradH[k + 1] = GH.col(k);
    }
    c.gradD[0] = -(c.gradD[1] + c.gradD[2] + c.gradD[3]);
    c.gradH[0] = -(c.gradH[1] + c.gradH[2] + c.gradH[3]);
    return true;
}

// What one step of length h asks of the constraints: their compliances alpha~ = 1 / (stiffness
// V h^2) per unit volume (divided by V per tetrahedron), gamma, and the damping factor.
struct TetStep {
    double complianceD = 0, complianceH = 0; // 1 / (mu h^2), 1 / (lambda h^2)
    float gamma = 1;
    double damping = 0; // XPBD's gamma = beta / h for Rayleigh damping beta K (Macklin et al. 2016, eq. 26)
};

// One tetrahedron by XPBD, both constraints at once. The linearized system in the multipliers
// (Macklin et al. 2016, eq. 18 with damping, eq. 26, one row per constraint, k, l in {D, H}):
//     (1 + g) sum_corners w grad C_k . grad C_l  dlambda_l + alpha~_k dlambda_k
//         = -C_k - alpha~_k lambda_k - g grad C_k . (x - x_start)
// solved exactly (2 x 2), then every corner moves by w (grad C_D dlambda_D + grad C_H dlambda_H).
static void solveTet(SoftTet& t, const TetStep& s, std::vector<Vector3>& p, const std::vector<float>& w,
                     const std::vector<Vector3>* start) {
    TetConstraints c;
    if (!evaluateTet(t, p, s.gamma, c)) return;
    double dd = 0, hh = 0, dh = 0, rateD = 0, rateH = 0;
    for (int k = 0; k < 4; ++k) {
        const double wk = w[size_t(t.v[k])];
        dd += wk * double(dot(c.gradD[k], c.gradD[k]));
        hh += wk * double(dot(c.gradH[k], c.gradH[k]));
        dh += wk * double(dot(c.gradD[k], c.gradH[k]));
        if (!start) continue;
        const Vector3 moved = p[size_t(t.v[k])] - (*start)[size_t(t.v[k])];
        rateD += double(dot(c.gradD[k], moved));
        rateH += double(dot(c.gradH[k], moved));
    }
    if (dd + hh == 0) return; // all four corners held
    const double g = start ? s.damping : 0.0, alphaD = s.complianceD / t.restVolume, alphaH = s.complianceH / t.restVolume;
    const double a11 = (1 + g) * dd + alphaD, a22 = (1 + g) * hh + alphaH, a12 = (1 + g) * dh;
    const double b1 = -(c.deviatoric + alphaD * t.lambdaD + g * rateD);
    const double b2 = -(c.hydrostatic + alphaH * t.lambdaH + g * rateH);
    const double det = a11 * a22 - a12 * a12;
    if (!(det > 0)) return;
    const double dD = (b1 * a22 - b2 * a12) / det, dH = (a11 * b2 - a12 * b1) / det;
    t.lambdaD += float(dD);
    t.lambdaH += float(dH);
    for (int k = 0; k < 4; ++k) {
        const float wk = w[size_t(t.v[k])];
        if (wk > 0) p[size_t(t.v[k])] += (c.gradD[k] * float(dD) + c.gradH[k] * float(dH)) * wk;
    }
}

void solveCubeTets(SoftBody& body, int cube, float h, std::vector<Vector3>& p, const std::vector<float>& invMass,
                   const std::vector<Vector3>* stepStart) {
    const Lame l = lameParameters(body.material);
    TetStep s;
    s.complianceD = 1.0 / (double(l.mu) * h * h);
    s.complianceH = 1.0 / (double(l.lambda) * h * h);
    s.gamma = 1.0f + l.mu / l.lambda;
    s.damping = double(body.material.damping) / h;
    for (int t = body.cubeStart[size_t(cube)]; t < body.cubeStart[size_t(cube) + 1]; ++t)
        solveTet(body.tets[size_t(t)], s, p, invMass, stepStart);
}

double elasticEnergy(const SoftBody& body, const std::vector<Vector3>& x) {
    // V Psi summed, Psi = mu/2 (tr F^T F - 3) + lambda/2 ((det F - gamma)^2 - (1 - gamma)^2): the
    // last term is the rest value of the second, so the rest shape stores exactly zero.
    const Lame l = lameParameters(body.material);
    const double gamma = 1.0 + double(l.mu) / l.lambda;
    double energy = 0;
    for (const SoftTet& t : body.tets) {
        const Matrix3x3 F = deformation(t, x);
        const double ic = double(length2(F.col(0))) + double(length2(F.col(1))) + double(length2(F.col(2)));
        const double J = F.determinant(), volume = (J - gamma) * (J - gamma) - (1.0 - gamma) * (1.0 - gamma);
        energy += t.restVolume * (0.5 * l.mu * (ic - 3.0) + 0.5 * l.lambda * volume);
    }
    return energy;
}

void turnNormalsWithTets(const SoftBody& body, const std::vector<Vector3>& x, const std::vector<Vector3>& restNormal,
                         std::vector<Vector3>& normal) {
    for (int i : body.particles) normal[size_t(i)] = Vector3(0.0f);
    for (const SoftTet& t : body.tets) {
        const Matrix3x3 C = cofactor(deformation(t, x));
        for (int k = 0; k < 4; ++k) normal[size_t(t.v[k])] += C * restNormal[size_t(t.v[k])];
    }
    for (int i : body.particles) {
        Vector3& n = normal[size_t(i)];
        n = length2(n) > 1e-12f ? normalize(n) : restNormal[size_t(i)];
    }
}

} // namespace rf
