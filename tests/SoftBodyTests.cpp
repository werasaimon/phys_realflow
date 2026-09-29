// Soft bodies against the standard problems of solid mechanics, as elastic solvers are checked
// (Smith, de Goes, Kim 2018; Macklin & Müller 2021; Chen et al. 2024), and against the principle
// of virtual work: statics against beam theory on three grids, the energy minimum of the rest
// state, the first bending frequency, a bar hanging under its own weight (extension, Poisson's
// contraction, no volume locking at nu 0.49), momentum and angular momentum of a spinning body,
// recovery from being dragged inside out, rest without jitter, independence of the step, and
// energy that never grows.
//
// The theory is evaluated for the body the solver's material fills: the tetrahedra span the
// particle centres, half a spacing in from every face, so a section n particles across is an
// elastic core of n - 1 spacings carrying the weight of n (the half-cell rim is mass without
// stiffness). The continuum answer is printed beside it; the two meet as n grows.
#include "TestRunner.h"
#include "Tests.h"

#include <chrono>
#include <cmath>
#include <memory>

namespace {

constexpr float kDt = 1.0f / 180.0f; // the particle substep of a 60 Hz frame (Simulation: 3 per frame)
constexpr float kBeamY = 1.0f;        // the height the beams are clamped at

// A particle system with particles of radius r in a large box and gravity g, and a soft box body
// of the given size centred at c.
struct Rig {
    ParticleSystem s;
    int body = -1;
    Rig(float r, const Vector3& c, const Vector3& size, const SoftMaterial& m, const Vector3& g, float courant = 0) {
        s.params.particleRadius = r;
        s.params.gravity = g;
        if (courant > 0) s.params.softCourant = courant;
        s.reset(AABB({-3, -3, -3}, {3, 3, 3}));
        TriMesh box = primitives::box(size * 0.5f);
        box.translate(c);
        body = s.addSoftBody(box, m, Vector3(1));
    }
    const SoftBody& b() const { return s.softBodies()[size_t(body)]; }
    const Vector3& x(int i) const { return s.positions()[size_t(i)]; }
    float mass() const {
        float m = 0;
        for (int i : b().particles) m += s.invMasses()[size_t(i)] > 0 ? 1.0f / s.invMasses()[size_t(i)] : 0.0f;
        return m;
    }
    float meanY(const std::vector<int>& ids) const {
        float y = 0;
        for (int i : ids) y += x(i).y;
        return y / float(ids.size());
    }
};

// The corotated energy the solver minimises (SoftBodySolver.cpp), in double: sum V [mu |F - R|^2 +
// lambda/2 (det F - 1)^2], R the rotation of F.
double elasticEnergy(const SoftBody& b, const std::vector<Vector3>& x) {
    const double mu = shearModulus(b.material), lambda = lameLambda(b.material);
    double U = 0;
    for (const SoftTet& t : b.tets) {
        const Matrix3x3 F = deformationGradient(t, x);
        const Matrix3x3 D = F - extractRotation(F, t.rotation, 30).toMatrix3x3();
        double n2 = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) n2 += double(D.m[i][j]) * D.m[i][j];
        U += t.restVolume * mu * n2;
    }
    for (const SoftNode& n : b.nodes) { // the volume per particle: lambda/2 V_i (J_i - 1)^2
        double J = 0;
        for (int s = n.starBegin; s < n.starEnd; ++s) {
            const SoftTet& t = b.tets[size_t(b.nodeStar[size_t(s)])];
            J += 0.25 * t.restVolume * deformationGradient(t, x).determinant() / n.restVolume;
        }
        U += 0.5 * lambda * n.restVolume * (J - 1) * (J - 1);
    }
    return U;
}

// The weight's potential, sum m g.x over the free particles.
double gravityPotential(const ParticleSystem& s, const SoftBody& b, const std::vector<Vector3>& x) {
    double V = 0;
    for (int i : b.particles)
        if (s.invMasses()[size_t(i)] > 0) V -= dot(s.params.gravity, x[size_t(i)]) / s.invMasses()[size_t(i)];
    return V;
}

float kineticEnergy(const ParticleSystem& s, const SoftBody& b) {
    double K = 0;
    for (int i : b.particles)
        if (s.invMasses()[size_t(i)] > 0) K += 0.5 * length2(s.velocities()[size_t(i)]) / s.invMasses()[size_t(i)];
    return float(K);
}

// The largest distance of a particle from where the best rigid motion of the rest shape puts it,
// over the body's size: 0 = the rest shape exactly (turned and moved anywhere).
float shapeError(const ParticleSystem& s, const SoftBody& b, const std::vector<Vector3>& rest) {
    Vector3 c(0.0f), c0(0.0f);
    for (size_t k = 0; k < rest.size(); ++k) c += s.positions()[size_t(b.particles[k])], c0 += rest[k];
    c /= float(rest.size());
    c0 /= float(rest.size());
    Matrix3x3 A = Matrix3x3::zero();
    float size = 0;
    for (size_t k = 0; k < rest.size(); ++k) {
        A += Matrix3x3::outer(s.positions()[size_t(b.particles[k])] - c, rest[k] - c0);
        size = std::max(size, length(rest[k] - c0));
    }
    const Matrix3x3 R = extractRotation(A, Quaternion(), 50).toMatrix3x3();
    float worst = 0;
    for (size_t k = 0; k < rest.size(); ++k) worst = std::max(worst, length(s.positions()[size_t(b.particles[k])] - c - R * (rest[k] - c0)));
    return worst / size;
}

bool allFinite(const ParticleSystem& s) {
    for (const std::vector<Vector3>* a : {&s.positions(), &s.velocities()})
        for (const Vector3& p : *a)
            if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) return false;
    return true;
}

// A cantilever of square section h x h and free length L along x, n particles across, its first
// two layers of particles clamped. Stepped for `seconds`; the tip's drop is recorded every step.
struct Cantilever {
    std::unique_ptr<Rig> rig;
    float s = 0, length = 0; // spacing; elastic length, from the clamp's last layer to the end layer
    int n = 0;
    std::vector<int> end;    // the free end's layer
    float rest = 0;          // its height at rest
    std::vector<float> tip;  // its drop, every step
    float deflection = 0;    // mean drop over the last half second
    float seconds = 0;       // of CPU
};
Cantilever cantilever(int n, float h, float L, const SoftMaterial& m, float g, float simulated, float courant = 0) {
    Cantilever c;
    c.n = n;
    c.s = h / float(n);
    const int layers = int(std::lround(L / c.s)) + 2;
    c.rig = std::make_unique<Rig>(0.5f * c.s, Vector3(0.5f * float(layers) * c.s, kBeamY, 0.0f), Vector3(float(layers) * c.s, h, h), m,
                                  Vector3(0, -g, 0), courant);
    c.rig->s.pinParticles([&](const Vector3& p) { return p.x < 2.0f * c.s; });
    for (int i : c.rig->b().particles)
        if (c.rig->x(i).x > (float(layers) - 1.0f) * c.s) c.end.push_back(i);
    c.rest = c.rig->meanY(c.end);
    c.length = (float(layers) - 2.0f) * c.s;
    const int steps = int(simulated / kDt), averaged = int(0.5f / kDt);
    const auto t0 = std::chrono::steady_clock::now();
    for (int k = 0; k < steps; ++k) {
        c.rig->s.step(kDt);
        c.tip.push_back(c.rest - c.rig->meanY(c.end));
        if (k >= steps - averaged) c.deflection += c.tip.back() / float(averaged);
    }
    c.seconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - t0).count();
    return c;
}

// Timoshenko's cantilever under its own weight q per length, elastic section b x b: bending
// qL^4 / (8 E I) plus shear qL^2 / (2 kappa G A), kappa = 5/6 for a rectangle.
float timoshenko(const SoftMaterial& m, float q, float b, float L) {
    const float A = b * b, I = b * b * b * b / 12.0f, G = shearModulus(m);
    return q * std::pow(L, 4.0f) / (8.0f * m.youngModulus * I) + q * L * L / (2.0f * (5.0f / 6.0f) * G * A);
}

// Soft rubber on a lowered gravity: the same small deflection as real rubber (5 MPa) under real
// gravity, but a shear wave slow enough that the step is cut in a handful, not fifteen.
const SoftMaterial kRubber{1000.0f, 2.5e5f, 0.3f, 10.0f};
constexpr float kBeamG = 0.5f, kBeamH = 0.1f, kBeamL = 0.4f;

} // namespace

// 1. The cantilever under its own weight against Timoshenko's beam, on three grids (4, 6, 8
//    particles across the 10 cm section): within 2 / (n - 1)^2 of the theory for the elastic core
//    the material fills - the bending error of linear tetrahedra, second order in the cells across
//    it (1.7 / (n - 1)^2 converged, so 3 cells cannot come within the 12 % once asked for) - and
//    closer to it on every finer grid.
void testSoftCantilever() {
    float previous = kInf;
    for (int n : {4, 6, 8}) {
        const Cantilever c = cantilever(n, kBeamH, kBeamL, kRubber, kBeamG, 2.0f);
        const float q = kRubber.density * kBeamG * kBeamH * kBeamH;
        const float core = timoshenko(kRubber, q, float(n - 1) * c.s, c.length), solid = timoshenko(kRubber, q, kBeamH, c.length);
        std::printf("  %d particles across: tip drops %.2f mm; Timoshenko for the elastic core %.2f mm (ratio %.3f), for the solid "
                    "section %.2f mm (ratio %.3f); %.1f s of CPU for 2 s\n",
                    n, 1000 * c.deflection, 1000 * core, c.deflection / core, 1000 * solid, c.deflection / solid, c.seconds);
        const float error = std::fabs(c.deflection / core - 1.0f);
        CHECK(error < 2.0f / sqr(float(n - 1)), "%d across: %.3f of the core's beam theory", n, c.deflection / core);
        CHECK(error < previous, "%d across: %.3f of the core's beam theory, no closer than on the coarser grid", n, c.deflection / core);
        previous = error;
    }
}

// 2. The principle of virtual work (Johann Bernoulli, 1717): at rest the elastic forces and the
//    weight together do no work on any admissible virtual displacement - seven smooth fields that
//    keep the clamp (bending, twisting, stretching, shearing and three polynomial ones) - and the
//    rest state is a minimum of the total potential energy: displaced either way along the first
//    four, the body's energy goes up. The forces are those of the energy the solver minimises.
void testSoftVirtualWork() {
    SoftMaterial rubber = kRubber;
    Cantilever c = cantilever(6, kBeamH, kBeamL, rubber, kBeamG, 2.5f);
    ParticleSystem& s = c.rig->s;
    const SoftBody& b = c.rig->b();
    const std::vector<Vector3> x = s.positions();
    // The forces at every particle: its weight, and minus the gradient of the elastic energy,
    // -sum V P b with the corotated stress P = 2 mu (F - R) + lambda (det F - 1) cof F.
    const float mu = shearModulus(rubber), lambda = lameLambda(rubber);
    std::vector<Vector3> force(x.size(), Vector3(0.0f));
    for (const SoftTet& t : b.tets) { // the shape: P = 2 mu (F - R)
        const Matrix3x3 F = deformationGradient(t, x);
        const Matrix3x3 R = extractRotation(F, t.rotation, 30).toMatrix3x3();
        const Matrix3x3 G = (F - R) * (2.0f * mu) * t.restInverse.transposed() * t.restVolume;
        const Vector3 g[4] = {-(G.col(0) + G.col(1) + G.col(2)), G.col(0), G.col(1), G.col(2)};
        for (int k = 0; k < 4; ++k) force[size_t(t.v[size_t(k)])] -= g[k];
    }
    for (const SoftNode& n : b.nodes) { // the volume per particle: -lambda (J_i - 1) sum_t V_t/4 cof F_t b_c
        double J = 0;
        for (int q = n.starBegin; q < n.starEnd; ++q) {
            const SoftTet& t = b.tets[size_t(b.nodeStar[size_t(q)])];
            J += 0.25 * t.restVolume * deformationGradient(t, x).determinant() / n.restVolume;
        }
        const float pressure = lambda * float(J - 1.0);
        for (int q = n.starBegin; q < n.starEnd; ++q) {
            const SoftTet& t = b.tets[size_t(b.nodeStar[size_t(q)])];
            const Matrix3x3 F = deformationGradient(t, x);
            const Vector3 f1 = F.col(0), f2 = F.col(1), f3 = F.col(2);
            const Matrix3x3 cof = Matrix3x3::fromColumns(cross(f2, f3), cross(f3, f1), cross(f1, f2));
            const Matrix3x3 G = cof * t.restInverse.transposed() * (0.25f * t.restVolume * pressure);
            const Vector3 g[4] = {-(G.col(0) + G.col(1) + G.col(2)), G.col(0), G.col(1), G.col(2)};
            for (int k = 0; k < 4; ++k) force[size_t(t.v[size_t(k)])] -= g[k];
        }
    }
    double residual = 0, weight = 0;
    for (int i : b.particles) {
        if (s.invMasses()[size_t(i)] == 0) continue;
        const Vector3 w = s.params.gravity / s.invMasses()[size_t(i)];
        residual += length2(force[size_t(i)] + w);
        weight += length2(w);
    }
    const float unbalanced = float(std::sqrt(residual / weight));
    // The second variation along fields that vanish at the clamp.
    const float clampX = 2.0f * c.s, L = c.length + 2.0f * c.s;
    auto field = [&](int kind, const Vector3& p) {
        const float u = std::max(0.0f, (p.x - clampX) / (L - clampX)); // 0 at the clamp, 1 at the tip
        switch (kind) {
        case 0: return Vector3(0, u * u, 0);                                        // bending
        case 1: return cross(Vector3(u, 0, 0), p - Vector3(p.x, kBeamY, 0.0f));       // twisting about the axis
        case 2: return Vector3(u, 0, 0);                                            // stretching
        default: return Vector3(0, 0, std::sin(3.0f * u) * (p.y - kBeamY) * 10.0f);   // shearing the section
        }
    };
    // The first variation (Bernoulli's principle of virtual work, the weak form of the balance):
    // along every admissible field phi the forces do no work, sum (f_i + w_i) . phi(x_i) = 0. It
    // averages the float rounding of single particles out, which the pointwise residual above
    // cannot (at 1 m from the origin a float position is good to 6e-8 m, and with gravity lowered
    // to 0.5 m/s^2 the weight of a particle is 1 / 20 000 of the forces the beam carries).
    auto extra = [&](int kind, const Vector3& p) { // smooth polynomial fields, zero at the clamp
        const float u = std::max(0.0f, (p.x - clampX) / (L - clampX)), y = (p.y - kBeamY) * 10.0f, z = p.z * 10.0f;
        switch (kind) {
        case 0: return Vector3(u * y, u * u * u, u * z);
        case 1: return Vector3(u * u * z, u * y * z, u * u);
        default: return Vector3(u * y * y, u * u * y, u * u * u * z);
        }
    };
    float worstWork = 0;
    std::printf("  first variation / weight's work scale:");
    for (int kind = 0; kind < 7; ++kind) {
        double work = 0, scaleWork = 0;
        for (int i : b.particles) {
            if (s.invMasses()[size_t(i)] == 0) continue;
            const Vector3 phi = kind < 4 ? field(kind, x[size_t(i)]) : extra(kind - 4, x[size_t(i)]);
            const Vector3 w = s.params.gravity / s.invMasses()[size_t(i)];
            work += dot(force[size_t(i)] + w, phi);
            scaleWork += length(w) * length(phi);
        }
        const float ratio = float(std::fabs(work) / scaleWork);
        worstWork = std::max(worstWork, ratio);
        std::printf(" %.2e", ratio);
    }
    std::printf("\n");
    const double eps = 1e-3, energy0 = elasticEnergy(b, x) + gravityPotential(s, b, x);
    bool minimum = true;
    std::printf("  virtual work at rest: the unbalanced force is %.2f %% of the weight; second variations", 100 * unbalanced);
    for (int kind = 0; kind < 4; ++kind) {
        std::vector<Vector3> plus = x, minus = x;
        for (int i : b.particles) {
            if (s.invMasses()[size_t(i)] == 0) continue;
            const Vector3 d = field(kind, x[size_t(i)]) * float(eps);
            plus[size_t(i)] += d;
            minus[size_t(i)] -= d;
        }
        const double second = elasticEnergy(b, plus) + gravityPotential(s, b, plus) + elasticEnergy(b, minus) + gravityPotential(s, b, minus) -
                              2 * energy0;
        std::printf(" %s %.2e J", kind == 0 ? "bend" : kind == 1 ? "twist" : kind == 2 ? "stretch" : "shear", second);
        minimum = minimum && second > 0;
    }
    std::printf("\n");
    CHECK(worstWork < 0.01f, "along a virtual displacement the forces at rest do %.2f %% of the weight's work", 100 * worstWork);
    CHECK(minimum, "the rest state is not a minimum of the potential energy");
}

// 3. The same beam without internal friction swings about its bent shape: its first bending
//    frequency against 3.516 sqrt(EI / (rho A L^4)) for the elastic core; and omega^2 delta / g =
//    3.516^2 / 8 for any stiffness - the dynamics against the statics, free of the grid. Four
//    seconds: the period is 0.9 s, and the period is read off four downward crossings.
void testSoftBeamFrequency() {
    SoftMaterial rubber = kRubber;
    rubber.damping = 0;
    const Cantilever c = cantilever(6, kBeamH, kBeamL, rubber, kBeamG, 4.0f);
    float mean = 0;
    for (float y : c.tip) mean += y / float(c.tip.size());
    std::vector<float> up; // times the tip passes its mean going down
    for (size_t k = 1; k < c.tip.size(); ++k)
        if (c.tip[k - 1] < mean && c.tip[k] >= mean) up.push_back((float(k) - (c.tip[k] - mean) / (c.tip[k] - c.tip[k - 1])) * kDt);
    const float period = up.size() >= 3 ? (up.back() - up.front()) / float(up.size() - 1) : 0;
    const float omega = period > 0 ? 2 * kPi / period : 0;
    const float core = float(c.n - 1) * c.s, I = core * core * core * core / 12.0f;
    const float theory = 3.516f * std::sqrt(rubber.youngModulus * I / (rubber.density * kBeamH * kBeamH * std::pow(c.length, 4.0f)));
    const float invariant = omega * omega * mean / kBeamG, invariantTheory = 3.516f * 3.516f / 8.0f;
    std::printf("  first bending mode: %.2f rad/s over %zu swings, Euler-Bernoulli for the core %.2f rad/s (ratio %.3f); "
                "w^2 d / g %.3f, theory %.3f\n",
                omega, up.size() > 0 ? up.size() - 1 : 0, theory, omega / theory, invariant, invariantTheory);
    CHECK(up.size() >= 4, "the beam does not swing (%zu crossings)", up.size());
    CHECK(std::fabs(omega / theory - 1.0f) < 0.15f, "the first bending frequency is %.3f of beam theory", omega / theory);
    CHECK(std::fabs(invariant / invariantTheory - 1.0f) < 0.1f, "w^2 delta / g = %.3f, beam theory %.3f", invariant, invariantTheory);
}

// A bar of section w x w, `across` particles, hanging from its top two layers. The bottom layer's
// drop from the clamp, the width change at mid-height and the volume change, over the last half second.
struct HangingBar {
    float extension = 0, elasticLength = 0, lateralStrain = 0, volumeChange = 0, coreStress = 0;
};
HangingBar hangingBar(const SoftMaterial& m, float w, float Lbar, int across) {
    const float s = w / float(across);
    const int layers = int(std::lround(Lbar / s));
    Rig rig(0.5f * s, Vector3(0, 1.0f, 0), Vector3(w, float(layers) * s, w), m, Vector3(0, -9.81f, 0));
    const float top = 1.0f + 0.5f * float(layers) * s;
    rig.s.pinParticles([&](const Vector3& p) { return p.y > top - 2.0f * s; });
    const float clampY = top - 1.5f * s, bottomY = top - (float(layers) - 0.5f) * s, midY = 0.5f * (clampY + bottomY);
    std::vector<int> bottom, middle;
    for (int i : rig.b().particles) {
        if (rig.x(i).y < bottomY + 0.5f * s) bottom.push_back(i);
        if (std::fabs(rig.x(i).y - midY) < 0.5f * s) middle.push_back(i);
    }
    auto width = [&] { // mean distance of the middle layer from the axis
        float sum = 0;
        for (int i : middle) sum += std::fabs(rig.x(i).x) + std::fabs(rig.x(i).z);
        return sum / float(2 * middle.size());
    };
    const float width0 = width(), bottom0 = rig.meanY(bottom);
    double volume0 = 0;
    for (const SoftTet& t : rig.b().tets) volume0 += t.restVolume;
    HangingBar r;
    const int steps = int(2.0f / kDt), averaged = int(0.5f / kDt);
    for (int k = 0; k < steps; ++k) {
        rig.s.step(kDt);
        if (k < steps - averaged) continue;
        r.extension += (bottom0 - rig.meanY(bottom)) / float(averaged);
        r.lateralStrain += (width() / width0 - 1.0f) / float(averaged);
        double volume = 0;
        for (const SoftTet& t : rig.b().tets) volume += t.restVolume * deformationGradient(t, rig.s.positions()).determinant();
        r.volumeChange += float(volume / volume0 - 1.0) / float(averaged);
    }
    r.elasticLength = clampY - bottomY;
    // The weight below mid-height (the full section) over the elastic core's area.
    r.coreStress = m.density * 9.81f * (midY - bottomY + 0.5f * s) * sqr(float(across) / float(across - 1));
    return r;
}

// 4. A bar hanging under its own weight stretches by rho g L^2 / (2E) (times the section over the
//    core) and narrows by Poisson's ratio: at mid-height the lateral strain is -nu times the axial.
// 5. The same bar nearly incompressible (nu = 0.49): it must still stretch as much - a tetrahedral
//    mesh with a volume term per tetrahedron can lock and hardly stretch - and keep its volume.
void testSoftHangingBar() {
    float extension[2] = {0, 0};
    int k = 0;
    for (const float nu : {0.3f, 0.49f}) {
        SoftMaterial jelly{1000.0f, 5e4f, nu, 10.0f};
        const HangingBar b = hangingBar(jelly, 0.1f, 0.5f, 5);
        const float area = sqr(5.0f / 4.0f);
        const float theory = jelly.density * 9.81f * b.elasticLength * b.elasticLength / (2.0f * jelly.youngModulus) * area;
        const float axial = b.coreStress / jelly.youngModulus, lateralTheory = -nu * axial;
        std::printf("  nu %.2f: bar of %.2f m stretches %.2f mm, theory for the core %.2f mm (ratio %.3f; solid section %.2f mm); "
                    "mid-height lateral strain %.4f, -nu x axial %.4f; volume change %.4f\n",
                    nu, b.elasticLength, 1000 * b.extension, 1000 * theory, b.extension / theory, 1000 * theory / area, b.lateralStrain,
                    lateralTheory, b.volumeChange);
        CHECK(std::fabs(b.extension / theory - 1.0f) < 0.15f, "nu %.2f: the bar stretches %.3f of theory", nu, b.extension / theory);
        CHECK(std::fabs(b.lateralStrain / lateralTheory - 1.0f) < 0.3f, "nu %.2f: Poisson contraction %.4f, theory %.4f", nu, b.lateralStrain,
              lateralTheory);
        extension[k++] = b.extension;
        if (nu > 0.4f) CHECK(std::fabs(b.volumeChange) < 0.004f, "nu 0.49: the volume changes by %.4f", b.volumeChange);
    }
    CHECK(extension[1] > 0.85f * extension[0], "nu 0.49 locks: it stretches %.2f mm against %.2f mm at nu 0.3", 1000 * extension[1],
          1000 * extension[0]);
}

// 6. A box spinning and tumbling in weightless space: nothing outside acts on it, so its momentum
//    and angular momentum stay (Noether) while it stretches under the spin.
void testSoftSpinMomentum() {
    SoftMaterial foam{200.0f, 1e5f, 0.3f, 3.0f};
    Rig rig(0.01f, Vector3(0.0f), Vector3(0.3f, 0.12f, 0.08f), foam, Vector3(0.0f));
    const Vector3 w(0.3f, 6.0f, 0.2f);
    for (int i : rig.b().particles) rig.s.addVelocity(i, cross(w, rig.x(i)));
    auto momenta = [&](Vector3& P, Vector3& L) {
        Vector3 c(0.0f);
        float M = 0;
        for (int i : rig.b().particles) c += rig.x(i) / rig.s.invMasses()[size_t(i)], M += 1.0f / rig.s.invMasses()[size_t(i)];
        c /= M;
        P = L = Vector3(0.0f);
        for (int i : rig.b().particles) {
            const float m = 1.0f / rig.s.invMasses()[size_t(i)];
            P += rig.s.velocities()[size_t(i)] * m;
            L += cross(rig.x(i) - c, rig.s.velocities()[size_t(i)] * m);
        }
    };
    Vector3 P0, L0, P1, L1;
    momenta(P0, L0);
    for (int k = 0; k < int(2.0f / kDt); ++k) rig.s.step(kDt);
    momenta(P1, L1);
    const float drift = length(L1 - L0) / length(L0);
    std::printf("  spinning box: |P| %.2e -> %.2e kg m/s, L (%.4f %.4f %.4f) -> (%.4f %.4f %.4f) kg m^2/s, change %.3f%% in 2 s\n", length(P0),
                length(P1), L0.x, L0.y, L0.z, L1.x, L1.y, L1.z, 100 * drift);
    CHECK(length(P1) < 1e-4f * rig.mass(), "a free spinning body drifts: |P| = %.2e", length(P1));
    CHECK(drift < 0.01f, "angular momentum changed by %.2f%%", 100 * drift);
}

// 7. A cube clamped by one half, the other half grabbed and dragged by the mouse right through the
//    clamped half and beyond (its tetrahedra squashed flat and turned inside out on the way), then
//    let go: it comes back to its rest shape.
void testSoftDragThrough() {
    SoftMaterial rubber{500.0f, 2e5f, 0.4f, 3.0f};
    Rig rig(0.01f, Vector3(0.0f), Vector3(0.16f), rubber, Vector3(0.0f));
    std::vector<Vector3> rest;
    for (int i : rig.b().particles) rest.push_back(rig.x(i));
    rig.s.pinParticles([](const Vector3& p) { return p.x < -0.02f; });
    CHECK(rig.s.grab(Vector3(0.07f, 0, 0)), "nothing to grab");
    int inverted = 0;
    for (int k = 0; k <= 90; ++k) { // 0.5 s: from x = 0.07 to -0.13, through the clamped half
        rig.s.setGrabTarget(Vector3(0.07f - 0.2f * float(k) / 90.0f, 0, 0));
        rig.s.step(kDt);
        int now = 0;
        for (const SoftTet& t : rig.b().tets) now += deformationGradient(t, rig.s.positions()).determinant() <= 0;
        inverted = std::max(inverted, now);
    }
    rig.s.releaseGrab();
    for (int k = 0; k < int(3.0f / kDt); ++k) rig.s.step(kDt);
    const float error = shapeError(rig.s, rig.b(), rest);
    float minJ = kInf;
    for (const SoftTet& t : rig.b().tets) minJ = std::min(minJ, deformationGradient(t, rig.s.positions()).determinant());
    std::printf("  dragged through its clamped half: up to %d of %zu tetrahedra inside out; 3 s after letting go shape error %.2f%%, "
                "smallest volume ratio %.3f, finite %d\n",
                inverted, rig.b().tets.size(), 100 * error, minJ, int(allFinite(rig.s)));
    CHECK(inverted > 0, "the drag did not turn any tetrahedron inside out");
    CHECK(allFinite(rig.s) && error < 0.02f && minJ > 0.9f, "the cube did not come back: shape error %.2f%%, volume ratio %.3f", 100 * error,
          minJ);
}

// 8. A jelly cube resting on the floor: after settling it is still - no jitter, no creep.
void testSoftRest() {
    SoftMaterial jelly{1000.0f, 3e4f, 0.45f, 3.0f};
    ParticleSystem s;
    s.params.particleRadius = 0.01f;
    s.reset(AABB({-0.5f, 0, -0.5f}, {0.5f, 1, 0.5f}));
    TriMesh box = primitives::box(Vector3(0.1f));
    box.translate({0, 0.1f, 0});
    const int b = s.addSoftBody(box, jelly, Vector3(1));
    for (int k = 0; k < int(2.0f / kDt); ++k) s.step(kDt);
    std::vector<Vector3> lo, hi;
    for (int i : s.softBodies()[size_t(b)].particles) lo.push_back(s.positions()[size_t(i)]), hi.push_back(s.positions()[size_t(i)]);
    for (int k = 0; k < int(3.0f / kDt); ++k) {
        s.step(kDt);
        const auto& ids = s.softBodies()[size_t(b)].particles;
        for (size_t n = 0; n < ids.size(); ++n) lo[n] = vmin(lo[n], s.positions()[size_t(ids[n])]), hi[n] = vmax(hi[n], s.positions()[size_t(ids[n])]);
    }
    float range = 0;
    for (size_t n = 0; n < lo.size(); ++n) range = std::max(range, maxComp(hi[n] - lo[n]));
    std::printf("  jelly cube at rest from 2 s to 5 s: every particle stays within %.3f mm\n", 1000 * range);
    CHECK(range < 0.5e-3f, "a resting jelly moves: %.3f mm", 1000 * range);
}

// 9. The static answer does not depend on the step: the cantilever's deflection with the step cut
//    so that a shear wave crosses at most 1, 1/2 and 1/4 of a cell per step.
// 10. Energy never grows: the same beam, undamped, swinging for 1.5 s - its total energy (kinetic,
//    elastic, the weight's) ends below where it began.
void testSoftTimeStepAndEnergy() {
    float d[3];
    const float courant[3] = {1.0f, 0.5f, 0.25f};
    for (int k = 0; k < 3; ++k) d[k] = cantilever(6, kBeamH, kBeamL, kRubber, kBeamG, 2.0f, courant[k]).deflection;
    const float spread = (std::max({d[0], d[1], d[2]}) - std::min({d[0], d[1], d[2]})) / d[2];
    std::printf("  Courant 1 / 0.5 / 0.25: tip drops %.3f / %.3f / %.3f mm (spread %.2f%%)\n", 1000 * d[0], 1000 * d[1], 1000 * d[2], 100 * spread);
    CHECK(spread < 0.03f, "the static deflection depends on the step: %.2f%%", 100 * spread);

    SoftMaterial undamped = kRubber;
    undamped.damping = 0;
    Cantilever c = cantilever(6, kBeamH, kBeamL, undamped, kBeamG, 0.0f);
    ParticleSystem& s = c.rig->s;
    const SoftBody& b = c.rig->b();
    auto total = [&] { return kineticEnergy(s, b) + elasticEnergy(b, s.positions()) + gravityPotential(s, b, s.positions()); };
    const double start = total();
    double swing = 0; // the largest kinetic energy of the swing: the scale of what it trades
    for (int k = 0; k < int(1.5f / kDt); ++k) {
        s.step(kDt);
        swing = std::max(swing, double(kineticEnergy(s, b)));
    }
    const double end = total();
    std::printf("  undamped beam: total energy %+.6f J after 1.5 s (the swing trades up to %.5f J of kinetic energy)\n", end - start, swing);
    CHECK(end - start < 0.01 * swing, "the energy grew by %.6f J", end - start);
}

namespace {

// A jelly asleep on a tilted static platform: neither side of their contact can move, and the
// contact is skipped - divided by its zero inverse mass, it turned every velocity to NaN.
bool sleepsFiniteOnStaticSlope() {
    ParticleSystem s;
    RigidWorld w;
    const AABB domain({-0.5f, 0, -0.5f}, {0.5f, 1.0f, 0.5f});
    w.setDomain(domain);
    s.setRigidWorld(&w);
    s.reset(domain);
    w.addBox({0, 0.1f, 0}, Vector3(0.3f, 0.1f, 0.3f), Quaternion::fromAxisAngle(normalize(Vector3(0.3f, 0, 1)), 0.05f), 0.0f, Vector3(1));
    TriMesh jelly = primitives::box(Vector3(0.08f));
    jelly.translate({0, 0.285f, 0});
    s.addSoftBody(jelly, SoftMaterial{400.0f, 3e4f, 0.45f, 3.0f}, Vector3(1));
    bool finite = true;
    for (int k = 0; k < int(2.0f / kDt); ++k) {
        w.step(kDt), s.step(kDt);
        finite = finite && allFinite(s);
    }
    return finite && s.sleepingSoftBodies() == 1;
}

} // namespace

// 10. Sleeping (as the rigid world's islands): a jelly lying still falls asleep and then does not
//     move by a single bit; a box dropped on it wakes it and comes to rest on it (a sleeping body
//     must not let a load through); both fall asleep again; the mouse wakes it. One asleep on a
//     static slope stays finite.
void testSoftSleep() {
    ParticleSystem s;
    RigidWorld w;
    const AABB domain({-0.5f, 0, -0.5f}, {0.5f, 1.0f, 0.5f});
    w.setDomain(domain);
    s.setRigidWorld(&w);
    s.reset(domain);
    TriMesh jelly = primitives::box(Vector3(0.1f));
    jelly.translate({0, 0.1f, 0});
    s.addSoftBody(jelly, SoftMaterial{400.0f, 3e4f, 0.45f, 3.0f}, Vector3(1));
    auto run = [&](float seconds) {
        for (int k = 0; k < int(seconds / kDt); ++k) w.step(kDt), s.step(kDt);
    };
    run(2.0f);
    const bool asleep = s.sleepingSoftBodies() == 1;
    const std::vector<Vector3> before = s.positions();
    run(0.5f);
    bool frozen = true;
    for (size_t i = 0; i < before.size(); ++i) {
        const Vector3& a = before[i], & b = s.positions()[i];
        frozen = frozen && a.x == b.x && a.y == b.y && a.z == b.z;
    }
    float top = -kInf;
    for (const Vector3& p : s.positions()) top = std::max(top, p.y + s.params.particleRadius);
    const int box = w.addBox({0, top + 0.05f + 0.3f, 0}, Vector3(0.05f), Quaternion(), 500.0f, Vector3(1));
    run(0.4f);
    const bool wokeByBox = s.sleepingSoftBodies() == 0;
    run(3.0f);
    float jellyTop = -kInf;
    for (const Vector3& p : s.positions())
        if (std::fabs(p.x) < 0.04f && std::fabs(p.z) < 0.04f) jellyTop = std::max(jellyTop, p.y + s.params.particleRadius);
    const float gap = (w.bodies()[size_t(box)].pos.y - 0.05f) - jellyTop;
    const bool asleepAgain = s.sleepingSoftBodies() == 1;
    s.grab(s.positions()[0]);
    run(0.1f);
    const bool wokeByHand = s.sleepingSoftBodies() == 0;
    s.releaseGrab();
    const bool slope = sleepsFiniteOnStaticSlope();
    std::printf("  jelly asleep after 2 s %d, still to the bit %d; a box dropped on it wakes it %d and rests on it: gap %.1f mm; "
                "asleep again %d; the mouse wakes it %d; asleep on a static slope and finite %d\n",
                int(asleep), int(frozen), int(wokeByBox), 1000 * gap, int(asleepAgain), int(wokeByHand), int(slope));
    CHECK(asleep && frozen, "the resting jelly does not sleep (asleep %d, frozen %d)", int(asleep), int(frozen));
    CHECK(wokeByBox && std::fabs(gap) < 0.5f * s.params.particleRadius, "the box on the sleeping jelly: woke %d, gap %.1f mm", int(wokeByBox),
          1000 * gap);
    CHECK(asleepAgain && wokeByHand, "asleep again %d, woken by the mouse %d", int(asleepAgain), int(wokeByHand));
    CHECK(slope, "a jelly asleep on a static slope: not asleep, or not finite");
}
