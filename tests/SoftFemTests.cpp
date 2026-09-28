// The Neo-Hookean soft body (particles/SoftTets.cpp) against the theory of elasticity: a cantilever
// against Euler-Bernoulli (tip deflection and first bending frequency, two resolutions), a nearly
// incompressible cube keeps its volume, Coulomb friction on an incline (sticks below the friction
// angle, slides with a = g (sin - mu cos) above it), energy never made, a soft ball rolls off the
// floor with no invisible wall, and what a frame of six soft barrels costs against shape matching.
#include "TestRunner.h"
#include "Tests.h"

#include "scene/SceneGraph.h"

#include <string>

namespace {

// The centre of mass and the mean velocity of a soft body's particles.
Vector3 bodyCentre(const ParticleSystem& s, const SoftBody& b) {
    Vector3 c(0.0f);
    for (int i : b.particles) c += s.positions()[size_t(i)];
    return c / float(b.particles.size());
}
Vector3 bodyVelocity(const ParticleSystem& s, const SoftBody& b) {
    Vector3 v(0.0f);
    for (int i : b.particles) v += s.velocities()[size_t(i)];
    return v / float(b.particles.size());
}

// The volume of a body's tetrahedra now, sum of det(D_s) / 6.
double tetVolume(const ParticleSystem& s, const SoftBody& b) {
    double v = 0;
    for (const SoftTet& t : b.tets) {
        const Vector3 x0 = s.positions()[size_t(t.v[0])];
        v += Matrix3x3::fromColumns(s.positions()[size_t(t.v[1])] - x0, s.positions()[size_t(t.v[2])] - x0,
                                    s.positions()[size_t(t.v[3])] - x0).determinant() / 6.0;
    }
    return v;
}

// A cantilever of square cross-section H x H and length L (the hull of the particle centres: the
// mesh is one spacing larger, its particles sit half a spacing in), clamped at z = 0 by pinning its
// first layer of particles, released straight under the load q = rho g' (H + s)^2 per metre - the
// particles' weight. The tip (the mean of the last layer) swings about the static deflection:
// averaged over whole periods it is that deflection, and the periods are the first mode's.
struct BeamResult {
    float deflection = 0, frequency = 0, theoryDeflection = 0, theoryFrequency = 0;
    int particles = 0, tets = 0, smallSteps = 0;
};

BeamResult swingBeam(float radius, float H, float L, const SoftMaterial& m, float g) {
    ParticleSystem s;
    s.params.particleRadius = radius;
    s.params.gravity = Vector3(0, -g, 0);
    s.params.softMaxSubsteps = 64;
    s.reset(AABB({-1, -1, -0.5f}, {1, 1, L + 0.5f}));
    const float sp = 2.0f * radius;
    TriMesh bar = primitives::box(Vector3(0.5f * (H + sp), 0.5f * (H + sp), 0.5f * (L + sp)));
    bar.translate({0, 0, 0.5f * L});
    const int b = s.addSoftBody(bar, m, Vector3(1));
    s.pinParticles([&](const Vector3& x) { return x.z < 0.25f * sp; });
    const SoftBody& body = s.softBodies()[size_t(b)];
    std::vector<int> tip;
    for (int i : body.particles)
        if (s.positions()[size_t(i)].z > L - 0.25f * sp) tip.push_back(i);
    auto tipY = [&]() {
        double y = 0;
        for (int i : tip) y += s.positions()[size_t(i)].y;
        return float(y / double(tip.size()));
    };
    const float y0 = tipY(), dt = 1.0f / 180.0f;
    std::vector<float> ys;
    for (int k = 0; k < 720; ++k) {
        s.step(dt);
        ys.push_back(y0 - tipY());
    }
    // Upward crossings of the mean: the periods; the mean over whole periods: the deflection.
    double mean = 0;
    for (float y : ys) mean += y;
    mean /= double(ys.size());
    std::vector<float> crossings;
    for (size_t k = 1; k < ys.size(); ++k)
        if (ys[k - 1] < mean && ys[k] >= mean) crossings.push_back((float(k - 1) + float((mean - ys[k - 1]) / (ys[k] - ys[k - 1]))) * dt);
    BeamResult r;
    if (crossings.size() >= 2) {
        const size_t first = size_t(crossings.front() / dt), last = size_t(crossings.back() / dt);
        double sum = 0;
        for (size_t k = first; k < last; ++k) sum += ys[k];
        r.deflection = float(sum / double(last - first));
        r.frequency = float(crossings.size() - 1) / (crossings.back() - crossings.front());
    }
    const Lame lame = lameParameters(m);
    const float E = 2.0f * lame.mu * (1.0f + m.poissonRatio), I = H * H * H * H / 12.0f;
    const float massPerLength = m.density * (H + sp) * (H + sp);
    r.theoryDeflection = massPerLength * g * L * L * L * L / (8.0f * E * I);                     // q L^4 / (8 E I)
    r.theoryFrequency = 1.875f * 1.875f / (2.0f * kPi * L * L) * std::sqrt(E * I / massPerLength); // Euler-Bernoulli, mode 1
    r.particles = int(body.particles.size());
    r.tets = int(body.tets.size());
    r.smallSteps = s.softSmallStepsLast();
    return r;
}

} // namespace

// Euler-Bernoulli (Timoshenko & Gere, "Mechanics of Materials"): a cantilever under a uniform
// load q deflects at its tip by q L^4 / (8 E I) and swings first at f1 = 1.875^2 / (2 pi L^2)
// sqrt(E I / (rho A)). L / H = 6 keeps the shear part of the deflection (Timoshenko's beam) near
// 3 %; the load g' = 1 m/s^2 keeps the deflection near 5 % of L, where the beam is linear.
void testSoftFemCantilever() {
    SoftMaterial m;
    m.density = 1000;
    m.youngModulus = 8e5f;
    m.poissonRatio = 0.3f;
    m.damping = 0;
    const float H = 0.1f, L = 0.6f;
    const BeamResult coarse = swingBeam(0.0125f, H, L, m, 1.0f), fine = swingBeam(0.00625f, H, L, m, 1.0f);
    auto error = [](float ours, float theory) { return std::fabs(ours / theory - 1.0f); };
    for (const BeamResult* r : {&coarse, &fine})
        std::printf("  cantilever %s (%d particles, %d tetrahedra, %d small steps): tip deflection %.2f mm (Euler-Bernoulli %.2f, "
                    "%+.1f%%), first frequency %.3f Hz (%.3f, %+.1f%%)\n",
                    r == &coarse ? "4 x 4 cells" : "8 x 8 cells", r->particles, r->tets, r->smallSteps, 1000 * r->deflection,
                    1000 * r->theoryDeflection, 100 * (r->deflection / r->theoryDeflection - 1), r->frequency, r->theoryFrequency,
                    100 * (r->frequency / r->theoryFrequency - 1));
    CHECK(error(fine.deflection, fine.theoryDeflection) < 0.10f, "tip deflection %.2f mm vs %.2f", 1000 * fine.deflection,
          1000 * fine.theoryDeflection);
    CHECK(error(fine.frequency, fine.theoryFrequency) < 0.10f, "first frequency %.3f Hz vs %.3f", fine.frequency, fine.theoryFrequency);
    CHECK(error(fine.deflection, fine.theoryDeflection) < error(coarse.deflection, coarse.theoryDeflection),
          "the deflection must converge with resolution");
}

// A cube of Poisson's ratio 0.49 (rubber-like) squashed by its weight on the floor: it bulges
// sideways and keeps its volume - the hydrostatic constraint holds det F near gamma = 1 + mu/lambda.
void testSoftFemVolume() {
    ParticleSystem s;
    s.params.particleRadius = 0.01f;
    s.reset(AABB({-1, 0, -1}, {1, 1, 1}));
    SoftMaterial m;
    m.density = 1000;
    m.youngModulus = 5e3f;
    m.poissonRatio = 0.49f;
    TriMesh cube = primitives::box(Vector3(0.1f));
    cube.translate({0, 0.1f, 0});
    const int b = s.addSoftBody(cube, m, Vector3(1));
    const double rest = tetVolume(s, s.softBodies()[size_t(b)]);
    const float top0 = bodyCentre(s, s.softBodies()[size_t(b)]).y;
    double worst = 0, settled = 0; // the largest change while it squashes and swings; the mean change over the last second
    for (int k = 0; k < 540; ++k) {
        s.step(1.0f / 180.0f);
        const double change = tetVolume(s, s.softBodies()[size_t(b)]) / rest - 1.0;
        worst = std::max(worst, std::fabs(change));
        if (k >= 360) settled += change / 180.0;
    }
    const float sag = top0 - bodyCentre(s, s.softBodies()[size_t(b)]).y;
    // Theory: the volume changes by -p / K, K = E / (3 (1 - 2 nu)) the bulk modulus and p the mean
    // pressure. sigma_yy grows from 0 at the top to rho g H at the floor: with free sides p = rho g
    // H / 6 on average; with the sides held (the floor's friction holds the bottom) up to rho g H / 2.
    const float K = m.youngModulus / (3.0f * (1.0f - 2.0f * m.poissonRatio)), H = 0.18f, p = m.density * 9.81f * H;
    std::printf("  cube nu = 0.49 (%zu tetrahedra) on the floor: centre sags %.1f mm; volume under the load %+.2f%% (theory -p/K from "
                "%.2f%% free to %.2f%% held), at most %.2f%% while it swings\n",
                s.softBodies()[size_t(b)].tets.size(), 1000 * sag, 100 * settled, -100 * p / 6.0f / K, -100 * p / 2.0f / K, 100 * worst);
    CHECK(sag > 0.002f, "the cube must squash under its weight (%.1f mm)", 1000 * sag);
    CHECK(std::fabs(settled) < 0.01, "the volume under load changed by %.2f%%", 100 * settled);
    CHECK(worst < 0.02, "the volume changed by %.2f%% while swinging", 100 * worst);
}

// A rubber cube on a slope (the gravity tilted by theta over a flat floor, the same thing), friction
// 0.5: at tan theta = 0.3 it must stay; at tan theta = 0.8 it slides with a = g (sin - mu cos).
static void slideOnSlope(float tanTheta, float& drift, float& acceleration) {
    ParticleSystem s;
    s.params.particleRadius = 0.01f;
    const float theta = std::atan(tanTheta);
    s.params.gravity = Vector3(9.81f * std::sin(theta), -9.81f * std::cos(theta), 0);
    s.reset(AABB({-1, 0, -1}, {8, 1, 1}));
    SoftMaterial m;
    m.density = 1000;
    m.youngModulus = 1e5f;
    m.poissonRatio = 0.4f;
    m.friction = 0.5f;
    TriMesh cube = primitives::box(Vector3(0.05f));
    cube.translate({0, 0.05f, 0});
    const int b = s.addSoftBody(cube, m, Vector3(1));
    const SoftBody& body = s.softBodies()[size_t(b)];
    auto bottom = [&]() { // the mean x of the layer on the floor
        double x = 0;
        int n = 0;
        for (int i : body.particles)
            if (s.positions()[size_t(i)].y < 0.02f) x += s.positions()[size_t(i)].x, ++n;
        return float(x / std::max(n, 1));
    };
    // The acceleration from the mean velocity at t = 0.3 s and 1 s: sliding all the while and slower
    // than the particles' speed limit (half a kernel radius per substep, 3.6 m/s here) - it is
    // reached at 1.6 s. The drift of the layer on the floor from t = 0.5 s (settled) to 2.5 s.
    const float dt = 1.0f / 180.0f;
    float x1 = 0, v1 = 0;
    for (int k = 1; k <= 450; ++k) {
        s.step(dt);
        if (k == 54) v1 = bodyVelocity(s, body).x;
        if (k == 180) acceleration = (bodyVelocity(s, body).x - v1) / (126 * dt);
        if (k == 90) x1 = bottom();
    }
    drift = bottom() - x1;
}

void testSoftFemFriction() {
    float stayDrift, stayA, slideDrift, slideA;
    slideOnSlope(0.3f, stayDrift, stayA);
    slideOnSlope(0.8f, slideDrift, slideA);
    const float theta = std::atan(0.8f), expected = 9.81f * (std::sin(theta) - 0.5f * std::cos(theta));
    std::printf("  rubber cube on a slope, mu 0.5: tan 0.3 drifts %.2f mm in 2 s; tan 0.8 slides at %.3f m/s^2 "
                "(g (sin - mu cos) = %.3f, %+.1f%%)\n",
                1000 * stayDrift, slideA, expected, 100 * (slideA / expected - 1));
    CHECK(std::fabs(stayDrift) < 0.002f, "below the friction angle the cube slid %.1f mm", 1000 * stayDrift);
    CHECK(std::fabs(slideA / expected - 1) < 0.1f, "sliding acceleration %.3f vs %.3f", slideA, expected);
}

// A soft ball with no damping dropped on the floor: kinetic + potential + elastic energy never
// rises above where it started (XPBD and the position-based contacts may only lose energy).
void testSoftFemEnergy() {
    ParticleSystem s;
    s.params.particleRadius = 0.01f;
    s.reset(AABB({-1, 0, -1}, {1, 2, 1}));
    SoftMaterial m;
    m.density = 1000;
    m.youngModulus = 5e4f;
    m.poissonRatio = 0.4f;
    m.damping = 0;
    TriMesh ball = primitives::sphere(0.1f, 24, 12);
    ball.translate({0, 0.5f, 0});
    s.addSoftBody(ball, m, Vector3(1));
    auto energy = [&]() {
        double e = s.softElasticEnergy();
        for (size_t i = 0; i < s.size(); ++i) {
            const double mass = 1.0 / s.invMasses()[i];
            e += 0.5 * mass * double(length2(s.velocities()[i])) + mass * 9.81 * s.positions()[i].y;
        }
        return e;
    };
    const double start = energy();
    double highest = start, lowest = start;
    for (int k = 0; k < 360; ++k) {
        s.step(1.0f / 180.0f);
        highest = std::max(highest, energy());
        lowest = std::min(lowest, energy());
    }
    std::printf("  soft ball dropped from 0.4 m, no damping: energy %.4f J at the start, highest %.4f J (%+.3f%%), after 2 s %.4f J\n",
                start, highest, 100 * (highest / start - 1), energy());
    CHECK(highest <= start * 1.001, "energy rose from %.4f to %.4f J", start, highest);
    CHECK(lowest > 0.2 * start, "the ball lost too much energy (%.4f of %.4f J)", lowest, start);
}

// The user's world without liquid: only the ground holds a soft body. A rubber ball given 3 m/s
// rolls off the 3 m floor plate onto the ground and on beyond the particles' box (x = 2 m).
void testSoftFemNoInvisibleWalls() {
    SceneGraph g;
    std::string error;
    const std::string scene = "world gravity 0 -9.81 0 size 4 4 4 gas 0 magneticGas 0\n"
                              "entity \"Пол\"\n  object id 1 visible 1 locked 1\n"
                              "  shape plane size 3 0.02 3 position 0 0.01 0 rotation 0 0 0 color 0.6 0.62 0.66\n"
                              "  rigid density 500 friction 0.5 restitution 0.2 fixed 1 velocity 0 0 0 spin 0 0 0\nend\n"
                              "entity \"Мяч\"\n  object id 2 visible 1 locked 0\n"
                              "  shape sphere size 0.3 0.3 0.3 position 0 0.17 0 rotation 0 0 0 color 0.9 0.3 0.3\n"
                              "  soft density 1100 young 1e6 poisson 0.47 friction 0.8\nend\n";
    CHECK(g.load(scene, error), "the scene does not load: %s", error.c_str());
    Simulation sim;
    sim.load(std::make_unique<GraphScene>(g));
    for (size_t i = 0; i < sim.particles.size(); ++i) sim.particles.addVelocity(int(i), Vector3(3, 0, 0));
    for (int f = 0; f < 180; ++f) sim.stepFrame();
    const Vector3 c = bodyCentre(sim.particles, sim.particles.softBodies()[0]);
    std::printf("  soft ball rolled at 3 m/s: after 3 s at x %.2f m, height %.2f m (particle box ends at x %.1f m)\n", c.x, c.y,
                sim.particles.domain().hi.x);
    CHECK(c.x > 2.5f && c.y < 0.3f, "the ball must roll on over the ground (x %.2f, y %.2f)", c.x, c.y);
}

// What a frame of the user's six soft barrels costs: the Neo-Hookean tetrahedra against the
// legacy shape matching they replace (the same scene, contacts and machine).
static double barrelFrameMs(bool shapeMatching) {
    std::string text = "world gravity 0 -9.81 0 size 4 4 4 gas 0 magneticGas 0\n"
                       "entity \"Пол\"\n  object id 1 visible 1 locked 1\n"
                       "  shape plane size 3 0.02 3 position 0 0.01 0 rotation 0 0 0 color 0.6 0.62 0.66\n"
                       "  rigid density 500 friction 0.5 restitution 0.2 fixed 1 velocity 0 0 0 spin 0 0 0\nend\n";
    for (int k = 0; k < 6; ++k) {
        char line[320];
        std::snprintf(line, sizeof line,
                      "entity \"Бочка %d\"\n  object id %d visible 1 locked 0\n"
                      "  shape cylinder size 0.3 0.3 0.3 position 0 %.4f 0 rotation 0 0 0 color 0.5 0.5 0.5\n"
                      "  soft density 500 %sstiffness 0.5\nend\n",
                      k + 1, k + 2, 0.17f + 0.302f * float(k), shapeMatching ? "model shape-matching " : "");
        text += line;
    }
    SceneGraph g;
    std::string error;
    g.load(text, error);
    Simulation sim;
    sim.load(std::make_unique<GraphScene>(g));
    double total = 0;
    for (int f = 0; f < 180; ++f) {
        sim.stepFrame();
        total += Probe::snapshot().value("frame/step ms");
    }
    return total / 180.0;
}

void testSoftFemPerformance() {
    const double legacy = barrelFrameMs(true), tets = barrelFrameMs(false);
    std::printf("  six soft barrels, 3 s: shape matching %.1f ms per frame, Neo-Hookean tetrahedra %.1f ms (x %.2f)\n", legacy, tets,
                tets / legacy);
    CHECK(tets < 1.5 * legacy, "the tetrahedra cost %.1f ms per frame, shape matching %.1f", tets, legacy);
}
