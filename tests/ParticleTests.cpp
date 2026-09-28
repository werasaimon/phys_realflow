// Particles (liquid, soft bodies, cloth): rest density and its error, buoyancy of bodies
// against Archimedes, shape recovery of soft bodies, tearing and burning of cloth, the walls in
// the density (no corner jets), momentum conservation in collisions.
#include "TestRunner.h"
#include "Tests.h"

#include "scene/SceneGraph.h"

#include <string>

void testSoftBodyAndCloth() {
    const float dt = 1.0f / 180.0f;
    // 1) Curtain pinned at two corners: XPBD stretch constraints + long range attachments.
    {
        ParticleSystem s;
        s.reset(AABB({-1, 0, -1}, {1, 2, 1}));
        ClothMaterial inextensible;
        inextensible.tensileStiffness = 0.0f;
        inextensible.strengthWarp = inextensible.strengthWeft = 0.0f;
        s.addCloth({-0.3f, 1.5f, 0}, {0.6f, 0, 0}, {0, -0.6f, 0}, inextensible, 1 | 2, Vector3(1));
        for (int k = 0; k < 360; ++k) s.step(dt);
        float worst = 0;
        for (const DistanceConstraint& d : s.cloths()[0].constraints)
            if (d.restLength < 1.5f * s.params.clothSpacing * s.params.particleRadius) // stretch / shear
                worst = std::max(worst, length(s.positions()[d.a] - s.positions()[d.b]) / d.restLength - 1.0f);
        std::printf("  curtain %dx%d: max stretch %.2f%%\n", s.cloths()[0].width, s.cloths()[0].height, 100 * worst);
        CHECK(worst < 0.01f, "cloth stretches %.1f%%", 100 * worst);
    }
    // 2) Soft cube dropped on the floor: squashes on impact, springs back.
    {
        ParticleSystem s;
        s.reset(AABB({-1, 0, -1}, {1, 2, 1}));
        TriMesh cube = primitives::box(Vector3(0.12f));
        cube.translate({0, 0.8f, 0});
        int sb = s.addSoftBody(cube, 400.0f, 0.1f, Vector3(1));
        std::vector<Vector3> rest;
        for (int i : s.softBodies()[sb].particles) rest.push_back(s.positions()[i]);
        float peak = 0;
        for (int k = 0; k < 540; ++k) {
            s.step(dt);
            peak = std::max(peak, softShapeError(s, s.softBodies()[sb], rest));
        }
        // Resting on the floor it stays squashed by its own weight (a jelly this soft sags by a few
        // per cent); elastic means it springs back once the load is gone: gravity off for a second.
        const float sag = softShapeError(s, s.softBodies()[sb], rest);
        s.params.gravity = Vector3(0.0f);
        for (int k = 0; k < 180; ++k) s.step(dt);
        const float unloaded = softShapeError(s, s.softBodies()[sb], rest);
        std::printf("  soft cube (%zu particles, %zu clusters): deformation at impact %.1f%%, sag under its weight %.2f%%, unloaded %.2f%%\n",
                    s.softBodies()[sb].particles.size(), s.softBodies()[sb].clusters.size(), 100 * peak, 100 * sag, 100 * unloaded);
        CHECK(peak > 0.03f, "a soft cube must squash on impact (%.1f%%)", 100 * peak);
        CHECK(sag < 0.15f, "a soft cube must carry its own weight (%.1f%% sag)", 100 * sag);
        CHECK(unloaded < 0.02f, "a soft cube must spring back once unloaded (%.1f%% left)", 100 * unloaded);
    }
    // 3) Canvas trampoline pinned at its corners holds a foam cube.
    {
        ParticleSystem s;
        s.reset(AABB({-1, 0, -1}, {1, 2, 1}));
        ClothMaterial canvas;
        canvas.areaDensity = 1.5f;
        canvas.tensileStiffness = 0.0f;
        canvas.bendCompliance = 1e-4f;
        canvas.strengthWarp = canvas.strengthWeft = 0.0f;
        s.addCloth({-0.4f, 0.8f, -0.4f}, {0.8f, 0, 0}, {0, 0, 0.8f}, canvas, 15, Vector3(1));
        TriMesh cube = primitives::box(Vector3(0.1f));
        cube.translate({0, 1.0f, 0});
        int sb = s.addSoftBody(cube, 150.0f, 0.5f, Vector3(1));
        for (int k = 0; k < 360; ++k) s.step(dt);
        float ylo = 1e9f, vmax = 0;
        for (int i : s.softBodies()[sb].particles) {
            ylo = std::min(ylo, s.positions()[i].y);
            vmax = std::max(vmax, length(s.velocities()[i]));
        }
        std::printf("  foam cube on a canvas trampoline: cube bottom y %.3f (cloth plane 0.8), max |v| %.3f m/s\n", ylo, vmax);
        CHECK(ylo > 0.65f, "cube fell through the cloth (bottom y %f)", ylo);
        CHECK(vmax < 0.2f, "cube not at rest on the cloth (%f m/s)", vmax);
    }
    // 4) Buoyancy: a soft cube of density 500 floats in water (about half submerged).
    {
        ParticleSystem s;
        s.reset(AABB({-0.4f, 0, -0.3f}, {0.4f, 1.2f, 0.3f}));
        s.addBlock(AABB({-0.4f, 0, -0.3f}, {0.4f, 0.3f, 0.3f}));
        TriMesh cube = primitives::box(Vector3(0.09f));
        cube.translate({0, 0.7f, 0});
        int sb = s.addSoftBody(cube, 500.0f, 0.6f, Vector3(1));
        for (int k = 0; k < 720; ++k) s.step(dt);
        Vector3 c(0.0f);
        for (int i : s.softBodies()[sb].particles) c += s.positions()[i];
        c /= float(s.softBodies()[sb].particles.size());
        std::vector<float> ys;
        for (size_t i = 0; i < s.size(); ++i)
            if (s.phases()[i] == uint8_t(ParticlePhase::Fluid)) ys.push_back(s.positions()[i].y);
        std::sort(ys.begin(), ys.end());
        const float surface = ys[size_t(0.98 * ys.size())];
        std::printf("  soft cube (500 kg/m3) in water: centre y %.3f, water surface %.3f\n", c.y, surface);
        CHECK(std::fabs(c.y - surface) < 0.06f, "a body half as dense as water must float at the surface (%f vs %f)", c.y, surface);
    }
    // 5) Tearing (thread tension vs strength, cracks along the weave) and mouse grab.
    //    a) a cotton curtain on a rod does not tear under its own weight, nor when waved gently;
    //    b) pulled down hard at the bottom, it rips across the load (horizontal crack);
    //    c) two panels sewn together across a hanging sheet, heavy enough that the seam (0.3 of the
    //       strength) carries more than it holds and the fabric less: they part along the seam;
    //    d) a grabbed soft body follows; its skinned surface matches the model at rest.
    auto brokenThreads = [](const Cloth& c, int& warp, int& weft, int seamColumn) {
        int onSeam = 0;
        warp = weft = 0;
        for (const DistanceConstraint& d : c.constraints) {
            if (!d.broken || d.strength <= 0) continue;
            if (d.kind == DistanceConstraint::Warp) { ++warp; onSeam += d.x == seamColumn; }
            else ++weft;
        }
        return onSeam;
    };
    {
        ParticleSystem s;
        s.reset(AABB({-1, 0, -1}, {1, 2, 1}));
        s.addCloth({-0.3f, 1.5f, 0}, {0.6f, 0, 0}, {0, -0.6f, 0}, ClothMaterial(), 16, Vector3(1)); // on a rod
        for (int k = 0; k < 180; ++k) s.step(dt);
        const int hanging = s.cloths()[0].tornThreads;
        const Cloth& c = s.cloths()[0];
        CHECK(s.grab(s.positions()[c.particle(c.width / 2, c.height - 1)]), "no cloth particle grabbed");
        Vector3 t = s.grabTarget();
        float wavedLoad = 0; // the largest thread tension while waved, as a part of the thread's strength
        for (int k = 0; k < 60; ++k) {
            t += Vector3(0, 0, 0.001f);
            s.setGrabTarget(t);
            s.step(dt);
            for (const DistanceConstraint& d : c.constraints)
                if (!d.broken && d.strength > 0) wavedLoad = std::max(wavedLoad, threadTension(d, s.positions(), dt) / d.strength);
        }
        const int gentle = s.cloths()[0].tornThreads;
        std::printf("  curtain waved gently: threads loaded to %.0f%% of their strength at most, cloth cut into %d small steps\n",
                    100 * wavedLoad, c.smallSteps);
        for (int k = 0; k < 150; ++k) { t += Vector3(0, -0.01f, 0.004f); s.setGrabTarget(t); s.step(dt); }
        int warp, weft;
        brokenThreads(s.cloths()[0], warp, weft, -1);
        s.releaseGrab();
        for (int k = 0; k < 90; ++k) s.step(dt);
        bool finite = true;
        for (const Vector3& p : s.positions()) finite &= std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
        std::printf("  curtain on a rod: torn while hanging %d, while waved %d; pulled down: %d weft (vertical) + %d warp threads%s",
                    hanging, gentle, weft, warp, "\n");
        CHECK(hanging == 0 && gentle == 0, "cloth tore without being loaded (%d / %d threads)", hanging, gentle);
        CHECK(weft > 20 && weft > warp && finite, "a downward pull must rip across (weft %d, warp %d)", weft, warp);
    }
    {
        // The seam runs across the sheet 10 rows under the rod. The weight per metre of width at the
        // rod is 0.6 of the weft strength, at the seam ~0.45 (what hangs below it): more than the
        // seam's 0.3, less than the fabric's 1. The weight comes on over 1.5 s (no dynamic overshoot).
        ParticleSystem s;
        s.reset(AABB({-1, -1, -1}, {1, 2, 1}));
        ClothMaterial sewn;
        sewn.seamRows = {10};
        sewn.areaDensity = 0.6f * sewn.strengthWeft / (9.81f * 0.6f);
        s.addCloth({-0.3f, 1.5f, 0}, {0.6f, 0, 0}, {0, -0.6f, 0}, sewn, 16, Vector3(1)); // on a rod
        const Vector3 g = s.params.gravity;
        for (int k = 0; k < 360; ++k) {
            s.params.gravity = g * std::min(1.0f, float(k) / 270.0f);
            s.step(dt);
        }
        const Cloth& c = s.cloths()[0];
        int onSeam = 0, other = 0;
        for (const DistanceConstraint& d : c.constraints) {
            if (!d.broken || d.strength <= 0) continue;
            if (d.kind == DistanceConstraint::Weft && d.y == 10) ++onSeam;
            else ++other;
        }
        std::printf("  sewn panels hanging, the seam loaded to 1.5 x what it holds: %d of %d seam stitches gave, %d other threads%s",
                    onSeam, c.width, other, "\n");
        CHECK(onSeam == c.width, "the seam must open along its whole length (%d of %d)", onSeam, c.width);
    }
    {
        ParticleSystem q;
        q.reset(AABB({-1, 0, -1}, {1, 2, 1}));
        TriMesh ball = primitives::sphere(0.08f, 16, 8);
        ball.translate({0, 0.1f, 0});
        int b = q.addSoftBody(ball, 150.0f, 0.4f, Vector3(1));
        std::vector<Vector3> surf;
        q.softBodySurface(b, surf);
        float skinErr = 0;
        for (size_t v = 0; v < surf.size(); ++v) skinErr = std::max(skinErr, length(surf[v] - q.softBodies()[b].surface.positions[v]));
        for (int k = 0; k < 60; ++k) q.step(dt);
        auto centre = [&]() {
            Vector3 m(0.0f);
            for (int i : q.softBodies()[b].particles) m += q.positions()[i];
            return m / float(q.softBodies()[b].particles.size());
        };
        const float y0 = centre().y;
        q.grab(q.positions()[q.softBodies()[b].particles.back()]);
        Vector3 g = q.grabTarget();
        for (int k = 0; k < 120; ++k) { g += Vector3(0, 0.003f, 0); q.setGrabTarget(g); q.step(dt); }
        const float y1 = centre().y;
        std::printf("  soft body grab: centre rises %.3f m (target 0.36 m); skin at rest error %.1e m%s", y1 - y0, skinErr, "\n");
        CHECK(y1 - y0 > 0.25f, "a grabbed soft body must follow (rose %f m)", y1 - y0);
        CHECK(skinErr < 1e-4f, "skinned surface differs from the model at rest by %f", skinErr);
    }
    // 6) Momentum: two soft bodies collide head-on in zero gravity (contacts split corrections by mass).
    {
        ParticleSystem s;
        s.params.gravity = Vector3(0.0f);
        s.reset(AABB({-1, -1, -1}, {1, 1, 1}));
        TriMesh a = primitives::box(Vector3(0.08f)), b = primitives::box(Vector3(0.06f));
        a.translate({-0.3f, 0, 0});
        b.translate({0.3f, 0, 0});
        int ia = s.addSoftBody(a, 300.0f, 0.5f, Vector3(1), {2.0f, 0, 0});
        int ib = s.addSoftBody(b, 600.0f, 0.5f, Vector3(1), {-1.0f, 0, 0});
        auto momentum = [&]() {
            Vector3 P(0.0f);
            for (int k : {ia, ib})
                for (int i : s.softBodies()[k].particles) P += s.velocities()[i] * (1.0f / s.invMasses()[i]);
            return P;
        };
        Vector3 before = momentum();
        for (int k = 0; k < 120; ++k) s.step(dt); // collide, not yet at the walls
        Vector3 after = momentum();
        std::printf("  soft bodies colliding in zero g: momentum x %.4f -> %.4f kg m/s\n", before.x, after.x);
        CHECK(std::fabs(after.x - before.x) < 0.02f * std::fabs(before.x) + 1e-3f, "momentum not conserved: %f -> %f", before.x,
              after.x);
    }
}

// Thread tensions of a cloth now, per unit width [N/m] (tension of the thread / the width it
// stands for), warp and weft threads that are still whole; rows y < nearRows also go to `nearPin`.
static std::vector<float> threadTensionsPerWidth(const Cloth& c, const std::vector<Vector3>& x, float dt, int nearRows,
                                                 float& nearPin) {
    std::vector<float> out;
    for (const DistanceConstraint& d : c.constraints) {
        if (d.broken || (d.kind != DistanceConstraint::Warp && d.kind != DistanceConstraint::Weft)) continue;
        const float t = threadTension(d, x, dt) / threadWidth(c, d);
        out.push_back(t);
        if (d.y < nearRows) nearPin = std::max(nearPin, t);
    }
    return out;
}

struct SwingResult {
    int torn = 0;
    float nearPin = 0;   // largest thread tension in the pinned quarter while the sheet swings down [N/m]
    float peak = 0;      // largest thread tension anywhere, any time [N/m]
    float p99 = 0;       // largest 99th percentile of the thread tensions of a step [N/m]
    float edgeSpeed = 0; // fastest the free edge moved [m/s]
};

// A 1 x 1 m sheet pinned along one edge, let go level, swings down on it and back for 3 s.
static SwingResult swingSheet(const ClothMaterial& material) {
    ParticleSystem s;
    s.reset(AABB({-1.2f, -0.5f, -1.2f}, {1.2f, 2.5f, 1.2f}));
    s.addCloth({-0.5f, 2.0f, -0.5f}, {1, 0, 0}, {0, 0, 1}, material, 16, Vector3(1));
    const float dt = 1.0f / 180.0f;
    SwingResult r;
    for (int k = 0; k < 540; ++k) {
        s.step(dt);
        const Cloth& c = s.cloths()[0];
        float nearPin = 0;
        std::vector<float> t = threadTensionsPerWidth(c, s.positions(), dt, c.height / 4, nearPin);
        std::sort(t.begin(), t.end());
        r.peak = std::max(r.peak, t.back());
        r.p99 = std::max(r.p99, t[size_t(0.99 * double(t.size() - 1))]);
        if (k * dt < 0.5f) r.nearPin = std::max(r.nearPin, nearPin); // swinging down, before the fold reaches the edge
        for (int x = 0; x < c.width; ++x) r.edgeSpeed = std::max(r.edgeSpeed, length(s.velocities()[c.particle(x, c.height - 1)]));
    }
    r.torn = s.cloths()[0].tornThreads;
    return r;
}

// A cotton sheet (0.3 kg/m^2, threads 30 kN/m, strength 4000 / 3000 N/m) swings on its pinned edge.
// What the swing really asks of the threads, per unit width, with M = 0.3 kg/m the sheet's mass
// per metre of the pinned edge:
//  - near the pinned edge, while the sheet swings down: the pull of a pendulum, at most (at the
//    bottom) M g (1 + 2) = 8.8 N/m for all the mass at the free end (a stiff plate: 5/2 M g) -
//    energy M g L = M v^2 / 2 at the end, centripetal M v^2 / L = 2 M g, plus the weight;
//  - at the free edge, when the fold that runs down the falling sheet reaches it: the edge is
//    stopped like the tip of a whip, T ~ v sqrt(k mu) (a strip of density mu and stiffness k
//    stopped from speed v: the stress wave's impedance), a few hundred N/m at 5 m/s.
// Both are far below the strength: not one thread may break. Two sheets: the default one, which
// shears and bends like cotton (ClothMaterial: Kawabata's KES-F values), must read the pendulum's
// pull near the pin; one like card (the old defaults: shear 3000 N/m, bendCompliance 1e-3, a
// hundred times stiffer) reads more there - its bending links, pushed together where the sheet
// folds over the rod, press on the threads - and must still lose nothing. Before the threads were
// solved as whole lines (Cloth.cpp, solveThreadLine) and the contact passes started from the
// threads' force (ParticleSystem::step), the solver's lag stretched them 3-8 % instead of 0.03 %:
// ~1000 N/m near the pin, 0.99 of the strength at the bottom (the scene graph's 0.6 m sheet tore
// ~200 threads).
void testClothSwingNoFalseTears() {
    const ClothMaterial cotton;
    ClothMaterial card;
    card.shearStiffness = 3000.0f;
    card.bendCompliance = 1e-3f;
    const SwingResult soft = swingSheet(cotton), stiff = swingSheet(card);
    const float pendulum = 3 * cotton.areaDensity * 1.0f * 9.81f;
    const float whip = soft.edgeSpeed * std::sqrt(cotton.tensileStiffness * cotton.areaDensity);
    std::printf("  sheet 1 x 1 m swinging 3 s (strength %.0f / %.0f N/m); pendulum pull 3 M g = %.1f N/m, whip v sqrt(k mu) = "
                "%.0f N/m at the edge's %.1f m/s\n", cotton.strengthWarp, cotton.strengthWeft, pendulum, whip, soft.edgeSpeed);
    for (const SwingResult* r : {&soft, &stiff})
        std::printf("    %s: %d threads torn; near the pin while swinging down %.1f N/m; anywhere max %.0f N/m, 99th percentile %.0f N/m\n",
                    r == &soft ? "cotton (default)" : "card            ", r->torn, r->nearPin, r->peak, r->p99);
    CHECK(soft.torn == 0 && stiff.torn == 0, "a swinging sheet tore threads (%d, %d)", soft.torn, stiff.torn);
    CHECK(soft.nearPin < pendulum, "near the pin the threads read %.1f N/m, more than the pendulum's 3 M g = %.1f N/m", soft.nearPin,
          pendulum);
    CHECK(soft.peak < whip, "the threads read %.0f N/m, more than the whip at the free edge gives (%.0f N/m)", soft.peak, whip);
    CHECK(stiff.peak < 0.5f * cotton.strengthWeft, "card-like sheet: %.0f N/m, within a factor 2 of tearing", stiff.peak);
}

struct HangResult {
    float weight = 0;     // what hangs below the rod row [N]
    float carried = 0;    // the vertical pull of every link across the gap under the rod row [N]
    float threads = 0;    // the part of it the weft threads carry [N]
    float topThread = 0;  // the middle column's top weft thread, per unit width [N/m]
    int torn = 0, tornTopRow = 0;
};

// Statics, then strength. A sheet 0.6 m tall hangs from a rod (its top row pinned); its weight is
// brought on slowly - gravity ramps up over 1.5 s, so there is no dynamic overshoot - and held
// 0.5 s. `loadOverStrength` sets its density so that its weight per metre of width, mu g L, is that
// many times the weft strength.
static HangResult hangSheet(float loadOverStrength) {
    ParticleSystem s;
    s.reset(AABB({-1, -1, -1}, {1, 2, 1}));
    ClothMaterial m;
    const float L = 0.6f, g = 9.81f, dt = 1.0f / 180.0f;
    m.areaDensity = loadOverStrength * m.strengthWeft / (g * L);
    s.addCloth({-0.3f, 1.5f, 0}, {L, 0, 0}, {0, -L, 0}, m, 16, Vector3(1));
    for (int k = 0; k < 360; ++k) {
        s.params.gravity = Vector3(0, -g * std::min(1.0f, float(k) / 270.0f), 0);
        s.step(dt);
    }
    const Cloth& c = s.cloths()[0];
    const auto& x = s.positions();
    HangResult r;
    r.weight = float(c.width * (c.height - 1)) * m.areaDensity * c.particleArea * g;
    auto row = [&](int particle) { return (particle - c.firstParticle) / c.width; };
    for (const DistanceConstraint& d : c.constraints) {
        if (d.broken || (row(d.a) == 0) == (row(d.b) == 0)) continue; // only links across the gap under the rod
        const Vector3 dir = normalize(x[d.b] - x[d.a]);                   // from the rod row down
        const float up = threadTension(d, x, dt) * -dir.y;               // the upward pull on what hangs
        r.carried += up;
        if (d.kind == DistanceConstraint::Weft) r.threads += up;
        if (d.kind == DistanceConstraint::Weft && d.x == c.width / 2) r.topThread = threadTension(d, x, dt) / threadWidth(c, d);
    }
    r.torn = c.tornThreads;
    for (const DistanceConstraint& d : c.constraints) r.tornTopRow += d.kind == DistanceConstraint::Weft && d.y == 0 && d.broken;
    return r;
}

// Newton at the rod: the links across the gap under it hold up exactly the weight below (the
// threads most of it, the shear diagonals and bending links the rest), so the thread tension the
// tearing reads is the real one. Then the strength: loaded to 0.8 of it the sheet holds; loaded
// to 1.3 (the threads' share of it still above the strength) it tears - across the top, where the
// load is largest.
void testClothTearsAtStrength() {
    const HangResult light = hangSheet(0.8f), heavy = hangSheet(1.3f);
    const float strength = ClothMaterial().strengthWeft;
    std::printf("  sheet hanging with 0.8 x its strength: links under the rod carry %.1f N of its %.1f N weight (%+.2f%%), "
                "threads %.0f%% of it; middle top thread %.0f N/m; %d torn\n", light.carried, light.weight,
                100 * (light.carried / light.weight - 1), 100 * light.threads / light.carried, light.topThread, light.torn);
    std::printf("  with 1.3 x (threads' share %.0f N/m vs strength %.0f N/m): %d threads torn, %d of them in the top row\n",
                1.3f * strength * light.threads / light.carried, strength, heavy.torn, heavy.tornTopRow);
    CHECK(std::fabs(light.carried / light.weight - 1) < 0.02f, "the links carry %.2f N, the weight is %.2f N", light.carried, light.weight);
    CHECK(light.torn == 0, "a sheet loaded below its strength tore %d threads", light.torn);
    CHECK(heavy.tornTopRow > 0, "a sheet loaded above its strength must tear at the top (%d torn, %d at the top)", heavy.torn,
          heavy.tornTopRow);
}

void testSPH() {
    ParticleSystem s;
    s.params.particleRadius = 0.025f;
    AABB dom({0, 0, 0}, {0.6f, 0.8f, 0.4f});
    s.reset(dom);
    s.addBlock(AABB({0, 0, 0}, {0.6f, 0.4f, 0.4f}));
    size_t n = s.size();
    CHECK(n > 500, "particles %zu", n);
    for (int f = 0; f < 90; ++f)
        for (int k = 0; k < s.params.substeps; ++k) s.step(1.0f / 60 / s.params.substeps);
    float maxY = 0;
    bool finite = true;
    for (const Vector3& p : s.positions()) {
        maxY = std::max(maxY, p.y);
        finite &= std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z);
        finite &= dom.contains(p);
    }
    CHECK(finite, "non-finite or escaped particles");
    CHECK(s.averageDensityError() < 0.05f, "density error %f", s.averageDensityError());
    // The column should neither collapse (compression) nor expand much.
    CHECK(maxY > 0.33f && maxY < 0.5f, "free surface height %f", maxY);
    CHECK(s.maxSpeed() < 0.5f, "fluid at rest, max speed %f", s.maxSpeed());
}

void testFloating() {
    ParticleSystem s;
    RigidWorld w;
    s.params.particleRadius = 0.02f;
    AABB dom({-0.4f, 0, -0.3f}, {0.4f, 0.8f, 0.3f});
    w.setDomain(dom);
    int light = w.addBox({-0.15f, 0.6f, 0}, {0.08f, 0.08f, 0.08f}, Quaternion(), 400, Vector3(1));
    int heavy = w.addBox({0.15f, 0.6f, 0}, {0.06f, 0.06f, 0.06f}, Quaternion(), 3000, Vector3(1));
    s.setRigidWorld(&w);
    s.reset(dom);
    s.addBlock(AABB(dom.lo, {dom.hi.x, 0.35f, dom.hi.z}));
    for (int f = 0; f < 150; ++f)
        for (int k = 0; k < s.params.substeps; ++k) {
            float dt = 1.0f / 60 / s.params.substeps;
            w.step(dt);
            s.step(dt);
        }
    float yl = w.bodies()[light].pos.y, yh = w.bodies()[heavy].pos.y;
    CHECK(yl > 0.25f, "light box should float, y=%f", yl);
    CHECK(yh < 0.12f, "heavy box should sink, y=%f", yh);
}

void testLiquidWalls() {
    // Walls in the density (Koschier & Bender 2017): a particle right at a wall sees half its
    // kernel beyond it; a tank of water sloshing into its corners never clumps particles on the
    // corner edges (before: 10x rest density there, jets along the edges at the speed limit).
    ParticleSystem ps;
    ps.reset(AABB({0, 0, 0}, {1, 1, 1}));
    Vector3 g;
    const float atWall = ps.wallVolume({0.5f, 0.0f, 0.5f}, g), inside = ps.wallVolume({0.5f, 0.5f, 0.5f}, g);
    std::printf("  wall part of the kernel: %.4f at the wall, %.4f inside\n", atWall, inside);
    CHECK(std::fabs(atWall - 0.5f) < 1e-3f && inside == 0.0f, "wall volume wrong (%f, %f)", atWall, inside);

    Simulation sim;
    loadSample(sim, Preset::Water);
    const ParticleSystem& P = sim.particles;
    // The bug's signature is the density: particles lined up on a corner edge reached 10x rest
    // density (10 279 kg/m^3 measured before the fix). Speed is no test - the dam-break front
    // legitimately runs at ~2 sqrt(g H) = 4.6 m/s and spray falls from a metre at 4.4 m/s.
    float rhoMax = 0;
    for (int f = 1; f <= 240; ++f) {
        sim.stepFrame();
        if (f < 30) continue;
        for (size_t i = 0; i < P.size(); ++i)
            if (P.phases()[i] == uint8_t(ParticlePhase::Fluid)) rhoMax = std::max(rhoMax, P.densities()[i]);
    }
    std::printf("  water tank over 4 s: max density %.0f kg/m^3 (rest 1000; 10 279 in the corners before the fix)\n", rhoMax);
    CHECK(rhoMax < 2000.0f, "particles clump (density %f)", rhoMax);
}

void testLightBodyInWater() {
    // A beach ball (80 kg/m^3, as light as 12 water particles) hit by the dam-break wave: many
    // particles at once must not kick it (before the Gauss-Seidel body contacts: 155 m/s and
    // 9000 rad/s, from summing impulses computed against an immovable body), and it must float.
    Simulation sim;
    loadSample(sim, Preset::Water);
    const RigidBody& ball = sim.rigid.bodies()[2];
    float vMax = 0, wMax = 0, yEnd = 0;
    for (int f = 1; f <= 360; ++f) {
        sim.stepFrame();
        vMax = std::max(vMax, length(ball.vel));
        wMax = std::max(wMax, length(ball.angVel));
        if (f > 300) yEnd += ball.pos.y / 60.0f;
    }
    // Calm water level: pool 0.2 m + the column's 0.45 x 0.55 m spread over the 2 m tank.
    const float level = 0.2f + 0.45f * 0.55f / 2.0f;
    std::printf("  beach ball: max |v| %.2f m/s, max |w| %.1f rad/s; centre at the end %.3f m, water level ~%.3f m (radius %.2f)\n",
                vMax, wMax, yEnd, level, ball.radius());
    CHECK(vMax < 5.0f && wMax < 50.0f, "the light ball was kicked (v %f, w %f)", vMax, wMax);
    CHECK(yEnd > level - 0.3f * ball.radius(), "the beach ball must float (centre %f)", yEnd);
}

// Validation against the classic experiment: Martin & Moyce 1952, a water column a wide and 2a
// tall (n^2 = 2) collapsing on a dry floor. The front position Z = x / a against the
// dimensionless time T = t sqrt(2 g / a); the reference table is the digitised curve for
// a = 2.25 in that VOF, SPH and MPS papers reproduce (scratchpad/benchmarks.md, item 1), the late
// slope dZ/dT ~ 1.7 (Ritter's frictionless shallow-water theory gives 2). Grid and particle
// codes sit within 3-8 % of the front; SPH/PBF run slightly ahead (no gate, no floor friction).
// The front is the bulk of the floor layer (the last 1 cm bin holding at least 5 particles), not
// the few leading splash particles. The column is 20 particles wide and the substeps keep the
// PBF velocity clamp (0.5 h / dt = 6 m/s) far above the front speed of ~3 m/s.
//
// What this test found (2026-09-27): the front ran 35 % ahead of the experiment and a column left
// at rest grew 80 % in height - the artificial pressure s_corr was scaled in absolute units while
// lambda has the units of h^2, so 5 mm particles were puffed up 9 times harder than the 15 mm ones
// of the demo scenes. With s_corr = -k h^2 (W / W_dq)^4 the resting column keeps its height and
// the front is within 12 % (literature 3-8 %; the rest is particle resolution and the PBF density
// solver at 4 iterations, 7 % over-dense under its own weight). The threshold is set at what the
// code does now with the 8 % of the literature as the target.
void testDamBreakMartinMoyce() {
    const float a = 0.2f, g = 9.81f, depth = 0.1f;
    auto makeColumn = [&](ParticleSystem& s) {
        s.params.particleRadius = 0.005f; // spacing 1 cm: 20 particles across the column
        s.params.substeps = 10;           // dt = 1/600 s: the velocity clamp is 6 m/s
        s.params.gravity = {0, -g, 0};
        s.reset(AABB({0, 0, 0}, {6 * a, 3 * a, depth}));
        s.addBlock(AABB({0, 0, 0}, {a, 2 * a, depth}));
    };
    // 1) The same column held by a wall at rest (a domain as wide as the column): how much it grows.
    float restGrowth = 0;
    {
        ParticleSystem s;
        s.params.particleRadius = 0.005f;
        s.params.substeps = 10;
        s.params.gravity = {0, -g, 0};
        s.reset(AABB({0, 0, 0}, {a, 3 * a, depth}));
        s.addBlock(AABB({0, 0, 0}, {a, 2 * a, depth}));
        for (int k = 0; k < 300; ++k) s.step(1.0f / 600); // 0.5 s
        float top = 0;
        for (const Vector3& p : s.positions()) top = std::max(top, p.y + s.params.particleRadius);
        restGrowth = top / (2 * a) - 1.0f;
    }
    // 2) The dam break.
    ParticleSystem s;
    makeColumn(s);
    const float r = s.params.particleRadius, dt = 1.0f / 60 / s.params.substeps, scale = std::sqrt(2 * g / a);
    std::printf("  dam break: %zu particles, column %.0f x %.0f particles; the column at rest grows by %.1f%% in height\n", s.size(), a / (2 * r),
                2 * a / (2 * r), 100 * restGrowth);
    // Reference (T, Z), n^2 = 2, digitised.
    const float refT[] = {0.41f, 0.84f, 1.19f, 1.43f, 1.63f, 1.83f, 2.02f, 2.20f, 2.37f, 2.53f, 2.69f, 2.85f, 3.00f};
    const float refZ[] = {1.11f, 1.22f, 1.44f, 1.67f, 1.89f, 2.11f, 2.33f, 2.56f, 2.78f, 3.00f, 3.22f, 3.44f, 3.67f};
    const int nRef = int(sizeof(refT) / sizeof(refT[0]));
    auto reference = [&](float T) { // linear interpolation, Z = 1 at T = 0
        if (T <= refT[0]) return 1.0f + (refZ[0] - 1.0f) * T / refT[0];
        for (int i = 1; i < nRef; ++i)
            if (T <= refT[i]) return refZ[i - 1] + (refZ[i] - refZ[i - 1]) * (T - refT[i - 1]) / (refT[i] - refT[i - 1]);
        return refZ[nRef - 1];
    };
    // The bulk front along the floor: the last 1 cm bin along x with at least 5 particles of the
    // floor layer (leading splash particles do not count); the column height at the back wall.
    auto measure = [&](float& Zfront, float& Hback) {
        std::vector<int> bins(int(6 * a / 0.01f) + 1, 0);
        float back = 0;
        for (const Vector3& p : s.positions()) {
            if (p.y < 2.5f * r) ++bins[std::min(int((p.x + r) / 0.01f), int(bins.size()) - 1)];
            if (p.x < 3 * r) back = std::max(back, p.y + r);
        }
        Zfront = 0;
        for (size_t b = 0; b < bins.size(); ++b)
            if (bins[b] >= 5) Zfront = 0.01f * float(b + 1) / a;
        Hback = back / (2 * a);
    };
    std::vector<std::pair<float, float>> series; // (T, Z)
    double devSum = 0;
    int devCount = 0;
    float hMax = 0;
    std::printf("      T    Z(sim)  Z(ref)  H(back wall)\n");
    for (int step = 0, next = 0; ; ++step) {
        const float t = step * dt, T = t * scale;
        if (T > 3.05f) break;
        if (step == next) { // every ~0.1 in T
            float Z, H;
            measure(Z, H);
            hMax = std::max(hMax, H);
            const float Zr = reference(T);
            series.push_back({T, Z});
            if (T >= 0.5f && T <= 3.0f) { devSum += std::fabs(Z - Zr) / Zr; ++devCount; }
            if (series.size() % 3 == 1) std::printf("    %5.2f  %6.2f  %6.2f  %6.2f\n", T, Z, Zr, H);
            next += 6;
        }
        s.step(dt);
    }
    // Late slope dZ/dT over T in [1.5, 3]: a least-squares line.
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int n = 0;
    for (auto [T, Z] : series)
        if (T >= 1.5f && T <= 3.0f) { sx += T; sy += Z; sxx += T * T; sxy += T * Z; ++n; }
    const float slope = n >= 2 ? float((n * sxy - sx * sy) / (n * sxx - sx * sx)) : 0.0f;
    const float deviation = devCount ? float(devSum / devCount) : 1.0f;
    std::printf("  dam break front: mean |Z - Z_ref| / Z_ref = %.1f%% over T in [0.5, 3] (literature 3-8%%), late slope dZ/dT = %.2f "
                "(experiment ~1.7, Ritter 2), back wall rises to %.2f of the height (experiment: falls), max speed %.2f m/s\n",
                100 * deviation, slope, hMax, s.maxSpeed());
    CHECK(deviation < 0.15f, "dam break front off the Martin & Moyce curve by %.1f%% (today 12%%; target 8%%)", 100 * deviation);
    CHECK(slope > 0.75f * 1.7f && slope < 1.25f * 1.7f, "late front slope %f vs ~1.7", slope);
    CHECK(s.maxSpeed() < 0.9f * 0.5f * 4 * r / dt, "the velocity clamp limits the front: %f m/s", s.maxSpeed());
}

// Many bodies in a pool: every particle looks for the bodies it may touch. Before the world tree
// it tested all 150 bodies for each of 30 000 particles in every pass; with the tree it asks for
// the few whose boxes overlap its own. Same contacts, less work.
namespace {
struct ManyBodiesScene : Scene {
    void configure(Simulation& sim) override { sim.particles.params.particleRadius = 0.0135f; } // ~30 000 particles
    void build(Simulation& sim) override {
        const AABB tank({-1.0f, 0.0f, -0.5f}, {1.0f, 1.0f, 0.5f});
        sim.useLiquidTank(tank);
        uint32_t seed = 3;
        auto rnd = [&] { seed = seed * 1664525u + 1013904223u; return (seed >> 8) * (1.0f / 16777216.0f); };
        for (int i = 0; i < 150; ++i) {
            const Vector3 p(-0.9f + 1.8f * rnd(), 0.5f + 0.4f * rnd(), -0.4f + 0.8f * rnd());
            const float s = 0.025f + 0.015f * rnd(); // half-size 2.5-4 cm: 5-8 cm bodies
            if (i % 2 == 0) sim.rigid.addBox(p, Vector3(s), Quaternion::fromAxisAngle({rnd(), rnd(), rnd() + 0.1f}, 6.28f * rnd()), 400.0f, Vector3(1));
            else sim.rigid.addSphere(p, s, 600.0f, Vector3(1));
        }
        sim.particles.addBlock(AABB(tank.lo, {tank.hi.x, 0.3f, tank.hi.z})); // the pool, 30 cm deep
    }
};
} // namespace

void testParticlesManyBodies() {
    Simulation sim;
    sim.load(std::make_unique<ManyBodiesScene>());
    std::printf("  many bodies: %zu particles, %zu bodies\n", sim.particles.size(), sim.rigid.bodies().size());
    double stepMs = 0;
    int frames = 0, contacts = 0;
    for (int f = 0; f < 60; ++f) { // 1 s: the bodies fall in and float or sink
        sim.stepFrame();
        const Probe::Snapshot s = Probe::snapshot();
        stepMs += s.value("frame/step ms");
        contacts = std::max(contacts, int(s.value("particles/body contacts")));
        ++frames;
    }
    stepMs /= frames;
    // Measured on the build machine: brute force over 150 bodies per particle 220 ms per frame;
    // with the world tree 69 ms (451 770 queries per frame, 15 passes x 30 118 particles), the
    // same 1566 contacts and the same numbers in every other particle test.
    std::printf("  many bodies: %.1f ms per frame, up to %d particle-body contacts in a pass, tree queries %.0f per frame\n",
                stepMs, contacts, Probe::snapshot().value("rigid/tree queries"));
    CHECK(std::isfinite(stepMs) && stepMs > 0, "no timing");
    CHECK(contacts > 100, "the bodies must touch the water (%d contacts)", contacts);
    CHECK(stepMs < 140.0f, "a frame of 30 000 particles and 150 bodies costs %.1f ms (tree: 69, brute force: 220)", stepMs);
}

// Removing one group of particles (the editor's meta-objects: a soft cube turned into something
// else): a pool of liquid, a soft cube in it and a cloth hanging on a rod. The cube's group goes
// mid-run, then the liquid; the cloth goes on hanging from the same rod, and a new soft body
// added afterwards works.
static bool allFinite(const ParticleSystem& s) {
    for (const Vector3& x : s.positions())
        if (!std::isfinite(x.x + x.y + x.z)) return false;
    return true;
}

void testRemoveParticleGroup() {
    const float dt = 1.0f / 180.0f;
    ParticleSystem s;
    s.reset(AABB({-0.5f, 0, -0.3f}, {0.5f, 1.2f, 0.3f}));
    const int liquid = s.addBlock(AABB({-0.5f, 0, -0.3f}, {0.5f, 0.15f, 0.3f}));
    TriMesh cube = primitives::box(Vector3(0.06f));
    cube.translate({-0.2f, 0.4f, 0});
    const int soft = s.softBodyGroup(s.addSoftBody(cube, 400.0f, 0.3f, Vector3(1)));
    const int cloth = s.clothGroup(s.addCloth({0.05f, 1.0f, -0.2f}, {0.4f, 0, 0}, {0, -0.4f, 0}, ClothMaterial(), 16, Vector3(1)));
    for (int k = 0; k < 90; ++k) s.step(dt);
    const size_t before = s.size(), softCount = s.groupSize(soft), liquidCount = s.groupSize(liquid);
    std::vector<Vector3> rod; // the pinned particles of the cloth
    const Cloth& c0 = s.cloths()[0];
    for (int x = 0; x < c0.width; ++x) rod.push_back(s.positions()[size_t(c0.particle(x, 0))]);
    s.removeGroup(soft);
    for (int k = 0; k < 90; ++k) s.step(dt);
    const Cloth& c1 = s.cloths()[0];
    float rodShift = 0;
    for (int x = 0; x < c1.width; ++x) rodShift = std::max(rodShift, length(s.positions()[size_t(c1.particle(x, 0))] - rod[size_t(x)]));
    std::printf("  remove the soft cube: %zu -> %zu particles (its %zu), soft bodies %zu, cloth rod moved %.1e m, finite %d\n", before,
                s.size(), softCount, s.softBodies().size(), rodShift, int(allFinite(s)));
    CHECK(s.size() == before - softCount, "%zu particles left, expected %zu", s.size(), before - softCount);
    CHECK(s.softBodies().empty() && s.cloths().size() == 1 && s.cloths()[0].group == cloth, "the cube must be gone, the cloth kept");
    CHECK(rodShift < 1e-6f, "the cloth's rod moved %e m", rodShift);
    CHECK(allFinite(s), "NaN after removing the soft cube");
    s.removeGroup(liquid);
    for (int k = 0; k < 90; ++k) s.step(dt);
    std::printf("  remove the liquid: %zu particles left (fluid %zu), cloth %d x %d\n", s.size(), s.fluidCount(), s.cloths()[0].width,
                s.cloths()[0].height);
    CHECK(s.fluidCount() == 0 && s.size() == before - softCount - liquidCount, "the liquid must be gone (%zu fluid)", s.fluidCount());
    CHECK(allFinite(s), "NaN after removing the liquid");
    TriMesh ball = primitives::sphere(0.05f, 12, 6);
    ball.translate({-0.2f, 0.3f, 0});
    const int body = s.addSoftBody(ball, 400.0f, 0.3f, Vector3(1));
    for (int k = 0; k < 180; ++k) s.step(dt);
    float lowest = 1e9f;
    for (int i : s.softBodies()[size_t(body)].particles) lowest = std::min(lowest, s.positions()[size_t(i)].y);
    std::printf("  a new soft ball afterwards: %zu particles, lowest at y %.3f (floor 0)\n", s.softBodies()[size_t(body)].particles.size(),
                lowest);
    CHECK(allFinite(s) && lowest > -0.001f && lowest < 0.05f, "the new soft ball must fall to the floor (lowest %f)", lowest);
}

// ------------------------------------------------------------------------------------------------
//  Soft bodies against each other
// ------------------------------------------------------------------------------------------------

namespace {

// The user's scene, as the editor saves it: six soft cylinders 0.3 m wide and tall, stacked on a
// rigid floor with 2 mm between them.
std::string softBarrelScene() {
    std::string s = "world gravity 0 -9.81 0 size 4 4 4 gas 0 magneticGas 0\n"
                    "entity \"Пол\"\n  object id 1 visible 1 locked 1\n"
                    "  shape plane size 3 0.02 3 position 0 0.01 0 rotation 0 0 0 color 0.6 0.62 0.66\n"
                    "  rigid density 500 friction 0.5 restitution 0.2 fixed 1 velocity 0 0 0 spin 0 0 0\nend\n";
    for (int k = 0; k < 6; ++k) {
        char line[256];
        std::snprintf(line, sizeof line,
                      "entity \"Бочка %d\"\n  object id %d visible 1 locked 0\n"
                      "  shape cylinder size 0.3 0.3 0.3 position 0 %.4f 0 rotation 0 0 0 color 0.5 0.5 0.5\n"
                      "  soft density 500 stiffness 0.5\nend\n",
                      k + 1, k + 2, 0.17f + 0.302f * float(k));
        s += line;
    }
    return s;
}

// The centre of mass of a soft body's particles.
Vector3 softCentre(const ParticleSystem& s, const SoftBody& b) {
    Vector3 c(0.0f);
    for (int i : b.particles) c += s.positions()[size_t(i)];
    return c / float(b.particles.size());
}

// How deep particles of DIFFERENT soft bodies sit in each other, in particle spacings d0: the
// largest d0 - |xi - xj| over all such pairs (0: no two touch closer than their diameter), and how
// many pairs overlap by more than `deeper` spacings.
struct SoftOverlap {
    float deepest = 0;
    int pairs = 0;
};
SoftOverlap softOverlap(const ParticleSystem& s, float deeper) {
    std::vector<int> bodyOf(s.size(), -1);
    for (size_t b = 0; b < s.softBodies().size(); ++b)
        for (int i : s.softBodies()[b].particles) bodyOf[size_t(i)] = int(b);
    const float d0 = s.spacing();
    SoftOverlap o;
    const auto& x = s.positions();
    for (size_t i = 0; i < x.size(); ++i) {
        if (bodyOf[i] < 0) continue;
        for (size_t j = i + 1; j < x.size(); ++j) {
            if (bodyOf[j] < 0 || bodyOf[j] == bodyOf[i]) continue;
            const float d2 = length2(x[i] - x[j]);
            if (d2 >= d0 * d0) continue;
            const float depth = (d0 - std::sqrt(d2)) / d0;
            o.deepest = std::max(o.deepest, depth);
            o.pairs += depth > deeper;
        }
    }
    return o;
}

// While the barrels stand one on the other (every centre within 1 cm of the bottom one's axis):
// how far neighbours cross [m] over the middle of their faces (within 8 cm of the axis) - the drawn
// skins (`skin`: top of each barrel's skin against the bottom of the next one's) and the particle
// clouds (`particles`: top particle of each against the bottom particle of the next, plus one
// spacing, so 0 = the layers just touch). Negative: a gap. False when the stack no longer stands.
bool stackCrossing(const ParticleSystem& s, float& skin, float& particles) {
    const auto& bodies = s.softBodies();
    const Vector3 axis = softCentre(s, bodies[0]);
    for (const SoftBody& b : bodies) {
        const Vector3 c = softCentre(s, b);
        if (sqr(c.x - axis.x) + sqr(c.z - axis.z) > sqr(0.01f)) return false;
    }
    std::vector<float> top(bodies.size(), -kInf), bottom(bodies.size(), kInf);
    std::vector<Vector3> surface;
    for (size_t b = 0; b < bodies.size(); ++b) {
        s.softBodySurface(b, surface);
        for (const Vector3& v : surface) {
            if (sqr(v.x - axis.x) + sqr(v.z - axis.z) > sqr(0.08f)) continue;
            top[b] = std::max(top[b], v.y);
            bottom[b] = std::min(bottom[b], v.y);
        }
    }
    skin = particles = -kInf;
    for (size_t b = 0; b + 1 < bodies.size(); ++b) skin = std::max(skin, top[b] - bottom[b + 1]);
    for (size_t b = 0; b + 1 < bodies.size(); ++b) {
        float hi = -kInf, lo = kInf;
        for (int i : bodies[b].particles) {
            const Vector3& x = s.positions()[size_t(i)];
            if (sqr(x.x - axis.x) + sqr(x.z - axis.z) <= sqr(0.08f)) hi = std::max(hi, x.y);
        }
        for (int i : bodies[b + 1].particles) {
            const Vector3& x = s.positions()[size_t(i)];
            if (sqr(x.x - axis.x) + sqr(x.z - axis.z) <= sqr(0.08f)) lo = std::min(lo, x.y);
        }
        particles = std::max(particles, hi - lo + s.spacing());
    }
    return true;
}

} // namespace

// Six soft barrels stacked on a rigid floor (the scene a user built): they may squash, bounce and
// topple - a column of six soft barrels on a narrow base is unstable - but two barrels must never
// sink into each other, and never stay stuck. Before FleX's stiff-stack mass scaling
// (ParticleParams::stackMassScaling) they did: within half a second particles of neighbouring
// barrels overlapped by 0.87 d0, the barrels merged into one sausage and stayed stuck.
void testSoftStackNoOverlap() {
    SceneGraph g;
    std::string error;
    CHECK(g.load(softBarrelScene(), error), "the barrel scene does not load: %s", error.c_str());
    Simulation sim;
    sim.load(std::make_unique<GraphScene>(g));
    const float d0 = sim.particles.spacing();
    std::vector<float> start;
    for (const SoftBody& b : sim.particles.softBodies()) start.push_back(softCentre(sim.particles, b).y);
    float stacked = 0, anyTime = 0, skin = -kInf, layers = -kInf, rise = 0;
    int stuckPairs = 0, straightFrames = 0;
    for (int frame = 1; frame <= 180; ++frame) {
        sim.stepFrame();
        for (size_t b = 0; b < start.size(); ++b)
            rise = std::max(rise, softCentre(sim.particles, sim.particles.softBodies()[b]).y - start[b]);
        if (frame % 5 != 0) continue;
        const SoftOverlap o = softOverlap(sim.particles, 0.05f);
        float skinCross, layerCross;
        if (stackCrossing(sim.particles, skinCross, layerCross)) {
            stacked = std::max(stacked, o.deepest);
            skin = std::max(skin, skinCross / d0);
            layers = std::max(layers, layerCross / d0);
            ++straightFrames;
        }
        anyTime = std::max(anyTime, o.deepest);
        if (frame > 150) stuckPairs = std::max(stuckPairs, o.pairs);
    }
    std::printf("  six soft barrels: deepest overlap of two barrels %.2f d0 while stacked (%d checks), %.2f d0 at any time; "
                "particle layers cross by %.2f d0 and skins by %.2f d0 while stacked; after 2.5 s %d pairs deeper than 0.05 d0; highest rise of a barrel %.0f mm\n",
                stacked, straightFrames, anyTime, layers, skin, stuckPairs, 1000 * rise);
    // Limits, measured and explained (docs/03, "Мягкие тела друг на друге"): pairs of particles
    // of two barrels overlap by 0.10 d0 while the stack stands and 0.17 d0 in the blows of its
    // fall (the iterations' residual under 53 kg); the squashed faces spread and the neighbour's
    // bottom layer settles into the widened gaps by 1.4 d0, the skin following it (1.1 d0) - it
    // was 10 d0, the barrels merged, before the mass scaling. A rise of more than 5 cm would be
    // energy made by the scaling (k = 1.5 lifted a barrel 0.4 m, k = 2 by 1.1 m).
    CHECK(straightFrames >= 3, "the stack fell before it could be measured (%d checks)", straightFrames);
    CHECK(stacked < 0.15f && anyTime < 0.25f, "barrels sink into each other: %.2f d0 stacked, %.2f d0 falling", stacked, anyTime);
    CHECK(layers < 2.0f && skin < 2.0f, "the faces nest too deep: layers %.2f d0, skins %.2f d0", layers, skin);
    CHECK(stuckPairs == 0, "%d pairs of particles stay stuck in another barrel", stuckPairs);
    CHECK(rise < 0.05f, "a barrel rose %.0f mm above its start: energy from nowhere", 1000 * rise);
}

// Two soft cubes thrown at each other in weightless space bounce apart - not a pair of lattices
// that slid into each other and locked - and a rigid box resting on a soft cube stays on top of it.
void testSoftPressedApartAndBox() {
    const float dt = 1.0f / 180.0f;
    {
        ParticleSystem s;
        s.params.gravity = Vector3(0.0f);
        s.reset(AABB({-1, -1, -1}, {1, 1, 1}));
        TriMesh a = primitives::box(Vector3(0.12f)), b = primitives::box(Vector3(0.12f));
        a.translate({-0.12f, 0, 0});
        b.translate({0.12f, 0, 0});
        s.addSoftBody(a, 400.0f, 0.3f, Vector3(1), {1.0f, 0, 0});
        s.addSoftBody(b, 400.0f, 0.3f, Vector3(1), {-1.0f, 0, 0});
        float deepest = 0;
        for (int k = 0; k < 180; ++k) {
            s.step(dt);
            if (k % 10 == 0) deepest = std::max(deepest, softOverlap(s, 0.05f).deepest);
        }
        const Vector3 ca = softCentre(s, s.softBodies()[0]), cb = softCentre(s, s.softBodies()[1]);
        Vector3 va(0.0f), vb(0.0f);
        for (int i : s.softBodies()[0].particles) va += s.velocities()[size_t(i)];
        for (int i : s.softBodies()[1].particles) vb += s.velocities()[size_t(i)];
        va /= float(s.softBodies()[0].particles.size());
        vb /= float(s.softBodies()[1].particles.size());
        const SoftOverlap end = softOverlap(s, 0.05f);
        std::printf("  two soft cubes thrown together at 2 m/s: deepest overlap %.2f d0, after 1 s centres %.3f m apart (size 0.120), "
                    "separating at %.2f m/s, %d pairs still deeper than 0.05 d0\n",
                    deepest, cb.x - ca.x, vb.x - va.x, end.pairs);
        CHECK(deepest < 0.1f && end.pairs == 0, "the cubes sank into each other (%.2f d0, %d pairs left)", deepest, end.pairs);
        CHECK(cb.x - ca.x > 0.12f && vb.x - va.x > 0.1f, "the cubes did not bounce apart (%.3f m apart, %.2f m/s)", cb.x - ca.x, vb.x - va.x);
    }
    {
        ParticleSystem s;
        RigidWorld w;
        const AABB domain({-0.5f, 0, -0.5f}, {0.5f, 1, 0.5f});
        w.setDomain(domain);
        s.setRigidWorld(&w);
        s.reset(domain);
        TriMesh jelly = primitives::box(Vector3(0.16f));
        jelly.translate({0, 0.08f, 0});
        s.addSoftBody(jelly, 400.0f, 0.3f, Vector3(1));
        const int box = w.addBox({0, 0.22f, 0}, Vector3(0.05f), Quaternion(), 500.0f, Vector3(1));
        for (int k = 0; k < 360; ++k) {
            w.step(dt);
            s.step(dt);
        }
        float top = -kInf; // the jelly's surface under the box
        for (int i : s.softBodies()[0].particles) {
            const Vector3& x = s.positions()[size_t(i)];
            if (std::fabs(x.x - w.bodies()[size_t(box)].pos.x) < 0.04f && std::fabs(x.z - w.bodies()[size_t(box)].pos.z) < 0.04f &&
                x.y < w.bodies()[size_t(box)].pos.y)
                top = std::max(top, x.y);
        }
        const float gap = (w.bodies()[size_t(box)].pos.y - 0.05f) - (top + s.params.particleRadius);
        std::printf("  rigid box on a soft cube after 2 s: box bottom %.1f mm above the jelly's surface, box at y %.3f m\n",
                    1000 * gap, w.bodies()[size_t(box)].pos.y);
        // The box rests on the particles' surface (their centres + one radius) to within half a radius.
        CHECK(std::isfinite(gap) && gap > -0.5f * s.params.particleRadius && gap < 0.01f, "the box is not resting on the jelly (%.1f mm)", 1000 * gap);
    }
}
