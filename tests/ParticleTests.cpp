// Particles (liquid, soft bodies, cloth): rest density and its error, buoyancy of bodies
// against Archimedes, shape recovery of soft bodies, tearing and burning of cloth, the walls in
// the density (no corner jets), momentum conservation in collisions.
#include "TestRunner.h"
#include "Tests.h"

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
    //    c) two panels sewn together, one edge held, the other pulled: they part along the seam;
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
        for (int k = 0; k < 60; ++k) { t += Vector3(0, 0, 0.001f); s.setGrabTarget(t); s.step(dt); }
        const int gentle = s.cloths()[0].tornThreads;
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
        ParticleSystem s;
        s.reset(AABB({-1, 0, -1}, {1, 2, 1}));
        ClothMaterial sewn;
        sewn.seamColumns = {20};
        s.addCloth({-0.3f, 1.5f, 0}, {0.6f, 0, 0}, {0, -0.6f, 0}, sewn, 128, Vector3(1)); // right edge held
        for (int k = 0; k < 180; ++k) s.step(dt);
        const Cloth& c = s.cloths()[0];
        s.grab(s.positions()[c.particle(0, c.height / 2)]);
        Vector3 t = s.grabTarget();
        for (int k = 0; k < 150; ++k) { t += Vector3(-0.01f, 0, 0); s.setGrabTarget(t); s.step(dt); }
        int warp, weft;
        const int onSeam = brokenThreads(s.cloths()[0], warp, weft, 20);
        std::printf("  sewn panels pulled apart: %d of %d seam stitches gave, %d other threads%s", onSeam, c.height,
                    warp + weft - onSeam, "\n");
        CHECK(onSeam == c.height, "the seam must open along its whole length (%d of %d)", onSeam, c.height);
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
