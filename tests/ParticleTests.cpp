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
    sim.loadPreset(Preset::Water);
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
    sim.loadPreset(Preset::Water);
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
