// RenderSnapshot of Simulation: everything the viewer draws and the panels show, extracted
// from the solvers - slices, velocity arrows, streamlines, magnetic field lines, the smoke
// volume, the bodies, the particles, the readouts and the plots.
#include "scene/Simulation.h"

#include "core/Format.h"
#include "core/Parallel.h"

namespace rf {

// ---------------------------------------------------------------------------
// Snapshot extraction
// ---------------------------------------------------------------------------
static const char* fieldLabel(GridField f) {
    switch (f) {
    case GridField::Speed: return "|V|, м/с";
    case GridField::Pressure: return "p, Па";
    case GridField::PressureCoeff: return "Cp";
    case GridField::Vorticity: return "|ω|, 1/с";
    case GridField::Smoke: return "дым";
    case GridField::Temperature: return "T (отн.)";
    case GridField::VelocityX: return "Vx, м/с";
    case GridField::VelocityY: return "Vy, м/с";
    case GridField::MagneticFlux: return "|B|, Тл";
    case GridField::CurrentDensity: return "|J|, А/м²";
    }
    return "";
}

void Simulation::extractVectors(RenderSnapshot& s) const {
    const int n[3] = {grid.nx(), grid.ny(), grid.nz()};
    const int axis = clampv(vis.sliceAxis, 0, 2);
    const int layer = sliceLayer();
    const int st = std::max(1, vis.vectorStride);
    const bool plane = vis.vectorDisplay == 1;
    const size_t kMaxArrows = 200000;
    const float dx = grid.dx();
    const Vector3 o = grid.origin();
    float vmax = 0;
    for (int k = 0; k < n[2]; ++k)
        for (int j = 0; j < n[1]; ++j)
            for (int i = 0; i < n[0]; ++i) {
                int c[3] = {i, j, k};
                bool keep = true;
                for (int a = 0; a < 3; ++a) {
                    if (plane && a == axis) keep &= c[a] == layer;
                    else keep &= c[a] % st == st / 2;
                }
                if (!keep || grid.solid(i, j, k)) continue;
                if (vis.vectorsWhereSmoke && grid.smoke().at(i, j, k) < 0.02f) continue;
                if (s.arrowPos.size() >= kMaxArrows) break;
                Vector3 v = grid.cellVelocity(i, j, k);
                s.arrowPos.push_back(o + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * dx);
                s.arrowVel.push_back(v);
                vmax = std::max(vmax, length(v));
            }
    s.arrowMax = std::max(vmax, 1e-6f);
    s.arrowLength = 0.9f * dx * st * vis.vectorScale;
}

void Simulation::extractSlice(RenderSnapshot& s) const {
    const int nx = grid.nx(), ny = grid.ny(), nz = grid.nz();
    const float dx = grid.dx();
    const Vector3 o = grid.origin();
    int axis = clampv(vis.sliceAxis, 0, 2);
    int n[3] = {nx, ny, nz};
    int layer = sliceLayer();
    int ua = axis == 0 ? 2 : 0;         // horizontal image axis
    int va = axis == 1 ? 2 : 1;         // vertical image axis
    s.sliceW = n[ua];
    s.sliceH = n[va];
    s.slice.assign(size_t(s.sliceW) * s.sliceH, 0.0f);
    s.sliceSolid.assign(size_t(s.sliceW) * s.sliceH, 0);
    Vector3 corner = o;
    corner[axis] += (layer + 0.5f) * dx;
    Vector3 U(0.0f), V(0.0f);
    U[ua] = n[ua] * dx;
    V[va] = n[va] * dx;
    s.sliceCorner = corner;
    s.sliceU = U;
    s.sliceV = V;
    parallelFor(s.sliceH, [&](int b) {
        for (int a = 0; a < s.sliceW; ++a) {
            int c[3];
            c[axis] = layer;
            c[ua] = a;
            c[va] = b;
            size_t id = size_t(b) * s.sliceW + a;
            if (grid.solid(c[0], c[1], c[2])) {
                s.sliceSolid[id] = 1;
                continue;
            }
            s.slice[id] = grid.cellValue(vis.sliceField, c[0], c[1], c[2]);
        }
    }, 1);
    s.hasSlice = true;
}

// Magnetic field lines: traced both ways from seeds on a shell around the obstacle (the magnet)
// or on a lattice through the domain, along B / |B| (midpoint rule), until they leave the
// domain, enter a solid or the field vanishes.
void Simulation::computeFieldLines(RenderSnapshot& s) const {
    std::vector<Vector3> seeds;
    fieldLineSeeds(seeds);
    s.fieldLines.clear();
    s.fieldLineStrength.clear();
    float bmax = 1e-12f;
    for (const Vector3& seed : seeds) {
        std::vector<Vector3> line;
        std::vector<float> strength;
        traceFieldLine(seed, line, strength);
        if (line.size() < 2) continue;
        for (float B : strength) bmax = std::max(bmax, B);
        s.fieldLines.push_back(std::move(line));
        s.fieldLineStrength.push_back(std::move(strength));
    }
    s.fieldLineMax = bmax;
}

// Where the lines start: the scene's own seeds if it has some (a tokamak's poloidal plane), else
// a shell around the obstacle (the magnet), else a lattice through the domain.
void Simulation::fieldLineSeeds(std::vector<Vector3>& seeds) const {
    const float dx = grid.dx();
    const AABB dom = grid.domain();
    if (scene_) scene_->fieldLineSeeds(*this, seeds); // the scene knows where its field is interesting
    if (!seeds.empty()) {
    } else if (grid.hasObstacle() && obstacleMesh_ && !obstacleMesh_->empty()) {
        const AABB ob = obstacleMesh_->bounds();
        const Vector3 c = ob.center();
        const float r = 0.5f * maxComp(ob.extent()) + 2.0f * dx;
        for (int a = 0; a < 12; ++a)
            for (int b = 1; b <= 4; ++b) {
                const float phi = 2 * kPi * a / 12, theta = kPi * b / 5; // latitudes between the poles
                seeds.push_back(c + Vector3(std::sin(theta) * std::cos(phi), std::cos(theta), std::sin(theta) * std::sin(phi)) * r);
            }
    } else {
        for (int i = 1; i <= 4; ++i)
            for (int k = 1; k <= 4; ++k) seeds.push_back(dom.lo + dom.extent() * Vector3(i / 5.0f, 0.5f, k / 5.0f));
    }
}

// One line through a seed: half a cell per step along B (midpoint rule), backwards first (that
// half is reversed into the line), then forwards; |B| along it is kept for the colouring.
void Simulation::traceFieldLine(const Vector3& seed, std::vector<Vector3>& line, std::vector<float>& strength) const {
    const MagneticField& m = grid.magnetic;
    const float dx = grid.dx();
    const AABB dom = grid.domain();
    auto blocked = [&](const Vector3& p) {
        if (!dom.contains(p)) return true;
        const Vector3 g = (p - dom.lo) / dx;
        const int i = int(g.x), j = int(g.y), k = int(g.z);
        return i >= 0 && j >= 0 && k >= 0 && i < grid.nx() && j < grid.ny() && k < grid.nz() && grid.solid(i, j, k);
    };
    const float h = 0.5f * dx;
    const int maxSteps = 3 * (grid.nx() + grid.ny() + grid.nz());
    for (int dir = -1; dir <= 1; dir += 2) { // backwards (reversed into the line), then forwards
        std::vector<Vector3> half;
        std::vector<float> halfB;
        Vector3 p = seed;
        for (int it = 0; it < maxSteps && !blocked(p); ++it) {
            const Vector3 b1 = m.fieldAt(p);
            const float B = length(b1);
            if (B < 1e-9f) break;
            half.push_back(p);
            halfB.push_back(B);
            const Vector3 mid = p + b1 * (0.5f * h * float(dir) / B);
            const Vector3 b2 = m.fieldAt(mid);
            if (length2(b2) < 1e-18f) break;
            p += normalize(b2) * (h * float(dir));
        }
        if (dir < 0) {
            line.assign(half.rbegin(), half.rend());
            strength.assign(halfB.rbegin(), halfB.rend());
        } else if (half.size() > 1) {
            line.insert(line.end(), half.begin() + 1, half.end());
            strength.insert(strength.end(), halfB.begin() + 1, halfB.end());
        }
    }
}

void Simulation::computeStreamlines(RenderSnapshot& s) const {
    const float dx = grid.dx();
    const AABB dom = grid.domain();
    std::vector<Vector3> seeds;
    int n = std::max(2, vis.streamlineSeeds);
    if (grid.params.bc[0] == BoundaryType::Inflow) {
        AABB region;
        if (grid.hasObstacle() && obstacleMesh_ && !obstacleMesh_->empty()) {
            AABB ob = obstacleMesh_->bounds();
            Vector3 c = ob.center(), e = ob.extent() * 0.8f + Vector3(2 * dx);
            region = AABB(c - e, c + e);
        } else {
            region = AABB(dom.lo + dom.extent() * 0.2f, dom.hi - dom.extent() * 0.2f);
        }
        float y0 = std::max(region.lo.y, dom.lo.y + dx), y1 = std::min(region.hi.y, dom.hi.y - dx);
        float z0 = std::max(region.lo.z, dom.lo.z + dx), z1 = std::min(region.hi.z, dom.hi.z - dx);
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                seeds.push_back({dom.lo.x + 1.5f * dx, y0 + (y1 - y0) * (j + 0.5f) / n, z0 + (z1 - z0) * (k + 0.5f) / n});
    } else {
        Vector3 c = grid.source.center;
        float r = std::max(grid.source.radius * 1.5f, 2 * dx);
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                seeds.push_back(c + Vector3(-r + 2 * r * (j + 0.5f) / n, 0.0f, -r + 2 * r * (k + 0.5f) / n));
    }
    s.streamlines.assign(seeds.size(), {});
    s.streamlineSpeed.assign(seeds.size(), {});
    const float h = 0.5f * dx;
    const int maxSteps = 4 * (grid.nx() + grid.ny() + grid.nz());
    float vmaxAll = parallelMax<float>(int(seeds.size()), 1e-6f, [&](int s0, int s1) {
        float vm = 0;
        for (int si = s0; si < s1; ++si) {
            Vector3 p = seeds[si];
            auto& line = s.streamlines[si];
            auto& spd = s.streamlineSpeed[si];
            for (int it = 0; it < maxSteps; ++it) {
                if (!dom.contains(p)) break;
                Vector3 g = (p - dom.lo) / dx;
                int ci = int(g.x), cj = int(g.y), ck = int(g.z);
                if (ci >= 0 && cj >= 0 && ck >= 0 && ci < grid.nx() && cj < grid.ny() && ck < grid.nz() &&
                    grid.solid(ci, cj, ck))
                    break;
                Vector3 v1 = grid.velocityAt(p);
                float sp = length(v1);
                line.push_back(p);
                spd.push_back(sp);
                vm = std::max(vm, sp);
                if (sp < 1e-4f) break;
                Vector3 mid = p + v1 * (0.5f * h / sp);
                Vector3 v2 = grid.velocityAt(mid);
                float sp2 = length(v2);
                if (sp2 < 1e-4f) break;
                p += v2 * (h / sp2);
            }
        }
        return vm;
    }, 1);
    s.streamlineMax = vmaxAll;
}

// Cloth sheets and soft-body surfaces, in any mode (and the liquid outside the liquid mode, which
// draws it itself coloured by speed / density).
void Simulation::fillParticleSolids(RenderSnapshot& s) const {
    const auto& x = particles.positions();
    if (mode_ != SimMode::Fluid && particles.fluidCount() > 0) {
        s.particleRadius = particles.params.particleRadius;
        for (size_t i = 0; i < x.size(); ++i)
            if (particles.phases()[i] == uint8_t(ParticlePhase::Fluid)) {
                if (vis.liquidSurface) {
                    s.liquid.push_back(x[i]);
                    continue;
                }
                s.particles.push_back(x[i]);
                s.particleScalar.push_back(RenderSnapshot::kOwnColor);
                s.particleColor.push_back({0.2f, 0.5f, 0.95f});
            }
    }
    for (const Cloth& c : particles.cloths()) {
        RenderSnapshot::ClothMesh cm;
        cm.width = c.width;
        cm.height = c.height;
        cm.color = c.color;
        cm.positions.assign(x.begin() + c.firstParticle, x.begin() + c.firstParticle + c.width * c.height);
        cm.cellIntact = c.cellIntact;
        if (c.material.flammable)
            for (float u : c.unburnt) cm.burnt.push_back(1.0f - u);
        s.cloths.push_back(std::move(cm));
    }
    for (size_t b = 0; b < particles.softBodies().size(); ++b) {
        RenderSnapshot::SoftMesh sm;
        particles.softBodySurface(b, sm.positions);
        sm.triangles = particles.softBodies()[b].surface.triangles;
        sm.color = particles.softBodies()[b].color;
        s.softMeshes.push_back(std::move(sm));
    }
}

// Everything the viewer draws and the panels show, from the solvers as they are now. The parts
// in order: the settings (so the panels can sync), the bodies, then what the mode shows - the
// liquid, the gas, or the rigid arena - the joints and the mouse, the scene's own readings, and
// the probe's channels.
void Simulation::fillSnapshot(RenderSnapshot& s) const {
    fillSettings(s);
    fillBodies(s);
    switch (mode_) {
    case SimMode::Fluid: fillLiquidView(s); break;
    case SimMode::WindTunnel: fillGasView(s); break;
    case SimMode::Rigid: fillRigidView(s); break;
    }
    fillJointsAndGrab(s);
    s.lights.clear();                       // a scene with lights of its own fills them in describe()
    if (scene_) scene_->describe(*this, s); // the scene's own readings, after the generic ones
    drawResearchLayers();
    fillContacts(s);
    s.probe = Probe::snapshot();            // every channel the engine reported this frame
    s.info.insert(s.info.begin(), {"Время", format("%.3f с", time_)});
    s.info.push_back({"Шаг расчёта", format("%.1f мс", lastStepMs_)});
}

// The research layers that depend on the view (see core/Probe.h): the gas on the viewer's slice
// plane and the field lines B coloured by |B|. The solvers draw their own layers inside the step;
// these are drawn here, once per Probe frame however often the viewer asks for a snapshot.
void Simulation::drawResearchLayers() const {
    if (Probe::layers() == 0 || researchDrawnFrame_ == Probe::frameIndex()) return;
    researchDrawnFrame_ = Probe::frameIndex();
    if (mode_ == SimMode::WindTunnel) grid.drawDebug(vis.sliceAxis, sliceLayer());
    if (!Probe::layerOn(DrawLayer::FieldLinesB) || !grid.magnetic.enabled) return;
    RenderSnapshot lines;
    computeFieldLines(lines);
    for (size_t n = 0; n < lines.fieldLines.size(); ++n) {
        const std::vector<Vector3>& p = lines.fieldLines[n];
        const std::vector<float>& B = lines.fieldLineStrength[n];
        for (size_t i = 1; i < p.size(); ++i)
            Probe::line(DrawLayer::FieldLinesB, p[i - 1], p[i], heatColor(B[i] / lines.fieldLineMax));
    }
    Probe::label(DrawLayer::FieldLinesB, grid.domain().hi, format("|B| макс. %.3g Тл", double(lines.fieldLineMax)));
}

// The contact points of the last rigid step, while the ContactPoints layer is on.
void Simulation::fillContacts(RenderSnapshot& s) const {
    s.contacts.clear();
    if (!Probe::layerOn(DrawLayer::ContactPoints)) return;
    for (const RigidWorld::DebugContact& c : rigid.debugContacts())
        s.contacts.push_back({c.position, c.normal, c.depth, c.impulse, c.frictionImpulse, c.a, c.b});
}

// The current parameters and the frame's counters, and every list of the snapshot cleared.
void Simulation::fillSettings(RenderSnapshot& s) const {
    s.frame = frame_;
    s.mode = mode_;
    s.time = time_;
    s.stepMs = lastStepMs_;
    s.sceneName = scene_ ? scene_->name : std::string();
    s.sceneParams = sceneParams();
    s.paramsVersion = paramsVersion_;
    s.particleParams = particles.params;
    s.gasParams = grid.params;
    s.rigid = rigid.params;
    s.obstacleSettings = obstacle;
    s.vis = vis;
    s.emitter = particles.emitter;
    s.heat = grid.source;
    s.gasPushesBodies = gasPushesBodies;
    s.obstacleVersion = obstacleVersion_;
    s.obstacle = obstacleMesh_;
    s.measurements = measurements_;
    s.info.clear();
    s.plots.clear();
    s.particles.clear();
    s.particleScalar.clear();
    s.bodies.clear();
    s.obstacleScalar.clear();
    s.particleColor.clear();
    s.liquid.clear();
    s.liquidRadius = particles.params.particleRadius;
    s.cloths.clear();
    s.softMeshes.clear();
    fillParticleSolids(s);
    s.hasSlice = s.hasVolume = false;
    s.streamlines.clear();
    s.streamlineSpeed.clear();
    s.fieldLines.clear();
    s.fieldLineStrength.clear();
    s.volumePlasma = false;
    s.arrowPos.clear();
    s.arrowVel.clear();
    s.gridNx = s.gridNy = s.gridNz = 0;
}

// Every rigid body with its pose, what it looks like and its colour. A body drawn as a mesh - a
// hull, a compound, a capsule, or any body with its own visual mesh (a horse on a capsule collider)
// - is reported as a mesh (ConvexHull, or Compound when the collider is one, whose parts the viewer
// may show); collisionShape is always the collider, for picking.
void Simulation::fillBodies(RenderSnapshot& s) const {
    for (const RigidBody& b : rigid.bodies()) {
        // A destroyed body keeps its slot so that the viewer's body indices stay the solver's
        // (picking reports them back): an entry of zero size that draws nothing and cannot be hit.
        if (!b.alive) {
            s.bodies.push_back({ShapeType::Sphere, b.pos, b.rot, Vector3(0.0f), 0.0f, b.color, nullptr, true, nullptr, false});
            continue;
        }
        ShapeType drawn = b.type();
        std::shared_ptr<const TriMesh> mesh;
        if (b.type() == ShapeType::ConvexHull) mesh = static_cast<const ConvexHullShape*>(b.shape.get())->mesh();
        else if (b.type() == ShapeType::Compound) mesh = static_cast<const CompoundShape*>(b.shape.get())->visualMesh();
        else if (b.type() == ShapeType::Capsule) mesh = static_cast<const CapsuleShape*>(b.shape.get())->mesh();
        if (b.visualMesh) mesh = b.visualMesh;
        if (mesh && drawn != ShapeType::Compound) drawn = ShapeType::ConvexHull;
        s.bodies.push_back({drawn, b.pos, b.rot, b.halfExtents(), b.radius(), b.color, mesh, b.sleeping, b.shape, b.invMass > 0});
    }
}

// The colour scale of the main scalar: the measured range, or the fixed one the user set.
void Simulation::setColorRange(RenderSnapshot& s, float lo, float hi) const {
    if (vis.autoRange) { s.colorMin = lo; s.colorMax = hi > lo ? hi : lo + 1e-3f; }
    else { s.colorMin = vis.rangeMin; s.colorMax = vis.rangeMax > vis.rangeMin ? vis.rangeMax : vis.rangeMin + 1e-3f; }
}

// The liquid: every fluid particle (or the surface's particles) coloured by speed or density,
// and the readings of the particle solver.
void Simulation::fillLiquidView(RenderSnapshot& s) const {
    s.domain = particles.domain();
    s.particleRadius = particles.params.particleRadius;
    const auto& x = particles.positions();
    const auto& v = particles.velocities();
    const auto& rho = particles.densities();
    const auto& phase = particles.phases();
    s.particles.clear();
    s.particles.reserve(x.size());
    s.particleScalar.clear();
    s.particleScalar.reserve(x.size());
    float lo = kInf, hi = -kInf;
    for (size_t i = 0; i < x.size(); ++i) {
        if (phase[i] != uint8_t(ParticlePhase::Fluid)) continue; // cloth / soft bodies: drawn as surfaces
        if (vis.liquidSurface) {
            s.liquid.push_back(x[i]);
            continue;
        }
        s.particles.push_back(x[i]);
        float val = 0;
        if (vis.particleColoring == ParticleColoring::Speed) val = length(v[i]);
        else if (vis.particleColoring == ParticleColoring::Density)
            val = i < rho.size() ? rho[i] / particles.params.restDensity : 1.0f;
        s.particleScalar.push_back(val);
        lo = std::min(lo, val);
        hi = std::max(hi, val);
    }
    if (lo > hi) lo = 0, hi = 1;
    if (vis.particleColoring == ParticleColoring::Speed) { lo = 0; s.colorLabel = "|V|, м/с"; }
    else if (vis.particleColoring == ParticleColoring::Density) s.colorLabel = "ρ/ρ0";
    else s.colorLabel = "";
    if (s.particles.empty()) { lo = 0; hi = 1; }
    setColorRange(s, lo, hi);
    s.info.push_back({"Частиц", format("%zu", particles.size())});
    s.info.push_back({"Радиус частицы", format("%.1f мм", particles.params.particleRadius * 1000)});
    s.info.push_back({"Масса частицы", format("%.3g кг", particles.particleMass())});
    s.info.push_back({"Ошибка плотности", format("%.2f %%", particles.averageDensityError() * 100)});
    s.info.push_back({"Макс. скорость", format("%.2f м/с", particles.maxSpeed())});
    s.info.push_back({"Тел", format("%zu", rigid.bodies().size())});
    if (!particles.softBodies().empty() || !particles.cloths().empty()) {
        s.info.push_back({"Мягких тел / тканей", format("%zu / %zu", particles.softBodies().size(), particles.cloths().size())});
        s.info.push_back({"Частиц жидкости / твёрдых", format("%zu / %zu", particles.fluidCount(), particles.size() - particles.fluidCount())});
        s.info.push_back({"Контактов частиц", format("%zu", particles.particleContactCount())});
    }
    s.plots.push_back({"Ошибка плотности, %", particles.averageDensityError() * 100});
    s.plots.push_back({"Макс. скорость, м/с", particles.maxSpeed()});
    Probe::set("particles/density error %", particles.averageDensityError() * 100);
}

// The gas: the slice of the chosen field with its colour scale, the smoke volume, streamlines
// and magnetic field lines, the grid, the velocity arrows, the pressure on the obstacle, and
// the readings.
void Simulation::fillGasView(RenderSnapshot& s) const {
    s.domain = grid.domain();
    if (vis.showSlice) extractSlice(s);
    float lo = kInf, hi = -kInf;
    for (size_t i = 0; i < s.slice.size(); ++i)
        if (!s.sliceSolid[i]) { lo = std::min(lo, s.slice[i]); hi = std::max(hi, s.slice[i]); }
    if (s.slice.empty()) { lo = 0; hi = 1; }
    if (vis.sliceField == GridField::Speed || vis.sliceField == GridField::Vorticity) lo = 0;
    setColorRange(s, lo, hi);
    s.colorLabel = fieldLabel(vis.sliceField);
    if (vis.showSmoke) fillSmokeVolume(s);
    if (vis.showStreamlines) computeStreamlines(s);
    if (grid.magnetic.enabled && vis.showFieldLines) computeFieldLines(s);
    if (grid.magnetic.enabled) fillMagneticInfo(s);
    s.gridNx = grid.nx();
    s.gridNy = grid.ny();
    s.gridNz = grid.nz();
    s.gridDx = grid.dx();
    s.gridOrigin = grid.origin();
    s.sliceLayer = sliceLayer();
    if (vis.vectorDisplay > 0) extractVectors(s);
    fillSurfacePressure(s);
    fillGasInfo(s);
    if (grid.hasObstacle()) fillAerodynamicInfo(s);
    fillGasCounts(s);
}

// The smoke (and, in a fire, the temperature) as bytes for the volume renderer.
void Simulation::fillSmokeVolume(RenderSnapshot& s) const {
    const Field3& sm = grid.smoke();
    s.volX = sm.nx; s.volY = sm.ny; s.volZ = sm.nz;
    auto toByte = [](float v) { return uint8_t(clampv(v, 0.0f, 1.0f) * 255.0f + 0.5f); };
    if (grid.combustion.enabled) { // fire: smoke and temperature interleaved
        const Field3& T = grid.temperature();
        s.volumeChannels = 2;
        s.ambientTemperature = grid.combustion.ambientTemperature;
        s.volume.resize(2 * sm.d.size());
        for (size_t i = 0; i < sm.d.size(); ++i) {
            s.volume[2 * i] = toByte(sm.d[i]);
            s.volume[2 * i + 1] = toByte(T.d[i] / s.volumeTemperatureScale);
        }
    } else {
        s.volumeChannels = 1;
        s.volume.resize(sm.d.size());
        for (size_t i = 0; i < sm.d.size(); ++i) s.volume[i] = toByte(sm.d[i]);
    }
    s.volumePlasma = grid.magnetic.enabled; // the tracer is glowing plasma
    s.hasVolume = true;
}

// The plasma's numbers: the field, the Alfven speed, the magnetic Reynolds number, the energy,
// and how well div B = 0 holds.
void Simulation::fillMagneticInfo(RenderSnapshot& s) const {
    const MagneticField& m = grid.magnetic;
    const float B = m.maxField(), rho = grid.params.fluidDensity;
    const float L = grid.hasObstacle() ? std::max(obstacle.size, grid.dx()) : grid.domain().extent().x;
    s.info.push_back({"Магнитное поле, макс.", format("%.1f мТл", B * 1000)});
    s.info.push_back({"Скорость Альфвена v_A", format("%.2f м/с (предел Бориса %.1f)", B / std::sqrt(MagneticField::kMu0 * rho),
                                                 m.speedLimit)});
    s.info.push_back({"Маг. число Рейнольдса Rm", format("%.3g", grid.params.inflowSpeed * L / m.resistivity())});
    s.info.push_back({"Энергия поля", format("%.3g Дж", m.energy())});
    s.info.push_back({"div B (отн.)", format("%.1e", m.maxDivergence())});
}

// Cp of the obstacle's triangles, area-weighted onto the vertices for smooth colouring.
void Simulation::fillSurfacePressure(RenderSnapshot& s) const {
    if (!(vis.surfacePressure && obstacleMesh_ && !obstacleMesh_->empty() &&
          surfaceLoads_.triangles.size() == obstacleMesh_->triangles.size()))
        return;
    const TriMesh& m = *obstacleMesh_;
    std::vector<float> sum(m.positions.size(), 0.0f), wsum(m.positions.size(), 0.0f);
    for (size_t t = 0; t < m.triangles.size(); ++t) {
        const TriangleLoad& L = surfaceLoads_.triangles[t];
        for (uint32_t v : m.triangles[t]) {
            sum[v] += L.cp * L.area;
            wsum[v] += L.area;
        }
    }
    s.obstacleScalar.resize(m.positions.size());
    for (size_t v = 0; v < sum.size(); ++v) s.obstacleScalar[v] = wsum[v] > 0 ? sum[v] / wsum[v] : 0.0f;
}

// The grid, the time step, the fire's power and temperature, the Reynolds number.
void Simulation::fillGasInfo(RenderSnapshot& s) const {
    float L = obstacle.size;
    float Re = grid.params.inflowSpeed * L / std::max(grid.params.kinematicViscosity, 1e-9f);
    s.info.push_back({"Сетка", format("%d × %d × %d  (%.1f тыс. ячеек)", grid.nx(), grid.ny(), grid.nz(),
                                   grid.nx() * grid.ny() * grid.nz() / 1000.0f)});
    s.info.push_back({"Шаг сетки dx", format("%.1f мм", grid.dx() * 1000)});
    s.info.push_back({"Шаг по времени", format("%.2f мс", lastGridDt_ * 1000)});
    if (grid.combustion.enabled) {
        float Tmax = 0;
        for (float t : grid.temperature().d) Tmax = std::max(Tmax, t);
        s.info.push_back({"Мощность пламени", format("%.1f кВт", grid.heatReleaseRate() / 1000)});
        s.info.push_back({"Макс. температура", format("%.0f K (%.0f °C)", Tmax + grid.combustion.ambientTemperature, Tmax + grid.combustion.ambientTemperature - 273.15f)});
        int burnt = 0;
        for (const Cloth& c : particles.cloths()) burnt += c.burntThreads;
        if (!particles.cloths().empty()) s.info.push_back({"Прогоревших нитей", format("%d", burnt)});
    }
    if (grid.hasObstacle()) s.info.push_back({"Число Рейнольдса", format("%.3g", Re)});
}

// The aerodynamic coefficients of the obstacle: drag, lift, their ratio, the force, the skin
// friction's share, and the same from the surface loads on the polygons.
void Simulation::fillAerodynamicInfo(RenderSnapshot& s) const {
    Vector3 F = grid.bodyForce();
    s.info.push_back({"Cd (сред.)", format("%.3f", grid.dragCoefficientAvg())});
    s.info.push_back({"Cl (сред.)", format("%.3f", grid.liftCoefficientAvg())});
    s.info.push_back({"Cd / Cl мгн.", format("%.3f / %.3f", grid.dragCoefficient(), grid.liftCoefficient())});
    if (std::fabs(grid.dragCoefficientAvg()) > 1e-4f)
        s.info.push_back({"L/D", format("%.2f", grid.liftCoefficientAvg() / grid.dragCoefficientAvg())});
    s.info.push_back({"Сила F", format("(%.2f, %.2f, %.2f) Н", F.x, F.y, F.z)});
    if (grid.params.wallFriction) {
        float cdf = grid.frictionForce().x / (grid.dynamicPressure() * std::max(grid.referenceArea(), 1e-9f));
        s.info.push_back({"из них трение (Cd тр.)", format("%.3f", cdf)});
    }
    const SurfaceLoads& SL = surfaceLoads_;
    if (!SL.triangles.empty()) {
        s.info.push_back({"По полигонам: Cd", format("%.3f (давл. %.3f + трение %.3f)", SL.cd, SL.cdPressure, SL.cdFriction)});
        s.info.push_back({"По полигонам: Cl / Cm", format("%.3f / %.3f", SL.cl, SL.cm)});
        s.info.push_back({"Треугольников, смоч. площадь", format("%zu, %.4f м²", SL.triangles.size(), SL.wettedArea)});
    }
    s.info.push_back({grid.params.usePlanformArea ? "Площадь (в плане)" : "Площадь (миделя)",
                      format("%.4f м²", grid.referenceArea())});
    s.plots.push_back({"Cd", grid.dragCoefficient()});
    s.plots.push_back({"Cl", grid.liftCoefficient()});
    Probe::set("gas/Cd", grid.dragCoefficient());
    Probe::set("gas/Cl", grid.liftCoefficient());
}

// The solver's counters: pressure iterations, speeds, smoke, the water, the solids in the gas.
void Simulation::fillGasCounts(RenderSnapshot& s) const {
    s.info.push_back({"Итераций давления", format("%d (невязка %.1e)", grid.lastPressureIterations(), grid.lastResidual())});
    s.info.push_back({"Макс. скорость", format("%.2f м/с", grid.maxVelocity())});
    s.info.push_back({"Макс. |div u| после проекции", format("%.2e 1/с", grid.maxDivergence())});
    s.info.push_back({"Дым в объёме", format("%.4f м³", grid.totalSmoke())});
    if (particles.fluidCount() > 0)
        s.info.push_back({"Вода: частиц / ячеек в воздухе", format("%zu / %d", particles.fluidCount(), grid.liquidCellCount())});
    if (particles.hasSolids()) {
        int torn = 0;
        for (const Cloth& c : particles.cloths()) torn += c.tornThreads;
        s.info.push_back({"Мягких тел / тканей", format("%zu / %zu (порвано нитей %d)", particles.softBodies().size(),
                                                         particles.cloths().size(), torn)});
    }
    if (!rigid.bodies().empty()) {
        s.info.push_back({"Тел в газе", format("%zu (%d ячеек, спят %zu)", rigid.bodies().size(), grid.movingSolidCells(),
                                            rigid.sleepingCount())});
        s.info.push_back({"Сила газа на тело (макс.)", format("%.3f Н", gasForceMax_)});
        if (rigid.anyHeld()) s.info.push_back({"Тела отпустятся через", format("%.1f с", std::max(0.0, double(releaseTime) - time_))});
    }
    if (!s.arrowPos.empty()) s.info.push_back({"Векторов скорости", format("%zu", s.arrowPos.size())});
    if (!grid.hasObstacle()) {
        s.plots.push_back({"Макс. скорость, м/с", grid.maxVelocity()});
        s.plots.push_back({"Дым, дм³", grid.totalSmoke() * 1000.0f});
        Probe::set("gas/smoke dm3", grid.totalSmoke() * 1000.0f);
    }
}

// The rigid arena: its box and the solver's counters.
void Simulation::fillRigidView(RenderSnapshot& s) const {
    s.domain = rigid.domain();
    s.colorMin = 0;
    s.colorMax = 1;
    s.info.push_back({"Тел", format("%zu", rigid.bodies().size())});
    s.info.push_back({"Контактов", format("%zu", rigid.contactCount())});
    s.info.push_back({"Спящих тел", format("%zu", rigid.sleepingCount())});
    s.info.push_back({"CCD: остановлено тел за шаг", format("%zu", rigid.ccdHits())});
    s.info.push_back({"Сочленений", format("%zu", rigid.joints().size())});
    s.info.push_back({"Кин. энергия", format("%.2f Дж", rigid.kineticEnergy())});
    s.plots.push_back({"Кин. энергия, Дж", rigid.kineticEnergy()});
    Probe::set("rigid/kinetic energy J", rigid.kineticEnergy());
}

// The joints' anchors and axes, and the mouse joint (of a body or of grabbed particles).
void Simulation::fillJointsAndGrab(RenderSnapshot& s) const {
    s.joints.clear();
    for (const auto& j : rigid.joints())
        s.joints.push_back({j->type(), j->worldAnchorA(rigid.bodies()), j->worldAnchorB(rigid.bodies()), j->worldAxis(rigid.bodies())});
    const auto& gj = rigid.grabJoint();
    s.grabActive = gj.active && gj.body >= 0 && gj.body < int(rigid.bodies().size());
    if (s.grabActive) {
        const RigidBody& gb = rigid.bodies()[gj.body];
        s.grabAnchor = gb.pos + gb.rotation() * gj.localAnchor;
        s.grabTarget = gj.target;
    } else if (particles.grabbing()) {
        s.grabActive = true;
        s.grabAnchor = particles.grabAnchor();
        s.grabTarget = particles.grabTarget();
    }
}

} // namespace rf
