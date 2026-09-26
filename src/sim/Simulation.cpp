#include "sim/Simulation.h"

#include "core/Parallel.h"
#include "rigid/ConvexDecomposition.h"

#include <chrono>
#include <cstdarg>
#include <cstdio>

namespace rf {

const char* presetName(Preset p) {
    switch (p) {
    case Preset::DamBreak: return "Жидкость: разрушение плотины";
    case Preset::FluidObstacle: return "Жидкость: поток на препятствие";
    case Preset::FloatingBodies: return "Жидкость: плавающие тела";
    case Preset::JetOnObject: return "Жидкость: струя на объект";
    case Preset::TunnelSphere: return "Аэротруба: сфера";
    case Preset::TunnelCylinder: return "Аэротруба: цилиндр";
    case Preset::TunnelWing: return "Аэротруба: крыло NACA";
    case Preset::TunnelStreamlined: return "Аэротруба: обтекаемое тело";
    case Preset::TunnelCube: return "Аэротруба: куб";
    case Preset::SmokePlume: return "Газ: тепловой шлейф дыма";
    case Preset::SmokeSphere: return "Дым: сферический источник (закрытый объём)";
    case Preset::RigidFalling: return "Твёрдые тела: падение на меш";
    case Preset::RigidGranular: return "Твёрдые частицы: сыпучая среда";
    case Preset::RigidPyramid: return "Твёрдые тела: пирамида и снаряд";
    case Preset::RigidConvex: return "Твёрдые тела: многогранники (SAT, GJK-EPA)";
    case Preset::RigidTower: return "Твёрдые тела: башня из 100 кубиков";
    case Preset::RigidJoints: return "Сочленения: 5 типов";
    case Preset::RigidCcd: return "CCD: пули и тонкая стена";
    case Preset::RigidTeapots: return "Невыпуклые: 100 чайников (выпуклая декомпозиция)";
    case Preset::SmokeBodies: return "Дым + твёрдые тела (двусторонняя связь)";
    case Preset::SoftCloth: return "Мягкие тела и ткань (единый решатель частиц)";
    case Preset::GasSoftCloth: return "Газ + мягкие тела + ткань + твёрдые тела";
    case Preset::Hydro: return "Гидродинамика: вода + воздух + тела";
    case Preset::Fire: return "Огонь: горелка, горящая штора, тела";
    case Preset::Water: return "Вода: волна в бассейне, плавающие тела (шейдер)";
    case Preset::Magnetosphere: return "Плазма: магнит отклоняет поток (магнитосфера)";
    default: return "?";
    }
}

SimMode presetMode(Preset p) {
    switch (p) {
    case Preset::DamBreak: case Preset::FluidObstacle: case Preset::FloatingBodies: case Preset::JetOnObject:
    case Preset::SoftCloth: case Preset::Water:
        return SimMode::Fluid;
    case Preset::TunnelSphere: case Preset::TunnelCylinder: case Preset::TunnelWing: case Preset::TunnelStreamlined:
    case Preset::TunnelCube: case Preset::SmokePlume: case Preset::SmokeSphere: case Preset::SmokeBodies:
    case Preset::GasSoftCloth: case Preset::Hydro: case Preset::Fire: case Preset::Magnetosphere:
        return SimMode::WindTunnel;
    default:
        return SimMode::Rigid;
    }
}

static std::string fmt(const char* f, ...) {
    char buf[256];
    va_list ap;
    va_start(ap, f);
    std::vsnprintf(buf, sizeof(buf), f, ap);
    va_end(ap);
    return buf;
}

static float deg(float d) { return d * kPi / 180.0f; }

Simulation::Simulation() {
    obstacleMesh_ = std::make_shared<TriMesh>();
    fluidDomain_ = AABB({-1.0f, 0.0f, -0.4f}, {1.0f, 1.2f, 0.4f});
    rigidDomain_ = AABB({-2.0f, 0.0f, -2.0f}, {2.0f, 5.0f, 2.0f});
    loadPreset(Preset::DamBreak);
}

// ---------------------------------------------------------------------------
// Presets
// ---------------------------------------------------------------------------
void Simulation::loadPreset(Preset p) {
    preset_ = p;
    mode_ = presetMode(p);
    obstacle = ObstacleSettings();
    vis = VisSettings();
    particles.params = ParticleParams();
    particles.emitter = ParticleEmitter();
    grid.params = NSParams();
    grid.source = HeatSource();
    grid.combustion = Combustion();
    grid.magnetic = MagneticField();
    rigid.params = RigidParams();

    switch (p) {
    case Preset::DamBreak:
        break;
    case Preset::FluidObstacle:
        obstacle.shape = ObstacleShape::Sphere;
        obstacle.size = 0.35f;
        obstacle.position = {0.3f, 0.12f, 0.0f};
        break;
    case Preset::FloatingBodies:
        particles.params.particleRadius = 0.016f;
        break;
    case Preset::JetOnObject:
        obstacle.shape = ObstacleShape::Wing;
        obstacle.size = 0.5f;
        obstacle.span = 0.6f;
        obstacle.angleOfAttackDeg = 25.0f;
        obstacle.position = {0.1f, 0.35f, 0.0f};
        particles.emitter.enabled = true;
        particles.emitter.position = {-0.9f, 0.75f, 0.0f};
        particles.emitter.direction = normalize(Vector3(1.0f, -0.25f, 0.0f));
        particles.emitter.radius = 0.07f;
        particles.emitter.speed = 3.0f;
        particles.params.maxParticles = 80000;
        break;
    case Preset::TunnelSphere:
        obstacle.shape = ObstacleShape::Sphere;
        obstacle.size = 0.5f;
        break;
    case Preset::TunnelCylinder:
        obstacle.shape = ObstacleShape::Cylinder;
        obstacle.size = 0.4f;
        obstacle.span = 1.4f;
        vis.sliceField = GridField::Vorticity;
        grid.params.vorticityConfinement = 0.5f;
        break;
    case Preset::TunnelWing:
        obstacle.shape = ObstacleShape::Wing;
        obstacle.size = 0.9f;
        obstacle.span = 1.2f;
        obstacle.angleOfAttackDeg = 6.0f;
        grid.params.usePlanformArea = true;
        vis.sliceField = GridField::PressureCoeff;
        break;
    case Preset::TunnelStreamlined:
        obstacle.shape = ObstacleShape::Streamlined;
        obstacle.size = 1.0f;
        obstacle.thickness = 0.3f;
        break;
    case Preset::TunnelCube:
        obstacle.shape = ObstacleShape::Cube;
        obstacle.size = 0.4f;
        break;
    case Preset::SmokePlume:
        obstacle.shape = ObstacleShape::Sphere;
        obstacle.size = 0.45f;
        obstacle.position = {0.0f, 0.1f, 0.0f};
        grid.params.domainSize = {2.0f, 3.0f, 2.0f};
        grid.params.resolutionX = 48;
        grid.params.inflowSpeed = 1.0f;
        grid.params.bc[0] = BoundaryType::Wall;
        grid.params.bc[1] = BoundaryType::Wall;
        grid.params.bc[3] = BoundaryType::Outflow;
        grid.params.smokeRake = false;
        grid.params.heatBuoyancy = 4.0f;
        grid.params.smokeBuoyancy = 0.3f;
        grid.params.vorticityConfinement = 1.5f;
        grid.params.smokeDissipation = 0.05f;
        grid.source.enabled = true;
        grid.source.center = {0.0f, -1.25f, 0.0f};
        grid.source.radius = 0.2f;
        vis.sliceField = GridField::Temperature;
        vis.showStreamlines = false;
        break;
    case Preset::SmokeSphere:
    case Preset::SmokeBodies:
    case Preset::GasSoftCloth:
        // Closed box of voxels, hot smoky sphere near the bottom, gas rises by buoyancy.
        grid.params.domainSize = {1.6f, 2.4f, 1.6f};
        grid.params.resolutionX = 32;               // dx = 5 cm -> 32 x 48 x 32 cells
        grid.params.inflowSpeed = 1.0f;             // only a reference scale for Cp here
        for (auto& b : grid.params.bc) b = BoundaryType::Wall;
        grid.params.smokeRake = false;
        grid.params.heatBuoyancy = 3.0f;
        grid.params.smokeBuoyancy = 0.2f;
        grid.params.temperatureDissipation = 0.4f;
        grid.params.smokeDissipation = 0.08f;       // tracer fades slowly so the closed box does not fog up
        grid.params.vorticityConfinement = 1.0f;
        grid.source.enabled = true;
        grid.source.center = {0.0f, -0.8f, 0.0f};
        grid.source.radius = 0.15f;
        vis.sliceField = GridField::Speed;
        vis.showSlice = false;
        vis.showStreamlines = false;
        // Scenes with bodies: smoke only (no voxel grid / vectors), 3/4 view, dense smoke that
        // leaves through the open ceiling.
        vis.gridDisplay = p == Preset::SmokeSphere ? 1 : 0;
        vis.vectorDisplay = p == Preset::SmokeSphere ? 1 : 0;
        // Cloth in the gas scene only flutters in the flow (no bodies thrown at it): the coarser
        // cloth spacing 2r - four times fewer particles - is enough.
        if (p == Preset::GasSoftCloth) particles.params.clothSpacing = 2.0f;
        if (p != Preset::SmokeSphere) { // dense smoke, open ceiling: the plume leaves, the box does not fill up
            grid.source.smoke = 2.0f;
            grid.source.radius = 0.2f;
            grid.params.smokeDissipation = 0.3f; // smoke lives ~3 s: the plume stays dense, the box clears
            grid.params.bc[3] = BoundaryType::Outflow;
        }
        break;
    case Preset::Magnetosphere:
        // A plasma wind blows at a magnetised sphere (Birkeland's terrella; the solar wind at the
        // Earth): the dipole's magnetic pressure B^2/2mu0 stops the flow where it matches the ram
        // pressure rho u^2 (the Chapman-Ferraro magnetopause, ~0.3 m upstream here); the plasma is
        // deflected around the magnetosphere and the field lines are swept back into a tail.
        // Conductivity 1e8 S/m: magnetic Reynolds number u L / eta ~ 100, the field is frozen in.
        obstacle.shape = ObstacleShape::Sphere;
        obstacle.size = 0.3f;
        obstacle.position = Vector3(0.0f);
        grid.params.domainSize = {2.4f, 1.4f, 1.4f}; // wide enough that the walls do not squeeze the flow
        grid.params.resolutionX = 60;       // dx = 4 cm
        grid.params.inflowSpeed = 1.5f;
        grid.params.fluidDensity = 1.0f;
        grid.params.pressureTolerance = 1e-3f;
        grid.params.smokeDissipation = 0.7f; // the streaks fade in ~1.5 s: fresh wind only, no fog
        grid.params.smokeRake = true;       // plasma streaks from the inflow show the deflection
        grid.magnetic.enabled = true;
        grid.magnetic.conductivity = 1e8f;
        grid.magnetic.speedLimit = 4.0f;    // Boris correction near the poles (v_A up to 80 m/s there)
        vis.sliceField = GridField::MagneticFlux;
        vis.planetSurface = true;           // the terrella as a little Earth
        vis.surfacePressure = false;
        vis.showSlice = false;
        vis.showStreamlines = false;
        vis.gridDisplay = 0;
        // Velocity arrows in the equatorial plane (the dipole points along y): the wind slowing at
        // the magnetopause, turning around the magnetosphere and closing behind it.
        vis.vectorDisplay = 1;
        vis.sliceAxis = 1;
        vis.slicePosition = 0.5f;
        vis.vectorStride = 2;
        vis.vectorScale = 1.3f;
        break;
    case Preset::Water:
        particles.params.particleRadius = 0.012f; // ~40 000 particles: a smooth enough surface
        vis.liquidSurface = true;
        break;
    case Preset::Fire:
        // A gas burner on the floor of a room-corner box open at the top; a cotton curtain hangs
        // next to it. The flame heats the curtain's lower edge until it ignites; the burning fabric
        // gives off fuel gas that burns in the air, so the fire climbs the curtain, which chars
        // through and falls apart. Temperatures in kelvin above ambient (Combustion).
        grid.params.domainSize = {1.2f, 1.6f, 1.0f};
        grid.params.resolutionX = 40;                // dx = 3 cm -> 40 x 53 x 33 cells
        grid.params.inflowSpeed = 1.0f;              // reference scale only
        for (auto& b : grid.params.bc) b = BoundaryType::Wall;
        grid.params.bc[3] = BoundaryType::Outflow;   // open top: the hot gas and the smoke leave
        grid.params.smokeRake = false;
        grid.params.smokeDissipation = 0.15f;
        grid.params.vorticityConfinement = 2.0f;     // small eddies of the flame lost on the grid
        grid.params.pressureTolerance = 1e-3f;
        grid.combustion.enabled = true;
        grid.source.enabled = true;                  // the burner: fuel gas, lit by a pilot flame
        grid.source.center = {0.0f, -0.75f, 0.03f}; // its flame licks the curtain's lower edge
        grid.source.radius = 0.06f;
        grid.source.fuel = 1.0f;
        grid.source.temperature = 400.0f;            // above ignition: burns as it leaves the burner
        grid.source.smoke = 0.0f;
        grid.source.velocity = {0.0f, 0.5f, 0.0f};
        particles.params.clothSpacing = 2.0f;        // 3 cm - the grid spacing
        vis.sliceField = GridField::Temperature;
        vis.showSlice = false;
        vis.showStreamlines = false;
        vis.gridDisplay = 0;
        vis.vectorDisplay = 0;
        break;
    case Preset::Hydro:
        // Water and air in one box: wind from the left over a pool, a water column that collapses
        // into a wave, floating and sinking bodies, a flag in the wind.
        grid.params.domainSize = {2.0f, 1.2f, 0.8f};
        // The air only has to carry the smoke streaks and push spray and flag: a coarser grid,
        // a looser pressure tolerance and larger (still stable, semi-Lagrangian) steps.
        grid.params.resolutionX = 48;                  // dx ~4 cm
        grid.params.pressureTolerance = 1e-3f;
        grid.params.cfl = 4.0f;
        grid.params.inflowSpeed = 3.0f;
        for (auto& b : grid.params.bc) b = BoundaryType::Wall;
        grid.params.bc[0] = BoundaryType::Inflow;       // wind
        grid.params.bc[1] = BoundaryType::Outflow;
        grid.params.bc[3] = BoundaryType::Outflow;      // open top
        grid.params.smokeRake = true;                   // streaks show the air flow
        grid.params.smokeDissipation = 0.15f;
        grid.params.vorticityConfinement = 0.5f;
        particles.params.particleRadius = 0.018f;
        particles.params.clothSpacing = 2.0f;
        vis.showSlice = false;
        vis.showStreamlines = false;
        vis.gridDisplay = 0;
        vis.vectorDisplay = 0;
        break;
    case Preset::RigidFalling:
        obstacle.shape = ObstacleShape::Sphere;
        obstacle.size = 1.6f;
        obstacle.position = {0.0f, 0.0f, 0.0f};
        break;
    case Preset::RigidGranular:
        obstacle.shape = ObstacleShape::Cone;
        obstacle.size = 1.2f;
        obstacle.rollDeg = 0.0f;
        obstacle.angleOfAttackDeg = 90.0f; // apex up
        obstacle.position = {0.0f, 0.6f, 0.0f};
        rigid.params.iterations = 8;
        rigid.params.substeps = 4; // loose spheres: no tall stacks, fewer substeps suffice
        break;
    case Preset::RigidPyramid:
        break;
    case Preset::RigidConvex:
        obstacle.shape = ObstacleShape::Cube;
        obstacle.size = 1.4f;
        obstacle.position = {0.0f, 0.1f, 0.0f};
        obstacle.angleOfAttackDeg = 20.0f;
        obstacle.yawDeg = 30.0f;
        break;
    default:
        break;
    }
    buildObstacleGeometry();
    reset();
    ++paramsVersion_;
}

void Simulation::buildObstacleGeometry() {
    TriMesh m;
    const float s = std::max(obstacle.size, 1e-3f);
    switch (obstacle.shape) {
    case ObstacleShape::None: break;
    case ObstacleShape::Sphere: m = primitives::sphere(0.5f * s); break;
    case ObstacleShape::Cube: m = primitives::box(Vector3(0.5f * s)); break;
    case ObstacleShape::Cylinder: m = primitives::cylinder(0.5f * s, obstacle.span); break;
    case ObstacleShape::Wing:
        m = primitives::nacaWing(obstacle.nacaCode, s, obstacle.span);
        m.translate({-0.35f * s, 0, 0}); // rotate about ~ the aerodynamic centre
        break;
    case ObstacleShape::Streamlined: m = primitives::streamlinedBody(s, obstacle.thickness); break;
    case ObstacleShape::Ellipsoid: m = primitives::ellipsoid({0.5f * s, 0.25f * s, 0.25f * s}); break;
    case ObstacleShape::Cone: m = primitives::cone(0.5f * s, s); break;
    case ObstacleShape::Custom:
        if (obstacle.customMesh) {
            m = *obstacle.customMesh;
            m.fitTo(Vector3(0.0f), s);
        }
        break;
    }
    if (!m.empty()) {
        Quaternion q = Quaternion::fromEuler(deg(obstacle.yawDeg), -deg(obstacle.angleOfAttackDeg), deg(obstacle.rollDeg));
        m.transform(q.toMatrix3x3(), Vector3(1.0f), obstacle.position);
    }
    obstacleMesh_ = std::make_shared<TriMesh>(std::move(m));
    obstacleBVH_.build(*obstacleMesh_);
    ++obstacleVersion_;
}

void Simulation::rebuildObstacle() {
    buildObstacleGeometry();
    reset();
    ++paramsVersion_;
}

bool Simulation::loadCustomMesh(const std::string& path, std::string& error) {
    auto m = std::make_shared<TriMesh>();
    if (!loadMesh(path, *m, error)) return false;
    obstacle.customMesh = m;
    auto slash = path.find_last_of("/\\");
    obstacle.customName = slash == std::string::npos ? path : path.substr(slash + 1);
    obstacle.shape = ObstacleShape::Custom;
    rebuildObstacle();
    return true;
}

void Simulation::reset() {
    time_ = 0;
    frame_ = 0;
    switch (mode_) {
    case SimMode::Fluid: setupFluidScene(); break;
    case SimMode::WindTunnel: setupTunnelScene(); break;
    case SimMode::Rigid: setupRigidScene(); break;
    }
}

void Simulation::setupFluidScene() {
    const AABB& d = fluidDomain_;
    rigid.clear();
    rigid.setDomain(d);
    rigid.setStaticMesh(&obstacleBVH_);
    particles.setStaticMesh(&obstacleBVH_);
    particles.setRigidWorld(&rigid);
    particles.reset(d);
    switch (preset_) {
    case Preset::DamBreak:
    case Preset::FluidObstacle:
        particles.addBlock(AABB(d.lo, {d.lo.x + 0.6f, 0.7f, d.hi.z}));
        break;
    case Preset::SoftCloth: {
        // Left: a canvas hammock pinned at its corners with a foam cube and a jelly ball dropped
        // on it. Right: a pool with a curtain hanging above it and a floating soft cube.
        particles.addBlock(AABB({0.1f, d.lo.y, d.lo.z}, {d.hi.x, 0.22f, d.hi.z}));
        ClothMaterial canvas;
        canvas.areaDensity = 1.5f;
        canvas.tensileStiffness = 0.0f; // inextensible canvas ...
        canvas.bendCompliance = 1e-4f;
        canvas.strengthWarp = canvas.strengthWeft = 0.0f; // ... that does not tear
        particles.addCloth({-0.9f, 0.55f, -0.3f}, {0.6f, 0, 0}, {0, 0, 0.6f}, canvas, 1 | 2 | 4 | 8, {0.85f, 0.72f, 0.4f});
        // Cotton curtain on a rod, sewn from two panels: the vertical seam in the middle holds 30 %
        // of the fabric strength - pull a panel hard sideways and it opens along the seam.
        ClothMaterial cotton;
        cotton.areaDensity = 0.3f;
        cotton.seamColumns = {20}; // 41 columns: the middle
        particles.addCloth({0.4f, 1.1f, -0.3f}, {0, 0, 0.6f}, {0, -0.55f, 0}, cotton, 16 /* on a rod */, {0.75f, 0.3f, 0.35f});
        TriMesh cube = primitives::box(Vector3(0.09f));
        cube.translate({-0.65f, 0.95f, 0.0f});
        particles.addSoftBody(cube, 150.0f, 0.4f, {0.3f, 0.75f, 0.95f});
        TriMesh ball = primitives::sphere(0.08f, 16, 8);
        ball.translate({-0.45f, 1.1f, 0.1f});
        particles.addSoftBody(ball, 150.0f, 0.15f, {0.55f, 0.9f, 0.35f});
        TriMesh floater = primitives::box(Vector3(0.07f));
        floater.translate({0.75f, 0.6f, 0.0f});
        particles.addSoftBody(floater, 500.0f, 0.6f, {0.95f, 0.6f, 0.2f});
        break;
    }
    case Preset::FloatingBodies: {
        rigid.addBox({-0.5f, 0.75f, 0.0f}, {0.12f, 0.06f, 0.1f}, Quaternion::fromAxisAngle({0, 0, 1}, 0.4f), 400.0f,
                     {0.85f, 0.55f, 0.25f});
        rigid.addBox({0.0f, 0.9f, 0.05f}, {0.1f, 0.1f, 0.1f}, Quaternion::fromAxisAngle({1, 1, 0}, 0.7f), 600.0f,
                     {0.3f, 0.75f, 0.35f});
        rigid.addSphere({0.5f, 0.8f, -0.1f}, 0.09f, 300.0f, {0.9f, 0.3f, 0.3f});
        rigid.addBox({0.4f, 1.05f, 0.1f}, {0.06f, 0.06f, 0.06f}, Quaternion(), 2500.0f, {0.5f, 0.5f, 0.55f});
        particles.addBlock(AABB(d.lo, {d.hi.x, 0.32f, d.hi.z}));
        break;
    }
    case Preset::JetOnObject:
        particles.addBlock(AABB(d.lo, {d.hi.x, 0.08f, d.hi.z}));
        break;
    case Preset::Water: {
        // A pool with floating things; a column of water at the left end collapses into a wave
        // that runs through them. Densities: pine 500, plank 600, beach ball 80, teapot (hollow
        // in reality, solid here) 700, steel-like ball 3000 (sinks), soft foam 300.
        rigid.addBox({-0.25f, 0.3f, -0.15f}, Vector3(0.07f), Quaternion::fromAxisAngle({0, 1, 0}, 0.4f), 500.0f,
                     {0.85f, 0.6f, 0.3f});
        rigid.addBox({0.1f, 0.3f, 0.15f}, {0.16f, 0.025f, 0.06f}, Quaternion::fromAxisAngle({0, 1, 0}, -0.3f), 600.0f,
                     {0.7f, 0.45f, 0.25f});
        rigid.addSphere({0.45f, 0.35f, -0.1f}, 0.08f, 80.0f, {0.95f, 0.3f, 0.25f});
        rigid.addSphere({0.7f, 0.5f, 0.2f}, 0.05f, 3000.0f, {0.6f, 0.62f, 0.66f});
        rigid.addCompound(teapotShape(), {0.3f, 0.4f, 0.2f}, Quaternion::fromAxisAngle({0, 1, 0}, 2.5f), 700.0f,
                          {0.95f, 0.95f, 0.92f});
        TriMesh foam = primitives::box(Vector3(0.06f));
        foam.translate({0.75f, 0.35f, -0.2f});
        particles.addSoftBody(foam, 300.0f, 0.4f, {0.4f, 0.85f, 0.4f});
        particles.addBlock(AABB(d.lo, {d.hi.x, 0.2f, d.hi.z}));                      // the pool
        particles.addBlock(AABB({d.lo.x, 0.2f, d.lo.z}, {d.lo.x + 0.45f, 0.75f, d.hi.z})); // the column
        break;
    }
    default:
        particles.addBlock(AABB(d.lo, {d.lo.x + 0.6f, 0.7f, d.hi.z}));
        break;
    }
}

// A solid given as overlapping closed parts -> minimal set of convex hulls fitted to its smooth surface.
static std::shared_ptr<const CompoundShape> decompose(const std::vector<TriMesh>& parts) {
    DecompositionParams prm;
    prm.resolution = 24;
    prm.maxHullVertices = 64;
    return std::make_shared<CompoundShape>(convexDecomposition(parts, prm), primitives::merge(parts));
}

std::shared_ptr<const CompoundShape> teapotShape() {
    static std::shared_ptr<const CompoundShape> shape = decompose(primitives::teapotParts(0.3f));
    return shape;
}

std::shared_ptr<const CompoundShape> bunnyShape() {
    static std::shared_ptr<const CompoundShape> shape = decompose(primitives::bunnyParts(0.3f));
    return shape;
}

void Simulation::setupTunnelScene() {
    Vector3 size = grid.params.domainSize;
    Vector3 origin;
    if (preset_ == Preset::SmokePlume || preset_ == Preset::SmokeSphere || preset_ == Preset::SmokeBodies ||
        preset_ == Preset::GasSoftCloth || preset_ == Preset::Fire)
        origin = Vector3(-0.5f * size.x, -0.5f * size.y, -0.5f * size.z);
    else if (preset_ == Preset::Hydro)
        origin = Vector3(-0.5f * size.x, 0.0f, -0.5f * size.z); // floor at y = 0
    else origin = Vector3(-0.3f * size.x, -0.5f * size.y, -0.5f * size.z); // body at 30% of the length
    grid.reset(origin, &obstacleBVH_);

    // Rigid bodies live in the gas box (walls = domain, the tunnel obstacle is static geometry).
    rigid.clear();
    rigid.setDomain(grid.domain());
    rigid.setStaticMesh(&obstacleBVH_);
    // The particle system (soft bodies, cloth) shares the gas box.
    particles.setStaticMesh(&obstacleBVH_);
    particles.setRigidWorld(&rigid);
    particles.reset(grid.domain());
    gasImpulse_.clear();
    gasAngularImpulse_.clear();
    softGasImpulse_.clear();
    gasForceMax_ = 0;
    held_.clear();
    if (preset_ == Preset::Hydro) {
        const AABB d = grid.domain();
        // Bodies first (the water is not placed inside them).
        rigid.addBox({0.05f, 0.3f, 0.1f}, Vector3(0.08f), Quaternion::fromAxisAngle({0, 1, 0}, 0.3f), 500.0f, {0.9f, 0.6f, 0.25f});
        rigid.addBox({0.4f, 0.3f, -0.15f}, {0.12f, 0.03f, 0.08f}, Quaternion(), 600.0f, {0.75f, 0.5f, 0.3f}); // plank
        rigid.addSphere({0.7f, 0.45f, 0.15f}, 0.06f, 3000.0f, {0.6f, 0.6f, 0.65f});                          // sinks
        rigid.addCompound(teapotShape(), {-0.15f, 0.35f, -0.15f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.5f), 400.0f,
                          {0.8f, 0.55f, 0.85f});
        TriMesh foam = primitives::box(Vector3(0.06f));
        foam.translate({0.3f, 0.45f, 0.2f});
        particles.addSoftBody(foam, 150.0f, 0.4f, {0.3f, 0.75f, 0.95f});
        // A flag on a pole downwind: its left edge is fixed, the wind makes it flutter.
        ClothMaterial flag;
        flag.areaDensity = 0.15f;
        flag.bendCompliance = 1e-2f;
        particles.addCloth({0.65f, 1.0f, 0.0f}, {0.3f, 0, 0}, {0, -0.22f, 0}, flag, 64, {0.9f, 0.25f, 0.25f});
        // Water: a pool over the whole floor and a column at the upwind end that collapses.
        particles.addBlock(AABB(d.lo, {d.hi.x, 0.16f, d.hi.z}));
        particles.addBlock(AABB({d.lo.x, 0.16f, d.lo.z}, {d.lo.x + 0.4f, 0.6f, d.hi.z}));
    }
    if (preset_ == Preset::Magnetosphere) {
        // The magnet: a dipole m = 600 A m^2 pointing up, from its vector potential
        // A = mu0/4pi (m x r) / r^3 (so div B = 0 exactly on the grid). Its equatorial field
        // mu0 m / 4pi r^3 balances the flow's ram pressure, B^2 / 2mu0 = rho u^2, at r = 0.29 m
        // (0.37 m counting the doubling of the field by the magnetopause currents).
        const Vector3 centre = obstacle.position, moment(0.0f, 600.0f, 0.0f);
        grid.magnetic.setBackgroundFromPotential([centre, moment](const Vector3& x) { // current-free B0
            const Vector3 r = x - centre;
            const float d = std::max(length(r), 0.08f); // inside the magnet (a conductor): no singularity
            return cross(moment, r) * (1e-7f / (d * d * d));
        });
    }
    if (preset_ == Preset::Fire) {
        const float floor = grid.domain().lo.y;
        // A cotton curtain on a rod, its lower edge a hand above the burner flame.
        ClothMaterial cotton;
        cotton.areaDensity = 0.2f;
        cotton.flammable = true; // cellulose: Arrhenius pyrolysis (ClothMaterial defaults)
        particles.addCloth({-0.3f, 0.55f, 0.08f}, {0.6f, 0, 0}, {0, -1.1f, 0}, cotton, 16, {0.85f, 0.8f, 0.7f});
        // Bodies beside the fire: a crate and a teapot on the floor, a ball thrown into the plume,
        // a jelly cube falling next to the curtain.
        rigid.addBox({0.38f, floor + 0.1f, -0.25f}, Vector3(0.1f), Quaternion::fromAxisAngle({0, 1, 0}, 0.4f), 500.0f,
                     {0.65f, 0.45f, 0.25f});
        rigid.addCompound(teapotShape(), {-0.35f, floor + 0.14f, -0.25f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.6f), 500.0f,
                          {0.8f, 0.55f, 0.85f});
        const int ball = rigid.addSphere({0.45f, 0.2f, 0.0f}, 0.06f, 300.0f, {0.9f, 0.35f, 0.3f});
        rigid.bodies()[ball].vel = {-1.2f, 1.0f, 0.0f};
        TriMesh jelly = primitives::box(Vector3(0.06f));
        jelly.translate({-0.4f, 0.3f, 0.2f});
        particles.addSoftBody(jelly, 150.0f, 0.3f, {0.55f, 0.9f, 0.35f});
    }
    if (preset_ == Preset::GasSoftCloth) {
        const float floor = grid.domain().lo.y;
        // Rigid: a teapot, a box and a ball on the floor beside the hot source.
        rigid.addCompound(teapotShape(), {-0.45f, floor + 0.16f, 0.3f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.6f), 500.0f,
                          {0.8f, 0.55f, 0.85f});
        rigid.addBox({0.45f, floor + 0.12f, 0.35f}, Vector3(0.12f), Quaternion::fromAxisAngle({0, 1, 0}, 0.4f), 400.0f,
                     {0.3f, 0.7f, 0.9f});
        rigid.addSphere({0.4f, floor + 0.1f, -0.3f}, 0.1f, 500.0f, {0.9f, 0.35f, 0.3f});
        // Soft: a foam cube and a jelly ball dropped from above.
        TriMesh cube = primitives::box(Vector3(0.08f));
        cube.translate({-0.2f, 0.7f, 0.35f});
        particles.addSoftBody(cube, 150.0f, 0.4f, {0.3f, 0.75f, 0.95f});
        TriMesh ball = primitives::sphere(0.08f, 16, 8);
        ball.translate({0.35f, 0.9f, -0.25f});
        particles.addSoftBody(ball, 150.0f, 0.15f, {0.55f, 0.9f, 0.35f});
        // Cloth: a silk handkerchief high above the hot source: it floats down into the plume, which
        // holds it up (terminal speed ~0.7 m/s, the plume rises at 2-3 m/s) ...
        ClothMaterial silk;
        silk.areaDensity = 0.04f;
        silk.bendCompliance = 1e-2f;
        particles.addCloth({-0.2f, 0.45f, -0.2f}, {0.4f, 0, 0}, {0, 0, 0.4f}, silk, 0, {0.95f, 0.85f, 0.5f});
        // ... and a cotton curtain on a rod beside it, sewn from two panels.
        ClothMaterial cotton;
        cotton.seamColumns = {13}; // 0.8 m at 3 cm: 28 columns, the middle
        particles.addCloth({-0.45f, 1.05f, -0.45f}, {0.8f, 0, 0}, {0, -0.9f, 0}, cotton, 16, {0.75f, 0.3f, 0.35f});
    }
    if (preset_ == Preset::SmokeBodies) {
        // A few bodies released above the hot plume: they fall through the smoke, push it aside
        // and shed vortices; the rising gas acts back on them.
        rigid.addBox({-0.35f, 0.5f, 0.1f}, Vector3(0.12f), Quaternion::fromAxisAngle({1, 0, 1}, 0.5f), 400.0f, {0.3f, 0.7f, 0.9f});
        rigid.addBox({0.3f, 0.9f, -0.15f}, {0.18f, 0.05f, 0.14f}, Quaternion::fromAxisAngle({0, 0, 1}, 0.3f), 300.0f,
                     {0.9f, 0.6f, 0.25f});
        rigid.addSphere({0.05f, 0.2f, 0.0f}, 0.12f, 500.0f, {0.9f, 0.35f, 0.3f});
        rigid.addSphere({-0.2f, 1.0f, -0.3f}, 0.09f, 500.0f, {0.95f, 0.85f, 0.3f});
        rigid.addConvex(primitives::cylinder(0.1f, 0.3f, 16), {0.35f, 0.3f, 0.35f}, Quaternion::fromAxisAngle({1, 0, 0}, 1.2f),
                        500.0f, {0.6f, 0.85f, 0.4f});
        rigid.addCompound(teapotShape(), {-0.4f, 0.85f, 0.35f}, Quaternion::fromAxisAngle({0, 1, 0}, 0.8f), 500.0f,
                          {0.8f, 0.55f, 0.85f});
        rigid.addCompound(teapotShape(), {0.4f, 0.55f, -0.35f}, Quaternion::fromAxisAngle({1, 0, 0.3f}, 2.2f), 500.0f,
                          {0.95f, 0.75f, 0.45f});
        rigid.addCompound(bunnyShape(), {0.0f, 0.95f, 0.0f}, Quaternion::fromAxisAngle({0, 1, 0}, -0.6f), 500.0f,
                          {0.92f, 0.9f, 0.86f}); // right above the plume: falls through it
        // They hang still (the smoke flows around them) until the plume has risen, then drop through it.
        releaseTime_ = 1.5f;
        for (int i = 0; i < int(rigid.bodies().size()); ++i) {
            RigidBody& b = rigid.bodies()[i];
            held_.push_back({i, b.invMass, b.invInertiaLocal});
            b.invMass = 0;
            b.invInertiaLocal = Vector3(0.0f);
            b.updateInertia();
        }
    }
}

void Simulation::updateSurfaceLoads() {
    if (!grid.hasObstacle() || !obstacleMesh_ || obstacleMesh_->empty()) {
        surfaceLoads_ = SurfaceLoads();
        return;
    }
    // Moment about the body's reference point (for the wing: ~ the aerodynamic centre), per
    // unit of the body size (diameter / chord).
    surfaceLoads_ = computeSurfaceLoads(grid, *obstacleMesh_, obstacle.position, grid.referenceArea(),
                                        std::max(obstacle.size, 1e-3f));
}

std::vector<MovingSolid> Simulation::movingSolids() const {
    std::vector<MovingSolid> out;
    out.reserve(rigid.bodies().size());
    for (const RigidBody& b : rigid.bodies()) {
        MovingSolid m;
        m.bounds = b.worldBounds();
        m.position = b.pos;
        m.velocity = b.vel;
        m.angularVelocity = b.angVel;
        m.length = 2.0f * b.boundingRadius();
        std::shared_ptr<const ConvexShape> shape = b.shape;
        const Matrix3x3 Rt = b.rotation().transposed();
        const Vector3 p = b.pos;
        m.inside = [shape, Rt, p](const Vector3& x) { return shape->contains(Rt * (x - p)); };
        m.resting = b.sleeping; // has not moved: the gas reuses its cells
        out.push_back(std::move(m));
    }
    // Soft bodies (after the rigid ones): cells within one particle spacing of a particle are
    // inside; the body moves with the mean velocity of its particles.
    const float s = particles.spacing();
    for (const SoftBody& sb : particles.softBodies()) {
        MovingSolid m;
        std::vector<Vector3> pts;
        Vector3 c(0.0f), v(0.0f);
        for (int i : sb.particles) {
            pts.push_back(particles.positions()[i]);
            m.bounds.expand(particles.positions()[i]);
            c += particles.positions()[i];
            v += particles.velocities()[i];
        }
        m.bounds.lo -= Vector3(s);
        m.bounds.hi += Vector3(s);
        m.position = c / float(pts.size());
        m.velocity = v / float(pts.size());
        m.length = maxComp(m.bounds.extent());
        const float r2 = s * s;
        m.inside = [pts = std::move(pts), r2](const Vector3& x) {
            for (const Vector3& p : pts)
                if (length2(p - x) < r2) return true;
            return false;
        };
        out.push_back(std::move(m));
    }
    return out;
}

void Simulation::applyGasDragOnCloth(float dt) {
    // Aerodynamic force of the gas on every cloth particle (its patch of area s^2): pressure drag
    // of a flat plate on the part of the relative wind normal to the sheet (Cd 1.2), skin friction
    // along it (Cf 0.02). Taken implicitly (the particle cannot overtake the wind in one step); the
    // reaction goes into the gas at the same point, so momentum is conserved.
    const float rho = grid.params.fluidDensity;
    const auto& x = particles.positions();
    const auto& v = particles.velocities();
    const auto& w = particles.invMasses();
    for (const Cloth& c : particles.cloths()) {
        const float area = c.spacing * c.spacing;
        for (int y = 0; y < c.height; ++y)
            for (int xg = 0; xg < c.width; ++xg) {
                const int i = c.particle(xg, y);
                if (w[i] == 0) continue;
                const Vector3 du = x[c.particle(std::min(xg + 1, c.width - 1), y)] - x[c.particle(std::max(xg - 1, 0), y)];
                const Vector3 dvv = x[c.particle(xg, std::min(y + 1, c.height - 1))] - x[c.particle(xg, std::max(y - 1, 0))];
                const Vector3 n = normalize(cross(du, dvv));
                const Vector3 rel = grid.velocityAt(x[i]) - v[i];
                const float wn = dot(rel, n);
                const Vector3 relT = rel - n * wn;
                const float kn = 0.5f * rho * 1.2f * area * std::fabs(wn) * dt * w[i];
                const float kt = 0.5f * rho * 0.02f * area * length(relT) * dt * w[i];
                const Vector3 dv = n * (wn * kn / (1.0f + kn)) + relT * (kt / (1.0f + kt));
                particles.addVelocity(i, dv);
                grid.addImpulse(x[i], -dv / w[i]);
            }
    }
}

void Simulation::applyGasDragOnLiquid(float dt) {
    // Drag of the air on the liquid particles that touch air (the surface, spray, drops): a sphere
    // of the particle's radius (Cd 0.47) in the relative wind, implicit; the reaction goes into the
    // air. The wind blows spray away and drags the surface; flying drops are slowed down.
    if (particles.fluidCount() == 0) return;
    const float rho = grid.params.fluidDensity, r = particles.params.particleRadius;
    const float area = kPi * r * r;
    const auto& x = particles.positions();
    const auto& v = particles.velocities();
    const auto& w = particles.invMasses();
    const auto& phase = particles.phases();
    for (size_t i = 0; i < x.size(); ++i) {
        if (phase[i] != uint8_t(ParticlePhase::Fluid) || w[i] == 0 || !grid.touchesGas(x[i])) continue;
        const Vector3 rel = grid.fluidVelocityAt(x[i]) - v[i];
        const float k = 0.5f * rho * 0.47f * area * length(rel) * dt * w[i];
        const Vector3 dv = rel * (k / (1.0f + k));
        particles.addVelocity(int(i), dv);
        grid.addImpulse(x[i], -dv / w[i]);
    }
}

void Simulation::passLiquidToGas() {
    if (particles.fluidCount() == 0) {
        grid.setLiquid({}, {});
        return;
    }
    std::vector<Vector3> pos, vel;
    pos.reserve(particles.fluidCount());
    vel.reserve(particles.fluidCount());
    for (size_t i = 0; i < particles.size(); ++i)
        if (particles.phases()[i] == uint8_t(ParticlePhase::Fluid)) {
            pos.push_back(particles.positions()[i]);
            vel.push_back(particles.velocities()[i]);
        }
    grid.setLiquid(std::move(pos), std::move(vel));
}

void Simulation::stepParticlesInGas() {
    // Soft bodies: the gas pressure impulses of the last gas steps (plus the buoyancy of the
    // displaced gas) as a velocity change of the whole body.
    const auto& bodies = particles.softBodies();
    if (gasPushesBodies)
        for (size_t b = 0; b < bodies.size() && b < softGasImpulse_.size(); ++b) {
            float mass = 0;
            for (int i : bodies[b].particles) mass += particles.invMasses()[i] > 0 ? 1.0f / particles.invMasses()[i] : 0.0f;
            if (mass <= 0) continue;
            const float volume = float(bodies[b].particles.size()) * std::pow(particles.spacing(), 3.0f);
            const Vector3 J = softGasImpulse_[b] - particles.params.gravity * (grid.params.fluidDensity * volume * frameDt);
            for (int i : bodies[b].particles)
                if (particles.invMasses()[i] > 0) particles.addVelocity(i, J / mass);
        }
    // Rigid bodies and particles interleaved (as in the liquid mode: floating needs it).
    const int n = std::max(1, particles.params.substeps);
    const float dt = frameDt / n;
    const int rs = std::max(1, rigid.params.substeps / n);
    for (int s = 0; s < n; ++s) {
        if (gasPushesBodies) {
            applyGasDragOnCloth(dt);
            applyGasDragOnLiquid(dt);
        }
        for (int r = 0; r < rs; ++r) rigid.step(dt / rs);
        particles.step(dt);
    }
}

void Simulation::stepGasWithBodies() {
    // Weak (staggered) two-way coupling, one exchange per frame:
    //  1) gas -> bodies: pressure impulses integrated over the previous frame's gas steps, plus the
    //     buoyancy of the displaced gas (the Boussinesq grid carries no hydrostatic pressure);
    //  2) rigid substeps;
    //  3) bodies -> gas: the bodies at their new poses and velocities are the moving boundaries of
    //     the gas steps that cover the same frame time.
    if (!held_.empty() && time_ >= releaseTime_) {
        for (const Held& h : held_) {
            if (h.body >= int(rigid.bodies().size())) continue;
            RigidBody& b = rigid.bodies()[h.body];
            b.invMass = h.invMass;
            b.invInertiaLocal = h.invInertiaLocal;
            b.updateInertia();
            rigid.wake(h.body);
        }
        held_.clear();
    }
    const int nb = int(rigid.bodies().size());
    if (gasPushesBodies) {
        const float rho = grid.params.fluidDensity;
        for (int i = 0; i < nb; ++i) {
            RigidBody& b = rigid.bodies()[i];
            if (b.invMass == 0) continue;
            Vector3 J = i < int(gasImpulse_.size()) ? gasImpulse_[i] : Vector3(0.0f);
            Vector3 L = i < int(gasAngularImpulse_.size()) ? gasAngularImpulse_[i] : Vector3(0.0f);
            J -= rigid.params.gravity * (rho * b.shape->volume() * frameDt); // Archimedes in the gas
            // A resting body is not woken by pressure noise far below its sleep threshold.
            if (b.sleeping && length(J) * b.invMass < 0.5f * rigid.params.sleepLinear) continue;
            rigid.applyExternalWrench(i, J, L);
        }
    }
    if (particles.size() == 0) {
        const int n = std::max(1, rigid.params.substeps);
        for (int s = 0; s < n; ++s) rigid.step(frameDt / n);
    } else {
        stepParticlesInGas(); // rigid bodies interleaved with the particles
    }

    // Fire: the cloths take heat from the gas (or burn and give it heat and fuel gas).
    if (grid.combustion.enabled) {
        std::vector<FireOutput> fire;
        auto gasHeat = [this](const Vector3& x) { return GasHeat{grid.temperatureAt(x), grid.irradianceAt(x)}; };
        particles.burnCloths(frameDt, gasHeat, grid.combustion.ambientTemperature, fire);
        for (const FireOutput& f : fire) grid.addEmission({f.position, f.fuel, f.heat, 0.0f});
    }

    gasImpulse_.assign(rigid.bodies().size(), Vector3(0.0f));
    gasAngularImpulse_.assign(rigid.bodies().size(), Vector3(0.0f));
    softGasImpulse_.assign(particles.softBodies().size(), Vector3(0.0f));
    float remaining = frameDt;
    for (int it = 0; it < 8 && remaining > 1e-6f; ++it) {
        grid.setMovingSolids(movingSolids());
        passLiquidToGas();
        float dt = grid.step(remaining);
        const auto& F = grid.movingForces();
        const auto& T = grid.movingTorques();
        for (size_t i = 0; i < F.size(); ++i) {
            if (i < gasImpulse_.size()) {
                gasImpulse_[i] += F[i] * dt;
                gasAngularImpulse_[i] += T[i] * dt;
            } else if (i - gasImpulse_.size() < softGasImpulse_.size()) {
                softGasImpulse_[i - gasImpulse_.size()] += F[i] * dt;
            }
        }
        remaining -= dt;
        lastGridDt_ = dt;
    }
    gasForceMax_ = 0;
    for (const Vector3& J : gasImpulse_) gasForceMax_ = std::max(gasForceMax_, length(J) / frameDt);
}

void Simulation::setupRigidScene() {
    rigid.clear();
    rigidDomain_ = preset_ == Preset::RigidTower ? AABB({-3.0f, 0.0f, -3.0f}, {3.0f, 24.0f, 3.0f})
                                                 : AABB({-2.0f, 0.0f, -2.0f}, {2.0f, 5.0f, 2.0f});
    rigid.setDomain(rigidDomain_);
    rigid.setStaticMesh(&obstacleBVH_);
    particles.setStaticMesh(&obstacleBVH_); // soft bodies / cloth can be added in this mode too
    particles.setRigidWorld(&rigid);
    particles.reset(rigidDomain_);
    uint32_t seed = 7;
    auto rnd = [&]() { seed = seed * 1664525u + 1013904223u; return (seed >> 8) * (1.0f / 16777216.0f); };
    auto randColor = [&]() { return Vector3(0.35f + 0.6f * rnd(), 0.35f + 0.6f * rnd(), 0.35f + 0.6f * rnd()); };
    switch (preset_) {
    case Preset::RigidFalling:
        for (int i = 0; i < 60; ++i) {
            Vector3 p(-1.2f + 2.4f * rnd(), 2.0f + 2.5f * rnd(), -1.2f + 2.4f * rnd());
            if (i % 2 == 0)
                rigid.addBox(p, Vector3(0.08f + 0.12f * rnd(), 0.08f + 0.12f * rnd(), 0.08f + 0.12f * rnd()),
                             Quaternion::fromAxisAngle({rnd(), rnd(), rnd() + 0.1f}, 6.28f * rnd()), 800.0f, randColor());
            else
                rigid.addSphere(p, 0.08f + 0.12f * rnd(), 800.0f, randColor());
        }
        break;
    case Preset::RigidGranular: {
        const float r = 0.05f;
        for (int k = 0; k < 12; ++k)
            for (int j = 0; j < 14; ++j)
                for (int i = 0; i < 12; ++i) {
                    Vector3 p(-0.6f + i * 2.1f * r + 0.01f * rnd(), 2.0f + j * 2.1f * r, -0.6f + k * 2.1f * r + 0.01f * rnd());
                    float t = float(j) / 14.0f;
                    rigid.addSphere(p, r * (0.85f + 0.3f * rnd()), 2000.0f, {0.85f - 0.3f * t, 0.65f, 0.35f + 0.4f * t});
                }
        break;
    }
    case Preset::RigidConvex: {
        // Mixed convex polyhedra dropped on a slope: exercises SAT (box-box) and GJK/EPA with
        // perturbation manifolds (hull-hull, hull-box, hull-mesh).
        std::vector<TriMesh> kinds = {primitives::cylinder(0.14f, 0.3f, 12), primitives::cone(0.16f, 0.35f, 12),
                                      primitives::sphere(0.16f, 8, 5), primitives::ellipsoid({0.22f, 0.1f, 0.14f}, 10, 6)};
        TriMesh tet;
        tet.positions = {{0.2f, 0.2f, 0.2f}, {-0.2f, -0.2f, 0.2f}, {-0.2f, 0.2f, -0.2f}, {0.2f, -0.2f, -0.2f}};
        tet.triangles = {{0, 1, 2}, {0, 3, 1}, {0, 2, 3}, {1, 3, 2}};
        tet.orientOutward();
        kinds.push_back(tet);
        for (int i = 0; i < 45; ++i) {
            Vector3 p(-1.3f + 2.6f * rnd(), 1.5f + 3.0f * rnd(), -1.3f + 2.6f * rnd());
            Quaternion q = Quaternion::fromAxisAngle({rnd() - 0.5f, rnd() - 0.5f, rnd() - 0.5f + 1e-3f}, 6.28f * rnd());
            if (i % 6 == 5)
                rigid.addBox(p, Vector3(0.1f + 0.1f * rnd(), 0.08f + 0.08f * rnd(), 0.1f + 0.1f * rnd()), q, 700.0f, randColor());
            else
                rigid.addConvex(kinds[i % kinds.size()], p, q, 700.0f, randColor());
        }
        break;
    }
    case Preset::RigidCcd: {
        // A 2 cm thick wall and bullets that travel 25 cm per substep: without CCD they tunnel.
        rigid.addBox({0.5f, 1.5f, 0.0f}, {0.01f, 1.5f, 1.5f}, Quaternion(), 0.0f, {0.6f, 0.62f, 0.66f});
        for (int i = 0; i < 12; ++i) {
            int b = rigid.addSphere({-1.8f, 0.4f + 0.2f * i, -1.0f + 0.18f * i}, 0.03f, 8000.0f, {1.0f, 0.9f, 0.3f});
            rigid.bodies()[b].vel = {150.0f, 5.0f, 0.0f};
        }
        // Moving targets: a free stack of boxes and thin plates hanging on hinges - bullets must hit
        // them (momentum transfer, bounce), never pass through or overlap.
        for (int i = 0; i < 6; ++i)
            rigid.addBox({-0.3f, 0.1f + 0.2f * i, -0.6f}, Vector3(0.1f), Quaternion(), 400.0f, {0.35f + 0.1f * i, 0.55f, 0.9f - 0.1f * i});
        for (int i = 0; i < 3; ++i) {
            float z = 0.2f + 0.3f * i;
            int plate = rigid.addBox({0.0f, 1.9f, z}, {0.01f, 0.4f, 0.12f}, Quaternion(), 1500.0f, {0.95f, 0.5f, 0.3f});
            rigid.addHingeJoint(plate, -1, {0.0f, 2.3f, z}, {0, 0, 1});
            int b = rigid.addSphere({-1.8f, 1.75f, z}, 0.03f, 8000.0f, {1.0f, 0.95f, 0.4f});
            rigid.bodies()[b].vel = {120.0f, 0.0f, 0.0f};
        }
        for (int i = 0; i < 6; ++i) {
            int b = rigid.addSphere({-1.8f, 0.1f + 0.2f * i, -0.6f}, 0.03f, 8000.0f, {1.0f, 0.95f, 0.4f});
            rigid.bodies()[b].vel = {150.0f, 0.0f, 0.0f};
        }
        // A thin fast-spinning plate next to a thin post.
        rigid.addBox({-1.0f, 1.5f, 1.2f}, {0.01f, 1.5f, 0.01f}, Quaternion(), 0.0f, {0.6f, 0.62f, 0.66f});
        int plate = rigid.addBox({-1.0f, 2.2f, 0.6f}, {0.02f, 0.02f, 0.5f}, Quaternion(), 800.0f, {0.4f, 0.8f, 1.0f});
        rigid.bodies()[plate].angVel = {0.0f, 120.0f, 0.0f};
        break;
    }
    case Preset::RigidJoints: {
        const Vector3 grey(0.7f, 0.72f, 0.76f);
        // 1) Ball joints: a chain of 10 links hanging from the ceiling, released horizontally.
        {
            Vector3 top(-1.6f, 4.4f, -1.2f);
            int prev = -1;
            for (int i = 0; i < 10; ++i) {
                Vector3 c = top + Vector3(0.12f + 0.24f * i, 0, 0);
                int id = rigid.addBox(c, {0.1f, 0.03f, 0.03f}, Quaternion(), 800.0f, {0.9f, 0.45f, 0.3f});
                rigid.addBallJoint(id, prev, c - Vector3(0.12f, 0, 0));
                prev = id;
            }
        }
        // 2) Hinges: a motor-driven paddle wheel and a swinging door with angle limits.
        {
            int paddle = rigid.addBox({0.0f, 1.0f, -1.2f}, {0.6f, 0.05f, 0.15f}, Quaternion(), 600.0f, {0.3f, 0.6f, 0.95f});
            HingeJoint& h = rigid.addHingeJoint(paddle, -1, {0.0f, 1.0f, -1.2f}, {0, 0, 1});
            h.motorEnabled = true;
            h.motorSpeed = 1.5f;
            h.maxMotorTorque = 400.0f;
            int door = rigid.addBox({-1.2f, 1.0f, 0.3f}, {0.35f, 0.6f, 0.03f}, Quaternion(), 400.0f, {0.95f, 0.8f, 0.3f});
            HingeJoint& d = rigid.addHingeJoint(door, -1, {-1.55f, 1.0f, 0.3f}, {0, 1, 0});
            d.limitEnabled = true;
            d.lower = -1.2f;
            d.upper = 1.2f;
            rigid.bodies()[door].angVel = {0, 4.0f, 0};
        }
        // 3) Slider: a carriage on a tilted rail with travel limits.
        {
            Vector3 axis = normalize(Vector3(1.0f, -0.4f, 0.0f));
            int cart = rigid.addBox({1.0f, 2.2f, 1.2f}, {0.15f, 0.1f, 0.1f}, Quaternion::fromAxisAngle({0, 0, 1}, std::atan2(axis.y, axis.x)),
                                    700.0f, {0.4f, 0.85f, 0.4f});
            SliderJoint& sl = rigid.addSliderJoint(cart, -1, axis);
            sl.limitEnabled = true;
            sl.lower = -0.9f;
            sl.upper = 0.9f;
        }
        // 4) Fixed: an L-shaped body welded from two boxes, dropped with a spin.
        {
            int a = rigid.addBox({1.2f, 3.0f, -0.3f}, {0.3f, 0.07f, 0.07f}, Quaternion(), 600.0f, {0.8f, 0.4f, 0.9f});
            int b = rigid.addBox({0.97f, 3.23f, -0.3f}, {0.07f, 0.3f, 0.07f}, Quaternion(), 600.0f, {0.8f, 0.4f, 0.9f});
            rigid.addFixedJoint(a, b);
            rigid.bodies()[a].angVel = {1.0f, 0.5f, 2.0f};
            rigid.bodies()[b].angVel = {1.0f, 0.5f, 2.0f};
        }
        // 5) Distance: a rope pendulum and a soft spring.
        {
            int bob = rigid.addSphere({1.5f, 3.4f, -1.4f}, 0.12f, 2000.0f, {0.95f, 0.3f, 0.35f});
            DistanceJoint& rope = rigid.addDistanceJoint(bob, -1, {1.5f, 3.4f, -1.4f}, {0.6f, 4.4f, -1.4f});
            rope.rope = true;
            int weight = rigid.addBox({0.4f, 3.0f, 1.3f}, Vector3(0.12f), Quaternion(), 900.0f, {0.3f, 0.8f, 0.85f});
            DistanceJoint& spring = rigid.addDistanceJoint(weight, -1, {0.4f, 3.12f, 1.3f}, {0.4f, 4.4f, 1.3f});
            spring.frequency = 1.2f;
            spring.dampingRatio = 0.1f;
        }
        // A few free bodies to hit with the paddle or grab with the mouse.
        for (int i = 0; i < 6; ++i)
            rigid.addBox({-0.3f + 0.12f * i, 2.0f + 0.3f * i, -1.2f}, Vector3(0.08f), Quaternion(), 500.0f, grey);
        break;
    }
    case Preset::RigidTeapots: {
        // Non-convex bodies: all 100 share one compound shape (see teapotShape()).
        const auto teapot = teapotShape();
        for (int i = 0; i < 100; ++i) {
            int layer = i / 25, cell = i % 25;
            Vector3 p(-1.2f + 0.6f * (cell % 5) + 0.08f * (rnd() - 0.5f), 0.4f + 0.5f * layer,
                   -1.2f + 0.6f * (cell / 5) + 0.08f * (rnd() - 0.5f));
            Quaternion q = Quaternion::fromAxisAngle({rnd() - 0.5f, rnd() - 0.5f, rnd() - 0.5f + 1e-3f}, 6.28f * rnd());
            float t = i / 99.0f;
            rigid.addCompound(teapot, p, q, 600.0f, {0.85f - 0.4f * t, 0.45f + 0.35f * std::sin(3.14f * t), 0.3f + 0.6f * t});
        }
        break;
    }
    case Preset::RigidTower: {
        // Standard stability test: 100 cubes, each released 1 cm above the one below.
        const float h = 0.1f, gap = 0.01f;
        for (int i = 0; i < 100; ++i) {
            float t = i / 99.0f;
            rigid.addBox({0.0f, h + i * (2 * h + gap), 0.0f}, Vector3(h), Quaternion(), 500.0f,
                         {0.25f + 0.65f * t, 0.55f + 0.2f * std::sin(6.28f * t), 0.9f - 0.6f * t});
        }
        break;
    }
    case Preset::RigidPyramid: {
        const float h = 0.15f;
        int levels = 7;
        for (int l = 0; l < levels; ++l)
            for (int i = 0; i < levels - l; ++i) {
                float x = (i - (levels - l - 1) * 0.5f) * (2 * h + 0.005f);
                rigid.addBox({x, h + l * 2 * h, 0.0f}, Vector3(h), Quaternion(), 500.0f,
                             {0.9f - 0.08f * l, 0.5f + 0.05f * l, 0.3f});
            }
        int ball = rigid.addSphere({-1.8f, 0.6f, 0.0f}, 0.2f, 3000.0f, {0.2f, 0.3f, 0.9f});
        rigid.bodies()[ball].vel = {9.0f, 1.5f, 0.0f};
        break;
    }
    default: break;
    }
}

// ---------------------------------------------------------------------------
// Stepping
// ---------------------------------------------------------------------------
void Simulation::stepFrame() {
    auto t0 = std::chrono::steady_clock::now();
    switch (mode_) {
    case SimMode::Fluid: {
        int n = std::max(1, particles.params.substeps);
        float dt = frameDt / n;
        const int rs = std::max(1, rigid.params.substeps / n); // rigid substeps per SPH substep
        for (int s = 0; s < n; ++s) {
            for (int r = 0; r < rs; ++r) rigid.step(dt / rs);
            particles.step(dt);
        }
        time_ += frameDt;
        break;
    }
    case SimMode::WindTunnel:
        if (rigid.bodies().empty() && particles.size() == 0) { // gas alone
            grid.setMovingSolids({});
            lastGridDt_ = grid.step(0.05f);
            // Plasma: the Alfven waves keep the steps short - take a few per frame (up to 1/60 s).
            for (int extra = 0; extra < 3 && grid.magnetic.enabled && grid.time() - time_ < frameDt; ++extra)
                lastGridDt_ = grid.step(frameDt - (grid.time() - time_));
            time_ = grid.time();
        } else {
            stepGasWithBodies();
            time_ += frameDt;
        }
        updateSurfaceLoads();
        break;
    case SimMode::Rigid: {
        if (particles.size() == 0) {
            int n = std::max(1, rigid.params.substeps);
            float dt = frameDt / n;
            for (int s = 0; s < n; ++s) rigid.step(dt);
        } else { // soft bodies / cloth added: interleave as in the liquid mode
            int n = std::max(1, particles.params.substeps);
            float dt = frameDt / n;
            const int rs = std::max(1, rigid.params.substeps / n);
            for (int s = 0; s < n; ++s) {
                for (int r = 0; r < rs; ++r) rigid.step(dt / rs);
                particles.step(dt);
            }
        }
        time_ += frameDt;
        break;
    }
    }
    ++frame_;
    lastStepMs_ = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

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

int Simulation::sliceLayer() const {
    int n[3] = {grid.nx(), grid.ny(), grid.nz()};
    int axis = clampv(vis.sliceAxis, 0, 2);
    return clampv(int(vis.slicePosition * n[axis]), 0, n[axis] - 1);
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

void Simulation::computeFieldLines(RenderSnapshot& s) const {
    // Magnetic field lines: traced both ways from seeds on a shell around the obstacle (the magnet)
    // or on a lattice through the domain, along B / |B| (midpoint rule), until they leave the
    // domain, enter a solid or the field vanishes.
    const MagneticField& m = grid.magnetic;
    const float dx = grid.dx();
    const AABB dom = grid.domain();
    std::vector<Vector3> seeds;
    if (grid.hasObstacle() && obstacleMesh_ && !obstacleMesh_->empty()) {
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
    auto blocked = [&](const Vector3& p) {
        if (!dom.contains(p)) return true;
        const Vector3 g = (p - dom.lo) / dx;
        const int i = int(g.x), j = int(g.y), k = int(g.z);
        return i >= 0 && j >= 0 && k >= 0 && i < grid.nx() && j < grid.ny() && k < grid.nz() && grid.solid(i, j, k);
    };
    const float h = 0.5f * dx;
    const int maxSteps = 3 * (grid.nx() + grid.ny() + grid.nz());
    s.fieldLines.clear();
    s.fieldLineStrength.clear();
    float bmax = 1e-12f;
    for (const Vector3& seed : seeds) {
        std::vector<Vector3> line;
        std::vector<float> strength;
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
        if (line.size() < 2) continue;
        for (float B : strength) bmax = std::max(bmax, B);
        s.fieldLines.push_back(std::move(line));
        s.fieldLineStrength.push_back(std::move(strength));
    }
    s.fieldLineMax = bmax;
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

void Simulation::fillSnapshot(RenderSnapshot& s) const {
    s.frame = frame_;
    s.mode = mode_;
    s.time = time_;
    s.stepMs = lastStepMs_;
    s.preset = preset_;
    s.paramsVersion = paramsVersion_;
    s.particleParams = particles.params;
    s.ns = grid.params;
    s.rigid = rigid.params;
    s.obstacleSettings = obstacle;
    s.vis = vis;
    s.emitter = particles.emitter;
    s.heat = grid.source;
    s.gasPushesBodies = gasPushesBodies;
    s.obstacleVersion = obstacleVersion_;
    s.obstacle = obstacleMesh_;
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

    for (const RigidBody& b : rigid.bodies())
    {
        std::shared_ptr<const TriMesh> mesh;
        if (b.type() == ShapeType::ConvexHull) mesh = static_cast<const ConvexHullShape*>(b.shape.get())->mesh();
        else if (b.type() == ShapeType::Compound) mesh = static_cast<const CompoundShape*>(b.shape.get())->visualMesh();
        s.bodies.push_back({b.type(), b.pos, b.rot, b.halfExtents(), b.radius(), b.color, mesh, b.sleeping, b.shape});
    }

    auto range = [&](float lo, float hi) {
        if (vis.autoRange) { s.colorMin = lo; s.colorMax = hi > lo ? hi : lo + 1e-3f; }
        else { s.colorMin = vis.rangeMin; s.colorMax = vis.rangeMax > vis.rangeMin ? vis.rangeMax : vis.rangeMin + 1e-3f; }
    };

    switch (mode_) {
    case SimMode::Fluid: {
        s.domain = particles.domain();
        s.particleRadius = particles.params.particleRadius;
        const auto& x = particles.positions();
        const auto& v = particles.velocities();
        const auto& rho = particles.densities();
        const auto& phase = particles.phases();
        const bool anySolid = particles.fluidCount() < particles.size();
        s.particles.clear();
        s.particles.reserve(x.size());
        s.particleScalar.clear();
        s.particleScalar.reserve(x.size());
        if (anySolid) s.particleColor.reserve(x.size());
        float lo = kInf, hi = -kInf;
        for (size_t i = 0; i < x.size(); ++i) {
            if (phase[i] != uint8_t(ParticlePhase::Fluid)) continue; // cloth / soft bodies: drawn as surfaces
            if (vis.liquidSurface) {
                s.liquid.push_back(x[i]);
                continue;
            }
            s.particles.push_back(x[i]);
            if (anySolid) s.particleColor.push_back(particles.particleColors()[i]);
            if (phase[i] == uint8_t(ParticlePhase::Soft)) {
                s.particleScalar.push_back(RenderSnapshot::kOwnColor);
                continue;
            }
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
        range(lo, hi);
        s.info.push_back({"Частиц", fmt("%zu", particles.size())});
        s.info.push_back({"Радиус частицы", fmt("%.1f мм", particles.params.particleRadius * 1000)});
        s.info.push_back({"Масса частицы", fmt("%.3g кг", particles.particleMass())});
        s.info.push_back({"Ошибка плотности", fmt("%.2f %%", particles.averageDensityError() * 100)});
        s.info.push_back({"Макс. скорость", fmt("%.2f м/с", particles.maxSpeed())});
        s.info.push_back({"Тел", fmt("%zu", rigid.bodies().size())});
        if (!particles.softBodies().empty() || !particles.cloths().empty()) {
            s.info.push_back({"Мягких тел / тканей", fmt("%zu / %zu", particles.softBodies().size(), particles.cloths().size())});
            s.info.push_back({"Частиц жидкости / твёрдых", fmt("%zu / %zu", particles.fluidCount(), particles.size() - particles.fluidCount())});
            s.info.push_back({"Контактов частиц", fmt("%zu", particles.particleContactCount())});
        }
        s.plots.push_back({"Ошибка плотности, %", particles.averageDensityError() * 100});
        s.plots.push_back({"Макс. скорость, м/с", particles.maxSpeed()});
        break;
    }
    case SimMode::WindTunnel: {
        s.domain = grid.domain();
        if (vis.showSlice) extractSlice(s);
        float lo = kInf, hi = -kInf;
        for (size_t i = 0; i < s.slice.size(); ++i)
            if (!s.sliceSolid[i]) { lo = std::min(lo, s.slice[i]); hi = std::max(hi, s.slice[i]); }
        if (s.slice.empty()) { lo = 0; hi = 1; }
        if (vis.sliceField == GridField::Speed || vis.sliceField == GridField::Vorticity) lo = 0;
        range(lo, hi);
        s.colorLabel = fieldLabel(vis.sliceField);

        if (vis.showSmoke) {
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
        if (vis.showStreamlines) computeStreamlines(s);
        if (grid.magnetic.enabled && vis.showFieldLines) computeFieldLines(s);
        if (grid.magnetic.enabled) {
            const MagneticField& m = grid.magnetic;
            const float B = m.maxField(), rho = grid.params.fluidDensity;
            const float L = grid.hasObstacle() ? std::max(obstacle.size, grid.dx()) : grid.domain().extent().x;
            s.info.push_back({"Магнитное поле, макс.", fmt("%.1f мТл", B * 1000)});
            s.info.push_back({"Скорость Альфвена v_A", fmt("%.2f м/с (предел Бориса %.1f)", B / std::sqrt(MagneticField::kMu0 * rho),
                                                         m.speedLimit)});
            s.info.push_back({"Маг. число Рейнольдса Rm", fmt("%.3g", grid.params.inflowSpeed * L / m.resistivity())});
            s.info.push_back({"Энергия поля", fmt("%.3g Дж", m.energy())});
            s.info.push_back({"div B (отн.)", fmt("%.1e", m.maxDivergence())});
        }
        s.gridNx = grid.nx();
        s.gridNy = grid.ny();
        s.gridNz = grid.nz();
        s.gridDx = grid.dx();
        s.gridOrigin = grid.origin();
        s.sliceLayer = sliceLayer();
        if (vis.vectorDisplay > 0) extractVectors(s);
        if (vis.surfacePressure && obstacleMesh_ && !obstacleMesh_->empty() &&
            surfaceLoads_.triangles.size() == obstacleMesh_->triangles.size()) {
            // Cp of the triangles, area-weighted onto the vertices for smooth colouring.
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
        float L = obstacle.size;
        float Re = grid.params.inflowSpeed * L / std::max(grid.params.kinematicViscosity, 1e-9f);
        s.info.push_back({"Сетка", fmt("%d × %d × %d  (%.1f тыс. ячеек)", grid.nx(), grid.ny(), grid.nz(),
                                       grid.nx() * grid.ny() * grid.nz() / 1000.0f)});
        s.info.push_back({"Шаг сетки dx", fmt("%.1f мм", grid.dx() * 1000)});
        s.info.push_back({"Шаг по времени", fmt("%.2f мс", lastGridDt_ * 1000)});
        if (grid.combustion.enabled) {
            float Tmax = 0;
            for (float t : grid.temperature().d) Tmax = std::max(Tmax, t);
            s.info.push_back({"Мощность пламени", fmt("%.1f кВт", grid.heatReleaseRate() / 1000)});
            s.info.push_back({"Макс. температура", fmt("%.0f K (%.0f °C)", Tmax + grid.combustion.ambientTemperature, Tmax + grid.combustion.ambientTemperature - 273.15f)});
            int burnt = 0;
            for (const Cloth& c : particles.cloths()) burnt += c.burntThreads;
            if (!particles.cloths().empty()) s.info.push_back({"Прогоревших нитей", fmt("%d", burnt)});
        }
        if (grid.hasObstacle()) s.info.push_back({"Число Рейнольдса", fmt("%.3g", Re)});
        if (grid.hasObstacle()) {
            Vector3 F = grid.bodyForce();
            s.info.push_back({"Cd (сред.)", fmt("%.3f", grid.dragCoefficientAvg())});
            s.info.push_back({"Cl (сред.)", fmt("%.3f", grid.liftCoefficientAvg())});
            s.info.push_back({"Cd / Cl мгн.", fmt("%.3f / %.3f", grid.dragCoefficient(), grid.liftCoefficient())});
            if (std::fabs(grid.dragCoefficientAvg()) > 1e-4f)
                s.info.push_back({"L/D", fmt("%.2f", grid.liftCoefficientAvg() / grid.dragCoefficientAvg())});
            s.info.push_back({"Сила F", fmt("(%.2f, %.2f, %.2f) Н", F.x, F.y, F.z)});
            if (grid.params.wallFriction) {
                float cdf = grid.frictionForce().x / (grid.dynamicPressure() * std::max(grid.referenceArea(), 1e-9f));
                s.info.push_back({"из них трение (Cd тр.)", fmt("%.3f", cdf)});
            }
            const SurfaceLoads& SL = surfaceLoads_;
            if (!SL.triangles.empty()) {
                s.info.push_back({"По полигонам: Cd", fmt("%.3f (давл. %.3f + трение %.3f)", SL.cd, SL.cdPressure, SL.cdFriction)});
                s.info.push_back({"По полигонам: Cl / Cm", fmt("%.3f / %.3f", SL.cl, SL.cm)});
                s.info.push_back({"Треугольников, смоч. площадь", fmt("%zu, %.4f м²", SL.triangles.size(), SL.wettedArea)});
            }
            s.info.push_back({grid.params.usePlanformArea ? "Площадь (в плане)" : "Площадь (миделя)",
                              fmt("%.4f м²", grid.referenceArea())});
            s.plots.push_back({"Cd", grid.dragCoefficient()});
            s.plots.push_back({"Cl", grid.liftCoefficient()});
        }
        s.info.push_back({"Итераций давления", fmt("%d (невязка %.1e)", grid.lastPressureIterations(), grid.lastResidual())});
        s.info.push_back({"Макс. скорость", fmt("%.2f м/с", grid.maxVelocity())});
        s.info.push_back({"Макс. |div u| после проекции", fmt("%.2e 1/с", grid.maxDivergence())});
        s.info.push_back({"Дым в объёме", fmt("%.4f м³", grid.totalSmoke())});
        if (particles.fluidCount() > 0)
            s.info.push_back({"Вода: частиц / ячеек в воздухе", fmt("%zu / %d", particles.fluidCount(), grid.liquidCellCount())});
        if (particles.hasSolids()) {
            int torn = 0;
            for (const Cloth& c : particles.cloths()) torn += c.tornThreads;
            s.info.push_back({"Мягких тел / тканей", fmt("%zu / %zu (порвано нитей %d)", particles.softBodies().size(),
                                                             particles.cloths().size(), torn)});
        }
        if (!rigid.bodies().empty()) {
            s.info.push_back({"Тел в газе", fmt("%zu (%d ячеек, спят %zu)", rigid.bodies().size(), grid.movingSolidCells(),
                                                rigid.sleepingCount())});
            s.info.push_back({"Сила газа на тело (макс.)", fmt("%.3f Н", gasForceMax_)});
            if (!held_.empty()) s.info.push_back({"Тела отпустятся через", fmt("%.1f с", std::max(0.0f, releaseTime_ - time_))});
        }
        if (!s.arrowPos.empty()) s.info.push_back({"Векторов скорости", fmt("%zu", s.arrowPos.size())});
        if (!grid.hasObstacle()) {
            s.plots.push_back({"Макс. скорость, м/с", grid.maxVelocity()});
            s.plots.push_back({"Дым, дм³", grid.totalSmoke() * 1000.0f});
        }
        break;
    }
    case SimMode::Rigid: {
        s.domain = rigidDomain_;
        s.colorMin = 0;
        s.colorMax = 1;
        s.info.push_back({"Тел", fmt("%zu", rigid.bodies().size())});
        s.info.push_back({"Контактов", fmt("%zu", rigid.contactCount())});
        s.info.push_back({"Спящих тел", fmt("%zu", rigid.sleepingCount())});
        s.info.push_back({"CCD: остановлено тел за шаг", fmt("%zu", rigid.ccdHits())});
        s.info.push_back({"Сочленений", fmt("%zu", rigid.joints().size())});
        s.info.push_back({"Кин. энергия", fmt("%.2f Дж", rigid.kineticEnergy())});
        s.plots.push_back({"Кин. энергия, Дж", rigid.kineticEnergy()});
        break;
    }
    }
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
    s.info.insert(s.info.begin(), {"Время", fmt("%.3f с", time_)});
    s.info.push_back({"Шаг расчёта", fmt("%.1f мс", lastStepMs_)});
}

} // namespace rf
