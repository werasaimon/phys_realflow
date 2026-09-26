// The scenes of Simulation: their names and modes, the parameters every preset sets, the
// obstacle geometry, and the bodies, particles, cloth and fields each scene starts with.
#include "scene/Simulation.h"

#include "rigid/ConvexDecomposition.h"

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
    case Preset::Tokamak: return "Плазма: токамак (кольцо плазмы в тороидальном поле)";
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
    case Preset::GasSoftCloth: case Preset::Hydro: case Preset::Fire: case Preset::Magnetosphere: case Preset::Tokamak:
        return SimMode::WindTunnel;
    default:
        return SimMode::Rigid;
    }
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
    grid.params = GasParams();
    grid.source = HeatSource();
    grid.combustion = Combustion();
    grid.magnetic = MagneticField();
    grid.vessel = nullptr;
    tokamak = Tokamak();
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
    case Preset::Tokamak: {
        // A ring of plasma in a toroidal vessel (Tokamak.h): the coils' toroidal field, the plasma
        // current's poloidal field and Shafranov's vertical field hold it; the field lines wind
        // with the safety factor q. With q_a = 0.7 (default) the column kinks into a helix within
        // a couple of seconds; set tokamak.safetyFactorEdge above 1 (or below the wall limit 0.4)
        // and reset to see it stay a ring. Scaled to millitesla: the Alfven speed is a few m/s, so
        // the kink takes seconds instead of microseconds.
        const Tokamak& t = tokamak; // R0 = 0.4 m, a = 0.12 m, b = 0.24 m, B0 = 2.5 mT, q_a = 0.7
        obstacle.shape = ObstacleShape::Custom;
        obstacle.customMesh = std::make_shared<TriMesh>(primitives::torus(t.majorRadius, t.vesselRadius));
        obstacle.customName = "торовый сосуд";
        obstacle.size = 2.0f * (t.majorRadius + t.vesselRadius); // the mesh as built
        grid.params.domainSize = {1.4f, 0.6f, 1.4f};
        grid.params.resolutionX = 56; // dx = 2.5 cm: 5 cells across the current channel's radius
        grid.params.inflowSpeed = 0.0f;
        for (BoundaryType& b : grid.params.bc) b = BoundaryType::Wall;
        grid.params.fluidDensity = 1.0f;
        grid.params.pressureTolerance = 1e-3f;
        grid.params.smokeRake = false;
        grid.params.wallFriction = false; // no boundary layer of a gas: the plasma slips along the wall
        grid.magnetic.enabled = true;
        grid.magnetic.conductivity = 1e10f;      // the current outlives the scene: tau = a^2 / (5.8 eta) ~ 30 s
        grid.magnetic.numericalDissipation = 0.0f; // only the Lax-Wendroff amount: f |u| dx would eat the current
        vis.sliceField = GridField::CurrentDensity;
        vis.sliceAxis = 2;
        vis.slicePosition = 0.5f; // the poloidal cross-section at z = 0
        vis.showSlice = false;
        vis.showStreamlines = false;
        vis.surfacePressure = false;
        vis.gridDisplay = 0;
        vis.vesselGlass = true;
        break;
    }
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
        Quaternion q = Quaternion::fromEuler(degToRad(obstacle.yawDeg), -degToRad(obstacle.angleOfAttackDeg), degToRad(obstacle.rollDeg));
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
        preset_ == Preset::GasSoftCloth || preset_ == Preset::Fire || preset_ == Preset::Tokamak)
        origin = Vector3(-0.5f * size.x, -0.5f * size.y, -0.5f * size.z);
    else if (preset_ == Preset::Hydro)
        origin = Vector3(-0.5f * size.x, 0.0f, -0.5f * size.z); // floor at y = 0
    else origin = Vector3(-0.3f * size.x, -0.5f * size.y, -0.5f * size.z); // body at 30% of the length
    if (preset_ == Preset::Tokamak) {
        // The torus is a vessel: the plasma fills it and everything outside is its conducting
        // wall. The mesh is only drawn (as glass), not voxelised as a body.
        grid.vessel = [t = tokamak](const Vector3& x) { return t.inside(x); };
        grid.reset(origin, nullptr);
    } else grid.reset(origin, &obstacleBVH_);

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
    if (preset_ == Preset::Tokamak) {
        // The coils' fields (toroidal B0 R0 / R and the vertical field) are current-free: the
        // background B0. The plasma current's poloidal field is the evolving part B1. The glowing
        // plasma sits where the current flows.
        const Tokamak t = tokamak;
        grid.magnetic.setBackgroundFromPotential([t](const Vector3& x) { return t.coilPotential(x); });
        grid.magnetic.addFromPotential([t](const Vector3& x) { return t.plasmaPotential(x); });
        // The resistive "vacuum" between the channel and the wall (see Tokamak.h).
        const float etaPlasma = grid.magnetic.resistivity();
        grid.magnetic.setResistivityMap([t, etaPlasma](const Vector3& x) { return t.resistivityAt(x, etaPlasma); });
        grid.setTracer([t](const Vector3& x) { return t.tracer(x); });
        tokamakBv_ = t.verticalFieldStrength();
        tokamakShift_ = 0;
        tokamakControlTime_ = 0;
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
        for (int i = 0; i < int(rigid.bodies().size()); ++i) rigid.hold(i);
    }
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

} // namespace rf
