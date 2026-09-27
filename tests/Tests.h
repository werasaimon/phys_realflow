#pragma once
// Every test, by module. Each is a numerical experiment with an exact answer or a hard physical
// criterion; see the files for what they check.

// math, meshes, spatial structures (MathTests.cpp)
void testMath();
void testPrimitives();
void testBVH();
void testAABBTree();
void testMassProperties();

// rigid bodies (RigidTests.cpp)
void testRigid();
void testGjkEpa();
void testBoxStack();
void testBroadPhase();
void testTallStack();
void testStack100();
void testStack200();
void testRaycastGrab();
void testJoints();
void testCcd();
void testCcdBodies();
void testCcdSpinningPlate();
void testGjkRandomThin();
void testConvexHullAndDecomposition();
void testTeapots();
void testBeamOverCubes();
void testConvexRest();
void testHardContacts();
void testHardContactDynamics();

// particles: liquid, soft bodies, cloth (ParticleTests.cpp)
void testSoftBodyAndCloth();
void testSPH();
void testFloating();
void testLiquidWalls();
void testLightBodyInWater();

// gas and fire (GasTests.cpp)
void testGasBodies();
void testSurfaceLoads();
void testGasParticles();
void testHydro();
void testGridUniform();
void testGridSphere();
void testGridWingLift();
void testSmokeClosedBox();
void testDisturbance();
void testCombustion();
void testHeatConduction();
void testPyrolysis();
void testFireScene();

// plasma (PlasmaTests.cpp)
void testMagneticField();
void testMagnetosphere();
void testTokamak();

// scenes (SceneTests.cpp)
void testSimulationPresets();
void testCoherence();
void testDeterminism();
// validation against published benchmarks (docs/00-vision.md, "Мировые эталоны")
void testCylinderStrouhal();
void testNoetherRigid();        // RigidTests.cpp: energy, momentum, angular momentum of the rigid solver
void testGridConvergence();     // GasTests.cpp: order of convergence of the gas solver on refined grids
void testDamBreakMartinMoyce(); // ParticleTests.cpp: dam-break front against Martin & Moyce 1952
