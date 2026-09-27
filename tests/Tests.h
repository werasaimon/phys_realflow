#pragma once
// Every test, by module. Each is a numerical experiment with an exact answer or a hard physical
// criterion; see the files for what they check.

// core utilities: the probe and the allocation census (CoreTests.cpp)
void testProbe();
void testAllocationsPerFrame();

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
void testNewtonCradle();        // RigidTests.cpp: a row of balls passes a hit on (restitution through the slop zone)
void testNoetherRigid();        // RigidTests.cpp: energy, momentum, angular momentum of the rigid solver
void testGridConvergence();     // GasTests.cpp: order of convergence of the gas solver on refined grids
void testDamBreakMartinMoyce(); // ParticleTests.cpp: dam-break front against Martin & Moyce 1952
void testVoronoiFracture();     // FractureTests.cpp: Voronoi cells of a convex body fill it, convex and watertight
void testTerrain();             // TerrainTests.cpp: 150 bodies on a 51 200-triangle static mesh (BVH)
void testParticlesManyBodies(); // ParticleTests.cpp: 30 000 particles find 150 bodies through the world tree
void testDestroyBody();         // RigidTests.cpp: one body removed, its slot reused, the rest untouched
void testCapsuleShape();        // RigidTests.cpp: capsule mass properties, lying / toppling on a box, raycasts
void testRemoveParticleGroup(); // ParticleTests.cpp: one soft body / liquid removed, the rest goes on

// relativity (RelativityTests.cpp): geodesics in Kerr, light bending, the shadow of a black hole
void testGeodesicInvariants();
void testLightDeflection();
void testPhotonSphere();
void testPerihelionPrecession();
void testShadow();
void testHorizonPenetration();

// scene graph of the editor (GraphTests.cpp): text file round trip, magnet forces, a graph scene runs
void testSceneGraphRoundTrip();
void testMagnetForce();
void testMagnetsAttract();
void testGraphSceneBuilds();
void testGraphCloth();
void testGraphEmitterFollows();
void testGraphFlammableCloth();
void testGraphMeshShape();
void testGraphGeometryOnly();
void testGraphColliderApart();  // GraphTests.cpp: a sphere on a box collider slides, a box on a sphere collider rolls
void testGraphColliderFit();    // GraphTests.cpp: a fitted capsule under a model, collider round trip, off-centre rebuild
void testGraphColliderOnly();   // GraphTests.cpp: a collider without the rigid role is a static obstacle
// meta-objects (MetaObjectTests.cpp): one entity changes between frames, the rest does not notice
void testMetaWaterSoftCycle();
void testMetaRigidToClothAndBack();
void testMetaMagnetToggle();
void testMetaNoGrowth();
