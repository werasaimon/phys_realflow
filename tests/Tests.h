#pragma once
// Every test, by module. Each is a numerical experiment with an exact answer or a hard physical
// criterion; see the files for what they check.

// core utilities: the probe and the allocation census (CoreTests.cpp)
void testProbe();
void testAllocationsPerFrame();
void testTimersCoverTheStep(); // the stage timers add up to the whole step: no grey gap in the profiler

// the research layers of the debug drawing (DebugDrawTests.cpp)
void testDebugLayersOff();
void testDebugContactLayers();
void testDebugWatchedPair();
void testDebugWorldTree();
void testDebugClothTension();
void testDebugGasLayers();
void testDebugLayerCap();

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
void testRigidPerfThousandCubes(); // PerfTests.cpp: 1000 cubes fall, ms per stage
void testRaycastGrab();
void testJoints();
void testCcd();
void testCcdBodies();
void testCcdSpinningPlate();
void testGjkRandomThin();
void testConvexHullAndDecomposition();
void testTeapots();
void testNonConvexDecomposition();
void testRingOnFloor();
void testChainOfRings();
void testRingOverPeg();
void testDropIntoCup();
void testRollingRace();
void testBowlHoldsPour();
void testStackOfTables();
void testNestedCups();
void testModelsOnTerrain();
void testBeamOverCubes();
void testConvexRest();
void testHardContacts();
void testHardContactDynamics();

// particles: liquid, soft bodies, cloth (ParticleTests.cpp)
void testSoftBodyAndCloth();
void testSoftStackNoOverlap();
void testSoftPressedApartAndBox();
void testClothSwingNoFalseTears();
void testClothTearsAtStrength();
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
void testAlfvenicState();
void testTokamak();

// scenes (SceneTests.cpp)
void testSimulationPresets();
void testCoherence();
void testDeterminism();
// determinism across machines (DeterminismTests.cpp)
void testDeterminismAcrossThreads(); // the same bits on 1, 2 and all threads of the pool
void testFireAcrossThreads();        // the whole fire scene on 1, 2 and all threads (minutes; on request)
void testRigidGoldenHash();          // the rigid benchmark scene after 10 s: one known hash (RF_STRICT_FP)
void testElementaryFunctionsAreStd(); // rf::sin and friends return std::'s bits (step 1 of 2)
// validation against published benchmarks (docs/00-vision.md, "Мировые эталоны")
void testCylinderStrouhal();
void testNewtonCradle();        // RigidTests.cpp: a row of balls passes a hit on (restitution through the slop zone)
void testNoetherRigid();        // RigidTests.cpp: energy, momentum, angular momentum of the rigid solver
void testMultigridPressure();   // GasTests.cpp: MGPCG against the Jacobi PCG, iterations on 32^3 .. 128^3
void testAdvectionReflection(); // GasTests.cpp: the inviscid Taylor-Green energy, projection vs reflection
void testPressureBenchmark();   // GasTests.cpp: pressure ms per frame on the gas scenes (on request)
void testGridConvergence();     // GasTests.cpp: order of convergence of the gas solver on refined grids
void testDamBreakMartinMoyce(); // ParticleTests.cpp: dam-break front against Martin & Moyce 1952
void testVoronoiFracture();     // FractureTests.cpp: Voronoi cells of a convex body fill it, convex and watertight
void testTerrain();             // TerrainTests.cpp: 150 bodies on a 51 200-triangle static mesh (BVH)
void testParticlesManyBodies(); // ParticleTests.cpp: 30 000 particles find 150 bodies through the world tree
void testDestroyBody();         // RigidTests.cpp: one body removed, its slot reused, the rest untouched
void testCapsuleShape();        // RigidTests.cpp: capsule mass properties, lying / toppling on a box, raycasts
void testDeadLanding();         // RigidTests.cpp: restitution 0 - a dropped body stops at the touch
void testBarrelTouchesFloor(); // PileTests.cpp: a barrel touching a wide floor gets its normal and true depth
void testBarrelPileSettles();  // PileTests.cpp: 50 barrels dropped into a pile settle, stand still and sleep
void testFlatLanding();         // RigidTests.cpp: a cube dropped flat lands without turning or sliding
void testRemoveParticleGroup(); // ParticleTests.cpp: one soft body / liquid removed, the rest goes on

// relativity (RelativityTests.cpp): geodesics in Kerr, light bending, the shadow of a black hole
void testGeodesicInvariants();
void testLightDeflection();
void testPhotonSphere();
void testPerihelionPrecession();
void testShadow();
void testHorizonPenetration();
void testCurvatureFlatAndSphere();  // CurvatureTests.cpp: coordinates are not curvature; the 2-sphere; loops vs einstein()
void testCurvatureBlackHoles();     // CurvatureTests.cpp: Schwarzschild, Kerr, Reissner-Nordstrom curvature
void testCurvatureCosmology();      // CurvatureTests.cpp: de Sitter, Friedmann, the wormhole's exotic matter
void testTidesAndGammaGeodesic();   // CurvatureTests.cpp: tidal tensor, the geodesic equation with Gamma vs Hamilton
void testKerrHamiltonianGradient(); // SymplecticTests.cpp: Carter's exact gradient vs the metric's derivatives
void testSymplecticLongOrbits();    // SymplecticTests.cpp: energy bounded over thousands of orbits, RK4 drifts
void testSymplecticOrder();         // SymplecticTests.cpp: orders 2 (midpoint, Tao-2) and 4 (Tao-4, RK4)
void testSymplecticReversibility(); // SymplecticTests.cpp: forward then back returns to the start
// tensors and the linear-algebra toolbox (TensorTests.cpp)
void testTensorEinstein();
void testMatrixToolbox();

// scene graph of the editor (GraphTests.cpp): text file round trip, magnet forces, a graph scene runs
void testSceneGraphRoundTrip();
void testSceneGraphLocale();
// SoftBodyTests.cpp
void testSoftCantilever();
void testSoftBeamFrequency();
void testSoftHangingBar();
void testSoftSpinMomentum();
void testSoftDragThrough();
void testSoftVirtualWork();
void testSoftRest();
void testSoftTimeStepAndEnergy();
void testSoftSleep();
void testMagnetForce();
void testMagnetsAttract();
void testGraphSceneBuilds();
void testGraphNoInvisibleWalls(); // GraphTests.cpp: a column taller than the world box falls freely; a ball rolls off the floor
void testGraphSoftOverlapNoFlight(); // GraphTests.cpp: soft bodies made inside each other are pushed apart, not thrown
void testGraphCloth();
void testGraphEmitterFollows();
void testGraphFlammableCloth();
void testGraphMeshShape();
void testGraphGeometryOnly();
void testGraphColliderApart();  // GraphTests.cpp: a sphere on a box collider slides, a box on a sphere collider rolls
void testGraphColliderFit();    // GraphTests.cpp: a fitted capsule under a model, collider round trip, off-centre rebuild
void testGraphColliderOnly();   // GraphTests.cpp: a collider without the rigid role is a static obstacle
void testGraphLightsCameras();  // GraphTests.cpp: lights and cameras round trip, spot direction, camera frame, snapshot lights
// meta-objects (MetaObjectTests.cpp): one entity changes between frames, the rest does not notice
void testMetaWaterSoftCycle();
void testMetaRigidToClothAndBack();
void testMetaMagnetToggle();
void testMetaNoGrowth();
void testHierarchyInstancesRoundTrip();
void testExpandArrayPatterns();
void testArrayPileAndGrow();
void testInstancesShareRoles();
void testGluedGroupTumbles();

// the action principle by energy balances (ActionTests.cpp, docs/11-action.md): what has an action
// keeps its invariants, what dissipates loses exactly the work of its Rayleigh function
void testActionFriction();
void testActionDampedPendulum();
void testActionViscousBalance();
void testActionJouleBalance();
void testActionFreeRotation();

// the channels (ChannelTests.cpp, scene/Channels.h): what a chart shows, checked against physics
void testChannelsFreeFall();
void testChannelsCollision();
void testChannelsFriction();
void testChannelsStackForces();
void testChannelsTaylorGreen();
void testChannelsPoiseuilleProbe();
void testChannelsHydrostaticProbe();
void testChannelsCostAndInnocence();
void testChannelsCatalog();

// the book "Язык природы" (MathBookTests.cpp, docs/math/): every picture of the book is a real run
void testMathBookRounding();
void testMathBookVectors();
void testMathBookDerivative();
void testMathBookSchemes();
void testMathBookLinearSystems();
void testMathBookFlip();
void testMathBookDivergence();
void testMathBookSigma();

// verification registry (VerificationTests.cpp, only with RF_BUILD_VERIFY): the formulas on exact
// sequences, two fast cases through the runner, the board written into a scratch copy
void testVerificationSmoke();
