cmake_minimum_required(VERSION 3.21)
project(PhysRealFlow VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(RF_BUILD_TESTS "Build the solver tests"       ON)
option(RF_WERROR      "Warnings are errors (CI)"     OFF)
option(RF_SANITIZE    "AddressSanitizer + UndefinedBehaviorSanitizer (GCC/Clang, CI)" OFF)
option(RF_STRICT_FP   "Bit-reproducible floating point across builds: no FMA contraction, no -march=native" OFF)
option(RF_LTO         "Link-time optimisation of the SDK (inlining across its files)" OFF)
option(RF_BUILD_VERIFY "Build the verification registry and rf_verify (ASME V&V, docs/10)" ON)

# ---------------------------------------------------------------------------
# Physics SDK: standard C++17 only - no Qt, no other libraries. The GUI program below may use Qt,
# the SDK never does (it must build and run headless, and embed in any other engine).
# ---------------------------------------------------------------------------
add_library(rfcore STATIC
  src/math/Math.h src/math/Scalar.h src/math/AABB.h src/math/ElementaryFunctions.h
  src/math/Vector2.h src/math/Vector3.h src/math/Vector4.h src/math/Quaternion.h
  src/math/Matrix3x3.h src/math/Matrix4x4.h
  src/math/MatrixNxN.h       src/math/MatrixNxN.cpp
  src/math/MatrixFunctions.cpp
  src/math/Tensor.h          src/math/Tensor.cpp
  src/core/Parallel.h        src/core/Parallel.cpp
  src/core/Mesh.h            src/core/Mesh.cpp
  src/core/Format.h
  src/core/Probe.h           src/core/Probe.cpp
  src/spatial/BVH.h          src/spatial/BVH.cpp
  src/spatial/AABBTree.h     src/spatial/AABBTree.cpp
  # rigid bodies
  src/rigid/RigidBody.h
  src/rigid/Shapes.h              src/rigid/Shapes.cpp
  src/rigid/GjkEpa.h              src/rigid/GjkEpa.cpp
  src/rigid/NarrowPhase.h         src/rigid/NarrowPhase.cpp
  src/rigid/BroadPhase.h          src/rigid/BroadPhase.cpp
  src/rigid/TimeOfImpact.h        src/rigid/TimeOfImpact.cpp
  src/rigid/Fracture.h            src/rigid/Fracture.cpp
  src/rigid/Joints.h              src/rigid/Joints.cpp
  src/rigid/ConvexDecomposition.h src/rigid/ConvexDecomposition.cpp
  src/rigid/RigidWorld.h          src/rigid/RigidWorld.cpp
  src/rigid/ContactSolver.cpp
  src/rigid/Islands.cpp
  src/rigid/ShockPropagation.cpp
  src/rigid/FreeRotation.cpp
  src/rigid/ContactLcp.h
  src/rigid/Grab.cpp
  src/rigid/XpbdSolver.cpp
  src/rigid/RigidDebugDraw.cpp
  src/rigid/WatchedPairDraw.cpp
  # particles: liquid, soft bodies, cloth
  src/particles/ParticleSystem.h  src/particles/ParticleSystem.cpp
  src/particles/ParticleDebugDraw.cpp
  src/particles/DensitySolver.cpp
  src/particles/ParticleContacts.cpp
  src/particles/SoftBody.h        src/particles/SoftBody.cpp
  src/particles/SoftTets.cpp      src/particles/ParticleSoftBodies.cpp
  src/particles/Cloth.h           src/particles/Cloth.cpp
  # gas on the MAC grid, fire
  src/gas/Field3.h                src/gas/Field3.cpp
  src/gas/GasSolver.h             src/gas/GasSolver.cpp
  src/gas/GasDebugDraw.cpp
  src/gas/GasSolverInternal.h
  src/gas/Advection.cpp
  src/gas/AdvectionReflection.cpp
  src/gas/Viscosity.cpp
  src/gas/PressureSolver.cpp
  src/gas/Multigrid.h             src/gas/Multigrid.cpp
  src/gas/MovingSolids.cpp
  src/gas/Heat.cpp
  src/gas/ExternalFields.cpp
  src/gas/Combustion.h            src/gas/Combustion.cpp
  src/gas/SurfaceLoads.h          src/gas/SurfaceLoads.cpp
  # plasma
  src/plasma/MagneticField.h      src/plasma/MagneticField.cpp
  # relativity: Kerr / Schwarzschild spacetime, geodesics, the reference ray tracer
  src/relativity/Metric.h         src/relativity/Metric.cpp
  src/relativity/Geodesic.h       src/relativity/Geodesic.cpp
  src/relativity/RayTracer.h      src/relativity/RayTracer.cpp
  src/relativity/Spacetime.h      src/relativity/Spacetime.cpp
  src/relativity/Curvature.h      src/relativity/Curvature.cpp
  src/relativity/KerrHamiltonian.h src/relativity/KerrHamiltonian.cpp
  src/relativity/SymplecticGeodesic.h src/relativity/SymplecticGeodesic.cpp
  # the scene facade: the solvers together; the scenes themselves are in samples/
  src/scene/Scene.h
  src/scene/Simulation.h          src/scene/Simulation.cpp
  src/scene/Obstacle.cpp
  src/scene/Coupling.cpp
  src/scene/Snapshot.cpp
  src/scene/SceneGraph.h          src/scene/SceneGraph.cpp
  src/scene/GraphScene.cpp
  src/scene/MetaObjects.cpp
  src/scene/EntityShapes.cpp
  src/scene/SceneHierarchy.cpp
  src/scene/GluedGroups.cpp
  src/scene/Magnets.h             src/scene/Magnets.cpp
  src/scene/Channels.h            src/scene/Channels.cpp
  src/scene/ChannelMath.h         src/scene/ObjectChannels.cpp
  src/scene/FieldProbes.cpp
)
target_include_directories(rfcore PUBLIC src)

# ---------------------------------------------------------------------------
# Sample scenes (as Box2D's samples/): solutions built on the engine, used by the app and the tests
# ---------------------------------------------------------------------------
add_library(rfsamples STATIC
  samples/Samples.h               samples/Samples.cpp
  samples/SamplesInternal.h
  samples/Models.h                samples/Models.cpp
  samples/LiquidScenes.cpp
  samples/SoftBodyScenes.cpp
  samples/RigidScenes.cpp
  samples/WindTunnelScenes.cpp
  samples/SmokeScenes.cpp
  samples/FireScene.cpp
  samples/HydroScene.cpp
  samples/plasma/Tokamak.h        samples/plasma/Tokamak.cpp
  samples/plasma/TokamakScene.h   samples/plasma/TokamakScene.cpp
  samples/plasma/MagnetosphereScene.cpp
  samples/LessonScenes.cpp
)
target_include_directories(rfsamples PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(rfsamples PUBLIC rfcore)
if(RF_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT RF_IPO_OK OUTPUT RF_IPO_MSG)
  if(RF_IPO_OK)
    set_property(TARGET rfcore PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
  else()
    message(WARNING "LTO not supported here: ${RF_IPO_MSG}")
  endif()
endif()
find_package(Threads REQUIRED)
target_link_libraries(rfcore PUBLIC Threads::Threads)
# Guard: fail the configure step if anything but the thread library gets linked into the SDK.
get_target_property(RF_CORE_LIBS rfcore LINK_LIBRARIES)
if(NOT RF_CORE_LIBS STREQUAL "Threads::Threads")
  message(FATAL_ERROR "rfcore (physics SDK) must depend on the C++ standard library only, got: ${RF_CORE_LIBS}")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(rfcore PRIVATE -Wall -Wextra -Wno-unused-parameter)
  if(RF_WERROR)
    target_compile_options(rfcore PRIVATE -Werror)
  endif()
  if(RF_STRICT_FP)
    # The same numbers whatever the compiler inlines: no fused multiply-add, no CPU-specific code.
    # RF_STRICT_FP tells the tests: the golden hash of the rigid benchmark scene must match.
    target_compile_options(rfcore PUBLIC $<$<CONFIG:Release>:-O3> -ffp-contract=off)
    target_compile_definitions(rfcore PUBLIC RF_STRICT_FP=1)
  else()
    target_compile_options(rfcore PUBLIC $<$<CONFIG:Release>:-O3 -march=native>)
  endif()
  if(RF_SANITIZE)
    target_compile_options(rfcore PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(rfcore PUBLIC -fsanitize=address,undefined)
  endif()
elseif(MSVC)
  target_compile_options(rfcore PRIVATE /W4 /wd4100 /utf-8)
  if(RF_WERROR)
    target_compile_options(rfcore PRIVATE /WX)
  endif()
endif()

# The Qt editor lives in its own repository (https://github.com/werasaimon/phys_realflow_editor)
# with this SDK as a git submodule; nothing here depends on Qt.

# ---------------------------------------------------------------------------
# Verification and validation (verification/, docs/10-verification.md): the benchmark registry
# as a library, rf_verify runs it and writes the results and the board (docs/09).
# ---------------------------------------------------------------------------
if(RF_BUILD_VERIFY)
  add_library(rfverify STATIC
    verification/Benchmark.h        verification/Benchmark.cpp
    verification/Convergence.h      verification/Convergence.cpp
    verification/ValidationMetric.h verification/ValidationMetric.cpp
    verification/Uncertainty.h      verification/Uncertainty.cpp
    verification/Mms.h              verification/Mms.cpp
    verification/Significance.h     verification/Significance.cpp
    verification/Sha256.h           verification/Sha256.cpp
    verification/Blinding.h         verification/Blinding.cpp
    verification/Report.h           verification/Report.cpp
    verification/Board.cpp
    verification/Charts.cpp
    verification/Runner.h           verification/Runner.cpp
    verification/Procedure.h        verification/Procedure.cpp
    verification/cases/RigidCases.cpp
    verification/cases/GasCases.cpp
    verification/cases/ParticleCases.cpp
    verification/cases/RelativityCases.cpp
    verification/cases/PlasmaCases.cpp
  )
  target_include_directories(rfverify PUBLIC verification)
  target_link_libraries(rfverify PUBLIC rfcore)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # The cases are written as brace lists that leave the trailing fields (split, commitment,
    # stochastic) at their defaults on purpose: no warning for that.
    target_compile_options(rfverify PRIVATE -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers $<$<BOOL:${RF_WERROR}>:-Werror>)
  elseif(MSVC)
    target_compile_options(rfverify PRIVATE /W4 /wd4100 /utf-8 $<$<BOOL:${RF_WERROR}>:/WX>)
  endif()
  add_executable(rf_verify verification/rf_verify.cpp)
  target_link_libraries(rf_verify PRIVATE rfverify)
  target_compile_definitions(rf_verify PRIVATE RF_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
  if(MSVC)
    target_compile_options(rf_verify PRIVATE /utf-8)
  endif()
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(RF_BUILD_TESTS)
  enable_testing()
  add_executable(rf_tests
    tests/main.cpp
    tests/TestRunner.h tests/Tests.h
    tests/CoreTests.cpp
    tests/DebugDrawTests.cpp
    tests/MathTests.cpp
    tests/RigidTests.cpp tests/HardContactTests.h
    tests/PileTests.cpp
    tests/ParticleTests.cpp
    tests/SoftFemTests.cpp
    tests/GasTests.cpp
    tests/PlasmaTests.cpp
    tests/SceneTests.cpp
    tests/DeterminismTests.cpp
    tests/TerrainTests.cpp
    tests/FractureTests.cpp
    tests/RelativityTests.cpp
    tests/GraphTests.cpp
    tests/MetaObjectTests.cpp
    tests/ArrayTests.cpp
    tests/ActionTests.cpp
    tests/ChannelTests.cpp
    tests/MathBookTests.cpp
    tests/TensorTests.cpp
    tests/CurvatureTests.cpp
    tests/SymplecticTests.cpp
    tests/PerfTests.cpp)
  target_link_libraries(rf_tests PRIVATE rfcore rfsamples)
  if(RF_BUILD_VERIFY) # the smoke test of the registry
    target_sources(rf_tests PRIVATE tests/VerificationTests.cpp)
    target_link_libraries(rf_tests PRIVATE rfverify)
    target_compile_definitions(rf_tests PRIVATE RF_WITH_VERIFY)
    if(MSVC)
      target_compile_options(rf_tests PRIVATE /utf-8)
    endif()
  endif()
  add_test(NAME rf_tests COMMAND rf_tests)
endif()
