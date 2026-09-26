cmake_minimum_required(VERSION 3.21)
project(PhysRealFlow VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(RF_BUILD_APP   "Build the Qt GUI application" ON)
option(RF_BUILD_TESTS "Build the solver tests"       ON)

# ---------------------------------------------------------------------------
# Physics SDK: standard C++17 only - no Qt, no other libraries. The GUI program below may use Qt,
# the SDK never does (it must build and run headless, and embed in any other engine).
# ---------------------------------------------------------------------------
add_library(rfcore STATIC
  src/math/Math.h src/math/Scalar.h src/math/AABB.h
  src/math/Vector2.h src/math/Vector3.h src/math/Vector4.h src/math/Quaternion.h
  src/math/Matrix3x3.h src/math/Matrix4x4.h
  src/math/MatrixNxN.h       src/math/MatrixNxN.cpp
  src/core/Parallel.h        src/core/Parallel.cpp
  src/core/Mesh.h            src/core/Mesh.cpp
  src/spatial/BVH.h          src/spatial/BVH.cpp
  src/spatial/AABBTree.h     src/spatial/AABBTree.cpp
  src/particles/ParticleSystem.h        src/particles/ParticleSystem.cpp
  src/particles/SoftBody.h         src/particles/SoftBody.cpp
  src/particles/Cloth.h            src/particles/Cloth.cpp
  src/grid/NSGridSolver.h    src/grid/NSGridSolver.cpp
  src/grid/SurfaceLoads.h    src/grid/SurfaceLoads.cpp
  src/grid/Combustion.h      src/grid/Combustion.cpp
  src/grid/Field3.h          src/grid/Field3.cpp
  src/grid/MagneticField.h   src/grid/MagneticField.cpp
  src/rigid/Shapes.h         src/rigid/Shapes.cpp
  src/rigid/GjkEpa.h         src/rigid/GjkEpa.cpp
  src/rigid/Narrowphase.h    src/rigid/Narrowphase.cpp
  src/rigid/BroadPhase.h     src/rigid/BroadPhase.cpp
  src/rigid/RigidWorld.h     src/rigid/RigidWorld.cpp
  src/rigid/XpbdSolver.cpp
  src/rigid/RigidBody.h
  src/rigid/Joints.h         src/rigid/Joints.cpp
  src/rigid/Ccd.h            src/rigid/Ccd.cpp
  src/rigid/ConvexDecomposition.h src/rigid/ConvexDecomposition.cpp
  src/sim/Simulation.h       src/sim/Simulation.cpp
)
target_include_directories(rfcore PUBLIC src)
find_package(Threads REQUIRED)
target_link_libraries(rfcore PUBLIC Threads::Threads)
# Guard: fail the configure step if anything but the thread library gets linked into the SDK.
get_target_property(RF_CORE_LIBS rfcore LINK_LIBRARIES)
if(NOT RF_CORE_LIBS STREQUAL "Threads::Threads")
  message(FATAL_ERROR "rfcore (physics SDK) must depend on the C++ standard library only, got: ${RF_CORE_LIBS}")
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(rfcore PRIVATE -Wall -Wextra -Wno-unused-parameter)
  target_compile_options(rfcore PUBLIC $<$<CONFIG:Release>:-O3 -march=native>)
endif()

# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------
if(RF_BUILD_APP)
  find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets OpenGL OpenGLWidgets Charts)
  qt_standard_project_setup()

  qt_add_executable(realflow WIN32
    app/main.cpp
    app/MainWindow.h      app/MainWindow.cpp
    app/Viewport.h        app/Viewport.cpp
    app/SimController.h   app/SimController.cpp
    app/ParamForm.h       app/ParamForm.cpp
    app/OrbitCamera.h     app/OrbitCamera.cpp
    app/ViewportTools.h   app/ViewportTools.cpp
    app/FluidSurfaceRenderer.h app/FluidSurfaceRenderer.cpp
    app/PlotPanel.h       app/PlotPanel.cpp
  )
  target_link_libraries(realflow PRIVATE rfcore
    Qt6::Core Qt6::Gui Qt6::Widgets Qt6::OpenGL Qt6::OpenGLWidgets Qt6::Charts)
  # Software OpenGL (Mesa llvmpipe, shipped with Qt) next to the program: with --software-gl the
  # same shaders run on the CPU - for PCs without a usable GPU / driver.
  if(WIN32)
    add_custom_command(TARGET realflow POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE_DIR:Qt6::Core>/opengl32sw.dll" "$<TARGET_FILE_DIR:realflow>")
  endif()
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(RF_BUILD_TESTS)
  enable_testing()
  add_executable(rf_tests tests/tests.cpp)
  target_link_libraries(rf_tests PRIVATE rfcore)
  add_test(NAME rf_tests COMMAND rf_tests)
endif()
