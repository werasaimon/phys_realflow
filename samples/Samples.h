#pragma once
// The ready-made scenes of PhysRealFlow, as Box2D keeps its demos in samples/: each one is a
// Scene (src/scene/Scene.h) built on the Simulation facade, and this registry lists them for the
// viewer's menu, the command line (--preset N) and the tests. The engine (rfcore) knows none of
// them; a program that wants its own scenes writes them the same way and needs no registry.
#include "scene/Scene.h"
#include "scene/Simulation.h"

#include <memory>
#include <vector>

namespace rf {

// The scenes by name; the value is the scene's number in the menu and on the command line
// (--preset N), kept stable so that the numbers in the documentation stay right.
enum class Preset {
    DamBreak, FluidObstacle, FloatingBodies, JetOnObject,
    TunnelSphere, TunnelCylinder, TunnelWing, TunnelStreamlined, TunnelCube, SmokePlume, SmokeSphere,
    RigidFalling, RigidGranular, RigidPyramid, RigidConvex, RigidTower, RigidJoints, RigidCcd, RigidTeapots,
    SmokeBodies, SoftCloth, GasSoftCloth, Hydro, Fire, Water, Magnetosphere, Tokamak, Terrain,
    Dzhanibekov, NewtonCradle, GalileoTower, // lessons: one law each, its number next to the theory
    RigidChains, RigidRace, RigidBowl, RigidTables, RigidCups, TerrainModels,
    RigidTower200, RigidSandbox, TorusChains,
    HangingTorus,
    Count
};

struct SampleEntry {
    Preset id;
    const char* category; // the viewer groups the menu by it
    const char* name;
    std::unique_ptr<Scene> (*create)();
};

// All samples, index = Preset. Every value of Preset has exactly one entry (SceneTests checks).
const std::vector<SampleEntry>& samples();

// Creates the sample and loads it into the simulation (Simulation::load).
void loadSample(Simulation& sim, Preset p);
inline void loadSample(Simulation& sim, int index) { loadSample(sim, Preset(index)); }

// Each file of samples/ adds its scenes to the list (any order; samples() sorts them by Preset).
void addLiquidSamples(std::vector<SampleEntry>& out);      // LiquidScenes.cpp
void addWindTunnelSamples(std::vector<SampleEntry>& out);  // WindTunnelScenes.cpp
void addSmokeSamples(std::vector<SampleEntry>& out);       // SmokeScenes.cpp
void addRigidSamples(std::vector<SampleEntry>& out);       // RigidScenes.cpp
void addSoftBodySamples(std::vector<SampleEntry>& out);    // SoftBodyScenes.cpp
void addHydroSamples(std::vector<SampleEntry>& out);       // HydroScene.cpp
void addFireSamples(std::vector<SampleEntry>& out);        // FireScene.cpp
void addPlasmaSamples(std::vector<SampleEntry>& out);      // plasma/MagnetosphereScene.cpp, plasma/TokamakScene.cpp
void addLessonSamples(std::vector<SampleEntry>& out);      // LessonScenes.cpp: Dzhanibekov, Newton's cradle, Galileo
void addNonConvexSamples(std::vector<SampleEntry>& out);   // NonConvexScenes.cpp: race, bowl, tables, cups, models on terrain
void addTorusChainSamples(std::vector<SampleEntry>& out); // TorusChainsScene.cpp: long contact-only chains
void addHangingTorusSample(std::vector<SampleEntry>& out); // HangingTorusScene.cpp: short interactive chain

} // namespace rf
