#pragma once
// Helpers shared between the files of samples/ (not part of the public sample registry).
#include "scene/Simulation.h"

namespace rf {

// A closed box of voxels with a hot smoky sphere near the bottom; the gas rises by buoyancy
// (the SmokeSphere, SmokeBodies and GasSoftCloth scenes). bodiesInside: the scene puts bodies
// into the box - then dense smoke that leaves through an open ceiling, and no voxel grid or
// velocity arrows drawn (smoke only, 3/4 view). Defined in SmokeScenes.cpp.
void configureClosedSmokeBox(Simulation& sim, bool bodiesInside);

// The terrain of the Terrain and TerrainModels scenes: a 12 x 12 m height field of 51 200
// triangles as the static obstacle, its lowest point at y = 0. Defined in RigidScenes.cpp.
void configureTerrain(Simulation& sim);

} // namespace rf
