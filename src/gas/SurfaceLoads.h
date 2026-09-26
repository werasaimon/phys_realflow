#pragma once
// Aerodynamic loads per triangle of a body surface, read from the gas grid.
//
// For every triangle of the body mesh: the gas pressure and velocity are read just outside the
// surface (interpolated over gas cells only), giving
//   Cp = p / q                         pressure coefficient (q = dynamic pressure of the free stream)
//   tau = 1/2 rho Cf |u_t| u_t         wall shear from the tangential velocity (wall function,
//                                      same flat-plate correlations as the solver)
//   F = (-p n + tau) A                 force on the triangle
// The sums give the pressure and friction forces, the moment about a reference point and the
// usual coefficients. The distributions (Cp, Cf per triangle) are what aerodynamic studies need.

#include "core/Mesh.h"
#include "gas/GasSolver.h"

#include <string>
#include <vector>

namespace rf {

struct TriangleLoad {
    Vector3 centroid;
    Vector3 normal; // outward, unit
    float area = 0;
    float cp = 0;   // pressure coefficient
    float cf = 0;   // skin-friction coefficient |tau| / q
    Vector3 force;  // pressure + friction [N]
};

struct SurfaceLoads {
    std::vector<TriangleLoad> triangles;
    Vector3 pressureForce, frictionForce; // [N]
    Vector3 moment;                       // [N m] about the reference point
    float cdPressure = 0, cdFriction = 0; // drag split: pressure (form) and friction
    float cd = 0, cl = 0, cs = 0, cm = 0; // drag (x), lift (y), side (z), pitching moment (about z)
    float wettedArea = 0;
};

// refArea / refLength: reference area and length of the coefficients (e.g. frontal area and
// diameter, or planform area and chord); momentRef: point about which the moment is taken.
SurfaceLoads computeSurfaceLoads(const GasSolver& gas, const TriMesh& surface, const Vector3& momentRef, float refArea,
                                 float refLength);

// One row per triangle: centroid, normal, area, Cp, Cf, force.
bool saveSurfaceLoadsCsv(const std::string& path, const SurfaceLoads& loads);

} // namespace rf
