#include "gas/SurfaceLoads.h"

#include <cstdio>
#include <filesystem>
#include <fstream>

namespace rf {

SurfaceLoads computeSurfaceLoads(const GasSolver& gas, const TriMesh& surface, const Vector3& momentRef, float refArea,
                                 float refLength) {
    SurfaceLoads out;
    const GasParams& prm = gas.params;
    const float q = gas.dynamicPressure();
    const float rho = prm.fluidDensity, nu = std::max(prm.kinematicViscosity, 1e-12f);
    // The gas is read at one and two cell spacings outside the surface (the voxelised body may
    // stick out of the true surface by up to half a cell; these points are in gas cells), and the
    // pressure is extrapolated linearly to the wall: p_wall = 2 p(1 dx) - p(2 dx). Reading at one
    // cell only would miss the pressure rise towards a stagnation point (Cp 0.7 instead of 1).
    const float probe = gas.dx();

    out.triangles.reserve(surface.triangles.size());
    for (size_t t = 0; t < surface.triangles.size(); ++t) {
        const auto& tri = surface.triangles[t];
        const Vector3 &a = surface.positions[tri[0]], &b = surface.positions[tri[1]], &c = surface.positions[tri[2]];
        Vector3 cr = cross(b - a, c - a);
        float area = 0.5f * length(cr);
        if (area <= 0) continue;
        TriangleLoad L;
        L.centroid = (a + b + c) / 3.0f;
        L.normal = cr / (2.0f * area);
        L.area = area;

        const Vector3 x = L.centroid + L.normal * probe;
        const float p = 2.0f * gas.fluidPressureAt(x) - gas.fluidPressureAt(L.centroid + L.normal * (2.0f * probe));
        L.cp = p / q;

        // Wall shear from the tangential gas velocity (the body is at rest in the tunnel).
        Vector3 u = gas.fluidVelocityAt(x);
        Vector3 ut = u - L.normal * dot(u, L.normal);
        float mag = length(ut);
        Vector3 tau(0.0f);
        if (prm.wallFriction && mag > 1e-6f) {
            // Same wall function and running length as the solver's own skin friction.
            const float cfLocal = skinFrictionCoefficient(mag * gas.wallLength() / nu);
            tau = ut * (0.5f * rho * cfLocal * mag);
        }
        L.cf = length(tau) / q;

        const Vector3 Fp = L.normal * (-p * area), Ff = tau * area;
        L.force = Fp + Ff;
        out.pressureForce += Fp;
        out.frictionForce += Ff;
        out.moment += cross(L.centroid - momentRef, L.force);
        out.wettedArea += area;
        out.triangles.push_back(L);
    }
    const float qa = q * std::max(refArea, 1e-12f);
    const Vector3 F = out.pressureForce + out.frictionForce;
    out.cdPressure = out.pressureForce.x / qa;
    out.cdFriction = out.frictionForce.x / qa;
    out.cd = F.x / qa;
    out.cl = F.y / qa;
    out.cs = F.z / qa;
    out.cm = out.moment.z / (qa * std::max(refLength, 1e-12f));
    return out;
}

bool saveSurfaceLoadsCsv(const std::string& path, const SurfaceLoads& loads) {
    // The path is UTF-8 (as loadMesh's): u8path keeps Cyrillic folder names working on Windows.
    std::ofstream f(std::filesystem::u8path(path));
    if (!f) return false;
    char line[512];
    std::snprintf(line, sizeof(line), "# Cd=%.5f (pressure %.5f, friction %.5f) Cl=%.5f Cs=%.5f Cm=%.5f wetted_area=%.6f\n",
                  loads.cd, loads.cdPressure, loads.cdFriction, loads.cl, loads.cs, loads.cm, loads.wettedArea);
    f << line << "triangle,cx,cy,cz,nx,ny,nz,area,Cp,Cf,Fx,Fy,Fz\n";
    for (size_t i = 0; i < loads.triangles.size(); ++i) {
        const TriangleLoad& t = loads.triangles[i];
        std::snprintf(line, sizeof(line), "%zu,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", i, t.centroid.x,
                      t.centroid.y, t.centroid.z, t.normal.x, t.normal.y, t.normal.z, t.area, t.cp, t.cf, t.force.x, t.force.y,
                      t.force.z);
        f << line;
    }
    f.close();
    return bool(f);
}

} // namespace rf
