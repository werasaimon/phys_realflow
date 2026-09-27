// The magnetosphere scene: a plasma wind against a magnetised sphere, built on the MHD solver of
// the engine; and the registry of the plasma samples.
#include "samples/Samples.h"
#include "samples/plasma/TokamakScene.h"

#include <algorithm>

namespace rf {

class MagnetosphereScene : public Scene {
public:
    void configure(Simulation& sim) override {
        // A plasma wind blows at a magnetised sphere (Birkeland's terrella; the solar wind at the
        // Earth): the dipole's magnetic pressure B^2/2mu0 stops the flow where it matches the ram
        // pressure rho u^2 (the Chapman-Ferraro magnetopause, ~0.3 m upstream here); the plasma is
        // deflected around the magnetosphere and the field lines are swept back into a tail.
        // Conductivity 1e8 S/m: magnetic Reynolds number u L / eta ~ 100, the field is frozen in.
        sim.obstacle.shape = ObstacleShape::Sphere;
        sim.obstacle.size = 0.3f;
        sim.obstacle.position = Vector3(0.0f);
        sim.grid.params.domainSize = {2.4f, 1.4f, 1.4f}; // wide enough that the walls do not squeeze the flow
        sim.grid.params.resolutionX = 60;       // dx = 4 cm
        sim.grid.params.inflowSpeed = 1.5f;
        sim.grid.params.fluidDensity = 1.0f;
        sim.grid.params.pressureTolerance = 1e-3f;
        sim.grid.params.smokeDissipation = 0.7f; // the streaks fade in ~1.5 s: fresh wind only, no fog
        sim.grid.params.smokeRake = true;       // plasma streaks from the inflow show the deflection
        sim.grid.magnetic.enabled = true;
        sim.grid.magnetic.conductivity = 1e8f;
        sim.grid.magnetic.speedLimit = 4.0f;    // Boris correction near the poles (v_A up to 80 m/s there)
        sim.vis.sliceField = GridField::MagneticFlux;
        sim.vis.planetSurface = true;           // the terrella as a little Earth
        sim.vis.surfacePressure = false;
        sim.vis.showSlice = false;
        sim.vis.showStreamlines = false;
        sim.vis.gridDisplay = 0;
        // Velocity arrows in the equatorial plane (the dipole points along y): the wind slowing at
        // the magnetopause, turning around the magnetosphere and closing behind it.
        sim.vis.vectorDisplay = 1;
        sim.vis.sliceAxis = 1;
        sim.vis.slicePosition = 0.5f;
        sim.vis.vectorStride = 2;
        sim.vis.vectorScale = 1.3f;
    }

    void build(Simulation& sim) override {
        sim.useGasBox({0.3f, 0.5f, 0.5f}); // the magnet at 30 % of the length
        // The magnet: a dipole m = 600 A m^2 pointing up, from its vector potential
        // A = mu0/4pi (m x r) / r^3 (so div B = 0 exactly on the grid). Its equatorial field
        // mu0 m / 4pi r^3 balances the flow's ram pressure, B^2 / 2mu0 = rho u^2, at r = 0.29 m
        // (0.37 m counting the doubling of the field by the magnetopause currents).
        const Vector3 centre = sim.obstacle.position, moment(0.0f, 600.0f, 0.0f);
        sim.grid.magnetic.setBackgroundFromPotential([centre, moment](const Vector3& x) { // current-free B0
            const Vector3 r = x - centre;
            const float d = std::max(length(r), 0.08f); // inside the magnet (a conductor): no singularity
            return cross(moment, r) * (1e-7f / (d * d * d));
        });
    }
};

void addPlasmaSamples(std::vector<SampleEntry>& out) {
    out.push_back({Preset::Magnetosphere, "Плазма", "Плазма: магнит отклоняет поток (магнитосфера)",
                   [] { return std::unique_ptr<Scene>(new MagnetosphereScene); }});
    out.push_back({Preset::Tokamak, "Плазма", "Плазма: токамак (кольцо плазмы в тороидальном поле)",
                   [] { return std::unique_ptr<Scene>(new TokamakScene); }});
}

} // namespace rf
