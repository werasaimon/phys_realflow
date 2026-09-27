// Wind tunnel scenes: a body in a uniform flow, the gas alone. The obstacle sits at the origin,
// 30 % down the length of the box; the tunnel measures its drag and lift (SurfaceLoads).
#include "samples/Samples.h"

namespace rf {

// The box every tunnel scene builds in: nothing but the gas around the obstacle.
static void buildTunnel(Simulation& sim) {
    sim.useGasBox({0.3f, 0.5f, 0.5f}); // body at 30% of the length
}

class TunnelSphereScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.obstacle.shape = ObstacleShape::Sphere;
        sim.obstacle.size = 0.5f;
    }
    void build(Simulation& sim) override { buildTunnel(sim); }
};

class TunnelCylinderScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.obstacle.shape = ObstacleShape::Cylinder;
        sim.obstacle.size = 0.4f;
        sim.obstacle.span = 1.4f;
        sim.vis.sliceField = GridField::Vorticity;
        sim.grid.params.vorticityConfinement = 0.5f;
    }
    void build(Simulation& sim) override { buildTunnel(sim); }
};

class TunnelWingScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.obstacle.shape = ObstacleShape::Wing;
        sim.obstacle.size = 0.9f;
        sim.obstacle.span = 1.2f;
        sim.obstacle.angleOfAttackDeg = 6.0f;
        sim.grid.params.usePlanformArea = true;
        sim.vis.sliceField = GridField::PressureCoeff;
    }
    void build(Simulation& sim) override { buildTunnel(sim); }
};

class TunnelStreamlinedScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.obstacle.shape = ObstacleShape::Streamlined;
        sim.obstacle.size = 1.0f;
        sim.obstacle.thickness = 0.3f;
    }
    void build(Simulation& sim) override { buildTunnel(sim); }
};

class TunnelCubeScene : public Scene {
public:
    void configure(Simulation& sim) override {
        sim.obstacle.shape = ObstacleShape::Cube;
        sim.obstacle.size = 0.4f;
    }
    void build(Simulation& sim) override { buildTunnel(sim); }
};

void addWindTunnelSamples(std::vector<SampleEntry>& out) {
    out.push_back({Preset::TunnelSphere, "Аэротруба", "Аэротруба: сфера", [] { return std::unique_ptr<Scene>(new TunnelSphereScene); }});
    out.push_back({Preset::TunnelCylinder, "Аэротруба", "Аэротруба: цилиндр", [] { return std::unique_ptr<Scene>(new TunnelCylinderScene); }});
    out.push_back({Preset::TunnelWing, "Аэротруба", "Аэротруба: крыло NACA", [] { return std::unique_ptr<Scene>(new TunnelWingScene); }});
    out.push_back({Preset::TunnelStreamlined, "Аэротруба", "Аэротруба: обтекаемое тело", [] { return std::unique_ptr<Scene>(new TunnelStreamlinedScene); }});
    out.push_back({Preset::TunnelCube, "Аэротруба", "Аэротруба: куб", [] { return std::unique_ptr<Scene>(new TunnelCubeScene); }});
}

} // namespace rf
