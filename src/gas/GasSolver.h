#pragma once
// Incompressible Navier-Stokes on a staggered (MAC) grid - "virtual wind tunnel" and gas/smoke.
//  * semi-Lagrangian RK2 advection, optional MacCormack (2nd order) with limiter
//  * implicit viscosity solved to a tolerance (conjugate gradient), no-slip on the wall itself
//    (ghost values), incremental pressure correction (Viscosity.cpp); vorticity confinement,
//    Boussinesq buoyancy
//  * pressure projection: preconditioned conjugate gradient, multigrid V-cycle preconditioner (MGPCG)
//  * obstacles voxelised from a triangle mesh via MeshBVH inside/outside queries
//  * pressure force integration on the body -> drag / lift coefficients
//  * moving solids (two-way coupling with rigid bodies): cells inside a body are solid and the
//    faces around them carry the body's velocity v + w x r (no-through condition); the pressure
//    solve then pushes the gas out of the way, and the pressure integrated over the body's faces
//    gives the force and torque of the gas on the body (Bridson, "Fluid Simulation for Computer
//    Graphics", ch. 5; Crane, Llamas, Tariq, GPU Gems 3 ch. 30)
//  * skin friction at solid surfaces by a wall function: the boundary layer (~mm) is far below the
//    grid spacing, so the wall shear comes from the flat-plate correlations (Schlichting: laminar
//    Cf = 1.328/sqrt(Re), turbulent Cf = 0.074 Re^-0.2) applied to the tangential velocity of the
//    gas next to the wall; the momentum taken from the gas is given to the body (action = reaction)
//  * fire: a fuel field that burns in hot gas, releasing heat, soot and expanding gas (Combustion)

#include "gas/Combustion.h"
#include "gas/Field3.h"
#include "gas/Multigrid.h"
#include "plasma/MagneticField.h"
#include "spatial/BVH.h"
#include "math/Math.h"

#include <cstdint>
#include <functional>
#include <vector>

namespace rf {

enum class BoundaryType { Wall, Inflow, Outflow };
enum class GridField { Speed, Pressure, PressureCoeff, Vorticity, Smoke, Temperature, VelocityX, VelocityY,
                       MagneticFlux, CurrentDensity };

struct GasParams {
    Vector3 domainSize{4.0f, 2.0f, 2.0f};  // [m]                          (*reset)
    int resolutionX = 96;                // cells along X                 (*reset)
    float inflowSpeed = 10.0f;           // [m/s]
    float fluidDensity = 1.225f;         // [kg/m^3] (air)
    float kinematicViscosity = 1.5e-5f;  // [m^2/s]
    float cfl = 2.0f;
    int maxPressureIterations = 400;
    float pressureTolerance = 1e-4f;
    // The implicit viscous step (Viscosity.cpp): conjugate gradient to this relative residual, at
    // most this many iterations per velocity component.
    float viscosityTolerance = 1e-6f;
    int maxViscosityIterations = 1000;
    PressurePreconditioner pressurePreconditioner = PressurePreconditioner::Multigrid; // of the PCG (Multigrid.h)
    bool maccormack = true;
    // Advection-reflection (AdvectionReflection.cpp): the projection that removes the divergent part
    // of the flow is done at mid-step as a REFLECTION, u <- 2 P(u) - u, which keeps the kinetic
    // energy, and once more at the end - instead of one projection at the end, which throws the
    // divergent part's energy away every step (Zehnder, Narain & Thomaszewski, SIGGRAPH 2018; second
    // order in time: Narain, Zehnder & Thomaszewski 2019). Numerical dissipation 13.8 % -> 1.6 % of
    // the loss, second order in time, at 1.5 ... 2.1 times the cost of a frame. OFF for now (user,
    // 2026-09-28), temporarily: with it the pressure's numerical boundary layer at a body grows, and
    // the surface-loads test's two pressure forces on the sphere drift 35.5 % apart (limit 30 %; the
    // classic step 21.5 %). It is switched on once that layer is fixed (rotational pressure correction).
    bool advectionReflection = false;
    float vorticityConfinement = 0.0f;
    bool smokeRake = true;               // smoke streaks injected at the inflow
    float smokeBuoyancy = 0.0f;          // [m/s^2] per unit smoke density (sinks)
    float heatBuoyancy = 0.0f;           // [m/s^2] per unit temperature (rises)
    float smokeDissipation = 0.0f;       // [1/s]
    float temperatureDissipation = 0.0f; // [1/s] cooling of the gas
    bool usePlanformArea = false;        // reference area for Cd/Cl: frontal (default) or planform
    bool wallFriction = true;            // skin friction on solid surfaces (wall function)
    BoundaryType bc[6] = {BoundaryType::Inflow, BoundaryType::Outflow, BoundaryType::Wall,
                          BoundaryType::Wall,   BoundaryType::Wall,    BoundaryType::Wall}; // -x +x -y +y -z +z (*reset)
};

struct HeatSource {
    bool enabled = false;
    Vector3 center;
    float radius = 0.1f;
    float temperature = 1.0f;
    float smoke = 1.0f;
    float fuel = 0.0f;      // fuel concentration kept in the source (a gas burner; see Combustion)
    Vector3 velocity{0.0f}; // [m/s] the gas leaves the source with this velocity (0 = no jet)
};

// Gas given off by something in the flow at a point (burning cloth: its pyrolysis gas and heat).
struct GasEmission {
    Vector3 position;
    float fuel = 0;  // [m^3] of fuel gas at unit concentration
    float heat = 0;  // [J]
    float smoke = 0; // [m^3] of smoke at unit density
};

// A local, one-shot perturbation of the flow (e.g. from the mouse): inside a sphere the velocity is
// blended towards `velocity` and smoke/heat are added, with a smooth (1 - r^2/R^2)^2 falloff.
// The next pressure projection makes the result divergence-free again.
struct Disturbance {
    Vector3 center;
    float radius = 0.1f;
    Vector3 velocity;
    float velocityBlend = 1.0f; // 0 = leave the velocity untouched (smoke/heat only)
    float smoke = 0.0f;
    float heat = 0.0f;
};

// Skin-friction coefficient of a flat plate at the Reynolds number Re of its running length
// (Schlichting): laminar 1.328 / sqrt(Re), turbulent 0.074 Re^-0.2 beyond Re = 5e5. Used by the
// solver's wall function and by SurfaceLoads, so both report the same friction.
inline float skinFrictionCoefficient(float Re) {
    Re = std::max(Re, 1.0f);
    return Re < 5e5f ? 1.328f / std::sqrt(Re) : 0.074f * std::pow(Re, -0.2f);
}

// A rigid body as a moving obstacle for the gas (world space).
struct MovingSolid {
    AABB bounds;
    Vector3 position, velocity, angularVelocity; // velocity of a point x: v + w x (x - position)
    std::function<bool(const Vector3&)> inside;  // x inside the body
    float length = 0.2f;                      // characteristic length for the wall Reynolds number [m]
    bool resting = false;                     // did not move since the last step: its cells are reused
    Vector3 pointVelocity(const Vector3& x) const { return velocity + cross(angularVelocity, x - position); }
};

class GasSolver {
public:
    GasParams params;
    HeatSource source;
    Combustion combustion; // fire (off by default); the temperature field is then in K above ambient
    // Electromagnetic field of a conducting gas (plasma): Lorentz force on the flow, induction by
    // the flow, Joule heating (off by default). The static obstacle cells are perfect conductors.
    MagneticField magnetic;

    // A closed vessel (a tokamak's torus): when set, every cell outside it is wall - a perfect
    // conductor for the magnetic field. The gas lives inside. Set before reset().
    std::function<bool(const Vector3&)> vessel;

    // Verification hooks (ExternalFields.cpp): the velocity (at the face centres) and the
    // temperature set from a function of the world point, and a heat source term [K/s] of the
    // point and the time, added every step after the conduction (fire on). For exact flows and
    // manufactured solutions (verification/); no scene uses them.
    void setVelocity(const std::function<Vector3(const Vector3&)>& velocity);
    void setTemperature(const std::function<float(const Vector3&)>& temperature);
    std::function<float(const Vector3& x, double t)> heatSource;

    // origin = lower corner of the domain in world space.
    void reset(const Vector3& origin, const MeshBVH* obstacle);
    // The smoke (the tracer: plasma, dye) set from a function of the world point, gas cells only.
    void setTracer(const std::function<float(const Vector3&)>& density);
    // Advances by an adaptive CFL time step <= maxDt; returns the step used.
    float step(float maxDt);
    void applyDisturbance(const Disturbance& d);
    // Moving obstacles for the following steps (replaces the previous set). After each step the
    // pressure force [N] and torque [N m, about MovingSolid::position] of the gas on every solid
    // are available in the same order.
    void setMovingSolids(std::vector<MovingSolid> solids) { moving_ = std::move(solids); }
    // Momentum put into the gas from outside [N s] at a world point (e.g. the reaction of the drag
    // on cloth particles); spread over the nearby faces, applied in the next step.
    void addImpulse(const Vector3& x, const Vector3& impulse) { pendingImpulses_.push_back({x, impulse}); }
    // Fuel, heat and smoke given to the gas at a point; added to its cell in the next step.
    void addEmission(const GasEmission& e) { pendingEmissions_.push_back(e); }
    // A liquid (particles of the particle system) as a moving obstacle for the gas: every cell that
    // holds a particle is solid and moves with its particles' mean velocity. Set before a step.
    void setLiquid(std::vector<Vector3> positions, std::vector<Vector3> velocities) {
        liquidPos_ = std::move(positions);
        liquidVelIn_ = std::move(velocities);
    }
    int liquidCellCount() const { return liquidCells_; }
    // Is there gas in the cell of x or next to it? (Liquid particles there feel the air.)
    bool touchesGas(const Vector3& x) const;
    const std::vector<Vector3>& movingForces() const { return movForce_; }
    const std::vector<Vector3>& movingTorques() const { return movTorque_; }
    int movingSolidCells() const { return movingCells_; }
    // Time of the last step spent on solids: voxelising moving bodies + their forces + wall friction [ms].
    float lastSolidCouplingMs() const { return solidMs_; }
    float lastPressureMs() const { return pressureMs_; }

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }
    float dx() const { return dx_; }
    Vector3 origin() const { return origin_; }
    AABB domain() const { return AABB(origin_, origin_ + Vector3(float(nx_), float(ny_), float(nz_)) * dx_); }

    bool solid(int i, int j, int k) const { return solid_[cidx(i, j, k)] != 0; }
    bool staticSolid(int i, int j, int k) const { return (solid_[cidx(i, j, k)] & kStatic) != 0; }
    Vector3 velocityAt(const Vector3& world) const;
    float pressureAt(const Vector3& world) const { return p_.sample((world - origin_) / dx_); }
    float temperatureAt(const Vector3& world) const { return temp_.sample((world - origin_) / dx_); }
    // Fire: thermal radiation of the flame arriving at a point [W/m^2]. The hot gas is an optically
    // thin emitter: every cell hotter than ~700 K radiates 4 kappa sigma (T^4 - T_air^4) V evenly in
    // all directions (Siegel & Howell, "Thermal Radiation Heat Transfer"), without shadowing.
    float irradianceAt(const Vector3& world) const;
    // Pressure / cell velocity interpolated over gas cells only (solid cells get no weight): the
    // gas state right at a solid surface, where plain interpolation would mix in the body.
    float fluidPressureAt(const Vector3& world) const;
    Vector3 fluidVelocityAt(const Vector3& world) const;
    float cellValue(GridField f, int i, int j, int k) const;
    const Field3& smoke() const { return smoke_; }
    const Field3& temperature() const { return temp_; }
    const Field3& fuel() const { return fuel_; }
    // Fire: heat released by the combustion in the last step [W] (the flame's power).
    float heatReleaseRate() const { return heatReleaseRate_; }

    // Diagnostics
    Vector3 bodyForce() const { return force_; }             // [N], pressure + friction force on obstacle
    Vector3 frictionForce() const { return staticFriction_; } // [N], skin-friction part of bodyForce()
    float referenceArea() const { return refArea_; }        // [m^2]
    float frontalArea() const { return frontalArea_; }
    float planformArea() const { return planformArea_; }
    float dragCoefficient() const { return cd_; }
    float liftCoefficient() const { return cl_; }
    float sideCoefficient() const { return cs_; }
    float dragCoefficientAvg() const { return cdAvg_; }
    float liftCoefficientAvg() const { return clAvg_; }
    int lastPressureIterations() const { return lastIters_; }
    float lastResidual() const { return lastResidual_; }
    float maxDivergence() const { return maxDivergence_; } // [1/s] after projection, should be ~0
    float totalSmoke() const { return totalSmoke_; }
    float lastDt() const { return lastDt_; }
    // Cell-centred velocity (average of the two face values per axis) [m/s].
    Vector3 cellVelocity(int i, int j, int k) const;
    // The research layers of the gas on one slice of cells (GasDebugDraw.cpp): the grid, the
    // velocity, -grad p / rho, div u, curl u and the current density - only the layers that are on.
    // axis: the slice's normal (0 x, 1 y, 2 z); layer: the cell index along it.
    void drawDebug(int axis, int layer) const;
    // div u of cell (i, j, k) from its six faces [1/s]: what the pressure projection drives to zero.
    float cellDivergence(int i, int j, int k) const;
    float maxVelocity() const { return maxVel_; }
    // Seconds simulated, summed in double: a float sum of 10^4 steps drifts by ~10^-4 s, enough to
    // show in a manufactured solution's error (verification/Mms.cpp).
    double time() const { return time_; }
    int lastViscousIterations() const { return lastViscousIters_; } // the most any component needed
    bool hasObstacle() const { return solidCount_ > 0; }
    // Running length of the boundary layer on the static obstacle (its extent along x) [m]: the
    // length of the wall-function Reynolds number.
    float wallLength() const { return staticLength_; }
    float dynamicPressure() const { return 0.5f * params.fluidDensity * sqr(std::max(params.inflowSpeed, 1e-3f)); }

private:
    size_t cidx(int i, int j, int k) const { return size_t(i) + size_t(nx_) * (size_t(j) + size_t(ny_) * size_t(k)); }
    Vector3 sampleVelGrid(const Vector3& gp) const; // velocity in cells/s
    Vector3 backtrace(const Vector3& gp, float dt) const;
    // scalar = true (smoke, heat): what comes in through an open (inflow/outflow) side is clean
    // ambient gas (value 0), not a copy of the boundary cell.
    void advect(const Field3& src, Field3& dst, float dt, bool scalar = false) const;
    void advectField(Field3& f, Field3& tmp1, Field3& tmp2, float dt, bool scalar = false);
    // Smoke, heat, fuel ...: cell-centred scalars advected together (one trace for all).
    void advectScalars(const std::vector<Field3*>& fields, float dt);
    std::vector<Vector3> departure_, arrival_;
    std::vector<uint8_t> fromOutside_;
    void applyVelocityBC();
    void addForces(float dt);
    void diffuse(float dt);
    // The steps of diffuse() (Viscosity.cpp), for one velocity component c on its field F.
    bool faceCells(int c, int i, int j, int k, size_t& before, size_t& after) const;
    void markViscousUnknowns(int c, const Field3& F);
    void buildViscousSystem(int c, const Field3& F, float a, float dt);
    void addViscousNeighbour(int c, const Field3& F, const int at[3], int d, int s, float a, double& diag, double& b) const;
    double oldPressureStep(int c, int i, int j, int k, float dt) const;
    void viscousMatVec(const Field3& F, float a, const std::vector<double>& x, std::vector<double>& out) const;
    int solveViscousPcg(Field3& F, float a);
    void project(float dt);
    void computeForces();
    void computeDiagnostics();
    void injectSources();
    void computeDiag();
    // The steps of reset() (GasSolver.cpp).
    void chooseGrid();
    void allocateFields();
    void voxelizeStaticObstacle(const MeshBVH* obstacle);
    void measureObstacle();
    void findStaticWallCells();
    // The steps of step() (GasSolver.cpp).
    float chooseTimeStep(float maxDt);
    void advanceFlowByProjection(float dt); // advection, forces, viscosity, one projection at the end
    void advanceFlowByReflection(float dt); // the same with a reflection at mid-step (AdvectionReflection.cpp)
    void advectAll(float dt);
    void advectVelocity(float dt);       // the three velocity components, carried by (u0_, v0_, w0_)
    void advectCarriedScalars(float dt); // smoke, temperature, fuel, products, carried by (u0_, v0_, w0_)
    void reflectVelocity();              // u <- 2 u - u_advected; the carrier u0 <- 2 u - u_start
    void rubWalls(float dt);             // applyWallFriction, timed with the solids
    void projectTimed(float dt);         // project, timed as the pressure
    void burn(float dt);
    void applyForces(float dt);
    void induceMagneticField(float dt);
    void reportStep(float dt) const;
    void dissipateScalars(float dt);
    void averageCoefficients(float dt);
    // The steps of project() (PressureSolver.cpp).
    void buildPressureRightHandSide(float dt);
    void removeMeanDivergence();
    void buildWeightedDiagonal(bool weighted);
    void solvePressurePcg(bool weighted);
    bool solvePressureMultigridPcg(bool weighted); // MGPCG; false: it broke down, Jacobi takes over
    void fillMultigridFaces(bool weighted);
    void applyPressureMatrix(const std::vector<double>& x, std::vector<double>& out, bool weighted) const;
    void subtractPressureGradient(bool weighted);
    double faceWeightU(int i, int j, int k, bool weighted) const;
    double faceWeightV(int i, int j, int k, bool weighted) const;
    double faceWeightW(int i, int j, int k, bool weighted) const;
    // The steps of voxelizeMovingSolids() and applyWallFriction() (MovingSolids.cpp).
    void collectSolidCells(std::vector<std::vector<size_t>>& cells);
    void markSolidCells(const std::vector<std::vector<size_t>>& cells, std::vector<size_t>& entered, bool& changed);
    void pushScalarsOutOfSolids(const std::vector<size_t>& entered);
    void collectWallCells(std::vector<size_t>& wallCells);
    struct WallFrictionChange { int comp; size_t face; float dv; }; // a velocity change of one face
    void rubWallCell(size_t c, float dt, std::vector<WallFrictionChange>& changes);
    // The 8 cell centres around a world point with trilinear weights, gas cells only (weights
    // renormalised). If none of them is gas, the nearest gas cell within 2 cells (weight 1).
    int fluidStencil(const Vector3& world, size_t cells[8], float weights[8]) const;
    void voxelizeMovingSolids();
    void computeMovingForces();
    void applyWallFriction(float dt);
    void applyPendingImpulses();
    std::vector<std::pair<Vector3, Vector3>> pendingImpulses_;
    std::vector<GasEmission> pendingEmissions_;
    void applyPendingEmissions();
    struct Radiator {
        Vector3 position;
        float power; // [W]
    };
    std::vector<Radiator> radiators_; // hot cells of the last step (fire)
    void collectRadiators();
    void conductHeat(float dt);
    void addHeatSource(float dt); // heatSource over one step (ExternalFields.cpp)
    std::vector<Vector3> liquidPos_, liquidVelIn_; // liquid particles (input)
    std::vector<Vector3> liquidVel_;               // per cell: mean velocity of its liquid particles
    int liquidCells_ = 0;
    bool isLiquidCell(size_t c) const { return (solid_[c] & kMoving) && owner_[c] >= int(moving_.size()); }
    // Boundary value of velocity component `comp` on the face at world point fw between cells ca
    // and cb (-1 = outside): false if neither is solid. Static solids: 0; moving: body velocity.
    bool solidFaceVelocity(long ca, long cb, const Vector3& fw, int comp, float& out) const;

    int nx_ = 0, ny_ = 0, nz_ = 0;
    float dx_ = 1;
    Vector3 origin_;
    Field3 u_, v_, w_, u0_, v0_, w0_, t1_, t2_;
    Field3 p_, smoke_, temp_, fuel_, products_; // products: burnt gas, oxygen used up (Combustion)
    std::vector<float> expansion_; // divergence of the burning gas per cell [1/s] (Combustion)
    Field3 weightU_, weightV_, weightW_; // pressure-equation face weights (Boris correction), see project()
    std::vector<double> diagW_;          // diagonal of the (weighted) pressure matrix
    float heatReleaseRate_ = 0;
    static constexpr uint8_t kStatic = 1, kMoving = 2;
    std::vector<uint8_t> solid_;           // kStatic | kMoving bits
    std::vector<int> owner_;               // moving solid index per cell, -1 if none
    std::vector<MovingSolid> moving_;
    std::vector<Vector3> movForce_, movTorque_;
    std::vector<Vector3> movFricForce_, movFricTorque_;
    std::vector<std::vector<size_t>> solidCells_; // cells of every moving solid (last voxelisation);
                                                  // the list after the bodies' is the liquid's
    std::vector<size_t> staticWallCells_;         // gas cells next to the static obstacle
    std::vector<uint32_t> stamp_;                 // per-cell marks for set operations
    uint32_t stampId_ = 0;
    Vector3 staticFriction_;
    float staticLength_ = 1.0f; // characteristic length of the static obstacle
    int movingCells_ = 0;
    float solidMs_ = 0, pressureMs_ = 0;
    bool hadMoving_ = false;
    std::vector<uint8_t> diag_;
    // Connected fluid regions (6-neighbourhood). A region without an open (outflow) boundary is a
    // pure Neumann problem: the pressure equation is solvable only if its net inflow is zero, so
    // the mean divergence is removed per region (a pocket trapped under a landed body is one).
    std::vector<int> region_;
    std::vector<char> regionOpen_;
    std::vector<double> regionMeanB_; // removed mean of -div per closed region (last projection)
    std::vector<double> q_, r_, z_, s_, As_, b_;
    PressureMultigrid multigrid_; // the levels of the multigrid preconditioner (Multigrid.h)
    // The viscous solve of one component (Viscosity.cpp), one value per face: 1 = unknown (between
    // two gas cells), 0 = kept (a domain side, next to a solid); the matrix diagonal, the right-hand
    // side, the old pressure's velocity step, and the conjugate gradient's vectors.
    std::vector<uint8_t> viscKind_;
    std::vector<double> viscDiag_, viscB_, viscG_, viscX_, viscR_, viscZ_, viscP_, viscAp_;
    int lastViscousIters_ = 0;
    // Advection-reflection's copies (AdvectionReflection.cpp): the flow at the start of the step, the
    // flow after the first half of advection, and the pressures of its two projections - the
    // mid-step one and the final one, each kept to start the same projection of the next step.
    Field3 uStart_, vStart_, wStart_, uAdvected_, vAdvected_, wAdvected_, pHalf_, pEnd_;
    int solidCount_ = 0;

    Vector3 force_;
    float frontalArea_ = 0, planformArea_ = 0, refArea_ = 0;
    float cd_ = 0, cl_ = 0, cs_ = 0, cdAvg_ = 0, clAvg_ = 0;
    int lastIters_ = 0;
    float lastResidual_ = 0, maxVel_ = 0;
    double time_ = 0;
    float maxDivergence_ = 0, totalSmoke_ = 0, lastDt_ = 0;
};

} // namespace rf
