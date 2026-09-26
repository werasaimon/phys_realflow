#pragma once
// Woven cloth as a grid of particles: the grid lines along x are the warp threads, along y the weft.
//
// Dynamics - XPBD (Macklin, Mueller, Chentanez 2016, "XPBD: Position-Based Simulation of Compliant
// Constrained Dynamics"): distance constraints for the threads (stretch), the cell diagonals
// (shear) and every second particle (bending). The threads are elastic with the fabric's real
// tensile stiffness (N/m per unit strain), so a thread's tension is simply k * elongation.
// Long range attachments (Kim, Chentanez, Mueller 2012) with a little slack keep a pinned cloth
// from sagging while leaving the threads free to carry - and concentrate - the load.
//
// Tearing: a thread breaks when its tension exceeds its strength (fabric strength in N/m times the
// width the thread stands for). Its load moves to the neighbours, so a crack runs along the weave,
// perpendicular to the load, thread by thread. A cell is cut through (its shear links removed, no
// longer drawn) once the crack has crossed it completely - until then it is the crack tip. Seams
// are stitch lines that hold only part of the fabric strength, so a garment opens there first.

#include "math/Math.h"

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace rf {

struct ClothMaterial {
    float areaDensity = 0.3f;          // kg/m^2 (cotton ~0.15-0.3, canvas ~0.5-1.5)
    float tensileStiffness = 30000.0f; // N/m: tension per unit width per unit strain (0 = inextensible)
    float shearStiffness = 3000.0f;    // N/m: the cells shear far more easily than the threads stretch
    float bendCompliance = 1e-3f;      // m/N
    float strengthWarp = 4000.0f;      // N/m: tension per unit width that breaks the threads along x (0 = never tears)
    float strengthWeft = 3000.0f;      // N/m: along y
    std::vector<int> seamColumns;      // seams: stitched lines between columns c and c+1 ...
    std::vector<int> seamRows;         // ... and between rows r and r+1
    float seamStrength = 0.3f;         // fraction of the fabric strength a seam holds

    // Fire. Heated fabric decomposes (pyrolysis) at the Arrhenius rate k = A exp(-E / (R T)) per
    // unit of what is left; the decomposition takes heat (endothermic) and gives off fuel gas that
    // burns in the air as the flame. There is no ignition temperature to set: the fabric starts to
    // decompose fast where the heat it gets outruns the heat the reaction takes, and the front then
    // sits at the temperature where they balance. Defaults: cellulose, single-step global kinetics
    // of Antal & Varhegyi 1995 ("Cellulose pyrolysis kinetics: the current state of knowledge").
    bool flammable = false;
    float pyrolysisPreExponential = 1.26e18f; // A [1/s]
    float pyrolysisActivationEnergy = 238e3f; // E [J/mol]
    float pyrolysisHeat = 4e5f;               // [J/kg] of fabric decomposed, taken from the fabric
    // m^3 of fuel gas (as stoichiometric mixture with air) per kg decomposed: the volatiles of
    // cellulose release ~12 MJ/kg burning, a stoichiometric mixture ~2.2 MJ/m^3 (1800 K x 1231 J/(m^3 K)).
    float fuelYield = 5.5f;
    float specificHeat = 1300.0f;      // J/(kg K)
    float heatTransfer = 40.0f;        // W/(m^2 K) per side: convection between the gas and the fabric
    float emissivity = 0.9f;           // grey body: absorbs this part of the flame's radiation, radiates e sigma T^4
    float charMassFraction = 0.2f;     // part of the mass left as char once burnt through
    float charStrength = 0.001f;       // part of the strength left once charred through (ash)
    float conductivity = 0.06f;        // W/(m K): heat conduction along the fabric (Fourier) ...
    float thickness = 0.5e-3f;         // m: ... through its cross-section
};

struct DistanceConstraint {
    enum Kind : uint8_t { Warp, Weft, Shear, Bend };
    int a = 0, b = 0;
    float restLength = 0;
    float compliance = 0;  // m/N
    float strength = 0;    // N: tension at which a thread breaks (0 = never)
    float lambda = 0;      // accumulated Lagrange multiplier of this substep
    Kind kind = Warp;
    int16_t x = 0, y = 0;  // grid position of the smaller corner
    bool broken = false;
};

// A particle may not get farther from a pinned particle than the cloth between them allows.
struct Tether {
    int particle = 0, anchor = 0;
    float maxLength = 0;
};

struct Cloth {
    int object = -1;          // particle object id
    int firstParticle = 0;    // particle (x, y) = firstParticle + x + width * y
    int width = 0, height = 0;
    float spacing = 0;
    ClothMaterial material;
    // Constraints sorted into independent batches (no shared particle within a batch), batch b =
    // constraints[batchStart[b] .. batchStart[b + 1]): a batch is solved in parallel, the batches
    // one after another - Gauss-Seidel convergence, all cores busy.
    std::vector<Vector3> rest;  // rest positions of the grid (x + width * y)
    std::vector<DistanceConstraint> constraints;
    std::vector<int> batchStart;
    std::unordered_map<uint64_t, int> constraintIndex; // (smaller, larger particle) -> constraint
    std::vector<Tether> tethers;                       // tethersPerParticle consecutive entries per particle
    int tethersPerParticle = 0;
    std::vector<int> anchors;                          // pinned particles (local indices)
    bool tethersDirty = false;                         // threads broke: re-measure the tethers
    // Cells (w-1)*(h-1): 1 while the cell holds together (drawn), 0 once a crack has cut it.
    std::vector<uint8_t> cellIntact;
    // Broken threads: warp (x, y)-(x+1, y), (w-1)*h, and weft (x, y)-(x, y+1), w*(h-1).
    std::vector<uint8_t> warpBroken, weftBroken;
    int tornThreads = 0;
    // Fire state per particle (x + width * y): temperature [K above ambient] and the part of the
    // fabric not burnt yet (1 fresh .. 0 charred through).
    std::vector<float> temperature, unburnt;
    int burntThreads = 0;
    Vector3 color{0.85f, 0.75f, 0.35f};
    int particle(int x, int y) const { return firstParticle + x + width * y; }
};

// Threads, shear and bending constraints of the grid (rest positions `rest`), in 16 independent
// batches (graph colouring of the grid), stiffnesses and strengths from the material.
void buildClothConstraints(Cloth& cloth, const std::vector<Vector3>& rest);

// Tethers from every particle to every pinned one, lengths measured along the cloth.
void buildTethers(Cloth& cloth, const std::vector<float>& invMass);

// One pass over the cloth's constraints (batch by batch), then the tethers.
void solveCloth(Cloth& cloth, std::vector<Vector3>& p, const std::vector<float>& invMass, float dt);

// Tension of a thread now [N].
float threadTension(const DistanceConstraint& c, const std::vector<Vector3>& p, float dt);

// Breaks the threads loaded beyond their strength and cuts the cells the cracks have crossed.
// Burning weakens the threads: strength x (charStrength + (1 - charStrength) x unburnt fraction).
// Returns the number of threads broken now (the tethers are then marked dirty).
int tearCloth(Cloth& cloth, const std::vector<Vector3>& p, float dt);

// Re-measures the tethers along what is left of a torn cloth (once per step, if dirty).
void updateTethers(Cloth& cloth);

// Fire, one step. Heat flows along the fabric by conduction (dT/dt = div(k grad T) / (rho c),
// through intact threads only - not across a burnt hole). Every patch of fabric (one particle) is
// heated by convection from the gas around
// it (gasTemperature per local particle, K above ambientTemperature [K]), absorbs the flame's
// radiation arriving at it (irradiance, W/m^2) and radiates itself (e sigma (T^4 - T_air^4), both
// faces), and decomposes at the Arrhenius rate of its temperature, which takes the heat of
// pyrolysis (the energy balance is solved implicitly - the kinetics are stiff) and gives off fuel
// gas. A thread whose two ends have both charred through falls apart. Adds per local particle the heat given to the gas [J] and the fuel
// gas released [m^3]. Returns the threads burnt through.
int burnCloth(Cloth& cloth, const std::vector<float>& gasTemperature, const std::vector<float>& irradiance,
              float ambientTemperature, float dt, std::vector<float>& heatToGas, std::vector<float>& fuelToGas);

} // namespace rf
