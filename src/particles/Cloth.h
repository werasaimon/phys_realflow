#pragma once
// Woven cloth as a grid of particles: the grid lines along x are the warp threads, along y the weft.
//
// Dynamics - XPBD (Macklin, Mueller, Chentanez 2016, "XPBD: Position-Based Simulation of Compliant
// Constrained Dynamics"): distance constraints for the threads (stretch), the cell diagonals
// (shear) and every second particle (bending). The threads are elastic with the fabric's real
// tensile stiffness (N/m per unit strain), so a thread's tension is simply k * elongation. For
// that to be the real tension the threads must be solved to convergence: every warp and weft line
// is a chain, solved exactly in one tridiagonal sweep (Thomas algorithm; the direct solve of
// chains of Servin & Lacoursiere 2008 and Deul et al. 2018) - one Gauss-Seidel pass would leave
// the solver's lag in the threads as stretch, read as tension a hundred times the real one.
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
    // Shear and bending as cotton has them. Shear: the two diagonal links of a square cell, each of
    // stiffness k, store k s^2 gamma^2 / 2 at a shear angle gamma - the energy of a shear rigidity
    // G = k per unit area, so shearStiffness is G [N/m per radian]. Kawabata's KES-F measures G
    // = 0.5-3 gf/(cm deg) for cotton fabrics (1 gf/(cm deg) = 56 N/m per radian): 28-170 N/m;
    // 100 is a medium cotton (S. Kawabata, "The Standardization and Analysis of Hand Evaluation",
    // Textile Machinery Society of Japan 1980; woven-fabric models fitted to such measurements:
    // Clyde, Teran, Tamstorf 2017, "Modeling and data-driven parameter estimation for woven fabrics").
    float shearStiffness = 100.0f;     // N/m per radian: shear rigidity G
    // Bending: a link over two cells (every second particle), compliance c. It resists only by
    // getting shorter, so its stiffness grows with the fold angle theta: it acts like a bending
    // rigidity B = s^2 theta^2 / (16 c) (s the particle spacing). With 0.1 and s = 15 mm that is
    // 1.3e-5 N m at a fold of 0.3 rad per particle - within KES-F's B = 0.01-0.2 gf cm^2/cm
    // (1e-6 - 2e-5 N m) for cotton. (1e-3, the old default, was a hundred times stiffer: card.)
    float bendCompliance = 0.1f;       // m/N
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
    float force = 0;       // N: the tension it ended the last solve with, -lambda / dt^2 (threads only)
    Kind kind = Warp;
    int16_t x = 0, y = 0;  // grid position of the smaller corner
    bool broken = false;
};

// One thread of the weave: the warp constraints of one row (or the weft constraints of one column)
// in order along it, entries [first, first + count) of Cloth::lineConstraints.
struct ThreadLine {
    int first = 0, count = 0;
};

// One row of a thread line's system (work space of solveCloth, see Cloth.cpp). Double: a stiff
// thread's alpha~ is ~1e-7 of the inverse masses beside it, below what a float resolves.
struct ThreadRow {
    Vector3 n;         // direction of the constraint, from its particle b to a (zero: not solved)
    double diag = 0;   // w_a + w_b + alpha~
    double upper = 0;  // coupling to the next constraint of the line: -w_shared n_k . n_k+1
    double rhs = 0;    // -C - alpha~ lambda; after the solve: the change of lambda
};

// A particle may not get farther from a pinned particle than the cloth between them allows.
struct Tether {
    int particle = 0, anchor = 0;
    float maxLength = 0;
};

struct Cloth {
    int object = -1;          // particle object id
    int group = -1;           // particle group (ParticleSystem::removeGroup removes it whole)
    int firstParticle = 0;    // particle (x, y) = firstParticle + x + width * y
    int width = 0, height = 0;
    float spacing = 0;
    float particleArea = 0;   // fabric one particle stands for [m^2]: its mass = areaDensity * this
    // The width of fabric one thread stands for [m], the same share as the particles' area: a warp
    // thread (a row) the sheet's height / rows, a weft thread (a column) its width / columns. Its
    // stiffness and strength are the fabric's per metre times this (buildClothConstraints).
    float warpWidth = 0, weftWidth = 0;
    ClothMaterial material;
    // Constraints sorted into independent batches (no shared particle within a batch), batch b =
    // constraints[batchStart[b] .. batchStart[b + 1]): a batch is solved in parallel, the batches
    // one after another - Gauss-Seidel convergence, all cores busy.
    std::vector<Vector3> rest;  // rest positions of the grid (x + width * y)
    std::vector<DistanceConstraint> constraints;
    std::vector<int> batchStart;
    // The threads (warp and weft constraints, batches 0-3) are solved line by line, exactly: the
    // warp lines (one per row) first, then the weft lines from weftLinesStart (one per column).
    std::vector<ThreadLine> threadLines;
    int weftLinesStart = 0;
    std::vector<int> lineConstraints;                  // constraint indices, line after line
    std::vector<ThreadRow> threadRows;                 // work space, one row per lineConstraints entry
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
    int smallSteps = 0;       // the small steps the last particle step was cut into for this cloth
    // Work space of ParticleSystem::stepClothsInSmallSteps, per particle: where the small step
    // starts, the velocity, and where a pinned or grabbed particle must end. Kept between steps,
    // so a step allocates nothing.
    std::vector<Vector3> stepStart, stepVelocity, stepEnd;
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

// The cloth's particles moved down by `offset` in the particle arrays (particles before it were
// removed; its own stay together and in order): every global index it keeps follows.
void shiftCloth(Cloth& cloth, int offset);

// Tethers from every particle to every pinned one, lengths measured along the cloth.
void buildTethers(Cloth& cloth, const std::vector<float>& invMass);

// How solveCloth treats the threads. ExactLines: every warp and weft line solved exactly (the
// cloth's own small steps, which integrate its motion). Projection: one Gauss-Seidel projection
// per thread, like shear and bending (the contact passes over the whole step, which only
// re-satisfy the cloth after contacts moved some of its particles; an exact solve there, along
// directions a whole step old, overshoots under tension).
enum class ThreadSolve { ExactLines, Projection };

// One pass over the cloth's constraints - the threads, then shear and bending batch by batch -
// then the tethers.
void solveCloth(Cloth& cloth, std::vector<Vector3>& p, const std::vector<float>& invMass, float dt, ThreadSolve threads);

// The longest cloth small step in which the threads' tension still pulls them straight stably.
// A tensioned thread is a string: bent sideways it springs back, the kink running along it as a
// wave at c = sqrt(T / rho_l) (rho_l = m / l, its mass per length). The solver applies that pull
// explicitly - from the positions at the start of the step - which is stable while the wave
// crosses at most one link per step, c dt <= l, i.e. dt <= sqrt(m l / T); faster, the correction
// overshoots, the thread zig-zags ever wider and reads a false tension. The tension is the one each
// thread ended the last solve with (DistanceConstraint::force); infinity while nothing pulls.
float stableClothStep(const Cloth& cloth, const std::vector<float>& invMass);

// Tension of a thread now [N].
float threadTension(const DistanceConstraint& c, const std::vector<Vector3>& p, float dt);

// The width of fabric a constraint stands for [m]: warpWidth or weftWidth for a thread, the
// spacing for the other links. Tension / this = tension per unit width [N/m].
float threadWidth(const Cloth& cloth, const DistanceConstraint& c);

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
