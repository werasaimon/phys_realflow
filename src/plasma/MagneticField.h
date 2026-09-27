#pragma once
// Electromagnetic field of a conducting gas (plasma) or liquid metal: resistive, incompressible
// magnetohydrodynamics (MHD) on the staggered grid of GasSolver. SI units.
//   Ampere (low-frequency / MHD limit, no displacement current):  J = curl B / mu0
//   Ohm's law in the moving conductor:                              E = -u x B + J / sigma
//   Faraday (induction):                                            dB/dt = -curl E
//   Lorentz force on the gas (per volume):                           f = J x B
// Together: the field lines are frozen into the flow (they are carried, stretched and twisted by
// it) and slip through it by resistive diffusion, eta = 1 / (mu0 sigma); the field pushes back by
// its magnetic pressure B^2 / 2mu0 and the tension of its lines B^2 / mu0 - Alfven waves travel
// along the lines at v_A = B / sqrt(mu0 rho). Joule heating J^2 / sigma warms the conductor.
//
// Numerics: constrained transport (Evans & Hawley 1988, ApJ 332). B lives on the cell faces (as
// the velocity), E and J on the cell edges; the flux through a face changes only by the
// circulation of E around its four edges (Faraday's law in integral form), so div B stays exactly
// zero, to rounding. The electric field on an edge gets an upwind-type dissipation proportional to
// the local flow speed |u| (as the averaged upwind EMFs of Balsara & Spicer 1999), which keeps the
// central differences of the advection stable; the Alfven waves need none (the Lorentz force and
// Faraday's law are stepped one after the other). The domain walls are perfect conductors: the tangential E
// on them is zero, so the flux through every wall stays what it was.

#include "gas/Field3.h"
#include "math/Math.h"

#include <cstdint>
#include <functional>
#include <vector>

namespace rf {

class MagneticField {
public:
    static constexpr float kMu0 = 1.25663706e-6f; // vacuum permeability [T m / A]

    bool enabled = false;
    // [T] uniform field of coils around the vessel. It is part of the fixed background B0, not of
    // the evolving B1: a weak wave on a strong field then keeps its own float precision (a 0.01 T
    // field stored in one float with a 3e-5 T wave left the wave three significant digits).
    Vector3 applied{0.0f};
    float conductivity = 1e6f;       // sigma [S/m]: resistivity eta = 1 / (mu0 sigma) [m^2/s]
    // Dissipation of the centred advection of B on an edge: eta_num = max(f |u| dx, |u|^2 h) with
    // this factor f and the substep h. The second term is twice the Lax-Wendroff amount
    // |u|^2 h / 2 that keeps a centred, forward-stepped advection stable - the least there is
    // (with f = 0 that is all); the first adds a margin independent of the step (f = 0.5: the
    // averaged upwind EMF of Balsara & Spicer 1999). Where the flow is slow next to the Alfven
    // speed the first is many times the second, and eats a current profile in seconds.
    float numericalDissipation = 0.5f;
    // Boris correction (Boris 1970; Gombosi et al. 2002, "Semirelativistic MHD and the Boris
    // correction", as in magnetosphere codes): the field's inertia caps the Alfven speed at this
    // reduced "speed of light", v_A' = v_A / sqrt(1 + v_A^2 / c^2), so a strong field (near a
    // magnet's poles) does not force tiny time steps. Force balances (steady states) do not change,
    // only how fast the strong-field regions respond. 0 = off.
    float speedLimit = 0.0f; // [m/s]

    void reset(int nx, int ny, int nz, float dx, const Vector3& origin);
    // Cells that are perfect conductors (the vessel wall): no tangential E on their edges.
    void setConductors(const std::vector<uint8_t>& conductor);
    // Adds curl A to the field (A [T m] at a world point): a divergence-free start field (e.g. the
    // field of a current column, from its vector potential).
    void addFromPotential(const std::function<Vector3(const Vector3&)>& A);
    // A fixed, current-free background field B0 = curl A (a magnet, coils) - the "B0 + B1 split" of
    // magnetosphere codes (Tanaka 1994): only the induced part B1 evolves and only its curl is the
    // current. On the grid the curl of a steep field such as a dipole's (~1/r^3) is not exactly
    // zero; taken as current it would drive spurious flows next to the magnet.
    void setBackgroundFromPotential(const std::function<Vector3(const Vector3&)>& A);
    // Resistivity per cell [m^2/s] (world point of the cell centre) instead of the uniform
    // 1 / (mu0 sigma): the resistive "vacuum" around a tokamak's plasma, a cold edge. On an edge
    // the mean of its 4 cells. Set after reset(); reset() clears it.
    void setResistivityMap(const std::function<float(const Vector3&)>& eta);
    // The largest resistivity of any cell (the uniform one without a map): the diffusion limit.
    float maxResistivity() const { return etaCell_.empty() ? resistivity() : etaMax_; }

    // Faraday's law over dt with the face velocities u, v, w [m/s] (the gas at `density`), in as
    // many substeps as the Alfven speed and the diffusion need.
    void induce(const Field3& u, const Field3& v, const Field3& w, float density, float dt);
    // Lorentz force J x B: velocity change of the inner faces (solid[c] cells skipped).
    void applyLorentzForce(Field3& u, Field3& v, Field3& w, const std::vector<uint8_t>& solid, float density, float dt);
    // Joule heating J^2 / sigma of every cell over dt [J/m^3]: the heat of the free edges (a wall's
    // edges do no work), a quarter of each edge's to each of its four cells - the exact partner of
    // the field energy the resistive step removes.
    void jouleHeating(std::vector<float>& heat, float dt);
    // Boris correction: 1 / (1 + v_A^2 / c^2) on every face (x, y, z faces as u, v, w) - the factor
    // by which the field's inertia slows every acceleration there (the pressure's too: see
    // GasSolver::project). 1 everywhere when speedLimit is 0.
    void borisWeights(Field3& wx, Field3& wy, Field3& wz, float density) const;
    // Largest stable time step of the gas for the Alfven waves (the fastest signal).
    float maxTimeStep(float maxFlowSpeed, float density) const;

    // Magnetic field B [T], current density J [A/m^2] and electric field E [V/m] (of the last step)
    // at a world point (trilinear) or a cell centre.
    Vector3 fieldAt(const Vector3& world) const;
    Vector3 currentAt(const Vector3& world) const;
    Vector3 electricFieldAt(const Vector3& world) const;
    Vector3 cellField(int i, int j, int k) const;
    Vector3 cellCurrent(int i, int j, int k) const;
    Vector3 cellElectricField(int i, int j, int k) const;

    // Diagnostics
    float maxField() const;       // [T]
    double energy() const;        // magnetic energy, integral of B^2 / 2mu0 [J]
    float maxDivergence() const;  // max |div B| dx / max |B| (0 up to rounding)
    float resistivity() const { return 1.0f / (kMu0 * conductivity); }

    // The evolving face fields (x, y, z components on the x, y, z faces): the whole field, or the
    // induced part B1 when a background B0 is set.
    Field3 bx, by, bz;

private:
    // Alfven speed of field B2 = |B|^2 in gas of `density`, capped by the Boris correction.
    float alfvenSpeed(float B2, float density) const {
        const float vA2 = B2 / (kMu0 * std::max(density, 1e-12f));
        return speedLimit > 0 ? std::sqrt(vA2 / (1.0f + vA2 / (speedLimit * speedLimit))) : std::sqrt(vA2);
    }
    void computeCurrent();                                                // J = curl B / mu0 on the edges
    void computeElectricField(const Field3& u, const Field3& v, const Field3& w, float substep);
    void applyFaraday(float dt);                                          // B -= dt curl E
    int nx_ = 0, ny_ = 0, nz_ = 0;
    float dx_ = 1;
    Vector3 origin_;
    Field3 ex_, ey_, ez_; // edges: x-edges at (i+1/2, j, k), y-edges at (i, j+1/2, k), z-edges at (i, j, k+1/2)
    Field3 jx_, jy_, jz_;
    Field3 b0x_, b0y_, b0z_;         // background B0 (empty = none)
    std::vector<float> etaCell_;     // resistivity map per cell (empty = uniform)
    Field3 etaX_, etaY_, etaZ_;      // the map on the edges
    float etaMax_ = 0;
    Field3 tbx_, tby_, tbz_;         // total B = B0 + B1, refreshed before it is used
    void updateTotal();
    void initBackground(); // B0 = the uniform `applied` field on every face
    void addCurl(const std::function<Vector3(const Vector3&)>& A, Field3& fx, Field3& fy, Field3& fz) const;
    std::vector<uint8_t> conductor_; // per cell, empty = none
    std::vector<uint8_t> fixedEx_, fixedEy_, fixedEz_; // edges on a conductor: E = 0
    void markConductorEdges();
    bool conductor(int i, int j, int k) const {
        if (i < 0 || j < 0 || k < 0 || i >= nx_ || j >= ny_ || k >= nz_) return true; // the domain walls
        return !conductor_.empty() && conductor_[size_t(i) + size_t(nx_) * (size_t(j) + size_t(ny_) * size_t(k))] != 0;
    }
};

} // namespace rf
