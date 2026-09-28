// Advection-reflection: the flow over one step with the projection turned into a mirror at mid-step
// (GasParams::advectionReflection).
//
// Why. The classic step (advanceFlowByProjection) carries the velocity along itself and then
// projects it, P(u) = u - grad p: the divergent part of u is thrown away. Carried by itself, a
// divergence-free flow picks up a divergent part every step, so every step some of its kinetic
// energy leaves with that part - a numerical viscosity on top of the real one, the 13.8 % of the
// energy loss the Taylor-Green case measured (docs/11-action.md). Zehnder, Narain & Thomaszewski
// ("An advection-reflection solver for detail-preserving fluid simulation", ACM TOG 37(4), SIGGRAPH
// 2018, doi 10.1145/3197517.3201324) noticed that P is an orthogonal projection in the kinetic-energy
// norm, so R = 2P - I is a reflection: it removes the divergent part not by cutting it off but by
// turning it around, and a reflection keeps the energy exactly. Done at mid-step, the divergent
// part picked up in the first half is turned back and largely cancelled by the second half; only
// the second, much smaller projection at the end removes energy.
//
// The recipe, second order in time (Narain, Zehnder & Thomaszewski, "A second-order advection-
// reflection solver", Proc. ACM Comput. Graph. Interact. Tech. 2(2), 2019, doi 10.1145/3340257,
// eq. 11 - the implicit midpoint rule with a projection at mid-step):
//   1. u_a      = A[u0; u0, dt/2]               carry the flow half a step by itself
//   2. u_half   = P[u_a + dt/2 f]               the mid-step forces, the mid-step projection
//   3. u_mirror = 2 u_half - u_a                the reflection (= R[u_a] + dt P[f])
//   4. u_end    = A[u_mirror; 2 u_half - u0, dt/2]   carried the second half by the estimate of
//                                                     the step's end flow, 2 u_half - u0
//   5. u1       = P[u_end]                      the final projection
// A[f; v, h] is our semi-Lagrangian (MacCormack) transport of f along v for time h.
//
// What else happens in the step, and where. The scalars (smoke, heat, fuel) ride the flow of the
// step's start for the whole step, as in the classic step, and the fire burns before the flow moves.
// The forces of step 2 are the smooth ones (buoyancy, confinement, Lorentz) over half a step: the
// reflection doubles their push, so the flow receives dt f. Viscosity, the bodies' impulses and the
// wall friction come after step 4, over the whole step, before the final projection: they are the
// dissipative physics and stay as they are.
//
// Two pressures. The mid-step projection's pressure p_half pushes the mirrored flow (by dt, through
// the reflection); the final projection's pressure p_end holds a steady flow against its viscosity
// (a channel's pressure drop). Each is kept to start the same projection of the next step. The
// viscous solve is incremental in p_end, as in the classic step (Viscosity.cpp, item 3): without it
// the final projection moves the gas on the wall by dt |grad p_end| / rho after the no-slip solve,
// and the walls fall back to first order - the Poiseuille case read an order of 0.93 until
// 2026-09-28. The physical pressure of the step, for the forces on the bodies, is p_half + p_end.
#include "gas/GasSolver.h"
#include "gas/GasSolverInternal.h"

#include "core/Parallel.h"
#include "core/Probe.h"


namespace rf {

void GasSolver::advanceFlowByReflection(float dt) {
    const float half = 0.5f * dt;
    // The two kept pressures; after a reset (or a switch from the classic step) the last total
    // pressure is the best guess for the final one.
    if (pEnd_.d.size() != p_.d.size()) {
        pEnd_ = p_;
        pHalf_ = p_;
        std::fill(pHalf_.d.begin(), pHalf_.d.end(), 0.0f);
    }
    // 0. The scalars ride the step's start flow for the whole step; the fire burns.
    {
        Probe::Timer timer("gas/advect ms");
        u0_ = u_; v0_ = v_; w0_ = w_;
        advectCarriedScalars(dt);
        uStart_ = u_; vStart_ = v_; wStart_ = w_;
        // 1. Half a step of the flow carried by itself.
        advectVelocity(half);
        applyVelocityBC();
        uAdvected_ = u_; vAdvected_ = v_; wAdvected_ = w_;
    }
    burn(dt);
    // 2. The smooth forces over half a step, then the mid-step projection.
    {
        Probe::Timer timer("gas/forces ms");
        addForces(half);
        if (magnetic.enabled) magnetic.applyLorentzForce(u_, v_, w_, solid_, params.fluidDensity, half);
    }
    applyVelocityBC();
    p_.d.swap(pHalf_.d); // the mid-step projection starts from last step's mid-step pressure
    projectTimed(half);
    p_.d.swap(pHalf_.d);
    // 3. The reflection, and the flow that will carry the second half.
    reflectVelocity();
    applyVelocityBC();
    // 4. The second half of the advection.
    {
        Probe::Timer timer("gas/advect ms");
        advectVelocity(half);
    }
    applyVelocityBC();
    // Viscosity (incremental in the final pressure, kept in p_ for it), the bodies' impulses and the
    // walls over the whole step.
    p_.d.swap(pEnd_.d);
    {
        Probe::Timer timer("gas/diffuse ms");
        diffuse(dt);
    }
    applyPendingImpulses();
    applyVelocityBC();
    rubWalls(dt);
    // 5. The final projection, from last step's final pressure; the step's pressure is the two
    //    projections' together.
    projectTimed(dt);
    pEnd_.d = p_.d;
    for (size_t c = 0; c < p_.d.size(); ++c) p_.d[c] += pHalf_.d[c];
}

// Step 3 of the recipe, face by face, on the three components: u_half is in (u_, v_, w_).
//   the mirrored flow      u_ <- 2 u_half - u_advected
//   the second-half carrier u0_ <- 2 u_half - u_start   (an estimate of the flow at the step's end)
void GasSolver::reflectVelocity() {
    Field3* now[3] = {&u_, &v_, &w_};
    Field3* carrier[3] = {&u0_, &v0_, &w0_};
    const Field3* start[3] = {&uStart_, &vStart_, &wStart_};
    const Field3* advected[3] = {&uAdvected_, &vAdvected_, &wAdvected_};
    for (int c = 0; c < 3; ++c) {
        std::vector<float>& u = now[c]->d;
        std::vector<float>& carry = carrier[c]->d;
        const std::vector<float>& u0 = start[c]->d;
        const std::vector<float>& ua = advected[c]->d;
        parallelFor(int(u.size()), [&](int f) {
            const float uHalf = u[size_t(f)];
            carry[size_t(f)] = 2.0f * uHalf - u0[size_t(f)];
            u[size_t(f)] = 2.0f * uHalf - ua[size_t(f)];
        }, 4096);
    }
}

} // namespace rf
