// Fields set from outside the gas solver, for verification: the velocity and the temperature set
// from a function of the world point, and a heat source term added every step. With them the
// solver can start from an exact flow (the Taylor-Green vortex) or run a manufactured solution:
// choose T(x, t), compute the source S = dT/dt - div(alpha grad T) by hand, add S every step and
// the solver must reproduce T to its order of accuracy (Roache 2002, "Code verification by the
// method of manufactured solutions", J. Fluids Eng. 124). Nothing here is used by the scenes.
#include "gas/GasSolver.h"

#include "core/Parallel.h"

namespace rf {

void GasSolver::setVelocity(const std::function<Vector3(const Vector3&)>& velocity) {
    // Each component at its own face centre (the staggered grid): u on the x faces, v on y, w on z.
    Field3* components[3] = {&u_, &v_, &w_};
    for (int c = 0; c < 3; ++c) {
        Field3& F = *components[c];
        parallelFor(F.nz, [&](int k) {
            for (int j = 0; j < F.ny; ++j)
                for (int i = 0; i < F.nx; ++i) {
                    const Vector3 x = origin_ + (Vector3(float(i), float(j), float(k)) + F.offset) * dx_;
                    F.at(i, j, k) = velocity(x)[c];
                }
        }, 1);
    }
    applyVelocityBC(); // the walls and the solid faces keep their own values
}

void GasSolver::setTemperature(const std::function<float(const Vector3&)>& temperature) {
    parallelFor(nz_, [&](int k) {
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                const size_t c = cidx(i, j, k);
                temp_.d[c] = solid_[c] ? 0.0f : temperature(origin_ + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * dx_);
            }
    }, 1);
}

// The source is added after the conduction of the step, at the middle of the step's time
// interval: T += dt S(x, t + dt/2) in every gas cell (operator splitting, first order in dt).
void GasSolver::addHeatSource(float dt) {
    if (!heatSource) return;
    const double tMid = time_ + 0.5 * double(dt);
    parallelFor(nz_, [&](int k) {
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                const size_t c = cidx(i, j, k);
                if (!solid_[c]) temp_.d[c] += dt * heatSource(origin_ + Vector3(i + 0.5f, j + 0.5f, k + 0.5f) * dx_, tMid);
            }
    }, 1);
}

} // namespace rf
