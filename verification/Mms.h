#pragma once
// The method of manufactured solutions (Salari & Knupp 2000, SAND2000-1444; Roache 2002, J. Fluids
// Eng. 124): choose a smooth field, put it into the equation, and whatever is left over is the
// source term that makes it an exact solution. Add that source to the solver and it must
// reproduce the chosen field, with an error that falls as h^p - p the order of the scheme. Any
// bug in a term of the discrete equation shows up as a wrong order.
//
// Here: the heat equation of the gas solver as it is coded, with the diffusivity of air that
// grows with temperature, alpha(T) = alpha0 ((T + Ta) / Ta)^1.75 (T is the excess over the
// ambient Ta), in a closed box with insulated walls (dT/dn = 0):
//     dT/dt = div(alpha(T) grad T) + S.
// The manufactured field breathes in time,
//     T = A (2 + phi g),   phi = cos(pi X) cos(pi Y),   g = 1 + sin(w t) / 2,   X = x / Lx, Y = y / Ly,
// its normal derivative vanishes on the walls, and the source is
//     S = dT/dt - alpha(T) lap T - alpha'(T) |grad T|^2
//       = A phi w cos(w t) / 2 + alpha(T) A g k^2 phi - alpha'(T) A^2 g^2 |grad phi|^2,
//     k^2 = pi^2 (1/Lx^2 + 1/Ly^2).
// The nonlinear term is verified along with the rest. The box is thin in z (4 cells) and the
// field does not depend on z: a 2D problem in the 3D code.
//
// Why it breathes: a steady manufactured field started at its exact value needs, on fine grids,
// corrections per conduction substep below one float ulp of T - the solver then does not move
// at all and the error looks smaller than it is (a finding of this case, docs/10). A field that
// changes every step keeps the rounding unbiased.

#include "Benchmark.h"

namespace rf::verify {

struct ManufacturedHeat {
    double Lx = 1, Ly = 1;      // the box [m] (its lower corner at the origin)
    double amplitude = 1.0;     // A [K]
    double alpha0 = 0.01;       // diffusivity at ambient [m^2/s]
    double ambient = 293.0;     // Ta [K], as Combustion::ambientTemperature
    double temperature(double x, double y, double t) const;
    double source(double x, double y, double t) const; // [K/s]
    double diffusionTime() const;                       // 1 / (alpha0 k^2) [s]
    double omega() const;                               // one period per diffusion time [1/s]
};

struct MmsLevel {
    int cells = 0;       // across x
    double dx = 0;
    double errorRms = 0; // RMS of T - T_exact over the cells [K], at the end
    double errorMax = 0;
    int steps = 0;
};

// One grid of the manufactured heat problem: `cells` across x, run for one diffusion time with
// steps small enough for the conduction to take one substep each (dt ~ dx^2).
MmsLevel runManufacturedHeat(const ManufacturedHeat& m, int cells);

} // namespace rf::verify
