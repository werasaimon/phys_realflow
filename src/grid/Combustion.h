#pragma once
// Combustion of a gaseous fuel in the grid gas - the flame model of Nguyen, Fedkiw, Jensen 2002
// ("Physically Based Modeling and Animation of Fire") in the simplified grid form of Bridson
// ("Fluid Simulation for Computer Graphics", 2nd ed., ch. "Fire"), the model behind Houdini Pyro:
//  * fuel: concentration of fuel gas in a cell (1 = stoichiometric mixture with the air), carried
//    by the flow like smoke;
//  * where the gas is hotter than the ignition temperature the fuel burns at a first-order rate,
//    limited by the oxygen: a cell holds `products` (burnt gas, 0 = fresh air .. 1 = no oxygen
//    left), and at most 1 - products of fuel can still burn there - so a fuel-rich cell burns
//    only as fast as fresh air mixes in and never gets hotter than the adiabatic flame temperature
//    (heatRelease above ambient). The reaction releases heat, soot
//    (smoke) and volume: the burnt gas expands - a divergence source of the pressure projection,
//    which pushes the flame outwards and makes it billow;
//  * heat conduction, dT/dt = div(alpha grad T) (Fourier), with the diffusivity of air;
//  * hot gas cools by radiation, dT/dt = -c (T / 1000 K)^4 (Nguyen et al. eq. 7; optically thin
//    gas, Stefan-Boltzmann), integrated exactly so that large steps stay stable;
//  * hot gas rises with the buoyancy of an ideal gas, g (T - T_air) / T (the Boussinesq form
//    g (T - T_air) / T_air overestimates it several times at flame temperatures).
// Temperatures are kelvin above the ambient air.

#include <cstdint>
#include <vector>

namespace rf {

class Combustion {
public:
    bool enabled = false;
    float ambientTemperature = 293.0f;  // [K] the air
    float ignitionTemperature = 300.0f; // [K above ambient] fuel burns in gas hotter than this
    float burnRate = 6.0f;              // [1/s] fraction of a cell's fuel burning per second
    float heatRelease = 1800.0f;        // [K] temperature rise per unit of fuel burnt
    float sootYield = 0.2f;             // smoke density per unit of fuel burnt
    float expansion = 1.0f;             // volume of burnt gas per volume of fuel burnt
    float radiativeCooling = 1500.0f;   // [K/s] cooling rate of gas 1000 K above ambient (~T^4)
    float absorptionCoefficient = 2.0f; // [1/m] of the sooty flame gas: how strongly it radiates
    float thermalDiffusivity = 2.2e-5f; // [m^2/s] of air at ambient; grows ~T^1.75 (kinetic theory)
    float gravity = 9.81f;              // [m/s^2]
    float specificHeat = 1005.0f;       // [J/(kg K)] of the gas (air)

    // One reaction step over all gas cells (solid[c] != 0 is skipped). expansionRate[c] receives
    // the divergence of the expanding burnt gas [1/s]. Returns the fuel burnt, summed over cells.
    double react(std::vector<float>& fuel, std::vector<float>& products, std::vector<float>& temperature,
                 std::vector<float>& smoke, std::vector<float>& expansionRate, const std::vector<uint8_t>& solid, float dt) const;

    // Upward acceleration of gas `temperature` K above ambient [m/s^2].
    float buoyancy(float temperature) const { return gravity * temperature / (ambientTemperature + temperature); }
};

} // namespace rf
