#!/usr/bin/env python3
"""The plots of docs/04-gas-navier-stokes.md, drawn from real runs of tests/GasTests.cpp.

Run the test with a folder for its curves, then this script turns the CSV into an SVG:

    RF_PLOT_DIR=plots RF_TEST="advection-reflection" ./rf_tests   (the test writes plots/*.csv)
    python tools/plot_gas.py plots docs/images/gas                 (this script writes the SVG)

It draws with the same plain-SVG helpers as tools/plot_action.py (Python standard library only,
its own light and dark palette, axes with units, a legend, a caption naming the test). Nothing is
drawn by hand: every point is a row of a CSV written by a test.
"""
import os
import sys

import plot_action as plot


def main():
    source, destination = sys.argv[1], sys.argv[2]
    os.makedirs(destination, exist_ok=True)
    run = "RF_PLOT_DIR=… rf_tests, tests/GasTests.cpp"
    # The inviscid Taylor-Green vortex: every joule it loses is lost by the numerics.
    d = plot.read_csv(source, "advection_reflection")
    lost = lambda kept: [100.0 * (1.0 - k) for k in kept]  # energy lost, in % of the start
    plot.line_chart(os.path.join(destination, "advection-reflection.svg"),
                    "Вихрь Тейлора–Грина без вязкости: сколько энергии теряет численная схема", "время, с",
                    "потеряно энергии, %",
                    f"32 ячейки, шаг dt = dx/2, вязкость 0: всё потерянное — численная диссипация. "
                    f"Тест «gas: advection-reflection…» ({run})",
                    d["t_s"], [("проекция в конце шага (Stam 1999)", lost(d["projection"])),
                               ("отражение на полушаге (Zehnder et al. 2018)", lost(d["reflection"]))])
    print("plots written to", destination)


if __name__ == "__main__":
    main()
