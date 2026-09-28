#!/usr/bin/env python3
"""The plot of the Alfvénic state in docs/06-mhd-plasma.md, drawn from a real run of
tests/PlasmaTests.cpp (testAlfvenicState).

Run the test with a folder for its curves, then this script turns the CSV into an SVG:

    RF_PLOT_DIR=plots RF_TEST="Alfvénic" ./rf_tests          (the test writes plots/alfvenic_state.csv)
    python tools/plot_mhd.py plots docs/img                   (this script writes the SVG)

The drawing itself (axes, legend, both colour themes) is tools/plot_action.py's line chart, so
every plot of the book looks the same. Python standard library only. Nothing is drawn by hand:
every point is a row the test wrote.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_action import line_chart, read_csv  # noqa: E402  (the shared drawing of the book's plots)


def main():
    src, dst = sys.argv[1], sys.argv[2]
    os.makedirs(dst, exist_ok=True)
    d = read_csv(src, "alfvenic_state")
    # The (2,2) cell is the telling one: the field keeps it on exp(-8 nu pi^2 t); without the
    # field the (1,1) cell stirs it until it turns over (its amplitude goes negative).
    line_chart(os.path.join(dst, "mhd-alfvenic-state.svg"),
               "Альфвеновское состояние: поле не даёт вихрям перемешаться", "время, с", "A(t) / A(0) ячейки (2,2)",
               "Две ячейки конвекции, 32 клетки, ν = η = 0,01 м²/с. Тест «MHD: an Alfvénic state v = b…» "
               "(RF_PLOT_DIR=… rf_tests, tests/PlasmaTests.cpp)",
               d["t_s"], [("МГД, v = b", d["mhd_22"]), ("то же течение без поля", d["hydro_22"]),
                          ("точно: e^(−8νπ²t)", d["exact_22"])])
    print("plot written to", dst)


if __name__ == "__main__":
    main()
