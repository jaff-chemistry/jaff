# ABOUTME: Driver of the generated solver: edit the block of inputs, then run
# ABOUTME: `python main.py`; prints the final state and saves the evolution
"""Integrate the network's chemistry from the inputs below.

Edit the block, then run ``python main.py``. A network that calls functions
``chemistry.py`` does not define (table interpolators and the like) needs them
defined there first; see its docstring.
"""

import numpy as np
from solver import solve

# ---- inputs: edit here ---------------------------------------------------------
N0 = {  # initial number densities [cm^-3]
    $$JAFF REPEAT idx, specie IN species
    "$specie$": 0.0,
    $$JAFF END
}
PARAMS = {  # placeholders: set physical values
    $$JAFF REPEAT idx, parameter IN parameters { THERMAL { MODE -> dtdt } }
    "$parameter$": 1.0,
    $$JAFF END
}
TEMPERATURE = 100.0  # initial gas temperature [K]
T_END = 3.15576e13  # end time [s] (1 Myr)
EVOLVE_TEMPERATURE = False  # True: evolve the temperature with the EOS
OUTPUT = "solution.dat"  # time, temperature and densities, one row per output time
# --------------------------------------------------------------------------------


def main():
    t_eval = np.concatenate(([0.0], np.geomspace(1.0, T_END, 100)))
    result = solve(
        N0,
        TEMPERATURE,
        T_END,
        PARAMS,
        evolve_temperature=EVOLVE_TEMPERATURE,
        t_eval=t_eval,
    )

    names = list(result.n)
    table = np.column_stack([result.t, result.temperature, *result.n.values()])
    np.savetxt(OUTPUT, table, header=" ".join(["t[s]", "T[K]", *names]))

    print(f"t = {result.t[-1]:.6e} s, T = {result.temperature[-1]:.6e} K")
    for name in names:
        print(f"  {name:>10s} {result.n[name][-1]:.6e} cm^-3")
    print(f"evolution saved to {OUTPUT}")


if __name__ == "__main__":
    main()
