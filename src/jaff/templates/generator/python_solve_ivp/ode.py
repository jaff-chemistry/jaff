import math

import numpy as np
from commons import *


def get_ode(y, tgas, crate, av, evolve_tgas=False):
    dy = np.zeros_like(y)

    # species number densities (rhs expressions index them as nden[i])
    nden = y[:nspecs]

    # full right-hand side: species dn_i/dt (0..nspecs-1) plus the gas-temperature
    # rate dT/dt (idx_tgas), generated from the network's EOS (heating/cooling and
    # the composition term from reactions changing the particle number)
    rhs = np.zeros(nvars)

    # $JAFF REPEAT idx, rhs, cse IN rhses $[THERMAL dtdt]$
    x$idx$ = $cse$
    rhs[$idx$] = $rhs$
    # $JAFF END

    dy[:nspecs] = rhs[:nspecs]

    # thermal coupling: evolve the gas temperature only when requested at runtime
    if evolve_tgas:
        dy[idx_tgas] = rhs[idx_tgas]

    return dy
