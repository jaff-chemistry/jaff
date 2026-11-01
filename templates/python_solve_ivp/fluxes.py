import numpy as np
from commons import *
from rates import get_rates


def get_fluxes(y, tgas, crate, av):
    k = get_rates(tgas, crate, av)

    flux = np.zeros_like(y)

    $$JAFF REPEAT idx, flux_expression IN flux_expressions

    flux[$idx$] = $flux_expression$

    $$JAFF END

    return flux
