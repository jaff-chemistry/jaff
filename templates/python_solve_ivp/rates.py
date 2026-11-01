import functools
import math

import numpy as np

from commons import nreactions


def get_rates(tgas, crate, av):
    k = np.zeros(nreactions)
    kphoto = np.zeros(nreactions)

    $$JAFF REPEAT idx, rate IN  rates

    k[$idx$] = $rate$

    $$JAFF END
    return k
