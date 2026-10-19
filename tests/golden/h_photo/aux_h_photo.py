# ABOUTME: jaffgen test template: every auxiliary function of the h_photo network
# ABOUTME: Rendered via GET aux_func; compared numerically against tests/golden
import math


def aux_functions():
    out = {}
    # $$JAFF GET aux_func IN deltae0
    out["deltae0"] = 6.40000000000000e-12
    # $$JAFF END
    # $$JAFF GET aux_func IN deltae1
    out["deltae1"] = -1.380649e-16*tgas*(0.684 - 0.0416*math.log(0.0001*tgas))
    # $$JAFF END
    return out
