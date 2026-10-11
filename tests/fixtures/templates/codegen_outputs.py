# ABOUTME: jaffgen test template exposing every expression generator as a Python function
# ABOUTME: Rendered per network, then evaluated against tests/golden with shared inputs
import math

# $JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign $[POS j NEG k REPLACE idx_ek idx_e]$
idx_$specie_with_normalized_sign$ = $idx$
# $JAFF END

SPECIE_MASSES = {}
# $JAFF REPEAT idx, specie_mass IN specie_masses
SPECIE_MASSES[$idx$] = $specie_mass$
# $JAFF END


def rates():
    k = {}
    # $JAFF REPEAT idx, rate, cse IN rates
    x$idx$ = $cse$
    k[$idx$] = $rate$
    # $JAFF END
    return k


def fluxes():
    k = rates()
    y = nden
    out = {}
    # $JAFF REPEAT idx, flux_expression IN flux_expressions
    out[$idx$] = $flux_expression$
    # $JAFF END
    return out


def ode_expressions():
    flux = fluxes()
    out = {}
    # $JAFF REPEAT idx, ode_expression IN ode_expressions
    out[$idx$] = $ode_expression$
    # $JAFF END
    return out


def odes():
    out = {}
    # $JAFF REPEAT idx, ode, cse IN odes
    cse$idx$ = $cse$
    out[$idx$] = $ode$
    # $JAFF END
    return out


def radodes():
    out = {}
    # $JAFF REPEAT idx, radode, cse IN radodes
    cse$idx$ = $cse$
    out[$idx$] = $radode$
    # $JAFF END
    return out


def dedt():
    # $JAFF SUB dedt
    return {0: $dedt$}
    # $JAFF END


def rhs_volumetric():
    out = {}
    # $JAFF REPEAT idx, rhs, cse IN rhses $[RADIATION True]$
    cse$idx$ = $cse$
    out[$idx$] = $rhs$
    # $JAFF END
    return out


def rhs_specific_mass():
    out = {}
    # $JAFF REPEAT idx, rhs, cse IN rhses $[RADIATION True DEDT_TYPE specific]$
    cse$idx$ = $cse$
    out[$idx$] = $rhs$
    # $JAFF END
    return out


def rhs_specific_number():
    out = {}
    # $JAFF REPEAT idx, rhs, cse IN rhses $[RADIATION True DEDT_TYPE per_particle]$
    cse$idx$ = $cse$
    out[$idx$] = $rhs$
    # $JAFF END
    return out


def jacobian_species():
    out = {}
    # $JAFF REPEAT idx, expr, cse IN jacobian $[RADIATION True]$
    cse$idx$ = $cse$
    out[($idx$, $idx$)] = $expr$
    # $JAFF END
    return out


def jacobian_volumetric():
    out = {}
    # $JAFF REPEAT idx, expr, cse IN jacobian $[RADIATION True THERMAL dedt]$
    cse$idx$ = $cse$
    out[($idx$, $idx$)] = $expr$
    # $JAFF END
    return out


def jacobian_specific_mass():
    out = {}
    # $JAFF REPEAT idx, expr, cse IN jacobian $[RADIATION True THERMAL dedt DEDT_TYPE specific]$
    cse$idx$ = $cse$
    out[($idx$, $idx$)] = $expr$
    # $JAFF END
    return out


def jacobian_specific_number():
    out = {}
    # $JAFF REPEAT idx, expr, cse IN jacobian $[RADIATION True THERMAL dedt DEDT_TYPE per_particle]$
    cse$idx$ = $cse$
    out[($idx$, $idx$)] = $expr$
    # $JAFF END
    return out
