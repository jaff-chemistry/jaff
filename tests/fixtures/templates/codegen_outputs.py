# ABOUTME: jaffgen test template exposing every expression generator as a Python function
# ABOUTME: Rendered per network, then evaluated against tests/golden with shared inputs
import math

$$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign {
    REPLACE { "idx_ek" -> "idx_e" }
    SIGN { POS -> j  NEG -> k }
}
idx_$specie_with_normalized_sign$ = $idx$
$$JAFF END

SPECIE_MASSES = {}
$$JAFF REPEAT idx, specie_mass IN specie_masses
SPECIE_MASSES[$idx$] = $specie_mass$
$$JAFF END


def rates():
    k = {}
    $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rates } }
    x$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, rate IN rates { CSE { GROUP -> rates } }
    k[$idx$] = $rate$
    $$JAFF END
    return k


def fluxes():
    k = rates()
    y = nden
    out = {}
    $$JAFF REPEAT idx, flux_expression IN flux_expressions
    out[$idx$] = $flux_expression$
    $$JAFF END
    return out


def ode_expressions():
    flux = fluxes()
    out = {}
    $$JAFF REPEAT idx, ode_expression IN ode_expressions
    out[$idx$] = $ode_expression$
    $$JAFF END
    return out


def odes():
    out = {}
    $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> odes } }
    cse$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, ode IN odes { CSE { GROUP -> odes } }
    out[$idx$] = $ode$
    $$JAFF END
    return out


def radodes():
    out = {}
    $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> radodes } }
    cse$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, radode IN radodes { CSE { GROUP -> radodes } }
    out[$idx$] = $radode$
    $$JAFF END
    return out


def dedt():
    $$JAFF SUB dedt
    return {0: $dedt$}
    $$JAFF END


def rhs_volumetric():
    out = {}
    $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rhs } }
    cse$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, rhs IN rhses { RADIATION { ENABLED -> True }  THERMAL { MODE -> dedt }  CSE { GROUP -> rhs } }
    out[$idx$] = $rhs$
    $$JAFF END
    return out


def rhs_specific_mass():
    out = {}
    $$JAFF REPEAT idx, cse IN cses {
        GROUP { NAME -> rhs_1 }
    }
    cse$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, rhs IN rhses {
        RADIATION { ENABLED -> True }
        THERMAL { MODE -> dedt  ENERGY -> specific }
        CSE { GROUP -> rhs_1 }
    }
    out[$idx$] = $rhs$
    $$JAFF END
    return out


def rhs_specific_number():
    out = {}
    $$JAFF REPEAT idx, cse IN cses {
        GROUP { NAME -> rhs_2 }
    }
    cse$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, rhs IN rhses {
        RADIATION { ENABLED -> True }
        THERMAL { MODE -> dedt  ENERGY -> per_particle }
        CSE { GROUP -> rhs_2 }
    }
    out[$idx$] = $rhs$
    $$JAFF END
    return out


def jacobian_species():
    out = {}
    $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> jacobian } }
    cse$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, jac IN jacobian { RADIATION { ENABLED -> True }  CSE { GROUP -> jacobian } }
    out[($idx$, $idx$)] = $jac$
    $$JAFF END
    return out


def jacobian_volumetric():
    out = {}
    $$JAFF REPEAT idx, cse IN cses {
        GROUP { NAME -> jacobian_1 }
    }
    cse$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, jac IN jacobian {
        RADIATION { ENABLED -> True }
        THERMAL { MODE -> dedt }
        CSE { GROUP -> jacobian_1 }
    }
    out[($idx$, $idx$)] = $jac$
    $$JAFF END
    return out


def jacobian_specific_mass():
    out = {}
    $$JAFF REPEAT idx, cse IN cses {
        GROUP { NAME -> jacobian_2 }
    }
    cse$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, jac IN jacobian {
        RADIATION { ENABLED -> True }
        THERMAL { MODE -> dedt  ENERGY -> specific }
        CSE { GROUP -> jacobian_2 }
    }
    out[($idx$, $idx$)] = $jac$
    $$JAFF END
    return out


def jacobian_specific_number():
    out = {}
    $$JAFF REPEAT idx, cse IN cses {
        GROUP { NAME -> jacobian_3 }
    }
    cse$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, jac IN jacobian {
        RADIATION { ENABLED -> True }
        THERMAL { MODE -> dedt  ENERGY -> per_particle }
        CSE { GROUP -> jacobian_3 }
    }
    out[($idx$, $idx$)] = $jac$
    $$JAFF END
    return out
