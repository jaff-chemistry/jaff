# ABOUTME: Template-engine regression fixture: REPEAT over codegen expressions (rates to Jacobian) with CSE
# ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

def rates():
    k = {}
    $$JAFF REPEAT idx, rate IN rates
    k[$idx$] = $rate$
    $$JAFF END
    return k


def rates_cse():
    k = {}
    $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rates } }
    cse$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, rate IN rates { CSE { GROUP -> rates } }
    k[$idx$] = $rate$
    $$JAFF END
    return k


def rates_inline():
    $$JAFF REPEAT rate IN rates
    return [$rate$, ]
    $$JAFF END


def fluxes():
    $$JAFF REPEAT idx, flux_expression IN flux_expressions
    flux[$idx$] = $flux_expression$
    $$JAFF END


def ode_expressions():
    $$JAFF REPEAT idx, ode_expression IN ode_expressions
    dy[$idx$] = $ode_expression$
    $$JAFF END


def odes():
    $$JAFF REPEAT idx, ode IN odes
    dy[$idx$] = $ode$
    $$JAFF END


def odes_cse():
    $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> odes } }
    tmp$idx$_v = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, ode IN odes { CSE { GROUP -> odes } }
    dy[$idx$] = $ode$
    $$JAFF END


def rhs():
    $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt } }
    f[$idx$] = $rhs$
    $$JAFF END


def rhs_cse_replace():
    $$JAFF REPEAT idx, cse IN cses {
        REPLACE {
            "nden\[(\d+)\]" -> "y[\1]"
            "tgas" -> "T"
        }
        GROUP { NAME -> rhs }
    }
    x$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, rhs IN rhses {
        REPLACE {
            "nden\[(\d+)\]" -> "y[\1]"
            "tgas" -> "T"
        }
        THERMAL { MODE -> dedt }
        CSE { GROUP -> rhs }
    }
    f[$idx+1$] = $rhs$
    $$JAFF END


def jacobian():
    $$JAFF REPEAT idx, jac IN jacobian
    J[$idx$][$idx$] = $jac$
    $$JAFF END


def jacobian_cse():
    $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> jacobian } }
    c$idx$ = $cse$
    $$JAFF END
    $$JAFF REPEAT idx, jac IN jacobian { CSE { GROUP -> jacobian } }
    J[$idx$, $idx$] = $jac$
    J1[$idx+1$, $idx+1$] = $jac$
    $$JAFF END
