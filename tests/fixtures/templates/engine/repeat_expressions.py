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
    $$JAFF REPEAT idx, rate, cse IN rates
    cse$idx$ = $cse$
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
    $$JAFF REPEAT idx, ode, cse IN odes
    tmp$idx$_v = $cse$
    dy[$idx$] = $ode$
    $$JAFF END


def rhs():
    $$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt } }
    f[$idx$] = $rhs$
    $$JAFF END


def rhs_cse_replace():
    $$JAFF REPEAT idx, rhs, cse IN rhses {
        REPLACE {
            "nden\[(\d+)\]" -> "y[\1]"
            "tgas" -> "T"
        }
        THERMAL { MODE -> dedt }
    }
    x$idx$ = $cse$
    f[$idx+1$] = $rhs$
    $$JAFF END


def jacobian():
    $$JAFF REPEAT idx, expr IN jacobian
    J[$idx$][$idx$] = $expr$
    $$JAFF END


def jacobian_cse():
    $$JAFF REPEAT idx, expr, cse IN jacobian
    c$idx$ = $cse$
    J[$idx$, $idx$] = $expr$
    J1[$idx+1$, $idx+1$] = $expr$
    $$JAFF END
