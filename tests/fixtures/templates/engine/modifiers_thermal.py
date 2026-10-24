# ABOUTME: Template-engine regression fixture: THERMAL options and the dedt/dtdt SUB tokens (energy network)
# ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

$$JAFF SUB dedt, dtdt
dedt = $dedt$
dtdt = $dtdt$
$$JAFF END

$$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> none } }
f_none[$idx$] = $rhs$
$$JAFF END

$$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt } }
f_dedt[$idx$] = $rhs$
$$JAFF END

$$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dtdt } }
f_dtdt[$idx$] = $rhs$
$$JAFF END

$$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt  ENERGY -> specific } }
f_specific[$idx$] = $rhs$
$$JAFF END

$$JAFF REPEAT idx, rhs IN rhses { THERMAL { MODE -> dedt  ENERGY -> per_particle } }
f_per_particle[$idx$] = $rhs$
$$JAFF END

$$JAFF REPEAT idx, expr, cse IN jacobian {
    REPLACE { "nden\[(\d+)\]" -> "y[\1]" }
    THERMAL { MODE -> dtdt  ENERGY -> specific }
}
c$idx$ = $cse$
J[$idx$, $idx$] = $expr$
$$JAFF END

$$JAFF REPEAT idx, expr IN jacobian { THERMAL { MODE -> none } }
J_none[$idx$, $idx$] = $expr$
$$JAFF END
