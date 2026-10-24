# ABOUTME: Template-engine regression fixture: radodes and RADIATION options (radiation network)
# ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

$$JAFF REPEAT idx, radode IN radodes
drad[$idx$] = $radode$
$$JAFF END

$$JAFF REPEAT idx, radode, cse IN radodes
r$idx$ = $cse$
drad[$idx$] = $radode$
$$JAFF END

$$JAFF REPEAT idx, rhs IN rhses { RADIATION { ENABLED -> True }  THERMAL { MODE -> dedt } }
f_rad[$idx$] = $rhs$
$$JAFF END

$$JAFF REPEAT idx, rhs IN rhses { RADIATION { ENABLED -> FALSE }  THERMAL { MODE -> dedt } }
f_norad[$idx$] = $rhs$
$$JAFF END

$$JAFF REPEAT idx, rhs IN rhses { RADIATION { ENABLED -> true  ORDER -> 1 }  THERMAL { MODE -> dedt } }
f_rad1[$idx$] = $rhs$
$$JAFF END

$$JAFF REPEAT idx, expr, cse IN jacobian {
    REPLACE {
        "photden\[\s*(\d+)\s*\]" -> "rn[2*\1]"
        "rflux\[\s*(\d+)\s*\]" -> "rn[2*\1+1]"
        "nden\[\s*(\d+)\s*\]" -> "xn[\1]"
        "tgas" -> "T"
        "c_hat" -> "chat"
    }
    RADIATION { ENABLED -> True }
}
c$idx$ = $cse$
J[$idx$, $idx$] = $expr$
$$JAFF END
