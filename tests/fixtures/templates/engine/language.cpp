// ABOUTME: Template-engine regression fixture: C++ smoke test: indented directives, brace arrays, CSE, REPLACE
// ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

// Plain C++ comment: passes through untouched.
#include <array>

$$JAFF SUB nspec, nreact, label
constexpr int NSPEC = $nspec$;
constexpr int NREACT = $nreact$;
constexpr int NSPEC_PLUS_ONE = $nspec+1$;
const char *LABEL = "$label$";
$$JAFF END

$$JAFF GET specie_idx IN H
constexpr int IDX_H = $specie_idx$;
$$JAFF END

$$JAFF HAS specie He+ { REPLACE { "1" -> "true"  "0" -> "false" } }
constexpr bool HAS_HEP = $specie$;
$$JAFF END

$$JAFF REPEAT specie IN species
const std::array<const char *, NSPEC> names = {"$specie$", };
$$JAFF END

$$JAFF REPEAT specie_mass IN specie_masses
const double masses[] = {$specie_mass$, };
$$JAFF END

$$JAFF REDUCE specie_charge IN specie_charges
constexpr double TOTAL_CHARGE = $($specie_charge$)$;
$$JAFF END

void rates(double *k, const double *nden, double tgas) {
    $$JAFF REPEAT idx, rate, cse IN rates
    const double cse$idx$ = $cse$;
    k[$idx$] = $rate$;
    $$JAFF END
}

void rhs(double *f, const double *nden, double tgas) {
    $$JAFF REPEAT idx, rhs, cse IN rhses {
        REPLACE {
            "nden\[(\d+)\]" -> "state.xn[\1]"
            "tgas" -> "state.T"
        }
        THERMAL { MODE -> dedt }
    }
    const amrex::Real cse$idx$ = $cse$;
    ydot($idx+1$) = $rhs$;
    $$JAFF END
}

void jacobian(double *J, const double *nden, double tgas) {
    $$JAFF REPEAT idx, expr, cse IN jacobian
    const double cse$idx$ = $cse$;
    J($idx+1$, $idx+1$) = $expr$;
    $$JAFF END
}
