// ABOUTME: Template-engine regression fixture: C++ smoke test: indented directives, brace arrays, CSE, REPLACE
// ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

// Plain C++ comment: passes through untouched.
#include <array>

// $$JAFF SUB nspec, nreact, label
constexpr int NSPEC = 3;
constexpr int NREACT = 1;
constexpr int NSPEC_PLUS_ONE = 4;
const char *LABEL = "test_jac_dedt";
// $$JAFF END

// $$JAFF GET specie_idx IN H
constexpr int IDX_H = 0;
// $$JAFF END

// $$JAFF HAS specie He+ { REPLACE { "1" -> "true"  "0" -> "false" } }
constexpr bool HAS_HEP = false;
// $$JAFF END

// $$JAFF REPEAT specie IN species
const std::array<const char *, NSPEC> names = {"H", "D", "He"};
// $$JAFF END

// $$JAFF REPEAT specie_mass IN specie_masses
const double masses[] = {1.673773e-24, 3.343586e-24, 6.646473e-24};
// $$JAFF END

// $$JAFF REDUCE specie_charge IN specie_charges
constexpr double TOTAL_CHARGE = ((0) + (0) + (0));
// $$JAFF END

void rates(double *k, const double *nden, double tgas) {
    // $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rates } }
    // $$JAFF END
    // $$JAFF REPEAT idx, rate IN rates { CSE { GROUP -> rates } }
    k[0] = nden[0];
    // $$JAFF END
}

void rhs(double *f, const double *nden, double tgas) {
    // $$JAFF REPEAT idx, cse IN cses {
        // REPLACE {
            // "nden\[(\d+)\]" -> "state.xn[\1]"
            // "tgas" -> "state.T"
        // }
        // GROUP { NAME -> rhs }
    // }
    const amrex::Real cse0 = state.xn[1]*(state.xn[0]*state.xn[0]);
    const amrex::Real cse1 = -cse0;
    // $$JAFF END
    // $$JAFF REPEAT idx, rhs IN rhses {
        // REPLACE {
            // "nden\[(\d+)\]" -> "state.xn[\1]"
            // "tgas" -> "state.T"
        // }
        // THERMAL { MODE -> dedt }
        // CSE { GROUP -> rhs }
    // }
    ydot(1) = cse1;
    ydot(2) = cse1;
    ydot(3) = cse0;
    ydot(4) = state.xn[1]*((state.xn[0]*state.xn[0])*state.xn[0]);
    // $$JAFF END
}

void jacobian(double *J, const double *nden, double tgas) {
    // $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> jacobian } }
    const double cse0 = 2*nden[0]*nden[1];
    const double cse1 = -cse0;
    const double cse2 = nden[0]*nden[0];
    const double cse3 = -cse2;
    // $$JAFF END
    // $$JAFF REPEAT idx, jac IN jacobian { CSE { GROUP -> jacobian } }
    J(1, 1) = cse1;
    J(1, 2) = cse3;
    J(2, 1) = cse1;
    J(2, 2) = cse3;
    J(3, 1) = cse0;
    J(3, 2) = cse2;
    // $$JAFF END
}
