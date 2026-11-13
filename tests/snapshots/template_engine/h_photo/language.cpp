// ABOUTME: Template-engine regression fixture: C++ smoke test: indented directives, brace arrays, CSE, REPLACE
// ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

// Plain C++ comment: passes through untouched.
#include <array>

// $$JAFF SUB nspec, nreact, label
constexpr int NSPEC = 3;
constexpr int NREACT = 2;
constexpr int NSPEC_PLUS_ONE = 4;
const char *LABEL = "h_photo";
// $$JAFF END

// $$JAFF GET specie_idx IN H
constexpr int IDX_H = 0;
// $$JAFF END

// $$JAFF HAS specie He+ { REPLACE { "1" -> "true"  "0" -> "false" } }
constexpr bool HAS_HEP = false;
// $$JAFF END

// $$JAFF REPEAT specie IN species
const std::array<const char *, NSPEC> names = {"H", "H+", "e-"};
// $$JAFF END

// $$JAFF REPEAT specie_mass IN specie_masses
const double masses[] = {1.673773e-24, 1.6728620616289998e-24, 9.10938371e-28};
// $$JAFF END

// $$JAFF REDUCE specie_charge IN specie_charges
constexpr double TOTAL_CHARGE = ((0) + (1) + (-1));
// $$JAFF END

void rates(double *k, const double *nden, double tgas) {
    // $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rates } }
    // $$JAFF END
    // $$JAFF REPEAT idx, rate IN rates { CSE { GROUP -> rates } }
    k[0] = 5.0572844545056314e-8*photden[0];
    k[1] = 1.6594178159829075e-10*std::pow(tgas, -0.69999999999999996);
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
    const amrex::Real cse0 = state.xn[1]*state.xn[2];
    const amrex::Real cse1 = 1.6594178159829075e-10*cse0*std::pow(state.T, -0.69999999999999996) - 5.0572844545056314e-8*state.xn[0]*photden[0];
    const amrex::Real cse2 = -cse1;
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
    ydot(2) = cse2;
    ydot(3) = cse2;
    ydot(4) = -2.2910735482189854e-26*cse0*std::pow(state.T, 0.30000000000000004)*(0.68400000000000005 - 0.041599999999999998*std::log(0.0001*state.T)) + 3.2366620508836041e-19*state.xn[0]*photden[0];
    // $$JAFF END
}

void jacobian(double *J, const double *nden, double tgas) {
    // $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> jacobian } }
    const double cse0 = 5.0572844545056314e-8*photden[0];
    const double cse1 = 1.6594178159829075e-10*std::pow(tgas, -0.69999999999999996);
    const double cse2 = cse1*nden[2];
    const double cse3 = cse1*nden[1];
    const double cse4 = -cse2;
    const double cse5 = -cse3;
    // $$JAFF END
    // $$JAFF REPEAT idx, jac IN jacobian { CSE { GROUP -> jacobian } }
    J(1, 1) = -cse0;
    J(1, 2) = cse2;
    J(1, 3) = cse3;
    J(2, 1) = cse0;
    J(2, 2) = cse4;
    J(2, 3) = cse5;
    J(3, 1) = cse0;
    J(3, 2) = cse4;
    J(3, 3) = cse5;
    // $$JAFF END
}
