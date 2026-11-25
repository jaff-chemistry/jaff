// ABOUTME: Template-engine regression fixture: C++ smoke test: indented directives, brace arrays, CSE, REPLACE
// ABOUTME: Rendered output is pinned in tests/snapshots/template_engine/

// Plain C++ comment: passes through untouched.
#include <array>

// $$JAFF SUB nspec, nreact, label
constexpr int NSPEC = 5;
constexpr int NREACT = 3;
constexpr int NSPEC_PLUS_ONE = 6;
const char *LABEL = "react_cie_hepp";
// $$JAFF END

// $$JAFF GET specie_idx IN H
constexpr int IDX_H = 0;
// $$JAFF END

// $$JAFF HAS specie He+ { REPLACE { "1" -> "true"  "0" -> "false" } }
constexpr bool HAS_HEP = true;
// $$JAFF END

// $$JAFF REPEAT specie IN species
const std::array<const char *, NSPEC> names = {"H", "e-", "H+", "He++", "He+"};
// $$JAFF END

// $$JAFF REPEAT specie_mass IN specie_masses
const double masses[] = {1.673773e-24, 9.10938371e-28, 1.6728620616289998e-24, 6.644651123258e-24, 6.645562061629e-24};
// $$JAFF END

// $$JAFF REDUCE specie_charge IN specie_charges
constexpr double TOTAL_CHARGE = ((0) + (-1) + (1) + (2) + (1));
// $$JAFF END

void rates(double *k, const double *nden, double tgas) {
    // $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rates } }
    // $$JAFF END
    // $$JAFF REPEAT idx, rate IN rates { CSE { GROUP -> rates } }
    k[0] = 5.8500000000000005e-11*std::sqrt(tgas)*std::exp(-157809.10000000001/tgas);
    k[1] = 4.1207223003988966e-10*std::pow(tgas, -0.80000000000000004);
    k[2] = nden[3];
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
    const amrex::Real cse0 = 1.0/state.T;
    const amrex::Real cse3 = std::log(state.T);
    const amrex::Real cse4 = std::exp(0.19999999999999996*cse3);
    const amrex::Real cse5 = std::sqrt(state.T);
    const amrex::Real cse1 = 4.1207223003988966e-10*cse0*cse4*state.xn[1]*state.xn[2] - 5.8500000000000005e-11*cse5*std::exp(-157809.10000000001*cse0)*state.xn[0]*state.xn[1];
    const amrex::Real cse2 = state.xn[1]*(state.xn[3]*state.xn[3]);
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
    ydot(2) = -cse1 - cse2;
    ydot(3) = -cse1;
    ydot(4) = -cse2;
    ydot(5) = cse2;
    ydot(6) = 0;
    // $$JAFF END
}

void jacobian(double *J, const double *nden, double tgas) {
    // $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> jacobian } }
    const double cse0 = 1.0/tgas;
    const double cse12 = std::sqrt(tgas);
    const double cse1 = 5.8500000000000005e-11*cse12*std::exp(-157809.10000000001*cse0);
    const double cse2 = cse1*nden[1];
    const double cse13 = std::log(tgas);
    const double cse14 = std::exp(0.19999999999999996*cse13);
    const double cse3 = 4.1207223003988966e-10*cse0*cse14;
    const double cse4 = cse3*nden[2];
    const double cse5 = cse1*nden[0];
    const double cse6 = -cse4 + cse5;
    const double cse7 = cse3*nden[1];
    const double cse8 = nden[3]*nden[3];
    const double cse9 = -cse7;
    const double cse10 = 2*nden[1]*nden[3];
    const double cse11 = -cse10;
    // $$JAFF END
    // $$JAFF REPEAT idx, jac IN jacobian { CSE { GROUP -> jacobian } }
    J(1, 1) = -cse2;
    J(1, 2) = -cse6;
    J(1, 3) = cse7;
    J(2, 1) = cse2;
    J(2, 2) = -cse4 + cse5 - cse8;
    J(2, 3) = cse9;
    J(2, 4) = cse11;
    J(3, 1) = cse2;
    J(3, 2) = cse6;
    J(3, 3) = cse9;
    J(4, 2) = -cse8;
    J(4, 4) = cse11;
    J(5, 2) = cse8;
    J(5, 4) = cse10;
    // $$JAFF END
}
