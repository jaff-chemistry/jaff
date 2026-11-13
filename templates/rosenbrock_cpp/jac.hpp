#ifndef ROSENBROCK_JAC_H
#define ROSENBROCK_JAC_H

#include <cmath>
#include <vector>

#include "integrator_data.hpp"
#include "parameters.hpp"
#include "state.hpp"

// Network symbols mapped onto the integrator's state.
$$JAFF DEFINE state {
    REPLACE { "nden\[\s*(\d+)\s*\]" -> "state.rho[\1]"  "\btgas\b" -> "T" }
}

// The RHS and the Jacobian at one state, sharing every rate coefficient and
// every common sub-expression. Writes only the Jacobian's non-zero entries.
template <typename StateT, typename VectorT>
inline void rhs_jac(const StateT &state, VectorT &ydot, DArray2D &jac) {
  [[maybe_unused]] const double T = state.get_T();

  $$JAFF REPEAT idx, cse IN cses {
    USE state
    GROUP { NAME -> rhs_jac }
  }
  const double cse$idx$ = $cse$;
  $$JAFF END

  $$JAFF REPEAT idx, rhs IN rhses {
    USE state
    THERMAL { MODE -> dedt }
    CSE { GROUP -> rhs_jac }
  }
  ydot[$idx$] = $rhs$;
  $$JAFF END

  $$JAFF REPEAT idx, jac IN jacobian {
    USE state
    THERMAL { MODE -> dedt }
    CSE { GROUP -> rhs_jac }
  }
  jac[$idx$][$idx$] = $jac$;
  $$JAFF END
}

#endif
