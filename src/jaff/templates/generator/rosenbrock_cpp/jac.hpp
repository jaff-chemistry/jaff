#ifndef ROSENBROCK_JAC_H
#define ROSENBROCK_JAC_H

#include <cmath>
#include <vector>

#include "integrator_data.hpp"
#include "parameters.hpp"
#include "state.hpp"

template <typename StateT>
inline void jac(const StateT &state, DArray2D &jac) {
  $$JAFF REPEAT idx, expr, cse IN jacobian {
      REPLACE { "nden\[\s*(\d+)\s*\]" -> "state.rho[\1]"  "tgas" -> "state.get_T()" }
      THERMAL { MODE -> dedt }
  }
  const double cse$idx$ = $cse$;
  jac[$idx$][$idx$] = $expr$;
  $$JAFF END
}

// The RHS and the Jacobian at one state, sharing every temperature-only rate
// coefficient and its derivative. Writes only the Jacobian's non-zero entries.
template <typename StateT, typename VectorT>
inline void rhs_jac(const StateT &state, VectorT &ydot, DArray2D &jac) {
  const double T = state.get_T();
  $$JAFF REPEAT idx, rhs, cse IN rhs_jacobian {
      REPLACE { "nden\[\s*(\d+)\s*\]" -> "state.rho[\1]"  "tgas" -> "T" }
      THERMAL { MODE -> dedt }
  }
  const double cse$idx$ = $cse$;
  ydot[$idx$] = $rhs$;
  $$JAFF END
  $$JAFF REPEAT idx, expr IN rhs_jacobian_entries {
      REPLACE { "nden\[\s*(\d+)\s*\]" -> "state.rho[\1]"  "tgas" -> "T" }
      THERMAL { MODE -> dedt }
  }
  jac[$idx$][$idx$] = $expr$;
  $$JAFF END
}

#endif
