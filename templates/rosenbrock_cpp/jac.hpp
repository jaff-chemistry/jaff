#ifndef ROSENBROCK_JAC_H
#define ROSENBROCK_JAC_H

#include <cmath>
#include <vector>

#include "integrator_data.hpp"
#include "parameters.hpp"
#include "state.hpp"

template <typename StateT>
inline void jac(const StateT &state, DArray2D &jac) {
  $$JAFF REPEAT idx, cse IN cses {
      REPLACE { "nden\[\s*(\d+)\s*\]" -> "state.rho[\1]"  "tgas" -> "state.get_T()" }
      GROUP { NAME -> jacobian }
  }
  const double cse$idx$ = $cse$;
  $$JAFF END
  $$JAFF REPEAT idx, jac IN jacobian {
      REPLACE { "nden\[\s*(\d+)\s*\]" -> "state.rho[\1]"  "tgas" -> "state.get_T()" }
      THERMAL { MODE -> dedt }
      CSE { GROUP -> jacobian }
  }
  jac[$idx$][$idx$] = $jac$;
  $$JAFF END
}

#endif
