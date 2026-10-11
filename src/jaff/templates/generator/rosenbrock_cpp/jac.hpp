#ifndef ROSENBROCK_JAC_H
#define ROSENBROCK_JAC_H

#include <cmath>
#include <vector>

#include "integrator_data.hpp"
#include "parameters.hpp"
#include "state.hpp"

template <typename StateT>
inline void jac(const StateT &state, DArray2D &jac) {
  // $JAFF REPEAT idx, expr, cse IN jacobian $[REPLACE nden\[\s*(\d+)\s*\] state.rho[\1] REPLACE tgas state.get_T() THERMAL dedt]$
  const double cse$idx$ = $cse$;
  jac[$idx$][$idx$] = $expr$;
  // $JAFF END
}

#endif
