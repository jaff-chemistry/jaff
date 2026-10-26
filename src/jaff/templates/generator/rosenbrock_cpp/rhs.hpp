#ifndef ROSENBROCK_RHS_H
#define ROSENBROCK_RHS_H

#include <cmath>

#include "integrator_data.hpp"
#include "parameters.hpp"
#include "state.hpp"

template <typename StateT, typename VectorT>
inline void rhs(const StateT &state, VectorT &ydot) {
  $$JAFF REPEAT idx, rhs, cse IN rhses {
      REPLACE { "nden\[\s*(\d+)\s*\]" -> "state.rho[\1]"  "tgas" -> "state.get_T()" }
  }
  const double cse$idx$ = $cse$;
  ydot[$idx$] = $rhs$;
  $$JAFF END
}

inline void rhs_init() {}

#endif
