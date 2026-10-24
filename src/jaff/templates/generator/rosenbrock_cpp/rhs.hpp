#ifndef ROSENBROCK_RHS_H
#define ROSENBROCK_RHS_H

#include <cmath>
#include <vector>

#include "integrator_data.hpp"
#include "parameters.hpp"
#include "state.hpp"

template <typename StateT>
inline void rhs(const StateT &state, std::vector<double> &ydot) {
  $$JAFF REPEAT idx, rhs, cse IN rhses {
      REPLACE { "nden\[\s*(\d+)\s*\]" -> "state.rho[\1]"  "tgas" -> "state.get_T()" }
      THERMAL { MODE -> dedt }
  }
  const double cse$idx$ = $cse$;
  ydot[$idx$] = $rhs$;
  $$JAFF END
}

inline void rhs_init() {}

#endif
