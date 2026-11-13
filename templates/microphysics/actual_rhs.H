#ifndef actual_rhs_H
#define actual_rhs_H

#include <AMReX_Array.H>
#include <AMReX_REAL.H>
#include <cmath>

#include <actual_network.H>
#include <burn_type.H>
#include <extern_parameters.H>
#include <jacobian_utilities.H>

using namespace ArrayUtil;
using namespace network_rp;

// Network symbols mapped onto burn_t, shared by the RHS and the Jacobian.
$$JAFF DEFINE microphysics {
    REPLACE {
        "radeden\[\s*(\d+)\s*\]" -> "state.rn[2*\1]"
        "photden\[\s*(\d+)\s*\]" -> "state.rn[2*\1]"
        "rflux\[\s*(\d+)\s*\]" -> "state.rn[2*\1+1]"
        "nden\[\s*(\d+)\s*\]" -> "state.xn[\1]"
        "tgas" -> "state.T"
        "c_hat" -> "state.c_hat"
    }
}

// The equations integrated, shared by the RHS and the Jacobian.
$$JAFF DEFINE microphysics_equations {
    THERMAL { MODE -> dedt  ENERGY -> specific }
    # RADIATION { ENABLED -> True }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE void
actual_rhs(const burn_t &state, amrex::Array1D<amrex::Real, 1, neqs> &ydot) {

  using namespace Rates;

  $$JAFF REPEAT idx, cse IN cses { USE microphysics  GROUP { NAME -> rhs } }
  const amrex::Real cse$idx$ = $cse$;
  $$JAFF END
  $$JAFF REPEAT idx, rhs IN rhses {
      USE microphysics  USE microphysics_equations  CSE { GROUP -> rhs }
  }
  ydot($idx+1$) = $rhs$;
  $$JAFF END
}

template <class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE void actual_jac(const burn_t &state,
                                                   MatrixType &jac) {
  jac.zero();

  $$JAFF REPEAT idx, cse IN cses { USE microphysics  GROUP { NAME -> jacobian } }
  const amrex::Real cse$idx$ = $cse$;
  $$JAFF END
  $$JAFF REPEAT idx, jac IN jacobian {
      USE microphysics  USE microphysics_equations  CSE { GROUP -> jacobian }
  }
  jac($idx+1$, $idx+1$) = $jac$;
  $$JAFF END
}

AMREX_INLINE
void actual_rhs_init() {}

#endif
