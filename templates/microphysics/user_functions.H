// ABOUTME: Functions the network calls but does not define (table interpolators, ...)
// ABOUTME: and their partial derivatives: replace each stub with an implementation
//
// actual_rhs.H includes this header. Each function below is a stub that aborts
// when called: replace the body of every one your run reaches. The Jacobian calls
// <name>_partial_<i>, the derivative of <name> with respect to its argument i
// (counted from 0). Stubs that are never called can stay as they are.
#ifndef JAFF_USER_FUNCTIONS_H
#define JAFF_USER_FUNCTIONS_H

#include <AMReX.H>
#include <AMReX_REAL.H>

$$JAFF REPEAT idx, function IN functions { REPLACE { "\b(a\d+)\b" -> "amrex::Real \1" } }
AMREX_GPU_HOST_DEVICE AMREX_INLINE amrex::Real $function$ { amrex::Abort("define $function$ in user_functions.H"); return 0.0; }
$$JAFF END

#endif
