// ABOUTME: Functions the network calls but does not define (table interpolators, ...)
// ABOUTME: and their partial derivatives: replace each stub with an implementation
//
// chemistry.hpp includes this header. Each function below is a stub that aborts
// the program when called: replace the body of every one your run reaches. The
// Jacobian calls <name>_partial_<i>, the derivative of <name> with respect to its
// argument i (counted from 0). Stubs that are never called can stay as they are.
#ifndef JAFF_USER_FUNCTIONS_HPP
#define JAFF_USER_FUNCTIONS_HPP

#include <cstdio>
#include <cstdlib>

$$JAFF REPEAT idx, function IN functions { REPLACE { "\b(a\d+)\b" -> "double \1" } }
inline double $function$ { std::fprintf(stderr, "define $function$ in user_functions.hpp\n"); std::abort(); }
$$JAFF END

#endif
