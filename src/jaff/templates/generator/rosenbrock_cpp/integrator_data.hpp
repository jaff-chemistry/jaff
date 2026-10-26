#ifndef INTEGRATOR_DATA_H
#define INTEGRATOR_DATA_H

#include "parameters.hpp"
#include <array>
#include <cstdint>
#include <limits>

// Define the size of the ODE system that will be integrated
// INT_NEQS = NumSpec (species) + 1 (energy)

constexpr int INT_NEQS = params::NumSpec + 1;

// for outputting the burn state failure, use enough precision
// to be reproducible (not guaranteed for std::fixed)
constexpr int OUTDIGITS = std::numeric_limits<double>::max_digits10;

enum integrator_errors : std::int8_t {
  IERR_SUCCESS = 1,
  IERR_BAD_INPUTS = -1,
  IERR_DT_UNDERFLOW = -2,
  IERR_SPRAD_CONVERGENCE = -3,
  IERR_TOO_MANY_STEPS = -4,
  IERR_TOO_MUCH_ACCURACY_REQUESTED = -5,
  IERR_CORRECTOR_CONVERGENCE = -6,
  IERR_LU_DECOMPOSITION_ERROR = -7,
  IERR_BAD_STATE_IN_CORRECTOR = -8,
  IERR_ENTERED_NSE = -100
};

// Fixed-size storage, so the integrator state, the Jacobian and its LU
// factors need no heap and can live in registers or private memory when
// the integrator runs inside a GPU kernel.
using DArray1D = std::array<double, INT_NEQS>;
using IArray1D = std::array<int, INT_NEQS>;
using DArray2D = std::array<std::array<double, INT_NEQS>, INT_NEQS>;

#endif
