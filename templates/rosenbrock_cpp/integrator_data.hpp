#ifndef INTEGRATOR_DATA_H
#define INTEGRATOR_DATA_H

#include "parameters.hpp"
#include <cstdint>
#include <limits>
#include <vector>

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

using DArray1D = std::vector<double>;
using SArray1D = std::vector<short>;
using IArray1D = std::vector<int>;
using DArray2D = std::vector<std::vector<double>>;

#endif
