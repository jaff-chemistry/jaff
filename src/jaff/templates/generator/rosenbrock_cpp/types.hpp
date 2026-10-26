#ifndef ROSENBROCK_TYPE_H
#define ROSENBROCK_TYPE_H

#include "integrator_data.hpp"

#include <vector>

// When checking the integration time to see if we're done,
// be careful with roundoff issues.

const double timestep_safety_factor = 1.0e-4;

template <int int_neqs> struct rosenbrock_t {
  static_assert(int_neqs == INT_NEQS, "storage is sized for INT_NEQS");

  double t;    // the starting time
  double tout; // the stopping time
  double dt;   // next internal timestep

  int n_step;
  int n_rhs;
  int n_jac;

  double atol_spec;
  double rtol_spec;

  double atol_echem;
  double rtol_echem;

  DArray1D y{};
  DArray1D ynew{};
  // stage vectors k_1 .. k_8, indexed from 0
  std::array<DArray1D, 8> ak{};
  DArray1D work{};

  // scratch buffers for the numerical Jacobian (avoid per-call allocation)
  DArray1D jac_ewt{};
  DArray1D jac_ybase{};

  DArray2D jac{};
  IArray1D pivot{};

  short jacobian_type;

  // optional trajectory record: one entry per accepted step (plus the
  // initial state). Filled by rosenbrock_integrator; consumed by the caller.
  bool record_history = false;
  std::vector<double> hist_t;
  std::vector<DArray1D> hist_y;
};

#endif
