#ifndef ROSENBROCK_TYPE_H
#define ROSENBROCK_TYPE_H

#include "integrator_data.hpp"

// When checking the integration time to see if we're done,
// be careful with roundoff issues.

const double timestep_safety_factor = 1.0e-4;

template <int int_neqs> struct rosenbrock_t {

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

  DArray1D y = DArray1D(int_neqs);
  DArray1D ynew = DArray1D(int_neqs);
  DArray1D ak1 = DArray1D(int_neqs);
  DArray1D ak2 = DArray1D(int_neqs);
  DArray1D ak3 = DArray1D(int_neqs);
  DArray1D ak4 = DArray1D(int_neqs);
  DArray1D ak5 = DArray1D(int_neqs);
  DArray1D ak6 = DArray1D(int_neqs);
  DArray1D ak7 = DArray1D(int_neqs);
  DArray1D ak8 = DArray1D(int_neqs);
  DArray1D work = DArray1D(int_neqs);

  // scratch buffers for the numerical Jacobian (avoid per-call allocation)
  DArray1D jac_ewt = DArray1D(int_neqs);
  DArray1D jac_ybase = DArray1D(int_neqs);

  DArray2D jac = DArray2D(int_neqs, DArray1D(int_neqs));
  IArray1D pivot = IArray1D(int_neqs);

  short jacobian_type;

  // optional trajectory record: one entry per accepted step (plus the
  // initial state). Filled by rosenbrock_integrator; consumed by the caller.
  bool record_history = false;
  std::vector<double> hist_t;
  std::vector<DArray1D> hist_y;
};

#endif
