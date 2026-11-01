#ifndef ROSENBROCK_INTEGRATOR_H
#define ROSENBROCK_INTEGRATOR_H

#include <algorithm>
#include <cmath>
#include <vector>

#include "integrator_data.hpp"
#include "linpack.hpp"
#include "parameters.hpp"
#include "tableau.hpp"
#include "types.hpp"

namespace rosenbrock {

template <int int_neqs>
inline double &stage(rosenbrock_t<int_neqs> &rstate, const int i, const int n) {
  switch (i) {
  case 0:
    return rstate.ak1[n];
  case 1:
    return rstate.ak2[n];
  case 2:
    return rstate.ak3[n];
  case 3:
    return rstate.ak4[n];
  case 4:
    return rstate.ak5[n];
  case 5:
    return rstate.ak6[n];
  case 6:
    return rstate.ak7[n];
  case 7:
    return rstate.ak8[n];
  default:
    return rstate.work[n];
  }
}

template <int int_neqs>
inline double stage(const rosenbrock_t<int_neqs> &rstate, const int i,
                    const int n) {
  switch (i) {
  case 0:
    return rstate.ak1[n];
  case 1:
    return rstate.ak2[n];
  case 2:
    return rstate.ak3[n];
  case 3:
    return rstate.ak4[n];
  case 4:
    return rstate.ak5[n];
  case 5:
    return rstate.ak6[n];
  case 6:
    return rstate.ak7[n];
  case 7:
    return rstate.ak8[n];
  default:
    return rstate.work[n];
  }
}

template <int int_neqs>
inline DArray1D &stage_vector(rosenbrock_t<int_neqs> &rstate, const int i) {
  switch (i) {
  case 0:
    return rstate.ak1;
  case 1:
    return rstate.ak2;
  case 2:
    return rstate.ak3;
  case 3:
    return rstate.ak4;
  case 4:
    return rstate.ak5;
  case 5:
    return rstate.ak6;
  case 6:
    return rstate.ak7;
  case 7:
    return rstate.ak8;
  default:
    return rstate.work;
  }
}

template <typename RosenbrockT>
inline double rtol_for(const RosenbrockT &rstate, const int n) {
  return n == params::IE ? rstate.rtol_echem : rstate.rtol_spec;
}

template <typename RosenbrockT>
inline double atol_for(const RosenbrockT &rstate, const int n) {
  return n == params::IE ? rstate.atol_echem : rstate.atol_spec;
}

template <int int_neqs>
inline double error_norm(const rosenbrock_t<int_neqs> &rstate) {
  double err = 0.0;
  for (int n = 0; n < int_neqs; ++n) {
    const double sk = atol_for(rstate, n) +
                      rtol_for(rstate, n) * std::max(std::abs(rstate.y[n]),
                                                     std::abs(rstate.ynew[n]));
    const double term = rstate.work[n] / sk;
    err += term * term;
  }
  return std::sqrt(err / static_cast<double>(int_neqs));
}

template <int int_neqs>
inline double yass_change_norm(const rosenbrock_t<int_neqs> &rstate) {
  double max_change = 0.0;
  double max_species = 0.0;

  for (int i = 0; i < params::NumSpec && i < int_neqs; ++i) {
    max_species = std::max(max_species, std::abs(rstate.y[i]));
  }

  const double active_floor = integrator::yass_floor * max_species;

  for (int i = 0; i < params::NumSpec && i < int_neqs; ++i) {
    const double yi = std::abs(rstate.y[i]);
    if (yi > active_floor) {
      max_change =
          std::max(max_change, std::abs(rstate.ynew[i] - rstate.y[i]) / yi);
    }
  }

  return max_change / integrator::yass_epsilon;
}

template <typename MatrixType, int int_neqs>
inline bool matrix_is_finite(const MatrixType &matrix) {
  for (int j = 0; j < int_neqs; ++j) {
    for (int i = 0; i < int_neqs; ++i) {
      const double value = matrix[i][j];
      if (std::isnan(value) || std::isinf(value)) {
        return false;
      }
    }
  }

  return true;
}

template <int int_neqs> inline bool valid_integrator_state(const DArray1D &y) {
  for (int n = 0; n < int_neqs; ++n) {
    const double yn = y[n];
    if (std::isnan(yn) || std::isinf(yn)) {
      return false;
    }
  }

  return true;
}

// Copy the integration vector y into the thermodynamic state that rhs()
// and jac() read from.
template <typename BurnT, int int_neqs>
inline void vector_to_state(BurnT &state, const DArray1D &y) {
  for (int i = 0; i < params::NumSpec; ++i) {
    state.rho[i] = y[i];
  }
  state.e = y[params::NumSpec];
}

template <typename BurnT, int int_neqs>
inline void evaluate_jacobian(BurnT &state, rosenbrock_t<int_neqs> &rstate,
                              const double time) {

  vector_to_state<BurnT, int_neqs>(state, rstate.y);
  rhs(state, rstate.ak1);
  rstate.n_rhs += 1;

  if (rstate.jacobian_type == 1) {
    // jac() only writes the structurally non-zero entries, and decompose()
    // overwrites rstate.jac in place with the LU factors. Zero the matrix
    // first so stale factors from the previous step do not leak into the
    // entries jac() leaves untouched.
    for (int i = 0; i < int_neqs; ++i) {
      for (int j = 0; j < int_neqs; ++j) {
        rstate.jac[i][j] = 0.0;
      }
    }

    jac(state, rstate.jac);

    if (matrix_is_finite<decltype(rstate.jac), int_neqs>(rstate.jac)) {
      rstate.n_jac += 1;
      return;
    }
  }

  constexpr double UROUND = std::numeric_limits<double>::epsilon();
  DArray1D &ewt = rstate.jac_ewt;
  DArray1D &ybase = rstate.jac_ybase;

  // rhs() may clean or renormalize all components of rstate.y.
  for (int i = 0; i < int_neqs; ++i) {
    ybase[i] = rstate.y[i];
  }

  double fac = 0.0;
  for (int i = 0; i < int_neqs; ++i) {
    ewt[i] = 1.0 / (rtol_for(rstate, i) * std::abs(rstate.y[i]) +
                    atol_for(rstate, i));
    fac += (rstate.ak1[i] * ewt[i]) * (rstate.ak1[i] * ewt[i]);
  }
  fac = std::sqrt(fac / static_cast<double>(int_neqs));

  double R0 = 1000.0 * std::abs(rstate.dt) * UROUND *
              static_cast<double>(int_neqs) * fac;
  if (R0 == 0.0) {
    R0 = 1.0;
  }

  for (int j = 0; j < int_neqs; ++j) {
    const double yj = ybase[j];
    const double R = std::max(std::sqrt(UROUND) * std::abs(yj), R0 / ewt[j]);
    for (int n = 0; n < int_neqs; ++n) {
      rstate.y[n] = ybase[n];
    }
    rstate.y[j] += R;

    vector_to_state<BurnT, int_neqs>(state, rstate.y);
    rhs(state, rstate.work);

    fac = 1.0 / R;
    for (int i = 0; i < int_neqs; ++i) {
      rstate.jac[i][j] = (rstate.work[i] - rstate.ak1[i]) * fac;
    }

    for (int n = 0; n < int_neqs; ++n) {
      rstate.y[n] = ybase[n];
    }
  }

  // restore the thermodynamic state to the unperturbed solution
  vector_to_state<BurnT, int_neqs>(state, rstate.y);

  rstate.n_rhs += int_neqs;
  rstate.n_jac += 1;
}

template <typename BurnT, int int_neqs>
inline void rhs_at(const double time, BurnT &state,
                   rosenbrock_t<int_neqs> &rstate, const std::vector<double> &y,
                   std::vector<double> &ydot) {
  // Evaluate the RHS at the trial vector y by syncing it into the
  // thermodynamic state, then restore the state to the current solution.
  vector_to_state<BurnT, int_neqs>(state, y);
  rhs(state, ydot);
  vector_to_state<BurnT, int_neqs>(state, rstate.y);
}

template <int int_neqs>
inline int decompose(rosenbrock_t<int_neqs> &rstate, const double fac) {
  for (int j = 0; j < int_neqs; ++j) {
    for (int i = 0; i < int_neqs; ++i) {
      rstate.jac[i][j] = -rstate.jac[i][j];
    }
    rstate.jac[j][j] += fac;
  }

  int ierr_linpack = 0;
  if (integrator::linalg_do_pivoting == 1) {
    constexpr bool allow_pivot{true};
    dgefa<int_neqs, allow_pivot>(rstate.jac, rstate.pivot, ierr_linpack);
  } else {
    constexpr bool allow_pivot{false};
    dgefa<int_neqs, allow_pivot>(rstate.jac, rstate.pivot, ierr_linpack);
  }

  return ierr_linpack;
}

template <int int_neqs>
inline void solve(rosenbrock_t<int_neqs> &rstate, DArray1D &b) {
  if (integrator::linalg_do_pivoting == 1) {
    constexpr bool allow_pivot{true};
    dgesl<int_neqs, allow_pivot>(rstate.jac, rstate.pivot, b);
  } else {
    constexpr bool allow_pivot{false};
    dgesl<int_neqs, allow_pivot>(rstate.jac, rstate.pivot, b);
  }
}

} // namespace rosenbrock

template <typename BurnT, typename Tableau>
inline int rosenbrock_integrator(BurnT &state, rosenbrock_t<INT_NEQS> &rstate) {
  using C = Tableau;
  constexpr int int_neqs = INT_NEQS;
  static_assert(C::stages <= 8, "Rosenbrock integrator stores up to 8 stages");

  if (rstate.tout == rstate.t) {
    return IERR_SUCCESS;
  }

  for (int n = 0; n < int_neqs; ++n) {
    if (rosenbrock::atol_for(rstate, n) <= 0.0 ||
        rosenbrock::rtol_for(rstate, n) <=
            10.0 * std::numeric_limits<double>::epsilon()) {
      return IERR_TOO_MUCH_ACCURACY_REQUESTED;
    }
  }

  if (integrator::rosenbrock_timestep_controller == 1 &&
      (integrator::yass_epsilon <= 0.0 || integrator::yass_floor < 0.0 ||
       integrator::yass_fac_min <= 0.0 ||
       integrator::yass_fac_max < integrator::yass_fac_min ||
       integrator::yass_reduction_fac <= 0.0)) {
    return IERR_BAD_INPUTS;
  }

  rstate.n_rhs = 0;
  rstate.n_jac = 0;
  rstate.n_step = 0;

  // start a fresh trajectory record (initial state first)
  if (rstate.record_history) {
    rstate.hist_t.clear();
    rstate.hist_y.clear();
    rstate.hist_t.push_back(rstate.t);
    rstate.hist_y.push_back(rstate.y);
  }

  // honor a caller-supplied initial timestep guess; otherwise span the
  // whole integration interval
  if (rstate.dt <= 0.0) {
    rstate.dt = rstate.tout - rstate.t;
  }

  const double posneg = rstate.tout >= rstate.t ? 1.0 : -1.0;
  const double hmax =
      std::min(integrator::ode_max_dt, std::abs(rstate.tout - rstate.t));
  double h = std::min(std::abs(rstate.dt), hmax) * posneg;
  if (std::abs(h) <= 10.0 * std::numeric_limits<double>::epsilon()) {
    h = 1.e-6 * posneg;
  }

  bool reject = false;
  bool last = false;
  int nsing = 0;
  int n_reject = 0;
  double facold = 1.0;
  double errold = 1.0;
  double hopt = h;
  double x = rstate.t;
  std::vector<double> rhs_tmp(int_neqs);

  // main evolution loop

  // within this loop, x will be the current time and h will be the
  // step size we are attempting.  The solution at the current time
  // is rstate.y().

  while (true) {

    if (rstate.n_step > integrator::ode_max_steps) {
      rstate.t = x;
      rstate.dt = h;
      return IERR_TOO_MANY_STEPS;
    }

    if (0.1 * std::abs(h) <=
        std::abs(x) * std::numeric_limits<double>::epsilon()) {
      rstate.t = x;
      rstate.dt = h;
      return IERR_DT_UNDERFLOW;
    }

    if (last) {
      rstate.t = x;
      rstate.dt = hopt;
      return IERR_SUCCESS;
    }

    hopt = h;
    if ((x + h * (1.0 + timestep_safety_factor) - rstate.tout) * posneg >=
        0.0) {
      h = rstate.tout - x;
      last = true;
    }

    // take a step -- this loop will keep trying until a step is
    // successful or a catastrophic error is encountered.

    while (true) {
      rstate.dt = h;

      if (0.1 * std::abs(h) <=
          std::abs(x) * std::numeric_limits<double>::epsilon()) {
        rstate.t = x;
        rstate.dt = h;
        return IERR_DT_UNDERFLOW;
      }

      const double fac = 1.0 / (h * C::gamma);
      rosenbrock::evaluate_jacobian(state, rstate, x);
      int ierr_linpack = rosenbrock::decompose(rstate, fac);

      if (ierr_linpack != 0) {
        nsing += 1;
        if (nsing >= 5) {
          rstate.t = x;
          rstate.dt = h;
          return IERR_LU_DECOMPOSITION_ERROR;
        }
        h *= 0.5;
        reject = true;
        last = false;
        continue;
      }

      rosenbrock::solve(rstate, rstate.ak1);

      bool valid_stage_state = true;
      for (int istage = 1; istage < C::stages; ++istage) {
        for (int n = 0; n < int_neqs; ++n) {
          rstate.ynew[n] = rstate.y[n];
          rosenbrock::stage(rstate, istage, n) = 0.0;

          for (int jstage = 0; jstage < istage; ++jstage) {
            const double akj = rosenbrock::stage(rstate, jstage, n);
            rstate.ynew[n] += C::a(istage + 1, jstage + 1) * akj;
            rosenbrock::stage(rstate, istage, n) +=
                (C::c(istage + 1, jstage + 1) / h) * akj;
          }
        }

        if (!rosenbrock::valid_integrator_state<int_neqs>(rstate.ynew)) {
          h *= 0.25;
          reject = true;
          n_reject += 1;
          last = false;
          valid_stage_state = false;
          break;
        }

        rosenbrock::rhs_at(x + C::ct(istage + 1) * h, state, rstate, rstate.ynew,
                           rhs_tmp);
        rstate.n_rhs += 1;

        for (int n = 0; n < int_neqs; ++n) {
          rosenbrock::stage(rstate, istage, n) += rhs_tmp[n];
        }
        rosenbrock::solve(rstate, rosenbrock::stage_vector(rstate, istage));
      }

      if (!valid_stage_state) {
        continue;
      }

      for (int n = 0; n < int_neqs; ++n) {
        double y_update = 0.0;
        double err = 0.0;

        for (int istage = 0; istage < C::stages; ++istage) {
          const double aki = rosenbrock::stage(rstate, istage, n);
          y_update += C::b(istage + 1) * aki;
          err += C::e(istage + 1) * aki;
        }

        rstate.ynew[n] = rstate.y[n] + y_update;
        rstate.work[n] = err;
      }

      rstate.n_step += 1;

      const bool valid_state =
          rosenbrock::valid_integrator_state<int_neqs>(rstate.ynew);
      if (!valid_state) {
        reject = true;
        n_reject += 1;
        last = false;
        h *= 0.25;
        continue;
      }

      constexpr double err_min = 1.e-10;
      double err;
      double controller_fac = 1.0;

      if (integrator::rosenbrock_timestep_controller == 1) {
        err = rosenbrock::yass_change_norm(rstate);
        controller_fac =
            std::clamp(1.0 / std::max(err, err_min), integrator::yass_fac_min,
                       integrator::yass_fac_max);
      } else {
        err = rosenbrock::error_norm(rstate);
        controller_fac = std::clamp(
            std::pow(1.0 / std::max(err, err_min),
                     1.0 / (integrator::h211b_b * integrator::h211b_k)) *
                std::pow(1.0 / std::max(errold, err_min),
                         1.0 / (integrator::h211b_b * integrator::h211b_k)) *
                std::pow(facold, -1.0 / integrator::h211b_b),
            integrator::h211b_fac_min, integrator::h211b_fac_max);
      }

      double hnew = h * controller_fac;

      // update the starting time and state for the next step
      // and predict the next timestep

      if (err <= 1.0) {
        // step accepted: advance the digital-filter history
        facold = controller_fac;
        errold = err;

        for (int n = 0; n < int_neqs; ++n) {
          rstate.y[n] = rstate.ynew[n];
        }
        x += h;

        // record the accepted step
        if (rstate.record_history) {
          rstate.hist_t.push_back(x);
          rstate.hist_y.push_back(rstate.y);
        }

        if (std::abs(hnew) > hmax) {
          hnew = posneg * hmax;
        }
        if (reject) {
          hnew = posneg * std::min(std::abs(hnew), std::abs(h));
        }
        reject = false;
        n_reject = 0;
        h = hnew;
        break;
      }

      // if we made it here, then we've failed
      // cut the timestep in prep for trying again

      reject = true;
      n_reject += 1;
      last = false;
      if (n_reject >= 2) {
        hnew *= integrator::rosenbrock_timestep_controller == 1
                    ? integrator::yass_reduction_fac
                    : integrator::h211b_reduction_fac;
      }
      const double reject_fac = integrator::rosenbrock_timestep_controller == 1
                                    ? integrator::yass_reduction_fac
                                    : integrator::h211b_reduction_fac;
      h = posneg * std::min(std::abs(hnew), reject_fac * std::abs(h));
    }
  }
}

template <typename BurnT>
inline int rosenbrock_integrator(BurnT &state, rosenbrock_t<INT_NEQS> &rstate) {
  if (integrator::rosenbrock_tableau == 1) {
    return rosenbrock_integrator<BurnT, rosenbrock::rodas4p_tableau>(state,
                                                                     rstate);
  }
  if (integrator::rosenbrock_tableau == 2) {
    return rosenbrock_integrator<BurnT, rosenbrock::rodas3p_tableau>(state,
                                                                     rstate);
  }
  if (integrator::rosenbrock_tableau == 3) {
    return rosenbrock_integrator<BurnT, rosenbrock::ros2s_tableau>(state,
                                                                   rstate);
  }
  if (integrator::rosenbrock_tableau == 4) {
    return rosenbrock_integrator<BurnT, rosenbrock::ros2_tableau>(state,
                                                                  rstate);
  }
  if (integrator::rosenbrock_tableau == 5) {
    return rosenbrock_integrator<BurnT, rosenbrock::rosenbrock_euler_tableau>(
        state, rstate);
  }

  return rosenbrock_integrator<BurnT, rosenbrock::rodas5p_tableau>(state,
                                                                   rstate);
}

#endif
