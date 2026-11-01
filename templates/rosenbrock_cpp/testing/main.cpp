/* This file is for testing purposes and is not meant for production code */

#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

#include "integrator_data.hpp"
#include "parameters.hpp"
#include "types.hpp"
#include "rosenbrock.hpp"
#include "jac.hpp"
#include "rhs.hpp"
#include "state.hpp"
#include "eos.hpp"

int main() {
  constexpr int neqs = INT_NEQS;

  state<double> burn_state;

  // Initialize species densities from the per-species initial conditions
  // defined in parameters.hpp
  for (int i = 0; i < params::NumSpec; ++i) {
    burn_state.rho[i] = params::n_init[i];
  }

  // Temperature and internal energy (energy density needs the total number density
  double T_init = params::T_init;
  burn_state.e = eos::energy_from_temperature(T_init, burn_state.n_tot());

  // Setup integrator state
  rosenbrock_t<neqs> rstate;

  // Initial conditions: species densities
  for (int i = 0; i < params::NumSpec; ++i) {
    rstate.y[i] = burn_state.rho[i];
  }

  // Energy equation
  rstate.y[params::NumSpec] = burn_state.e;

  rstate.t = 0.0;
  rstate.tout = 1.0e-6;   // Integrate to 1 microsecond
  rstate.dt = 1.0e-12;

  // Tolerances
  rstate.rtol_spec = integrator::rtol_spec;
  rstate.atol_spec = integrator::atol_spec;
  rstate.rtol_echem = integrator::rtol_echem;
  rstate.atol_echem = integrator::atol_echem;

  // Jacobian type: 1 = analytical
  rstate.jacobian_type = integrator::jacobian;

  // record the trajectory (one entry per accepted step) for output
  rstate.record_history = true;

  std::cout << "Stiff Chemistry Network Integration (Rosenbrock Method)" << std::endl;
  std::cout << "=========================================================\n" << std::endl;

  std::cout << "Initial state:" << std::endl;
  std::cout << "  NumSpec: " << params::NumSpec << std::endl;
  std::cout << "  Temperature: " << std::scientific << std::setprecision(3) << burn_state.get_T() << " K" << std::endl;
  std::cout << "  Energy: " << burn_state.e << " erg/cm^3" << std::endl;
  std::cout << "  Integration: " << rstate.t << " to " << rstate.tout << " s\n" << std::endl;

  // Integrate
  int ierr = rosenbrock_integrator(burn_state, rstate);

  // Update state
  for (int i = 0; i < params::NumSpec; ++i) {
    burn_state.rho[i] = rstate.y[i];
  }
  burn_state.e = rstate.y[params::NumSpec];

  std::cout << "Final state (t = " << std::scientific << rstate.t << " s):" << std::endl;
  std::cout << "  Temperature: " << burn_state.get_T() << " K" << std::endl;
  std::cout << "  Energy: " << burn_state.e << " erg/cm^3" << std::endl;

  std::cout << "\nIntegration statistics:" << std::endl;
  std::cout << "  Error code: " << static_cast<int>(ierr) << std::endl;
  std::cout << "  RHS calls:  " << rstate.n_rhs << std::endl;
  std::cout << "  Jac evals:  " << rstate.n_jac << std::endl;
  std::cout << "  Steps:      " << rstate.n_step << std::endl;

  // Write the recorded trajectory: time, species densities, temperature
  {
    std::ofstream out("evolution.dat");
    out << "# t";
    for (int i = 0; i < params::NumSpec; ++i) out << " n" << i;
    out << " T\n";
    out << std::scientific << std::setprecision(8);
    for (std::size_t k = 0; k < rstate.hist_t.size(); ++k) {
      out << rstate.hist_t[k];
      double n_tot = 0.0;
      for (int i = 0; i < params::NumSpec; ++i) {
        out << " " << rstate.hist_y[k][i];
        n_tot += rstate.hist_y[k][i];
      }
      out << " "
          << eos::temperature_from_energy(rstate.hist_y[k][params::NumSpec], n_tot)
          << "\n";
    }
    std::cout << "  Trajectory: " << rstate.hist_t.size()
              << " points -> evolution.dat" << std::endl;
  }

  if (ierr == IERR_SUCCESS) {
    std::cout << "\n✓ Integration successful!" << std::endl;
    return 0;
  } else {
    std::cerr << "\n✗ Integration failed" << std::endl;
    return 1;
  }
}
