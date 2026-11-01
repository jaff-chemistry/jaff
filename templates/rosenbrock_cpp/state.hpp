#ifndef ROSENBROCK_STATE_H
#define ROSENBROCK_STATE_H

#include "integrator_data.hpp"
#include "parameters.hpp"
#include "eos.hpp"

template <typename T> struct state {
  DArray1D rho = DArray1D(params::NumSpec);  // Species number densities [cm^-3]
  double e = 0.0;                            // Thermal energy density [erg/cm^3]

  // Total number density of all species
  double n_tot() const {
    double n = 0.0;
    for (int i = 0; i < params::NumSpec; ++i) {
      n += rho[i];
    }
    return n;
  }

  // Temperature derived from the energy density via the EOS
  double get_T() const {
    return eos::temperature_from_energy(e, n_tot());
  }
};

#endif
