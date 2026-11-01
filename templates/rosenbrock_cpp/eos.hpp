#ifndef ROSENBROCK_EOS_H
#define ROSENBROCK_EOS_H

namespace eos {


constexpr double k_B = 1.380649e-16;  // cgs

inline double temperature_from_energy(double e, double n_tot) {
  if (e <= 0.0 || n_tot <= 0.0) {
    return 1.0e-10;  // Floor temperature
  }
  return e / (1.5 * n_tot * k_B);
}

inline double energy_from_temperature(double T, double n_tot) {
  if (T <= 0.0 || n_tot <= 0.0) {
    return 1.0e-30;  // Floor energy
  }
  return 1.5 * n_tot * k_B * T;
}

}  // namespace eos

#endif
