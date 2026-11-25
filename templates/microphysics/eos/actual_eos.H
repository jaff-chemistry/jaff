// ABOUTME: Multigamma ideal-gas EOS of the JAFF-generated network for AMReX Microphysics:
// ABOUTME: E = sum_i n_i k T / (gamma_i - 1) with the network's gammas and masses
//
// The same internal energy as JAFF's EOS, so the energy equation in actual_rhs.H and
// the temperature found here agree. state.xn holds number densities [cm^-3]
// (integrator.use_number_densities = 1); the mass density is sum_i n_i m_i. Point
// the host's EOS directory at this folder (Microphysics: EOS_PATH=<run>/eos).
#ifndef ACTUAL_EOS_H
#define ACTUAL_EOS_H

#include <AMReX.H>
#include <cmath>
#include <string_view>

#include <actual_eos_data.H>
#include <extern_parameters.H>
#include <fundamental_constants.H>
#include <network.H>

constexpr std::string_view eos_name = "multigamma";

// Adiabatic index and mass [g] of every species, from the network
inline void actual_eos_init() {
  $$JAFF REPEAT idx, specie_gamma IN specie_gammas
  gammas[$idx$] = $specie_gamma$;
  $$JAFF END
  $$JAFF REPEAT idx, specie_mass IN specie_masses
  spmasses[$idx$] = $specie_mass$;
  $$JAFF END
}

template <typename I>
AMREX_GPU_HOST_DEVICE AMREX_INLINE bool is_input_valid(I input) {
  static_assert(std::is_same_v<I, eos_input_t>, "input must be an eos_input_t");

  return input == eos_input_rt || input == eos_input_re || input == eos_input_rp;
}

template <typename I, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE void actual_eos(I input, T &state) {
  static_assert(std::is_same_v<I, eos_input_t>, "input must be an eos_input_t");

  // n_tot [cm^-3], mass density [g cm^-3] and heat capacity per volume
  // c_v = sum_i n_i k / (gamma_i - 1) [erg cm^-3 K^-1], from the composition
  amrex::Real ntot = 0.0_rt, rhotot = 0.0_rt, cv = 0.0_rt;
  for (int n = 0; n < NumSpec; ++n) {
    ntot += state.xn[n];
    rhotot += state.xn[n] * spmasses[n];
    cv += state.xn[n] * C::k_B / (gammas[n] - 1.0_rt);
  }
  state.mu = rhotot / (ntot * C::m_u);

  amrex::Real temp = NAN;
  switch (input) {
  case eos_input_rt:
    temp = state.T;
    break;
  case eos_input_re:
    if constexpr (has_energy<T>) {
      AMREX_ASSERT(state.e > 0.);
      temp = state.e * rhotot / cv;
    }
    break;
  case eos_input_rp:
    if constexpr (has_pressure<T>) {
      AMREX_ASSERT(state.p > 0.);
      temp = state.p / (ntot * C::k_B);
    }
    break;
  default:
#ifndef AMREX_USE_GPU
    amrex::Error("EOS: input mode not supported by the JAFF multigamma EOS");
#endif
    break;
  }

  // specific internal energy e = E / rho [erg g^-1]
  state.T = temp;
  if constexpr (has_energy<T>) {
    state.e = cv * temp / rhotot;
    state.dedT = cv / rhotot;
    state.dedr = 0.0_rt;
  }
  if constexpr (has_pressure<T>) {
    const amrex::Real pressure = ntot * C::k_B * temp;
    const amrex::Real gamma1 = 1.0_rt + ntot * C::k_B / cv;
    state.p = pressure;
    state.dpdT = pressure / temp;
    state.dpdr = pressure / state.rho;
    state.cs = std::sqrt(gamma1 * pressure / state.rho);
    if constexpr (has_G<T>) {
      state.G = 0.5_rt * (1.0_rt + gamma1);
    }
    if constexpr (has_energy<T>) {
      state.dpde = state.dpdT / state.dedT;
    }
  }
}

inline void actual_eos_finalize() {}

#endif
