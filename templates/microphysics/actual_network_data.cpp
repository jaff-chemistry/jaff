// ABOUTME: Network setup of the JAFF-generated chemistry for AMReX Microphysics:
// ABOUTME: unset-parameter check at init, electron density from charge balance
#include <string>

#include <AMReX.H>
#include <actual_network.H>
#include <extern_parameters.H>

namespace network {}

void actual_network_init() {
  // -1.e300 is the "unset" default of every network parameter (_parameters)
  constexpr amrex::Real unset = -1.e300;
  std::string missing;
  $$JAFF REPEAT idx, parameter IN parameters { THERMAL { MODE -> dedt } }
  if (network_rp::$parameter$ == unset) missing += " network.$parameter$";
  $$JAFF END
  if (!missing.empty()) {
    amrex::Abort("network parameters not set in the inputs file:" + missing);
  }
}

void balance_charge(burn_t &state) {
  // update the number density of electrons due to charge conservation
  $$JAFF SUB e_idx
  const int e_idx = $e_idx$;
  $$JAFF END

  $$JAFF REDUCE charged_specie_index_ne, charged_specie_charge_ne IN charged_specie_indices_ne, charged_specie_charges_ne
  state.xn[e_idx] = $($charged_specie_charge_ne$*state.xn[$charged_specie_index_ne$])$;
  $$JAFF END
}
