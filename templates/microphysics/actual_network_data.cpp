#include <actual_network.H>

namespace network {}

void actual_network_init() {}

void balance_charge(burn_t &state) {
  // update the number density of electrons due to charge conservation
  $$JAFF SUB e_idx
  const int e_idx = $e_idx$;
  $$JAFF END

  $$JAFF REDUCE charged_specie_index_ne, charged_specie_charge_ne IN charged_specie_indices_ne, charged_specie_charges_ne
  state.xn[e_idx] = $($charged_specie_charge_ne$*state.xn[$charged_specie_index_ne$])$;
  $$JAFF END
}
