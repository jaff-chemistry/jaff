// ABOUTME: One-zone test of a JAFF-generated Microphysics network: reads initial number
// ABOUTME: densities by species name, burns once to tmax, prints success, T and densities
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <AMReX.H>
#include <burner.H>
#include <eos.H>
#include <extern_parameters.H>
#include <network.H>
#include <unit_test.H>

int main(int argc, char *argv[]) {
  amrex::Initialize(argc, argv);
  int status = 0;
  {
    init_unit_test();
    eos_init(unit_test_rp::small_temp, unit_test_rp::small_dens);
    network_init();

    burn_t state;
    for (int n = 0; n < NumSpec; ++n) {
      state.xn[n] = 0.0_rt;
    }
    std::ifstream initial(unit_test_rp::initial_file);
    std::string name;
    amrex::Real value;
    while (initial >> name >> value) {
      const int n = network_spec_index(name);
      if (n < 0) {
        amrex::Abort("unknown species in the initial file: " + name);
      }
      state.xn[n] = value;
    }

    // constant density: the mass of the initial composition
    state.rho = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
      state.rho += state.xn[n] * spmasses[n];
    }
    state.T = unit_test_rp::temperature;
    eos(eos_input_rt, state);

    burner(state, unit_test_rp::tmax);

    std::cout << std::setprecision(17) << "success " << state.success << "\n"
              << "T " << state.T << "\n";
    for (int n = 0; n < NumSpec; ++n) {
      std::cout << spec_names_cxx[n] << " " << state.xn[n] << "\n";
    }
    status = state.success ? 0 : 1;
  }
  amrex::Finalize();

  return status;
}
