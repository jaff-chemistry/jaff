// ABOUTME: Driver solving the generated chemistry with header-only-integrators VODE:
// ABOUTME: edit the block of inputs, then build with CMake and run build/chemistry_ode
//
// A network that calls functions it does not define needs them implemented in
// user_functions.hpp first.
#include <cmath>
#include <cstdio>
#include <string>

#include "chemistry_ode.hpp"

using namespace integrators;

int main() {
  ChemistryODE::state_type y{};

  // ---- inputs: edit here ---------------------------------------------------------
  // initial number densities [cm^-3]
  $$JAFF REPEAT idx, specie IN species
  y[$idx$] = 0.0;  // $specie$
  $$JAFF END
  // parameters: placeholders, set physical values
  $$JAFF REPEAT idx, parameter IN parameters { THERMAL { MODE -> dtdt } }
  ChemistryODE::params.$parameter$ = 1.0;
  $$JAFF END
  const double temperature = 100.0;          // initial gas temperature [K]
  const double t_end = 3.15576e13;           // end time [s] (1 Myr)
  ChemistryODE::evolve_temperature = false;  // true: evolve T with the EOS
  // VODE tolerances. This VODE fails to converge on stiff networks below about
  // rtol = 1e-10 (it drops to order 1 and stalls); 1e-8 is the tightest safe value.
  const double rtol = 1e-8, atol = 1e-20;
  const char *output = "solution.dat";
  // --------------------------------------------------------------------------------

  for (auto name : ChemistryODE::unset_parameters()) {
    std::fprintf(stderr, "parameter not set: %s\n", std::string(name).c_str());
    return 1;
  }
  y[ChemistryODE::nspecs] = temperature;
  const ChemistryODE::state_type y0 = y;

  std::FILE *out = std::fopen(output, "w");
  std::fprintf(out, "# t[s] T[K]");
  for (const auto name : ChemistryODE::species) {
    std::fprintf(out, " %s", std::string(name).c_str());
  }
  std::fprintf(out, "\n");

  // The state at t = 0, then at 100 log-spaced times from 1 s to t_end
  VODE<ChemistryODE> integrator{};
  VODEState<ChemistryODE::neqs> state{};
  state.y = y;
  state.t = 0.0;
  state.rtol = rtol;
  state.atol = atol;
  state.max_steps = 1000000;
  state.jacobian_analytic = true;
  ChemistryODE::state_type problem{};  // unused by the library
  constexpr int noutput = 100;
  for (int k = 0; k <= noutput; ++k) {
    const double t_next = k == 0 ? 0.0
                          : k == noutput ? t_end
                                         : std::pow(t_end, double(k - 1) / (noutput - 1));
    if (t_next > state.t) {
      state.tout = t_next;
      const IntegratorResult result = integrator.integrate(problem, state);
      if (result != IntegratorResult::SUCCESS) {
        std::fprintf(stderr, "VODE failed with result %d at t = %g s\n",
                     static_cast<int>(result), state.t);
        return 2;
      }
    }
    std::fprintf(out, "%.16e %.16e", t_next, state.y[ChemistryODE::nspecs]);
    for (int i = 0; i < ChemistryODE::nspecs; ++i) {
      std::fprintf(out, " %.16e", state.y[i]);
    }
    std::fprintf(out, "\n");
  }
  std::fclose(out);

  std::printf("t = %.6e s, T = %.6e K\n", t_end, state.y[ChemistryODE::nspecs]);
  for (int i = 0; i < ChemistryODE::nspecs; ++i) {
    std::printf("  %10s %.6e cm^-3\n", std::string(ChemistryODE::species[i]).c_str(),
                state.y[i]);
  }

  // Element and charge conservation report
  $$JAFF SUB nelem
  constexpr int n_elements = $nelem$;
  $$JAFF END
  $$JAFF REPEAT element IN elements
  constexpr const char *element_names[n_elements] = {"$element$", };
  $$JAFF END
  $$JAFF REPEAT specie_charge IN specie_charges
  constexpr int species_charge[ChemistryODE::nspecs] = {$specie_charge$, };
  $$JAFF END
  $$JAFF REPEAT element IN element_density_matrix
  constexpr int elem_matrix[n_elements][ChemistryODE::nspecs] = {$element$, };
  $$JAFF END
  std::printf("\nElement/charge conservation:\n");
  for (int e = 0; e < n_elements; ++e) {
    long double init = 0.0L, fin = 0.0L;
    for (int j = 0; j < ChemistryODE::nspecs; ++j) {
      init += static_cast<long double>(elem_matrix[e][j]) * y0[j];
      fin += static_cast<long double>(elem_matrix[e][j]) * state.y[j];
    }
    const long double rel = init != 0.0L ? (fin - init) / init : 0.0L;
    std::printf("  %s: init=%.10e, final=%.10e, rel_delta=%.10e\n", element_names[e],
                static_cast<double>(init), static_cast<double>(fin), static_cast<double>(rel));
  }
  long double q_init = 0.0L, q_fin = 0.0L;
  for (int j = 0; j < ChemistryODE::nspecs; ++j) {
    q_init += static_cast<long double>(species_charge[j]) * y0[j];
    q_fin += static_cast<long double>(species_charge[j]) * state.y[j];
  }
  const long double q_rel = q_init != 0.0L ? (q_fin - q_init) / q_init : 0.0L;
  std::printf("  net_charge: init=%.10e, final=%.10e, rel_delta=%.10e\n",
              static_cast<double>(q_init), static_cast<double>(q_fin),
              static_cast<double>(q_rel));

  std::printf("\nevolution saved to %s\n", output);

  return 0;
}
