// ABOUTME: Example driver of the generated integrator: edit the block of inputs, then
// ABOUTME: `make` and `./test`; prints the final state and saves the evolution
//
// It shows the calls a host code makes: set the parameters, convert the starting
// temperature to a specific energy, and call chem::solve once per time interval.
// A network that calls functions it does not define needs them implemented in
// user_functions.hpp first.
#include <cmath>
#include <cstdio>
#include <string>

#include "solver.hpp"

int main() {
  double n[chem::nspecs] = {};
  chem::Parameters p;

  // ---- inputs: edit here ---------------------------------------------------------
  // initial number densities [cm^-3]: placeholders, set physical values (the
  // specific energy needs some mass: not every density may be 0)
  $$JAFF REPEAT idx, specie IN species
  n[$idx$] = 1.0;  // $specie$
  $$JAFF END
  // parameters: placeholders, set physical values
  $$JAFF REPEAT idx, parameter IN parameters { THERMAL { MODE -> dedt  ENERGY -> specific } }
  p.$parameter$ = 1.0;
  $$JAFF END
  const double temperature = 100.0;  // initial gas temperature [K]
  const double t_end = 3.15576e13;   // end time [s] (1 Myr)
  chem::Options options;             // tolerances, linear solver, ...
  options.evolve_energy = true;      // false: hold the specific energy fixed
  const char *output = "solution.dat";
  // --------------------------------------------------------------------------------

  double e = chem::energy(n, temperature);
  std::FILE *out = std::fopen(output, "w");
  std::fprintf(out, "# t[s] T[K]");
  for (const auto name : chem::species) {
    std::fprintf(out, " %s", std::string(name).c_str());
  }
  std::fprintf(out, "\n");

  // The state at t = 0, then at 100 log-spaced times from 1 s to t_end
  constexpr int noutput = 100;
  double t = 0.0;
  for (int k = 0; k <= noutput; ++k) {
    const double t_next = k == 0 ? 0.0
                          : k == noutput ? t_end
                                         : std::pow(t_end, double(k - 1) / (noutput - 1));
    const chem::Status status = chem::solve(n, e, t_next - t, p, options);
    if (status != chem::Status::success) {
      std::fprintf(stderr, "chem::solve failed with status %d at t = %g s\n",
                   static_cast<int>(status), t);
      for (const auto name : chem::unset_parameters(p)) {
        std::fprintf(stderr, "parameter not set: %s\n", std::string(name).c_str());
      }
      return 1;
    }
    t = t_next;
    std::fprintf(out, "%.16e %.16e", t, chem::temperature(n, e));
    for (int i = 0; i < chem::nspecs; ++i) {
      std::fprintf(out, " %.16e", n[i]);
    }
    std::fprintf(out, "\n");
  }
  std::fclose(out);

  std::printf("t = %.6e s, T = %.6e K\n", t, chem::temperature(n, e));
  for (int i = 0; i < chem::nspecs; ++i) {
    std::printf("  %10s %.6e cm^-3\n", std::string(chem::species[i]).c_str(), n[i]);
  }
  std::printf("evolution saved to %s\n", output);

  return 0;
}
