// ABOUTME: Generated chemistry of one network: species, parameters, EOS, and the rhs and
// ABOUTME: Jacobian of the densities and the specific internal energy
//
// The state y is the species number densities n [cm^-3] followed by the specific
// internal energy e [erg g^-1]. The generated expressions are written in the
// densities (nden[i]) and the gas temperature (tgas), which is recovered from e
// through the EOS. Functions the network calls but does not define come from
// user_functions.hpp.
#ifndef JAFF_CHEMISTRY_HPP
#define JAFF_CHEMISTRY_HPP

#include <array>
#include <cmath>
#include <limits>
#include <string_view>
#include <vector>

#include "user_functions.hpp"

// The ODE system: densities and specific internal energy
$$JAFF DEFINE system {
    THERMAL { MODE -> dedt  ENERGY -> specific }
}

namespace chem {

$$JAFF SUB nspec
constexpr int nspecs = $nspec$;  // number of species
$$JAFF END
constexpr int neq = nspecs + 1;  // equations: densities, then e

using Matrix = std::array<std::array<double, neq>, neq>;

// Species names, in state order
constexpr std::array<std::string_view, nspecs> species = {
    $$JAFF REPEAT idx, specie IN species
    "$specie$",
    $$JAFF END
};

// Species masses [g], in state order
constexpr std::array<double, nspecs> masses = {
    $$JAFF REPEAT idx, specie_mass IN specie_masses
    $specie_mass$,
    $$JAFF END
};

// Index of a species by name, or -1
inline int species_index(std::string_view name) {
  for (int i = 0; i < nspecs; ++i) {
    if (species[i] == name) {
      return i;
    }
  }

  return -1;
}

constexpr double unset = std::numeric_limits<double>::quiet_NaN();

// Parameters the rhs and Jacobian read; set every one before solving
struct Parameters {
  $$JAFF REPEAT idx, parameter IN parameters { USE system }
  double $parameter$ = unset;
  $$JAFF END
};

// Names of the parameters still unset
inline std::vector<std::string_view> unset_parameters(const Parameters &p) {
  std::vector<std::string_view> names;
  $$JAFF REPEAT idx, parameter IN parameters { USE system }
  if (std::isnan(p.$parameter$)) names.push_back("$parameter$");
  $$JAFF END

  return names;
}

// Mass density [g cm^-3]
inline double density(const double *nden) {
  double rho = 0.0;
  for (int i = 0; i < nspecs; ++i) {
    rho += masses[i] * nden[i];
  }

  return rho;
}

// Internal energy per volume E(T, n) [erg cm^-3] from the EOS
inline double internal_energy(const double *nden, double tgas) {
  $$JAFF SUB eos
  return $eos$;
  $$JAFF END
}

// dE/dT at fixed densities [erg cm^-3 K^-1]
inline double heat_capacity(const double *nden, double tgas) {
  $$JAFF SUB eos_dt
  return $eos_dt$;
  $$JAFF END
}

// Specific internal energy e [erg g^-1] at temperature T [K]
inline double energy(const double *nden, double tgas) {
  return internal_energy(nden, tgas) / density(nden);
}

// Temperature [K] with specific internal energy e: Newton on E(T, n) = rho e.
// Returns NaN if it does not converge to a positive temperature.
inline double temperature(const double *nden, double e) {
  const double target = density(nden) * e;
  double tgas = target / heat_capacity(nden, 1.0);  // exact for E linear in T
  for (int iter = 0; iter < 100; ++iter) {
    if (!(tgas > 0.0) || !std::isfinite(tgas)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const double step = (internal_energy(nden, tgas) - target) / heat_capacity(nden, tgas);
    tgas -= step;
    if (std::abs(step) <= 1e-14 * std::abs(tgas)) {
      return tgas > 0.0 ? tgas : std::numeric_limits<double>::quiet_NaN();
    }
  }

  return std::numeric_limits<double>::quiet_NaN();
}

// dy/dt at state y (densities, then e). Returns false if the temperature cannot
// be recovered from e.
inline bool rhs(const double *y, const Parameters &p, double *ydot) {
  const double *nden = y;
  const double tgas = temperature(nden, y[nspecs]);
  if (std::isnan(tgas)) {
    return false;
  }
  $$JAFF REPEAT idx, parameter IN parameters { USE system }
  [[maybe_unused]] const double $parameter$ = p.$parameter$;
  $$JAFF END

  $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rhs } }
  const double cse$idx$ = $cse$;
  $$JAFF END
  $$JAFF REPEAT idx, rhs IN rhses { USE system  CSE { GROUP -> rhs } }
  ydot[$idx$] = $rhs$;
  $$JAFF END

  return true;
}

// dy/dt and the Jacobian split as J = A + a c^T, sharing every subexpression:
// A = d(dy/dt)/dy at fixed temperature (sparse; its energy column is empty),
// a = d(dy/dt)/dT and c = dT/dy (dT/dn_j at fixed e, then dT/de). The rank-one
// term is what makes J dense (T depends on every density); solve() keeps it apart.
// Writes only nonzero entries: zero A, a and c first. Returns false if the
// temperature cannot be recovered from e.
inline bool rhs_jac_split(const double *y, const Parameters &p, double *ydot, Matrix &A,
                          double *a, double *c) {
  const double *nden = y;
  const double tgas = temperature(nden, y[nspecs]);
  if (std::isnan(tgas)) {
    return false;
  }
  $$JAFF REPEAT idx, parameter IN parameters { USE system }
  [[maybe_unused]] const double $parameter$ = p.$parameter$;
  $$JAFF END

  $$JAFF REPEAT idx, cse IN cses { GROUP { NAME -> rhs_jac } }
  const double cse$idx$ = $cse$;
  $$JAFF END
  $$JAFF REPEAT idx, rhs IN rhses { USE system  CSE { GROUP -> rhs_jac } }
  ydot[$idx$] = $rhs$;
  $$JAFF END
  $$JAFF REPEAT idx, jac IN jacobian_fixed_t { USE system  CSE { GROUP -> rhs_jac } }
  A[$idx$][$idx$] = $jac$;
  $$JAFF END
  $$JAFF REPEAT idx, rhs_dt IN rhs_dt { USE system  CSE { GROUP -> rhs_jac } }
  a[$idx$] = $rhs_dt$;
  $$JAFF END
  $$JAFF REPEAT idx, dtdy IN temperature_gradient { USE system  CSE { GROUP -> rhs_jac } }
  c[$idx$] = $dtdy$;
  $$JAFF END

  return true;
}

// dy/dt and the full Jacobian J = A + a c^T at state y. Returns false if the
// temperature cannot be recovered from e.
inline bool rhs_jac(const double *y, const Parameters &p, double *ydot, Matrix &jac) {
  double a[neq] = {}, c[neq] = {};
  for (auto &row : jac) {
    row.fill(0.0);
  }
  if (!rhs_jac_split(y, p, ydot, jac, a, c)) {
    return false;
  }
  for (int i = 0; i < neq; ++i) {
    if (a[i] != 0.0) {
      for (int j = 0; j < neq; ++j) {
        jac[i][j] += a[i] * c[j];
      }
    }
  }

  return true;
}

// The structurally nonzero entries {row, column} of A, as rhs_jac_split writes them
inline constexpr int jacobian_pattern[][2] = {
    $$JAFF REPEAT idx, jac IN jacobian_fixed_t { USE system  CSE { GROUP -> rhs_jac } }
    {$idx$, $idx$},
    $$JAFF END
};

}  // namespace chem

#endif
