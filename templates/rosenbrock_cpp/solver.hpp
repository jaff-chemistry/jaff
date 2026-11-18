// ABOUTME: Header-only Rosenbrock integrator (RODAS5P) for the generated chemistry:
// ABOUTME: chem::solve advances one cell's densities and specific energy by dt
//
// Include this header from the host code:
//
//   chem::Parameters p;                // set every member (chem::unset_parameters)
//   double n[chem::nspecs] = {...};    // number densities [cm^-3]
//   double e = chem::energy(n, T);     // specific internal energy [erg g^-1]
//   chem::Status s = chem::solve(n, e, dt, p);
//
// solve() works in place and leaves n and e unchanged when it fails. It keeps no
// state between calls and allocates nothing, so cells can be solved in parallel.
#ifndef JAFF_SOLVER_HPP
#define JAFF_SOLVER_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "chemistry.hpp"
#include "tableau.hpp"

namespace chem {

enum class Status : int {
  success = 0,
  unset_parameters = 1,       // a parameter is unset (see unset_parameters())
  bad_input = 2,              // dt < 0, non-finite state, or no temperature for e
  too_many_steps = 3,         // Options::max_steps reached before dt
  step_underflow = 4,         // the step size fell below round-off
  singular_matrix = 5,        // the linear system stayed singular after retries
  bad_tolerances = 6,         // rtol or atol not positive
};

// sparse: LU on the Jacobian's sparsity pattern with a fixed pivot order (no
//   pivoting; a zero pivot rejects the step). dense: LU with partial pivoting.
enum class LinearSolver { sparse, dense };

struct Options {
  double rtol = 1e-12;          // relative tolerance
  double atol = 1e-20;          // absolute tolerance [cm^-3, erg g^-1]
  long max_steps = 1000000;     // accepted + rejected steps per call
  bool evolve_energy = true;    // false: e is held fixed (densities only)
  LinearSolver lu = LinearSolver::sparse;
  double h_init = 0.0;          // first step [s]; <= 0: dt
};

struct Stats {
  long steps = 0;     // accepted steps
  long rejected = 0;  // rejected steps
  long rhs = 0;       // rhs evaluations
  long jacobians = 0; // Jacobian evaluations
};

namespace detail {

using Vector = std::array<double, neq>;
using Tableau = rosenbrock::rodas5p_tableau;

// Dense LU with partial pivoting, in place. Returns false if singular.
inline bool lu_decompose(Matrix &a, std::array<int, neq> &pivot) {
  for (int k = 0; k < neq; ++k) {
    int p = k;
    for (int i = k + 1; i < neq; ++i) {
      if (std::abs(a[i][k]) > std::abs(a[p][k])) {
        p = i;
      }
    }
    pivot[k] = p;
    if (a[p][k] == 0.0 || !std::isfinite(a[p][k])) {
      return false;
    }
    if (p != k) {
      std::swap(a[p], a[k]);
    }
    const double inv = 1.0 / a[k][k];
    for (int i = k + 1; i < neq; ++i) {
      a[i][k] *= inv;
      const double f = a[i][k];
      if (f != 0.0) {
        for (int j = k + 1; j < neq; ++j) {
          a[i][j] -= f * a[k][j];
        }
      }
    }
  }

  return true;
}

inline void lu_solve(const Matrix &a, const std::array<int, neq> &pivot, Vector &b) {
  for (int k = 0; k < neq; ++k) {
    std::swap(b[k], b[pivot[k]]);
    for (int i = k + 1; i < neq; ++i) {
      b[i] -= a[i][k] * b[k];
    }
  }
  for (int i = neq - 1; i >= 0; --i) {
    for (int j = i + 1; j < neq; ++j) {
      b[i] -= a[i][j] * b[j];
    }
    b[i] /= a[i][i];
  }
}

// Elimination schedule of the sparse LU, built once from jacobian_pattern and
// read-only afterwards. The matrix W = I/(gamma h) - J has the Jacobian's pattern
// plus the diagonal; pivots follow a minimum-degree order of that pattern
// (symmetrised), and for pivot k the schedule lists the rows below it and the
// columns right of it that are nonzero once fill-in is included. Factorisation
// and solves work in the dense Matrix but touch only these entries.
struct SparseSchedule {
  std::vector<int> order;                     // pivot k -> state index
  std::vector<int> lower_start{0}, lower;     // rows i with L(i, order[k]) != 0
  std::vector<int> upper_start{0}, upper;     // columns j with U(order[k], j) != 0
};

inline SparseSchedule build_schedule() {
  std::vector<std::vector<char>> nz(neq, std::vector<char>(neq, 0));
  for (const auto &entry : jacobian_pattern) {
    nz[entry[0]][entry[1]] = 1;
  }
  for (int i = 0; i < neq; ++i) {
    nz[i][i] = 1;
  }

  // Minimum-degree order on the symmetrised pattern (ties: lowest index)
  auto sym = nz;
  for (int i = 0; i < neq; ++i) {
    for (int j = 0; j < neq; ++j) {
      sym[i][j] = static_cast<char>(nz[i][j] || nz[j][i]);
    }
  }
  SparseSchedule s;
  std::vector<char> done(neq, 0);
  for (int step = 0; step < neq; ++step) {
    int best = -1, best_degree = neq + 1;
    for (int i = 0; i < neq; ++i) {
      if (done[i]) {
        continue;
      }
      int degree = 0;
      for (int j = 0; j < neq; ++j) {
        degree += (j != i && !done[j] && sym[i][j]) ? 1 : 0;
      }
      if (degree < best_degree) {
        best = i;
        best_degree = degree;
      }
    }
    s.order.push_back(best);
    done[best] = 1;
    std::vector<int> neighbours;
    for (int j = 0; j < neq; ++j) {
      if (!done[j] && sym[best][j]) {
        neighbours.push_back(j);
      }
    }
    for (int a : neighbours) {
      for (int b : neighbours) {
        sym[a][b] = 1;
      }
    }
  }

  // Symbolic factorisation in that order, recording the fill
  std::vector<int> position(neq);
  for (int k = 0; k < neq; ++k) {
    position[s.order[k]] = k;
  }
  for (int k = 0; k < neq; ++k) {
    const int p = s.order[k];
    std::vector<int> rows, cols;
    for (int i = 0; i < neq; ++i) {
      if (position[i] > k && nz[i][p]) {
        rows.push_back(i);
      }
      if (position[i] > k && nz[p][i]) {
        cols.push_back(i);
      }
    }
    for (int i : rows) {
      for (int j : cols) {
        nz[i][j] = 1;
      }
    }
    s.lower.insert(s.lower.end(), rows.begin(), rows.end());
    s.upper.insert(s.upper.end(), cols.begin(), cols.end());
    s.lower_start.push_back(static_cast<int>(s.lower.size()));
    s.upper_start.push_back(static_cast<int>(s.upper.size()));
  }

  return s;
}

inline const SparseSchedule &schedule() {
  static const SparseSchedule s = build_schedule();  // thread-safe, built once

  return s;
}

// Sparse LU of a, in place, without pivoting. Returns false on a zero pivot.
inline bool sparse_decompose(Matrix &a) {
  const SparseSchedule &s = schedule();
  for (int k = 0; k < neq; ++k) {
    const int p = s.order[k];
    const double pivot = a[p][p];
    if (pivot == 0.0 || !std::isfinite(pivot)) {
      return false;
    }
    for (int r = s.lower_start[k]; r < s.lower_start[k + 1]; ++r) {
      const int i = s.lower[r];
      const double l = (a[i][p] /= pivot);
      if (l == 0.0) {
        continue;
      }
      for (int c = s.upper_start[k]; c < s.upper_start[k + 1]; ++c) {
        const int j = s.upper[c];
        a[i][j] -= l * a[p][j];
      }
    }
  }

  return true;
}

inline void sparse_solve(const Matrix &a, Vector &b) {
  const SparseSchedule &s = schedule();
  for (int k = 0; k < neq; ++k) {
    const int p = s.order[k];
    for (int r = s.lower_start[k]; r < s.lower_start[k + 1]; ++r) {
      b[s.lower[r]] -= a[s.lower[r]][p] * b[p];
    }
  }
  for (int k = neq - 1; k >= 0; --k) {
    const int p = s.order[k];
    for (int c = s.upper_start[k]; c < s.upper_start[k + 1]; ++c) {
      b[p] -= a[p][s.upper[c]] * b[s.upper[c]];
    }
    b[p] /= a[p][p];
  }
}

inline bool finite(const Vector &y) {
  return std::all_of(y.begin(), y.end(), [](double v) { return std::isfinite(v); });
}

// The working storage of one solve; lives on the stack.
struct Workspace {
  Vector y{}, ynew{}, err{}, f{};
  std::array<Vector, Tableau::stages> k{};
  Matrix jac{};
  std::array<int, neq> pivot{};
};

// dy/dt with the energy row zeroed when e is held fixed
inline bool rate(const Vector &y, const Parameters &p, const Options &o, Vector &f) {
  if (!rhs(y.data(), p, f.data())) {
    return false;
  }
  if (!o.evolve_energy) {
    f[nspecs] = 0.0;
  }

  return finite(f);
}

// Jacobian of rate(); the energy row and column are zeroed when e is held fixed
inline bool jacobian(const Vector &y, const Parameters &p, const Options &o,
                     Vector &f, Matrix &jac) {
  for (auto &row : jac) {
    row.fill(0.0);
  }
  if (!rhs_jac(y.data(), p, f.data(), jac)) {
    return false;
  }
  if (!o.evolve_energy) {
    f[nspecs] = 0.0;
    for (int i = 0; i < neq; ++i) {
      jac[nspecs][i] = 0.0;
      jac[i][nspecs] = 0.0;
    }
  }

  return true;
}

inline double error_norm(const Workspace &w, const Options &o) {
  double sum = 0.0;
  for (int i = 0; i < neq; ++i) {
    const double scale = o.atol + o.rtol * std::max(std::abs(w.y[i]), std::abs(w.ynew[i]));
    sum += (w.err[i] / scale) * (w.err[i] / scale);
  }

  return std::sqrt(sum / neq);
}

// One Rosenbrock step of size h from w.y into w.ynew (error estimate in w.err).
// Returns false on a singular matrix or a state the rhs cannot evaluate.
inline bool step(Workspace &w, const Parameters &p, const Options &o, double h,
                 Stats &stats, bool &singular) {
  using C = Tableau;
  singular = false;
  if (!jacobian(w.y, p, o, w.f, w.jac)) {
    return false;
  }
  stats.jacobians += 1;
  stats.rhs += 1;

  // W = I / (gamma h) - J
  const double fac = 1.0 / (h * C::gamma);
  for (int i = 0; i < neq; ++i) {
    for (int j = 0; j < neq; ++j) {
      w.jac[i][j] = -w.jac[i][j];
    }
    w.jac[i][i] += fac;
  }
  const bool sparse = o.lu == LinearSolver::sparse;
  if (!(sparse ? sparse_decompose(w.jac) : lu_decompose(w.jac, w.pivot))) {
    singular = true;
    return false;
  }
  auto linear_solve = [&](Vector &b) {
    if (sparse) {
      sparse_solve(w.jac, b);
    } else {
      lu_solve(w.jac, w.pivot, b);
    }
  };

  w.k[0] = w.f;
  linear_solve(w.k[0]);
  for (int s = 1; s < C::stages; ++s) {
    Vector stage{};
    for (int i = 0; i < neq; ++i) {
      w.ynew[i] = w.y[i];
      for (int j = 0; j < s; ++j) {
        w.ynew[i] += C::a(s + 1, j + 1) * w.k[j][i];
        stage[i] += (C::c(s + 1, j + 1) / h) * w.k[j][i];
      }
    }
    if (!finite(w.ynew) || !rate(w.ynew, p, o, w.f)) {
      return false;
    }
    stats.rhs += 1;
    for (int i = 0; i < neq; ++i) {
      w.k[s][i] = stage[i] + w.f[i];
    }
    linear_solve(w.k[s]);
  }

  for (int i = 0; i < neq; ++i) {
    double update = 0.0;
    double err = 0.0;
    for (int s = 0; s < C::stages; ++s) {
      update += C::b(s + 1) * w.k[s][i];
      err += C::e(s + 1) * w.k[s][i];
    }
    w.ynew[i] = w.y[i] + update;
    w.err[i] = err;
  }

  return finite(w.ynew);
}

}  // namespace detail

// Advance the densities n [cm^-3] and specific internal energy e [erg g^-1] of
// one cell by dt [s]. On success n and e hold the new state; on failure they
// are unchanged and the status says why. stats, if given, receives counters.
inline Status solve(double *n, double &e, double dt, const Parameters &p,
                    const Options &o = {}, Stats *stats = nullptr) {
  // H211b step-size controller (as in the original integrator)
  constexpr double b = 4.0, kc = 2.5, fac_min = 0.2, fac_max = 6.0, reduce = 0.5;
  constexpr double err_min = 1e-10, eps = std::numeric_limits<double>::epsilon();

  Stats local;
  Stats &st = stats ? *stats : local;
  st = Stats{};
  if (!unset_parameters(p).empty()) {
    return Status::unset_parameters;
  }
  if (!(o.rtol > 10 * eps) || !(o.atol > 0.0)) {
    return Status::bad_tolerances;
  }
  if (!(dt >= 0.0)) {
    return Status::bad_input;
  }

  detail::Workspace w;
  for (int i = 0; i < nspecs; ++i) {
    w.y[i] = n[i];
  }
  w.y[nspecs] = e;
  if (!detail::finite(w.y) || std::isnan(temperature(n, e))) {
    return Status::bad_input;
  }
  if (dt == 0.0) {
    return Status::success;
  }

  double t = 0.0;
  double h = o.h_init > 0.0 ? std::min(o.h_init, dt) : dt;
  double facold = 1.0, errold = 1.0;
  bool reject = false;
  int singular_tries = 0;

  while (t < dt) {
    if (st.steps + st.rejected >= o.max_steps) {
      return Status::too_many_steps;
    }
    const bool last = t + h * (1.0 + 1e-4) >= dt;
    if (last) {
      h = dt - t;
    }
    if (0.1 * h <= t * eps || h <= 0.0) {
      return Status::step_underflow;
    }

    bool singular = false;
    if (!detail::step(w, p, o, h, st, singular)) {
      st.rejected += 1;
      if (singular && ++singular_tries >= 5) {
        return Status::singular_matrix;
      }
      h *= singular ? 0.5 : 0.25;
      reject = true;
      continue;
    }

    const double err = detail::error_norm(w, o);
    const double fac = std::clamp(std::pow(1.0 / std::max(err, err_min), 1.0 / (b * kc)) *
                                      std::pow(1.0 / std::max(errold, err_min), 1.0 / (b * kc)) *
                                      std::pow(facold, -1.0 / b),
                                  fac_min, fac_max);
    double hnew = h * fac;

    if (err <= 1.0) {
      facold = fac;
      errold = err;
      w.y = w.ynew;
      t = last ? dt : t + h;
      st.steps += 1;
      singular_tries = 0;
      if (reject) {
        hnew = std::min(hnew, h);
      }
      reject = false;
      h = std::min(hnew, dt);
    } else {
      st.rejected += 1;
      reject = true;
      h = std::min(hnew, reduce * h);
    }
  }

  for (int i = 0; i < nspecs; ++i) {
    n[i] = w.y[i];
  }
  e = w.y[nspecs];

  return Status::success;
}

}  // namespace chem

#endif
