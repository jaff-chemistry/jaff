// This file has been Adapted from AMReX Microphysics / AMReX-Astrophysics:
//   https://github.com/AMReX-Astro/Microphysics
// Original LINPACK routines dgefa and dgesl have
// been modified for standard C++.
//
// lu_decompose_shifted/lu_solve restructure them as a left-looking (Crout)
// factorization with LAPACK-style whole-row pivoting. For a nonsingular
// matrix each entry receives the same sequence of floating-point operations
// as in dgefa/dgesl, so the factors and solutions are bitwise identical;
// only the loop order and the storage of the row interchanges differ. The
// right-looking dgefa vectorizes better on CPUs; the left-looking form
// finishes each column before reading the next, which keeps a GPU lane's
// stores and reloads off the critical path.

#ifndef LINPACK_H
#define LINPACK_H

#include "integrator_data.hpp"
#include <cmath>

template <int num_eqs, bool allow_pivot>
inline void dgesl(const DArray2D &a, const IArray1D &pivot, DArray1D &b) {

  constexpr int nm1 = num_eqs - 1;

  // solve a * x = b
  // first solve l * y = b
  if constexpr (nm1 >= 1) {
    for (int k = 0; k < nm1; ++k) {

      double t{};
      if constexpr (allow_pivot) {
        int l = pivot[k];
        t = b[l];
        if (l != k) {
          b[l] = b[k];
          b[k] = t;
        }
      } else {
        t = b[k];
      }

      for (int j = k + 1; j < num_eqs; ++j) {
        b[j] += t * a[j][k];
      }
    }
  }

  // now solve u * x = y
  for (int kb = 0; kb < num_eqs; ++kb) {

    const int k = num_eqs - 1 - kb;
    b[k] = b[k] / a[k][k];
    double t = -b[k];
    for (int j = 0; j < k; ++j) {
      b[j] += t * a[j][k];
    }
  }
}

template <int num_eqs, bool allow_pivot>
inline void dgefa(DArray2D &a, IArray1D &pivot, int &info) {

  // dgefa factors a matrix by gaussian elimination.
  // a is returned in the form a = l * u where
  // l is a product of permutation and unit lower
  // triangular matrices and u is upper triangular.

  // gaussian elimination with partial pivoting

  info = 0;
  constexpr int nm1 = num_eqs - 1;

  double t;

  if constexpr (nm1 >= 1) {

    for (int k = 0; k < nm1; ++k) {

      // find l = pivot index
      int l = k;

      if constexpr (allow_pivot) {
        double dmax = std::abs(a[k][k]);
        for (int i = k + 1; i < num_eqs; ++i) {
          double ai = std::abs(a[i][k]);
          if (ai > dmax) {
            l = i;
            dmax = ai;
          }
        }

        pivot[k] = l;
      }

      // zero pivot implies this column already triangularized
      if (a[l][k] != 0.0e0) {

        if constexpr (allow_pivot) {
          // interchange if necessary
          if (l != k) {
            t = a[l][k];
            a[l][k] = a[k][k];
            a[k][k] = t;
          }
        }

        // compute multipliers
        t = -1.0e0 / a[k][k];
        for (int j = k + 1; j < num_eqs; ++j) {
          a[j][k] *= t;
        }

        // row elimination with column indexing
        for (int j = k + 1; j < num_eqs; ++j) {
          t = a[l][j];

          if constexpr (allow_pivot) {
            if (l != k) {
              a[l][j] = a[k][j];
              a[k][j] = t;
            }
          }

          for (int i = k + 1; i < num_eqs; ++i) {
            a[i][j] += t * a[i][k];
          }
        }

      } else {
        info = k;
      }
    }
  }

  if constexpr (allow_pivot) {
    pivot[num_eqs - 1] = num_eqs - 1;
  }

  if (a[num_eqs - 1][num_eqs - 1] == 0.0e0) {
    info = num_eqs - 1;
  }
}

// Factor M = fac * I - J in place, where a holds J on entry.
//
// Column j is formed from J, brought up to date with the interchanges and
// elimination steps k < j, then pivoted and scaled, so each column is
// finished before the next one is read and the matrix shift needs no
// separate pass. A pivot swaps whole rows of the finished columns (the
// multipliers included), so lu_solve can permute b once and then run a plain
// forward substitution.
//
// On return a holds U on and above the diagonal and the negated multipliers
// below it, pivot[k] is the row interchanged with row k at step k, and info
// is the index of the last zero pivot (as in dgefa), or 0. A zero pivot
// records no interchange; dgefa recorded one that it did not apply.
template <int num_eqs, bool allow_pivot>
inline void lu_decompose_shifted(DArray2D &a, const double fac, IArray1D &pivot,
                                 int &info) {
  info = 0;

  for (int j = 0; j < num_eqs; ++j) {
    // column j of M
    for (int i = 0; i < num_eqs; ++i) {
      a[i][j] = -a[i][j];
    }
    a[j][j] += fac;

    // Bring the column into the row order of the finished columns: their
    // multipliers have moved with every interchange so far.
    if constexpr (allow_pivot) {
      for (int k = 0; k < j; ++k) {
        const int l = pivot[k];
        if (l != k) {
          const double t = a[l][j];
          a[l][j] = a[k][j];
          a[k][j] = t;
        }
      }
    }

    // elimination steps k < j, in order
    for (int k = 0; k < j; ++k) {
      const double t = a[k][j];
      for (int i = k + 1; i < num_eqs; ++i) {
        a[i][j] += t * a[i][k];
      }
    }

    if (j == num_eqs - 1) {
      break;
    }

    // find l = pivot index
    int l = j;
    if constexpr (allow_pivot) {
      double dmax = std::abs(a[j][j]);
      for (int i = j + 1; i < num_eqs; ++i) {
        const double ai = std::abs(a[i][j]);
        if (ai > dmax) {
          l = i;
          dmax = ai;
        }
      }
    }

    // zero pivot implies this column is already triangularized
    if (a[l][j] == 0.0e0) {
      info = j;
      if constexpr (allow_pivot) {
        pivot[j] = j;
      }
      continue;
    }

    if constexpr (allow_pivot) {
      pivot[j] = l;
      if (l != j) {
        // whole-row interchange of the finished columns 0..j
        for (int c = 0; c <= j; ++c) {
          const double t = a[l][c];
          a[l][c] = a[j][c];
          a[j][c] = t;
        }
      }
    }

    // compute multipliers
    const double t = -1.0e0 / a[j][j];
    for (int i = j + 1; i < num_eqs; ++i) {
      a[i][j] *= t;
    }
  }

  if constexpr (allow_pivot) {
    pivot[num_eqs - 1] = num_eqs - 1;
  }

  if (a[num_eqs - 1][num_eqs - 1] == 0.0e0) {
    info = num_eqs - 1;
  }
}

// Solve (fac * I - J) x = b with the factors from lu_decompose_shifted.
template <int num_eqs, bool allow_pivot>
inline void lu_solve(const DArray2D &a, const IArray1D &pivot, DArray1D &b) {

  // apply the row interchanges, then solve l * y = b
  if constexpr (allow_pivot) {
    for (int k = 0; k < num_eqs - 1; ++k) {
      const int l = pivot[k];
      if (l != k) {
        const double t = b[l];
        b[l] = b[k];
        b[k] = t;
      }
    }
  }

  for (int k = 0; k < num_eqs - 1; ++k) {
    const double t = b[k];
    for (int j = k + 1; j < num_eqs; ++j) {
      b[j] += t * a[j][k];
    }
  }

  // now solve u * x = y
  for (int kb = 0; kb < num_eqs; ++kb) {
    const int k = num_eqs - 1 - kb;
    b[k] = b[k] / a[k][k];
    const double t = -b[k];
    for (int j = 0; j < k; ++j) {
      b[j] += t * a[j][k];
    }
  }
}

#endif
