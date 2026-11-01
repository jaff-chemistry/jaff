// This file has been Adapted from AMReX Microphysics / AMReX-Astrophysics:
//   https://github.com/AMReX-Astro/Microphysics
// Original LINPACK routines dgefa and dgesl have
// been modified for standard C++

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

#endif
