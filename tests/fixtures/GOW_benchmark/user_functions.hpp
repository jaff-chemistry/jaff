// ABOUTME: GOW's CO-cooling and dust tables (from networks/GOW/GOW.hdf5) for the generated
// ABOUTME: rosenbrock_cpp and kokkos_ode code: replaces their user_functions.hpp (GOW benchmark)
//
// The C++ counterpart of tests/gow_tables.py (same arrays, same lookups): values
// are log-encoded (L0CO, LLTECO and the dust logps as -log10, nhalfCO as +log10,
// alphaCO direct); the CO tables are stored [Neff][T] with axes T [K] and log10
// N_CO,eff, the dust table [nH][Tg] with axes log10 Tg and log10 nH. Lookups
// interpolate linearly in log10 of the arguments and clamp to the table edges;
// the CO second argument is floored at 1e-13 before the log. CO self-shielding
// is 1 and every partial derivative is 0 (approximate Jacobian).
#ifndef JAFF_USER_FUNCTIONS_HPP
#define JAFF_USER_FUNCTIONS_HPP

#include <algorithm>
#include <cmath>

namespace gow_tables {

inline constexpr double co_t[11] = {
    10.0, 20.0, 30.0, 50.0,
    80.0, 100.0, 300.0, 600.0,
    1000.0, 1500.0, 2000.0};
inline constexpr double co_log_neff[11] = {
    14.0, 14.5, 15.0, 15.5,
    16.0, 16.5, 17.0, 17.5,
    18.0, 18.5, 19.0};
inline constexpr double co_l0[11] = {
    24.77, 24.38, 24.21, 24.03,
    23.89, 23.82, 23.34238089, 22.99832519,
    22.75384686, 22.56640625, 22.43740866};
inline constexpr double co_llte[121] = {
    21.08, 21.09, 21.11, 21.18,
    21.37, 21.67, 22.04, 22.44,
    22.87, 23.3, 23.76, 20.35,
    20.35, 20.37, 20.4, 20.51,
    20.73, 21.05, 21.42, 21.82,
    22.23, 22.66, 19.94, 19.95,
    19.96, 19.98, 20.05, 20.23,
    20.52, 20.86, 21.24, 21.65,
    22.06, 19.45, 19.45, 19.46,
    19.47, 19.52, 19.64, 19.87,
    20.19, 20.55, 20.94, 21.35,
    19.01, 19.01, 19.01, 19.02,
    19.05, 19.13, 19.32, 19.6,
    19.95, 20.32, 20.71, 18.8,
    18.8, 18.8, 18.81, 18.83,
    18.9, 19.06, 19.33, 19.66,
    20.03, 20.42, 17.81, 17.81,
    17.81, 17.82, 17.82, 17.85,
    17.92, 18.08, 18.34, 18.67,
    19.03, 17.23, 17.23, 17.23,
    17.23, 17.23, 17.25, 17.28,
    17.38, 17.59, 17.89, 18.26,
    16.86, 16.86, 16.86, 16.87,
    16.87, 16.88, 16.9, 16.97,
    17.15, 17.48, 17.93, 16.66,
    16.66, 16.66, 16.66, 16.66,
    16.67, 16.69, 16.75, 16.91,
    17.26, 17.74, 16.55, 16.55,
    16.55, 16.55, 16.55, 16.56,
    16.58, 16.63, 16.78, 17.12,
    17.61};
inline constexpr double co_nhalf[121] = {
    3.29, 3.27, 3.22, 3.07,
    2.72, 2.24, 1.74, 1.24,
    0.742, 0.242, -0.258, 3.49,
    3.48, 3.45, 3.34, 3.09,
    2.65, 2.15, 1.65, 1.15,
    0.652, 0.152, 3.67, 3.66,
    3.64, 3.56, 3.35, 2.95,
    2.47, 1.97, 1.47, 0.966,
    0.466, 3.97, 3.96, 3.94,
    3.89, 3.74, 3.42, 2.95,
    2.45, 1.95, 1.45, 0.95,
    4.3, 4.3, 4.29, 4.26,
    4.16, 3.92, 3.49, 3.0,
    2.5, 2.0, 1.5, 4.46,
    4.45, 4.45, 4.42, 4.34,
    4.14, 3.74, 3.25, 2.75,
    2.25, 1.75, 5.17, 5.16,
    5.16, 5.15, 5.13, 5.06,
    4.86, 4.47, 3.98, 3.48,
    2.98, 5.47, 5.47, 5.47,
    5.46, 5.45, 5.41, 5.3,
    5.02, 4.57, 4.07, 3.57,
    5.53, 5.53, 5.53, 5.52,
    5.51, 5.48, 5.39, 5.16,
    4.73, 4.24, 3.74, 5.3,
    5.3, 5.3, 5.3, 5.29,
    5.26, 5.17, 4.94, 4.52,
    4.03, 3.53, 4.7, 4.7,
    4.7, 4.7, 4.68, 4.64,
    4.53, 4.27, 3.84, 3.35,
    2.85};
inline constexpr double co_alpha[121] = {
    0.439, 0.436, 0.428, 0.416,
    0.416, 0.45, 0.492, 0.529,
    0.555, 0.582, 0.596, 0.409,
    0.407, 0.401, 0.388, 0.378,
    0.396, 0.435, 0.473, 0.503,
    0.528, 0.546, 0.392, 0.391,
    0.385, 0.373, 0.36, 0.367,
    0.403, 0.441, 0.473, 0.499,
    0.519, 0.37, 0.368, 0.364,
    0.353, 0.338, 0.334, 0.362,
    0.404, 0.44, 0.469, 0.492,
    0.361, 0.359, 0.356, 0.347,
    0.332, 0.322, 0.339, 0.381,
    0.423, 0.457, 0.483, 0.357,
    0.356, 0.352, 0.345, 0.33,
    0.317, 0.329, 0.37, 0.414,
    0.451, 0.479, 0.385, 0.385,
    0.383, 0.38, 0.371, 0.355,
    0.343, 0.362, 0.418, 0.47,
    0.51, 0.437, 0.437, 0.436,
    0.434, 0.429, 0.419, 0.406,
    0.41, 0.446, 0.487, 0.516,
    0.428, 0.427, 0.427, 0.425,
    0.421, 0.414, 0.401, 0.392,
    0.404, 0.432, 0.448, 0.354,
    0.354, 0.352, 0.349, 0.341,
    0.329, 0.317, 0.316, 0.335,
    0.364, 0.372, 0.322, 0.322,
    0.32, 0.316, 0.307, 0.292,
    0.276, 0.272, 0.289, 0.31,
    0.313};
inline constexpr double dust_log_t[10] = {
    0.5, 0.88888889, 1.27777778, 1.66666667,
    2.05555556, 2.44444444, 2.83333333, 3.22222222,
    3.61111111, 4.0};
inline constexpr double dust_log_nh[15] = {
    0.0, 0.42857143, 0.85714286, 1.28571429,
    1.71428571, 2.14285714, 2.57142857, 3.0,
    3.42857143, 3.85714286, 4.28571429, 4.71428571,
    5.14285714, 5.57142857, 6.0};
inline constexpr double dust_ps[150] = {
    33.60923439, 32.35048647, 31.6458604, 31.02132235,
    30.42222289, 29.83261, 29.24673384, 28.66235604,
    28.0785789, 27.49504472, 33.18091039, 31.92216147,
    31.21753302, 30.59298934, 29.99387726, 29.40423904,
    28.8183211, 28.23389204, 27.65007024, 27.06650687,
    32.75300149, 31.49424423, 30.78959638, 30.16501048,
    29.5658181, 28.97605855, 28.39000596, 27.80546791,
    27.22157705, 26.6379757, 32.32619855, 31.06738031,
    30.36260198, 29.73777592, 29.13823862, 28.54811853,
    27.96178914, 27.37708183, 26.79309991, 26.20945204,
    31.90230918, 30.64308622, 29.93758718, 29.3117745,
    28.71125923, 28.12042115, 27.53366605, 26.94873467,
    26.36464058, 25.78093706, 31.48587783, 30.22433325,
    29.51586354, 28.88730605, 28.28488476, 27.69295575,
    27.10563836, 26.52043033, 25.93620173, 25.35243225,
    31.0872853, 29.81525158, 29.09829079, 28.46439304,
    27.85909326, 27.26572717, 26.67771529, 26.09217504,
    25.50778679, 24.92393938, 30.72591396, 29.41932635,
    28.68506438, 28.0430074, 27.4339034, 26.83875916,
    26.24991203, 25.66397699, 25.07939997, 24.49546057,
    30.42733764, 29.03861249, 28.27630163, 27.62323376,
    27.00938204, 26.41209045, 25.82224861, 25.23584618,
    24.65104624, 24.06699834, 30.21146872, 28.67535786,
    27.87248634, 27.20529029, 26.58563564, 25.9857724,
    25.39474974, 24.80779463, 24.22273154, 23.63855567,
    30.08004755, 28.33375429, 27.47457781, 26.78951505,
    26.16280582, 25.55986866, 24.96744513, 24.37983659,
    23.79446288, 23.21013604, 30.0136033, 28.02095875,
    27.08405993, 26.37636621, 25.74107062, 25.13445646,
    24.54037027, 23.95198897, 23.36624852, 22.7817436,
    29.98466705, 27.74775444, 26.70307122, 25.96644194,
    25.32065037, 24.70962891, 24.1135674, 23.52427179,
    22.93809825, 22.35338321, 29.97312226, 27.52792637,
    26.3346648, 25.56052005, 24.90181737, 24.28549826,
    23.68708687, 23.09670875, 22.51002362, 21.92506064,
    29.96870043, 27.37312455, 25.98324711, 25.15962125,
    24.48490972, 23.86220014, 23.26098873, 22.66932799,
    22.08203828, 21.49678269};

inline double safe_log10(double x) { return std::log10(std::max(x, 1e-300)); }

// Interval i (0-based, 0..n-2) holding x and its fraction f, x clamped to the axis
inline void locate(double x, const double *axis, int n, int &i, double &f) {
  const double xc = std::min(std::max(x, axis[0]), axis[n - 1]);
  i = 0;
  while (i < n - 2 && axis[i + 1] < xc) {
    ++i;
  }
  f = (xc - axis[i]) / (axis[i + 1] - axis[i]);
}

inline double linear(double x, const double *axis, const double *values, int n) {
  int i;
  double f;
  locate(x, axis, n, i, f);
  return (1 - f) * values[i] + f * values[i + 1];
}

// z stored row-major [y][x]
inline double bilinear(double x, double y, const double *xs, int nx, const double *ys,
                       int ny, const double *z) {
  int i, j;
  double fx, fy;
  locate(x, xs, nx, i, fx);
  locate(y, ys, ny, j, fy);
  const double low = (1 - fx) * z[j * nx + i] + fx * z[j * nx + i + 1];
  const double high = (1 - fx) * z[(j + 1) * nx + i] + fx * z[(j + 1) * nx + i + 1];
  return (1 - fy) * low + fy * high;
}

inline const double *log_co_t() {
  static const auto axis = [] {
    struct { double v[11]; } a{};
    for (int k = 0; k < 11; ++k) a.v[k] = std::log10(co_t[k]);
    return a;
  }();
  return axis.v;
}

inline double co_lookup(const double *table, double tgas, double lvg) {
  return bilinear(safe_log10(tgas), safe_log10(lvg + 1e-13), log_co_t(), 11, co_log_neff,
                  11, table);
}

}  // namespace gow_tables

inline double l0_co_interp1d(double a0) {
  return std::pow(10.0, -gow_tables::linear(gow_tables::safe_log10(a0),
                                            gow_tables::log_co_t(), gow_tables::co_l0, 11));
}
inline double llte_co_interp2d(double a0, double a1) {
  return std::pow(10.0, -gow_tables::co_lookup(gow_tables::co_llte, a0, a1));
}
inline double nhalf_co_interp2d(double a0, double a1) {
  return std::pow(10.0, gow_tables::co_lookup(gow_tables::co_nhalf, a0, a1));
}
inline double alpha_co_interp2d(double a0, double a1) {
  return gow_tables::co_lookup(gow_tables::co_alpha, a0, a1);
}
inline double psigd_coll_interp2d(double a0, double a1) {  // n_H, tgas
  return std::pow(10.0, -gow_tables::bilinear(gow_tables::safe_log10(a1),
                                               gow_tables::safe_log10(a0),
                                               gow_tables::dust_log_t, 10,
                                               gow_tables::dust_log_nh, 15, gow_tables::dust_ps));
}
inline double fshield_co_interp2d(double, double) { return 1.0; }

// Partial derivatives: 0 (the Jacobian is then approximate; the solution is not)
inline double l0_co_interp1d_partial_0(double) { return 0.0; }
inline double llte_co_interp2d_partial_0(double, double) { return 0.0; }
inline double llte_co_interp2d_partial_1(double, double) { return 0.0; }
inline double nhalf_co_interp2d_partial_0(double, double) { return 0.0; }
inline double nhalf_co_interp2d_partial_1(double, double) { return 0.0; }
inline double alpha_co_interp2d_partial_0(double, double) { return 0.0; }
inline double alpha_co_interp2d_partial_1(double, double) { return 0.0; }
inline double psigd_coll_interp2d_partial_0(double, double) { return 0.0; }
inline double psigd_coll_interp2d_partial_1(double, double) { return 0.0; }
inline double fshield_co_interp2d_partial_0(double, double) { return 0.0; }
inline double fshield_co_interp2d_partial_1(double, double) { return 0.0; }

#endif
