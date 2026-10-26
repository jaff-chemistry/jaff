#ifndef PARAMETERS_H
#define PARAMETERS_H

namespace params {

$$JAFF SUB nspec
constexpr int NumSpec = $nspec$;
constexpr int IE = $nspec$;
$$JAFF END

constexpr double T_init = 1.0e6; // Initial temperature

inline constexpr double n_init[NumSpec] = {
    $$JAFF REPEAT idx, specie IN species
    0.0, // [$idx$] $specie$
    $$JAFF END
};

} // namespace params

namespace integrator {

// these are the parameters for the integration module

// Whether to use an analytical or numerical Jacobian.
// 1 == Analytical
// 2 == Numerical
constexpr int jacobian = 1;

// Tolerances for the solver (relative and absolute), for the
// species and energy equations.

// relative tolerance for species
constexpr double rtol_spec = 1.e-6;

// relative tolerance for energy
constexpr double rtol_echem= 1.e-6;

// absolute tolerance for species
constexpr double atol_spec = 1.e-20;

// absolute tolerance for energy
constexpr double atol_echem= 1.e-20;

// The absolute cutoff for species -- note that this might be larger
// than ``small_x``, but the issue is that we need to prevent underflow
// issues and keep mass fractions positive in the integrator.  You may
// have to increase the floor to, e.g. 1.e-20 if your rates are large.
constexpr double SMALL_X_SAFE = 1.0e-30;

// The maximum temperature for reactions in the integration.
constexpr double MAX_TEMP = 1.0e11;

// maximum number of timesteps for the integrator
constexpr int ode_max_steps = 150000;

// maximum timestep for the integrator
constexpr double ode_max_dt = 1.e30;

// flag for turning on the use of number densities for all species
constexpr int use_number_densities = 0;

// for the linear algebra, do we allow pivoting?
constexpr int linalg_do_pivoting = 1;

// LU factorization order (both give bitwise-identical results):
//   0 = right-looking LINPACK dgefa/dgesl (vectorizes well on CPUs)
//   1 = left-looking Crout with the matrix shift fused in (for GPU lanes)
constexpr int linalg_left_looking = 0;

// for the step rejection logic on mass fractions, we only consider
// species that are > X_reject_buffer * atol_spec
constexpr double X_reject_buffer = 1.0;

// Rosenbrock tableau selector:
//   0 = Rodas5P 8-stage method from DifferentialEquations.jl
//   1 = Rodas4P 6-stage method from DifferentialEquations.jl
//   2 = Rodas3P 5-stage method from DifferentialEquations.jl
//   3 = existing 3-stage ROS2S tableau
//   4 = Verwer et al. (1999) 2-stage ROS2 tableau
//   5 = Rosenbrock-Euler 1-stage method
constexpr int rosenbrock_tableau = 0;

// Timestep controller selector:
//   0 = H211b error-history controller
//   1 = Khokhlov/YASS concentration-change controller
constexpr int rosenbrock_timestep_controller = 0;

// H211b timestep controller parameters
constexpr double h211b_b = 4.0;
constexpr double h211b_k = 2.5;
constexpr double h211b_fac_min = 0.2;
constexpr double h211b_fac_max = 6.0;
constexpr double h211b_reduction_fac = 0.5;

// Khokhlov/YASS timestep controller parameters. This controller accepts a
// step when each active species changes by no more than yass_epsilon
// relative to its starting value. Species with concentration <= yass_floor
// times the most abundant species do not constrain the step.
constexpr double yass_epsilon = 0.1;
constexpr double yass_floor = 1.e-3;
constexpr double yass_fac_min = 0.2;
constexpr double yass_fac_max = 6.0;
constexpr double yass_reduction_fac = 0.5;

} // namespace integrator

#endif
