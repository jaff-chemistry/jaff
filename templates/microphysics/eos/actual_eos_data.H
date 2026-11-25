// ABOUTME: Per-species data of the JAFF-generated multigamma EOS: adiabatic index and
// ABOUTME: mass of every species, filled by actual_eos_init()
#ifndef actual_eos_data_H
#define actual_eos_data_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <network.H>

inline AMREX_GPU_MANAGED amrex::Real gammas[NumSpec];
inline AMREX_GPU_MANAGED amrex::Real spmasses[NumSpec];

#endif
