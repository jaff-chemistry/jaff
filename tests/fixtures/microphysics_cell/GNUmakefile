# ABOUTME: One-zone Microphysics test of a JAFF-generated network (tests/e2e): build with
# ABOUTME: make NETWORK_CUSTOM_PATH=<run> EOS_PATH=<run>/eos MICROPHYSICS_HOME=... AMREX_HOME=...
PRECISION  = DOUBLE
PROFILE    = FALSE
DEBUG      = FALSE
DIM        = 1
COMP       = gnu
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_CUDA   = FALSE
USE_REACT  = TRUE

EBASE = main

# The generated network and its EOS come from the command line:
# NETWORK_CUSTOM_PATH=<jaffgen run folder>, EOS_PATH=<run folder>/eos
NETWORK_DIR := jaff
EOS_DIR     := jaff
CONDUCTIVITY_DIR := stellar
INTEGRATOR_DIR = VODE

EXTERN_SEARCH += .

Bpack   := ./Make.package
Blocs   := .

include $(MICROPHYSICS_HOME)/unit_test/Make.unit_test
