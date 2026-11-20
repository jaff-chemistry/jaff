! ABOUTME: Data of the one-zone test: initial temperature and number densities, the
! ABOUTME: output file, and the network's species names and masses in network order
module Simulation_data

#include "Flash.h"
#include "constants.h"

  implicit none

  real, save :: sim_initTemp
  character(len=MAX_STRING_LENGTH), save :: sim_initialFile, sim_outputFile
  ! initial number density [cm^-3] of every species, in network order
  real, save, dimension(NSPECIES) :: sim_initialDensity
  ! network name, FLASH unk name and mass [g] of every species, in network order
  character(len=16), save, dimension(NSPECIES) :: sim_speciesName, sim_unkName
  real, save, dimension(NSPECIES) :: sim_speciesMass

end module Simulation_data
