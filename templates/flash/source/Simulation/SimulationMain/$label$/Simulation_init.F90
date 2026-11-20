! ABOUTME: Reads the one-zone test's settings and its file of initial number densities
! ABOUTME: (one "name value" pair per line, network species names)
subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"

  character(len=64) :: name
  real :: value
  integer :: unit, ios, n
  logical :: found

  call RuntimeParameters_get("sim_initTemp", sim_initTemp)
  call RuntimeParameters_get("sim_initialFile", sim_initialFile)
  call RuntimeParameters_get("sim_outputFile", sim_outputFile)

  $$JAFF REPEAT idx, specie IN species
  sim_speciesName($idx+1$) = "$specie$"
  $$JAFF END
  $$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign {
      SIGN { POS -> j  NEG -> k }
  }
  sim_unkName($idx+1$) = "$specie_with_normalized_sign$"
  $$JAFF END
  $$JAFF REPEAT idx, specie_mass IN specie_masses
  sim_speciesMass($idx+1$) = $specie_mass$
  $$JAFF END

  sim_initialDensity = 0.0
  inquire(file=trim(sim_initialFile), exist=found)
  if (.not. found) call Driver_abortFlash("missing sim_initialFile: " // trim(sim_initialFile))
  open(newunit=unit, file=trim(sim_initialFile), status="old", action="read")
  do
     read(unit, *, iostat=ios) name, value
     if (ios /= 0) exit
     n = findloc(sim_speciesName, trim(name), dim=1)
     if (n == 0) call Driver_abortFlash("unknown species in sim_initialFile: " // trim(name))
     sim_initialDensity(n) = value
  end do
  close(unit)

end subroutine Simulation_init
