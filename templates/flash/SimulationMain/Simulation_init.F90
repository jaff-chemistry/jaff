subroutine Simulation_init()

  use Simulation_data, ONLY : sim_initDens, sim_initTemp
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

  call RuntimeParameters_get("sim_initDens", sim_initDens)
  call RuntimeParameters_get("sim_initTemp", sim_initTemp)

end subroutine Simulation_init
