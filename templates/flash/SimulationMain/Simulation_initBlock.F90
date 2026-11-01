subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY : sim_initDens, sim_initTemp
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Simulation_interface, ONLY : Simulation_mapStrToInt

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockID

  integer :: i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, pointer, dimension(:,:,:,:) :: solnData

  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
  call Grid_getBlkPtr(blockID, solnData)

  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)

           solnData(DENS_VAR,i,j,k) = sim_initDens
           solnData(TEMP_VAR,i,j,k) = sim_initTemp
           solnData(VELX_VAR,i,j,k) = 0.0
           solnData(VELY_VAR,i,j,k) = 0.0
           solnData(VELZ_VAR,i,j,k) = 0.0

           ! Initial mass fractions, one line per species, resolved by name
           ! (FLASH orders SPECIES alphabetically, not in network order).
           ! Fix the values by hand.
           $$JAFF REPEAT idx, specie_with_normalized_sign IN species_with_normalized_sign {
               SIGN { POS -> j  NEG -> k }
           }
           solnData(spidx("$specie_with_normalized_sign$"),i,j,k) = 0.0   ! $idx$
           $$JAFF END

        enddo
     enddo
  enddo

  call Grid_releaseBlkPtr(blockID, solnData)

contains

  integer function spidx(nm)
    character(len=*), intent(in) :: nm
    call Simulation_mapStrToInt(nm, spidx, MAPBLOCK_UNK)
  end function spidx

end subroutine Simulation_initBlock
