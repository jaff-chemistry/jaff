! ABOUTME: Fills the one-zone test's cells: density sum_i n_i m_i, the initial temperature,
! ABOUTME: zero velocity and the species mass fractions n_i m_i / rho
subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Simulation_interface, ONLY : Simulation_mapStrToInt

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: blockID

  integer :: i, j, k, n, unk
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: rho
  real, pointer, dimension(:,:,:,:) :: solnData

  rho = sum(sim_initialDensity * sim_speciesMass)

  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)
  call Grid_getBlkPtr(blockID, solnData)

  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
           solnData(DENS_VAR,i,j,k) = rho
           solnData(TEMP_VAR,i,j,k) = sim_initTemp
           solnData(VELX_VAR,i,j,k) = 0.0
           solnData(VELY_VAR,i,j,k) = 0.0
           solnData(VELZ_VAR,i,j,k) = 0.0
           do n = 1, NSPECIES
              call Simulation_mapStrToInt(trim(sim_unkName(n)), unk, MAPBLOCK_UNK)
              solnData(unk,i,j,k) = sim_initialDensity(n) * sim_speciesMass(n) / rho
           end do
        end do
     end do
  end do

  call Grid_releaseBlkPtr(blockID, solnData)

end subroutine Simulation_initBlock
