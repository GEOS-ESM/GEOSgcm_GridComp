#include "MAPL_Generic.h"
program CreateInterpWeights

  !--------------------------------------------------------------------!
  ! purpose: driver for MPI-IO test module                             !
  !--------------------------------------------------------------------!
  use MAPL_Mod
  use CreateInterpWeights_GridCompMod,      only: SetServices

  implicit none

  integer :: status

   call MAPL_CAP(SetServices, rc=STATUS)

end program CreateInterpWeights

