#define I_AM_MAIN

#include "MAPL_Generic.h"

program StandAlone_AdvCore
   use MAPL_Mod
!   use AdvCore_GridCompMod, only: SetServices
   use MPI

   implicit none

!EOP

!EOC

   character(*), parameter :: IAM = __FILE__

   type (MAPL_Cap) :: cap
   type (MAPL_FlapCapOptions) :: cap_options
   integer :: status

   cap_options = MAPL_FlapCapOptions( &
        description = 'FV Standalone dvCore', &
        authors     = 'S.J. Lin, R. Rood, W. Putman')

   ! comment this out and configure in command line
   cap_options%shared_obj='libFVdycoreCubed_GridComp.so'
   cap_options%proc_name="advcore_setservices"
   !cap = MAPL_Cap('Standalone FV3 AdvCore', SetServices, cap_options=cap_options)
   cap = MAPL_Cap('Standalone FV3 AdvCore', cap_options=cap_options)
   call cap%run(_RC)

end program StandAlone_AdvCore

!EOC


