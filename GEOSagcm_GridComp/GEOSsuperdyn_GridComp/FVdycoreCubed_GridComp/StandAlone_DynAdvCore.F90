#define I_AM_MAIN

#include "MAPL_Generic.h"

program StandAlone_DynAdvCore
   use MAPL_Mod
!  use StandAlone_DynAdvCore_GridCompMod, only: SetServices
   use MPI

   implicit none

!EOP

!EOC

   character(*), parameter :: IAM = __FILE__

   type (MAPL_Cap) :: cap
   type (MAPL_FlapCapOptions) :: cap_options
   integer :: status

   cap_options = MAPL_FlapCapOptions( &
        description = 'FV Standalone DyAdvCore', &
        authors     = 'S.J. Lin, R. Rood, W. Putman')

   !cap = MAPL_Cap('Standalone FV3 DynAdvCore', SetServices, cap_options = cap_options)

   ! comment this out and configure in command line
   cap_options%shared_obj='libFVdycoreCubed_GridComp.so'
   cap_options%proc_name="dynadv_setservices"

   cap = MAPL_Cap('Standalone FV3 DynAdvCore', cap_options = cap_options)
   call cap%run(_RC)

 end Program StandAlone_DynAdvCore

!EOC




