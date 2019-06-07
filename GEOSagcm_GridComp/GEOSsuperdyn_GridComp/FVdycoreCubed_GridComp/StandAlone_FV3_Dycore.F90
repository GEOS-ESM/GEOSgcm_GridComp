#define I_AM_MAIN

#include "MAPL_Generic.h"
#define _RC rc=status); _VERIFY(status

program StandAlone_FV3_Dycore
   use MAPL_Mod
   use FVdycoreCubed_GridComp,      only: SetServices
   use MAPL_CapOptionsMod
   use MAPL_FlapCapOptionsMod
   implicit none

!EOP

!EOC

   character(*), parameter :: IAM = __FILE__

   type (MAPL_Cap) :: cap
   type (MAPL_FlapCapOptions), allocatable :: cap_options
   integer :: status

   cap_options = MAPL_FlapCapOptions(description = 'FV Standalone Dycore',&
                                     authors      =  'S.J. Lin, R. Rood, W. Putman')
   cap = MAPL_Cap('GCM', SetServices, cap_options = cap_options)
   call cap%run(_RC)

 end Program StandAlone_FV3_Dycore

!EOC

