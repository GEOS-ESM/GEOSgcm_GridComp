!===============================================================================
! SVN $Id$
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/branch_tags/ccsm4_0_rel_tags/ccsm4_0_rel_02_share3_100228/shr/shr_const_mod.F90 $
!===============================================================================

MODULE shr_const_mod
   
   use MAPL_ConstantsMod ! use GEOS5 constants

   !----------------------------------------------------------------------------
   ! physical constants (all data public)
   !----------------------------------------------------------------------------
   public

   real,parameter :: SHR_CONST_PI      = MAPL_PI      ! pi
   real,parameter :: SHR_CONST_CDAY    = 86400.0      ! sec in calendar day ~ sec
   real,parameter :: SHR_CONST_G       = MAPL_GRAV    ! acceleration of gravity ~ m/s^2

   real,parameter :: SHR_CONST_RGAS    = MAPL_RUNIV   ! Universal gas constant ~ J/K/kmole
 
   real,parameter :: SHR_CONST_TKFRZ   = MAPL_TICE    ! freezing T of fresh water          ~ K 

   real,parameter :: SHR_CONST_RHOFW   = MAPL_RHOWTR  ! density of fresh water     ~ kg/m^3

END MODULE shr_const_mod
