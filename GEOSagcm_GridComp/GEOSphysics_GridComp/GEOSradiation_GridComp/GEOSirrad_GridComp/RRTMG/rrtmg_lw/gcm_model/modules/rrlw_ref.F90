#include "_gpudef.inc"    
     
      module rrlw_ref

       !use parkind, only : im => kind , rb => kind 

      implicit none
      

!------------------------------------------------------------------
! rrtmg_lw reference atmosphere 
! Based on standard mid-latitude summer profile
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! pref   :  real   : Reference pressure levels
! preflog:  real   : Reference pressure levels, ln(pref)
! tref   :  real   : Reference temperature levels for MLS profile
! chi_mls:  real   : 
!------------------------------------------------------------------

      real , dimension(59) :: pref
      real , dimension(59) :: preflog
      real , dimension(59) :: tref
      real :: chi_mls(7,59)

      ! (dmb 2012) These GPU arrays are defined as constant so that they are cached.
      ! This is really needed because they accessed in quite a scattered pattern.
      real _gpucon :: chi_mlsd(7,59)
      real _gpucon :: preflogd(59)
      real _gpucon :: trefd(59)

      contains
      ! (dmb 2012) Copy the reference arrays over to the GPU
      subroutine copyToGPUref()

        chi_mlsd = chi_mls
        preflogd = preflog
        trefd = tref

      end subroutine 

      end module rrlw_ref
