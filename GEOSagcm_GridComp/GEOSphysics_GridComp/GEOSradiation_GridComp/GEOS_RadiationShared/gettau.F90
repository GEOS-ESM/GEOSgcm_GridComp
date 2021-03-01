! $Id$
! $Name$

!=============================================================================
!BOP

! !MODULE: gettau -- A module to calculate cloud optical depth

! !INTERFACE:

module gettau

! !USES:

   use MAPL_ConstantsMod, only: MAPL_GRAV

   implicit none
   private

! !PUBLIC ROUTINES:

   public :: getvistau, getnirtau, getirtau

!EOP

contains

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: getvistau
!
! !INTERFACE:
   subroutine getvistau(nlevs,cosz,dp,fcld,reff,hydromets,ict,icb,&
                        taubeam,taudiff,asycl)

! !USES:
      use rad_constants, only: caib, caif, &
                               aib_uv, awb_uv, arb_uv, &
                               aig_uv, awg_uv, arg_uv

      implicit none

! !INPUT PARAMETERS:
      integer, intent(IN ) :: nlevs          !  Number of levels
      real,    intent(IN ) :: cosz           !  Cosine of solar zenith angle
      real,    intent(IN ) :: dp(:)          !  Delta pressure (Pa)
      real,    intent(IN ) :: fcld(:)        !  Cloud fraction (used sometimes)
      real,    intent(IN ) :: reff(:,:)      !  Effective radius (microns)
      real,    intent(IN ) :: hydromets(:,:) !  Hydrometeors (kg/kg)
      integer, intent(IN ) :: ict, icb       !  Flags for various uses 
!                 ict  = 0   Indicates that in-cloud values have been given
!                            and are expected
!                 ict != 0   Indicates that overlap computation is needed, and:
!                               ict is the level of the mid-high boundary
!                               icb is the level of the low-mid  boundary
!                
! !OUTPUT PARAMETERS:
      real,    intent(OUT) :: taubeam(:,:)   !  Optical Depth for Beam Radiation
      real,    intent(OUT) :: taudiff(:,:)   !  Optical Depth for Diffuse Radiation
      real,    intent(OUT) ::   asycl(:  )   !  Cloud Asymmetry Factor
! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for visible wavelengths
!  In general will compute in-cloud - will do grid mean when called
!  for diagnostic use only. ict flag will indicate which to do.
!  Slots for reff, hydrometeors, taubeam, taudiff, and asycl are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
!
!  In the below calculations, the constants used in the tau calculation are in 
!  m$^2$ g$^{-1}$ and m$^2$ g$^{-1}$ $\mu$m. Thus, we must convert the kg contained in the 
!  pressure (Pa = kg m$^{-1}$ s$^{-2}$) to grams.
!
! !REVISION HISTORY: 
!    2011.10.27   Molod moved to Radiation_Shared and revised arg list, units
!    2011.11.16   MAT: Generalized to a call that is per-column
!
!EOP
!------------------------------------------------------------------------------
!BOC

      integer            :: k,in,im,it,ia,kk
      real               :: fm,ft,fa,xai,tauc,asyclt
      real               :: cc(3)
      real               :: taucld1,taucld2,taucld3,taucld4
      real               :: g1,g2,g3,g4

      real               :: reff_snow

#include "getvistau.code"

      return

!EOC
   end subroutine getvistau
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: getnirtau
!
! !INTERFACE:
   subroutine getnirtau(ib,nlevs,cosz,dp,fcld,reff,hydromets,ict,icb,&
                        taubeam,taudiff,asycl,ssacl)

! !USES:
      use rad_constants, only: caib, caif, &
                               aib_nir, awb_nir, arb_nir, &
                               aia_nir, awa_nir, ara_nir, &
                               aig_nir, awg_nir, arg_nir
      implicit none

! !INPUT PARAMETERS:
      integer, intent(IN ) :: ib             !  Band number
      integer, intent(IN ) :: nlevs          !  Number of levels
      real,    intent(IN ) :: cosz           !  Cosine of solar zenith angle
      real,    intent(IN ) :: dp(:)          !  Delta pressure in Pa
      real,    intent(IN ) :: fcld(:)        !  Cloud fraction (used sometimes)
      real,    intent(IN ) :: reff(:,:)      !  Effective radius (microns)
      real,    intent(IN ) :: hydromets(:,:) !  Hydrometeors (kg/kg)
      integer, intent(IN ) :: ict, icb           !  Flags for various uses 
!                 ict  = 0   Indicates that in-cloud values have been given
!                            and are expected
!                 ict != 0   Indicates that overlap computation is needed, and:
!                               ict is the level of the mid-high boundary
!                               icb is the level of the low-mid  boundary
!                
! !OUTPUT PARAMETERS:
      real,    intent(OUT) :: taubeam(:,:)   !  Optical depth for beam radiation
      real,    intent(OUT) :: taudiff(:,:)   !  Optical depth for diffuse radiation
      real,    intent(OUT) ::   ssacl(:  )   !  Cloud single scattering albedo
      real,    intent(OUT) ::   asycl(:  )   !  Cloud asymmetry factor
! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for near-infrared wavelengths
!  In general will compute in-cloud - will do grid mean when called
!  for diagnostic use only. ict flag will indicate which to do.
!  Slots for reff, hydrometeors and tauall are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
!
!  In the below calculations, the constants used in the tau calculation are in 
!  m$^2$ g$^{-1}$ and m$^2$ g$^{-1}$ $\mu$m. Thus, we must convert the kg contained in the 
!  pressure (Pa = kg m$^{-1}$ s$^{-2}$) to grams.
!
! !REVISION HISTORY: 
!    2011.10.27   Molod moved to Radiation_Shared and revised arg list, units
!    2011.11.16   MAT: Generalized to a call that is per-column
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer            :: k,in,im,it,ia,kk
      real               :: fm,ft,fa,xai,tauc,asyclt,ssaclt
      real               :: cc(3)
      real               :: taucld1,taucld2,taucld3,taucld4
      real               :: g1,g2,g3,g4
      real               :: w1,w2,w3,w4

      real               :: reff_snow

#include "getnirtau.code"

      return

!EOC
   end subroutine getnirtau
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: getirtau
!
! !INTERFACE:
   subroutine getirtau(ib,nlevs,dp,fcld,reff,hydromets,&
                       taudiag,tcldlyr,enn)

! !USES:
      use rad_constants, only: aib_ir, awb_ir, &
                               aiw_ir, aww_ir, &
                               aig_ir, awg_ir

      implicit none

! !INPUT PARAMETERS:
      integer, intent(IN ) :: ib             !  Band number
      integer, intent(IN ) :: nlevs          !  Number of levels
      real,    intent(IN ) :: dp(:)          !  Delta pressure in Pa
      real,    intent(IN ) :: fcld(:)        !  Cloud fraction (used sometimes)
      real,    intent(IN ) :: reff(:,:)      !  Effective radius (microns)
      real,    intent(IN ) :: hydromets(:,:) !  Hydrometeors (kg/kg)

! !OUTPUT PARAMETERS:
      real,    intent(OUT) :: taudiag( :,:)  !  Optical depth for beam radiation
      real,    intent(OUT) :: tcldlyr(0:  )  !  Flux transmissivity?
      real,    intent(OUT) ::     enn(0:  )  !  Flux transmissivity of a cloud layer?

! !DESCRIPTION:
!  Compute in-cloud or grid mean optical depths for infrared wavelengths
!  Slots for reff, hydrometeors and tauall are as follows:
!                 1         Cloud Ice
!                 2         Cloud Liquid
!                 3         Falling Liquid (Rain)
!                 4         Falling Ice (Snow)
!
!  In the below calculations, the constants used in the tau calculation are in 
!  m$^2$ g$^{-1}$ and m$^2$ g$^{-1}$ $\mu$m. Thus, we must convert the kg contained in the 
!  pressure (Pa = kg m$^{-1}$ s$^{-2}$) to grams.
!
! !REVISION HISTORY: 
!    2011.11.18   MAT moved to Radiation_Shared and revised arg list, units
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer            :: k
      real               :: taucld1,taucld2,taucld3,taucld4
      real               :: g1,g2,g3,g4,gg
      real               :: w1,w2,w3,w4,ww
      real               :: ff,tauc

      real               :: reff_snow

#include "getirtau.code"

   end subroutine getirtau

end module gettau

