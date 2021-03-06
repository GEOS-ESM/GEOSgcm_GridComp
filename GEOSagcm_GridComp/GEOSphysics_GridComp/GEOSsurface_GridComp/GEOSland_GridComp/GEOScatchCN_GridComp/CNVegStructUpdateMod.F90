module CNVegStructUpdateMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNVegStructUpdateMod
!
! !DESCRIPTION:
! Module for vegetation structure updates (LAI, SAI, htop, hbot)
!
! !USES:
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public :: CNVegStructUpdate
!
! !REVISION HISTORY:
! 4/23/2004: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNVegStructUpdate
!
! !INTERFACE:
subroutine CNVegStructUpdate(num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, use C state variables and epc to diagnose
! vegetation structure (LAI, SAI, height)
!
! !USES:
   use clmtype
   use pftvarcon    , only: noveg, ncrop, nbrdlf_evr_shrub, nbrdlf_dcd_brl_shrub
   use shr_const_mod, only: SHR_CONST_PI
   use clm_time_manager , only : get_rad_step_size => get_step_size
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilp                 ! number of column soil points in pft filter
   integer, intent(in) :: filter_soilp(:)   ! pft filter for soil points
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 10/28/03: Created by Peter Thornton
! 2/29/08, David Lawrence: revised snow burial fraction for short vegetation
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)        ! pft vegetation type
   integer , pointer :: pcolumn(:)    ! column index associated with each pft
   integer , pointer :: pgridcell(:)  ! pft's gridcell index
   real, pointer :: snowdp(:)     ! snow height (m)
   real, pointer :: leafc(:)      ! (kgC/m2) leaf C
   real, pointer :: deadstemc(:)  ! (kgC/m2) dead stem C
   real, pointer :: woody(:)      !binary flag for woody lifeform (1=woody, 0=not woody)
   real, pointer :: slatop(:)     !specific leaf area at top of canopy, projected area basis [m^2/gC]
   real, pointer :: dsladlai(:)   !dSLA/dLAI, projected area basis [m^2/gC]
   real, pointer :: z0mr(:)       !ratio of momentum roughness length to canopy top height (-)
   real, pointer :: displar(:)    !ratio of displacement height to canopy top height (-)
   real, pointer :: forc_hgt_u_pft(:) ! observational height of wind at pft-level [m]
   real, pointer :: dwood(:)      ! density of wood (kgC/m^3)
!
! local pointers to implicit in/out scalars
!
   integer , pointer :: frac_veg_nosno_alb(:) ! frac of vegetation not covered by snow [-]
   real, pointer :: tlai(:) !one-sided leaf area index, no burying by snow
   real, pointer :: tsai(:) !one-sided stem area index, no burying by snow
   real, pointer :: htop(:) !canopy top (m)
   real, pointer :: hbot(:) !canopy bottom (m)
   real, pointer :: elai(:)     ! one-sided leaf area index with burying by snow
   real, pointer :: esai(:)     ! one-sided stem area index with burying by snow
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: p,c,g        !indices
   integer :: fp           !lake filter indices
   real:: taper        ! ratio of height:radius_breast_height (tree allometry)
   real:: stocking     ! #stems / ha (stocking density)
   real:: ol           ! thickness of canopy layer covered by snow (m)
   real:: fb           ! fraction of canopy layer covered by snow
   real :: tlai_old    ! for use in Zeng tsai formula
   real :: tsai_old    ! for use in Zeng tsai formula
   real :: tsai_min    ! PFT derived minimum tsai
   real :: tsai_alpha  ! monthly decay rate of tsai
   real dt             ! radiation time step (sec)

   real, parameter :: dtsmonth = 2592000. ! number of seconds in a 30 day month (60x60x24x30)
!EOP
!-----------------------------------------------------------------------
! tsai formula from Zeng et. al. 2002, Journal of Climate, p1835
!
! tsai(p) = max( tsai_alpha(ivt(p))*tsai_old + max(tlai_old-tlai(p),0), tsai_min(ivt(p)) )
! notes:
! * RHS tsai & tlai are from previous timestep
! * should create tsai_alpha(ivt(p)) & tsai_min(ivt(p)) in pftvarcon.F90 - slevis
! * all non-crop pfts use same values:
!   crop    tsai_alpha,tsai_min = 0.0,0.1
!   noncrop tsai_alpha,tsai_min = 0.5,1.0  (includes bare soil and urban)
!-------------------------------------------------------------------------------

   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%l%c%p%itype
    pcolumn                        => clm3%g%l%c%p%column
    pgridcell                      => clm3%g%l%c%p%gridcell
    leafc                          => clm3%g%l%c%p%pcs%leafc
    deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
    snowdp                         => clm3%g%l%c%cps%snowdp
    woody                          => pftcon%woody
    slatop                         => pftcon%slatop
    dsladlai                       => pftcon%dsladlai
    z0mr                           => pftcon%z0mr
    displar                        => pftcon%displar
    dwood                          => pftcon%dwood

   ! assign local pointers to derived type arrays (out)
    tlai                           => clm3%g%l%c%p%pps%tlai
    tsai                           => clm3%g%l%c%p%pps%tsai
    htop                           => clm3%g%l%c%p%pps%htop
    hbot                           => clm3%g%l%c%p%pps%hbot
    elai                           => clm3%g%l%c%p%pps%elai
    esai                           => clm3%g%l%c%p%pps%esai
    frac_veg_nosno_alb             => clm3%g%l%c%p%pps%frac_veg_nosno_alb
    forc_hgt_u_pft                 => clm3%g%l%c%p%pps%forc_hgt_u_pft

   dt = real( get_rad_step_size() )

   ! constant allometric parameters
   taper = 200.
   stocking = 1000.

   ! convert from stems/ha -> stems/m^2
   stocking = stocking / 10000.

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)
      g = pgridcell(p)

      if (ivt(p) /= noveg) then

          tlai_old = tlai(p) ! n-1 value
          tsai_old = tsai(p) ! n-1 value

          ! update the leaf area index based on leafC and SLA
          ! Eq 3 from Thornton and Zimmerman, 2007, J Clim, 20, 3902-3923. 
          if (dsladlai(ivt(p)) > 0.) then
             tlai(p) = (slatop(ivt(p))*(exp(leafc(p)*dsladlai(ivt(p))) - 1.))/dsladlai(ivt(p))
          else
             tlai(p) = slatop(ivt(p)) * leafc(p)
          end if
          tlai(p) = max(0., tlai(p))

          ! update the stem area index and height based on LAI, stem mass, and veg type.
          ! With the exception of htop for woody vegetation, this follows the DGVM logic.

          ! tsai formula from Zeng et. al. 2002, Journal of Climate, p1835 (see notes)
          ! Assumes doalb time step .eq. CLM time step, SAI min and monthly decay factor
          ! alpha are set by PFT, and alpha is scaled to CLM time step by multiplying by
          ! dt and dividing by dtsmonth (seconds in average 30 day month)
          ! tsai_min scaled by 0.5 to match MODIS satellite derived values
          if (ivt(p) >= ncrop ) then    ! crops (corn, wheat in CLM)

             tsai_alpha = 1.0-1.0*dt/dtsmonth
             tsai_min = 0.1
          else
             tsai_alpha = 1.0-0.5*dt/dtsmonth
             tsai_min = 1.0
          end if
          tsai_min = tsai_min * 0.5
          tsai(p) = max(tsai_alpha*tsai_old+max(tlai_old-tlai(p),0.),tsai_min)

          if (woody(ivt(p)) == 1.) then

             ! trees and shrubs

             ! if shrubs have a squat taper 
             if (ivt(p) >= nbrdlf_evr_shrub .and. ivt(p) <= nbrdlf_dcd_brl_shrub) then
                taper = 10.
             ! otherwise have a tall taper
             else
                taper = 200.
             end if

             ! trees and shrubs for now have a very simple allometry, with hard-wired
             ! stem taper (height:radius) and hard-wired stocking density (#individuals/area)
                htop(p) = ((3. * deadstemc(p) * taper * taper)/ &
                  (SHR_CONST_PI * stocking * dwood(ivt(p))))**(1./3.)

             ! Peter Thornton, 5/3/2004
             ! Adding test to keep htop from getting too close to forcing height for windspeed
             ! Also added for grass, below, although it is not likely to ever be an issue.
             htop(p) = min(htop(p),(forc_hgt_u_pft(p)/(displar(ivt(p))+z0mr(ivt(p))))-3.)

             ! Peter Thornton, 8/11/2004
             ! Adding constraint to keep htop from going to 0.0.
             ! This becomes an issue when fire mortality is pushing deadstemc
             ! to 0.0.
             htop(p) = max(htop(p), 0.01)

             hbot(p) = max(0., min(3., htop(p)-1.))

          else
             ! grasses

             ! height for grasses depends only on LAI
             htop(p) = max(0.25, tlai(p) * 0.25)

             htop(p) = min(htop(p),(forc_hgt_u_pft(p)/(displar(ivt(p))+z0mr(ivt(p))))-3.)

             ! Peter Thornton, 8/11/2004
             ! Adding constraint to keep htop from going to 0.0.
             htop(p) = max(htop(p), 0.01)

             hbot(p) = max(0.0, min(0.05, htop(p)-0.20))
          end if

      else
          tlai(p) = 0.
          tsai(p) = 0.
          htop(p) = 0.
          hbot(p) = 0.
      end if
      
      ! adjust lai and sai for burying by snow. 

      ! snow burial fraction for short vegetation (e.g. grasses) as in
      ! Wang and Zeng, 2007.
      if (ivt(p) > noveg .and. ivt(p) <= nbrdlf_dcd_brl_shrub ) then
         ol = min( max(snowdp(c)-hbot(p), 0.), htop(p)-hbot(p))
         fb = 1. - ol / max(1.e-06, htop(p)-hbot(p))
      else
         fb = 1. - max(min(snowdp(c),0.2),0.)/0.2   ! 0.2m is assumed
              !depth of snow required for complete burial of grasses
      endif

      elai(p) = max(tlai(p)*fb, 0.0)
      esai(p) = max(tsai(p)*fb, 0.0)

      ! Fraction of vegetation free of snow
      if ((elai(p) + esai(p)) > 0.) then
         frac_veg_nosno_alb(p) = 1
      else
         frac_veg_nosno_alb(p) = 0
      end if

   end do

end subroutine CNVegStructUpdate
!-----------------------------------------------------------------------

end module CNVegStructUpdateMod
