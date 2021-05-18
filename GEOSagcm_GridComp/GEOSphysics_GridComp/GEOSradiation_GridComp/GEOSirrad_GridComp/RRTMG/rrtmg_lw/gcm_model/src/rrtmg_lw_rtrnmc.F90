!pmn: recombine?
! (dmb 2012) This is the GPU version of the rtrnmc subroutine. This has been greatly
! modified to be efficiently run on the GPU. Originally, there was a g-point loop within
! this subroutine to perform the summation of the fluxes over the g-points. This has been
! modified so that this subroutine can be run in parallel across the g-points. This was
! absolutely critical because of two reasons.
! 1. For a relatively low number of profiles, there wouldn't be enough threads to keep
!    the GPU busy enough to run at full potential. As a result of this, this subroutine
!    would end up being a bottleneck.
! 2. The memory access for the GPU arrays would be innefient because there would be very
!    little coalescing which is critical for obtaining optimal performance.

module rrtmg_lw_rtrnmc

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! --------- Modules ----------

   use parrrtm, only : nbndlw, ngptlw
   use rrlw_tbl, only: bpade, tblint, tau_tbl, exp_tbl, tfn_tbl
   use rrtmg_lw_setcoef, only : pwvcm, planklay, planklev, plankbnd, dplankbnd_dt
   use rrlw_con, only: fluxfac
   use rrlw_wvn, only: delwave

   implicit none 
      
contains

   !-----------------------------------------------------
   subroutine rtrnmc (ncol, nlay, idrv, ngb, &
      semiss, taug, pfracs, cloudy, cldfmc, taucmc, &
      totuflux, totdflux, totuclfl, totdclfl, &
      dtotuflux_dt, dtotuclfl_dt, &
      olrb06, olrb09, olrb10, olrb11, &
      dolrb06_dt, dolrb09_dt, dolrb10_dt, dolrb11_dt)
   !-----------------------------------------------------
   !
   !  Original version:   E. J. Mlawer, et al. RRTM_V3.0
   !  Revision for GCMs:  Michael J. Iacono; October, 2002
   !  Revision for F90:  Michael J. Iacono; June, 2006
   !  Revision for dFdT option: M. J. Iacono and E. J. Mlawer, Nov 2009
   !
   !  This program calculates the upward fluxes and downward fluxes
   !  for an arbitrary clear or cloudy atmosphere. The input is the
   !  atmospheric profile, all Planck function information, and the
   !  McICA cloud mask. A variable diffusivity angle (secdiff)
   !  is used for the angle integration. Bands 2-3 and 5-9 use a
   !  value for SECDIFF that varies from 1.50 to 1.80 as a function 
   !  of the column water vapor, and other bands use a value of 1.66.
   !  The Gaussian weight appropriate to this angle (wtdiff=0.5) is
   !  applied here. Note that use of the emissivity angle for the
   !  flux integration can cause errors of 1 to 4 W/m2 within cloudy
   !  layers. Clouds are treated with the McICA stochastic approach.
   !
   !  Also provides the optional capability to calculate the deriv-
   !  ative of upward flux respect to surface temperature using the
   !  pre-tabulated derivative of the Planck function with respect 
   !  to temperature integrated over each spectral band.
   !-------------------------------------------------------

      integer, intent(in) :: ncol  ! number of columns
      integer, intent(in) :: nlay  ! number of layers
      integer, intent(in) :: idrv  ! do d/dTsurf calcs if == 1

      integer, intent(in) :: ngb(ngptlw)  ! band index of g-point

      real,    intent(in) :: semiss  (     nbndlw,ncol)  ! surface emissivity
      real,    intent(in) :: taug    (nlay,ngptlw,ncol)  ! gas optical depth
      real,    intent(in) :: pfracs  (nlay,ngptlw,ncol)  ! Planck fractions
      logical, intent(in) :: cloudy  (nlay,       ncol)  ! cloudy for ANY g-point
      real,    intent(in) :: cldfmc  (nlay,ngptlw,ncol)  ! cloud fraction
      real,    intent(in) :: taucmc  (nlay,ngptlw,ncol)  ! cloud optical thickness
     
      ! spectrally summed fluxes and upward flux derivatives wrt Tsurf
      real, intent(out) :: totuflux     (0:nlay,ncol)  ! upward longwave flux (W/m2)
      real, intent(out) :: totdflux     (0:nlay,ncol)  ! downward longwave flux (W/m2)
      real, intent(out) :: totuclfl     (0:nlay,ncol)  ! clrsky upward lw flux (W/m2)
      real, intent(out) :: totdclfl     (0:nlay,ncol)  ! clrsky downward lw flux (W/m2)
      real, intent(out) :: dtotuflux_dt (0:nlay,ncol)  ! d/d(Tsurf) (W/m2/K)
      real, intent(out) :: dtotuclfl_dt (0:nlay,ncol)  ! d/d(Tsurf) (W/m2/K)

      ! TOA OLR in bands 6 & 9-11 and their derivatives wrt Tsurf
      real, intent(out) :: olrb06     (ncol)  ! (W/m2)
      real, intent(out) :: olrb09     (ncol)  ! (W/m2)
      real, intent(out) :: olrb10     (ncol)  ! (W/m2)
      real, intent(out) :: olrb11     (ncol)  ! (W/m2)
      real, intent(out) :: dolrb06_dt (ncol)  ! (W/m2/K)
      real, intent(out) :: dolrb09_dt (ncol)  ! (W/m2/K)
      real, intent(out) :: dolrb10_dt (ncol)  ! (W/m2/K)
      real, intent(out) :: dolrb11_dt (ncol)  ! (W/m2/K)

      ! ----- Local -----
   
      ! g-point fluxes
      real, dimension (0:nlay,ngptlw,ncol) :: &
         gurad,          & ! upward longwave flux (W/m2)
         gdrad,          & ! downward longwave flux (W/m2)
         gclrurad,       & ! clear sky upward longwave flux (W/m2)
         gclrdrad,       & ! clear sky downward longwave flux (W/m2)
         gdtotuflux_dt,  & ! d(upward longwave flux)/d(Tsurf) (W/m2/K)
         gdtotuclfl_dt     ! d(clrsky upward lw flux)/d(Tsurf) (W/m2/K)
  
!pmn: better for all _dt -> _dTs so not confused with time

      real :: agas(nlay)    ! gas absorptivity
      real :: atot(nlay)    ! gas and cloud absorptivity
      real :: bbugas(nlay)  ! gas Planck function for upward rt
      real :: bbutot(nlay)  ! gas and cloud Planck function for upward calc.
     
      real :: secdiff    ! diffusivity angle
      real :: odepth     ! gas optical depth
      real :: odclds     ! cloud optical depth
      real :: absclds    ! cloud absorption
      real :: efclfracs  ! effective cloud fraction
      real :: odtot      ! optical depth of gas and cloud
      real :: tfacgas    ! gas pade factor, used for planck fn
      real :: tfactot    ! gas and cloud pade factor, used for planck fn
      real :: tausfac    ! general pade factor, used for planck fn
      real :: bbdgas     ! gas planck function for downward rt
      real :: gassrc     ! source radiance due to gas only
      real :: bbdtot     ! gas and cloud planck function for downward rt
      real :: tblind     ! real lookup table index
      real :: radld      ! downward radiance
      real :: radclrd    ! downward radiance for clear column
      real :: radlu      ! upward radiance
      real :: radclru    ! upward radiance for clear column
      real :: rad0       ! surface emitted radiance
      real :: reflect    ! surface reflectance

      real :: plfrac, blay, dplankup, dplankdn

      ! derivatives with respect to surface temperature
      real :: d_rad0_dt, d_radlu_dt, d_radclru_dt

      integer :: icol          ! column index
      integer :: lev           ! level index
      integer :: ig            ! g-point index
      integer :: ibnd          ! band index
      integer :: iclddn        ! cloud reached flag in down path
      integer :: ittot, itgas  ! lookup table indices
   
      ! This weight corresponds to the standard sec(diffusivity angle) of 1.66.
      real, parameter :: wtdiff = 0.5  ! weight for radiance to flux conversion

      ! diffusivity factors
      real, parameter, dimension(16) :: a0 = &
         [1.66 , 1.55 , 1.58 , 1.66 , &
          1.54 , 1.454, 1.89 , 1.33 , &
          1.668, 1.66 , 1.66 , 1.66 , &
          1.66 , 1.66 , 1.66 , 1.66 ]
      real, parameter, dimension(16) :: a1 = &
         [0.00 , 0.25 , 0.22 , 0.00 , &
          0.13 , 0.446,-0.10 , 0.40 , &
         -0.006, 0.00 , 0.00 , 0.00 , &
          0.00 , 0.00 , 0.00 , 0.00 ]
      real, parameter, dimension(16) :: a2 = &
         [0.00 ,-12.0 ,-11.7 , 0.00 , &
         -0.72 ,-0.243, 0.19 ,-0.062, &
          0.414, 0.00 , 0.00 , 0.00 , &
          0.00 , 0.00 , 0.00 , 0.00 ]

      ! zero g-point fluxes
! pmn: seems that only top down needs to be zeroed and that all the 
! pmn: += assignments below could be simply assignments
      gurad = 0. 
      gdrad = 0. 
      gclrurad = 0. 
      gclrdrad = 0. 
      if (idrv == 1) then
         gdtotuflux_dt = 0. 
         gdtotuclfl_dt = 0. 
      endif

! pmn: check new ordering is better ie, icol last ilay first
      ! main loop
      ! g-points are "monochromatic" in CKD
      do icol = 1,ncol
      do ig = 1,ngptlw 

         ibnd = ngb(ig)
      
         ! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
         ! and 1.80) as a function of total column water vapor. The function
         ! has been defined to minimize flux and cooling rate errors in these
         ! bands over a wide range of precipitable water values.

         if (ibnd == 1 .or. ibnd == 4 .or. ibnd >= 10) then
            secdiff = 1.66 
         else
            secdiff = a0(ibnd) + a1(ibnd) * exp(a2(ibnd)*pwvcm(icol))
            if (secdiff > 1.80 ) secdiff = 1.80 
            if (secdiff < 1.50 ) secdiff = 1.50 
         endif
    
         ! downward LW at TOA is zero
         radld = 0. 
         radclrd = 0. 

         ! column clear until cloud layer found
         iclddn = 0

         ! Downward radiative transfer loop
         do lev = nlay,1,-1

            plfrac = pfracs(lev,ig,icol)
            blay = planklay(ibnd,lev,icol)
            dplankup = planklev(ibnd,lev,  icol) - blay
            dplankdn = planklev(ibnd,lev-1,icol) - blay

            odepth = secdiff * taug(lev,ig,icol)
            if (odepth .lt. 0.) odepth = 0. 

            ! cloudy for any g-point
            if (cloudy(lev,icol)) then

               iclddn = 1  ! flag any cloud in column

               odclds = secdiff * taucmc(lev,ig,icol)
               absclds = 1. - exp(-odclds)
               efclfracs = absclds * cldfmc(lev,ig,icol)

               ! table lookup based on optical depth:
               ! bpade    1/(pade constant)
               ! tau_tbl  clear sky optical depth look-up table
               !          (used in cloudy radiative transfer)
               ! exp_tbl  exponential look-up table for transmittance
               ! tfn_tbl  tau transition function look-up table
               !          (i.e. the transition of the Planck function
               !          from that for the mean layer temperature to
               !          that for the layer boundary temperature as
               !          a function of optical depth. The "linear in
               !          tau" method is used to make the table.)

               ! gas only
               tblind = odepth/(bpade+odepth)
               itgas = tblint * tblind + 0.5 
               agas(lev) = 1. - exp_tbl(itgas)
               tfacgas = tfn_tbl(itgas)
               odepth = tau_tbl(itgas)

               ! gas and cloud
               odtot = odepth + odclds
               tblind = odtot/(bpade+odtot)
               ittot = tblint * tblind + 0.5 
               atot(lev) = 1. - exp_tbl(ittot)
               tfactot = tfn_tbl(ittot)

               bbdgas      = plfrac * (blay + tfacgas * dplankdn)
               bbugas(lev) = plfrac * (blay + tfacgas * dplankup)
               bbdtot      = plfrac * (blay + tfactot * dplankdn)
               bbutot(lev) = plfrac * (blay + tfactot * dplankup)

               ! Update of downward LW, radld, through a mixed
               ! gas and cloud layer. See simpler clear (gas
               ! only) case in else-branch first.
               !    If f is cloud fraction, then effective
               ! layer absorbtivity is
               !   alayer = (1-f) * agas + f * aboth
               ! where
               !   aboth = 1 - Tboth = 1 - Tgas * Tcld
               !     = 1 - (1-agas) * (1-acld)
               !     = agas + acld - agas * acld,
               ! so
               !   alayer = agas + f * acld * (1-agas)
               !          = agas + efclfracs * (1-agas),
               ! hence attenuation term, -radld * alayer, below.
               !   For source terms,
               ! slayer = f * bbdtot * atot + (1-f) * bbdgas * agas
               !        = f * bbdtot * atot + (1-f) * gassrc
               !        = gassrc + f * (bbdtot * atot - gassrc),
               ! as below.

               gassrc = agas(lev) * bbdgas
               radld = radld - radld * (agas(lev) + &
                  efclfracs * (1. - agas(lev))) + &
                  gassrc + cldfmc(lev,ig,icol) * &
                  (bbdtot * atot(lev) - gassrc)
               gdrad(lev-1,ig,icol) = gdrad(lev-1,ig,icol) + radld

            else  ! clear layer

               tblind = odepth/(bpade+odepth)
               itgas = tblint * tblind + 0.5 
               agas(lev) = 1. - exp_tbl(itgas)
               tausfac = tfn_tbl(itgas)

               bbdgas      = plfrac * (blay + tausfac * dplankdn)
               bbugas(lev) = plfrac * (blay + tausfac * dplankup)

               ! Update of downward LW, radld, through a
               ! layer of transmissivity T:
               !   radld -> radld * T + bbd * (1-T)
               !         -> radld + (bbd - radld) * (1-T)
               ! where bbd is effective Planck func for layer.

               radld = radld + (bbdgas - radld) * agas(lev)
               gdrad(lev-1,ig,icol) = gdrad(lev-1,ig,icol) + radld

            endif

            ! Set clear sky stream to total sky stream as long as layers
            ! remain clear. Streams diverge when a cloud is reached (iclddn=1),
            ! and clear sky stream must be computed separately from that point.

            if (iclddn == 1) then
               radclrd = radclrd + (bbdgas-radclrd) * agas(lev) 
               gclrdrad(lev-1,ig,icol) = gclrdrad(lev-1,ig,icol) + radclrd
            else
               radclrd = radld
               gclrdrad(lev-1,ig,icol) = gdrad(lev-1,ig,icol)
            endif

         enddo  ! downward loop

         ! Include the contribution of spectrally varying longwave emissivity
         ! and reflection from the surface to the upward radiative transfer.
         ! Note: Spectral and Lambertian reflection are identical for the
         ! diffusivity angle flux integration used here.
         ! Note: The emissivity is applied to plankbnd and dplankbnd_dt when 
         ! they are defined in subroutine setcoef. 
    
         ! surface emission
         rad0 = pfracs(1,ig,icol) * plankbnd(ibnd,icol)
         if (idrv == 1) then
            d_rad0_dt = pfracs(1,ig,icol) * dplankbnd_dt(ibnd,icol)
         endif

         ! Add in specular reflection of surface downward radiance
         reflect = 1. - semiss(ibnd,icol)
         radlu   = rad0 + reflect * radld
         radclru = rad0 + reflect * radclrd
         gurad   (0,ig,icol) = gurad   (0,ig,icol) + radlu
         gclrurad(0,ig,icol) = gclrurad(0,ig,icol) + radclru
         if (idrv == 1) then
            d_radlu_dt   = d_rad0_dt
            d_radclru_dt = d_rad0_dt
            gdtotuflux_dt(0,ig,icol) = d_radlu_dt
            gdtotuclfl_dt(0,ig,icol) = d_radclru_dt
         endif

         ! Upward radiative transfer loop
         do lev = 1,nlay

            ! cloudy for any g-point
            if (cloudy(lev,icol)) then

               gassrc = bbugas(lev) * agas(lev)
               odclds = secdiff * taucmc(lev,ig,icol)
               absclds = 1. - exp(-odclds)
               efclfracs = absclds * cldfmc(lev,ig,icol)

! pmn: can we eliminate some repeated calcs from down loop? at least for clds?

               radlu = radlu - radlu * (agas(lev) + &
                   efclfracs * (1. - agas(lev))) + &
                   gassrc + cldfmc(lev,ig,icol) * &
                   (bbutot(lev) * atot(lev) - gassrc)
               gurad(lev,ig,icol) = gurad(lev,ig,icol) + radlu
               if (idrv == 1) then
                  d_radlu_dt = d_radlu_dt * cldfmc(lev,ig,icol)        * (1. - atot(lev)) + &
                               d_radlu_dt * (1. - cldfmc(lev,ig,icol)) * (1. - agas(lev))
! pmn: -----
! pmn: this makes sense as the sum of the transmitted sensitivities for the cloudy
! pmn: and clear portions, but why use atot here vs. aboth for the non-derivative?
! pmn: -----
                  gdtotuflux_dt(lev,ig,icol) = gdtotuflux_dt(lev,ig,icol) + d_radlu_dt
               endif

            else  ! clear layer

               radlu = radlu + (bbugas(lev) - radlu) * agas(lev)
               gurad(lev,ig,icol) = gurad(lev,ig,icol) + radlu
               if (idrv == 1) then
                  d_radlu_dt = d_radlu_dt * (1. - agas(lev))
                  gdtotuflux_dt(lev,ig,icol) = gdtotuflux_dt(lev,ig,icol) + d_radlu_dt
               endif

            endif

            ! Set clear sky stream to total sky stream as long as all layers
            ! are clear (iclddn=0). Streams must be calculated separately at 
            ! all layers when a cloud is present (iclddn=1), because surface 
            ! reflectance is different for each stream.

            if (iclddn == 1) then
               radclru = radclru + (bbugas(lev) - radclru) * agas(lev) 
               gclrurad(lev,ig,icol) = gclrurad(lev,ig,icol) + radclru
            else
               radclru = radlu
               gclrurad(lev,ig,icol) = gurad(lev,ig,icol)
            endif
            if (idrv == 1) then
               if (iclddn == 1) then
                  d_radclru_dt = d_radclru_dt * (1. - agas(lev))
                  gdtotuclfl_dt(lev,ig,icol) = gdtotuclfl_dt(lev,ig,icol) + d_radclru_dt
               else
                  d_radclru_dt = d_radlu_dt
                  gdtotuclfl_dt(lev,ig,icol) = gdtotuflux_dt(lev,ig,icol)
               endif
            endif

         enddo  ! layer
          
         ! modify g-points values so later summation is simpler 
         tblind = wtdiff * delwave(ibnd) * fluxfac
         do lev = 0,nlay  
           gurad   (lev,ig,icol) = gurad   (lev,ig,icol) * tblind
           gdrad   (lev,ig,icol) = gdrad   (lev,ig,icol) * tblind
           gclrurad(lev,ig,icol) = gclrurad(lev,ig,icol) * tblind
           gclrdrad(lev,ig,icol) = gclrdrad(lev,ig,icol) * tblind
         end do
         if (idrv == 1) then
           do lev = 0,nlay
             gdtotuflux_dt(lev,ig,icol) = gdtotuflux_dt(lev,ig,icol) * tblind
             gdtotuclfl_dt(lev,ig,icol) = gdtotuclfl_dt(lev,ig,icol) * tblind
           end do
         endif

      end do  ! g-point
      end do  ! column

      ! zero spectral sum accumulators
      totuflux = 0.
      totdflux = 0.
      totuclfl = 0.
      totdclfl = 0.
      olrb06   = 0.
      olrb09   = 0.
      olrb10   = 0.
      olrb11   = 0.
      if (idrv == 1) then
         dtotuflux_dt = 0.
         dtotuclfl_dt = 0.
         dolrb06_dt   = 0.
         dolrb09_dt   = 0.
         dolrb10_dt   = 0.
         dolrb11_dt   = 0.
      end if

      ! Adds up the indivial g-point fluxes to arrive at a final
      ! upward and downward flux value for each column and layer.
! pmn test different loop order

      do icol = 1,ncol
      do lev = 0,nlay
         do ig = 1,ngptlw
            totuflux(lev,icol) = totuflux(lev,icol) + gurad   (lev,ig,icol)
            totdflux(lev,icol) = totdflux(lev,icol) + gdrad   (lev,ig,icol)
            totuclfl(lev,icol) = totuclfl(lev,icol) + gclrurad(lev,ig,icol)
            totdclfl(lev,icol) = totdclfl(lev,icol) + gclrdrad(lev,ig,icol)
         end do
      end do
      end do

! pmn: seperate index for ig's in an band?
      do icol = 1,ncol
         do ig = 1,ngptlw
            if (ngb(ig) ==  6) olrb06(icol) = olrb06(icol) + gurad(nlay,ig,icol)
            if (ngb(ig) ==  9) olrb09(icol) = olrb09(icol) + gurad(nlay,ig,icol)
            if (ngb(ig) == 10) olrb10(icol) = olrb10(icol) + gurad(nlay,ig,icol)
            if (ngb(ig) == 11) olrb11(icol) = olrb11(icol) + gurad(nlay,ig,icol)
         end do
      end do

      if (idrv == 1) then

         do icol = 1,ncol
         do lev = 0,nlay
            do ig = 1,ngptlw
               dtotuflux_dt(lev,icol) = dtotuflux_dt(lev,icol) + gdtotuflux_dt(lev,ig,icol)
               dtotuclfl_dt(lev,icol) = dtotuclfl_dt(lev,icol) + gdtotuclfl_dt(lev,ig,icol)
            end do
         end do
         end do

! pmn: seperate index for ig's in an band?
         do icol = 1,ncol
            do ig = 1,ngptlw
               if (ngb(ig) ==  6) dolrb06_dt(icol) = dolrb06_dt(icol) + gdtotuflux_dt(nlay,ig,icol)
               if (ngb(ig) ==  9) dolrb09_dt(icol) = dolrb09_dt(icol) + gdtotuflux_dt(nlay,ig,icol)
               if (ngb(ig) == 10) dolrb10_dt(icol) = dolrb10_dt(icol) + gdtotuflux_dt(nlay,ig,icol)
               if (ngb(ig) == 11) dolrb11_dt(icol) = dolrb11_dt(icol) + gdtotuflux_dt(nlay,ig,icol)
            end do
         end do

      end if

   end subroutine rtrnmc

end module rrtmg_lw_rtrnmc
