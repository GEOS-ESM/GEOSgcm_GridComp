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
   use rrtmg_lw_setcoef, only : &
     pwvcm, planklay, planklev, plankbnd, dplankbnd_dTs
   use rrlw_con, only: fluxfac
   use rrlw_wvn, only: ngb, delwave

   implicit none 
      
contains

   !-----------------------------------------------------
   subroutine rtrnmc (ncol, nlay, dudTs, &
      semiss, taug, pfracs, cloudy, cldfmc, taucmc, &
      totuflux, totdflux, totuclfl, totdclfl, &
      dtotuflux_dTs, dtotuclfl_dTs, &
      band_output, olrb, dolrb_dTs)
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

      integer, intent(in) :: ncol   ! number of columns
      integer, intent(in) :: nlay   ! number of layers
      logical, intent(in) :: dudTs  ! do d(upflux)/d(Tsurf) calcs

      real,    intent(in) :: semiss  (     nbndlw,ncol)  ! surface emissivity
      real,    intent(in) :: taug    (nlay,ngptlw,ncol)  ! gas optical depth
      real,    intent(in) :: pfracs  (nlay,ngptlw,ncol)  ! Planck fractions
      logical, intent(in) :: cloudy  (nlay,       ncol)  ! cloudy for ANY g-point
      real,    intent(in) :: cldfmc  (nlay,ngptlw,ncol)  ! cloud fraction
      real,    intent(in) :: taucmc  (nlay,ngptlw,ncol)  ! cloud optical thickness
     
      ! spectrally summed fluxes and upward flux derivatives wrt Tsurf
      real, intent(out) :: totuflux      (0:nlay,ncol)  ! upward longwave flux [W/m2]
      real, intent(out) :: totdflux      (0:nlay,ncol)  ! downward longwave flux [W/m2]
      real, intent(out) :: totuclfl      (0:nlay,ncol)  ! clrsky upward lw flux [W/m2]
      real, intent(out) :: totdclfl      (0:nlay,ncol)  ! clrsky downward lw flux [W/m2]
      real, intent(out) :: dtotuflux_dTs (0:nlay,ncol)  ! d/d(Tsurf) [W/m2/K]
      real, intent(out) :: dtotuclfl_dTs (0:nlay,ncol)  ! d/d(Tsurf) [W/m2/K]

      ! which band OLRs to calculate?
      logical, intent(in) :: band_output (nbndlw)

      ! band OLRs and d/dTs
      real, intent(out) :: olrb      (nbndlw,ncol)  ! [W/m2]
      real, intent(out) :: dolrb_dTs (nbndlw,ncol)  ! [W/m2/K]

      ! ----- Local -----
   
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

      real :: plfrac, blay, dplankup, dplankdn, sumfac, duflx

      ! derivatives with respect to surface temperature
      real :: d_rad0_dTs, d_radlu_dTs, d_radclru_dTs

      integer :: icol          ! column index
      integer :: lev           ! level index
      integer :: ig            ! g-point index
      integer :: ibnd          ! band index
      integer :: ittot, itgas  ! lookup table indices
   
      ! flags when clear-/total-sky downward streams diverge
      logical :: down_streams_diverge

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

      ! zero spectral sum accumulators
      totuflux = 0.
      totdflux = 0.
      totuclfl = 0.
      totdclfl = 0.
      if (dudTs) then
         dtotuflux_dTs = 0.
         dtotuclfl_dTs = 0.
      end if

      ! zero band output accumulators
      if (any(band_output)) then
         olrb = 0.
         if (dudTs) dolrb_dTs = 0.
      end if 

      ! main loop
      ! g-points are "monochromatic" in CKD
      do icol = 1,ncol
      do ig = 1,ngptlw 

         ibnd = ngb(ig)

         ! factor for accumulating fluxes across g-points
         sumfac = wtdiff * delwave(ibnd) * fluxfac
      
         ! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
         ! and 1.80) as a function of total column water vapor. The function
         ! has been defined to minimize flux and cooling rate errors in these
         ! bands over a wide range of precipitable water values.

         if (ibnd == 1 .or. ibnd == 4 .or. ibnd >= 10) then
            secdiff = 1.66 
         else
            secdiff = a0(ibnd) + a1(ibnd) * exp(a2(ibnd)*pwvcm(icol))
            if (secdiff > 1.80) then
               secdiff = 1.80 
            else if (secdiff < 1.50) then
               secdiff = 1.50 
            end if
         endif
    
         ! TOA downward LW is zero
         ! and column starts "clear"
         radld = 0.; radclrd = 0. 
         down_streams_diverge = .false.

         ! Downward radiative transfer loop
         do lev = nlay,1,-1

            plfrac = pfracs(lev,ig,icol)
            blay = planklay(ibnd,lev,icol)
            dplankup = planklev(ibnd,lev,  icol) - blay
            dplankdn = planklev(ibnd,lev-1,icol) - blay

            odepth = secdiff * taug(lev,ig,icol)
            if (odepth .lt. 0.) odepth = 0. 

            ! if layer is cloudy for any gpoint/subcol then
            ! use cloudy RT (even if this subcol is clear)
!pmn: have sub-branch for clear subcol?? to save calcs

            if (cloudy(lev,icol)) then

               ! Have now reached a layer that is cloudy (in any of
               ! its gpoint/subcols) on the way down. Above here the
               ! clear- and total-sky downward streams were identical.
               ! Here and below they must be separately calculated.
               down_streams_diverge = .true.

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

            else  ! clear layer (in every gpoint/subcolumn)

               tblind = odepth/(bpade+odepth)
               itgas = tblint * tblind + 0.5 
               agas(lev) = 1. - exp_tbl(itgas)
               tausfac = tfn_tbl(itgas)

               bbdgas      = plfrac * (blay + tausfac * dplankdn)
               bbugas(lev) = plfrac * (blay + tausfac * dplankup)

               ! Update of downward LW, radld, through a layer
               ! of transmissivity T:
               !   radld -> radld * T + bbd * (1-T)
               !         -> radld + (bbd - radld) * (1-T)
               ! where bbd is effective Planck func for layer.

               radld = radld + (bbdgas - radld) * agas(lev)

            endif
            totdflux(lev-1,icol) = totdflux(lev-1,icol) + sumfac * radld

            ! Set clear-sky stream to total-sky stream as long as layers
            ! remain clear. Streams diverge when a cloud is reached, and
            ! clear-sky stream must be computed separately from that point.

            if (down_streams_diverge) then
               radclrd = radclrd + (bbdgas-radclrd) * agas(lev) 
            else
               radclrd = radld
            endif

            totdclfl(lev-1,icol) = totdclfl(lev-1,icol) + sumfac * radclrd

         enddo  ! downward loop

         ! Include the contribution of spectrally varying longwave emissivity
         ! and reflection from the surface to the upward radiative transfer.
         ! Note: Spectral and Lambertian reflection are identical for the
         ! diffusivity angle flux integration used here.
         ! Note: The emissivity is applied to plankbnd and dplankbnd_dTs when 
         ! they are defined in subroutine setcoef. 
    
         ! surface emission
         rad0 = pfracs(1,ig,icol) * plankbnd(ibnd,icol)
         if (dudTs) d_rad0_dTs = pfracs(1,ig,icol) * dplankbnd_dTs(ibnd,icol)

         ! Add in specular reflection of surface downward radiance
         reflect = 1. - semiss(ibnd,icol)
         radlu   = rad0 + reflect * radld
         radclru = rad0 + reflect * radclrd
         totuflux(0,icol) = totuflux(0,icol) + sumfac * radlu
         totuclfl(0,icol) = totuclfl(0,icol) + sumfac * radclru
         if (dudTs) then
            d_radlu_dTs   = d_rad0_dTs
            d_radclru_dTs = d_rad0_dTs
            dtotuflux_dTs(0,icol) = dtotuflux_dTs(0,icol) + sumfac * d_radlu_dTs
            dtotuclfl_dTs(0,icol) = dtotuclfl_dTs(0,icol) + sumfac * d_radclru_dTs
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

               if (dudTs) then
                  d_radlu_dTs = d_radlu_dTs *       cldfmc(lev,ig,icol)  * (1. - atot(lev)) + &
                                d_radlu_dTs * (1. - cldfmc(lev,ig,icol)) * (1. - agas(lev))
! pmn: -----
! pmn: make above more efficient
! pmn: this makes sense as the sum of the transmitted sensitivities for the cloudy
! pmn: and clear portions, but why use atot here vs. aboth for the non-derivative?
! pmn: -----
               endif

            else  ! clear layer

               radlu = radlu + (bbugas(lev) - radlu) * agas(lev)
               if (dudTs) d_radlu_dTs = d_radlu_dTs * (1. - agas(lev))

            endif

            duflx = sumfac * radlu
            totuflux(lev,icol) = totuflux(lev,icol) + duflx
            if (dudTs) dtotuflux_dTs(lev,icol) = &
               dtotuflux_dTs(lev,icol) + sumfac * d_radlu_dTs

            ! Set clear sky stream to total sky stream as long as all layers
            ! are clear (iclddn=0). Streams must be calculated separately at 
            ! all layers when a cloud is present (iclddn=1), because surface 
            ! reflectance is different for each stream.

            if (down_streams_diverge) then
               radclru = radclru + (bbugas(lev) - radclru) * agas(lev) 
            else
               radclru = radlu
            endif
            totuclfl(lev,icol) = totuclfl(lev,icol) + sumfac * radclru
            if (dudTs) then
               if (down_streams_diverge) then
                  d_radclru_dTs = d_radclru_dTs * (1. - agas(lev))
               else
                  d_radclru_dTs = d_radlu_dTs
               endif
               dtotuclfl_dTs(lev,icol) = &
                  dtotuclfl_dTs(lev,icol) + sumfac * d_radclru_dTs
            endif

         enddo  ! layer
          
         ! OLR band output (we are at TOA)
         if (band_output(ibnd)) then
            olrb(ibnd,icol) = olrb(ibnd,icol) + duflx
            if (dudTs) dolrb_dTs(ibnd,icol) = &
               dolrb_dTs(ibnd,icol) + sumfac * d_radlu_dTs
         end if

      end do  ! g-point
      end do  ! column

   end subroutine rtrnmc

end module rrtmg_lw_rtrnmc
