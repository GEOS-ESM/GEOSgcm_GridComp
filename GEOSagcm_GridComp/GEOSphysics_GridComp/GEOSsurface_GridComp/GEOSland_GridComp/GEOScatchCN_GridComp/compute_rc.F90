   module compute_rc_mod

   use clmtype

   implicit none

   private
   public compute_rc

   contains

   subroutine compute_rc(nch,nveg,tc,qa,tm,pbot,coszen,pardir,pardif, &
                         elai,esai,ityp,fveg,btran,fwet,              &
                         rc,rcdtc,rcdea,psnsun,psnsha,laisun,laisha,  &
                         dayl_fac,co2v,dtc,dea,parabs,sifsun,sifsha)

   use clm_varcon, only: tfrz
   use MAPL_SatVaporMod

   implicit none

   integer, intent(in) :: nch        ! vector length
   integer, intent(in) :: nveg       ! number of PFTs (currently 2 for catchment model)
   real, intent(in) :: tc(nch)       ! canopy temperature (K)
   real, intent(in) :: qa(nch)       ! canopy air specific humidity (kg/kg)
   real, intent(in) :: tm(nch)       ! air temperature at agcm reference height (K)
   real, intent(in) :: pbot(nch)     ! surface pressure (Pa)
   real, intent(in) :: coszen(nch)   ! cosine solar zenith angle
   real, intent(in) :: pardir(nch)   ! direct PAR (W/m2)
   real, intent(in) :: pardif(nch)   ! diffuse PAR (W/m2)
   real, intent(in) :: elai(nch,nveg)! one-sided leaf area index
   real, intent(in) :: esai(nch,nveg)! one-sided stem area index
   integer, intent(in) :: ityp(nch,nveg) ! canopy vegetation index (PFT)
   real, intent(in) :: fveg(nch,nveg)    ! canopy vegetation fractions
   real, intent(in) :: btran(nch)    ! soil water transpiration factor (0-1)
   real, intent(in) :: fwet(nch)     ! fraction of canopy that is wet (0-1)
   real, intent(in) :: co2v(nch)     ! atmospheric carbon dioxide concentration
   real, intent(in) :: dayl_fac(nch) ! daylength factor (0-1)

   real, intent(out) :: rc(nch)      ! canopy stomatal resistance (s/m)
   real, intent(out) :: rcdtc(nch)   ! canopy stomatal resistance (s/m) for Tc+d(Tc)
   real, intent(out) :: rcdea(nch)   ! canopy stomatal resistance (s/m) for ea+d(ea)
   real, intent(out) :: psnsun(nch,nveg)    ! sunlit foliage photosynthesis (umol co2 /m**2/ s) [always +]
   real, intent(out) :: psnsha(nch,nveg)    ! shaded foliage photosynthesis (umol co2 /m**2/ s) [always +]
   real, intent(out) :: sifsun(nch,nveg)    ! sunlit foliage fluorescence
   real, intent(out) :: sifsha(nch,nveg)    ! shaded foliage fluorescence
   real, intent(out) :: laisun(nch,nveg)    ! sunlit projected leaf area index
   real, intent(out) :: laisha(nch,nveg)    ! shaded projected leaf area index
   real, intent(out) :: parabs(nch)         ! total absorbed PAR

!  local
  
   real ei(nch)        ! vapor pressure inside leaf (sat vapor press at tc) (Pa)
   real ea(nch)        ! vapor pressure of canopy air (Pa)
   real o2(nch)        ! atmospheric o2 concentration (Pa)
   real co2(nch)       ! atmospheric co2 concentration (Pa)
   real rb(nch)        ! boundary layer resistance (s/m)
   real tl(nch)        ! canopy temperature (K)
   real parsun(nch,nveg)  ! par absorbed per unit lai (w/m**2)
   real parsha(nch,nveg)  ! par absorbed per unit lai (w/m**2)
   real slasun(nch,nveg)  ! specific leaf area for sunlit canopy, projected area basis (m^2/gC)
   real slasha(nch,nveg)  ! specific leaf area for shaded canopy, projected area basis (m^2/gC)
   real  rssun(nch,nveg)  ! sunlit stomatal resistance (s/m)
   real  rssha(nch,nveg)  ! shaded stomatal resistance (s/m)
   real psn(nch,nveg)     ! foliage photosynthesis (umol co2 /m**2/ s) [always +]
   real sif(nch,nveg)     ! foliage fluorescence
   real deldT(nch)     ! d(es)/d(T)

   real :: gdir                      ! leaf projection in solar direction (0-1)

   integer n, nv, ivt
   real qsatl, qsatldT, par
   real*8 cosz, ext, t1, t2, rcs, fvs, rs, rcd, fsun
   real*8 vai, wl, ws, rho, tau, omegal, omega, sun_add, tot_aid, sun_aid, sun_aii, sha_aid, sha_aii 
   real*8 chil, phi1, phi2, temp, fabd, fabi, avmu, temp0, temp1, temp2, asu, betadl, betail
   real*8 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
   real*8 betad, betai, b, c1, d, d1, d2, f, h, sigma, p1, p2, p3, p4, s1, s2, u1, u2, u3
   real*8 h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, ftdd, ftid, albd, ftii, albi

   real, intent(in) :: dtc ! canopy temperature perturbation (K) [approx 1:10000]
   real, intent(in) :: dea ! vapor pressure perturbation (Pa) [approx 1:10000]

   real, pointer :: slatop(:)      ! specific leaf area at top of canopy, projected area basis [m^2/gC]
   real, pointer :: dsladlai(:)    ! dSLA/dLAI, projected area basis [m^2/gC]
   real, pointer :: xl(:)          ! leaf/stem orientation index
   real, pointer :: rhol(:)        ! leaf reflectance (visible)
   real, pointer :: rhos(:)        ! stem reflectance (visible)
   real, pointer :: taul(:)        ! leaf transmittance (visible)
   real, pointer :: taus(:)        ! stem transmittance (visible)

   real, parameter :: omegas = 0.80 ! two-stream parameter omega for snow (visible)
   real, parameter :: betads = 0.50 ! two-stream parameter betad for snow
   real, parameter :: betais = 0.50 ! two-stream parameter betai for snow
   real, parameter :: albgrd = 0.10 ! ground albedo direct
   real, parameter :: albgri = 0.10 ! ground albedo indirect

! assign local pointers to derived type arrays
! --------------------------------------------
   slatop   => pftcon%slatop
   dsladlai => pftcon%dsladlai
   xl       => pftcon%xl
   rhol     => pftcon%rhol
   rhos     => pftcon%rhos
   taul     => pftcon%taul
   taus     => pftcon%taus

      o2(:) = 0.20946*pbot   !  O2 partial pressure constant ratio
     co2(:) = co2v(:)*pbot   ! C02 partial pressure constant ratio [internal leaf CO2 partial pressure]
      rb(:) = 10.            ! gkw: for now, assume a small value for rb (see 8/3/10 email) 
      ea(:) = pbot(:) * qa(:) / (0.622 + qa(:))   ! canopy air vapor pressure (Pa)

     parabs(:) = 0.          ! initialize absorbed PAR to zero

! compute saturation vapor pressure
! ---------------------------------
   do n = 1,nch
     ei(n) = MAPL_EQsat(tc(n),DQ=deldT(n))
   end do

! set solar fluxes
! ----------------
   do n = 1,nch

     par = pardir(n) + pardif(n)

! partition PAR into sunlit and shaded canopy for each vegetation type
! --------------------------------------------------------------------
     do nv = 1,nveg

       if(par > 0. .and. elai(n,nv) > 0. .and. fveg(n,nv) > 1.e-4 .and. ityp(n,nv) > 0) then

         ivt = ityp(n,nv) ! mapped PFT

         cosz = max(0.001, coszen(n))

! Calculate two-stream parameters omega, betad, betai, avmu, gdir, ext.
! Omega, betad, betai are adjusted for snow. Values for omega*betad
! and omega*betai are calculated and then divided by the new omega
! because the product omega*betai, omega*betad is used in solution.
! Also, the transmittances and reflectances (tau, rho) are linear
! weights of leaf and stem values.

! leaf projection in solar direction (0-1)
! ----------------------------------------
         chil = min(max(xl(ivt),-0.4),0.6)
         if(abs(chil) <= 0.01) chil = 0.01
         phi1 = 0.5 - 0.633*chil - 0.330*chil*chil
         phi2 = 0.877 * (1. - 2.*phi1)
         gdir = phi1 + phi2*cosz ! 3.3

         ext = gdir/cosz
         t1 = min(ext*elai(n,nv), 40.0)
         t2 = exp(-t1)
         fsun = (1.-t2)/t1

         avmu = (1. - phi1/phi2 * log((phi1+phi2)/phi1)) / phi2 ! average diffuse optical depth 3.4
         temp0 = gdir + phi2*cosz
         temp1 = phi1*cosz
         temp2 = ( 1. - temp1/temp0 * log((temp1+temp0)/temp1) )

         vai = elai(n,nv) + esai(n,nv) ! elai+esai

         wl = elai(n,nv) / vai ! LAI weight
         ws = esai(n,nv) / vai ! SAI weight

         rho = rhol(ivt)*wl + rhos(ivt)*ws ! weighted reflectance
         tau = taul(ivt)*wl + taus(ivt)*ws ! weighted transmittance

         omegal = rho + tau ! snow-free visible scattering coefficient (3.8, 3.11, 3.12)

         asu = 0.5*omegal*gdir/temp0 *temp2 ! single scattering albedo
         betadl = (1.+avmu*ext)/(omegal*avmu*ext)*asu                        ! betad for leaves
         betail = 0.5 * ((rho+tau) + (rho-tau) * ((1.+chil)/2.)**2) / omegal ! betai for leaves

	 if (tc(n) > tfrz) then	 ! no snow
 	   tmp0 = omegal
           tmp1 = betadl
           tmp2 = betail
	  else
           tmp0 =  (1.-fwet(n))*omegal	      + fwet(n)*omegas
           tmp1 = ((1.-fwet(n))*omegal*betadl + fwet(n)*omegas*betads) / tmp0
           tmp2 = ((1.-fwet(n))*omegal*betail + fwet(n)*omegas*betais) / tmp0
	 end if
         omega = tmp0	      
         betad = tmp1 ! upscatter parameter for direct beam radiation
         betai = tmp2 ! upscatter parameter for diffuse radiation

! Absorbed, reflected, transmitted fluxes per unit incoming radiation
! -------------------------------------------------------------------
         b = 1. - omega + omega*betai ! 3.21
         c1 = omega*betai             ! 3.22
         tmp0 = avmu*ext
         d = tmp0 * omega*betad       ! 3.23
         f = tmp0 * omega*(1.-betad)  ! 3.24
         tmp1 = b*b - c1*c1
         h = sqrt(tmp1) / avmu        ! 3.25
         sigma = tmp0*tmp0 - tmp1     ! 3.26
         p1 = b + avmu*h              ! 3.32
         p2 = b - avmu*h              ! 3.33
         p3 = b + tmp0                ! 3.34
         p4 = b - tmp0                ! 3.35
         
         t1 = min(h*vai, 40.)
         s1 = exp(-t1)                ! 3.30
         t1 = min(ext*vai, 40.)
         s2 = exp(-t1)                ! 3.31

! unit incoming direct flux
! -------------------------
         u1 = b - c1/albgrd           ! 3.27
         u2 = b - c1*albgrd           ! 3.28
         u3 = f + c1*albgrd           ! 3.29

         tmp2 = u1 - avmu*h
         tmp3 = u1 + avmu*h
         d1 = p1*tmp2/s1 - p2*tmp3*s1 ! 3.36
         tmp4 = u2 + avmu*h
         tmp5 = u2 - avmu*h
         d2 = tmp4/s1 - tmp5*s1       ! 3.37
         h1 = -d*p4 - c1*f            ! 3.38

         tmp6 = d - h1*p3/sigma
         tmp7 = (d - c1 - h1/sigma*(u1+tmp0)) * s2
         h2 = (tmp6*tmp2/s1 - p2*tmp7) / d1   ! 3.39
         h3 = - (tmp6*tmp3*s1 - p1*tmp7) / d1 ! 3.40
         h4 = -f*p3 - c1*d                    ! 3.41
         tmp8 = h4/sigma
         tmp9 = (u3 - tmp8*(u2-tmp0)) * s2
         h5 = - (tmp8*tmp4/s1 + tmp9) / d2    ! 3.42
         h6 = (tmp8*tmp5*s1 + tmp9) / d2      ! 3.43
         h7 = (c1*tmp2) / (d1*s1)             ! 3.44
         h8 = (-c1*tmp3*s1) / d1              ! 3.45
         h9 = tmp4 / (d2*s1)                  ! 3.46
         h10 = (-tmp5*s1) / d2                ! 3.47

! Downward direct and diffuse fluxes below vegetation
! ---------------------------------------------------
         ftdd = s2
         ftid = h4*s2/sigma + h5*s1 + h6/s1

! Flux reflected by vegetation (direct)
! -------------------------------------
         albd = h1/sigma + h2 + h3

! Flux absorbed by vegetation (direct)
! ------------------------------------
         fabd = 1. - albd - (1.-albgrd)*ftdd - (1.-albgri)*ftid

! unit incoming diffuse flux
! --------------------------
         u1 = b - c1/albgri
         u2 = b - c1*albgri
         u3 = f + c1*albgri

         tmp2 = u1 - avmu*h
         tmp3 = u1 + avmu*h
         d1 = p1*tmp2/s1 - p2*tmp3*s1
         tmp4 = u2 + avmu*h
         tmp5 = u2 - avmu*h
         d2 = tmp4/s1 - tmp5*s1
         h1 = -d*p4 - c1*f
         tmp6 = d - h1*p3/sigma
         tmp7 = (d - c1 - h1/sigma*(u1+tmp0)) * s2
         h2 = (tmp6*tmp2/s1 - p2*tmp7) / d1
         h3 = - (tmp6*tmp3*s1 - p1*tmp7) / d1
         h4 = -f*p3 - c1*d
         tmp8 = h4/sigma
         tmp9 = (u3 - tmp8*(u2-tmp0)) * s2
         h5 = - (tmp8*tmp4/s1 + tmp9) / d2
         h6 = (tmp8*tmp5*s1 + tmp9) / d2
         h7 = (c1*tmp2) / (d1*s1)
         h8 = (-c1*tmp3*s1) / d1
         h9 = tmp4 / (d2*s1)
         h10 = (-tmp5*s1) / d2

! Downward direct and diffuse fluxes below vegetation
! ---------------------------------------------------
         ftii = h9*s1 + h10/s1

! Flux reflected by vegetation (indirect)
! ----------------------------
         albi = h7 + h8

! Flux absorbed by vegetation (indirect)
! --------------------------------------
         fabi = 1. - albi - (1.-albgri)*ftii

! partition PAR
! -------------         
         sun_add = pardir(n) * (1.-ftdd) * (1.-omega)
         tot_aid = max((pardir(n) * fabd) - sun_add , 0.)

         sun_aid = tot_aid * fsun
         sun_aii = pardif(n)*fabi*fsun
         sha_aid = tot_aid * (1.-fsun)
         sha_aii = pardif(n)*fabi*(1.-fsun)

         parsun(n,nv) = sun_add + sun_aid + sun_aii
         parsha(n,nv) = sha_aid + sha_aii

         parsun(n,nv) = parsun(n,nv) * wl ! sunlit canopy PAR for leaves
         parsha(n,nv) = parsha(n,nv) * wl ! shaded canopy PAR for leaves

         parabs(n) = parabs(n) + (parsun(n,nv) + parsha(n,nv))*fveg(n,nv) ! save absorbed PAR for FPAR calculation

         if(elai(n,nv) .gt. 0.01) then
           laisun(n,nv) = elai(n,nv)*fsun
           laisha(n,nv) = elai(n,nv)*(1.-fsun)

           parsun(n,nv) = parsun(n,nv) / laisun(n,nv) ! express as PAR per unit LAI 
           parsha(n,nv) = parsha(n,nv) / laisha(n,nv) ! express as PAR per unit LAI

           slasun(n,nv) = (t2*dsladlai(ivt)*ext*elai(n,nv) + t2*dsladlai(ivt) + &
                              t2*slatop(ivt)*ext - dsladlai(ivt) - &
                              slatop(ivt)*ext) / (ext*(t2-1.))

           slasha(n,nv) = ((slatop(ivt) + &
                             (dsladlai(ivt) * elai(n,nv)/2.0)) * elai(n,nv) - &
                             laisun(n,nv)*slasun(n,nv)) / laisha(n,nv)
          else
                 ! special case for low elai
           laisun(n,nv) = elai(n,nv)
           laisha(n,nv) = 0.
           slasun(n,nv) = slatop(ivt)
           slasha(n,nv) = 0.
           parsun(n,nv) = parsun(n,nv) / laisun(n,nv)
           parsha(n,nv) = 0.
         endif

        else

         laisun(n,nv) = 0.
         laisha(n,nv) = elai(n,nv)
         slasun(n,nv) = 0.
         slasha(n,nv) = 0.
         parsun(n,nv) = 0.
         parsha(n,nv) = 0.
       endif
       
     end do ! end PFT loop

   end do   ! end column loop


!  compute stomatal resistance using CLM routine; also compute photosynthesis

!  obtain stomatal resistance and photosynthesis

   call stomata(nch, nveg, ei, ea, o2, co2, rb, dayl_fac, pbot, ityp, tm, tc, &
                btran, parsun, slasun, rssun, psnsun, sifsun) ! sunlit
   call stomata(nch, nveg, ei, ea, o2, co2, rb, dayl_fac, pbot, ityp, tm, tc, &
                btran, parsha, slasha, rssha, psnsha, sifsha) ! shaded


!  combine resistance as reciprocal of vegetation weighted conductance (weighted harmonic sum)

!DIR$ NOVECTOR
   do n = 1,nch
     rcs = 0.
     do nv = 1,nveg
       rs = laisun(n,nv)/rssun(n,nv) + laisha(n,nv)/rssha(n,nv)
       rcs = rcs + fveg(n,nv)*rs
     end do 
     rc(n) = 1./max(rcs,5.e-5) + rb(n)
   end do


!  compute resistance with small delta ea

   ea(:) = ea(:) + dea

   call stomata(nch, nveg, ei, ea, o2, co2, rb, dayl_fac, pbot, ityp, tm, tc, &
                btran, parsun, slasun, rssun, psn, sif) ! sunlit
   call stomata(nch, nveg, ei, ea, o2, co2, rb, dayl_fac, pbot, ityp, tm, tc, &
                btran, parsha, slasha, rssha, psn, sif) ! shaded

!DIR$ NOVECTOR
   do n = 1,nch
     rcs = 0.
     do nv = 1,nveg
       rs = laisun(n,nv)/rssun(n,nv) + laisha(n,nv)/rssha(n,nv)
       rcs = rcs + fveg(n,nv)*rs
     end do 
     rcdea(n) = 1./max(rcs,5.e-5) + rb(n)
   end do


!  compute resistance with small delta Tc

   tl(:) = tc(:) + dtc
   ei(:) = ei(:) + deldT(:)*dtc ! ei=esat(Tc)+[d(esat)/d(Tc)]dTc

   ea(:) = pbot(:) * qa(:) / (0.622 + qa(:))  ! reset input canopy air vapor pressure (Pa)

   call stomata(nch, nveg, ei, ea, o2, co2, rb, dayl_fac, pbot, ityp, tm, tl, &
                btran, parsun, slasun, rssun, psn, sif) ! sunlit
   call stomata(nch, nveg, ei, ea, o2, co2, rb, dayl_fac, pbot, ityp, tm, tl, &
                btran, parsha, slasha, rssha, psn, sif) ! shaded

!DIR$ NOVECTOR
   do n = 1,nch
     rcs = 0.
     do nv = 1,nveg
       rs = laisun(n,nv)/rssun(n,nv) + laisha(n,nv)/rssha(n,nv)
       rcs = rcs + fveg(n,nv)*rs
     end do 
     rcdtc(n) = 1./max(rcs,5.e-5) + rb(n)
   end do

   end subroutine compute_rc 

!*******************************************************************************
!  gkw: this is how "btran" is calculated in CanopyFluxes

!!!  btran(:) = 0.
!!!
!!!  do j = 1,nlevgrnd
!!!    p = filterp(f)
!!!    c = pcolumn(p)
!!!    l = plandunit(p)
!!!
!!!! Root resistance factors
!!!
!!!    vol_ice = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
!!!    eff_porosity = watsat(c,j)-vol_ice
!!!    vol_liq = min(eff_porosity, h2osoi_liq(c,j)/(dz(c,j)*denh2o))
!!!    if (vol_liq .le. 0. .or. t_soisno(c,j) .le. tfrz-2.) then
!!!     rootr(p,j) = 0.
!!!      else
!!!      s_node = max(vol_liq/eff_porosity,0.01)
!!!      smp_node = max(smpsc(ivt(p)), -sucsat(c,j)*s_node**(-bsw(c,j)))
!!!
!!!      rresis(p,j) = min( (eff_porosity/watsat(c,j))* &
!!!     	  (smp_node - smpsc(ivt(p))) / (smpso(ivt(p)) - smpsc(ivt(p))), 1.)
!!!      rootr(p,j) = rootfr(p,j)*rresis(p,j)
!!!      btran(p) = btran(p) + rootr(p,j)
!!!    endif 
!!!  end do
!*******************************************************************************

!*******************************************************************************
! gkw: this is how rb is calculated in CanopyFluxes

! Determine aerodynamic resistances & Bulk boundary layer resistance of leaves
! gkw: code imported from CanopyFluxes for later use 8/3/2010

!!!  do n = 1,nch
!!!    ram1 = 1./(ustar(n)*ustar(n)/um(n))            ! gkw: 5.55
!!!    uaf  = um(p)*sqrt( 1./(ram1*um(n)) )           ! gkw: 5.100
!!!    rb(n) = 100.*sqrt(dleaf(ivt))/sqrt(uaf)        ! gkw: 5.109
!!!  end do
!*******************************************************************************

   subroutine Stomata (nch, nveg, ei, ea, o2, co2, rb, dayl_factor, forc_pbot, ityp, tgcm, tl, &
                       btran, apar, sla, rs, psn, sif)
!
! !DESCRIPTION: 
! Leaf stomatal resistance and leaf photosynthesis. Modifications for CN code.

! !REVISION HISTORY:
! 22 January 2004: Created by Peter Thornton
! 4/14/05: Peter Thornton: Converted Ci from local variable to pps struct member
!    now returns cisun or cisha per pft as implicit output argument.
!    Also sets alphapsnsun and alphapsnsha.
! 4/25/05, Peter Thornton: Adopted as the default code for CLM, together with
!   modifications for sun/shade canopy.  Renamed from StomataCN to Stomata,
!   and eliminating the older Stomata subroutine
! 3/6/09: Peter Thornton; added dayl_factor control on Vcmax, from Bill Bauerle

! !USES:
     use shr_const_mod, only: SHR_CONST_TKFRZ,SHR_CONST_RGAS
!
! !ARGUMENTS:
     implicit none
     integer , intent(in)    :: nch          ! number of land tiles
     integer , intent(in)    :: nveg         ! number of PFTs
     real, intent(in)    :: ei(nch)          ! vapor pressure inside leaf (sat vapor press at tl) (pa)
     real, intent(in)    :: ea(nch)          ! vapor pressure of canopy air (pa)
     real, intent(in)    :: o2(nch)          ! atmospheric o2 concentration (pa)
     real, intent(in)    :: co2(nch)         ! atmospheric co2 concentration (pa)
     real, intent(inout) :: rb(nch)          ! boundary layer resistance (s/m)
     real, intent(in)    :: dayl_factor(nch) ! scalar (0-1) for daylength
     real, intent(in)    :: forc_pbot(nch)   ! atmospheric pressure (Pa)
     integer, intent(in) :: ityp(nch,nveg)   ! vegetation type
     real, intent(in) :: tgcm(nch)           ! air temperature at agcm reference height (kelvin)
     real, intent(in) :: tl(nch)             ! leaf temperature (Kelvin)
     real, intent(in) :: btran(nch)          ! soil water transpiration factor (0 to 1)

     real, intent(in) :: apar(nch,nveg)      ! par absorbed per unit lai (w/m**2)
     real, intent(in) :: sla(nch,nveg)       ! specific leaf area, projected area basis (m^2/gC)
!
! !CALLED FROM:
! subroutine CanopyFluxes in this module
!
! !LOCAL VARIABLES:
!
!  these are function of PFT only; map into Catchment type

     real, pointer :: qe25(:)	   ! quantum efficiency at 25C (umol CO2 / umol photon)
     real, pointer :: c3psn(:)	   ! photosynthetic pathway: 0. = c4, 1. = c3
     real, pointer :: mp(:)	   ! slope of conductance-to-photosynthesis relationship
     real, pointer :: leafcn(:)	   ! leaf C:N (gC/gN)
     real, pointer :: flnr(:)	   ! fraction of leaf N in the Rubisco enzyme (gN Rubisco / gN leaf)
     real, pointer :: fnitr(:)	   ! foliage nitrogen limitation factor (-)
     real, pointer :: slatop(:)	   ! specific leaf area at top of canopy, projected area basis [m^2/gC]
     real, pointer :: dsladlai(:)  ! dSLA/dLAI, projected area basis [m^2/gC]
!
! local pointers to implicit out variables
!
     real, intent(out) :: rs(nch,nveg)          ! leaf stomatal resistance (s/m)
     real, intent(out) :: psn(nch,nveg)         ! foliage photosynthesis (umol co2 /m**2/ s) [always +]
     real, intent(out) :: sif(nch,nveg)         ! foliage fluorescence
!    real, intent(out) :: ci(nch)          ! intracellular leaf CO2 (Pa)
!    real, intent(out) :: lnc(nch)         ! leaf N concentration per unit projected LAI (gN leaf/m^2)
!    real, intent(out) :: vcmx(nch)        ! maximum rate of carboxylation (umol co2/m**2/s)
!
!
! !LOCAL VARIABLES:
!EOP
!
     real, parameter :: mpe = 1.e-6   ! prevents overflow error if division by zero
     integer , parameter :: niter = 3     ! number of iterations
     integer  :: n,nv ! indices
     integer  :: iter    ! iteration index
     integer  :: ivt    ! mapped PFT
     real :: ab      ! used in statement functions
     real :: bc      ! used in statement functions
     real :: f1      ! generic temperature response (statement function)
     real :: f2      ! generic temperature inhibition (statement function)
     real :: tc      ! leaf temperature (degree celsius)
     real :: cs      ! co2 concentration at leaf surface (pa)
     real :: kc      ! co2 michaelis-menten constant (pa)
     real :: ko      ! o2 michaelis-menten constant (pa)
     real :: atmp    ! intermediate calculations for rs
     real :: btmp    ! intermediate calculations for rs
     real :: ctmp    ! intermediate calculations for rs
     real :: q       ! intermediate calculations for rs
     real :: r1,r2   ! roots for rs
     real :: ppf     ! absorb photosynthetic photon flux (umol photons/m**2/s)
     real :: wc      ! rubisco limited photosynthesis (umol co2/m**2/s)
     real :: wj      ! light limited photosynthesis (umol co2/m**2/s)
     real :: we      ! export limited photosynthesis (umol co2/m**2/s)
     real :: cp      ! co2 compensation point (pa)
     real :: awc     ! intermediate calcuation for wc
     real :: j       ! electron transport (umol co2/m**2/s)
     real :: cea     ! constrain ea or else model blows up
     real :: cf      ! s m**2/umol -> s/m
     real :: rsmax0  ! maximum stomatal resistance [s/m]
     real :: kc25    ! co2 michaelis-menten constant at 25c (pa)
     real :: akc     ! q10 for kc25
     real :: ko25    ! o2 michaelis-menten constant at 25c (pa)
     real :: ako     ! q10 for ko25
     real :: avcmx   ! q10 for vcmx25
     real :: bp      ! minimum leaf conductance (umol/m**2/s)
     ! additional variables for new treatment of Vcmax, Peter Thornton, 1/26/04
     real :: act25   ! (umol/mgRubisco/min) Rubisco activity at 25 C
     real :: act     ! (umol/mgRubisco/min) Rubisco activity
     real :: q10act  ! (DIM) Q_10 for Rubisco activity
     real :: fnr     ! (gRubisco/gN in Rubisco)
     real :: slasun  ! specific leaf area for sunlit canopy, projected area basis (m^2/gC)
     real :: slasha  ! specific leaf area for shaded canopy, projected area basis (m^2/gC)
     real :: ci      ! intracellular leaf CO2 (Pa)
     real :: lnc     ! leaf N concentration per unit projected LAI (gN leaf/m^2)
     real :: vcmx    ! maximum rate of carboxylation (umol co2/m**2/s)
     real :: je      ! actual electron transport
     real :: xn      ! je/j
     real :: fs      ! fluorescnce yield at Fs

!------------------------------------------------------------------------------

     ! Set statement functions

     f1(ab,bc) = ab**((bc-25.)/10.)
     f2(ab) = 1. + exp((-2.2e05+710.*(ab+SHR_CONST_TKFRZ))/(SHR_CONST_RGAS*0.001*(ab+SHR_CONST_TKFRZ)))

     ! assign local pointers to derived type arrays

     qe25     => pftcon%qe25
     c3psn    => pftcon%c3psn
     mp       => pftcon%mp
     leafcn   => pftcon%leafcn
     flnr     => pftcon%flnr
     fnitr    => pftcon%fnitr
     slatop   => pftcon%slatop
     dsladlai => pftcon%dsladlai
 
     ! Set constant values

     kc25  = 30.
     akc   = 2.1
     ko25  = 30000.
     ako   = 1.2
     avcmx = 2.4
     bp    = 2000.

     ! New constants for CN code, added 1/26/04

     act25 = 3.6
     q10act = 2.4
     fnr = 7.16
     
     ! Convert rubisco activity units from umol/mgRubisco/min -> umol/gRubisco/s

     act25 = act25 * 1000.0 / 60.0

     do n = 1, nch

       do nv = 1,nveg

        ivt = ityp(n,nv) !  mapped vegetation type into CLM PFT

        ! Initialize rs=rsmax and psn=0 because calculations are performed only
        ! when apar > 0, in which case rs <= rsmax and psn >= 0
        ! Set constants

        rsmax0 = 2.e4
        cf = forc_pbot(n)/(SHR_CONST_RGAS*0.001*tgcm(n))*1.e06

        if (apar(n,nv) <= 0.) then          ! night time
           rs(n,nv) = min(rsmax0, 1./bp * cf)
           psn(n,nv) = 0.
           sif(n,nv) = 0.
           lnc = 0.
           vcmx = 0.
        else                             ! day time

           tc = tl(n) - SHR_CONST_TKFRZ
           ppf = 4.6 * apar(n,nv)                            ! gkw: used in 8.3
           j = ppf * qe25(ivt)                            ! gkw: used in 8.3
           kc = kc25 * f1(akc,tc)                         ! gkw: 8.5  
           ko = ko25 * f1(ako,tc)                         ! gkw: 8.6
           awc = kc * (1.+o2(n)/ko)
           cp = 0.5*kc/ko*o2(n)*0.21                      ! gkw: 8.7

           ! Modification for shrubs proposed by X.D.Z
           ! Why does he prefer this line here instead of in subr.
           ! CanopyFluxes? (slevis)

           ! new calculations for vcmax, 1/26/04
           lnc = 1. / (sla(n,nv) * leafcn(ivt))
		   act = act25 * f1(q10act,tc)
           vcmx = lnc * flnr(ivt) * fnr * act / f2(tc) * btran(n) * dayl_factor(n) !* fnitr(ivt) ! gkw: 8.13
           
           ! First guess ci

           ci = 0.7*co2(n)*c3psn(ivt) + 0.4*co2(n)*(1.-c3psn(ivt))

           ! rb: s/m -> s m**2 / umol

           rb(n) = rb(n)/cf 

           ! Constrain ea

           cea = max(0.25*ei(n)*c3psn(ivt)+0.40*ei(n)*(1.-c3psn(ivt)), min(ea(n),ei(n)) ) 

           ! ci iteration for 'actual' photosynthesis

           do iter = 1, niter
              wj = max(ci-cp,0.)*j/(ci+2.*cp)*c3psn(ivt) + j*(1.-c3psn(ivt))           ! gkw: 8.3
              wc = max(ci-cp,0.)*vcmx/(ci+awc)*c3psn(ivt) + vcmx*(1.-c3psn(ivt)) ! gkw: 8.2
              we = 0.5*vcmx*c3psn(ivt) + 4000.*vcmx*ci/forc_pbot(n)*(1.-c3psn(ivt)) ! gkw: 8.4
              psn(n,nv) = min(wj,wc,we)                                             ! leaf photosynthesis
              cs = max( co2(n)-1.37*rb(n)*forc_pbot(n)*psn(n,nv), mpe )             ! gkw: 8.25
              atmp = mp(ivt)*psn(n,nv)*forc_pbot(n)*cea / (cs*ei(n)) + bp           ! gkw: term of 8.27
              btmp = ( mp(ivt)*psn(n,nv)*forc_pbot(n)/cs + bp ) * rb(n) - 1.        ! gkw: term of 8.27 
              ctmp = -rb(n)
              if (btmp >= 0.) then
                 q = -0.5*( btmp + sqrt(btmp*btmp-4.*atmp*ctmp) )                ! gkw: solving 8.27 -> 8.1
              else
                 q = -0.5*( btmp - sqrt(btmp*btmp-4.*atmp*ctmp) )
              end if
              r1 = q/atmp
              r2 = ctmp/q
              rs(n,nv) = max(r1,r2)
              ci = max( cs-psn(n,nv)*forc_pbot(n)*1.65*rs(n,nv), 0. )               ! gkw: 8.28
           end do

           ! rs, rb:  s m**2 / umol -> s/m 

           rs(n,nv) = min(rsmax0, rs(n,nv)*cf)
           rb(n) = rb(n) * cf 

! fluorescence; code from Jung-Eun Lee, implemented & modified by gkw 1/28/14
! ------------
           je = max(psn(n,nv)*(ci+2.*cp)/max(ci+2.*cp-3.*c3psn(ivt)*cp,1.e-8) , 0.)

!!!        xn=1.-je/j/0.8
           xn=1.-je/j      ! gkw: 0.8 factor removed 20141108 (email from Jung-Eun)
           xn=max(xn,0.)   ! gkw: added 1/31/14

           if (wj <= 0.)  xn=0.
           call fluorescence(xn,fs)
           sif(n,nv) = fs*ppf

        end if

       end do

     end do

  end subroutine Stomata

!------------------------------------------------------------------------------
!
! !IROUTINE: Fluorescence
!
! !INTERFACE:
   subroutine fluorescence(x,fs)
!
! !DESCRIPTION: 
! Chlorophyll fluorescence
! writen by Jung-Eun Lee using van der Tol and Berry (2012)

! !USES:
     implicit none
     real, intent(in)    :: x       ! degree of light saturation
     real, intent(out)   :: fs      ! fluorescence yield
     real :: Kn      ! rate constant for non-photochemical quenching
     real :: Kf      ! rate constant for fluorescence
     real :: Kd      ! rate constant for thermal deactivation at Fm
     real :: Kp      ! rate constant for photochemisty
     real :: po0
     real :: ps
     real :: fo0
     real :: fo      ! fluorescnce yield at Fo
     real :: fm      ! fluorescnce yield at Fm
     real :: fm0
     real :: eta
     real :: qQ 
     real :: qE 

     Kf          = 0.05                 ! rate constant for fluorescence
     Kd          = 0.95                 ! rate constant for thermal deactivation at Fm
     Kp          = 4.0                  ! rate constant for photochemisty

     po0         = Kp/(Kf+Kd+Kp)        ! dark photochemistry fraction (Genty et al., 1989)
     ps          = po0*(1.-x)           ! photochemical yield
     Kn          = (6.2473 * x - 0.5944)*x ! empirical fit to Flexas' data
!    Kn          = (3.9867 * x - 1.0589)*x ! empirical fit to Flexas, Daumard, Rascher, Berry data

     fo0         = Kf/(Kf+Kp+Kd)        ! dark adapted fluorescence yield Fo
     fo          = Kf/(Kf+Kp+Kd+Kn)     ! dark adapted fluorescence yield Fo
     fm          = Kf/(Kf   +Kd+Kn)     ! light adapted fluorescence yield Fm
     fm0         = Kf/(Kf   +Kd)        ! light adapted fluorescence yield Fm
     fs          = fm*(1.-ps)           ! fluorescence as fraction of PAR
     eta         = fs/fo0               ! fluorescence as fraction of dark adapted

     qQ          = 1.-(fs-fo)/(fm-fo)   ! photochemical quenching
     qE          = 1.-(fm-fo)/(fm0-fo0) !non-photochemical quenching

  end subroutine fluorescence

  end module compute_rc_mod

