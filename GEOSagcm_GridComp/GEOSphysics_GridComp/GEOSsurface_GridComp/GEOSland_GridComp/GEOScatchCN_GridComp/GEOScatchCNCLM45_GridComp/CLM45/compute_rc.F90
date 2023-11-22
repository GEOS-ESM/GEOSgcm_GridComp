#define CN
#undef  CND
   module compute_rc_mod

   use clmtype
   use update_model_para4cn, only : LocalTileID

   implicit none

   private
   public compute_rc

   contains

   subroutine compute_rc(nch, nveg, tc, qa, t10, tm, pbot, coszen, pardir, pardif, &
                         albgrd, albgri, elai, esai, ityp, fveg, btran, fwet,    &
                         rc, rcdtc, rcdea, psnsun, psnsha, laisun, laisha,      &
                         dayl_fac, co2v, dtc, dea, parabs, sifsun, sifsha,      &
                         lmrsun, lmrsha, fpar_sf)

   use clm_varcon, only: tfrz
   use MAPL_SatVaporMod

   implicit none

   integer, intent(in):: nch               ! vector length
   integer, intent(in):: nveg              ! number of PFTs (currently 2 for catchment model)
   real, intent(in):: tc(nch)              ! canopy temperature (K)
   real, intent(in):: qa(nch)              ! canopy air specific humidity (kg/kg)
   real, intent(in):: t10(nch)             ! 10-day "running mean" of the 2 m temperature (K)
   real, intent(in):: tm(nch)              ! air temperature at agcm reference height (K)
   real, intent(in):: pbot(nch)            ! surface pressure (Pa)
   real, intent(in):: coszen(nch)          ! cosine solar zenith angle
   real, intent(in):: pardir(nch)          ! direct PAR (W/m2)
   real, intent(in):: pardif(nch)          ! diffuse PAR (W/m2)
   real, intent(in):: albgrd(nch, nveg)     ! ground albedo visible direct 
   real, intent(in):: albgri(nch, nveg)     ! ground albedo visible diffuse 
   real, intent(in):: elai(nch, nveg)       ! one-sided leaf area index
   real, intent(in):: esai(nch, nveg)       ! one-sided stem area index
   integer, intent(in):: ityp(nch, nveg)    ! canopy vegetation index (PFT)
   real, intent(in):: fveg(nch, nveg)       ! canopy vegetation fractions
   real, intent(in):: btran(nch)           ! soil water transpiration factor (0-1)
   real, intent(in):: fwet(nch)            ! fraction of canopy that is wet (0-1)
   real, intent(in):: co2v(nch)            ! atmospheric carbon dioxide concentration
   real, intent(in):: dayl_fac(nch)        ! daylength factor (0-1)
   real, intent(in), optional:: fpar_sf(nch, nveg)   ! FPAR Scale factor = SCALED_FPAR/CLM4_FPAR

   real, intent(out):: rc(nch)             ! canopy stomatal resistance (s/m). NOTE: it's resistance, not conductance, fzeng, 20 Feb 2018
   real, intent(out):: rcdtc(nch)          ! canopy stomatal resistance (s/m) for Tc+d(Tc)
   real, intent(out):: rcdea(nch)          ! canopy stomatal resistance (s/m) for ea+d(ea)
   real, intent(out):: psnsun(nch, nveg)    ! sunlit foliage photosynthesis (umol co2/m**2/s) [always +]
   real, intent(out):: psnsha(nch, nveg)    ! shaded foliage photosynthesis (umol co2/m**2/s) [always +]
   real, intent(out):: sifsun(nch, nveg)    ! sunlit foliage fluorescence
   real, intent(out):: sifsha(nch, nveg)    ! shaded foliage fluorescence
   real, intent(out):: laisun(nch, nveg)    ! sunlit projected leaf area index
   real, intent(out):: laisha(nch, nveg)    ! shaded projected leaf area index
   real, intent(out):: lmrsun(nch, nveg)    ! sunlit leaf maintenance respiration rate (umol CO2/m**2/s)
   real, intent(out):: lmrsha(nch, nveg)    ! shaded leaf maintenance respiration rate (umol CO2/m**2/s)
   real, intent(out):: parabs(nch, nveg)    ! total absorbed PAR
!  local
   integer:: stomatal_model_choice 
   real ei(nch)                             ! vapor pressure inside leaf (sat vapor press at tc) (Pa)
   real ea(nch)                             ! vapor pressure of canopy air (Pa)
   real o2(nch)                             ! atmospheric o2 concentration (Pa)
   real co2(nch)                            ! atmospheric co2 concentration (Pa)
   real rb(nch)                             ! boundary layer resistance (s/m)
   real tl(nch)                             ! canopy temperature (K)
   real parsun(nch, nveg)                    ! par absorbed per unit lai (w/m**2)
   real parsha(nch, nveg)                    ! par absorbed per unit lai (w/m**2)
   real rssun(nch, nveg)                     ! sunlit stomatal resistance (s/m)
   real rssha(nch, nveg)                     ! shaded stomatal resistance (s/m)
   real vcmaxcintsun(nch, nveg)              ! leaf to canopy scaling coefficient, sunlit leaf vcmax
   real vcmaxcintsha(nch, nveg)              ! leaf to canopy scaling coefficient, shaded leaf vcmax
   real psn(nch, nveg)                       ! foliage photosynthesis (umol co2/m**2/s) [always +]
   real sif(nch, nveg)                       ! foliage fluorescence
   real lmr(nch, nveg)                       ! leaf maintenance respiration rate (umol CO2/m**2/s)
   real deldT(nch)                          ! d(es)/d(T)

   integer n, nv
   integer ivt                              ! pft vegetation type
   real*8 par                               ! photosynthetically active radiation
   real*8 cosz                              ! cosine of solar zenith angle, 0.001 <= cosz <= 1.000
   real*8 extkb                             ! direct beam extinction coefficient
   real*8 chil                              ! -0.4 <= xl <= 0.6
   real*8 gdir                              ! leaf projection in solar direction (0-1)
   real*8 omega                             ! fraction of intercepted radiation that is scattered (0 to 1)
   real*8 omegal                            ! omega for leaves
   real*8 fsun                              ! sunlit fraction of canopy
   real*8 fabd                              ! flux absorbed by canopy per unit direct flux
   real*8 fabd_sun                          ! flux absorbed by sunlit canopy per unit direct flux
   real*8 fabd_sha                          ! flux absorbed by shaded canopy per unit direct flux
   real*8 fabi                              ! flux absorbed by canopy per unit diffuse flux
   real*8 fabi_sun                          ! flux absorbed by sunlit canopy per unit diffuse flux
   real*8 fabi_sha                          ! flux absorbed by shaded canopy per unit diffuse flux
   real*8 ftdd                              ! down direct flux below canopy per unit direct flux
   real*8 ftid                              ! down diffuse flux below canopy per unit direct flux
   real*8 ftii                              ! down diffuse flux below canopy per unit diffuse flu
   real*8 rho                               ! leaf/stem reflectance weighted by fraction LAI and SAI 
   real*8 tau                               ! leaf/stem transmittance weighted by fraction LAI and SAI 
   real*8 albd                              ! surface albedo (direct)
   real*8 albi                              ! surface albedo (diffuse)
   real*8 vai                               ! LAI+SAI
   real*8 wl                                ! fraction of LAI+SAI that is LAI
   real*8 ws                                ! fraction of LAI+SAI that is SAI
   real*8 avmu                              ! average diffuse optical depth
   real*8 asu                               ! single scattering albedo
   real*8 betai                             ! upscatter parameter for diffuse radiation
   real*8 betail                            ! betai for leaves
   real*8 betad                             ! upscatter parameter for direct beam radiation
   real*8 betadl                            ! betad for leaves
   real*8 t1, rcs, rs                       ! temporary
   real*8 phi1, phi2, temp0, temp1, temp2   ! temporary
   real*8 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9  ! temporary
   real*8 b, c1, d, d1, d2, f, h, sigma, p1, p2, p3, p4, s1, s2, u1, u2, u3, a1, a2  ! temporary
   real*8 h1, h2, h3, h4, h5, h6, h7, h8, h9, h10  ! temporary

   real, intent(in):: dtc  ! canopy temperature perturbation (K) [approx 1:10000]
   real, intent(in):: dea  ! vapor pressure perturbation (Pa) [approx 1:10000]

   real, pointer:: xl(:)                   ! ecophys const-leaf/stem orientation index
   real, pointer:: rhol(:)                 ! leaf reflectance (visible)
   real, pointer:: rhos(:)                 ! stem reflectance (visible)
   real, pointer:: taul(:)                 ! leaf transmittance (visible)
   real, pointer:: taus(:)                 ! stem transmittance (visible)

   real, parameter:: omegas = 0.80         ! two-stream parameter omega for snow (visible)
   real, parameter:: betads = 0.50         ! two-stream parameter betad for snow
   real, parameter:: betais = 0.50         ! two-stream parameter betai for snow
   
   real, parameter:: mpe = 1.e-06          ! prevents overflow for division by zero
   
!  real, parameter:: extkn = 0.30          ! nitrogen allocation coefficient
   real, parameter:: extkn = 0.11          ! nitrogen allocation coefficient, fzeng changed it back to 0.11 (see Bonan et al., 2012) on 1/15/2019 


!assign the type of stomatal model that you will use. 0 indicates Ball-Berry, 1 indicates Medlyn
   stomatal_model_choice = 1

! assign local pointers to derived type arrays
! --------------------------------------------
   xl       => pftcon%xl
   rhol     => pftcon%rhol
   rhos     => pftcon%rhos
   taul     => pftcon%taul
   taus     => pftcon%taus

    o2(:) = 0.20946*pbot		       ! O2 partial pressure constant ratio
   co2(:) = co2v(:)*pbot		       ! CO2 partial pressure constant ratio [internal leaf CO2 partial pressure]
    rb(:) = 10. 			       ! gkw: for now, assume a small value for rb (see 8/3/10 email) 
    ea(:) = pbot(:) * qa(:) / (0.622+qa(:))  ! canopy air vapor pressure (Pa)

   parabs(:,:) = 0.			       ! initialize absorbed PAR to zero    

! compute saturation vapor pressure
! ---------------------------------
   do n = 1, nch
     ei(n) = MAPL_EQsat(tc(n), DQ = deldT(n))
   end do

! set solar fluxes
! ----------------
   do n = 1, nch

     par = pardir(n) + pardif(n)

! partition PAR into sunlit and shaded canopy for each vegetation type
! --------------------------------------------------------------------
     do nv = 1, nveg

       if(par > 0. .and. elai(n, nv) > 0. .and. fveg(n, nv) > 1.e-4 .and. ityp(n, nv) > 0) then

         ivt = ityp(n, nv)  ! mapped PFT

! Sun/shade big leaf code uses only one layer, with canopy integrated values from above and also canopy-integrated scaling coefficients
! -------------------------------------------------------------------------------------------------------------------------------------

! Calculate two-stream parameters omega, betad, betai, avmu, gdir, ext.
! Omega, betad, betai are adjusted for snow. Values for omega*betad
! and omega*betai are calculated and then divided by the new omega
! because the product omega*betai, omega*betad is used in solution.
! Also, the transmittances and reflectances (tau, rho) are linear
! weights of leaf and stem values.

! Weight reflectance/transmittance by lai and sai
! Only perform on vegetated pfts where coszen > 0
! -----------------------------------------------

         vai = elai(n, nv) + esai(n, nv)  ! elai+esai

         wl = elai(n, nv) / max(vai, mpe)  ! fraction of LAI+SAI that is LAI
         ws = esai(n, nv) / max(vai, mpe)  ! fraction of LAI+SAI that is SAI

         rho = max( rhol(ivt)*wl+rhos(ivt)*ws, mpe )  ! weighted reflectance, 3.11
         tau = max( taul(ivt)*wl+taus(ivt)*ws, mpe )  ! weighted transmittance, 3.12

! leaf projection in solar direction (0-1)
! CLM4.5 SurfaceAlbedoMod.F90
! ----------------------------------------

         cosz = max(0.001, coszen(n))

         chil = min(max(xl(ivt), -0.4), 0.6)
         if(abs(chil) <= 0.01) chil = 0.01
         phi1 = 0.5-0.633*chil-0.330*chil*chil
         phi2 = 0.877 * (1. - 2.*phi1)
         gdir = phi1+phi2*cosz  ! 3.3
         extkb = gdir/cosz       ! same as twostext in SurfacdAlbedoMod.F90
         avmu = (1. - phi1/phi2*log((phi1+phi2)/phi1)) / phi2  ! average diffuse optical depth, 3.4 
         temp0 = gdir+phi2*cosz
         temp1 = phi1*cosz
         temp2 = ( 1. - temp1/temp0*log((temp1+temp0)/temp1) )

         omegal = rho+tau                                  ! snow-free visible scattering coefficient 
         asu = 0.5*omegal*gdir/temp0*temp2                                  ! single scattering albedo, 3.16
         betadl = (1.+avmu*extkb)/(omegal*avmu*extkb)*asu                    ! betad for leaves, 3.15
         betail = 0.5 * ((rho+tau) + (rho-tau) * ((1.+chil)/2.)**2) / omegal  ! betai for leaves, 3.13, 3.14

! Adjust omega, betad, and betai for intercepted snow
! CLM4.5 SurfaceAlbedoMod.F90
! ---------------------------------------------------
	 
	 if (tc(n) > tfrz) then	                       ! no snow
 	   tmp0 = omegal       ! 3.8 
           tmp1 = betadl       ! 3.10
           tmp2 = betail       ! 3.9
	  else
           tmp0 =  (1.-fwet(n))*omegal	      + fwet(n)*omegas                 ! 3.5
           tmp1 = ((1.-fwet(n))*omegal*betadl+fwet(n)*omegas*betads) / tmp0  ! 3.7
           tmp2 = ((1.-fwet(n))*omegal*betail+fwet(n)*omegas*betais) / tmp0  ! 3.6
	 end if
         omega = tmp0	      
         betad = tmp1                                  ! upscatter parameter for direct beam radiation
         betai = tmp2                                  ! upscatter parameter for diffuse radiation

! Common terms
! CLM4.5 SurfaceAlbedoMod.F90
! ---------------------------

         b = 1. - omega+omega*betai  ! 3.31
         c1 = omega*betai             ! 3.32
         tmp0 = avmu*extkb
         d = tmp0*omega*betad       ! 3.33
         f = tmp0*omega*(1.-betad)  ! 3.34
         tmp1 = b*b - c1*c1
         h = sqrt(tmp1) / avmu        ! 3.35
         sigma = tmp0*tmp0-tmp1     ! 3.36
         p1 = b+avmu*h              ! 3.42
         p2 = b-avmu*h              ! 3.43
         p3 = b+tmp0                ! 3.44
         p4 = b-tmp0                ! 3.45

! Absorbed, reflected, transmitted fluxes per unit incoming radiation for full canopy
! CLM4.5 SurfaceAlbedoMod.F90
! -----------------------------------------------------------------------------------
         
         t1 = min(h*vai, 40.)
         s1 = exp(-t1)                ! 3.40
         t1 = min(extkb*vai, 40.)
         s2 = exp(-t1)                ! 3.41

! Direct beam
! CLM4.5 SurfaceAlbedoMod.F90 
! ---------------------------
         
         u1 = b-c1/albgrd(n, nv)     ! 3.37
         u2 = b-c1*albgrd(n, nv)     ! 3.38
         u3 = f+c1*albgrd(n, nv)     ! 3.39

         tmp2 = u1-avmu*h
         tmp3 = u1+avmu*h
         d1 = p1*tmp2/s1-p2*tmp3*s1  ! 3.46
         tmp4 = u2+avmu*h
         tmp5 = u2-avmu*h
         d2 = tmp4/s1-tmp5*s1       ! 3.47
         h1 = -d*p4-c1*f            ! 3.48

         tmp6 = d-h1*p3/sigma
         tmp7 = (d-c1-h1/sigma*(u1+tmp0)) * s2
         h2 = (tmp6*tmp2/s1-p2*tmp7) / d1   ! 3.49
         h3 = -(tmp6*tmp3*s1-p1*tmp7) / d1  ! 3.50
         h4 = -f*p3-c1*d                    ! 3.51
         tmp8 = h4/sigma
         tmp9 = (u3-tmp8*(u2-tmp0)) * s2
         h5 = -(tmp8*tmp4/s1+tmp9) / d2     ! 3.52
         h6 = (tmp8*tmp5*s1+tmp9) / d2      ! 3.53
	 
	 albd = h1/sigma+h2+h3            ! Flux reflected by vegetation (direct), 3.17
         ftid = h4*s2/sigma+h5*s1+h6/s1   ! Downward diffuse flux below vegetation, 3.19
         ftdd = s2                            ! Downward direct flux below vegetation
         fabd = 1. - albd - (1.-albgrd(n, nv))*ftdd - (1.-albgri(n, nv))*ftid   ! Flux absorbed by vegetation (direct), 3.21

         a1 = h1/sigma * (1. - s2*s2) / (2. * extkb) &                ! 3.25 
            + h2	 * (1. - s2*s1) / (extkb+h) &
            + h3	 * (1. - s2/s1) / (extkb-h)

         a2 = h4/sigma * (1. - s2*s2) / (2. * extkb) &                ! 3.26 
            + h5	 * (1. - s2*s1) / (extkb+h) &
            + h6	 * (1. - s2/s1) / (extkb-h)

         fabd_sun = (1. - omega) * ( 1. - s2+1. / avmu * (a1+a2) )  ! Flux absorbed by vegetation (direct), sunlit, 3.23
         fabd_sha = fabd-fabd_sun                                     ! Flux absorbed by vegetation (direct), shaded, 3.24
	 
! Diffuse
! CLM4.5 SurfaceAlbedoMod.F90
! ---------------------------

	 u1 = b-c1/albgri(n, nv)
	 u2 = b-c1*albgri(n, nv)
	 tmp2 = u1-avmu*h
	 tmp3 = u1+avmu*h
	 d1 = p1*tmp2/s1-p2*tmp3*s1
	 tmp4 = u2+avmu*h
	 tmp5 = u2-avmu*h
	 d2 = tmp4/s1-tmp5*s1	 
         h7 = (c1*tmp2) / (d1*s1)             ! 3.54
         h8 = (-c1*tmp3*s1) / d1              ! 3.55
         h9 = tmp4 / (d2*s1)                  ! 3.56
         h10 = (-tmp5*s1) / d2                ! 3.57
	 
	 albi = h7+h8                             ! Flux reflected by vegetation (indirect), 3.18
	 ftii = h9*s1+h10/s1                      ! Downward direct and diffuse fluxes below vegetation, 3.20
	 fabi = 1. - albi - (1.-albgri(n, nv))*ftii  ! Flux absorbed by vegetation (indirect), 3.22

	 a1 = h7 * (1. - s2*s1) / (extkb+h) +  h8 * (1. - s2/s1) / (extkb-h)   ! 3.29
	 a2 = h9 * (1. - s2*s1) / (extkb+h) + h10 * (1. - s2/s1) / (extkb-h)   ! 3.30

	 fabi_sun = (1. - omega) / avmu * (a1+a2)  ! 3.27
	 fabi_sha = fabi-fabi_sun                 ! 3.28

! PAR absorbed by vegetation
! --------------------------
!        parabs(n) = parabs(n) + (pardir(n)*fabd+pardif(n)*fabi)*fveg(n, nv)*wl  ! save leaf-absorbed PAR for FPAR calculation

! sunlit fraction of canopy
! -------------------------
         fsun = (1. - s2) / t1
	 
!         if(fsun <= 0. .or. fsun >= 1.) then
!           print *, 'elai =',elai(n, nv), 'fsun =',fsun
!           stop 'compute_rc: fsun out of bound!!'
!         end if

! leaf to canopy scaling coefficients. Need to separate for nlevcan == 1 and nlevcan > 1 cases!!
! See L807-818 in SurfaceAlbedoMod.F90. This is default though.   
! -------------------------------------------------------------

         vcmaxcintsun(n, nv) = (1. - exp(-(extkn+extkb)*elai(n, nv))) / (extkn+extkb)
         vcmaxcintsha(n, nv) = (1. - exp(-extkn*elai(n, nv))) / extkn-vcmaxcintsun(n, nv)

         if(elai(n, nv) .gt. 0.01) then      

           ! absorbed PAR (per unit sun/shade lai+sai)
           ! -----------------------------------------	 
	   fabd_sun = fabd_sun / (fsun*vai)
	   fabi_sun = fabi_sun / (fsun*vai)
	   fabd_sha = fabd_sha / ((1. - fsun)*vai)
	   fabi_sha = fabi_sha / ((1. - fsun)*vai)	      

           ! sunlit & shaded leaf area
           ! -------------------------
           laisun(n, nv) = elai(n, nv)*fsun
	   laisha(n, nv) = elai(n, nv)*(1.-fsun)

           vcmaxcintsun(n, nv) = vcmaxcintsun(n, nv) / laisun(n, nv)
           vcmaxcintsha(n, nv) = vcmaxcintsha(n, nv) / laisha(n, nv)
         
         else   ! special case for low elai, when fsun can be > 1.0

           ! absorbed PAR (per unit sun/shade lai+sai)
           ! -----------------------------------------	 
	   fabd_sun = fabd_sun/vai
	   fabi_sun = fabi_sun/vai
	   fabd_sha = 0.
	   fabi_sha = 0.	      

           ! sunlit & shaded leaf area
           ! -------------------------
           laisun(n, nv) = elai(n, nv)
	   laisha(n, nv) = 0.

           vcmaxcintsun(n, nv) = vcmaxcintsun(n, nv) / laisun(n, nv)
           vcmaxcintsha(n, nv) = 0.
         
         endif
         
         ! CLM4.5 SurfaceRadiationMod.F90, L453-458
         parsun(n, nv) = pardir(n)*fabd_sun+pardif(n)*fabi_sun   ! sunlit canopy PAR for leaves per vai
         parsha(n, nv) = pardir(n)*fabd_sha+pardif(n)*fabi_sha   ! shaded canopy PAR for leaves per vai

        else

! elai = 0, no vegetation, or PAR = 0 (night)
! ---------------------------------------
         laisun(n, nv) = 0.
         laisha(n, nv) = elai(n, nv)
         parsun(n, nv) = 0.
         parsha(n, nv) = 0.
         vcmaxcintsun(n, nv) = 0.
         if(ityp(n, nv) > 0) then
           vcmaxcintsha(n, nv) = (1. - exp(-extkn*elai(n, nv))) / extkn
           if(elai(n, nv) > 0.) then
             vcmaxcintsha(n, nv) = vcmaxcintsha(n, nv) / elai(n, nv)
            else
             vcmaxcintsha(n, nv) = 0.
           end if
          else
           vcmaxcintsha(n, nv) = 0.
         endif
       endif

       if(present (fpar_sf)) then 

          ! scaling to match MODIS FPAR

          parsun(n, nv) = parsun(n, nv) * fpar_sf(n, nv)  ! sunlit canopy PAR for leaves per vai
          parsha(n, nv) = parsha(n, nv) * fpar_sf(n, nv)  ! shaded canopy PAR for leaves per vai

       endif

!      save leaf-absorbed PAR for FPAR calculation
       parabs(n, nv) =  parsun(n, nv) * laisun(n, nv) + parsha(n, nv) * laisha(n, nv)
       
     end do  ! end PFT loop

   end do   ! end column loop

!  compute stomatal resistance using CLM routine; also compute photosynthesis

!  obtain stomatal resistance and photosynthesis

   call Photosynthesis(nch, nveg, ei, ea, o2, co2, rb, dayl_fac, pbot, ityp, tm, t10, tc,  &
		      btran, elai, laisun, parsun, vcmaxcintsun, rssun, psnsun, sifsun, lmrsun, stomatal_model_choice)  ! sunlit
   call Photosynthesis(nch, nveg, ei, ea, o2, co2, rb, dayl_fac, pbot, ityp, tm, t10, tc,  &
		      btran, elai, laisha, parsha, vcmaxcintsha, rssha, psnsha, sifsha, lmrsha, stomatal_model_choice)  ! shaded

!  combine resistance as reciprocal of vegetation weighted conductance (weighted harmonic sum)

!DIR$ NOVECTOR
   do n = 1, nch
     rcs = 0.
     do nv = 1, nveg
       rs = laisun(n, nv)/rssun(n, nv) + laisha(n, nv)/rssha(n, nv)  ! rssun and rssha: stomatal resistance; rs: stomatal conductance, fzeng, 20 Feb 2018
       rcs = rcs+fveg(n, nv)*rs
     end do 
     rc(n) = 1./max(rcs, 5.e-5) + rb(n)    ! rc: stomatal resistance, fzeng, 20 Feb 2018
   end do


!  compute resistance with small delta ea

   ea(:) = ea(:) + dea

   call Photosynthesis(nch, nveg, ei, ea, o2, co2, rb, dayl_fac, pbot, ityp, tm, t10, tc,  &
		      btran, elai, laisun, parsun, vcmaxcintsun, rssun, psn, sif, lmr, stomatal_model_choice)  ! sunlit
   call Photosynthesis(nch, nveg, ei, ea, o2, co2, rb, dayl_fac, pbot, ityp, tm, t10, tc,  &
		      btran, elai, laisha, parsha, vcmaxcintsha, rssha, psn, sif, lmr, stomatal_model_choice)  ! shaded

!DIR$ NOVECTOR
   do n = 1, nch
     rcs = 0.
     do nv = 1, nveg
       rs = laisun(n, nv)/rssun(n, nv) + laisha(n, nv)/rssha(n, nv)
       rcs = rcs+fveg(n, nv)*rs
     end do 
     rcdea(n) = 1./max(rcs, 5.e-5) + rb(n)
   end do


!  compute resistance with small delta Tc

   tl(:) = tc(:) + dtc
   ei(:) = ei(:) + deldT(:)*dtc  ! ei = esat(Tc)+[d(esat)/d(Tc)]dTc

   ea(:) = pbot(:) * qa(:) / (0.622+qa(:))  ! reset input canopy air vapor pressure (Pa)

   call Photosynthesis(nch, nveg, ei, ea, o2, co2, rb, dayl_fac, pbot, ityp, tm, t10, tl,  &
		      btran, elai, laisun, parsun, vcmaxcintsun, rssun, psn, sif, lmr, stomatal_model_choice)  ! sunlit
   call Photosynthesis(nch, nveg, ei, ea, o2, co2, rb, dayl_fac, pbot, ityp, tm, t10, tl,  &
		      btran, elai, laisha, parsha, vcmaxcintsha, rssha, psn, sif, lmr, stomatal_model_choice)  ! shaded

!DIR$ NOVECTOR
   do n = 1, nch
     rcs = 0.
     do nv = 1, nveg
       rs = laisun(n, nv)/rssun(n, nv) + laisha(n, nv)/rssha(n, nv)
       rcs = rcs+fveg(n, nv)*rs
     end do 
     rcdtc(n) = 1./max(rcs, 5.e-5) + rb(n)
   end do

   end subroutine compute_rc
   
!*******************************************************************************
!  fzeng: this is how "btran" and "btran2" are calculated in CanopyFluxes

!!!  real(r8), parameter:: btran0 = 0.0_r8  ! initial value

!!!   ! Initialize
!!!
!!!   do f = 1, fn
!!!      p = filterp(f)
!!!      btran(p)  = btran0
!!!      btran2(p)  = btran0
!!!   end do
!!!
!!!   ! Effective porosity of soil, partial volume of ice and liquid (needed for btran)
!!!   ! and root resistance factors
!!!
!!!   do j = 1, nlevgrnd
!!!      do f = 1, fn
!!!         p = filterp(f)
!!!         c = pcolumn(p)
!!!         l = plandunit(p)
!!!
!!!         ! Root resistance factors
!!!
!!!         vol_ice = min(watsat(c, j), h2osoi_ice(c, j)/(dz(c, j)*denice))
!!!         eff_porosity = watsat(c, j)-vol_ice
!!!         vol_liq = min(eff_porosity, h2osoi_liq(c, j)/(dz(c, j)*denh2o))
!!!         if (vol_liq .le. 0._r8 .or. t_soisno(c, j) .le. tfrz-2._r8) then
!!!            rootr(p, j) = 0._r8
!!!         else
!!!            s_node = max(vol_liq/eff_porosity, 0.01_r8)
!!!            smp_node = max(smpsc(ivt(p)), -sucsat(c, j)*s_node**(-bsw(c, j)))
!!!
!!!            rresis(p, j) = min( (eff_porosity/watsat(c, j))* &
!!!                          (smp_node-smpsc(ivt(p))) / (smpso(ivt(p)) - smpsc(ivt(p))), 1._r8)
!!!            if (.not. (perchroot .or. perchroot_alt) ) then
!!!                rootr(p, j) = rootfr(p, j)*rresis(p, j)
!!!            else
!!!               rootr(p, j) = rootfr_unf(p, j)*rresis(p, j)
!!!            end if
!!!            btran(p)    = btran(p) + rootr(p, j)
!!!            smp_node_lf = max(smpsc(ivt(p)), -sucsat(c, j)*(h2osoi_vol(c, j)/watsat(c, j))**(-bsw(c, j))) 
!!!            btran2(p)   = btran2(p) +rootfr(p, j)*min((smp_node_lf-smpsc(ivt(p))) / (smpso(ivt(p)) - smpsc(ivt(p))), 1._r8)
!!!         endif 
!!!      end do
!!!   end do
!*******************************************************************************

!*******************************************************************************
! gkw: this is how rb is calculated in CanopyFluxes
! fzeng: no change from CLM4 to CLM4.5

! Determine aerodynamic resistances & Bulk boundary layer resistance of leaves
! gkw: code imported from CanopyFluxes for later use 8/3/2010

!!!  do n = 1, nch
!!!    ram1 = 1./(ustar(n)*ustar(n)/um(n))            ! gkw: 5.55
!!!    uaf  = um(p)*sqrt( 1./(ram1*um(n)) )           ! gkw: 5.100
!!!    rb(n) = 100.*sqrt(dleaf(ivt))/sqrt(uaf)        ! gkw: 5.109
!!!  end do
!******************************************************************************* 

   subroutine Photosynthesis (nch, nveg, esat_tv, eair, oair, cair, rb, dayl_factor, forc_pbot, ityp, tgcm, t10, t_veg,  &
			      btran_in, tlai, lai, apar, vcmaxcint, rs, psn, sif, lmr, stomatal_model_choice)
!
! !DESCRIPTION:
! Leaf photosynthesis and stomatal conductance calculation as described by
! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 and extended to
! a multi-layer canopy
!
! !REVISION HISTORY:

! !USES:
   use MAPL
   use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
   use clmtype
   use clm_varcon  , only : rgas, tfrz
   use clm_varpar  , only : nlevcan
   use pftvarcon   , only : nbrdlf_dcd_tmp_shrub
   use pftvarcon   , only : nsoybean, nsoybeanirrig, npcropmin
!#if (defined CN)
   use CNAllocationMod, only : CNAllocation_Carbon_only
   use pso_params
   use compute_rc_functions, only: pft_clm_to_pso
   use read_env, only : read_env_data
!#endif
   implicit none

! !ARGUMENTS:
   integer, intent(in)    :: nch	           ! number of land tiles
   integer, intent(in)    :: nveg	           ! number of PFTs   
   real(r8), intent(in)    :: esat_tv(nch)         ! saturation vapor pressure at t_veg (Pa)
   real(r8), intent(in)    :: eair(nch)            ! vapor pressure of canopy air (Pa)
   real(r8), intent(in)    :: oair(nch)            ! Atmospheric O2 partial pressure (Pa)
   real(r8), intent(in)    :: cair(nch)            ! Atmospheric CO2 partial pressure (Pa)
   real(r8), intent(in)    :: rb(nch)              ! boundary layer resistance (s/m)
   real(r8), intent(in)    :: dayl_factor(nch)     ! scalar (0-1) for daylength
   real(r8), intent(in)    :: forc_pbot(nch)       ! atmospheric pressure (Pa)
   integer,  intent(in)    :: ityp(nch, nveg)       ! vegetation type
   real(r8), intent(in)    :: tgcm(nch) 	   ! air temperature at agcm reference height (kelvin)
!KO
   real(r8), intent(in)    :: t10(nch)             ! 10-day running mean of the 2 m temperature (K)
!KO     
   real(r8), intent(in)    :: t_veg(nch)	   ! vegetation temperature (Kelvin)   
   
   real(r8), intent(in)    :: btran_in(nch)        ! soil water transpiration factor (0 to 1)

! gkw & fzeng:   
   real(r8), intent(in)    :: tlai(nch, nveg)       ! total leaf area index
   real(r8), intent(in)    :: lai(nch, nveg)        ! leaf area index for canopy layer, sunlit or shaded
   real(r8), intent(in)    :: apar(nch, nveg)	   ! par absorbed per unit lai (w/m**2)
   real(r8), intent(in)    :: vcmaxcint(nch, nveg)  ! leaf to canopy scaling coefficient   
   integer,  intent(in)    :: stomatal_model_choice  ! 0 for Ball-Berry Model, 1 for Medlyn Model
   real(r8), intent(out)   :: rs(nch, nveg)	   ! leaf stomatal resistance (s/m)
   real(r8), intent(out)   :: psn(nch, nveg)	   ! foliage photosynthesis (umol co2/m**2/s) [always +]
   real(r8), intent(out)   :: sif(nch, nveg)	   ! foliage fluorescence
   real(r8), intent(out)   :: lmr(nch, nveg)	   ! leaf maintenance respiration rate (umol CO2/m**2/s)
   
! !CALLED FROM:
! subroutine CanopyFluxes in this module

! !LOCAL VARIABLES:
!
! local pointers to implicit in variables
   real(r8), pointer:: c3psn(:)                   ! photosynthetic pathway: 0. = c4, 1. = c3
   real(r8), pointer:: slatop(:)                  ! specific leaf area at top of canopy, projected area basis [m^2/gC]
   real(r8), pointer:: flnr(:)                    ! fraction of leaf N in the Rubisco enzyme (gN Rubisco/gN leaf)
   real(r8), pointer:: fnitr(:)                   ! foliage nitrogen limitation factor (-)
   real(r8), pointer:: leafcn(:)                  ! leaf C:N (gC/gN)

   integer  :: nrad(nch)	                   ! number of canopy layers, above snow for radiative transfer
   real(r8):: tlai_z(nch, nlevcan)                 ! total leaf area index for canopy layer
   real(r8):: lai_z(nch, nlevcan)                  ! leaf area index for canopy layer, sunlit or shaded
   real(r8):: par_z(nch, nlevcan)                  ! par absorbed per unit lai for canopy layer (w/m**2)

!! fzeng: comment out C13 for now. Add it back later!!
   !!! C13
!   real(r8), pointer:: alphapsn(nch)              ! 13C fractionation factor for PSN ()

! local pointers to implicit out variables
   real(r8):: psn_z(nch, nlevcan)                  ! canopy layer: foliage photosynthesis (umol co2/m**2/s) [always +]
   real(r8):: lmr_z(nch, nlevcan)                  ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
   real(r8):: rs_z(nch, nlevcan)                   ! canopy layer: leaf stomatal resistance (s/m)
   real(r8):: ci_z(nch, nlevcan)                   ! intracellular leaf CO2 (Pa)
   real(r8):: psn_wc(nch)	                   ! Rubisco-limited foliage photosynthesis (umol co2/m**2/s) [always +]
   real(r8):: psn_wj(nch)	                   ! RuBP-limited foliage photosynthesis (umol co2/m**2/s) [always +]
   real(r8):: psn_wp(nch)	                   ! product-limited foliage photosynthesis (umol co2/m**2/s) [always +]

!KO
   real(r8):: rh_leaf(nch)                        ! fractional humidity at leaf surface (dimensionless)
!KO

! Leaf photosynthesis parameters
   real(r8):: vcmax_z(nch, nlevcan)                ! maximum rate of carboxylation (umol co2/m**2/s)
   real(r8):: jmax_z(nch, nlevcan)                 ! maximum electron transport rate (umol electrons/m**2/s)
   real(r8):: tpu_z(nch, nlevcan)                  ! triose phosphate utilization rate (umol CO2/m**2/s)
   real(r8):: kp_z(nch, nlevcan)                   ! initial slope of CO2 response curve (C4 plants)

   logical  :: c3flag(nch)                         ! true if C3 and false if C4
   real(r8):: lnc(nch)                            ! leaf N concentration (gN leaf/m^2)
   real(r8):: kc(nch)                             ! Michaelis-Menten constant for CO2 (Pa)
   real(r8):: ko(nch)                             ! Michaelis-Menten constant for O2 (Pa)
   real(r8):: cp(nch)                             ! CO2 compensation point (Pa)
   real(r8):: bbbopt(nch)                         ! Ball-Berry minimum leaf conductance, unstressed (umol H2O/m**2/s)
   real(r8):: bbb(nch)                            ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
   real(r8):: mbbopt(nch)                         ! Ball-Berry slope of conductance-photosynthesis relationship, unstressed
   real(r8):: mbb(nch)                            ! Ball-Berry slope of conductance-photosynthesis relationship
   real(r8):: kn(nch)                             ! leaf nitrogen decay coefficient
   real(r8):: btran(nch)                          ! transpiration coefficient
   real(r8):: vcmax25top                          ! canopy top: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
   real(r8):: jmax25top                           ! canopy top: maximum electron transport rate at 25C (umol electrons/m**2/s)
   real(r8):: tpu25top                            ! canopy top: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
   real(r8):: lmr25top                            ! canopy top: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
   real(r8):: kp25top                             ! canopy top: initial slope of CO2 response curve (C4 plants) at 25C

   real(r8):: vcmax25                             ! leaf layer: maximum rate of carboxylation at 25C (umol CO2/m**2/s)
   real(r8):: jmax25                              ! leaf layer: maximum electron transport rate at 25C (umol electrons/m**2/s)
   real(r8):: tpu25                               ! leaf layer: triose phosphate utilization rate at 25C (umol CO2/m**2/s)
   real(r8):: lmr25                               ! leaf layer: leaf maintenance respiration rate at 25C (umol CO2/m**2/s)
   real(r8):: kp25                                ! leaf layer: Initial slope of CO2 response curve (C4 plants) at 25C
   real(r8):: kc25                                ! Michaelis-Menten constant for CO2 at 25C (Pa)
   real(r8):: ko25                                ! Michaelis-Menten constant for O2 at 25C (Pa)
   real(r8):: cp25                                ! CO2 compensation point at 25C (Pa)

   real(r8):: vcmaxha                             ! activation energy for vcmax (J/mol)
   real(r8):: jmaxha                              ! activation energy for jmax (J/mol)
   real(r8):: tpuha                               ! activation energy for tpu (J/mol)
   real(r8):: lmrha                               ! activation energy for lmr (J/mol)
   real(r8):: kcha                                ! activation energy for kc (J/mol)
   real(r8):: koha                                ! activation energy for ko (J/mol)
   real(r8):: cpha                                ! activation energy for cp (J/mol)

   real(r8):: vcmaxhd                             ! deactivation energy for vcmax (J/mol)
   real(r8):: jmaxhd                              ! deactivation energy for jmax (J/mol)
   real(r8):: tpuhd                               ! deactivation energy for tpu (J/mol)
   real(r8):: lmrhd                               ! deactivation energy for lmr (J/mol)

   real(r8):: vcmaxse                             ! entropy term for vcmax (J/mol/K)
   real(r8):: jmaxse                              ! entropy term for jmax (J/mol/K)
   real(r8):: tpuse                               ! entropy term for tpu (J/mol/K)
   real(r8):: lmrse                               ! entropy term for lmr (J/mol/K)

   real(r8):: vcmaxc                              ! scaling factor for high temperature inhibition (25 C = 1.0)
   real(r8):: jmaxc                               ! scaling factor for high temperature inhibition (25 C = 1.0)
   real(r8):: tpuc                                ! scaling factor for high temperature inhibition (25 C = 1.0)
   real(r8):: lmrc                                ! scaling factor for high temperature inhibition (25 C = 1.0)

   real(r8):: qe(nch)                             ! quantum efficiency, used only for C4 (mol CO2/mol photons)
   real(r8):: fnps                                ! fraction of light absorbed by non-photosynthetic pigments
   real(r8):: theta_psii                          ! empirical curvature parameter for electron transport rate

   real(r8):: theta_cj(nch)                       ! empirical curvature parameter for ac, aj photosynthesis co-limitation
   real(r8):: theta_ip                            ! empirical curvature parameter for ap photosynthesis co-limitation

! Other
   integer  :: n, p, g, iv, nv, ivt                     ! indices
   real(r8):: cf                                  ! s m**2/umol -> s/m
   real(r8):: rsmax0                              ! maximum stomatal resistance [s/m]
   real(r8):: gb                                  ! leaf boundary layer conductance (m/s)
   real(r8):: cs                                  ! CO2 partial pressure at leaf surface (Pa)
   real(r8):: gs                                  ! leaf stomatal conductance (m/s)
   real(r8):: hs                                  ! fractional humidity at leaf surface (dimensionless)
   real(r8):: sco                                 ! relative specificity of rubisco
   real(r8):: ft				   ! photosynthesis temperature response (statement function)
   real(r8):: fth				   ! photosynthesis temperature inhibition (statement function)
   real(r8):: fth25				   ! ccaling factor for photosynthesis temperature inhibition (statement function)
   real(r8):: tl                                  ! leaf temperature in photosynthesis temperature function (K)
   real(r8):: ha                                  ! activation energy in photosynthesis temperature function (J/mol)
   real(r8):: hd                                  ! deactivation energy in photosynthesis temperature function (J/mol)
   real(r8):: se                                  ! entropy term in photosynthesis temperature function (J/mol/K)
   real(r8):: cc                                  ! scaling factor for high temperature inhibition (25 C = 1.0)
   real(r8):: ciold                               ! previous value of Ci for convergence check
   real(r8):: gs_mol_err                          ! gs_mol for error check
   real(r8):: je(nch)                             ! electron transport rate (umol electrons/m**2/s)
   real(r8):: qabs                                ! PAR absorbed by PS II (umol photons/m**2/s)
   real(r8):: aquad, bquad, cquad                   ! terms for quadratic equations
   real(r8):: r1, r2                               ! roots of quadratic equation
   real(r8):: ceair                               ! vapor pressure of air, constrained (Pa)
   real(r8):: fnr                                 ! (gRubisco/gN in Rubisco)
   real(r8):: act25                               ! (umol/mgRubisco/min) Rubisco activity at 25 C
   integer  :: niter                               ! iteration loop index
   real(r8):: nscaler                             ! leaf nitrogen scaling coefficient

   real(r8):: ac(nch, nlevcan)                     ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
   real(r8):: aj(nch, nlevcan)                     ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
   real(r8):: ap(nch, nlevcan)                     ! product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
   real(r8):: ag(nch, nlevcan)                     ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
   real(r8):: an(nch, nlevcan)                     ! net leaf photosynthesis (umol CO2/m**2/s)
   real(r8):: gs_mol(nch, nlevcan)                 ! leaf stomatal conductance (umol H2O/m**2/s)
   real(r8):: gb_mol(nch)	                   ! leaf boundary layer conductance (umol H2O/m**2/s)

   real(r8):: psn_wc_z(nch, nlevcan)               ! Rubisco-limited contribution to psn_z (umol CO2/m**2/s)
   real(r8):: psn_wj_z(nch, nlevcan)               ! RuBP-limited contribution to psn_z (umol CO2/m**2/s)
   real(r8):: psn_wp_z(nch, nlevcan)               ! product-limited contribution to psn_z (umol CO2/m**2/s)

   real(r8):: psncan                              ! canopy sum of psn_z
   real(r8):: psncan_wc                           ! canopy sum of psn_wc_z
   real(r8):: psncan_wj                           ! canopy sum of psn_wj_z
   real(r8):: psncan_wp                           ! canopy sum of psn_wp_z
   real(r8):: lmrcan                              ! canopy sum of lmr_z
   real(r8):: gscan                               ! canopy sum of leaf conductance
   real(r8):: laican                              ! canopy sum of lai_z
   real(r8):: rh_can

! gkw & fzeng: added for SIF; taken from stomata, may not correct here
! --------------------------------------------------------------------
   real(r8):: cican                               ! canopy mean intracellular leaf CO2 (Pa)
   real(r8):: ppf                                 ! absorb photosynthetic photon flux (umol photons/m**2/s)
   real(r8):: j                                   ! electron transport (umol co2/m**2/s)
   real(r8):: je_sif                              ! actual electron transport
   real(r8):: xn                                  ! je/j
   real(r8):: fs                                  ! fluorescnce yield at Fs

   integer:: iulog = 6
   integer:: this_particle  ! current ensemble number for running PSO optimization
   integer:: g1_ef_choice
   logical:: pso_alloc_stat
   integer:: num_params
   integer:: fake_num_ensemble
   integer:: reason
   integer:: pso_pft
   integer:: g0_pso_pft
   real   :: pso_val, offline_const, offline_int, g1_pso_val, g0_pso_val
   real   :: pso_intercept, pso_slope
   real,dimension(1) :: map_val_tile
   real :: map_val_tile_real
   logical:: pso_exists
   character (len = 100):: met_tag_map
   integer :: EXP_ID
   integer :: STATUS
   integer :: this_tile
   integer, dimension(1) :: this_map_idx


!------------------------------------------------------------------------------

   ! Temperature and soil water response functions

   ft(tl, ha) = exp( ha / (rgas*1.e-3_r8*(tfrz+25._r8)) * (1._r8 - (tfrz+25._r8)/tl) )
   fth(tl, hd, se, cc) = cc / ( 1._r8+exp( (-hd+se*tl) / (rgas*1.e-3_r8*tl) ) )
   fth25(hd, se) = 1._r8+exp( (-hd+se*(tfrz+25._r8)) / (rgas*1.e-3_r8*(tfrz+25._r8)) )

   ! fzeng: in CLM4.5 BiogeophysRestMod.F90 nrad is set to nlevcan (set to 1 in CLM4.5 clm_varpar.F90) if it's not in the restart file.
   
   if (nlevcan == 1) then 
     nrad(:) = nlevcan   ! gkw: nlevcan, unless it's buried by snow; could compare elai and tlai. fzeng: need modification if nlevcan > 1
   else 
     stop 'compute_rc: nlevcan not equals 1'
   endif
   btran = btran_in    ! gkw & fzeng: make local copy, because btran is modified in this routine. Confirmed by Randy. 

   ! Assign local pointers to pft constants
   
   c3psn     => pftcon%c3psn
   leafcn    => pftcon%leafcn
   flnr      => pftcon%flnr
   fnitr     => pftcon%fnitr
   slatop    => pftcon%slatop
   
   !==============================================================================!
   ! Photosynthesis and stomatal conductance parameters, from:
   ! Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
   !==============================================================================!

   ! vcmax25 parameters, from CN

   fnr = 7.16_r8
   act25 = 3.6_r8   !umol/mgRubisco/min
   ! Convert rubisco activity units from umol/mgRubisco/min -> umol/gRubisco/s
   act25 = act25*1000.0_r8/60.0_r8

   ! Activation energy, from:
   ! Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
   ! Bernacchi et al (2003) Plant, Cell and Environment 26:1419-1430
   ! except TPU from: Harley et al (1992) Plant, Cell and Environment 15:271-282

   kcha    = 79430._r8
   koha    = 36380._r8
   cpha    = 37830._r8
!KO   vcmaxha = 65330._r8
!KO   jmaxha  = 43540._r8
!KO   tpuha   = 53100._r8
!KO
   vcmaxha = 72000._r8
   jmaxha  = 50000._r8
   tpuha   = 72000._r8
!KO
   lmrha   = 46390._r8

   ! High temperature deactivation, from:
   ! Leuning (2002) Plant, Cell and Environment 25:1205-1210
   ! The factor "c" scales the deactivation to a value of 1.0 at 25C

!KO   vcmaxhd = 149250._r8
!KO   jmaxhd  = 152040._r8
!KO   tpuhd   = 150650._r8
!KO
   vcmaxhd = 200000._r8
   jmaxhd  = 200000._r8
   tpuhd   = 200000._r8
!KO
   lmrhd   = 150650._r8

!KO   vcmaxse = 485._r8
!KO   jmaxse  = 495._r8
!KO   tpuse   = 490._r8
   lmrse   = 490._r8

!KO   vcmaxc = fth25 (vcmaxhd, vcmaxse)
!KO   jmaxc  = fth25 (jmaxhd, jmaxse)
!KO   tpuc   = fth25 (tpuhd, tpuse)
   lmrc   = fth25 (lmrhd, lmrse)

   ! Miscellaneous parameters, from Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593

   fnps = 0.15_r8
   theta_psii = 0.7_r8
   theta_ip = 0.95_r8
   
   do nv = 1, nveg      ! gkw: loop over the four vegetation types
   
   do n = 1, nlevcan    ! fzeng: need modification if nlevcan > 1
     tlai_z(:,n) = tlai(:,nv)  ! gkw: tlai
     lai_z(:,n)  = lai(:,nv)  ! gkw: lai
     par_z(:,n)  = apar(:,nv)  ! gkw: par
   end do
   
   do p = 1, nch        ! gkw: loop over tiles

      if(ityp(p, nv) > 0) then
      
      g = p            ! fzeng: do so to allow fewer modifications from the original CLM4.5 code
      ivt = ityp(p, nv)  ! mapped vegetation type into CLM PFT

      ! Modification for shrubs proposed by X.D.Z
      ! Why does he prefer this line here instead of in subr.
      ! CanopyFluxes? (slevis)
      ! Equivalent modification for soy following AgroIBIS
#if (defined CNDV)
      if (ivt == nbrdlf_dcd_tmp_shrub .or. ivt == nbrdlf_dcd_tmp_shrub2)  btran(p) = min(1._r8, btran(p) * 3.33_r8)  ! gkw: should we do this for the seasonal deciduous split type?
#endif
      if (ivt == nsoybean .or. ivt == nsoybeanirrig) btran(p) = min(1._r8, btran(p) * 1.25_r8)
           
      ! C3 or C4 photosynthesis logical variable

      if (nint(c3psn(ivt)) == 1) then
         c3flag(p) = .true. 
      else if (nint(c3psn(ivt)) == 0) then
         c3flag(p) = .false.
      end if

      ! C3 and C4 dependent parameters
      g1_ef_choice = 0
      if (g1_ef_choice == 0) then
         if (c3flag(p)) then
            qe(p) = 0._r8
            theta_cj(p) = 0.98_r8
            bbbopt(p) = 10000._r8
            mbbopt(p) = 4
         else
            qe(p) = 0.05_r8
            theta_cj(p) = 0.80_r8
            bbbopt(p) = 40000._r8
            mbbopt(p) = 4
         endif
      elseif (g1_ef_choice == 1) then
         if (c3flag(p)) then
            qe(p) = 0._r8
            theta_cj(p) = 0.98_r8
            bbbopt(p) = 10000._r8
         else
            qe(p) = 0.05_r8
            theta_cj(p) = 0.80_r8
            bbbopt(p) = 40000._r8
         end if
          
         this_particle = pso_vals%particle_num + 1
         this_tile = pso_vals%local_tile_nums(p)
         this_map_idx = findloc(pso_vals%all_tile_nums,VALUE=this_tile)
         map_val_tile = pso_vals%map_vals(this_map_idx)
         map_val_tile_real = map_val_tile(1)

         !!!!! UNCOMMENT BELOW FOR PFT-BASED PSO OPTIMIZATION !!!!!!
         call pft_clm_to_pso(ityp(p, nv), pso_pft)
         g1_pso_val = pso_vals%param_vals(pso_pft, this_particle)
         g0_pso_pft = pso_pft + 5
         g0_pso_val = pso_vals%param_vals(g0_pso_pft, this_particle)
         !offline_const = 0.025
         !offline_int = -0.163747
         !mbbopt(p) = pso_val*(offline_int + map_val_tile_real*offline_const)
         ! for experiment to just set g1 equal to ai
         mbbopt(p) = g1_pso_val
         ! for experiment where both g0 and g1 are optimized
         !mbbopt(p) = g0_pso_val + g1_pso_val*map_val_tile_real
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!! UNCOMMENT BELOW FOR C1 C2 BASED PSO OPTIMIZATION !!!!!
         !pso_intercept = pso_vals%param_vals(1,this_particle)
         !pso_slope = pso_vals%param_vals(2,this_particle)
         !mbbopt(p) = pso_intercept + pso_slope*map_val_tile_real
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (mbbopt(p) < 0.5) then
            mbbopt(p) = 0.5
         endif
         if (p == 1 .and. nv == 1) then
            !write(*,*) 'pso_pft'
            !write(*,*) pso_pft
            !write(*,*) 'this_particle'
            !write(*,*) this_particle
            !write(*,*) 'pso_vals%param_vals'
            !write(*,*) pso_vals%param_vals
            !write(*,*) 'pso_val'
            !write(*,*) pso_val
            !write(*,*) 'mbbopt(p)'
            !write(*,*) mbbopt(p)
            continue
         endif
      endif
      ! Soil water stress applied to Ball-Berry parameters

      bbb(p) = max (bbbopt(p)*btran(p), 1._r8)
      mbb(p) = mbbopt(p)
      ! kc, ko, cp, from: Bernacchi et al (2001) Plant, Cell and Environment 24:253-259
      !
      !       kc25 = 404.9 umol/mol
      !       ko25 = 278.4 mmol/mol
      !       cp25 = 42.75 umol/mol
      !
      ! Derive sco from cp and O2 using present-day O2 (0.209 mol/mol) and re-calculate
      ! cp to account for variation in O2 using cp = 0.5 O2/sco
      !

      kc25 = (404.9_r8/1.e06_r8) * forc_pbot(g)
      ko25 = (278.4_r8/1.e03_r8) * forc_pbot(g)
      sco  = 0.5_r8*0.209_r8 / (42.75_r8/1.e06_r8)
      cp25 = 0.5_r8*oair(p) / sco

      kc(p) = kc25*ft(t_veg(p), kcha)
      ko(p) = ko25*ft(t_veg(p), koha)
      cp(p) = cp25*ft(t_veg(p), cpha)
      
      endif

   end do

   ! Multi-layer parameters scaled by leaf nitrogen profile.
   ! Loop through each canopy layer to calculate nitrogen profile using
   ! cumulative lai at the midpoint of the layer

   do p = 1, nch        ! gkw: loop over tiles
   
      if(ityp(p, nv) > 0) then
      
      g = p
      ivt = ityp(p, nv) !  mapped vegetation type into CLM PFT

      ! Leaf nitrogen concentration at the top of the canopy (g N leaf/m**2 leaf)

      lnc(p) = 1._r8 / (slatop(ivt) * leafcn(ivt))

      ! vcmax25 at canopy top, as in CN but using lnc at top of the canopy

      vcmax25top = lnc(p) * flnr(ivt) * fnr*act25*dayl_factor(p)
!#ifndef CN
!      vcmax25top = vcmax25top*fnitr(ivt)
!#else
      if ( CNAllocation_Carbon_only() ) vcmax25top = vcmax25top*fnitr(ivt) 
!#endif

      ! Parameters derived from vcmax25top. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
      ! used jmax25 = 1.97 vcmax25, from Wullschleger (1993) Journal of Experimental Botany 44:907-920.

!KO      jmax25top = 1.97_r8*vcmax25top
!KO
      jmax25top = (2.59_r8-0.035_r8*min(max((t10(p)-tfrz), 11._r8), 35._r8)) * vcmax25top
!KO
      tpu25top = 0.167_r8*vcmax25top
      kp25top = 20000._r8*vcmax25top

      ! Nitrogen scaling factor. Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593 used
      ! kn = 0.11. Here, derive kn from vcmax25 as in Lloyd et al (2010) Biogeosciences, 7, 1833-1859
      ! Remove daylength factor from vcmax25 so that kn is based on maximum vcmax25
      ! But not used as defined here if using sun/shade big leaf code. Instead, 
      ! will use canopy integrated scaling factors from SurfaceAlbedo.

      if (dayl_factor(p) .eq. 0._r8) then
         kn(p) =  0._r8
      else
         kn(p) = exp(0.00963_r8*vcmax25top/dayl_factor(p) - 2.43_r8)
      end if

#if (defined CN)
      ! Leaf maintenance respiration to match the base rate used in CN
      ! but with the new temperature functions for C3 and C4 plants.
      !
      ! Base rate for maintenance respiration is from:
      ! M. Ryan, 1991. Effects of climate change on plant respiration.
      ! Ecological Applications, 1(2), 157-167.
      ! Original expression is br = 0.0106 molC/(molN h)
      ! Conversion by molecular weights of C and N gives 2.525e-6 gC/(gN s)
      !
      ! Base rate is at 20C. Adjust to 25C using the CN Q10 = 1.5
      !
      ! CN respiration has units:  g C/g N [leaf] / s. This needs to be
      ! converted from g C/g N [leaf] / s to umol CO2/m**2 [leaf] / s
      !
      ! Then scale this value at the top of the canopy for canopy depth

      lmr25top = 2.525e-6_r8 * (1.5_r8 ** ((25._r8-20._r8)/10._r8))
      lmr25top = lmr25top*lnc(p) / 12.e-06_r8

#else
      ! Leaf maintenance respiration in proportion to vcmax25top

      if (c3flag(p)) then
         lmr25top = vcmax25top*0.015_r8
      else
         lmr25top = vcmax25top*0.025_r8
      end if

#endif

      ! Loop through canopy layers (above snow). Respiration needs to be
      ! calculated every timestep. Others are calculated only if daytime

      laican = 0._r8
      do iv = 1, nrad(p)  ! gkw: only one canopy layer, this will be 0 if vegetation covered by snow; use 1 for now

         ! Cumulative lai at middle of layer

         if (iv == 1) then
            laican = 0.5_r8*tlai_z(p, iv)
         else
            laican = laican+0.5_r8 * (tlai_z(p, iv-1)+tlai_z(p, iv))
         end if

         ! Scale for leaf nitrogen profile. If multi-layer code, use explicit
         ! profile. If sun/shade big leaf code, use canopy integrated factor.

         if (nlevcan == 1) then
            nscaler = vcmaxcint(p, nv)
         else if (nlevcan > 1) then
            nscaler = exp(-kn(p) * laican)
         end if

         ! Maintenance respiration

         lmr25 = lmr25top*nscaler
         if (c3flag(p)) then
            lmr_z(p, iv) = lmr25*ft(t_veg(p), lmrha) * fth(t_veg(p), lmrhd, lmrse, lmrc)
         else
            lmr_z(p, iv) = lmr25*2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
            lmr_z(p, iv) = lmr_z(p, iv) / (1._r8+exp( 1.3_r8*(t_veg(p)-(tfrz+55._r8)) ))
         end if

         if (par_z(p, iv) <= 0._r8) then           ! night time

            vcmax_z(p, iv) = 0._r8
            jmax_z(p, iv) = 0._r8
            tpu_z(p, iv) = 0._r8
            kp_z(p, iv) = 0._r8

!! fzeng: comment out C13 for now. Add it back later!!
!            if ( use_c13 ) then
!               alphapsn(p) = 1._r8
!            end if

         else                                     ! day time

            vcmax25 = vcmax25top*nscaler
            jmax25 = jmax25top*nscaler
            tpu25 = tpu25top*nscaler
            kp25 = kp25top*nscaler

            ! Adjust for temperature

!KO
            vcmaxse = 668.39_r8-1.07_r8*min(max((t10(p)-tfrz), 11._r8), 35._r8)
            jmaxse  = 659.70_r8-0.75_r8*min(max((t10(p)-tfrz), 11._r8), 35._r8)
            tpuse = vcmaxse
            vcmaxc = fth25 (vcmaxhd, vcmaxse)
            jmaxc  = fth25 (jmaxhd, jmaxse)
            tpuc   = fth25 (tpuhd, tpuse)
!KO
            vcmax_z(p, iv) = vcmax25*ft(t_veg(p), vcmaxha) * fth(t_veg(p), vcmaxhd, vcmaxse, vcmaxc)
            jmax_z(p, iv) = jmax25*ft(t_veg(p), jmaxha) * fth(t_veg(p), jmaxhd, jmaxse, jmaxc)
            tpu_z(p, iv) = tpu25*ft(t_veg(p), tpuha) * fth(t_veg(p), tpuhd, tpuse, tpuc)

            if (.not. c3flag(p)) then
               vcmax_z(p, iv) = vcmax25*2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)
               vcmax_z(p, iv) = vcmax_z(p, iv) / (1._r8+exp( 0.2_r8*((tfrz+15._r8)-t_veg(p)) ))
               vcmax_z(p, iv) = vcmax_z(p, iv) / (1._r8+exp( 0.3_r8*(t_veg(p)-(tfrz+40._r8)) ))
            end if

            kp_z(p, iv) = kp25*2._r8**((t_veg(p)-(tfrz+25._r8))/10._r8)

         end if

         ! Adjust for soil water

         vcmax_z(p, iv) = vcmax_z(p, iv) * btran(p)
         lmr_z(p, iv) = lmr_z(p, iv) * btran(p)

      end do       ! canopy layer loop
      endif
   end do          ! fzeng: tile loop

   !==============================================================================!
   ! Leaf-level photosynthesis and stomatal conductance
   !==============================================================================!

   rsmax0 = 2.e4_r8

   do p = 1, nch        ! gkw: loop over tiles
   
      if(ityp(p, nv) > 0) then

      g = p
      ivt = ityp(p, nv) !  mapped vegetation type into CLM PFT

      ! Leaf boundary layer conductance, umol/m**2/s

      cf = forc_pbot(g)/(rgas*1.e-3_r8*tgcm(p))*1.e06_r8
      gb = 1._r8/rb(p)
      gb_mol(p) = gb*cf

      ! Loop through canopy layers (above snow). Only do calculations if daytime

      do iv = 1, nrad(p)

         if (par_z(p, iv) <= 0._r8) then           ! night time

            ac(p, iv) = 0._r8
            aj(p, iv) = 0._r8
            ap(p, iv) = 0._r8
            ag(p, iv) = 0._r8
            an(p, iv) = ag(p, iv) - lmr_z(p, iv)
            psn_z(p, iv) = 0._r8
            psn_wc_z(p, iv) = 0._r8
            psn_wj_z(p, iv) = 0._r8
            psn_wp_z(p, iv) = 0._r8
            rs_z(p, iv) = min(rsmax0, 1._r8/bbb(p) * cf)
            ci_z(p, iv) = 0._r8
!KO
            rh_leaf(p) = 0._r8
!KO

         else                                     ! day time

            !now the constraint is no longer needed, Jinyun Tang
            ceair = min( eair(p),  esat_tv(p) )
            rh_can = ceair/esat_tv(p)

            ! Electron transport rate for C3 plants. Convert par from W/m2 to 
            ! umol photons/m**2/s using the factor 4.6

            qabs  = 0.5_r8 * (1._r8-fnps) * par_z(p, iv) * 4.6_r8
            aquad = theta_psii
            bquad = -(qabs+jmax_z(p, iv))
            cquad = qabs*jmax_z(p, iv)
            call quadratic (aquad, bquad, cquad, r1, r2)
            je(p) = min(r1, r2)

            ! Iterative loop for ci beginning with initial guess

            if (c3flag(p)) then
               ci_z(p, iv) = 0.7_r8*cair(p)
            else
               ci_z(p, iv) = 0.4_r8*cair(p)
            end if

            niter = 0

            ! Increment iteration counter. Stop if too many iterations

            niter = niter+1

            ! Save old ci

            ciold = ci_z(p, iv)
	    
            !find ci and stomatal conductance

            ! gkw & fzeng: modified the arguments
            ! twr: added btran, t_veg as argument	    
            call hybrid(ciold, gb_mol(p), je(p), cair(p), oair(p),                                 &
               lmr_z(p, iv), par_z(p, iv), rh_can, gs_mol(p, iv), niter,                              &
	       c3flag(p), vcmax_z(p, iv), cp(p), kc(p), ko(p), qe(p), tpu_z(p, iv), kp_z(p, iv),      &
	       theta_cj(p), forc_pbot(g), bbb(p), mbb(p), ac(p, iv), aj(p, iv), ap(p, iv), ag(p, iv), an(p, iv), btran(p), t_veg(p), stomatal_model_choice)

            ! End of ci iteration.  Check for an < 0, in which case gs_mol = bbb
            
            ! Modified by Jinyun Tang, Jan 2017
            ! if (an(p, iv) < 0._r8) gs_mol(p, iv) = bbb(p)
            ! Brutely force an to zero if gs_mol is at its minimal value
            if(abs(gs_mol(p, iv)-bbb(p))<1.e-14_r8) then
              an(p, iv)=1.e-20  ! 0._r8
              ag(p, iv)=max(lmr_z(p, iv), 1.e-20)  ! lmr_z(p, iv)
            end if

            ! Final estimates for cs and ci (needed for early exit of ci iteration when an < 0)

            cs = cair(p) - 1.4_r8/gb_mol(p) * an(p, iv) * forc_pbot(g)
            cs = max(cs, 1.e-06_r8)

!            ci_z(p, iv) = cair(p) - an(p, iv) * forc_pbot(g) * (1.4_r8*gs_mol(p, iv)+1.6_r8*gb_mol(p)) / (gb_mol(p)*gs_mol(p, iv))

            ! suggested by Jinyun Tang, Jan 2017: only update ci_z when gs_mol(p, iv) is very close to bbb
            if(abs(gs_mol(p, iv)-bbb(p))<1.e-14_r8) then
              ci_z(p, iv) = cair(p) - an(p, iv) * forc_pbot(g) * (1.4_r8*gs_mol(p, iv)+1.6_r8*gb_mol(p)) / (gb_mol(p)*gs_mol(p, iv))
            else
              ci_z(p, iv) = ciold 
            end if 

            ! Make sure ci is correct. fzeng, 24 Jan 2017
            if(an(p, iv)<0) then
              print *, 'negative an:',an(p, iv)
              stop 'Photosynthesis: negative an'
            end if
            if(ci_z(p, iv)<0. .or. ci_z(p, iv)>cair(p)) then
              print *, 'ci out of bound:',ci_z(p, iv), 'cair:',cair(p)
              print *, 'an:',an(p, iv), 'forc_pbot:',forc_pbot(g), 'gs_mol:',gs_mol(p, iv), 'gb_mol:',gb_mol(p), 'bbb:',bbb(p)
              stop 'Photosynthesis: ci out of bound'
            end if

            ! Convert gs_mol (umol H2O/m**2/s) to gs (m/s) and then to rs (s/m)

            gs = gs_mol(p, iv) / cf
            rs_z(p, iv) = min(1._r8/gs, rsmax0)

            ! Photosynthesis. Save rate-limiting photosynthesis

            psn_z(p, iv) = ag(p, iv)

            psn_wc_z(p, iv) = 0._r8
            psn_wj_z(p, iv) = 0._r8
            psn_wp_z(p, iv) = 0._r8
            if (ac(p, iv) <= aj(p, iv) .and. ac(p, iv) <= ap(p, iv)) then
               psn_wc_z(p, iv) =  psn_z(p, iv)
            else if (aj(p, iv) < ac(p, iv) .and. aj(p, iv) <= ap(p, iv)) then
               psn_wj_z(p, iv) =  psn_z(p, iv)
            else if (ap(p, iv) < ac(p, iv) .and. ap(p, iv) < aj(p, iv)) then
               psn_wp_z(p, iv) =  psn_z(p, iv)
            end if

            ! Make sure iterative solution is correct

            if (gs_mol(p, iv) < 0._r8) then
               write (iulog, *) 'Negative stomatal conductance:'
               write (iulog, *) gs_mol(p, iv)
               stop  ! call endrun()
            end if

            ! Compare with Ball-Berry model: gs_mol = m*an*hs/cs p+b

            hs = (gb_mol(p)*ceair+gs_mol(p, iv)*esat_tv(p)) / ((gb_mol(p)+gs_mol(p, iv))*esat_tv(p))
!KO
            rh_leaf(p) = hs
!KO
            gs_mol_err = mbb(p)*max(an(p, iv), 0._r8)*hs/cs*forc_pbot(g) + bbb(p)

!!          if (abs(gs_mol(p, iv)-gs_mol_err) > 1.e-01_r8) then        ! gkw: too stringent for 32-bit real
	    if (abs(gs_mol(p, iv)-gs_mol_err) > 1.e-6*gs_mol_err) then
               !write (iulog, *) 'Ball-Berry Value:'
               !write (iulog, *) gs_mol_err
               !write (iulog, *) 'Medlyn value:'
               !write (iulog, *) gs_mol(p, iv)
            end if

         end if    ! night or day if branch
      end do       ! canopy layer loop
     end if        ! vegetated or not if branch
   end do          ! tile loop
   
   !==============================================================================!
   ! Canopy photosynthesis and stomatal conductance
   !==============================================================================!

   ! Sum canopy layer fluxes and then derive effective leaf-level fluxes (per
   ! unit leaf area), which are used in other parts of the model. Here, laican
   ! sums to either laisun or laisha.

   do p = 1, nch
   
     if(ityp(p, nv) > 0) then
     
      psncan = 0._r8
      psncan_wc = 0._r8
      psncan_wj = 0._r8
      psncan_wp = 0._r8
      lmrcan = 0._r8
      gscan = 0._r8
      cican = 0._r8        ! gkw: for sif calculation
      laican = 0._r8
      do iv = 1, nrad(p)
         psncan = psncan+psn_z(p, iv) * lai_z(p, iv)
         psncan_wc = psncan_wc+psn_wc_z(p, iv) * lai_z(p, iv)
         psncan_wj = psncan_wj+psn_wj_z(p, iv) * lai_z(p, iv)
         psncan_wp = psncan_wp+psn_wp_z(p, iv) * lai_z(p, iv)
         lmrcan = lmrcan+lmr_z(p, iv) * lai_z(p, iv)
         gscan = gscan+lai_z(p, iv) / (rb(p)+rs_z(p, iv))
	 cican = cican+ci_z(p, iv) * lai_z(p, iv)
         laican = laican+lai_z(p, iv)
      end do
      if (laican > 0._r8) then
         psn(p, nv) = psncan/laican
         psn_wc(p) = psncan_wc/laican
         psn_wj(p) = psncan_wj/laican
         psn_wp(p) = psncan_wp/laican
         lmr(p, nv) = lmrcan/laican
	 cican = cican/laican
         rs(p, nv) = laican/gscan-rb(p)
	 
        ! fluorescence; code from Jung-Eun Lee, implemented & modified by gkw 1/28/14; adapted from stomata 8/15/15; needs work

         if(apar(p, nv) > 0.) then
           ivt = ityp(p, nv)
           
           ! fzeng: assume pftcon%qe25 values in CN_Driver just to get the code compiling. 
           ! Need to ask Jung-Eun Lee for the CLM4.5 version of fluorescence calculation!!
           ppf = 4.6*apar(p, nv)                            ! gkw: taken from stomata, may not be correct usage
           j = ppf*pftcon%qe25(ivt)                        ! gkw: taken from stomata; may not be correct usage
 
           je_sif = max(psn(p, nv)*(cican+2.*cp(p))/max(cican+2.*cp(p)-3.*c3psn(ivt)*cp(p), 1.e-8), 0.)  ! gkw: may not be correct here

           xn = 1.-je_sif/j      ! gkw: 0.8 factor removed 20141108 (email from Jung-Eun)  ! gkw: could use je(p) here...
           xn = max(xn, 0.)       ! gkw: added 1/31/14

           if (psn_wj(p) <= 0.)  xn = 0.
           call fluorescence(xn, fs)
           sif(p, nv) = fs*ppf
          else
           sif(p, nv) = 0.
         endif
	 
      else                      ! vegetated point but LAI = 0
         psn(p, nv) =  0._r8
         psn_wc(p) =  0._r8
         psn_wj(p) =  0._r8
         psn_wp(p) =  0._r8
         lmr(p, nv) =  0._r8
         rs(p, nv)  =  rsmax0
	 sif(p, nv) =  0._r8
      end if
      
     else
     
! noveg: ITYP = 0
! -------------
         psn(p, nv) =  0.
         psn_wc(p) =  0.
         psn_wj(p) =  0.
         psn_wp(p) =  0.
         lmr(p, nv) =  0.
         rs(p, nv)  =  rsmax0
         sif(p, nv) =  0.
     
     endif
   end do
   
   end do  ! end FVEG loop

   end subroutine Photosynthesis

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ci_func
!
! !INTERFACE:
   subroutine ci_func(ci, fval, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
      c3flag, vcmax_z, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, forc_pbot, bbb, mbb, ac, aj, ap, ag, an, btran, t_veg, stomatal_model_choice)
   
   !
   !! DESCRIPTION:
   ! evaluate the function
   ! f(ci)=ci - (ca - (1.37rb+1.65rs))*patm*an
   
   ! remark:  I am attempting to maintain the original code structure, also
   ! considering one may be interested to output relevant variables for the
   ! photosynthesis model, I have decided to add these relevant variables to
   ! the clmtype structure.
   
   ! !REVISION HISTORY:
   ! Dec 14, 2012: Created by Jinyun Tang
   ! March 2022: twr changed to include Medlyn option
   !!USES
   use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
   use clmtype
   use pso_params
   
   !
   !!ARGUMENTS:
   implicit none
      
   real(r8), intent(in):: ci                ! intracellular leaf CO2 (Pa)
   real(r8), intent(in):: lmr_z             ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
   real(r8), intent(in):: par_z             ! par absorbed per unit lai for canopy layer (w/m**2)
   real(r8), intent(in):: gb_mol            ! leaf boundary layer conductance (umol H2O/m**2/s)
   real(r8), intent(in):: je                ! electron transport rate (umol electrons/m**2/s)
   real(r8), intent(in):: cair              ! Atmospheric CO2 partial pressure (Pa)
   real(r8), intent(in):: oair              ! Atmospheric O2 partial pressure (Pa)
   real(r8), intent(in):: rh_can            ! canopy air realtive humidity
   real(r8), intent(out):: fval             ! return function of the value f(ci)
   real(r8), intent(out):: gs_mol           ! leaf stomatal conductance (umol H2O/m**2/s)

! gkw & fzeng: make the following argument list variables; modify this subroutine, hybrid, and brent
   logical,  intent(in):: c3flag            ! true if C3 and false if C4
   integer,  intent(in):: stomatal_model_choice  ! Ball-Berry model used if 0, Medlyn model used if 1
   real(r8), intent(in):: vcmax_z           ! maximum rate of carboxylation (umol co2/m**2/s)
   real(r8), intent(in):: cp		     ! CO2 compensation point (Pa)
   real(r8), intent(in):: kc		     ! Michaelis-Menten constant for CO2 (Pa)
   real(r8), intent(in):: ko		     ! Michaelis-Menten constant for O2 (Pa)
   real(r8), intent(in):: qe		     ! quantum efficiency, used only for C4 (mol CO2/mol photons)
   real(r8), intent(in):: tpu_z             ! triose phosphate utilization rate (umol CO2/m**2/s)
   real(r8), intent(in):: kp_z              ! initial slope of CO2 response curve (C4 plants)
   real(r8), intent(in):: theta_cj          ! empirical curvature parameter for ac, aj photosynthesis co-limitation
   real(r8), intent(in):: forc_pbot         ! atmospheric pressure (Pa)
   real(r8), intent(in):: bbb  	     ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
   real(r8), intent(in):: mbb  	     ! Ball-Berry slope of conductance-photosynthesis relationship   
   real(r8), intent(out):: ac  	     ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
   real(r8), intent(out):: aj  	     ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
   real(r8), intent(out):: ap  	     ! product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
   real(r8), intent(out):: ag  	     ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
   real(r8), intent(out):: an  	     ! net leaf photosynthesis (umol CO2/m**2/s)
   real(r8), intent(in):: btran         !water stress factor
   real(r8), intent(in):: t_veg         !Canopy temperature (K)
!!CALLED FROM:
! subroutine hybrid and brent in this module
   
   !local variables   
   real(r8):: ai                   ! intermediate co-limited photosynthesis (umol CO2/m**2/s)
   real(r8):: cs                   ! CO2 partial pressure at leaf surface (Pa)
   
   real(r8):: aquad, bquad, cquad  ! terms for quadratic equations
   real(r8):: r1, r2               ! roots of quadratic equation
   real(r8):: fnps                 ! fraction of light absorbed by non-photosynthetic pigments
   real(r8):: theta_psii           ! empirical curvature parameter for electron transport rate
   real(r8):: theta_ip             ! empirical curvature parameter for ap photosynthesis co-limitation
   real(r8):: d                    ! used in VPD calculation according to FATES numberical implementation of Medlyn
   real(r8):: Da                   ! used in VPD calculation according to FATES numberical implementation of Medlyn 
   real(r8):: svp                  ! saturation vapor pressure
   real(r8):: t_veg_C              ! vegitation temperature in Celcius
   integer :: this_particle

   ! Miscellaneous parameters, from Bonan et al (2011) JGR, 116, doi:10.1029/2010JG001593
   fnps = 0.15_r8
   theta_psii = 0.7_r8
   theta_ip = 0.95_r8

   if (c3flag) then
   
      ! C3: Rubisco-limited photosynthesis
      ac = vcmax_z*max(ci-cp, 0._r8) / (ci+kc*(1._r8+oair/ko))

      ! C3: RuBP-limited photosynthesis
      aj = je*max(ci-cp, 0._r8) / (4._r8*ci+8._r8*cp)

      ! C3: Product-limited photosynthesis 
      ap = 3._r8*tpu_z

   else

      ! C4: Rubisco-limited photosynthesis
      ac = vcmax_z

      ! C4: RuBP-limited photosynthesis
      aj = qe*par_z*4.6_r8

      ! C4: PEP carboxylase-limited (CO2-limited)
      ap = kp_z*max(ci, 0._r8) / forc_pbot

   end if

   ! Gross photosynthesis. First co-limit ac and aj. Then co-limit ap

   aquad = theta_cj
   bquad = -(ac+aj)
   cquad = ac*aj
   call quadratic (aquad, bquad, cquad, r1, r2)
   ai = min(r1, r2)

   aquad = theta_ip
   bquad = -(ai+ap)
   cquad = ai*ap
   call quadratic (aquad, bquad, cquad, r1, r2)
   ag = min(r1, r2)

   ! Net photosynthesis. Exit iteration if an < 0

   an = ag-lmr_z
!   if (an < 0._r8) then
!      fval = 0._r8
!      return
!   endif
   
   ! Quadratic gs_mol calculation with an known. Valid for an >= 0.
   ! With an <= 0, then gs_mol = bbb
   !write(*,*) 'mbb when calculating'
   !write(*,*) mbb
   !write(*,*) 'bbb when calculating' 
   !write(*,*) bbb

   if(an <= 0.0)then
     gs_mol = bbb
     if(aj <= 1.e-20_r8)then
       fval = 0._r8
     else
       fval = ci-cair
     endif
   else if (stomatal_model_choice == 0)then  ! WHERE twr BEGINS EDITING!!!!
     cs = cair-1.4_r8/gb_mol*an*forc_pbot  ! divide by gb_mol because 1/gb_mol = rb, an is net photynthesis, forc_pbot is atoms. pressure (pa)
     cs = max(cs, 1.e-06_r8)  ! make sure cs is greater than 1e-06
     aquad = cs
     bquad = cs*(gb_mol-bbb) - mbb*an*forc_pbot
     cquad = -gb_mol*(cs*bbb+mbb*an*forc_pbot*rh_can)
     call quadratic (aquad, bquad, cquad, r1, r2)
     gs_mol = max(r1, r2)
     ! Derive new estimate for ci
     fval = ci-cair+an*forc_pbot * (1.4_r8*gs_mol+1.6_r8*gb_mol) / (gb_mol*gs_mol)
   else if (stomatal_model_choice == 1)then
    !calculate g1 from environmental filtering (here refered to as mbb, in line with original CLM45 terminology)
    !for EF, use mean annual temperature, mean annual precipitation, and Precip/PET (How to import these into the model????)
    !mbb = 4.1  ! TWR TESTING G1 ESTIMATE FROM Li 2022
    cs = cair-1.4_r8/gb_mol*an*forc_pbot  ! divide by gb_mol because 1/gb_mol = rb, an is net photynthesis, forc_pbot is atoms. pressure (pa)
    cs = max(cs, 1.e-06_r8)  ! make sure cs is greater than 1e-06
    !aquad = cs
    aquad = 1.0_r8
    !bquad = cs*(gb_mol-bbb) - mbb*an*forc_pbot
    d = (1.6_r8*an)/(cs/forc_pbot)
    t_veg_C = t_veg-273.15  ! convert canopy temp from K to C
    svp = 610.78*exp(t_veg_C/((t_veg_C+238.3)*17.2694))  ! calculate saturated vapor pressure
    Da = svp * (1 - (rh_can/100))  ! this is VPD
    Da = min(Da, 5000._r8)  ! constraint on VPD
    bquad = -(2._r8 * (bbb*btran+d) + (((mbb**2)*(d)**2)/(gb_mol*Da)))
    !cquad = -gb_mol*(cs*bbb+mbb*an*forc_pbot*rh_can)
    cquad = (bbb*btran)**2 + (2._r8*bbb*btran+d*(1 - (mbb**2/Da))) * d
    call quadratic (aquad, bquad, cquad, r1, r2)
    gs_mol = max(r1, r2)
    
    this_particle = pso_vals%particle_num
    !write(*,*) 'inside medlyn'
    !write(*,*) 'this_particle'
    !write(*,*) this_particle
    !write(*,*) 'mbb'
    !write(*,*) mbb
    !write(*,*) 'aquad'
    !write(*,*) aquad
    !write(*,*) ' '
    !write(*,*) 'rh_can'
    !write(*,*) rh_can
    !write(*,*) 't_veg'
    !write(*,*) t_veg
    !write(*,*) 'an'
    !write(*,*) an
    !write(*,*) 'cs'
    !write(*,*) cs
    !write(*,*) 'forc_pbot'
    !write(*,*) forc_pbot
    !write(*,*) ' '
    !write(*,*) 'bbb'
    !write(*,*) bbb
    !write(*,*) 'btran'
    !write(*,*) btran
    !write(*,*) 'd'
    !write(*,*) d
    !write(*,*) 'gb_mol'
    !write(*,*) gb_mol
    !write(*,*) 'Da'
    !write(*,*) Da
    !write(*,*) ' '
    !write(*,*) 'bquad'
    !write(*,*) bquad
    !write(*,*) ' '
    !write(*,*) 'cquad'
    !write(*,*) cquad
    !write(*,*) ' ' 
    !write(*,*) 'r1'
    !write(*,*) r1
    !write(*,*) 'r2'
    !write(*,*) r2
    !write(*,*) 'gs_mol'
    !write(*,*) gs_mol
    ! Derive new estimate for ci
    fval = ci-cair+an*forc_pbot * (1.4_r8*gs_mol+1.6_r8*gb_mol) / (gb_mol*gs_mol)
   endif
   
   end subroutine ci_func

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: quadratic
!
! !INTERFACE:
   subroutine quadratic (a, b, c, r1, r2)
!
!
! !DESCRIPTION:
!==============================================================================!
!----------------- Solve quadratic equation for its two roots-----------------!
!==============================================================================!
! Solution from Press et al (1986) Numerical Recipes: The Art of Scientific
! Computing (Cambridge University Press, Cambridge), pp. 145.
!
! !CALLED FROM:
! subroutine Photosynthesis in this module
!
! !REVISION HISTORY:
! 4/5/10: Adapted from/home/bonan/ecm/psn/An_gs_iterative.f90 by Keith Oleson
!
! !USES:
   use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
   implicit none
!
! !ARGUMENTS:
   real(r8), intent(in)  :: a, b, c       ! Terms for quadratic equation
   real(r8), intent(out):: r1, r2       ! Roots of quadratic equation
!
! !LOCAL VARIABLES:
   real(r8):: q                        ! Temporary term for quadratic solution
   
   integer:: iulog = 6                 ! gkw
!------------------------------------------------------------------------------

   if (a == 0._r8) then
      write (iulog, *) 'Quadratic solution error: a = ',a
      stop  ! call endrun()
   end if

   if (b >= 0._r8) then
      q = -0.5_r8 * (b+sqrt(b*b - 4._r8*a*c))
   else
      q = -0.5_r8 * (b-sqrt(b*b - 4._r8*a*c))
   end if

   r1 = q/a
   if (q /= 0._r8) then
      r2 = c/q
   else
      r2 = 1.e36_r8
   end if

   end subroutine quadratic

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: brent
!
! !INTERFACE:
   subroutine brent(x, x1, x2, f1, f2, tol, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
      c3flag, vcmax_z, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, forc_pbot, bbb, mbb, ac, aj, ap, ag, an, btran, t_veg, stomatal_model_choice)
      
   !
   !!DESCRIPTION:
   !Use Brent's method to find the root of a single variable function ci_func, which is known to exist between x1 and x2.
   !The found root will be updated until its accuracy is tol.
   
   !!REVISION HISTORY:
   !Dec 14/2012: Jinyun Tang, modified from numerical recipes in F90 by press et al. 1188-1189
   !
   !!USES:
   use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4
   
   !
   !!ARGUMENTS:
   implicit none
      
   real(r8), intent(out):: x                !indepedent variable of the single value function ci_func(x)
   real(r8), intent(in):: x1, x2, f1, f2    !minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
   real(r8), intent(in):: tol               !the error tolerance

   real(r8), intent(in):: lmr_z             ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
   real(r8), intent(in):: par_z             ! par absorbed per unit lai for canopy layer (w/m**2)
   real(r8), intent(in):: gb_mol            ! leaf boundary layer conductance (umol H2O/m**2/s)
   real(r8), intent(in):: je                ! electron transport rate (umol electrons/m**2/s)
   real(r8), intent(in):: cair              ! Atmospheric CO2 partial pressure (Pa)
   real(r8), intent(in):: oair              ! Atmospheric O2 partial pressure (Pa)
   real(r8), intent(in):: rh_can            ! inside canopy relative humidity 
   real(r8), intent(out):: gs_mol           ! leaf stomatal conductance (umol H2O/m**2/s)
   
! gkw & fzeng: make the following argument list variables; modify this subroutine, hybrid, and brent
   logical,  intent(in):: c3flag            ! true if C3 and false if C4
   real(r8), intent(in):: vcmax_z           ! maximum rate of carboxylation (umol co2/m**2/s)
   real(r8), intent(in):: cp		     ! CO2 compensation point (Pa)
   real(r8), intent(in):: kc		     ! Michaelis-Menten constant for CO2 (Pa)
   real(r8), intent(in):: ko		     ! Michaelis-Menten constant for O2 (Pa)
   real(r8), intent(in):: qe		     ! quantum efficiency, used only for C4 (mol CO2/mol photons)
   real(r8), intent(in):: tpu_z             ! triose phosphate utilization rate (umol CO2/m**2/s)
   real(r8), intent(in):: kp_z              ! initial slope of CO2 response curve (C4 plants)
   real(r8), intent(in):: theta_cj          ! empirical curvature parameter for ac, aj photosynthesis co-limitation
   real(r8), intent(in):: forc_pbot         ! atmospheric pressure (Pa)
   real(r8), intent(in):: bbb  	     ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
   real(r8), intent(in):: mbb  	     ! Ball-Berry slope of conductance-photosynthesis relationship   
   real(r8), intent(out):: ac  	     ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
   real(r8), intent(out):: aj  	     ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
   real(r8), intent(out):: ap  	     ! product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
   real(r8), intent(out):: ag  	     ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
   real(r8), intent(out):: an  	     ! net leaf photosynthesis (umol CO2/m**2/s)

!trent added these parameters for Medlyn implementation
   real(r8), intent(in):: btran
   real(r8), intent(in):: t_veg
   integer, intent(in):: stomatal_model_choice

! !CALLED FROM:
! subroutine hybrid in this module

   integer, parameter:: ITMAX = 30            !maximum number of iterations, increased from 20 to 30 by Jinyun Tang, Jan 2017
   real(r8), parameter:: EPS = 1.e-2_r8       !relative error tolerance
   
   integer:: iter
   real(r8)  :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
   
   integer:: iulog = 6                      ! gkw
   
   a = x1
   b = x2
   fa = f1
   fb = f2
   if((fa > 0._r8 .and. fb > 0._r8).or.(fa < 0._r8 .and. fb < 0._r8))then
      write(iulog, *) 'root must be bracketed for brent'
      stop 'brent' ! call endrun()
   endif 
   c = b
   fc = fb
   iter = 0
   do
      if(iter == ITMAX)exit
      iter = iter+1
      if((fb > 0._r8 .and. fc > 0._r8) .or. (fb < 0._r8 .and. fc < 0._r8))then
         c = a   !Rename a, b, c and adjust bounding interval d.
         fc = fa
         d = b-a
         e = d
      endif
      if( abs(fc) < abs(fb)) then
         a = b
         b = c
         c = a
         fa = fb
         fb = fc
         fc = fa
      endif
      tol1 = 2._r8*EPS*abs(b)+0.5_r8*tol  !Convergence check.   
      xm = 0.5_r8*(c-b)
      if(abs(xm) <= tol1 .or. fb == 0.)then
         x = b
         return
      endif
      if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
         s = fb/fa  ! Attempt inverse quadratic interpolation.
         if(a == c) then
            p = 2._r8*xm*s
            q = 1._r8-s
         else
            q = fa/fc
            r = fb/fc
            p = s*(2._r8*xm*q*(q-r)-(b-a)*(r-1._r8))
            q=(q-1._r8)*(r-1._r8)*(s-1._r8)
         endif
         if(p > 0._r8) q = -q  ! Check whether in bounds.
         p = abs(p)
         if(2._r8*p < min(3._r8*xm*q-abs(tol1*q), abs(e*q))) then
            e = d  ! Accept interpolation.
            d = p/q
         else
            d = xm  !Interpolation failed, use bisection.
            e = d
         endif
      else  ! Bounds decreasing too slowly, use bisection.
         d = xm
         e = d
      endif
      a = b  ! Move last best guess to a.
      fa = fb
      if(abs(d) > tol1) then  ! Evaluate new trial root.
         b = b+d
      else
         b = b+sign(tol1, xm)
      endif
      
      call ci_func(b, fb, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
                   c3flag, vcmax_z, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, forc_pbot, bbb, mbb, ac, aj, ap, ag, an, btran, t_veg, stomatal_model_choice)
      if(abs(fb)<1.e-5_r8) then         
         exit
      endif
   enddo
!   if(iter == ITMAX)write(iulog, *) 'brent exceeding maximum iterations', b, fb
   x = b
   return
   end subroutine brent
   
!-------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hybrid
!
! !INTERFACE:

   subroutine hybrid(x0, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, iter, &
     c3flag, vcmax_z, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, forc_pbot, bbb, mbb, ac, aj, ap, ag, an, btran, t_veg, stomatal_model_choice)
   !
   !! DESCRIPTION:
   ! use a hybrid solver to find the root of equation  
   ! f(x) = x-h(x), 
   !we want to find x, s.t. f(x) = 0.
   !the hybrid approach combines the strength of the newton secant approach (find the solution domain)
   !and the bisection approach implemented with the Brent's method to guarrantee convergence.
   
   !! REVISION HISTORY:
   !Dec 14/2012: created by Jinyun Tang
   !Jan 2017: modified by Jinyun Tang
   !March 2022: modified by twr to contain btran
   
   !
   !!USES:   
   use MAPL_ConstantsMod, ONLY: r8 => MAPL_R4  
   
   !
   !! ARGUMENTS:
   implicit none
      
   real(r8), intent(inout):: x0              !initial guess and final value of the solution
   real(r8), intent(in):: lmr_z              ! canopy layer: leaf maintenance respiration rate (umol CO2/m**2/s)
   real(r8), intent(in):: par_z              ! par absorbed per unit lai for canopy layer (w/m**2)
   real(r8), intent(in):: rh_can             ! canopy air relative humidity
   real(r8), intent(in):: gb_mol             ! leaf boundary layer conductance (umol H2O/m**2/s)
   real(r8), intent(in):: je                 ! electron transport rate (umol electrons/m**2/s)
   real(r8), intent(in):: cair               ! Atmospheric CO2 partial pressure (Pa)
   real(r8), intent(in):: oair               ! Atmospheric O2 partial pressure (Pa)
   real(r8), intent(out):: gs_mol            ! leaf stomatal conductance (umol H2O/m**2/s)
   integer,  intent(out):: iter              !number of iterations used, for record only 

! gkw & fzeng: make the following argument list variables; modify this subroutine, ci_func, and brent   
   logical,  intent(in):: c3flag             ! true if C3 and false if C4
   real(r8), intent(in):: vcmax_z            ! maximum rate of carboxylation (umol co2/m**2/s)
   real(r8), intent(in):: cp		      ! CO2 compensation point (Pa)
   real(r8), intent(in):: kc		      ! Michaelis-Menten constant for CO2 (Pa)
   real(r8), intent(in):: ko		      ! Michaelis-Menten constant for O2 (Pa)
   real(r8), intent(in):: qe		      ! quantum efficiency, used only for C4 (mol CO2/mol photons)
   real(r8), intent(in):: tpu_z              ! triose phosphate utilization rate (umol CO2/m**2/s)
   real(r8), intent(in):: kp_z               ! initial slope of CO2 response curve (C4 plants)
   real(r8), intent(in):: theta_cj           ! empirical curvature parameter for ac, aj photosynthesis co-limitation
   real(r8), intent(in):: forc_pbot          ! atmospheric pressure (Pa)
   real(r8), intent(in):: bbb  	      ! Ball-Berry minimum leaf conductance (umol H2O/m**2/s)
   real(r8), intent(in):: mbb  	      ! Ball-Berry slope of conductance-photosynthesis relationship
   real(r8), intent(in):: btran          !soil water stress factor
   real(r8), intent(in):: t_veg          ! vegitation temperature (K)
   integer, intent(in)  :: stomatal_model_choice  ! 0 for Ball-Berry, 1 for Medlyn
   real(r8), intent(out):: ac                ! Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
   real(r8), intent(out):: aj                ! RuBP-limited gross photosynthesis (umol CO2/m**2/s)
   real(r8), intent(out):: ap                ! product-limited (C3) or CO2-limited (C4) gross photosynthesis (umol CO2/m**2/s)
   real(r8), intent(out):: ag                ! co-limited gross leaf photosynthesis (umol CO2/m**2/s)
   real(r8), intent(out):: an                ! net leaf photosynthesis (umol CO2/m**2/s)  

! !CALLED FROM:
! subroutine photosynthesis in this module

   !local variables
   real(r8):: x1, f0, f1
   real(r8):: x, dx
   real(r8), parameter:: eps = 1.e-2_r8      !relative accuracy
   real(r8), parameter:: eps1 = 1.e-4_r8
   integer,  parameter:: itmax = 40          !maximum number of iterations
   real(r8):: tol, minx, minf

   real(r8):: ci_val(5)
   real(r8):: fi_val(5)
   integer  :: ii, mi
   
   integer:: iulog = 6
   
   iter = 0
   call ci_func(x0, f0, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
                c3flag, vcmax_z, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, forc_pbot, bbb, mbb, ac, aj, ap, ag, an, btran, t_veg, stomatal_model_choice)
   if(abs(f0) < 1.e-14_r8)return
   ci_val(3)=x0
   fi_val(3)=f0
   
   !compute the minimum ci value
   if(c3flag)then
     ci_val(1)=cp+1.e-6_r8
   else
     ci_val(1)=1.e-6_r8
   endif
   call ci_func(ci_val(1), fi_val(1), gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
                c3flag, vcmax_z, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, forc_pbot, bbb, mbb, ac, aj, ap, ag, an, btran, t_veg, stomatal_model_choice)
   
   ci_val(2)=(ci_val(1)+ci_val(3))*0.5
   call ci_func(ci_val(2), fi_val(2), gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
                c3flag, vcmax_z, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, forc_pbot, bbb, mbb, ac, aj, ap, ag, an, btran, t_veg, stomatal_model_choice)
   
   ci_val(4)=(cair+ci_val(3))*0.5
   call ci_func(ci_val(4), fi_val(4), gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
                c3flag, vcmax_z, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, forc_pbot, bbb, mbb, ac, aj, ap, ag, an, btran, t_veg, stomatal_model_choice)
   
   !compute the maximum ci value
   ci_val(5)=cair*0.999_r8
   call ci_func(ci_val(5), fi_val(5), gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
                c3flag, vcmax_z, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, forc_pbot, bbb, mbb, ac, aj, ap, ag, an, btran, t_veg, stomatal_model_choice)
   
   mi = -1
   do ii = 1, 4
     if(fi_val(ii)*fi_val(ii+1)<0._r8)then
        mi = ii
     endif
   end do
   if(mi > 0)then

      x0 = ci_val(mi)
      f0 = fi_val(mi)
      x1 = ci_val(mi+1)
      f1 = fi_val(mi+1)
      tol = 0.5 * (x0+x1) * eps    ! This is missing in Jinyun's modifications. fzeng added following CanopyFluxesMod.F90, 9 Nov 2017
      call brent(x, x0, x1, f0, f1, tol, gb_mol, je, cair, oair, lmr_z, par_z, rh_can, gs_mol, &
                    c3flag, vcmax_z, cp, kc, ko, qe, tpu_z, kp_z, theta_cj, forc_pbot, bbb, mbb, &
                    ac, aj, ap, ag, an, btran, t_veg, stomatal_model_choice)
      
      ! fzeng added for debugging, 24 Jan 2017
!      if(an <= 0) then
!        write(iulog, '(L8, 16(X, E15.8), X)')c3flag, ci_val(1), fi_val(1), ci_val(2), fi_val(2), &
!        ci_val(3), fi_val(3), ci_val(4), fi_val(4), ci_val(5), fi_val(5), x, &
!        ag, an, lmr_z, gs_mol, bbb
!      endif 

      x0 = x
      
   else
   
     ! write(iulog, '(I8, 13(X, E15.8), X, L8)')p, ci_val(1), fi_val(1), ci_val(2), fi_val(2), &
     !   ci_val(3), fi_val(3), ci_val(4), fi_val(4), ci_val(5), fi_val(5), photosyns_vars%aj_patch(p, iv), &
     !   photosyns_vars%ag_patch(p, iv), photosyns_vars%an_patch(p, iv), all(fi_val < 0._r8)
     ! write(iulog, *)'no solution for ci and gs_mol is forced to be the minimum for pft',p
     ! write(iulog, *)'no solution found for ci'
     ! call endrun(decomp_index = p, clmlevel = namep, msg = errmsg(__FILE__, __LINE__))

     ! Jinyun commented this print statements out. 
     ! It seems that when the if condition above is not met, all(fi_val < 0) is true. 
     ! So the aj, ag, an etc. from call ci_func(ci_val(5), fi_val(5) ...) are the final output of hybrid. 
     ! fzeng, 9 Nov 2017

!      write(iulog, *)'no solution found for ci'     
!      write(iulog, '(L8, 16(X, E15.8), X, L8)')c3flag, ci_val(1), fi_val(1), ci_val(2), fi_val(2), &
!        ci_val(3), fi_val(3), ci_val(4), fi_val(4), ci_val(5), fi_val(5), aj, &
!        ag, an, lmr_z, gs_mol, bbb, all(fi_val < 0._r8)

      x0 = ci_val(5)  ! fzeng, 9 Nov 2017
     
   endif
      
   end subroutine hybrid
   
!------------------------------------------------------------------------------
!
! !IROUTINE: Fluorescence
!
! !INTERFACE:
   subroutine fluorescence(x, fs)
!
! !DESCRIPTION: 
! Chlorophyll fluorescence
! writen by Jung-Eun Lee using van der Tol and Berry (2012)

! !USES:
     implicit none
     real, intent(in)    :: x       ! degree of light saturation
     real, intent(out)   :: fs      ! fluorescence yield
     real:: Kn      ! rate constant for non-photochemical quenching
     real:: Kf      ! rate constant for fluorescence
     real:: Kd      ! rate constant for thermal deactivation at Fm
     real:: Kp      ! rate constant for photochemisty
     real:: po0
     real:: ps
     real:: fo0
     real:: fo      ! fluorescnce yield at Fo
     real:: fm      ! fluorescnce yield at Fm
     real:: fm0
     real:: eta
     real:: qQ 
     real:: qE 

     Kf          = 0.05                 ! rate constant for fluorescence
     Kd          = 0.95                 ! rate constant for thermal deactivation at Fm
     Kp          = 4.0                  ! rate constant for photochemisty

     po0         = Kp/(Kf+Kd+Kp)        ! dark photochemistry fraction (Genty et al., 1989)
     ps          = po0*(1.-x)           ! photochemical yield
     Kn          = (6.2473*x - 0.5944)*x  ! empirical fit to Flexas' data
!    Kn          = (3.9867*x - 1.0589)*x  ! empirical fit to Flexas, Daumard, Rascher, Berry data

     fo0         = Kf/(Kf+Kp+Kd)        ! dark adapted fluorescence yield Fo
     fo          = Kf/(Kf+Kp+Kd+Kn)     ! dark adapted fluorescence yield Fo
     fm          = Kf/(Kf   +Kd+Kn)     ! light adapted fluorescence yield Fm
     fm0         = Kf/(Kf   +Kd)        ! light adapted fluorescence yield Fm
     fs          = fm*(1.-ps)           ! fluorescence as fraction of PAR
     eta         = fs/fo0               ! fluorescence as fraction of dark adapted

     qQ          = 1.-(fs-fo)/(fm-fo)   ! photochemical quenching
     qE          = 1.-(fm-fo)/(fm0-fo0)  ! non-photochemical quenching

  end subroutine fluorescence

  end module compute_rc_mod
