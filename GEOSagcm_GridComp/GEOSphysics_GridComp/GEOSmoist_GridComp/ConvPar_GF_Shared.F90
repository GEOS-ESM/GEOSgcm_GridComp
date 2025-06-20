MODULE ConvPar_GF_SharedParams

USE MAPL
USE GEOSmoist_Process_Library, only : CNV_Tracers 

 IMPLICIT NONE

 !- plume spectral size
 INTEGER ,PARAMETER  :: maxiens=3 , deep=1 , shal=2 , mid=3
 CHARACTER(LEN=10),PARAMETER,DIMENSION(maxiens)  :: cumulus_type = (/ &
                                                   'deep      ' &
                                                  ,'shallow   ' &
                                                  ,'mid       ' &
                                                                    /)
 !------------------- namelist variables
 !-- plume to be activated (1 true, 0 false): deep, shallow, congestus
 INTEGER, DIMENSION(maxiens) :: icumulus_gf = (/1,0,1/)

 !-- choice for the closures:
 !--  deep   : 0 ensemble (all)          , 1 GR, 4 ll omega, 7 moist conv, 10 PB
 !--  shallow: 0 ensemble (Wstar/BLQE)   , 1 Wstar, 4 heat-engine or 7 BLQE
 !--  mid    : 0 ensemble (Wstar/BLQE/PB), 1 Wstar, 2 BLQE, 3 PB, 4 PB_BL
 INTEGER, DIMENSION(maxiens) :: closure_choice = (/0,  7,  3/) ! deep, shallow, congestus

 !-- gross entraiment rate: deep, shallow, congestus
 REAL,       DIMENSION(maxiens) :: cum_entr_rate = (/&
                                          1.00e-4  & !deep
                                         ,2.00e-3  & !shallow
                                         ,9.00e-4  & !mid
                                                   /)

 INTEGER :: USE_TRACER_TRANSP = 1 != 0/1     - default 1

 INTEGER :: USE_TRACER_SCAVEN = 1 != 0/1/2/3 - default 2

 INTEGER :: USE_FLUX_FORM     = 1 != 1/2/3   - default 1

 INTEGER :: USE_FCT           = 1 != 0/1     - default 1 (only for USE_FLUX_FORM     = 2)

 INTEGER :: USE_TRACER_EVAP   = 1 != 0/1     - default 1 (only for USE_TRACER_SCAVEN > 0)

 INTEGER :: USE_SCALE_DEP     = 1 != 0/1:  scale dependence flag, default = 1

 INTEGER :: DICYCLE           = 1 != 0/1:  diurnal cycle closure, default = 1

 REAL    :: ALP1              = 1 != 0/0.5/1: apply subsidence transport of LS/anvil cloud fraction using
                                  !=          time implicit discretization

                                  != boundary condition determination for the plumes
 INTEGER :: BC_METH           = 0 ! 0: simple arithmetic mean around k22
                                  ! 1: mass weighted mean around k22

 REAL,   DIMENSION(maxiens) :: CUM_AVE_LAYER     =(/50.,   30.,   50. /)!= layer depth for average the properties
                                                                        != of source air parcels (mbar)
 REAL    ::  AVE_LAYER         != layer depth for average the properties of source air parcels (mbar)

 REAL    ::  TAU_DEEP         = 5400.  != deep      convective timescale
 REAL    ::  TAU_MID          = 3600.  != congestus convective timescale

 REAL    ::  C0_DEEP          = 2.e-3 != default= 3.e-3   conversion rate (cloud to rain, m-1) - for deep      plume
 REAL    ::  C0_MID           = 2.e-3 != default= 2.e-3   conversion rate (cloud to rain, m-1) - for congestus plume
 REAL    ::  C0_SHAL          = 0.    != default= 0.e-3   conversion rate (cloud to rain, m-1) - for shallow   plume
 REAL    ::  QRC_CRIT         = 2.e-4 != default= 2.e-4   kg/kg
 REAL    ::  QRC_CRIT_LND     = 3.e-4 != default= 2.e-4   kg/kg
 REAL    ::  QRC_CRIT_OCN     = 3.e-4 != default= 2.e-4   kg/kg
 REAL    ::  C1               = 0.0   != default= 1.e-3   conversion rate (cloud to rain, m-1) - for the 'C1d' detrainment approach

 !- physical constants
 REAL, PARAMETER ::  &
  rgas    = 287.,    & ! J K-1 kg-1
  cp      = 1004.,   & ! J K-1 kg-1
  rv      = 461.,    & ! J K-1 kg-1
  p00     = 1.e5,    & ! hPa
  tcrit   = 258.,    & ! K
  g       = MAPL_GRAV,&! m s-2
  cpor    = cp/rgas, &
  xlv     = 2.5e6,   & ! J kg-1
  akmin   = 1.0,     & ! #
  tkmin   = 1.e-5,   & ! m+2 s-2
  ccnclean= 250.,    & ! # cm-3
  T_0     = 273.16,  & ! K
  T_ice   = 235.16,  & ! K
  xlf     = 0.333e6, & ! latent heat of freezing (J K-1 kg-1)
  max_qsat= 0.5,     & ! kg/kg
  mx_buoy = cp*5. + xlv*2.e-3 ! temp exc=5 K, q deficit=2 g/kg (=> mx_buoy ~ 10 kJ/kg)

  ! Default autoconversion parameter for GEOS-Chem species [s-1]
  REAL, PARAMETER       :: KC_DEFAULT_GCC = 5.e-3

 LOGICAL :: CNV_2MOM = .FALSE.

 CHARACTER(len=10) :: GF_ENV_SETTING = 'DYNAMICS' ! or 'CURRENT'

 CONTAINS

  subroutine get_incloud_sc_chem_up(cumulus,mtp,se,se_cup,sc_up,pw_up,tot_pw_up_chem&
                                   ,z_cup,rho,po,po_cup,qco  &
                                   ,qrco,tempco,pwo,zuo,up_massentro,up_massdetro,vvel2d,vvel1d&
                                   ,start_level,k22,kbcon,ktop,klcl,ierr,xland,itf,ktf,its,ite, kts,kte)
     IMPLICIT NONE
     !-inputs
     integer                               ,intent (in)  :: itf,ktf, its,ite, kts,kte
     integer                               ,intent (in)  :: mtp
     integer, dimension (its:ite)          ,intent (in)  :: ierr ,kbcon,ktop,k22,klcl,start_level
     character *(*)                        ,intent (in)  :: cumulus
     real, dimension (mtp ,its:ite,kts:kte),intent (in)  :: se,se_cup
     real, dimension (its:ite,kts:kte)     ,intent (in)  :: z_cup,rho,po_cup,qco,qrco,tempco,pwo,zuo &
                                                           ,up_massentro,up_massdetro,po


     real,    dimension (its:ite,kts:kte)  ,intent (in)  :: vvel2d
     real,    dimension (its:ite        )  ,intent (in)  :: vvel1d,xland

     !-outputs
     real, dimension (mtp ,its:ite,kts:kte),intent (out) :: sc_up,pw_up
     real, dimension (mtp ,its:ite        ),intent (out) :: tot_pw_up_chem

     !-locals
     real, parameter :: scav_eff = 0.6  ! for smoke : Chuang et al. (1992) J. Atmos. Sci.
     real   , dimension (mtp ,its:ite) ::  sc_b
     real   , dimension (mtp) :: conc_mxr
     real :: x_add,dz,XZZ,XZD,XZE,denom,henry_coef,w_upd,fliq,dp
     integer :: i,k,ispc
     real, parameter :: cte_w_upd = 10. ! m/s
!    real, parameter :: kc = 5.e-3  ! s-1
     real, parameter :: kc = 2.e-3  ! s-1        !!! autoconversion parameter in GF is lower than what is used in GOCART
     real, dimension (mtp ,its:ite,kts:kte) ::  factor_temp
     ! GEOS-Chem update
     real            :: this_w_upd
     real            :: fsol
     real            :: kc_scaled, ftemp, l2g
     logical         :: is_gcc

     !--initialization
     sc_up          = se_cup
     pw_up          = 0.0
     tot_pw_up_chem = 0.0

     do i=its,itf
       if(ierr(i) /= 0) cycle

       do ispc=1,mtp
             call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),se_cup(ispc,i,kts:kte),sc_b(ispc,i),k22(i))
       enddo
       do k=kts,k22(i)
            sc_up   (:,i,k) = sc_b  (:,i)
       enddo
     enddo

     do i=its,itf
          if(ierr(i) /= 0) cycle
loopk:      do k=k22(i)+1,ktop(i)+1

            !-- entr,detr, mass flux ...
            XZZ=             zuo(i,k-1)
            XZD=0.5*up_massdetro(i,k-1)
            XZE=    up_massentro(i,k-1)
            denom =  (XZZ-XZD+XZE)

            !-- transport + mixing
            if(denom > 0.) then
               sc_up(:,i,k) = (sc_up(:,i,k-1)*XZZ - sc_up(:,i,k-1)*XZD + se(:,i,k-1)*XZE) / denom
            else
               sc_up(:,i,k) = sc_up(:,i,k-1)
            endif

            !-- scavenging section
            if(USE_TRACER_SCAVEN==0 .or. trim(cumulus) == 'shallow') cycle loopk
            dz=z_cup(i,k)-z_cup(i,k-1)

            !-- in-cloud vert velocity for scavenging formulation 2
            w_upd = vvel2d(i,k)

            do ispc = 1,mtp
                !--use GEOS-Chem washout parameterization?
                is_gcc = CNV_Tracers(ispc)%use_gcc_washout
                IF(CNV_Tracers(ispc)%fscav > 1.e-6) THEN ! aerosol scavenging

                    !--formulation 1 as in GOCART with RAS conv_par
                    if(USE_TRACER_SCAVEN==1) &
                    pw_up(ispc,i,k) = max(0.,sc_up(ispc,i,k)*(1.-exp(- CNV_Tracers(ispc)%fscav  * (dz/1000.))))

                    !--formulation 2 as in GOCART
                    if(USE_TRACER_SCAVEN==2) then

                       ! Use GEOS-Chem formulation for GEOS-Chem species
                       if ( is_gcc ) then
                          if ( CNV_Tracers(ispc)%use_gocart ) then
                             kc_scaled  = kc
                             this_w_upd = w_upd
                             ftemp = 1.0
                             if ( tempco(i,k) < CNV_Tracers(ispc)%ftemp_threshold ) ftemp = 0.0
                          else
                             call compute_ki_gcc_aerosol ( ispc, tempco(i,k), kc_scaled )
                             ftemp      = CNV_Tracers(ispc)%fscav  ! apply aerosol scavenging efficiency 
                             this_w_upd = get_w_upd_gcc( vvel2d(i,k), xland(i), CNV_Tracers(ispc)%online_vud )
                          endif
                          ! if it's not a wetdep species, force kc_scaled to 0. This ensures no washout
                          if ( .not. CNV_Tracers(ispc)%is_wetdep ) kc_scaled = 0.0 
                          ! calculate soluble fraction and apply to tracer
                          fsol = min(1.,max(0.,(1.-exp(- kc_scaled * (dz/this_w_upd)))*ftemp))
                          pw_up(ispc,i,k) = sc_up(ispc,i,k)*fsol
 
                       ! Original formulation
                       else
                          pw_up(ispc,i,k) = max(0.,sc_up(ispc,i,k)*(1.-exp(- kc * (dz/w_upd)))*factor_temp(ispc,i,k))
                       endif
                    endif

                    !--formulation 3 - orignal GF conv_par
                    if(USE_TRACER_SCAVEN==3) then
                       !--- cloud liquid water tracer concentration
                       conc_mxr(ispc)  =  scav_eff* sc_up(ispc,i,k) !unit [kg(aq)/kg(air)]  for aerosol/smoke
                       !---   aqueous-phase concentration in rain water
                       pw_up(ispc,i,k) = conc_mxr(ispc)*pwo(i,k)/(1.e-8+qrco(i,k))
                    endif

                    !---(in cloud) total mixing ratio in gas and aqueous phases
                    sc_up(ispc,i,k) = sc_up(ispc,i,k) - pw_up(ispc,i,k)

                    !
                ELSEIF(CNV_Tracers(ispc)%Vect_Hcts(1)>1.e-6) THEN ! tracer gas phase scavenging

                    !--- equilibrium tracer concentration - Henry's law
                    if ( is_gcc ) then
                       henry_coef=henry_gcc(ispc,tempco(i,k),rho(i,k))
                    else
                       henry_coef=henry(ispc,tempco(i,k),rho(i,k))
                    endif

                    if(USE_TRACER_SCAVEN==3) then
                      !--- cloud liquid water tracer concentration
                      conc_mxr(ispc) = (henry_coef*qrco(i,k) /(1.+henry_coef*qrco(i,k)) )* sc_up(ispc,i,k)
                      !
                      !---   aqueous-phase concentration in rain water
                      pw_up(ispc,i,k) = conc_mxr(ispc)*pwo(i,k)/(1.e-8+qrco(i,k))

                    else

                       ! Use GEOS-Chem definitions if it's a GEOS-Chem tracer 
                       if ( is_gcc ) then
                          call compute_ki_gcc_gas ( ispc, tempco(i,k), po_cup(i,k), qco(i,k), qrco(i,k), henry_coef, &
                             kc_scaled, l2g )
                          this_w_upd = get_w_upd_gcc( vvel2d(i,k), xland(i), CNV_Tracers(ispc)%online_vud )
                          ! if it's not a wetdep species, force kc_scaled to 0. This ensures no washout
                          if ( .not. CNV_Tracers(ispc)%is_wetdep ) kc_scaled = 0.0
                          ! calculate soluble fraction and apply to tracer
                          fsol = min(1.,max(0.,(1.-exp(-kc_scaled*dz/this_w_upd)))) !*factor_temp(ispc,i,k))
                          pw_up(ispc,i,k) = sc_up(ispc,i,k)*fsol

                       !-- this the 'alpha' parameter in Eq 8 of Mari et al (2000 JGR) = X_aq/X_total
                       else
                          fliq       = henry_coef*qrco(i,k) /(1.+henry_coef*qrco(i,k))
                          !---   aqueous-phase concentration in rain water
                          pw_up(ispc,i,k) = max(0.,sc_up(ispc,i,k)*(1.-exp(-fliq*kc*dz/w_upd)))!*factor_temp(ispc,i,k))
                       endif

                    endif

                    !---(in cloud) total mixing ratio in gas and aqueous phases
                    sc_up(ispc,i,k) = sc_up(ispc,i,k) - pw_up(ispc,i,k)

                    !
                    !---(in cloud)  mixing ratio in aqueous phase
                    !sc_up_aq(ispc,i,k) = conc_mxr(ispc) !if using set to zero at the begin.
                endIF
            enddo
            !
            !-- total aerosol/gas in the rain water
            dp=100.*(po_cup(i,k)-po_cup(i,k+1))

            tot_pw_up_chem(:,i) = tot_pw_up_chem(:,i) + pw_up(:,i,k)*dp/g
          enddo loopk
           !
     enddo
  end subroutine get_incloud_sc_chem_up
!---------------------------------------------------------------------------------------------------
  FUNCTION henry(ispc,temp,rhoair) RESULT(henry_coef)
  !--- calculate Henry's constant for solubility of gases into cloud water
  !--- inputs : ak0(ispc), dak(ispc),  hstar(ispc), dhr(ispc)
  implicit none
  integer, intent(in) :: ispc
  real   , intent(in) :: temp,rhoair
  real :: henry_coef
  real :: fct ,tcorr, corrh

  !--- define some constants!
  real, parameter:: rgas  =  8.205e-2 ! atm M^-1 K^-1 ! 8.314 gas constant [J/(mol*K)]
  real, parameter:: temp0i= 1./298.15! inverse of standard temperature [K]
  real, parameter:: hplus = 1.175E-4          !  for cloud water. pH is asuumed to be 3.93: pH=3.93 =>hplus=10**(-pH)

  ! aqueous-phase concentrations XXXa [mol/m3(air)]!
  ! gas-phase concentrations XXXg [mol/m3(air)]!
  ! Henry constants XXXh for scavenging [mol/(l*atm)]!
  ! converted to [(mol(aq)/m3(aq))/(mol(g)/m3(air))], i.e. dimensionless!
  ! in equilibrium XXXa = XXXh * LWC * XXXg!
  tcorr = 1./temp - temp0i
  fct   = rgas * temp ! - for henry_coef dimensioless

  !-taking into account the acid dissociation constant
  ! ak=ak0*exp(dak*(1/t-1/298))
  corrh=1.+CNV_Tracers(ispc)%Vect_Hcts(3) * exp(CNV_Tracers(ispc)%Vect_Hcts(4) * tcorr)/hplus

  !-- for concentration in mol[specie]/mol[air] - Eq 5 in 'Compilation of Henry's law constants (version 4.0) for
  !-- water as solvent, R. Sander, ACP 2015'.
  henry_coef =  CNV_Tracers(ispc)%Vect_Hcts(1) * exp(CNV_Tracers(ispc)%Vect_Hcts(2)*tcorr) * fct * corrh

  end FUNCTION henry

!---------------------------------------------------------------------------------------------------
  FUNCTION henry_gcc(ispc,temp,rhoair) RESULT(henry_coef)
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Return henry coefficient (liquid to gas) as defined by GEOS-Chem
    !EOP  
    !=====================================================================================
    implicit none
    integer, intent(in) :: ispc
    real   , intent(in) :: temp,rhoair
    real                :: henry_coef          ! effective gas/aq constant [-] 
    ! parameter
    real*8, parameter :: pH   = 4.5d0
    REAL*8, PARAMETER :: TREF = 298.15d0        ! [K          ]
    REAL*8, PARAMETER :: R    = 8.3144598d0     ! [J K-1 mol-1]
    REAL*8, PARAMETER :: ATM  = 101.325d0       ! [mPa (!)    ]
    ! local variables
    real*8          :: hstar8, dhr8, ak08, temp8, h8 !, dak8
    ! cast all variables to r*8 locally to prevent overflows
    hstar8 = CNV_Tracers(ispc)%Vect_Hcts(1) ! Henry coefficient [M/atm]]
    dhr8   = CNV_Tracers(ispc)%Vect_Hcts(2) ! temperature dependency of hstar [-d ln(kH) / d(1/T)] 
    ak08   = CNV_Tracers(ispc)%Vect_Hcts(3) ! pKa value [-]
    !dak8   = CNV_Tracers(ispc)%Vect_Hcts(4) ! temperature dependency of ak0, currently not used
    temp8  = temp
    ! calculate henry coefficient
    h8 = hstar8 * exp ( dhr8 * (1./temp8 - 1./TREF) ) * R * temp8 / ATM
    if ( ak08 > 0.0d0 ) then
       h8 = h8 * ( 1.0 + 10.0**(pH-ak08) )
    endif       
    ! limit henry coefficient to 1.0e30
    henry_coef = real( min(h8,1.0d+30) )

  END FUNCTION henry_gcc
!---------------------------------------------------------------------------------------------------
  subroutine get_incloud_sc_chem_dd(cumulus,mtp,se,se_cup,sc_dn,pw_dn ,pw_up,sc_up                          &
                                   ,tot_pw_up_chem,tot_pw_dn_chem                                                 &
                                   ,z_cup,rho,po_cup,qrcdo,pwdo,pwevo,edto,zdo,dd_massentro,dd_massdetro,pwavo,pwo &
                                   ,jmin,ierr,itf,ktf,its,ite, kts,kte)
     IMPLICIT NONE
     !-inputs
     integer                               ,intent (in)  :: itf,ktf, its,ite, kts,kte
     integer                               ,intent (in)  :: mtp
     integer, dimension (its:ite)          ,intent (in)  :: ierr ,jmin
     character *(*)                        ,intent (in)  :: cumulus
     real, dimension (mtp ,its:ite,kts:kte),intent (in)  :: se,se_cup,pw_up,sc_up
     real, dimension (its:ite)             ,intent (in)  :: edto,pwavo,pwevo
     real, dimension (its:ite,kts:kte)     ,intent (in)  :: z_cup,rho,po_cup&
                       ,qrcdo,pwdo,zdo,dd_massentro,dd_massdetro,pwo
     real, dimension (mtp ,its:ite        ),intent (in)  :: tot_pw_up_chem

     !-outputs
     real, dimension (mtp ,its:ite,kts:kte),intent (out) :: sc_dn,pw_dn
     real, dimension (mtp ,its:ite        ),intent (out) :: tot_pw_dn_chem

     !-locals
     real   , dimension (mtp) :: conc_mxr
     real :: x_add,dz,XZZ,XZD,XZE,denom, evaporate,pwdper,x1,frac_evap,dp,xkk
     integer :: i,k,ispc

     sc_dn          = 0.0
     pw_dn          = 0.0
     tot_pw_dn_chem = 0.0
     if(trim(cumulus) == 'shallow') return

     do i=its,itf
       if(ierr(i) /= 0) cycle

       !--- fration of the total rain that was evaporated
       frac_evap = - pwevo(i)/(1.e-16+pwavo(i))

       !--- scalar concentration in-cloud - downdraft

       !--- at k=jmim
       k=jmin(i)
       pwdper = pwdo(i,k)/(1.e-16+pwevo(i)) *frac_evap  ! > 0
       if(USE_TRACER_EVAP == 0 ) pwdper = 0.0

       dp= 100.*(po_cup(i,k)-po_cup(i,k+1))

       do ispc=1,mtp
        !--downdrafts will be initiate with a mixture of 50% environmental and in-cloud concentrations
              sc_dn(ispc,i,k) = se_cup(ispc,i,k)
              pw_dn(ispc,i,k) = - pwdper * tot_pw_up_chem(ispc,i)*g/dp
              sc_dn(ispc,i,k) = sc_dn(ispc,i,k) - pw_dn(ispc,i,k)
              tot_pw_dn_chem(ispc,i) = tot_pw_dn_chem(ispc,i) + pw_dn(ispc,i,k)*dp/g
       enddo
       !
       !--- calculate downdraft mass terms
       do k=jmin(i)-1,kts,-1
         XZZ=              zdo(i,k+1)
         XZD= 0.5*dd_massdetro(i,k  )
         XZE=     dd_massentro(i,k  )

         denom =  (XZZ-XZD+XZE)
         !-- transport + mixing
         if(denom > 0.) then
            sc_dn(:,i,k) = (sc_dn(:,i,k+1)*XZZ - sc_dn(:,i,k+1)*XZD + se(:,i,k)*XZE) / denom
         else
            sc_dn(:,i,k) = sc_dn(:,i,k+1)
         endif
        !
        !-- evaporation term
        if(USE_TRACER_EVAP == 0 )cycle

           dp= 100.*(po_cup(i,k)-po_cup(i,k+1))

           !-- fraction of evaporated precip per layer
           pwdper   = pwdo(i,k)/(1.e-16+pwevo(i))! > 0

           !-- fraction of the total precip that was actually evaporated at layer k
           pwdper   = pwdper * frac_evap

           !-- sanity check
           pwdper   = min(1.,max(pwdper,0.))

           do ispc=1,mtp
              !-- amount evaporated by the downdraft from the precipitation
              pw_dn(ispc,i,k) = - pwdper * tot_pw_up_chem (ispc,i)*g/dp ! < 0. => source term for the downdraft tracer concentration

              !-- final tracer in the downdraft
              sc_dn(ispc,i,k) = sc_dn(ispc,i,k) - pw_dn(ispc,i,k) ! observe that -pw_dn is > 0.

              !-- total evaporated tracer
              tot_pw_dn_chem(ispc,i) = tot_pw_dn_chem(ispc,i) + pw_dn(ispc,i,k)*dp/g
           enddo
        enddo
        !
     enddo
   end subroutine get_incloud_sc_chem_dd
!---------------------------------------------------------------------------------------------------
   subroutine get_cloud_bc(cumulus,kts,kte,ktf,xland,po,array,x_aver,k22,add,Tpert)
    implicit none
    character *(*)   ,intent (in) :: cumulus
    integer,intent(in)            :: kts,kte,ktf,k22
    real   ,intent(in)            :: array(kts:kte),po(kts:kte),xland
    real   ,optional ,intent(in)  :: add
    real   ,optional ,intent(in)  :: Tpert(kts:kte)
    real   ,intent(out)           :: x_aver
    integer                       :: i,local_order_aver,order_aver, i_beg,i_end,ic
    real,    parameter            :: frac_ave_layer_ocean= 0.3
    real                          :: count,dp,dp_layer,effec_frac,x_ave_layer

    !-- dimensions of the average:
    !-- a) to pick the value at k22 level, instead of an average between
    !--    (k22-order_aver, ..., k22-1, k22) set order_aver=kts
    !-- b) to average between kts and k22 => set order_aver = k22
    !order_aver = 4    !=> bc_meth 0: average between k22, k22-1, k22-2 ...
                       !=> bc_meth 1: average between ... k22+1,k22, k22-1 ...
    !-- order_aver = kts !=> average between k22, k22-1 and k22-2

     if(bc_meth == 0) then

      order_aver = 3
      local_order_aver=min(k22,order_aver)

      x_aver=0.
      do i = kts,local_order_aver
        x_aver = x_aver + array(k22-i+1)
      enddo
      x_aver = x_aver/float(local_order_aver)

     elseif(bc_meth == 1) then
      effec_frac  = (1.-xland) +xland*frac_ave_layer_ocean
      if(trim(cumulus) == 'deep'   ) x_ave_layer = cum_ave_layer(deep)*effec_frac
      if(trim(cumulus) == 'shallow') x_ave_layer = cum_ave_layer(shal)*effec_frac
      if(trim(cumulus) == 'mid'    ) x_ave_layer = cum_ave_layer(mid )*effec_frac

      i_beg = minloc(abs(po(kts:ktf)-(po(k22)+0.5*x_ave_layer)),1)
      i_end = minloc(abs(po(kts:ktf)-(po(k22)-0.5*x_ave_layer)),1)
      i_beg = min(ktf,max(i_beg,kts))
      i_end = min(ktf,max(i_end,kts))

      if(i_beg >= i_end) then
         x_aver   = array(k22)
         dp_layer = 0.
         ic       = i_beg

      else
         dp_layer = 1.e-06
         x_aver   = 0.
         ic       = 0
         do i = i_beg,ktf
              dp = -(po(i+1)-po(i))
              if(dp_layer + dp <= x_ave_layer)  then
                  dp_layer =  dp_layer  + dp
                  x_aver   =  x_aver    + array(i)*dp

              else
                  dp       =  x_ave_layer - dp_layer
                  dp_layer =  dp_layer    + dp
                  x_aver   =  x_aver      + array(i)*dp

                  exit
              endif
      enddo
      x_aver = x_aver/dp_layer
         ic  = max(i_beg,i)
      endif
      !print*,"xaver1=",real(x_aver,4),real(dp_layer,4)

      !-- this perturbation is included only for MSE
      if(present(Tpert)) x_aver = x_aver + cp*maxval(Tpert(i_beg:ic))  ! version 2 - maxval in the layer

    endif
    IF(present(add)) x_aver = x_aver + add

   end subroutine get_cloud_bc
!-----------------------------------------------------------------------------------------
  FUNCTION get_w_upd_gcc( vud, xland, online_vud ) RESULT( w_upd )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Return updraft vertical velocity to be used for GEOS-Chem convective washout. 
    !EOP  
    !=====================================================================================
    implicit none
    real, intent(in) :: vud        ! online updraft velocity [m/s]
    real, intent(in) :: xland      ! land flag (1.-FRLAND): greater value means more water 
    real, intent(in) :: online_vud ! use online vud (1.0) or set vud based on land/water (0.0)
    real             :: w_upd      ! updraft velocity to use 
    ! use environment vud if specified so 
    if ( online_vud == 1.0 ) then
       w_upd = vud
    ! use parameterization otherwise: 10m/s over land, 5m/s over water.
    else
       ! over water
       if ( xland > 0.9 ) then
          w_upd = 5.0
       ! over land
       else
          w_upd = 10.0
       endif
    endif

  END FUNCTION get_w_upd_gcc
!---------------------------------------------------------------------------------------------------
  FUNCTION gcc_e_ice( temp ) RESULT( vpress )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Calculate saturation vapor pressure over ice at given temperature - adapted from GEOS-Chem
    !  Marti & Mauersberber (GRL '93) formulation of saturation
    !  vapor pressure of ice [Pa] is: log P = A/TK + B
    !EOP
    !=====================================================================================
    implicit none
    real, intent(in) :: temp   ! Temperature [K]
    real             :: vpress ! Saturation vapor pressure [hPa]
    ! parameter
    real, parameter  :: A = -2663.5
    real, parameter  :: B =  12.537
    ! Saturation vap press of Ice [Pa]
    if ( temp <= 1.e-5 ) then
       vpress = 0.0
    else
       vpress = ( 10.**( A/temp + B ) )
    endif

  END FUNCTION gcc_e_ice
!---------------------------------------------------------------------------------------------------
  SUBROUTINE compute_ki_gcc_aerosol( ispc, temp, kc_scaled )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    ! Compute the loss rate of a GEOS-Chem aerosol. This follows the parameterization described 
    ! in Jacob et al. 2000, as implemented in GEOS-Chem.
    !EOP
    !=====================================================================================
    implicit none

    integer, intent(in)      :: ispc
    real,    intent(in)      :: temp               ! temperature [K]
    real,    intent(out)     :: kc_scaled          ! loss rate [s-1]
    ! local variables
    real                     :: kcscal1            ! temperature-dependent scale factor for temperature range 1 
    real                     :: kcscal2            ! temperature-dependent scale factor for temperature range 2 
    real                     :: kcscal3            ! temperature-dependent scale factor for temperature range 3 
    ! parameter
    real, parameter          :: TEMP1  = 237.0        ! K
    real, parameter          :: TEMP2  = 258.0        ! K

    ! Get scale factors
    kcscal1 = CNV_Tracers(ispc)%KcScal(1) 
    kcscal2 = CNV_Tracers(ispc)%KcScal(2) 
    kcscal3 = CNV_Tracers(ispc)%KcScal(3) 
    ! start with default kc, then scale based on temperature and aerosol-specific scale factor
    kc_scaled = KC_DEFAULT_GCC
    if ( temp < TEMP1 ) then
       kc_scaled = kc_scaled * kcscal1
    else if ( (temp>=TEMP1) .and. (temp<TEMP2) ) then
       kc_scaled = kc_scaled * kcscal2
    else
       kc_scaled = kc_scaled * kcscal3
    endif

  END SUBROUTINE compute_ki_gcc_aerosol
!---------------------------------------------------------------------------------------------------
   SUBROUTINE compute_ki_gcc_gas( ispc, temp, press, q, cldh2o, Heff, kc_scaled, l2g )
    !=====================================================================================
    !BOP
    ! !DESCRIPTION:
    !  Compute the loss rate of a GEOS-Chem gas. This follows the parameterization described 
    !  in Jacob et al. 2000, as implemented in GEOS-Chem.
    !EOP
    !=====================================================================================
     implicit none

     integer, intent(in)      :: ispc
     real,    intent(in)      :: temp               ! temperature [K]
     real,    intent(in)      :: press              ! pressure [Pa] 
     real,    intent(in)      :: q                  ! water vapor mixing ratio [kg/kg] 
     real,    intent(in)      :: cldh2o             ! cloud total water [kg/kg] 
     real,    intent(in)      :: Heff               ! effective gas/aq Henry constant [-]
     real,    intent(out)     :: kc_scaled          ! loss rate [s-1]
     real,    intent(out)     :: l2g                ! liquid to gas ratio 

     ! parameter
     real, parameter       :: T_zero = 273.16  ! K, as in ConvPar_GF_GEOS5 
     real, parameter       :: T_ice  = 250.16  ! K, as in ConvPar_GF_GEOS5
     real, parameter       :: TEMP3  = 248.0   ! K
     real, parameter       :: TEMP4  = 268.0   ! K

     ! local variables
     real            :: fract_liq_f
     real            :: cldliq, cldice, c_h2o
     real            :: i2g
     real            :: c_tot, f_l, f_i
     real            :: airdens

     ! Compute environmental variables: cloud water and H2O mixing ratio.
     ! This corresponds to the computations done in SETUP_WETSCAV in module
     ! wetscav_mod.F90 in GEOS-Chem.

     ! compute cloud liquid water content and cloud ice water. 
     ! Compute either based on environmental variables or use original GEOS-Chem formulation
     if ( CNV_Tracers(ispc)%online_cldliq == 1.0 ) then
        ! compute from cloud total water, using formulation as suggested by Saulo Freitas
        fract_liq_f = min(1., (max(0.,(temp-T_ice))/(T_zero-T_ice))**2)
        ! liquid and ice water in kg/kg 
        cldliq = cldh2o * fract_liq_f       ! kg/kg
        cldice = cldh2o * (1.-fract_liq_f)  ! kg/kg
        ! to convert to cm3/cm3, need air density
        airdens = 100.*press/(287.04*temp*(1.+0.608*q))
        cldliq  = cldliq*airdens*1.e-3      ! cm3/cm3
        cldice  = cldice*airdens*1.e-3      ! cm3/cm3

     else
        ! original GEOS-Chem formulation
        IF ( temp >= TEMP4 ) THEN
           cldliq = 1e-6
        ELSE IF ( temp > TEMP3 .and. temp < TEMP4 ) THEN
           cldliq = 1e-6 * ((temp-TEMP3)/(TEMP4-TEMP3))
        ELSE
           cldliq = 0.0
        ENDIF
        cldliq = max(cldliq,0.0)      ! cm3 H2O/cm3 air
        cldice = max(1e-6-cldliq,0.0) ! cm3 ice/cm3 air
     endif

     ! mixing ratio of H2O [v/v]: compute using Dalton's law
     c_h2o = gcc_e_ice(temp) / press

     ! ice to gas ratio 
     i2g = 0.0
     if ( (CNV_Tracers(ispc)%liq_and_gas==1.) .and. (c_h2o>0.0) ) then
        i2g = ( cldice / c_h2o ) * CNV_Tracers(ispc)%convfaci2g
     endif

     ! liquid to gas ratio
     l2g = Heff * cldliq

     ! fraction of species in liquid & ice phases (Eqs. 4, 5, 6, Jacob et al, 2000)
     c_tot = 1.0 + l2g + i2g
     f_l   = l2g / c_tot
     f_i   = i2g / c_tot

     ! compute the rate constant Ki for loss of species from
     ! convective updraft scavenging (Eq. 1, Jacob et al, 2000)
     if ( temp >= TEMP4 ) then
        kc_scaled = KC_DEFAULT_GCC * ( f_l + f_l )
     else if ( temp > TEMP3 .and. temp < TEMP4 ) THEN
        kc_scaled = KC_DEFAULT_GCC * ( ( CNV_Tracers(ispc)%retfactor * f_l ) + f_i )
     else
        kc_scaled = KC_DEFAULT_GCC * f_i
     endif

   END SUBROUTINE compute_ki_gcc_gas

end module ConvPar_GF_SharedParams
