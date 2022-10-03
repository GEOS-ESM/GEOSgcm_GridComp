MODULE ConvPar_GF_SharedParams

USE MAPL
USE GEOSmoist_Process_Library

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

 LOGICAL :: CNV_2MOM = .FALSE.

 CHARACTER(len=10) :: GF_ENV_SETTING = 'DYNAMICS' ! or 'CURRENT'

 INTEGER, PARAMETER :: MAX_NSPEC=200
 CHARACTER(len=100),DIMENSION(MAX_NSPEC)    ::  CHEM_NAME
 INTEGER           ,DIMENSION(MAX_NSPEC)    ::  CHEM_NAME_MASK,CHEM_NAME_MASK_EVAP
 REAL              ,DIMENSION(MAX_NSPEC)    ::  CHEM_ADJ_AUTOC

 CONTAINS

  subroutine get_incloud_sc_chem_up(cumulus,mtp,se,se_cup,sc_up,pw_up,tot_pw_up_chem&
                                   ,z_cup,rho,po,po_cup  &
                                   ,qrco,tempco,pwo,zuo,up_massentro,up_massdetro,vvel2d,vvel1d  &
                                   ,start_level,k22,kbcon,ktop,klcl,ierr,xland,itf,ktf,its,ite, kts,kte)
     IMPLICIT NONE
     !-inputs
     integer                               ,intent (in)  :: itf,ktf, its,ite, kts,kte
     integer                               ,intent (in)  :: mtp
     integer, dimension (its:ite)          ,intent (in)  :: ierr ,kbcon,ktop,k22,klcl,start_level
     character *(*)                        ,intent (in)  :: cumulus
     real, dimension (mtp ,its:ite,kts:kte),intent (in)  :: se,se_cup
     real, dimension (its:ite,kts:kte)     ,intent (in)  :: z_cup,rho,po_cup,qrco,tempco,pwo,zuo &
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

     !--initialization
     sc_up          = se_cup
     pw_up          = 0.0
     tot_pw_up_chem = 0.0

     IF(USE_TRACER_SCAVEN==2 .and. cumulus /= 'shallow') THEN
        factor_temp = 1.
        do i=its,itf
            if(ierr(i) /= 0) cycle
            do ispc = 1,mtp
               ! - if tracer is type "carbon" then set coefficient to 0 for hydrophobic
               if( TRIM(CHEM_name (ispc)(1:len_trim('OCphobic') )) == 'OCphobic') factor_temp(ispc,:,:) = 0.0

               ! - suppress scavenging most aerosols at cold T except BCn1 (hydrophobic), dust, and HNO3
               if( TRIM(CHEM_name (ispc)(1:len_trim('BCphobic') )) == 'BCphobic') then
                  where(tempco < 258.) factor_temp(ispc,:,:) = 0.0
               endif

               if( TRIM(CHEM_name (ispc)) == 'sulfur'   .or. &

                   TRIM(CHEM_name (ispc)(1:len_trim('ss') )) == 'ss'  .or. & ! 'seasalt'

                   TRIM(CHEM_name (ispc)) == 'SO2'      .or. &
                   TRIM(CHEM_name (ispc)) == 'SO4'      .or. &

                   TRIM(CHEM_name (ispc)) == 'nitrate'  .or. &
                   TRIM(CHEM_name (ispc)) == 'bromine'  .or. &
                   TRIM(CHEM_name (ispc)) == 'NH3'      .or. &
                   TRIM(CHEM_name (ispc)) == 'NH4a'          ) then

                   where(tempco < 258.) factor_temp(ispc,:,:) = 0.0
               endif

             enddo
          enddo
     endIF
     do i=its,itf
       if(ierr(i) /= 0) cycle
      !start_level(i) = klcl(i)
      !start_level(i) = k22(i)

       do ispc=1,mtp
             call get_cloud_bc(cumulus,kts,kte,ktf,xland(i),po(i,kts:kte),se_cup(ispc,i,kts:kte),sc_b(ispc,i),k22(i))
       enddo
       do k=kts,start_level(i)
            sc_up   (:,i,k) = sc_b  (:,i)
           !sc_up   (:,i,k) = se_cup(:,i,k)
       enddo
     enddo

     do i=its,itf
          if(ierr(i) /= 0) cycle
loopk:      do k=start_level(i)+1,ktop(i)+1

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
            if(USE_TRACER_SCAVEN==0 .or. cumulus == 'shallow') cycle loopk
            dz=z_cup(i,k)-z_cup(i,k-1)

            !-- in-cloud vert velocity for scavenging formulation 2
!           w_upd = cte_w_upd
!           w_upd = vvel1d(i)
            w_upd = vvel2d(i,k)

            do ispc = 1,mtp
                IF(CNV_Tracers(ispc)%fscav > 1.e-6) THEN ! aerosol scavenging

                    !--formulation 1 as in GOCART with RAS conv_par
                    if(USE_TRACER_SCAVEN==1) &
                    pw_up(ispc,i,k) = max(0.,sc_up(ispc,i,k)*(1.-exp(- CNV_Tracers(ispc)%fscav  * (dz/1000.))))

                    !--formulation 2 as in GOCART
                    if(USE_TRACER_SCAVEN==2) &
                    pw_up(ispc,i,k) = max(0.,sc_up(ispc,i,k)*(1.-exp(- CHEM_ADJ_AUTOC(ispc) * kc * (dz/w_upd)))*factor_temp(ispc,i,k))

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
                    henry_coef=henry(ispc,tempco(i,k),rho(i,k))

                    if(USE_TRACER_SCAVEN==3) then
                      !--- cloud liquid water tracer concentration
                      conc_mxr(ispc) = (henry_coef*qrco(i,k) /(1.+henry_coef*qrco(i,k)) )* sc_up(ispc,i,k)
                      !
                      !---   aqueous-phase concentration in rain water
                      pw_up(ispc,i,k) = conc_mxr(ispc)*pwo(i,k)/(1.e-8+qrco(i,k))

                    else

                      !-- this the 'alpha' parameter in Eq 8 of Mari et al (2000 JGR) = X_aq/X_total
                       fliq = henry_coef*qrco(i,k) /(1.+henry_coef*qrco(i,k))

                      !---   aqueous-phase concentration in rain water
                      pw_up(ispc,i,k) = max(0.,sc_up(ispc,i,k)*(1.-exp(-fliq* CHEM_ADJ_AUTOC(ispc) *kc*dz/w_upd)))!*factor_temp(ispc,i,k))

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
           !----- get back the in-cloud updraft gas-phase mixing ratio : sc_up(ispc,k)
!          do k=start_level(i)+1,ktop(i)+1
!            do ispc = 1,mtp
!             sc_up(ispc,i,k) = sc_up(ispc,i,k) - sc_up_aq(ispc,i,k)
!            enddo
!          enddo
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
  real, parameter:: avogad=  6.022e23! Avogadro constant [1/mol]
  real, parameter:: rhoH2O=  999.9668! density of water [kg/m3]
  real, parameter:: temp0 =    298.15! standard temperature [K]
  real, parameter:: temp0i= 1./298.15! inverse of standard temperature [K]
  real, parameter:: MWH2O =        18.02! molecular mass of water [kg/kmol]
  real, parameter:: MWAIR =        28.97! effective molecular mass of air [kg/kmol]
  real, parameter:: conv3 = avogad / 1.0e6!  [mol(g)/m3(air)]  to [molec(g)/cm3(air)]
  real, parameter:: conv4 = 100.          !  [m]          to [cm]
  real, parameter:: conv5 = 1000.          !  [m^3]            to [l]
  real, parameter:: conv7 = 1/conv5          !  [l]    to [m^3]
  real, parameter:: conv6 = 1. / 101325.  !  [Pa]            to [atm]
  real, parameter:: hplus = 1.175E-4          !  for cloud water. pH is asuumed to be 3.93: pH=3.93 =>hplus=10**(-pH)

  ! aqueous-phase concentrations XXXa [mol/m3(air)]!
  ! gas-phase concentrations XXXg [mol/m3(air)]!
  ! Henry constants XXXh for scavenging [mol/(l*atm)]!
  ! converted to [(mol(aq)/m3(aq))/(mol(g)/m3(air))], i.e. dimensionless!
  ! in equilibrium XXXa = XXXh * LWC * XXXg!
  tcorr = 1./temp - temp0i

  !-P. Colarco corrected the expression below
  !fct   = conv7 * rgas * temp ! - for henry_coef in units 1/m3
   fct   =         rgas * temp ! - for henry_coef dimensioless


  !-taking into account the acid dissociation constant
  ! ak=ak0*exp(dak*(1/t-1/298))
  corrh=1.+CNV_Tracers(ispc)%Vect_Hcts(3) * exp(CNV_Tracers(ispc)%Vect_Hcts(4) * tcorr)/hplus

  !-- for concentration in mol[specie]/mol[air] - Eq 5 in 'Compilation of Henry's law constants (version 4.0) for
  !-- water as solvent, R. Sander, ACP 2015'.
  henry_coef =  CNV_Tracers(ispc)%Vect_Hcts(1) * exp(CNV_Tracers(ispc)%Vect_Hcts(2)*tcorr) * fct * corrh

  end FUNCTION henry
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
     if(cumulus == 'shallow') return

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
             !sc_dn(ispc,i,k) = 0.9*se_cup(ispc,i,k)+0.1*sc_up(ispc,i,k)

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

              !if(ispc==1) print*,"pw=",pwdper,tot_pw_up_chem (ispc,i),pwevo(i)/pwavo(i),pwdo(i,k)/(1.e-16+pwo(i,k))

              !-- final tracer in the downdraft
              sc_dn(ispc,i,k) = sc_dn(ispc,i,k) - pw_dn(ispc,i,k) ! observe that -pw_dn is > 0.

              !-- total evaporated tracer
              tot_pw_dn_chem(ispc,i) = tot_pw_dn_chem(ispc,i) + pw_dn(ispc,i,k)*dp/g

              !print*,"to=",k,tot_pw_dn_chem(ispc,i),pwdo(i,k)/(1.e-16+pwevo(i)),frac_evap,tot_pw_dn_chem(ispc,i)/tot_pw_up_chem (ispc,i)

           enddo
        enddo
        !
     enddo
    end subroutine get_incloud_sc_chem_dd
!---------------------------------------------------------------------------------------------------
    subroutine interface_aerchem(mtp,aer_chem_mech)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mtp
      CHARACTER(len=*), INTENT(IN) :: AER_CHEM_MECH
      !-local vars
      INTEGER :: ispc, len_ACM, len_spc, irun=0
      CHARACTER(len=100) :: TMP_AER_NAME
      CHEM_NAME_MASK      = 1
      CHEM_NAME_MASK_EVAP = 1
      CHEM_ADJ_AUTOC      = 1.0

      !- GOCART + PCHEM section
      IF(AER_CHEM_MECH=='GOCART') then
       len_ACM=len(TRIM(AER_CHEM_MECH))
       do ispc=1,mtp
           len_spc           = len(TRIM(CNV_Tracers(ispc)%cname))
           CHEM_name (ispc)  = TRIM(CNV_Tracers(ispc)%qname)
       enddo
      endIF
   end subroutine interface_aerchem
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
      x_ave_layer = ave_layer*effec_frac


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

end module ConvPar_GF_SharedParams
