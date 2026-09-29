module edmf_mod

  !
  ! Mass flux updraft parameterization implemented by Kay Suselj (JPL).
  ! Reference: Suselj et al 2021, DOI: 10.1175/MWR-D-20-0183.1
  ! Additional development by Nathan Arnold and David New (GMAO).
  !
  
  use MAPL_Constants, only: mapl_epsilon, mapl_grav, mapl_cp,  &
                               mapl_alhl, mapl_p00, mapl_vireps,  &
                               mapl_alhs, mapl_alhf, mapl_kappa,  &
                               mapl_pi, mapl_celsius_to_kelvin
  
  use MAPL_Mod,          only: mapl_undef
  
  use GEOS_Mod
  
  implicit none
  
  type EDMFPARAMS_TYPE
      logical :: DOTRACERS
      integer :: DISCRETE
      integer :: IMPLICIT
      integer :: ENTRAIN
      integer :: NUP
      integer :: ET
      integer :: UPABUOYDEP
      real    :: L0
      real    :: L0fac
      real    :: STOCHFRAC
      real    :: ENTUFAC
      real    :: ENT0
      real    :: ALPHATH
      real    :: ALPHAQT
      real    :: ALPHAW
      real    :: PWMAX
      real    :: PWMIN
      real    :: WA
      real    :: WB
      real    :: WC
      real    :: WCTHRESH
      real    :: MFLIMFAC
      real    :: ICE_RAMP
      real    :: PRCPCRIT
      real    :: TREFF
  endtype EDMFPARAMS_TYPE
  type (EDMFPARAMS_TYPE) :: MFPARAMS
  
  public run_edmf, mfparams
  
  contains
                                                       !==== Inputs ===================================================
  SUBROUTINE RUN_EDMF(its,ite, jts,jte, kts,kte, dt, & ! Index limits and timestep (s)
                      phis,                          & ! Surface geopotential (m2 s-2)
                      zlo3,                          & ! Surface-relative heights (m)
                      zw3,                           & ! Surface-relative edge heights (m)
                      pw3,                           & ! Edge pressures (Pa)
                      rhoe3,                         & ! Edge densities (kg m-3)
                      tke3,                          & ! Turbulent kinetic energy (m2 s-2)
                      u3,                            & ! U wind component (m s-1)
                      v3,                            & ! V wind component (m s-1)
                      t3,                            & ! Temperature (K)
                      thl3,                          & ! Liquid water potential temperature (K)
                      thv3,                          & ! Virtual potential temperature (K)
                      qv3,                           & ! Specific humidity (kg kg-1)
                      ql3,                           & ! Specific large-scale cloud liquid (kg kg-1)
                      qi3,                           & ! Specific large-scale cloud ice (kg kg-1)
                      exf3,                          & ! Full-level Exner function
                      exfh3,                         & ! Edge-level Exner function
                      wthl2,                         & ! Surface sensible heat flux (W m-2)
                      wqt2,                          & ! Surface evaporation (kg m-2 s-1)
                      frland,                        & ! Land fraction
                      pblh2,                         & ! PBL height (m)
                                                       !==== Outputs - variables needed for solver ====================
                      ae3,                           & ! Environmental area fraction
                      aw3,                           & ! Area-weighted updraft vertical velocity (m s-1)
                      aws3,                          & ! Dry static energy flux (J m s-1 kg-1)
                      awqv3,                         & ! Specific humidity flux (kg kg-1 m s-1)
                      awql3,                         & ! Specific liquid flux (kg kg-1 m s-1)
                      awqi3,                         & ! Specific ice flux (kg kg-1 m s-1)
                      awu3,                          & ! Kinematic U momentum flux (m2 s-2)
                      awv3,                          & ! Kinematic V momentum flux (m2 s-2)
                      YS,                            & ! Dry static energy increment for trisolver (J kg-1)
                      YQV,                           & ! Specific humidity increment for trisolver (kg kg-1)
                      YQL,                           & ! Liquid water increment for trisolver (kg kg-1)
                      YQI,                           & ! Ice water increment (kg kg-1)
                      YU,                            & ! U wind increment (m s-1)
                      YV,                            & ! V wind increment (m s-1)
                                                       !==== Outputs required for SHOC and ADG PDF ====================
                      mfw2,                          & ! Area-weighted vertical velocity squared (m2 s-2)
                      mfw3,                          & ! Area-weighted vertical velocity cubed (m3 s-3)
                      mfqt3,                         & ! Area-weighted updraft total water anomaly cubed (kg3 kg-3)
                      mfhl3,                         & ! Area-weighted updraft liquid static energy anomaly cubed (K3)
                      mfwqt,                         & ! Kinematic total water flux on midlevels (kg kg-1 m s-1)
                      mfhlqt,                        & ! Updraft total water-temperature anomaly covariance (kg kg-1 K)
                      mfwhl,                         & ! Kinematic static energy flux on midlevels (J kg-1 m s-1)
                      mftke,                         & ! Updraft kinetic energy (m2 s-2)
                      buoyf,                         & ! Updraft buoyancy flux (K m s-1)
                      edmfmf,                        & ! Updraft mass flux (kg m s-1)
                      dry_a3,                        & ! Dry updraft area fraction
                      moist_a3,                      & ! Cloudy updraft area fraction
                      dqrdt,                         & ! Liquid precipitation tendency
                      dqsdt,                         & ! Frozen precipitation tendency
                      edmf_dmf,                      & ! Detrained mass flux (kg m-2 s-1)
                                                       !==== Diagnostic outputs - updraft properties ==================
                      dry_w3,                        & ! Dry updraft mean vertical velocity (m s-1)
                      moist_w3,                      & ! Cloudy updraft mean vertical velocity (m s-1)
                      dry_qt3,                       & ! Dry updraft mean total water (kg kg-1)
                      moist_qt3,                     & ! Cloudy updraft mean total water (kg kg-1)
                      dry_thl3,                      & ! Dry updraft liquid water potential temperature (K)
                      moist_thl3,                    & ! Cloudy updraft liquid water potential temperature (K)
                      dry_u3,                        & ! Dry updraft mean u wind (m s-1)
                      moist_u3,                      & ! Cloudy updraft mean u wind (m s-1)
                      dry_v3,                        & ! Dry updraft mean v wind (m s-1)
                      moist_v3,                      & ! Cloudy updraft mean v wind (m s-1)
                      moist_qc3,                     & ! Cloudy updraft mean total condensate (kg kg-1)
                      entx,                          & ! Mean fractional lateral mixing rate
                      mfdepth,                       & ! Diagnostic updraft depth used in lateral mixing (m-1)
                      edmf_plumes_w,                 & ! Individual updraft vertical velocities (m s-1)
                      edmf_plumes_thl,               & ! Individual updraft liquid water potential temperatures (K)
                      edmf_plumes_qt )                 ! Individual updraft total waters (kg kg-1)
  
  
     INTEGER, INTENT(IN) :: ITS,ITE,JTS,JTE,KTS,KTE
     REAL,    INTENT(IN) :: DT
  
     REAL,DIMENSION(ITS:ITE,JTS:JTE,KTS:KTE), INTENT(IN) :: U3,    &
                                                            V3,    &
                                                            T3,    &
                                                            THL3,  &
                                                            THV3,  &
                                                            QV3,   &
                                                            QL3,   &
                                                            QI3,   &
                                                            ZLO3,  &
                                                            TKE3,  &
                                                            exf3
  
     REAL,DIMENSION(ITS:ITE,JTS:JTE,KTS-1:KTE), INTENT(IN) :: ZW3, PW3, rhoe3, exfh3
  
     REAL,DIMENSION(ITS:ITE,JTS:JTE), INTENT(IN) :: WTHL2,   &
                                                    WQT2,    &
                                                    PBLH2,   &
                                                    FRLAND,  &
                                                    PHIS
  
     ! Required outputs
     REAL,DIMENSION(ITS:ITE,JTS:JTE,KTS-1:KTE), INTENT(OUT) :: dry_a3,   &
                                                               moist_a3, &
                                                               ae3,      &
                                                               aw3,      &
                                                               aws3,     &
                                                               awqv3,    &
                                                               awql3,    &
                                                               awqi3,    &
                                                               awu3,     &
                                                               awv3,     &
                                                               edmfmf,   &
                                                               mfwhl,    &
                                                               mfwqt,    &
                                                               mftke
  
     REAL,DIMENSION(ITS:ITE,JTS:JTE,KTS:KTE), INTENT(OUT) :: buoyf, mfw2, mfw3,    &
                                                             mfqt3, mfhl3, mfhlqt, &
                                                             dqrdt, dqsdt
  
     REAL,DIMENSION(ITS:ITE,JTS:JTE,KTS:KTE), INTENT(INOUT) :: YS, YQV, YQL, &
                                                               YQI, YU, YV
  
    ! Diagnostic outputs
     REAL, DIMENSION(:,:),     POINTER :: mfdepth
  
     REAL, DIMENSION(:,:,:),   POINTER :: dry_w3,   moist_w3,   &
                                          dry_qt3,  moist_qt3,  &
                                          dry_thl3, moist_thl3, &
                                          dry_u3,   moist_u3,   &
                                          dry_v3,   moist_v3,   &
                                          moist_qc3, entx,      &
                                          edmf_dmf
  
     REAL, DIMENSION(:,:,:,:), POINTER :: EDMF_PLUMES_W,   &
                                          EDMF_PLUMES_THL, &
                                          EDMF_PLUMES_QT
  
  
  !============= Local variables =============
  
     ! updraft properties
     REAL,DIMENSION(KTS-1:KTE,1:MFPARAMS%NUP) :: UPW, UPTHL, UPQT, &
                                                 UPQL, UPQI, UPA,  &
                                                 UPU, UPV, UPTHV
     ! entrainment variables
     REAl,DIMENSION(KTS:KTE,1:MFPARAMS%NUP) :: ENT, ENTf
     INTEGER,DIMENSION(KTS:KTE,1:MFPARAMS%NUP) :: ENTi
     
     INTEGER :: K,KTMP,I,IH,JH,NUP2
     REAL :: wthv,wstar,qstar,thstar, &
             sigmaW,sigmaQT,sigmaTH,  &
             wmin,wmax,wlv,wtv,wp,    &
             B,QTn,THLn,THVn,QCn,QP,  &
             Un,Vn,Wn2,EntEXP,EntEXPU,&
             EntW,wf, WTHL, WQT, PBLH
  
     ! internal flipped variables (GEOS)
     REAL,DIMENSION(KTS:KTE)   :: U,V,THL,QT,THV,QV,QL,QI,ZLO,QR,QS
     REAL,DIMENSION(KTS-1:KTE) :: ZW,P,THLI,QTI
     REAL,DIMENSION(KTS-1:KTE) :: UI, VI, QVI, QLI, QII
  
     ! Updraft diagnostics
     REAL,DIMENSION(KTS-1:KTE) :: dry_a, moist_a,dry_w,moist_w,          &
                                  dry_qt,moist_qt,dry_thl,moist_thl,     &
                                  dry_u,moist_u,dry_v,moist_v, moist_qc
  
     REAL,DIMENSION(KTS-1:KTE) :: s_aw,s_aws,s_awqv,s_awql,s_awqi,s_awu,s_awv
     REAL,DIMENSION(KTS:KTE)   :: s_buoyf
     REAL,DIMENSION(KTS-1:KTE) :: s_aw2,s_aw3,s_aqt3,s_ahl3,s_aqt2,  &
                                  s_ahlqt,s_awqt,s_ahl2,s_awhl,qte
     REAL,DIMENSION(KTS:KTE)   :: exf,dp,pmid,wcfac
     REAL,DIMENSION(KTS-1:KTE) :: exfh,rhoe
  
     ! temporary/dummy variables
     REAL :: tmp, tmp2, sqrt2_sigmaW, tanh_term
  
     REAL :: L0,ztop,ltm,QTsrfF,THVsrfF,mft,mfthvt,mf,factor
     INTEGER, DIMENSION(2) :: the_seed
  
     LOGICAL :: calc_avg_diag
  
  
  ! min values to avoid singularities
     REAL, PARAMETER ::    &
         WSTARmin = 1.e-3, &
         PBLHmin = 100.
  
     ! If any average diagnostics requested, perform required calculations,
     ! otherwise skip for efficiency
     calc_avg_diag = (associated(dry_w3)   .or. associated(moist_w3)   .or. &
                      associated(dry_qt3)  .or. associated(moist_qt3)  .or. &
                      associated(dry_thl3) .or. associated(moist_thl3) .or. &
                      associated(dry_u3)   .or. associated(moist_u3)   .or. &
                      associated(dry_v3)   .or. associated(moist_v3)   .or. &
                      associated(moist_qc3) )
  
     ! set updraft properties to zero/undef
     dry_a3   = 0.
     moist_a3 = 0.
     if (calc_avg_diag) then
        if (associated(dry_w3))     dry_w3     = mapl_undef
        if (associated(moist_w3))   moist_w3   = mapl_undef
        if (associated(dry_qt3))    dry_qt3    = mapl_undef
        if (associated(moist_qt3))  moist_qt3  = mapl_undef
        if (associated(dry_thl3))   dry_thl3   = mapl_undef
        if (associated(moist_thl3)) moist_thl3 = mapl_undef
        if (associated(dry_u3))     dry_u3     = mapl_undef
        if (associated(moist_u3))   moist_u3   = mapl_undef
        if (associated(dry_v3))     dry_v3     = mapl_undef
        if (associated(moist_v3))   moist_v3   = mapl_undef
        if (associated(moist_qc3))  moist_qc3  = mapl_undef
     end if
  
     ! outputs - variables needed for solver
     aw3   =0.
     aws3  =0.
     awqv3 =0.
     awql3 =0.
     awqi3 =0.
     awu3  =0.
     awv3  =0.
     buoyf =0.
     mfw2  =0.
     mfw3  =0.
     mfqt3 =0.
     mfhl3 =0.
     mfwqt =0.
     mfhlqt=0.
     mfwhl =0.
     mftke =0.
     edmfmf=0.
     dqrdt =0.
     dqsdt =0.
  
     if (associated(entx)) entx = mapl_undef
  
     ! Initialize the environmental area. Updraft area will be subtracted below.
     ae3=1.
  
     ! OPTIMIZATION: Swap JH and IH loop order for stride-1 memory access!
     DO JH=JTS,JTE
      DO IH=ITS,ITE 
  
        wthl=wthl2(IH,JH)/mapl_cp
        wqt=wqt2(IH,JH)
        pblh=pblh2(IH,JH)
  
        pblh=max(pblh,pblhmin)
        wthv=wthl+mapl_epsilon*thv3(IH,JH,kte)*wqt
  
        ! calc average TKE below 100 m
        tmp = 0.
        tmp2 = 0.
        k = kte
        do  while (zlo3(IH,JH,k)<100. .and. k>1)
           tmp2 = tmp2 + (zw3(IH,JH,k-1)-zw3(IH,JH,k))
           tmp = tmp+tke3(IH,JH,k)*(zw3(IH,JH,k-1)-zw3(IH,JH,k))
           k = k-1
        end do
        tmp = tmp/tmp2  ! avg TKE
  
        ! Activate mass-flux only if positive surface buoyancy flux
        ! and mean TKE below 100m exceeds threshold (for stability)
        IF (wthv > 0.0 .and. tmp>0.05 .and. phis(IH,JH).lt.3e4) then
  
        nup2 = MFPARAMS%NUP
  
        UPW=0.
        UPTHL=0.
        UPTHV=0.
        UPQT=0.
        UPA=0.
        UPU=0.
        UPV=0.
        UPQI=0.
        UPQL=0.
        ENT=0.
        QR = 0.
        QS = 0.
  
        ! Estimate scale height for entrainment calculation
        if (mfparams%ET == 2 ) then
           pmid = 0.5*(pw3(IH,JH,kts-1:kte-1)+pw3(IH,JH,kts:kte))
           call calc_mf_depth(kts,kte,t3(IH,JH,:),zlo3(IH,JH,:)-zw3(IH,JH,kte),qv3(IH,JH,:),pmid,ztop,wthv,wqt)
           L0 = max(min(ztop,2500.),500.) / mfparams%L0fac
           if (associated(mfdepth)) mfdepth(IH,JH) = ztop
        else ! if mfparams%ET not 2
           L0 = mfparams%L0
        end if
  
        if (ztop.gt.100.) then
  
        !
        ! flipping variables
        !
        DO k=kts,kte
          exf(k) = exf3(IH,JH,kte-k+kts)
          zlo(k)=zlo3(IH,JH,kte-k+kts)-zw3(IH,JH,kte)
          u(k)=u3(IH,JH,kte-k+kts)
          v(k)=v3(IH,JH,kte-k+kts)
          thl(k)=thl3(IH,JH,kte-k+kts)
          thv(k)=thv3(IH,JH,kte-k+kts)
          qv(k)=qv3(IH,JH,kte-k+kts)
          ql(k)=ql3(IH,JH,kte-k+kts)
          qi(k)=qi3(IH,JH,kte-k+kts)
        END DO
        if (MFPARAMS%DISCRETE == 0) then
           DO k=kts,kte-1
              ui(k)   = 0.5*( u3(IH,JH,kte-k+kts)   + u3(IH,JH,kte-k+kts-1) )
              vi(k)   = 0.5*( v3(IH,JH,kte-k+kts)   + v3(IH,JH,kte-k+kts-1) )
              thli(k) = 0.5*( thl3(IH,JH,kte-k+kts) + thl3(IH,JH,kte-k+kts-1) )
              qvi(k)  = 0.5*( qv3(IH,JH,kte-k+kts)  + qv3(IH,JH,kte-k+kts-1) )
              qli(k)  = 0.5*( ql3(IH,JH,kte-k+kts)  + ql3(IH,JH,kte-k+kts-1) )
	      qii(k)  = 0.5*( qi3(IH,JH,kte-k+kts)  + qi3(IH,JH,kte-k+kts-1) )              
           END DO
        else
           DO k=kts,kte-1
              ui(k)   = u3(IH,JH,kte-k+kts-1)
              vi(k)   = v3(IH,JH,kte-k+kts-1)
              thli(k) = thl3(IH,JH,kte-k+kts-1)
              qvi(k)  = qv3(IH,JH,kte-k+kts-1)
              qli(k)  = ql3(IH,JH,kte-k+kts-1)
              qii(k)  = qi3(IH,JH,kte-k+kts-1)
           END DO
        end if
        ui(kte)     = u(kte)
        vi(kte)     = v(kte)
        thli(kte)   = thl(kte)
        qvi(kte)    = qv(kte)
        qli(kte)    = ql(kte)
        qii(kte)    = qi(kte)
        ui(kts-1)   = u(kts)
        vi(kts-1)   = v(kts)
        thli(kts-1) = thl(kts)  ! approximate
        qvi(kts-1)  = qv(kts)
        qli(kts-1)  = ql(kts)
        qii(kts-1)  = qi(kts)
        qt  = qv+ql+qi
        qti = qvi+qli+qii
  
        DO k=kts-1,kte
          exfh(k) = exfh3(IH,JH,kte-k+kts-1)
          rhoe(k) = rhoe3(IH,JH,kte-k+kts-1)
          zw(k)   = zw3(IH,JH,kte-k+kts-1)-zw3(IH,JH,kte)
          p(k)    = pw3(IH,JH,kte-k+kts-1)
        ENDDO
  
        dp = p(kts-1:kte-1)-p(kts:kte)
  
        !
        ! compute entrainment coefficient
        !
  
        ! get dz/L0
        do i=1,Nup2
          do k=kts,kte
            ENTf(k,i)=((ZW(k)-ZW(k-1))/L0)
          enddo
        enddo
  
        ! get Poisson P(dz/L0)
        THE_SEED(1) = 1000000 * ( 100*thl(kts) - INT(100*thl(kts)))
        THE_SEED(2) = 1000000 * ( 100*thl(kts+1) - INT(100*thl(kts+1)))
        if(THE_SEED(1) == 0) THE_SEED(1) =  5
        if(THE_SEED(2) == 0) THE_SEED(2) = -5
  
        if (L0 .gt. 0. ) then
          ! entrainent: Ent=Ent0/dz*P(dz/L0)
          if (MFPARAMS%ENTRAIN==0 .or. MFPARAMS%ENTRAIN==4) then
            call Poisson(kts,kte,1,Nup2,ENTf,ENTi,the_seed)
            do i=1,Nup2
              do k=kts,kte
                 ENT(k,i) = (1.-MFPARAMS%STOCHFRAC) * MFPARAMS%Ent0/L0 &
                          + MFPARAMS%STOCHFRAC * real(ENTi(k,i))*MFPARAMS%Ent0/(ZW(k)-ZW(k-1))
                 ! Increase ent above 2500m to limit deepest plumes
                 ENT(k,i) = ENT(k,i)*(1.+MAX(0.,ZW(k)-2500.)/500.)
              enddo
            enddo
          else if (MFPARAMS%ENTRAIN==1 ) then
            call Poisson(kts,kte,1,Nup2,ENTf,ENTi,the_seed)
            do i=1,Nup2   ! Vary entrainment across updrafts, 0.75-1.25x
              do k=kts,kte
                ENT(k,i) = ((FLOAT(Nup2-i)*0.5/FLOAT(Nup2))+0.75)*( (1.-MFPARAMS%STOCHFRAC) * MFPARAMS%Ent0/L0 &
                        + MFPARAMS%STOCHFRAC * real(ENTi(k,i))*MFPARAMS%Ent0/(ZW(k)-ZW(k-1)) ) !&
              enddo
            enddo
          end if
  
        else ! if L0 <= 0
          ENT=0.
        end if
  
  
        !
        ! surface conditions
        !
        wstar=max(wstarmin,(mapl_grav*wthv*pblh/300.)**(1./3.))  ! convective velocity scale
        qstar=max(0.,wqt)/wstar
        thstar=max(0.,wthl)/wstar
  
        sigmaW=MFPARAMS%AlphaW*wstar
        sigmaQT=MFPARAMS%AlphaQT*qstar
        sigmaTH=MFPARAMS%AlphaTH*thstar
  
        wmin=sigmaW*MFPARAMS%pwmin
        wmax=sigmaW*MFPARAMS%pwmax
  
        ! Identify inversions below 1.5km, calculate stability in overlying 1km to define
        ! a dynamic pressure deceleration factor in the updraft w equation below.
        ! OPTIMIZATION: Used local 1D t(k) instead of strided 3D t3()
        wcfac = 0.
        tmp = 0.
        k = kts+1
        do while (zlo(k).lt.1500.)
           if ( t3(IH,JH,kte-k) > t3(IH,JH,kte-k+1) ) then
              tmp = thv(k)   ! THV at inversion
              exit
           end if
           k = k+1
        end do
        if (tmp.ne.0.) then
           ktmp = k
           do while (zlo(ktmp).lt.zlo(k)+1e3)
              ktmp = ktmp+1
           end do
           wcfac(1:k) = min(10.,max(0.,thv(ktmp)-thv(k)-MFPARAMS%WCTHRESH))*exp(-(zlo(k)-zlo(1:k))/200. )
        end if
  
        ! precalculate loop invariant variables
        sqrt2_sigmaW = sqrt(2.)*sigmaW
        if (MFPARAMS%UPABUOYDEP/=0) then
           tanh_term = 0.5 + 0.5*TANH((wthv-0.02)/0.09)
        else
           tanh_term = 1.0
        end if
  
        ! define surface conditions
        DO I=1,NUP2
  
          wlv=wmin+(wmax-wmin)/(real(NUP2))*(real(i)-1.)
          wtv=wmin+(wmax-wmin)/(real(NUP2))*real(i)
  
          UPW(kts-1,I)=min(0.5*(wlv+wtv), 5.)
          if (MFPARAMS%UPABUOYDEP/=0) then
             UPA(kts-1,I)=tanh_term*(0.5*ERF(wtv/sqrt2_sigmaW)-0.5*ERF(wlv/sqrt2_sigmaW))
          else
             UPA(kts-1,I)=(0.5*ERF(wtv/sqrt2_sigmaW)-0.5*ERF(wlv/sqrt2_sigmaW))
          end if
  
          UPU(kts-1,I)=U(kts)
          UPV(kts-1,I)=V(kts)
  
          UPQT(kts-1,I)=QT(kts)+0.32*UPW(kts-1,I)*sigmaQT/sigmaW
          UPTHV(kts-1,I)=THV(kts)+0.58*UPW(kts-1,I)*sigmaTH/sigmaW
  
        ENDDO ! NUP
  
        !
        ! If needed, rescale UPW to ensure that the mass-flux does not exceed layer mass
        !
  
        mf = SUM(RHOE(kts-1)*UPA(kts-1,:)*UPW(kts-1,:))
        factor = dp(kts)/(mf*MAPL_GRAV*dt)
        if (factor .lt. 1.0) then
          UPW(kts-1,:) = UPW(kts-1,:)*factor
        end if
  
        !
        ! Make sure that the thv and qt fluxes are not more than
        ! their values computed from the surface scheme.
        !
  
        ! Total updraft flux from lowest level
        QTsrfF=0.
        THVsrfF=0.
        DO I=1,NUP2
          QTsrfF  = QTsrfF +UPW(kts-1,I)*UPA(kts-1,I)*(UPQT(kts-1,I)-QT(kts))
          THVsrfF = THVsrfF+UPW(kts-1,I)*UPA(kts-1,I)*(UPTHV(kts-1,I)-THV(kts))
        ENDDO
  
        ! Adjust updraft THV so updraft flux is <90% of surface flux
        if (THVsrfF .gt. 0.9*wthv .and. THVsrfF .gt. 0.1) then
          UPTHV(kts-1,:)=(UPTHV(kts-1,:)-THV(kts))*0.9*wthv/THVsrfF+THV(kts)
        endif
  
        ! Adjust updraft QT so updraft flux is <90% of surface flux
        IF ( (QTsrfF .gt. 0.9*wqt) .and. (wqt .gt. 0.) )  then
          UPQT(kts-1,:)=(UPQT(kts-1,:)-QT(kts))*0.9*wqt/QTsrfF+QT(kts)
        ENDIF
  
        ! Compute condensation and initial updraft THL, QL, QI
        DO I=1,NUP2
          call condensation_edmfA(UPTHV(kts-1,I),UPQT(kts-1,I),P(kts-1), EXFH(kts-1), &
                                  UPTHL(kts-1,I),UPQL(kts-1,I),UPQI(kts-1,I), &
                                  mfparams%ice_ramp)
        ENDDO
  
        !========================
        !  Integrate updrafts
        !========================
  
        DO I=1,NUP2  ! updraft loop
          vertint: DO k=KTS,KTE ! vertical loop
  
            if (UPW(K-1,I).le.0.) exit vertint
  
              EntExp  = exp(-ENT(K,I)*(ZW(k)-ZW(k-1)))
              EntExpU = exp(-ENT(K,I)*(ZW(k)-ZW(k-1))*MFPARAMS%EntUFac)
  
              ! Effect of mixing on thermodynamic variables in updraft
              QTn  = QT(K)*(1-EntExp)+UPQT(K-1,I)*EntExp
              THLn = THL(K)*(1-EntExp)+UPTHL(K-1,I)*EntExp
              Un   = U(K)*(1-EntExpU)+UPU(K-1,I)*EntExpU
              Vn   = V(K)*(1-EntExpU)+UPV(K-1,I)*EntExpU
  
              ! Calculate condensation
              call condensation_edmf(QTn,THLn,P(K),EXFH(K),THVn,QCn,wf,mfparams%ice_ramp)
  
              ! Calculate and remove precipitation
              if (MFPARAMS%PRCPCRIT.gt.0.) then
                QP = max(0.,QCn-MFPARAMS%PRCPCRIT)
                QCn = QCn - QP
                QTn = QTn - QP
                THLn = THLn + (MAPL_ALHL*wf+(1.-wf)*MAPL_ALHS)/mapl_cp*QP/EXFH(k)
                QR(K) = QR(K) + UPA(K-1,I)*QP*wf
                QS(K) = QS(K) + UPA(K-1,I)*QP*(1.-wf)
              end if
  
              ! vertical velocity
              B=mapl_grav*(0.5*(THVn+UPTHV(k-1,I))/THV(k)-1.)
              ! represent deceleration from dynamic pressure approaching inversion
              WP=MFPARAMS%WB*ENT(K,I)+MFPARAMS%WC*wcfac(k)
              IF (WP==0.) THEN
                Wn2=UPW(K-1,I)**2+2.*MFPARAMS%WA*B*(ZW(k)-ZW(k-1))
              ELSE
                EntW=exp(-2.*WP*(ZW(k)-ZW(k-1)))
                Wn2=EntW*UPW(k-1,I)**2+(1.-EntW)*MFPARAMS%WA*B/WP
              END IF
  
  
              IF (Wn2>0.) THEN
                 UPW(K,I)=sqrt(Wn2)
                 UPTHV(K,I)=THVn
                 UPTHL(K,I)=THLn
                 UPQT(K,I)=QTn
                 UPQL(K,I)=QCn*wf
                 UPQI(K,I)=QCn*(1.-wf)
                 UPU(K,I)=Un
                 UPV(K,I)=Vn
                 UPA(K,I)=UPA(K-1,I)
              ELSE
                UPW(K,I) = 0.
                UPA(K,I) = 0.
                exit vertint
              END IF

          ENDDO vertint  ! loop over vertical
        ENDDO ! I: loop over updrafts
  
        if (associated(entx)) then
          do k=kts,kte
            tmp = sum(UPA(k,:))
            if (tmp .gt. 0.) then   ! weighted avg of lateral entrainment rate
              entx(IH,JH,KTE-k+KTS) = sum(UPA(k,:)*ENT(k,:))/tmp
            else
              entx(IH,JH,KTE-k+KTS) = MAPL_UNDEF
            end if
          end do
        end if
  
  
    ! CFL condition: Check that mass flux does not exceed layer mass at any level
    ! If it does, rescale updraft area.
    ! See discussion in Beljaars et al 2018 [ECMWF Tech Memo]
  
        factor = 1.0
        DO k=KTS,KTE
          mf = SUM(RHOE(K)*UPA(K,:)*UPW(K,:))
          if (mf .gt. MFPARAMS%MFLIMFAC*dp(K)/(MAPL_GRAV*dt)) then
             factor = min(factor,MFPARAMS%MFLIMFAC*dp(K)/(mf*MAPL_GRAV*dt) )
          end if
        ENDDO
        UPA = factor*UPA
        QR  = factor*QR
        QS  = factor*QS
  
    ! Rescale UPA if MF TKE more than half of prognostic TKE near surface
    ! Prevents instability due to MF without KH
        K = KTS
        tmp = 0.
        tmp2 = 0.
        factor = 1.
        DO WHILE (ZW(K).lt.200. .and. SUM(UPW(K-1,:)).gt.0.)
           tmp = 0.25*SUM(UPA(K,:)*UPW(K,:)*UPW(K,:)+UPA(K-1,:)*UPW(K-1,:)*UPW(K-1,:))
           tmp2 = TKE3(IH,JH,KTE-K+KTS)
           factor = min(factor,0.5*tmp2/tmp)
           K = K+1
        END DO
        if (factor.lt.1.) then
          UPA    = factor*UPA
          QR     = factor*QR
          QS     = factor*QS
        end if
  
        DO k=KTS,KTE
          edmfmf(IH,JH,KTE-k+KTS-1) = rhoe(K)*SUM(upa(K,:)*upw(K,:))
        ENDDO
        DQRDT(IH,JH,KTS:KTE) = QR(KTE:KTS:-1)/DT
        DQSDT(IH,JH,KTS:KTE) = QS(KTE:KTS:-1)/DT
  
  
        !
        ! writing updraft properties for output
        ! all variables, except Areas are now multipled by the area
        !
        dry_a     = 0.
        moist_a   = 0.

        DO I=1,NUP2 ! first sum over all i-updrafts
          DO k=KTS-1,KTE  ! loop in vertical
            IF ((UPQL(K,I)>0.) .OR. UPQI(K,I)>0.)  THEN
              moist_a(K) = moist_a(K)+UPA(K,I)
            ELSE
              dry_a(K) = dry_a(K)+UPA(K,I)
            ENDIF
          ENDDO
        END DO
  
        if (calc_avg_diag) then
          dry_w     = 0.
          moist_w   = 0.
          dry_qt    = 0.
          moist_qt  = 0.
          dry_thl   = 0.
          moist_thl = 0.
          dry_u     = 0.
          moist_u   = 0.
          dry_v     = 0.
          moist_v   = 0.
          moist_qc  = 0.

          DO I=1,NUP2   ! sum over all i-updrafts
            DO k=KTS-1,KTE  ! loop over vertical
              IF ((UPQL(K,I)>0.) .OR. UPQI(K,I)>0.)  THEN
                moist_w(K)   = moist_w(K)   + UPA(K,I)*UPW(K,I)
                moist_qt(K)  = moist_qt(K)  + UPA(K,I)*UPQT(K,I)
                moist_thl(K) = moist_thl(K) + UPA(K,I)*UPTHL(K,I)
                moist_u(K)   = moist_u(K)   + UPA(K,I)*UPU(K,I)
                moist_v(K)   = moist_v(K)   + UPA(K,I)*UPV(K,I)
                moist_qc(K)  = moist_qc(K)  + UPA(K,I)*(UPQL(K,I)+UPQI(K,I))
              ELSE
                dry_w(K)   = dry_w(K)   + UPA(K,I)*UPW(K,I)
                dry_qt(K)  = dry_qt(K)  + UPA(K,I)*UPQT(K,I)
                dry_thl(K) = dry_thl(K) + UPA(K,I)*UPTHL(K,I)
                dry_u(K)   = dry_u(K)   + UPA(K,I)*UPU(K,I)
                dry_v(K)   = dry_v(K)   + UPA(K,I)*UPV(K,I)
              ENDIF
            ENDDO ! vertical loop
          END DO ! updraft loop
         
          DO k = KTS-1,KTE
            IF (dry_a(k)>0.) THEN  ! divide by area for average
              dry_w(k)   = dry_w(k)  /dry_a(k)
              dry_qt(k)  = dry_qt(k) /dry_a(k)
              dry_thl(k) = dry_thl(k)/dry_a(k)
              dry_u(k)   = dry_u(k)  /dry_a(k)
              dry_v(k)   = dry_v(k)  /dry_a(k)
            ELSE
              dry_w(k)   = mapl_undef
              dry_qt(k)  = mapl_undef
              dry_thl(k) = mapl_undef
              dry_u(k)   = mapl_undef
              dry_v(k)   = mapl_undef
            ENDIF
            IF (moist_a(k)>0.) THEN
              moist_w(k)   = moist_w(k)  / moist_a(k)
              moist_qt(k)  = moist_qt(k) / moist_a(k)
              moist_thl(k) = moist_thl(k)/ moist_a(k)
              moist_u(k)   = moist_u(k)  / moist_a(k)
              moist_v(k)   = moist_v(k)  / moist_a(k)
              moist_qc(k)  = moist_qc(k) / moist_a(k)
            ELSE
              moist_w(k)   = mapl_undef
              moist_qt(k)  = mapl_undef
              moist_thl(k) = mapl_undef
              moist_u(k)   = mapl_undef
              moist_v(k)   = mapl_undef
              moist_qc(k)  = mapl_undef
            ENDIF
          ENDDO     ! loop in vertical
        end if
  
    !
    ! computing variables needed for solver
    !
        s_aw   = 0.
        s_aws  = 0.
        s_awqv = 0.
        s_awql = 0.
        s_awqi = 0.
        s_awu  = 0.
        s_awv  = 0.
  
        s_buoyf = 0.
        s_aqt2  = 0.
        s_awqt  = 0.
        s_aw2   = 0.
        s_aw3   = 0.
        s_aqt3  = 0.
        s_ahl3  = 0.
        s_ahl2  = 0.
        s_awhl  = 0.
        s_ahlqt = 0.
  
        qte = (QTI(:)-SUM(UPA(:,:)*UPQT(:,:),DIM=2))/(1.-SUM(UPA(:,:),DIM=2))
        s_aqt3(:) = (1.-SUM(UPA,DIM=2))*(QTE-QTI)**3

        if (MFPARAMS%implicit == 1) then
           do i=1,nup2
              do k=kts-1,kte
                 s_aw(K)=s_aw(K)+UPA(K,I)*UPW(K,I)
                 s_aw2(K)=s_aw2(K)+UPA(K,I)*UPW(K,I)*UPW(K,I)
                 s_aw3(K)=s_aw3(K)+UPA(K,I)*UPW(K,I)*UPW(K,I)*UPW(K,I)
                 s_aqt2(K)=s_aqt2(K)+UPA(K,I)*(UPQT(K,I)-QTI(K))*(UPQT(K,I)-QTI(K))
                 s_aqt3(K)=s_aqt3(K)+UPA(K,I)*(UPQT(K,I)-QTI(K))**3
                 s_ahlqt(K)=s_ahlqt(K)+exfh(k)*UPA(K,I)*(UPQT(K,I)-QTI(K))*(UPTHL(K,i)-THLI(K))
                 tmp = mapl_cp*exfh(k)*UPTHL(K,i) + mapl_grav*zw(k) + phis(IH,JH) &
                     + mapl_alhl*UPQL(K,i) + UPQI(K,I)*mapl_alhs
                 ltm=exfh(k)*(UPTHL(K,i)-THLI(K))
                 s_aws(k)  = s_aws(K)+UPA(K,i)*UPW(K,i)*tmp 
                 s_ahl2(k) = s_ahl2(K)+UPA(K,i)*ltm*ltm
                 s_ahl3(k) = s_ahl3(K)+UPA(K,i)*ltm*ltm*ltm
                 s_awhl(k) = s_awhl(K)+UPA(K,i)*UPW(K,I)*ltm
                 s_awu(k)  = s_awu(K)  + UPA(K,i)*UPW(K,I)*UPU(K,I)
                 s_awv(k)  = s_awv(K)  + UPA(K,i)*UPW(K,I)*UPV(K,I)
                 s_awqv(k) = s_awqv(K) + UPA(K,i)*UPW(K,I)*(UPQT(K,I) - UPQI(K,I) - UPQL(K,I))
                 s_awql(k) = s_awql(K) + UPA(K,i)*UPW(K,I)*UPQL(K,I)
                 s_awqi(k) = s_awqi(K) + UPA(K,i)*UPW(K,I)*UPQI(K,I)
                 s_awqt(k)  = s_awqt(K)  + UPA(K,i)*UPW(K,I)*(UPQT(K,I) - QTI(K))
                 mftke(IH,JH,k) = mftke(IH,JH,k) + UPA(KTE+KTS-K-1,i)*0.5*UPW(KTE+KTS-K-1,I)*UPW(KTE+KTS-K-1,I)
              end do
              do k=kts,kte
                 mfthvt=0.5*(UPA(k-1,I)*UPW(k-1,I)*UPTHV(k-1,I)+UPA(k,I)*UPW(k,I)*UPTHV(k,I))
                 mft=0.5*(UPA(k-1,I)*UPW(k-1,I)+UPA(k,I)*UPW(k,I))
                 s_buoyf(k)=s_buoyf(k)+(mfthvt-mft*THV(k))*exf(k)
              end do
           end do
        else  ! if not implicit approach
           do i=1,nup2
              do k=kts-1,kte
                 s_aw(K)=s_aw(K)+UPA(K,I)*UPW(K,I)
                 s_aw2(K)=s_aw2(K)+UPA(K,I)*UPW(K,I)*UPW(K,I)
                 s_aw3(K)=s_aw3(K)+UPA(K,I)*UPW(K,I)*UPW(K,I)*UPW(K,I)
                 s_aqt2(K)=s_aqt2(K)+UPA(K,I)*(UPQT(K,I)-QTI(K))*(UPQT(K,I)-QTI(K))
                 s_aqt3(K)=s_aqt3(K)+UPA(K,I)*(UPQT(K,I)-QTI(K))**3
                 s_ahlqt(K)=s_ahlqt(K)+exfh(k)*UPA(K,I)*(UPQT(K,I)-QTI(K))*(UPTHL(K,i)-THLI(K))
                 tmp = mapl_cp*exfh(k)*( UPTHL(K,i) - THLI(K) ) &
                     + mapl_alhl*( UPQL(K,i) - QLI(K) )   &
                     + mapl_alhs*( UPQI(K,I) - QII(K) )
                 ltm=exfh(k)*(UPTHL(K,i)-THLI(K))
                 s_aws(k)  = s_aws(K)+UPA(K,i)*UPW(K,i)*tmp 
                 s_ahl2(k) = s_ahl2(K)+UPA(K,i)*ltm*ltm
                 s_ahl3(k) = s_ahl3(K)+UPA(K,i)*ltm*ltm*ltm
                 s_awhl(k) = s_awhl(K)+UPA(K,i)*UPW(K,I)*ltm
                 s_awu(k)  = s_awu(K)  + UPA(K,i)*UPW(K,I)*(UPU(K,I) - UI(K))
                 s_awv(k)  = s_awv(K)  + UPA(K,i)*UPW(K,I)*(UPV(K,I) - VI(K))
                 s_awqv(k) = s_awqv(K) + UPA(K,i)*UPW(K,I)*(UPQT(K,I) - UPQI(K,I) - UPQL(K,I) - QVI(K))
                 s_awql(k) = s_awql(K) + UPA(K,i)*UPW(K,I)*(UPQL(K,I) - QLI(K))
                 s_awqi(k) = s_awqi(K) + UPA(K,i)*UPW(K,I)*(UPQI(K,I) - QII(K))
                 s_awqt(k)  = s_awqt(K)  + UPA(K,i)*UPW(K,I)*(UPQT(K,I) - QTI(K))
                 mftke(IH,JH,k) = mftke(IH,JH,k) + UPA(KTE+KTS-K-1,i)*0.5*UPW(KTE+KTS-K-1,I)*UPW(KTE+KTS-K-1,I)
              end do
              do k=kts,kte
                 mfthvt=0.5*(UPA(k-1,I)*UPW(k-1,I)*UPTHV(k-1,I)+UPA(k,I)*UPW(k,I)*UPTHV(k,I))
                 mft=0.5*(UPA(k-1,I)*UPW(k-1,I)+UPA(k,I)*UPW(k,I))
                 s_buoyf(k)=s_buoyf(k)+(mfthvt-mft*THV(k))*exf(k)
              end do
           end do
        end if
  
        !
        ! turn around the outputs and fill them in the 3d fields
        !
        dry_a3(IH,JH,KTS-1:KTE)    = dry_a(KTE:KTS-1:-1)
        moist_a3(IH,JH,KTS-1:KTE)  = moist_a(KTE:KTS-1:-1)
        if (associated(dry_w3))     dry_w3(IH,JH,KTS-1:KTE)     = dry_w(KTE:KTS-1:-1)
        if (associated(moist_w3))   moist_w3(IH,JH,KTS-1:KTE)   = moist_w(KTE:KTS-1:-1)
        if (associated(dry_qt3))    dry_qt3(IH,JH,KTS-1:KTE)    = dry_qt(KTE:KTS-1:-1)
        if (associated(moist_qt3))  moist_qt3(IH,JH,KTS-1:KTE)  = moist_qt(KTE:KTS-1:-1)
        if (associated(dry_thl3))   dry_thl3(IH,JH,KTS-1:KTE)   = dry_thl(KTE:KTS-1:-1)
        if (associated(moist_thl3)) moist_thl3(IH,JH,KTS-1:KTE) = moist_thl(KTE:KTS-1:-1)
        if (associated(dry_u3))     dry_u3(IH,JH,KTS-1:KTE)     = dry_u(KTE:KTS-1:-1)
        if (associated(moist_u3))   moist_u3(IH,JH,KTS-1:KTE)   = moist_u(KTE:KTS-1:-1)
        if (associated(dry_v3))     dry_v3(IH,JH,KTS-1:KTE)     = dry_v(KTE:KTS-1:-1)
        if (associated(moist_v3))   moist_v3(IH,JH,KTS-1:KTE)   = moist_v(KTE:KTS-1:-1)
        if (associated(moist_qc3))  moist_qc3(IH,JH,KTS-1:KTE)  = moist_qc(KTE:KTS-1:-1)
  
  
        ! Note values were initialized to zero above.
        ! Ending loop at interface above lowest layer.
        DO K=KTS-1,KTE-1
          ! outputs - variables needed for solver
          aw3(IH,JH,K)   = s_aw(KTE+KTS-K-1)
          aws3(IH,JH,K)  = s_aws(KTE+KTS-K-1)
          awqv3(IH,JH,K) = s_awqv(KTE+KTS-K-1)
          awql3(IH,JH,K) = s_awql(KTE+KTS-K-1)
          awqi3(IH,JH,K) = s_awqi(KTE+KTS-K-1)
          awu3(IH,JH,K)  = s_awu(KTE+KTS-K-1)
          awv3(IH,JH,K)  = s_awv(KTE+KTS-K-1)
          ae3(IH,JH,K)   = (1.-dry_a(KTE+KTS-K-1)-moist_a(KTE+KTS-K-1))
          mfwhl(IH,JH,K) = s_awhl(KTE+KTS-K-1)
          mfwqt(IH,JH,K) = s_awqt(KTE+KTS-K-1)
        ENDDO
  
        s_awqv(KTS-1) = 0.
        s_awql(KTS-1) = 0.
        s_awqi(KTS-1) = 0.
  
  
        ! buoyancy is defined on full levels
        DO k=kts,kte
          buoyf(IH,JH,K)  = s_buoyf(KTE+KTS-K)    ! can be used in SHOC
          mfw2(IH,JH,K)   = 0.5*(s_aw2(KTE+KTS-K-1)+s_aw2(KTE+KTS-K))
          mfw3(IH,JH,K)   = 0.5*(s_aw3(KTE+KTS-K-1)+s_aw3(KTE+KTS-K))
          mfhlqt(IH,JH,K) = 0.5*(s_ahlqt(KTE+KTS-K-1)+s_ahlqt(KTE+KTS-K))
          if (SUM(moist_a(KTS-1:KTE+KTS-K)).le.1e-4) then
            mfqt3(IH,JH,K) = 0.
            mfhl3(IH,JH,K) = 0.
          else
            mfqt3(IH,JH,K)  = 0.5*(s_aqt3(KTE+KTS-K-1)+s_aqt3(KTE+KTS-K))
            mfhl3(IH,JH,K)  = 0.5*(s_ahl3(KTE+KTS-K-1)+s_ahl3(KTE+KTS-K))
          end if
        ENDDO
  
        where (UPA.eq.0.)
          UPW   = MAPL_UNDEF
          UPTHL = MAPL_UNDEF
          UPQT  = MAPL_UNDEF
        end where
        if (associated(EDMF_PLUMES_W))   EDMF_PLUMES_W(IH,JH,KTS-1:KTE,:)   = upw(KTE:KTS-1:-1,:)
        if (associated(EDMF_PLUMES_THL)) EDMF_PLUMES_THL(IH,JH,KTS-1:KTE,:) = upthl(KTE:KTS-1:-1,:)
        if (associated(EDMF_PLUMES_QT))  EDMF_PLUMES_QT(IH,JH,KTS-1:KTE,:)  = upqt(KTE:KTS-1:-1,:)
  

        ! OPTIMIZATION: Moved Trisolver 3D updates INSIDE the column loop!
        ! Completely avoids reading back all the memory.
        
        YS(IH,JH,kte)  = -rhoe3(IH,JH,kte-1) * aws3(IH,JH,kte-1)
        YQV(IH,JH,kte) = -rhoe3(IH,JH,kte-1) * awqv3(IH,JH,kte-1)
        YQL(IH,JH,kte) = -rhoe3(IH,JH,kte-1) * awql3(IH,JH,kte-1)
        YQI(IH,JH,kte) = -rhoe3(IH,JH,kte-1) * awqi3(IH,JH,kte-1)
        YU(IH,JH,kte)  = -rhoe3(IH,JH,kte-1) * awu3(IH,JH,kte-1)
        YV(IH,JH,kte)  = -rhoe3(IH,JH,kte-1) * awv3(IH,JH,kte-1)
      
        DO k = kts, kte-1
           YS(IH,JH,k)  = rhoe3(IH,JH,k)*aws3(IH,JH,k)   - rhoe3(IH,JH,k-1)*aws3(IH,JH,k-1)
           YQV(IH,JH,k) = rhoe3(IH,JH,k)*awqv3(IH,JH,k)  - rhoe3(IH,JH,k-1)*awqv3(IH,JH,k-1)
           YQL(IH,JH,k) = rhoe3(IH,JH,k)*awql3(IH,JH,k)  - rhoe3(IH,JH,k-1)*awql3(IH,JH,k-1)
           YQI(IH,JH,k) = rhoe3(IH,JH,k)*awqi3(IH,JH,k)  - rhoe3(IH,JH,k-1)*awqi3(IH,JH,k-1)
           YU(IH,JH,k)  = rhoe3(IH,JH,k)*awu3(IH,JH,k)   - rhoe3(IH,JH,k-1)*awu3(IH,JH,k-1)
           YV(IH,JH,k)  = rhoe3(IH,JH,k)*awv3(IH,JH,k)   - rhoe3(IH,JH,k-1)*awv3(IH,JH,k-1)
        END DO
      
        ! Deal with implied condensation
        DO k = kts, kte
           if (YQI(IH,JH,k) < 0. .and. YQL(IH,JH,k) > 0.) then
              tmp = min(YQL(IH,JH,k), -YQI(IH,JH,k))
              YQL(IH,JH,k) = YQL(IH,JH,k) - tmp
              YQI(IH,JH,k) = YQI(IH,JH,k) + tmp
              YS(IH,JH,k)  = YS(IH,JH,k) + tmp*MAPL_ALHF
           end if
           if (YQI(IH,JH,k) < 0.) then
              YQV(IH,JH,k) = YQV(IH,JH,k) + YQI(IH,JH,k)
              YS(IH,JH,k)  = YS(IH,JH,k) - YQI(IH,JH,k)*MAPL_ALHS
              YQI(IH,JH,k) = 0.
           end if
           if (YQL(IH,JH,k) < 0.) then
              YQV(IH,JH,k) = YQV(IH,JH,k) + YQL(IH,JH,k)
              YS(IH,JH,k)  = YS(IH,JH,k) - YQL(IH,JH,k)*MAPL_ALHL
              YQL(IH,JH,k) = 0.
           end if
           
           tmp = mapl_grav * dt / ( pw3(IH,JH,k)-pw3(IH,JH,k-1) )
           YS(IH,JH,k)  = tmp * YS(IH,JH,k)
           YQV(IH,JH,k) = tmp * YQV(IH,JH,k)
           YQL(IH,JH,k) = tmp * YQL(IH,JH,k)
           YQI(IH,JH,k) = tmp * YQI(IH,JH,k)
           YU(IH,JH,k)  = tmp * YU(IH,JH,k)
           YV(IH,JH,k)  = tmp * YV(IH,JH,k)
        END DO
  
        ! Detrained mass flux calculation within the column loop
        if (associated(edmf_dmf)) then
           DO k = kts, kte
              edmf_dmf(IH,JH,k) = max(0., edmfmf(IH,JH,k-1) - edmfmf(IH,JH,k))
              edmf_dmf(IH,JH,k) = edmf_dmf(IH,JH,k) + edmfmf(IH,JH,k) * &
                                  (1. - exp( -entx(IH,JH,k)*(zw3(IH,JH,k-1)-zw3(IH,JH,k)) ))
              if (moist_a3(IH,JH,k) <= 0.) edmf_dmf(IH,JH,k) = 0.
           END DO
        end if

       end if  !  IF ( mfdepth>100m)
      END IF   !  IF ( wthv > 0.0 )
  
      ENDDO ! IH loop (INNER LOOP)
    ENDDO ! JH loop (OUTER LOOP)
    
  END SUBROUTINE run_edmf
  
  
  subroutine calc_mf_depth(kts,kte,t,z,q,p,ztop,wthv,wqt)
  
    integer, intent(in   )                     :: kts, kte
    real,    intent(in   ), dimension(kts:kte) :: t, z, q, p
    real,    intent(in   )                     :: wthv, wqt
    real,    intent(  out)                     :: ztop
  
    real     :: tep,z1,z2,t1,t2,qp,pp,qsp,dqp,dqsp,wstar,qstar,thstar,sigmaQT,sigmaTH
    integer  :: k
  
     wstar=max(0.1,(mapl_grav*wthv*1e3/t(kte))**(1./3.))  ! convective velocity scale
     qstar=max(0.,wqt)/wstar
     thstar=max(0.,wthv)/wstar
  
     sigmaQT=2.0*qstar
     sigmaTH=2.0*thstar
  
    tep  = t(kte)+max(0.1,sigmaTH) ! parcel values
    qp   = q(kte)+sigmaQT
  
    t1   = t(kte)
    z1   = z(kte)
    ztop = z(kte)
  
    do k = kte-1 , kts+1, -1
      z2 = z(k)
      t2 = t(k)
      pp = p(k)
  
      tep   = tep - MAPL_GRAV*( z2-z1 )/MAPL_CP
  
      qp    = qp  + (0.7/1000.)*(z2-z1)*(q(k)-qp)   ! assume fractional entrainment rate of 0.75/km
      tep   = tep + (0.7/1000.)*(z2-z1)*(t(k)-tep)
  
      dqsp  = GEOS_DQSAT(tep , pp , qsat=qsp,  pascals=.true. )
  
      dqp   = max( qp - qsp, 0. )/(1.+(MAPL_ALHL/MAPL_CP)*dqsp )
      qp    = qp - dqp
      tep   = tep  + MAPL_ALHL * dqp/MAPL_CP
  
  
      ! compare Tv env vs parcel
      if ( t2*(1.+MAPL_VIREPS*q(k)) .ge. tep*(1.+MAPL_VIREPS*qp)+0.2 ) then
        ztop = 0.5*(z2+z1)
        exit
      end if
  
      z1 = z2
      t1 = t2
    enddo  ! k loop
  
    return
  
  end subroutine calc_mf_depth
  
  
  subroutine condensation_edmf(QT,THL,P,EXN,THV,QC,wf,ice_ramp)
  !
  ! Newton-Raphson condensation solver
  ! Calculates THV and QC
  !
  use GEOS_UtilsMod, only : GEOS_Qsat
  
  real,intent(in)  :: QT,THL,P,EXN,ice_ramp
  real,intent(out) :: THV,QC,wf
  
  
  integer :: niter,i
  real    :: diff,t,qs,qcold
  real    :: Rv, L, dQS_dT, dT_dQC, f_prime, QC_new
  
  ! mapl_vireps = Rv/Rd - 1 => Rv = Rd*(1+mapl_vireps)
  ! mapl_kappa = Rd/Cp => Rd = mapl_kappa * mapl_cp
  Rv = mapl_kappa * mapl_cp * (1.0 + mapl_vireps)
  
  niter=5
  diff=1.e-5
  QC=0.
  
  do i=1,NITER
    L = get_alhl(EXN*THL, ice_ramp) ! First guess for L
    T = EXN*THL + L/mapl_cp*QC
    
    L = get_alhl(T, ice_ramp)       ! Re-evaluate L at current T
    QS = geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
    
    dQS_dT = QS * L / (Rv * max(T*T, 1.0))
    dT_dQC = L / mapl_cp
    f_prime = 1.0 + dQS_dT * dT_dQC
    
    QCOLD=QC
    QC_new = QC - (QC - (QT - QS)) / f_prime
    QC = max(0.0, QC_new)
    
    if (abs(QC-QCOLD)<Diff) exit
  enddo
  
  T = EXN*THL + get_alhl(T,ice_ramp)/mapl_cp*QC
  QS = geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
  QC = max(QT-QS,0.)
  THV = (THL+get_alhl(T,ice_ramp)/mapl_cp*QC/EXN)*(1.+MAPL_VIREPS*(QT-QC)-QC)
  wf = water_f(T,ice_ramp)
  
  end subroutine condensation_edmf
  
  
  subroutine condensation_edmfA(THV,QT,P,EXN,THL,QL,QI,ice_ramp)
  !
  ! Newton-Raphson condensation solver 
  ! Calculates QL,QI from THV and QT
  !
  
  use GEOS_UtilsMod, only : GEOS_Qsat
  
  real,intent(in)  :: THV,QT,P,EXN,ice_ramp
  real,intent(out) :: THL,QL,QI
  
  
  integer :: niter,i
  real    :: diff,t,qs,qcold,wf,qc
  real    :: Rv, L, dQS_dT, dT_dQC, f_prime, QC_new, Denom
  
  Rv = mapl_kappa * mapl_cp * (1.0 + mapl_vireps)
  
  niter=5
  diff=1.e-5
  QC=0.
  
  do i=1,NITER
     Denom = 1.0 + MAPL_VIREPS*(QT-QC) - QC
     T = EXN*THV/Denom
     QS = geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
     
     L = get_alhl(T,ice_ramp)
     dQS_dT = QS * L / (Rv * max(T*T, 1.0))
     dT_dQC = T * (1.0 + MAPL_VIREPS) / Denom
     f_prime = 1.0 + dQS_dT * dT_dQC
     
     QCOLD=QC
     QC_new = QC - (QC - (QT - QS)) / f_prime
     QC = max(0.0, QC_new)
     
     if (abs(QC-QCOLD)<Diff) exit
  enddo
  
   THL=(T-QC*get_alhl(T,ice_ramp)/mapl_cp)/EXN
   wf=water_f(T,ice_ramp)
   QL=QC*wf
   QI=QC*(1.-wf)
  
  end subroutine condensation_edmfA
  
  
  function  get_alhl3(T,IM,JM,LM,iceramp)
  
  real,dimension(IM,JM,LM) ::  T,get_alhl3
  real :: iceramp
  integer :: IM,JM,LM
  integer :: IH,JH,L
  
  do jh=1,jm
    do ih=1,im
      do l=1,lm
          get_alhl3(IH,JH,l)=get_alhl(T(IH,JH,l),iceramp)
      enddo
    enddo
  enddo
  
  end function get_alhl3
  
  function get_alhl(T,iceramp)
     real :: T,get_alhl,iceramp,wf
      wf=water_f(T,iceramp)
      get_alhl=wf*mapl_alhl+(1.-wf)*mapl_alhs
  end function get_alhl
  
  
  ! OPTIMIZATION: Removed if/else branches for pure math compilation
  function water_f(T,iceramp)
    real :: T, iceramp, water_f, Tw, Tmin, Tmax
    Tmax = 0.
    Tmin = -abs(iceramp)
    Tw = T - mapl_celsius_to_kelvin
    water_f = MIN(MAX((Tw - Tmin) / (Tmax - Tmin), 0.0), 1.0)
  end function water_f
  
  
  subroutine Poisson(kstart,kend,istart,iend,mu,POI,seed)
  
  integer, intent(in) :: istart,iend,kstart,kend
  real,dimension(kstart:kend,istart:iend),intent(in) :: MU
  integer, dimension(kstart:kend,istart:iend), intent(out) :: POI
  integer,dimension(2),  intent(in) :: seed
  
  integer :: i,k
  integer(8) :: rng_state
  
  ! Initialize deterministic state from the passed seed
  rng_state = int(seed(1), 8) * 2147483647_8 + int(seed(2), 8)
  if (rng_state == 0_8) rng_state = 123456789_8
  
  do i=istart,iend
    do k=kstart,kend
      poi(k,i)=poidev(mu(k,i), rng_state)
    enddo
  enddo
  
  end subroutine Poisson
  
  
        ! OPTIMIZATION: Thread-safety achieved by removing SAVE statements
        FUNCTION poidev(xm, rng_state)
        REAL :: poidev,xm
        INTEGER(8), INTENT(INOUT) :: rng_state
        
        REAL :: alxm,em,g,sq,t,y
        
        if (xm.lt.12.)then
          g=exp(-xm)
          em=-1
          t=1.
  2       em=em+1.
          t=t*ran1_lcg(rng_state)
          if (t.gt.g) goto 2
        else
          sq=sqrt(2.*xm)
          alxm=log(xm)
          g=xm*alxm-gammln(xm+1.)
  1       y=tan(MAPL_PI*ran1_lcg(rng_state))
          em=sq*y+xm
          if (em.lt.0.) goto 1
          em=int(em)
          t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
          if (ran1_lcg(rng_state).gt.t) goto 1
        endif
        poidev=em
        return
        END FUNCTION poidev
  
        ! OPTIMIZATION: Thread-safety achieved by using PARAMETER arrays
        FUNCTION gammln(xx)
        REAL gammln,xx
        INTEGER j
        DOUBLE PRECISION ser,tmp,x,y
        DOUBLE PRECISION, PARAMETER :: stp_val = 2.5066282746310005d0
        DOUBLE PRECISION, PARAMETER :: cof(6) = [76.18009172947146d0, -86.50532032941677d0, &
                                                 24.01409824083091d0, -1.231739572450155d0,  &
                                                  0.1208650973866179d-2, -0.5395239384953d-5]

        x=xx
        y=x
        tmp=x+5.5d0
        tmp=(x+0.5d0)*log(tmp)-tmp
        ser=1.000000000190015d0
        do 11 j=1,6
          y=y+1.d0
          ser=ser+cof(j)/y
  11    continue
        gammln=tmp+log(stp_val*ser/x)
        return
        END FUNCTION gammln
  
        ! A fast inline 32-bit Linear Congruential Generator (LCG)
        FUNCTION ran1_lcg(state)
          INTEGER(8), INTENT(INOUT) :: state
          REAL :: ran1_lcg
          state = mod(state * 1103515245_8 + 12345_8, 2147483648_8)
          ran1_lcg = real(state) / 2147483648.0
        END FUNCTION ran1_lcg
  
   end module edmf_mod
      
