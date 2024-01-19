module edmf_mod

!
! Mass flux updraft parameterization implemented by Kay Suselj (JPL).
! Reference: Suselj et al 2021, DOI: 10.1175/MWR-D-20-0183.1
! Additional development by Nathan Arnold and David New (GMAO).
!

use MAPL_ConstantsMod, only: mapl_epsilon, mapl_grav, mapl_cp,  &
                             mapl_alhl, mapl_p00, mapl_vireps,  &
                             mapl_alhs, mapl_kappa, mapl_rgas

use MAPL_Mod,          only: mapl_undef

use GEOS_Mod

implicit none

real, parameter ::     &
     WSTARmin = 1.e-3, &
     zpblmin  = 100.,  &
     onethird = 1./3., &
     r        = 2.

 type EDMFPARAMS_TYPE
    integer :: DISCRETE
    integer :: IMPLICIT
    integer :: ENTRAIN
    integer :: DOCLASP
    integer :: NUP
    integer :: THERMAL_PLUME
    integer :: TEST
    integer :: DEBUG
    integer :: ET
    real    :: L0
    real    :: L0fac
    real    :: STOCHFRAC
    real    :: ENTWFAC
    real    :: EDFAC
    real    :: ENT0
    real    :: ENT0LTS
    real    :: ALPHATH
    real    :: ALPHAQT
    real    :: ALPHAW
    real    :: PWMAX
    real    :: PWMIN
    real    :: WA
    real    :: WB
    real    :: AU0
    real    :: CTH1
    real    :: CTH2
    real    :: RH0_QB
    real    :: C_KH_MF
    real    :: MFLIMFAC
    real    :: ICE_RAMP
 endtype EDMFPARAMS_TYPE
 type (EDMFPARAMS_TYPE) :: MFPARAMS

public run_edmf, mfparams

contains

SUBROUTINE RUN_EDMF(its,ite, jts,jte, kts,kte, dt, &   ! Inputs
                    phis,                          &
                    zlo3,                          &
                    zw3,                           &
                    pw3,                           &
                    rhoe3,                         &
                    tke3,                          &
                    u3,                            &
                    v3,                            &
                    t3,                            &
                    thl3,                          &
                    thv3,                          &
                    qt3,                           &
                    qv3,                           &
                    ql3,                           &
                    qi3,                           &
                    wthl2,                         &
                    wqt2,                          &
                    frland,                        &
                    pblh2,                         &
                    ! CLASP inputs for surface heterogeneity
     !              mfsrcthl, mfsrcqt, mfw, mfarea, &
                    ! outputs - variables needed for solver
                    ae3,                           &
                    aw3,                           &
                    aws3,                          &
                    awqv3,                         &
                    awql3,                         &
                    awqi3,                         &
                    awu3,                          &
                    awv3,                          &
                    ! Outputs required for SHOC and ADG PDF
                    mfw2,                          &
                    mfw3,                          &
                    mfqt3,                         &
                    mfhl3,                         &
                    mfwqt,                         &
                    mfhlqt,                        &
                    mfwhl,                         &
                    mftke,                         &
                    buoyf,                         &
                    edmfmf,                        &
                    dry_a3,                        &
                    moist_a3,                      &
                    ! Diagnostic outputs - updraft properties
                    dry_w3,                        &
                    moist_w3,                      &
                    dry_qt3,                       &
                    moist_qt3,                     &
                    dry_thl3,                      &
                    moist_thl3,                    &
                    dry_u3,                        &
                    moist_u3,                      &
                    dry_v3,                        &
                    moist_v3,                      &
                    moist_qc3,                     &
                    entx,                          &
                    mfdepth,                       &
                    edmf_plumes_w,                 &
                    edmf_plumes_thl,               &
                    edmf_plumes_qt )


   INTEGER, INTENT(IN) :: ITS,ITE,JTS,JTE,KTS,KTE
   REAL,    INTENT(IN) :: DT

   REAL,DIMENSION(ITS:ITE,JTS:JTE,KTS:KTE), INTENT(IN) :: U3,    &
                                                          V3,    &
                                                          T3,    &
                                                          THL3,  &
                                                          QT3,   &
                                                          THV3,  &
                                                          QV3,   &
                                                          QL3,   &
                                                          QI3,   &
                                                          ZLO3,  &
                                                          TKE3

   REAL,DIMENSION(ITS:ITE,JTS:JTE,KTS-1:KTE), INTENT(IN) :: ZW3, PW3, rhoe3

   REAL,DIMENSION(ITS:ITE,JTS:JTE), INTENT(IN) :: WTHL2,   &
                                                  WQT2,    &
                                                  PBLH2,   &
                                                  FRLAND,  &
                                                  PHIS

   ! CLASP inputs
   REAL,DIMENSION(ITS:ITE,JTS:JTE,KTS:KTE) :: mfsrcqt,mfsrcthl,mfw,mfarea

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

   REAL,DIMENSION(ITS:ITE,JTS:JTE,KTS:KTE), INTENT(OUT) :: buoyf,mfw2,mfw3,mfqt3,mfhl3,&!mfqt2,mfhl2,&
                                                        mfhlqt


  ! Diagnostic outputs
   REAL, DIMENSION(:,:), POINTER   :: mfdepth
   REAL, DIMENSION(:,:,:), POINTER :: dry_w3,   moist_w3,   &
                                      dry_qt3,  moist_qt3,  &
                                      dry_thl3, moist_thl3, &
                                      dry_u3,   moist_u3,   &
                                      dry_v3,   moist_v3,   &
                                      moist_qc3, entx

   REAL, DIMENSION(:,:,:,:), POINTER :: EDMF_PLUMES_W,   &
                                        EDMF_PLUMES_THL, &
                                        EDMF_PLUMES_QT


!============= Local variables =============

   ! updraft properties
   REAL,DIMENSION(KTS-1:KTE,1:MFPARAMS%NUP) :: UPW,UPTHL,UPQT,UPQL,UPQI, &
                                               UPA,UPU,UPV,UPTHV
   ! entrainment variables
   REAl,DIMENSION(KTS:KTE,1:MFPARAMS%NUP) :: ENT,ENTf
   INTEGER,DIMENSION(KTS:KTE,1:MFPARAMS%NUP) :: ENTi

   INTEGER :: K,I,IH,JH,NUP2
   REAL :: wthv,wstar,qstar,thstar,sigmaW,sigmaQT,sigmaTH,z0, &
           wmin,wmax,wlv,wtv,wp
   REAL :: B,QTn,THLn,THVn,QCn,Un,Vn,Wn2,EntEXP,EntEXPU,EntW,wf

! internal flipped variables (GEOS5)
   REAL,DIMENSION(KTS:KTE) :: U,V,THL,QT,THV,QV,QL,QI,ZLO
   REAL,DIMENSION(KTS-1:KTE)  :: ZW,P,THLI,QTI
   REAL,DIMENSION(KTS-1:KTE) :: UI, VI, QVI, QLI, QII

! internal surface cont
   REAL :: UST,WTHL,WQT,PBLH
   REAL,DIMENSION(KTS-1:KTE) :: dry_a, moist_a,dry_w,moist_w,          &
                                dry_qt,moist_qt,dry_thl,moist_thl,     &
                                dry_u,moist_u,dry_v,moist_v, moist_qc
   REAL,DIMENSION(KTS-1:KTE) :: s_aw,s_aws,s_awqv,s_awql,s_awqi,s_awu,s_awv
   REAL,DIMENSION(KTS:KTE)   :: s_buoyf
   REAL,DIMENSION(KTS-1:KTE) :: s_aw2,s_aw3,s_aqt3,s_ahl3,s_aqt2,  &
                                s_ahlqt,s_awqt,s_ahl2,s_awhl,qte
! exner function
   REAL,DIMENSION(KTS:KTE)   :: exf,dp,pmid
   REAL,DIMENSION(KTS-1:KTE) :: exfh,tmp1d
   REAL,DIMENSION(KTS-1:KTE) :: rhoe

   REAL :: L0,ztop,tmp,ltm,MFsrf,QTsrfF,THVsrfF,mft,mfthvt,mf,factor,lts
   INTEGER, DIMENSION(2) :: seedmf,the_seed

   LOGICAL :: calc_avg_diag

! velocity equation parameters
   REAL,PARAMETER :: Wa=1., &    ! buoyancy term
                     Wb=1.5      ! entrainment term
!                    Wa=1., &    ! original
!                    Wb=1.5

! min values to avoid singularities
   REAL,PARAMETER :: &
       WSTARmin=1.e-3, &
       PBLHmin=100.

   ! temporary, set
   mfsrcthl = -999.
   mfsrcqt  = -999.
   mfw      = -999.
   mfarea   = -999.

   tmp1d(:) = 1e-3

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
   !      mfqt2 =0.
   !      mfhl2 =0.

   if (associated(entx)) entx = mapl_undef

   ! this is the environmental area - by default 1.
   ae3=MFPARAMS%EDfac


   DO IH=ITS,ITE ! loop over the horizontal dimensions
    DO JH=JTS,JTE

      wthl=wthl2(IH,JH)/mapl_cp
      wqt=wqt2(IH,JH)
      !ust=ust2(IH,JH)
      pblh=pblh2(IH,JH)

      pblh=max(pblh,pblhmin)
      wthv=wthl+mapl_epsilon*thv3(IH,JH,kte)*wqt

      ! if CLASP enabled: mass flux is input
      ! if CLASP disabled: mass-flux if positive surface buoyancy flux and
      !                    TKE at 2nd model level above threshold
      IF ( (wthv > 0.0 .and. TKE3(IH,JH,kte-1)>0.01 .and. MFPARAMS%doclasp==0 .and. phis(IH,JH).lt.2e4)      &
      .or. (any(mfsrcthl(IH,JH,1:MFPARAMS%NUP) >= -2.0) .and. MFPARAMS%doclasp/=0)) then

!     print *,'wthv=',wthv,' wqt=',wqt,' wthl=',wthl,' edmfdepth=',edmfdepth(IH,JH)

      if (MFPARAMS%doclasp/=0) then
       nup2 = count(mfsrcthl(IH,JH,1:MFPARAMS%NUP)>=-2.0)
      else
       nup2 = MFPARAMS%NUP
      end if

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

      ! Estimate scale height for entrainment calculation
      if (mfparams%ET == 2 ) then
        pmid = 0.5*(pw3(IH,JH,kts-1:kte-1)+pw3(IH,JH,kts:kte))
        call calc_mf_depth(kts,kte,t3(IH,JH,:),zlo3(IH,JH,:)-zw3(IH,JH,kte),qv3(IH,JH,:),pmid,ztop,wthv,wqt)
        L0 = max(min(ztop,2500.),500.) / mfparams%L0fac
        if (associated(mfdepth)) mfdepth(IH,JH) = ztop

        ! Reduce L0 over ocean where LTS is large to encourage StCu formation
        lts =  0.0
        if (FRLAND(IH,JH)<0.5) then
           do k = kte-1,kts+1,-1
              if (zlo3(IH,JH,k)-zw3(IH,JH,kte).gt.3000.0) then
                 lts = thv3(IH,JH,k+1)
                 exit
              end if
           end do
           lts = lts - thv3(IH,JH,kte)
           L0 = L0/( 1.0 + (mfparams%ent0lts/mfparams%ent0-1.)*(0.5+0.5*tanh(0.3*(lts-19.))) )
        end if
      else ! if mfparams%ET not 2
        L0 = mfparams%L0
      end if

      !
      ! flipping variables
      !
      DO k=kts,kte
        zlo(k)=zlo3(IH,JH,kte-k+kts)-zw3(IH,JH,kte)
        u(k)=u3(IH,JH,kte-k+kts)
        v(k)=v3(IH,JH,kte-k+kts)
        thl(k)=thl3(IH,JH,kte-k+kts)
        thv(k)=thv3(IH,JH,kte-k+kts)
        qt(k)=qt3(IH,JH,kte-k+kts)
        qv(k)=qv3(IH,JH,kte-k+kts)
        ql(k)=ql3(IH,JH,kte-k+kts)
        qi(k)=qi3(IH,JH,kte-k+kts)
        if (k<kte) then
           if (MFPARAMS%DISCRETE == 0) then
              ui(k)   = 0.5*( u3(IH,JH,kte-k+kts)   + u3(IH,JH,kte-k+kts-1) )
              vi(k)   = 0.5*( v3(IH,JH,kte-k+kts)   + v3(IH,JH,kte-k+kts-1) )
              thli(k) = 0.5*( thl3(IH,JH,kte-k+kts) + thl3(IH,JH,kte-k+kts-1) )
              qti(k)  = 0.5*( qt3(IH,JH,kte-k+kts)  + qt3(IH,JH,kte-k+kts-1) )
              qvi(k)  = 0.5*( qv3(IH,JH,kte-k+kts)  + qv3(IH,JH,kte-k+kts-1) )
              qli(k)  = 0.5*( ql3(IH,JH,kte-k+kts)  + ql3(IH,JH,kte-k+kts-1) )
              qii(k)  = 0.5*( qi3(IH,JH,kte-k+kts)  + qi3(IH,JH,kte-k+kts-1) )
           else
              ui(k)   = u3(IH,JH,kte-k+kts-1)
              vi(k)   = v3(IH,JH,kte-k+kts-1)
              thli(k) = thl3(IH,JH,kte-k+kts-1)
              qti(k)  = qt3(IH,JH,kte-k+kts-1)
              qvi(k)  = qv3(IH,JH,kte-k+kts-1)
              qli(k)  = ql3(IH,JH,kte-k+kts-1)
              qii(k)  = qi3(IH,JH,kte-k+kts-1)
           end if
        end if
      ENDDO
      ui(kte)     = u(kte)
      vi(kte)     = v(kte)
      thli(kte)   = thl(kte)
      qti(kte)    = qt(kte)
      qvi(kte)    = qv(kte)
      qli(kte)    = ql(kte)
      qii(kte)    = qi(kte)
      ui(kts-1)   = u(kts)
      vi(kts-1)   = v(kts)
      thli(kts-1) = thl(kts)  ! approximate
      qti(kts-1)  = qt(kts)
      qvi(kts-1)  = qv(kts)
      qli(kts-1)  = ql(kts)
      qii(kts-1)  = qi(kts)

      DO k=kts-1,kte
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
          call Poisson(1,Nup2,kts,kte,ENTf,ENTi,the_seed)
          do i=1,Nup2
            do k=kts,kte
            ENT(k,i) = (1.-MFPARAMS%STOCHFRAC) * MFPARAMS%Ent0/L0 &
                    + MFPARAMS%STOCHFRAC * real(ENTi(k,i))*MFPARAMS%Ent0/(ZW(k)-ZW(k-1))
            enddo
          enddo
        else if (MFPARAMS%ENTRAIN==1 ) then
          call Poisson(1,Nup2,kts,kte,ENTf,ENTi,the_seed)
          do i=1,Nup2   ! Vary entrainment across updrafts, 0.75-1.25x
            do k=kts,kte
              ENT(k,i) = ((FLOAT(Nup2-i)*0.5/FLOAT(Nup2))+0.75)*( (1.-MFPARAMS%STOCHFRAC) * MFPARAMS%Ent0/L0 &
                      + MFPARAMS%STOCHFRAC * real(ENTi(k,i))*MFPARAMS%Ent0/(ZW(k)-ZW(k-1)) ) !&
            enddo
          enddo
        else if (MFPARAMS%ENTRAIN==2) then
          do i=1,Nup2   ! alternate approach from Soares et al 2004
            do k=kts,kte
              ENT(k,i) = MFPARAMS%Ent0*(1./(ZW(k)+ZW(k)-ZW(k-1))+1./(max(0.,L0-ZW(k))+ZW(k)-ZW(k-1)))
            enddo
          enddo
        end if

      else ! if L0 <= 0
        ENT=0.
      end if

      ! exner function
      exfh=(p/mapl_p00)**mapl_kappa
      exf=(0.5*(p(1:kte)+p(0:kte-1))/mapl_p00)**mapl_kappa

      !
      ! surface conditions
      !
      wstar=max(wstarmin,(mapl_grav*wthv*pblh/300.)**(1./3.))  ! convective velocity scale
      qstar=max(0.,wqt)/wstar
      thstar=max(0.01,wthv)/wstar

      sigmaW=MFPARAMS%AlphaW*wstar
      sigmaQT=MFPARAMS%AlphaQT*qstar
      sigmaTH=MFPARAMS%AlphaTH*thstar

      if (MFPARAMS%doclasp/= 0) then
        wmin=2.*sigmaW
        wmax=2.*sigmaW
      else
        wmin=sigmaW*MFPARAMS%pwmin
        wmax=sigmaW*MFPARAMS%pwmax
      end if

      ! define surface conditions
      DO I=1,NUP2

        wlv=wmin+(wmax-wmin)/(real(NUP2))*(real(i)-1.)
        wtv=wmin+(wmax-wmin)/(real(NUP2))*real(i)

        if (MFPARAMS%doclasp/=0) then
          UPW(kts-1,I) = MFW(IH,JH,I)
          UPA(kts-1,I)=MFAREA(IH,JH,I) !0.5*(ERF(3.0/sqrt(2.))-ERF(1.0/sqrt(2.)))/real(NUP)  ! assume equal size for now
        else
          UPW(kts-1,I)=min(0.5*(wlv+wtv), 5.)
          UPA(kts-1,I)=MIN(1.0,0.5+wthv/0.2)*(0.5*ERF(wtv/(sqrt(2.)*sigmaW))-0.5*ERF(wlv/(sqrt(2.)*sigmaW)))
        end if

        UPU(kts-1,I)=U(kts)
        UPV(kts-1,I)=V(kts)

        if (MFPARAMS%doclasp/=0) then   ! if CLASP, use tile-based perturbations
          UPQT(kts-1,I)=QT(kts)+MFSRCQT(IH,JH,I)
          UPTHV(kts-1,I)=THV(kts)+MFSRCTHL(IH,JH,I)
        else
          UPQT(kts-1,I)=QT(kts)-(-1.**I)*0.32*UPW(kts-1,I)*sigmaQT/sigmaW
!          UPQT(kts-1,I)=QT(kts)+0.32*UPW(kts-1,I)*sigmaQT/sigmaW
          UPTHV(kts-1,I)=THV(kts)+0.58*UPW(kts-1,I)*sigmaTH/sigmaW
        end if

      ENDDO

      !
      ! If needed, rescale UPW to ensure that the mass-flux does not exceed layer mass
      !

      mf = SUM(RHOE(kts-1)*UPA(kts-1,:)*UPW(kts-1,:))
      factor = dp(kts)/(mf*MAPL_GRAV*dt)
      !   print *,'factor=',factor
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

      ! Adjust updraft THV so updraft flux is 90% of surface flux
      if (THVsrfF .gt. 0.9*wthv .and. THVsrfF .gt. 0.1) then
        UPTHV(kts-1,:)=(UPTHV(kts-1,:)-THV(kts))*0.9*wthv/THVsrfF+THV(kts)
      !  print *,'adjusting surface THV perturbation by a factor',0.9*wthv/THVsrfF
      endif

      ! Adjust updraft QT so updraft flux is 90% of surface flux
      IF ( (QTsrfF .gt. 0.9*wqt) .and. (wqt .gt. 0.) )  then
        UPQT(kts-1,:)=(UPQT(kts-1,:)-QT(kts))*0.9*wqt/QTsrfF+QT(kts)
      !  print *,'adjusting surface QT perturbation by a factor',0.9*wqt/QTsrfF
      ENDIF

      ! Compute condensation and initial updraft THL, QL, QI
      DO I=1,NUP2
        call condensation_edmfA(UPTHV(kts-1,I),UPQT(kts-1,I),P(kts-1),      &
                                UPTHL(kts-1,I),UPQL(kts-1,I),UPQI(kts-1,I), &
                                mfparams%ice_ramp)
      ENDDO

      !========================
      !  Integrate updrafts
      !========================

      ! loop over vertical
      vertint: DO k=KTS,KTE
        DO I=1,NUP2  ! loop over updrafts

          if (UPW(K-1,I).gt.0.) then
            if (MFPARAMS%ENTRAIN==3) then  ! dynamic entrainment rates
              ENT(K,I) = MFPARAMS%ENT0*max(1e-4,B)/max(0.1,UPW(K,I)**2)
            elseif (MFPARAMS%ENTRAIN==4) then
            !   ENT(K,I) = (1.-MFPARAMS%STOCHFRAC)*MFPARAMS%Ent0/L0 &
            !              + MFPARAMS%STOCHFRAC*MFPARAMS%ENT0*0.0032/max(0.1,UPW(K-1,I))
            !   ENT(K,I) = 1e-3*(MFPARAMS%Ent0/(max(min(UPW(K-1,I),2.0),0.5))-0.5)
               ENT(K,I) = ENT(K,I)*(1.+3.0*sqrt(TKE3(IH,JH,I))/max(0.1,UPW(K-1,I)))
            end if

            EntExp  = exp(-ENT(K,I)*(ZW(k)-ZW(k-1)))
            EntExpU = exp(-ENT(K,I)*(ZW(k)-ZW(k-1))*MFPARAMS%EntWFac)

            ! Effect of mixing on thermodynamic variables in updraft
            QTn  = QT(K)*(1-EntExp)+UPQT(K-1,I)*EntExp
            THLn = THL(K)*(1-EntExp)+UPTHL(K-1,I)*EntExp
            Un   = U(K)*(1-EntExpU)+UPU(K-1,I)*EntExpU
            Vn   = V(K)*(1-EntExpU)+UPV(K-1,I)*EntExpU

            ! Calculate condensation
            call condensation_edmf(QTn,THLn,P(K),THVn,QCn,wf,mfparams%ice_ramp)

            ! vertical velocity
            B=mapl_grav*(0.5*(THVn+UPTHV(k-1,I))/THV(k)-1.)
            WP=Wb*ENT(K,I)
            IF (WP==0.) THEN
              Wn2=UPW(K-1,I)**2+2.*Wa*B*(ZW(k)-ZW(k-1))
            ELSE
              EntW=exp(-2.*WP*(ZW(k)-ZW(k-1)))
              Wn2=EntW*UPW(k-1,I)**2+Wa*B/WP*(1.-EntW)
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
!                  EXIT vertint
            END IF
          end if ! check if updraft still rising
        ENDDO   ! I: loop over updrafts

        ! Near-surface CFL: To prevent instability, rescale updraft velocities
        ! if mass flux exceeds MFLIMFAC times the layer mass
        if (ZW(k)<300.) then
          mf = SUM(RHOE(k)*UPA(k,:)*UPW(k,:))
          factor = (1.+(MFPARAMS%MFLIMFAC-1.)*(ZW(k)/300.))*dp(K)/(1e-8+mf*MAPL_GRAV*dt)
          if (factor .lt. 1.0) then
             UPW(k,:) = UPW(k,:)*factor
        !                  print *,'rescaling UPW by factor: ',factor
          end if
        end if

      ! loop over vertical
      ENDDO vertint

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
      ! if (factor.ne.1.0) print *,'*** CFL rescale by factor: ',factor
      UPA = factor*UPA

      DO k=KTS,KTE
        edmfmf(IH,JH,KTE-k+KTS-1) = rhoe(K)*SUM(upa(K,:)*upw(K,:))
      ENDDO

      !
      ! writing updraft properties for output
      ! all variables, except Areas are now multipled by the area
      !

      if (associated(EDMF_PLUMES_W))   EDMF_PLUMES_W(IH,JH,KTS-1:KTE,:)   = upw(KTE:KTS-1:-1,:)
      if (associated(EDMF_PLUMES_THL)) EDMF_PLUMES_THL(IH,JH,KTS-1:KTE,:) = upthl(KTE:KTS-1:-1,:)
      if (associated(EDMF_PLUMES_QT))  EDMF_PLUMES_QT(IH,JH,KTS-1:KTE,:)  = upqt(KTE:KTS-1:-1,:)

      dry_a     = 0.
      moist_a   = 0.

      DO k=KTS-1,KTE  ! loop in vertical
        DO I=1,NUP2 ! first sum over all i-updrafts
          IF ((UPQL(K,I)>0.) .OR. UPQI(K,I)>0.)  THEN
            moist_a(K) = moist_a(K)+UPA(K,I)
          ELSE
            dry_a(K) = dry_a(K)+UPA(K,I)
          ENDIF
        ENDDO  ! first sum over all i-updrafts
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

        DO k=KTS-1,KTE  ! loop over vertical
          DO I=1,NUP2   ! sum over all i-updrafts
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
          ENDDO
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

      DO I=1,NUP2
        DO k=KTS-1,KTE
          s_aw(K)=s_aw(K)+UPA(K,I)*UPW(K,I)
          s_aw2(K)=s_aw2(K)+UPA(K,I)*UPW(K,I)*UPW(K,I)
          s_aw3(K)=s_aw3(K)+UPA(K,I)*UPW(K,I)*UPW(K,I)*UPW(K,I)
          s_aqt2(K)=s_aqt2(K)+UPA(K,I)*(UPQT(K,I)-QTI(K))*(UPQT(K,I)-QTI(K))
          s_aqt3(K)=s_aqt3(K)+UPA(K,I)*(UPQT(K,I)-QTI(K))**3
          s_ahlqt(K)=s_ahlqt(K)+exfh(k)*UPA(K,I)*(UPQT(K,I)-QTI(K))*(UPTHL(K,i)-THLI(K))
          ! tmp is dry static energy of updraft, including surface geopotential
          if (MFPARAMS%IMPLICIT == 1) then
             tmp = mapl_cp*exfh(k)*UPTHL(K,i) + mapl_grav*zw(k) + phis(IH,JH) + mapl_alhl*UPQL(K,i) + UPQI(K,I)*mapl_alhs ! updraft S
          else   ! cp*T - Lv*ql + g*z + g*z0 + Lv*Ql + Ls*Qi
             tmp =   mapl_cp*exfh(k)*( UPTHL(K,i) - THLI(K) ) &
                         + mapl_alhl*( UPQL(K,i) - QLI(K) )   &
                         + mapl_alhs*( UPQI(K,I) - QII(K) )
          end if
          ltm=exfh(k)*(UPTHL(K,i)-THLI(K)) !+mapl_grav*zw(k)/mapl_cp
          s_aws(k)  = s_aws(K)+UPA(K,i)*UPW(K,i)*tmp ! for trisolver
          s_ahl2(k) = s_ahl2(K)+UPA(K,i)*ltm*ltm
          s_ahl3(k) = s_ahl3(K)+UPA(K,i)*ltm*ltm*ltm
          s_awhl(k) = s_awhl(K)+UPA(K,i)*UPW(K,I)*ltm
          if (MFPARAMS%IMPLICIT == 1) then
             s_awu(k)  = s_awu(K)  + UPA(K,i)*UPW(K,I)*UPU(K,I)
             s_awv(k)  = s_awv(K)  + UPA(K,i)*UPW(K,I)*UPV(K,I)
             s_awqv(k) = s_awqv(K) + UPA(K,i)*UPW(K,I)*(UPQT(K,I) - UPQI(K,I) - UPQL(K,I))
             s_awql(k) = s_awql(K) + UPA(K,i)*UPW(K,I)*UPQL(K,I)
             s_awqi(k) = s_awqi(K) + UPA(K,i)*UPW(K,I)*UPQI(K,I)
          else
             s_awu(k)  = s_awu(K)  + UPA(K,i)*UPW(K,I)*(UPU(K,I) - UI(K))
             s_awv(k)  = s_awv(K)  + UPA(K,i)*UPW(K,I)*(UPV(K,I) - VI(K))
             s_awqv(k) = s_awqv(K) + UPA(K,i)*UPW(K,I)*(UPQT(K,I) - UPQI(K,I) - UPQL(K,I) - QVI(K))
             s_awql(k) = s_awql(K) + UPA(K,i)*UPW(K,I)*(UPQL(K,I) - QLI(K))
             s_awqi(k) = s_awqi(K) + UPA(K,i)*UPW(K,I)*(UPQI(K,I) - QII(K))
          end if
          s_awqt(k)  = s_awqt(K)  + UPA(K,i)*UPW(K,I)*(UPQT(K,I) - QTI(K))
          mftke(IH,JH,k) = mftke(IH,JH,k) + UPA(KTE+KTS-K-1,i)*0.5*UPW(KTE+KTS-K-1,I)*UPW(KTE+KTS-K-1,I)
        ENDDO

        DO k=KTS,KTE
        ! mass-flux on half levels
        ! need to be careful to treat properly zeros in the UPTHV
          mfthvt=0.5*(UPA(k-1,I)*UPW(k-1,I)*UPTHV(k-1,I)+UPA(k,I)*UPW(k,I)*UPTHV(k,I))
          mft=0.5*(UPA(k-1,I)*UPW(k-1,I)+UPA(k,I)*UPW(k,I))

          s_buoyf(k)=s_buoyf(k)+(mfthvt-mft*THV(k))*exf(k)
        ENDDO
      ENDDO

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


      ! Note values were initialized to zero above
      DO K=KTS-1,KTE-1
        ! outputs - variables needed for solver
        aw3(IH,JH,K)   = s_aw(KTE+KTS-K-1)
        aws3(IH,JH,K)  = s_aws(KTE+KTS-K-1)
        awqv3(IH,JH,K) = s_awqv(KTE+KTS-K-1)
        awql3(IH,JH,K) = s_awql(KTE+KTS-K-1)
        awqi3(IH,JH,K) = s_awqi(KTE+KTS-K-1)
        awu3(IH,JH,K)  = s_awu(KTE+KTS-K-1)
        awv3(IH,JH,K)  = s_awv(KTE+KTS-K-1)
        ae3(IH,JH,K)   = (1.-dry_a(KTE+KTS-K-1)-moist_a(KTE+KTS-K-1))*MFPARAMS%EDfac
        mfwhl(IH,JH,K) = s_awhl(KTE+KTS-K-1)
        mfwqt(IH,JH,K) = s_awqt(KTE+KTS-K-1)
      ENDDO
      mfwhl(IH,JH,KTE) = s_awhl(KTS-1)
      mfwqt(IH,JH,KTE) = s_awqt(KTS-1)

    

      ! buoyancy is defined on full levels
      DO k=kts,kte
        buoyf(IH,JH,K)  = s_buoyf(KTE+KTS-K)    ! can be used in SHOC
        mfw2(IH,JH,K)   = 0.5*(s_aw2(KTE+KTS-K-1)+s_aw2(KTE+KTS-K))
        mfw3(IH,JH,K)   = 0.5*(s_aw3(KTE+KTS-K-1)+s_aw3(KTE+KTS-K))
        mfqt3(IH,JH,K)  = 0.5*(s_aqt3(KTE+KTS-K-1)+s_aqt3(KTE+KTS-K))
        mfhl3(IH,JH,K)  = 0.5*(s_ahl3(KTE+KTS-K-1)+s_ahl3(KTE+KTS-K))
        mfhlqt(IH,JH,K) = 0.5*(s_ahlqt(KTE+KTS-K-1)+s_ahlqt(KTE+KTS-K))
  !      mfhl2(IH,JH,K)=0.5*(s_ahl2(KTE+KTS-K-1)+s_ahl2(KTE+KTS-K))  ! no longer needed
  !      mfqt2(IH,JH,K)=0.5*(s_aqt2(KTE+KTS-K-1)+s_aqt2(KTE+KTS-K))  ! no longer needed
      ENDDO

    END IF   !  IF ( wthv > 0.0 )

  ENDDO ! JH loop over horizontal area
  ENDDO ! IH


END SUBROUTINE run_edmf


subroutine calc_mf_depth(kts,kte,t,z,q,p,ztop,wthv,wqt)

  integer, intent(in   )                     :: kts, kte
!  real,    intent(in   )                     :: ent0
  real,    intent(in   ), dimension(kts:kte) :: t, z, q, p
  real,    intent(in   )                     :: wthv, wqt
  real,    intent(  out)                     :: ztop

  real     :: tep,z1,z2,t1,t2,qp,pp,qsp,dqp,dqsp,wstar,qstar,thstar,sigmaQT,sigmaTH
  integer  :: k

   wstar=max(0.1,(mapl_grav/300.*wthv*1e3)**(1./3.))  ! convective velocity scale
   qstar=max(0.,wqt)/wstar
   thstar=max(0.,wthv)/wstar

!   sigmaW=MFPARAMS%AlphaW*wstar
   sigmaQT=1.0*qstar
   sigmaTH=2.0*thstar

! print *,'sigQT=',sigmaQT,'  sigTH=',sigmaTH,'  wstar=',wstar

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

    qp    = qp  + (0.25/300.)*(z2-z1)*(q(k)-qp)
    tep   = tep + (0.25/300.)*(z2-z1)*(t(k)-tep)

!    print *,'mfdepth: tep=',tep,' pp=',pp
    dqsp  = GEOS_DQSAT(tep , pp , qsat=qsp,  pascals=.true. )

    dqp   = max( qp - qsp, 0. )/(1.+(MAPL_ALHL/MAPL_CP)*dqsp )
    qp    = qp - dqp
    tep   = tep  + MAPL_ALHL * dqp/MAPL_CP


    ! compare Tv env vs parcel
    if ( t2*(1.+MAPL_VIREPS*q(k)) .ge. tep*(1.+MAPL_VIREPS*qp)+0.1 ) then
      ztop = 0.5*(z2+z1)
      exit
    end if

    z1 = z2
    t1 = t2
  enddo  ! k loop

  return

end subroutine calc_mf_depth


subroutine condensation_edmf(QT,THL,P,THV,QC,wf,ice_ramp)
!
! zero or one condensation for edmf: calculates THV and QC
!
use GEOS_UtilsMod, only : GEOS_Qsat

real,intent(in) :: QT,THL,P
real,intent(in) :: ice_ramp
real,intent(out):: THV,QC,wf


integer :: niter,i
real :: diff,exn,t,qs,qcold

! max number of iterations
niter=50
! minimum difference
diff=2.e-5

EXN=(P/mapl_p00)**mapl_kappa
QC=0.

T=EXN*THL

do i=1,NITER
  T=EXN*THL+get_alhl(T,ice_ramp)/mapl_cp*QC
! qsat, p is in pascal
  QS=geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
  QCOLD=QC
  QC=max(0.5*QC+0.5*(QT-QS),0.)
if (abs(QC-QCOLD)<Diff) exit
enddo

T=EXN*THL+get_alhl(T,ice_ramp)/mapl_cp*QC
QS=geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
QC=max(QT-QS,0.)
THV=(THL+get_alhl(T,ice_ramp)/mapl_cp*QC/EXN)*(1.+MAPL_VIREPS*(QT-QC)-QC)
!THV=(THL+get_alhl(T,ice_ramp)/mapl_cp*QC/EXN)*(1.+(mapl_epsilon)*(QT-QC)-QC)
wf=water_f(T,ice_ramp)

end subroutine condensation_edmf


subroutine condensation_edmfA(THV,QT,P,THL,QL,QI,ice_ramp)
!
! zero or one condensation for edmf: calculates QL,QI from THV and QT
!

use GEOS_UtilsMod, only : GEOS_Qsat

real,intent(in) :: THV,QT,P
real,intent(in) :: ice_ramp
real,intent(out):: THL,QL,QI


integer :: niter,i
real :: diff,exn,t,qs,qcold,wf,qc

! max number of iterations
niter=50
! minimum difference
diff=2.e-5

EXN=(P/mapl_p00)**mapl_kappa
QC=0.

!T=EXN*THL

do i=1,NITER
   T=EXN*THV/(1.+MAPL_VIREPS*(QT-QC)-QC)
!   T=EXN*THV/(1.+mapl_epsilon*(QT-QC)-QC)
   QS=geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
   QCOLD=QC
   QC=max(0.5*QC+0.5*(QT-QS),0.)
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

do ih=1,im
  do jh=1,jm
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


function water_f(T,iceramp)
!
! computes water fraction
!
real ::T,iceramp,water_f,Tw
real :: Tmax,Tmin

  Tmax=0.
  Tmin=Tmax-abs(iceramp)
  Tw=T-273.16

! water fraction
  IF (Tw>Tmax) THEN
    water_f=1.
  ELSE IF (Tw<Tmin) THEN
    water_f=0.
  ELSE
    water_f=(Tw-Tmin)/(Tmax-Tmin);
  END IF

end function water_f


subroutine Poisson(istart,iend,jstart,jend,mu,POI,seed)

integer, intent(in) :: istart,iend,jstart,jend
real,dimension(istart:iend,jstart:jend),intent(in) :: MU
integer, dimension(istart:iend,jstart:jend), intent(out) :: POI
integer :: sed_len
integer,dimension(2),  intent(in) :: seed

integer :: seed_len
integer :: IH,JH,idum,p
integer,allocatable :: theseed(:)

call random_seed(SIZE=seed_len)
allocate(theseed(seed_len))

theseed(1:2)=seed
! Gfortran uses longer seeds, so fill the rest with zero
if (seed_len > 2) theseed(3:) = 1+seed(2)


call random_seed(put=theseed)


do ih=istart,iend
  do jh=jstart,jend
!    poi(IH,JH)=poidev(mu(IH,JH),idum)
    poi(IH,JH)=poidev(mu(IH,JH))
  enddo
enddo

end subroutine Poisson



!      FUNCTION poidev(xm,idum)
      FUNCTION poidev(xm)
      INTEGER idum
      REAL poidev,xm,PI
      PARAMETER (PI=3.141592654)
!CU    USES gammln,ran1
      REAL alxm,em,g,oldm,sq,t,y
      SAVE alxm,g,oldm,sq
      DATA oldm /-1./
      if (xm.lt.12.)then
        if (xm.ne.oldm) then
          oldm=xm
          g=exp(-xm)
        endif
        em=-1
        t=1.
2       em=em+1.
!        t=t*ran1(idum)
        t=t*ran1()
        if (t.gt.g) goto 2
      else
        if (xm.ne.oldm) then
          oldm=xm
          sq=sqrt(2.*xm)
          alxm=log(xm)
          g=xm*alxm-gammln(xm+1.)
        endif
!1       y=tan(PI*ran1(idum))
1       y=tan(PI*ran1())
        em=sq*y+xm
        if (em.lt.0.) goto 1
        em=int(em)
        t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
!        if (ran1(idum).gt.t) goto 1
        if (ran1().gt.t) goto 1
      endif
      poidev=em
      return
      END FUNCTION poidev

      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
     24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
     -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END FUNCTION gammln

!      FUNCTION ran1(idum)
      FUNCTION ran1()
      real ran1
!      integer idum

      call random_number(ran1)
      END FUNCTION ran1

end module edmf_mod
