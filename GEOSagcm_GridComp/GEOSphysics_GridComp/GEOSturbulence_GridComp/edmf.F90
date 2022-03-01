!#define EDMF_DIAG 1
module edmf_mod

!
! Mass flux updraft parameterization implemented by Kay Suselj (JPL). Additional
! development by Nathan Arnold and David New (GMAO).
!

use MAPL_ConstantsMod, only: mapl_epsilon, mapl_grav, mapl_cp,  &
                             mapl_alhl, mapl_p00, mapl_vireps,  &
                             mapl_alhs, mapl_kappa, mapl_rgas

use MAPL_Mod,          only: mapl_undef

use GEOS_Mod

use edmfparams

implicit none

real, parameter ::     &
     WSTARmin = 1.e-3, &
     zpblmin  = 100.,  &
     onethird = 1./3., &
     r        = 2.

public run_edmf

contains

SUBROUTINE RUN_EDMF(its,ite,kts,kte,dt,phis, &
              zlo3,zw3,pw3,rhoe3,nup,&
              u3,v3,t3,thl3,thv3,qt3,qv3,ql3,qi3,&
              ust2,wthl2,wqt2,frland,pblh2, &
!              mfsrcthl, mfsrcqt, mfw, mfarea, &
            ! outputs - variables needed for solver
             ae3,aw3,aws3,awqv3,awql3,awqi3,awu3,awv3, &
             mfw2,mfw3,mfqt3,mfhl3,mfwqt,mfqt2,mfhl2,mfhlqt,mfwhl, &
             ! outputs - updraft properties
             dry_a3,moist_a3, &
              dry_w3,moist_w3, &
             dry_qt3,moist_qt3, &
             dry_thl3,moist_thl3, &
             dry_u3,moist_u3, &
             dry_v3,moist_v3, &
             moist_qc3, &
             buoyf, entx, edmfmf, &
#ifdef EDMF_DIAG
             w_plume1,w_plume2,w_plume3,w_plume4,w_plume5, &
             w_plume6,w_plume7,w_plume8,w_plume9,w_plume10, &
             qt_plume1,qt_plume2,qt_plume3,qt_plume4,qt_plume5, &
             qt_plume6,qt_plume7,qt_plume8,qt_plume9,qt_plume10, &
             thl_plume1,thl_plume2,thl_plume3,thl_plume4,thl_plume5, &
             thl_plume6,thl_plume7,thl_plume8,thl_plume9,thl_plume10, &
#endif
             params)

! Variables needed for solver:
! ae = sum_i (1-a_i)
! aw3 = sum (a_i w_i)
! aws3 = sum (a_i w_i*s_i); s=thl*cp
! aws3,awqt3,awu3,awv3 similar as above except for different variables
!
!
!Mass flux variables - diagnostic outputs (on edges):
!   UPA,UPW,UPQT,... KTS:KTE+1
!  dry_a,moist_a,dry_w,moist_w, ... KTS:KTE+1
!
!Higher-order moments (needed for SHOC, on mid-points):
!  buoyf= sum_i  a_i*w_i*(thv_i-<thv)*pi
!  mfw2=sum_i a_i w_i^2
!  mfw3=sum_i a_i w_i^3
!  mfwqt=sum_i a_i w_i qt_i
!  mfqt2=sum_i a_i qt_i^2
!  mfhl2=sum_i a_i h_i^2
!
! three dimensional outputs are on edges as well, but turned around
! dry_a3,moist_a3,dry_thl3, ... (ITS:ITE,KTS-1:KTE)
! s_aw3,s_awthl3 ... (ITS:ITE,KTS-1:KTE)

       type (EDMFPARAMS_TYPE), INTENT(IN) :: PARAMS
       INTEGER, INTENT(IN) :: ITS,ITE,KTS,KTE,NUP!,DOCLASP
       REAL,DIMENSION(ITS:ITE,KTS:KTE), INTENT(IN) :: U3,V3,T3,THL3,QT3,THV3,QV3,QL3,QI3,ZLO3
       REAL,DIMENSION(ITS:ITE,KTS-1:KTE), INTENT(IN) :: ZW3,PW3, rhoe3
       REAL,DIMENSION(ITS:ITE,KTS:KTE) :: mfsrcqt,mfsrcthl,mfw,mfarea
       REAL,DIMENSION(ITS:ITE), INTENT(IN) :: UST2,WTHL2,WQT2,PBLH2,FRLAND,PHIS
       REAL :: DT
       INTEGER :: NUP2

! outputs
       REAL,DIMENSION(ITS:ITE,KTS-1:KTE), INTENT(OUT) :: dry_a3, moist_a3,dry_w3,moist_w3, &
               dry_qt3,moist_qt3,dry_thl3,moist_thl3,dry_u3,moist_u3,dry_v3,moist_v3,moist_qc3

#ifdef EDMF_DIAG
       REAL,DIMENSION(ITS:ITE,KTS-1:KTE), INTENT(OUT) :: w_plume1,w_plume2,w_plume3,w_plume4, &
                                                         w_plume5,w_plume6,w_plume7, &
                                                         w_plume8,w_plume9,w_plume10
       REAL,DIMENSION(ITS:ITE,KTS-1:KTE), INTENT(OUT) :: qt_plume1,qt_plume2,qt_plume3,qt_plume4, &
                                                         qt_plume5,qt_plume6,qt_plume7, &
                                                         qt_plume8,qt_plume9,qt_plume10
       REAL,DIMENSION(ITS:ITE,KTS-1:KTE), INTENT(OUT) :: thl_plume1,thl_plume2,thl_plume3,thl_plume4, &
                                                         thl_plume5,thl_plume6,thl_plume7, &
                                                         thl_plume8,thl_plume9,thl_plume10
#endif

  ! outputs - variables needed for solver (s_aw - sum ai*wi, s_awphi - sum ai*wi*phii)
        REAL,DIMENSION(ITS:ITE,KTS-1:KTE), INTENT(OUT) :: ae3,aw3,aws3,awqv3,awql3,awqi3,awu3,awv3
   ! output - buoyancy flux: sum_i a_i*w_i*(thv_i-<thv>) ... for TKE equation
         REAL,DIMENSION(ITS:ITE,KTS:KTE), INTENT(OUT) :: buoyf,mfw2,mfw3,mfqt3,mfhl3,mfqt2,mfwqt,mfhl2,&
                                                         mfhlqt,mfwhl,entx
      REAL, DIMENSION(ITS:ITE,KTS-1:KTE), INTENT(OUT) :: edmfmf
! updraft properties
      REAL,DIMENSION(KTS-1:KTE,1:NUP) :: UPW,UPTHL,UPQT,UPQL,UPQI,UPA,UPU,UPV,UPTHV
 ! entrainment variables
      REAl,DIMENSION(KTS:KTE,1:NUP) :: ENT,ENTf
      INTEGER,DIMENSION(KTS:KTE,1:NUP) :: ENTi
! internal variables
       INTEGER :: K,I,IH
       REAL :: wthv,wstar,qstar,thstar,sigmaW,sigmaQT,sigmaTH,z0, &
               wmin,wmax,wlv,wtv,wp
       REAL :: B,QTn,THLn,THVn,QCn,Un,Vn,Wn2,EntEXP,EntEXPU,EntW,wf

! internal flipped variables (GEOS5)

       REAL,DIMENSION(KTS:KTE) :: U,V,THL,QT,THV,QV,QL,QI,ZLO
       REAL,DIMENSION(KTS-1:KTE)  :: ZW,P,THLI,QTI
       REAL,DIMENSION(KTS-1:KTE) :: UI, VI, QVI, QLI, QII

! internal surface cont
      REAL :: UST,WTHL,WQT,PBLH
       REAL,DIMENSION(KTS-1:KTE) :: dry_a, moist_a,dry_w,moist_w, &
               dry_qt,moist_qt,dry_thl,moist_thl,dry_u,moist_u,dry_v,moist_v, moist_qc
        REAL,DIMENSION(KTS-1:KTE) :: s_aw,s_aws,s_awqv,s_awql,s_awqi,s_awu,s_awv
        REAL,DIMENSION(KTS:KTE) ::  s_buoyf
        REAL,DIMENSION(KTS-1:KTE) :: s_aw2,s_aw3,s_aqt3,s_ahl3,s_aqt2,s_ahlqt,s_awqt,s_ahl2,s_awhl
! exner function
        REAL,DIMENSION(KTS:KTE) :: exf,dp,pmid
        REAL,DIMENSION(KTS-1:KTE) :: exfh
        REAL,DIMENSION(KTS-1:KTE) :: rhoe

        REAL :: L0,ztop,stmp,ltm,MFsrf,QTsrfF,THVsrfF,mft,mfthvt,mf,factor
        INTEGER, DIMENSION(2) :: seedmf,the_seed


! w parameters
 REAL,PARAMETER :: &
        Wa=1., &
        Wb=1.5

! min values to avoid singularities
  REAL,PARAMETER :: &
     WSTARmin=1.e-3, &
     PBLHmin=100.

     ! temporary, set 
      mfsrcthl = 0.
      mfsrcqt  = 0.
      mfw      = 0.
      mfarea   = 0.

     ! set updraft properties to zero/undef
      dry_a3=0.
      moist_a3=0.
      dry_w3=mapl_undef
      moist_w3=mapl_undef
      dry_qt3=mapl_undef
      moist_qt3=mapl_undef
      dry_thl3=mapl_undef
      moist_thl3=mapl_undef
      dry_u3=mapl_undef
      moist_u3=mapl_undef
      dry_v3=mapl_undef
      moist_v3=mapl_undef
      moist_qc3=mapl_undef
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
      mfqt2 =0.
      mfwqt =0.
      mfhl2 =0.
      mfhlqt=0.
      mfwhl =0.
      entx = 0.
      edmfmf=0.

   ! this is the environmental area - by default 1.

     ae3=PARAMS%EDfac

#ifdef EDMF_DIAG
      w_plume1 = mapl_undef
      w_plume2 = mapl_undef
      w_plume3 = mapl_undef
      w_plume4 = mapl_undef
      w_plume5 = mapl_undef
      w_plume6 = mapl_undef
      w_plume7 = mapl_undef
      w_plume8 = mapl_undef
      w_plume9 = mapl_undef
      w_plume10 = mapl_undef
      qt_plume1 = mapl_undef
      qt_plume2 = mapl_undef
      qt_plume3 = mapl_undef
      qt_plume4 = mapl_undef
      qt_plume5 = mapl_undef
      qt_plume6 = mapl_undef
      qt_plume7 = mapl_undef
      qt_plume8 = mapl_undef
      qt_plume9 = mapl_undef
      qt_plume10 = mapl_undef
      thl_plume1 = mapl_undef
      thl_plume2 = mapl_undef
      thl_plume3 = mapl_undef
      thl_plume4 = mapl_undef
      thl_plume5 = mapl_undef
      thl_plume6 = mapl_undef
      thl_plume7 = mapl_undef
      thl_plume8 = mapl_undef
      thl_plume9 = mapl_undef
      thl_plume10 = mapl_undef
#endif

DO IH=ITS,ITE ! loop over the horizontal dimension


wthl=wthl2(IH)/mapl_cp
wqt=wqt2(IH)
ust=ust2(IH)
pblh=pblh2(IH)

pblh=max(pblh,pblhmin)
wthv=wthl+mapl_epsilon*thv3(IH,kte)*wqt

! if surface buoyancy is positive then mass-flux, otherwise not
  IF ( (wthv > 0.0 .and. PARAMS%doclasp==0) .or. (any(mfsrcthl(IH,1:nup) >= -2.0) .and. PARAMS%doclasp/=0)) then

     if (PARAMS%doclasp/=0) then
       nup2 = count(mfsrcthl(IH,1:nup)>=-2.0)
     else
       nup2 = nup
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
 if (params%ET == 2 ) then
    pmid = 0.5*(pw3(IH,kts-1:kte-1)+pw3(IH,kts:kte))
    call calc_mf_depth(kts,kte,t3(IH,:),zlo3(IH,:),qv3(IH,:),pmid,ztop)
    L0 = max(min(ztop,3000.),500.) / params%L0fac
 else
    L0 = params%L0
 end if   
!
! flipping variables (GEOS5)
!


  DO k=kts,kte
      zlo(k)=zlo3(IH,kte-k+kts)
      u(k)=u3(IH,kte-k+kts)
      v(k)=v3(IH,kte-k+kts)
      thl(k)=thl3(IH,kte-k+kts)
      thv(k)=thv3(IH,kte-k+kts)
      qt(k)=qt3(IH,kte-k+kts)
      qv(k)=qv3(IH,kte-k+kts)
      ql(k)=ql3(IH,kte-k+kts)
      qi(k)=qi3(IH,kte-k+kts)
      if (k<kte) then
         if (PARAMS%DISCRETE == 0) then
            ui(k)   = 0.5*( u3(IH,kte-k+kts)   + u3(IH,kte-k+kts-1) )
            vi(k)   = 0.5*( v3(IH,kte-k+kts)   + v3(IH,kte-k+kts-1) )
            thli(k) = 0.5*( thl3(IH,kte-k+kts) + thl3(IH,kte-k+kts-1) )
            qti(k)  = 0.5*( qt3(IH,kte-k+kts)  + qt3(IH,kte-k+kts-1) )
            qvi(k)  = 0.5*( qv3(IH,kte-k+kts)  + qv3(IH,kte-k+kts-1) )
            qli(k)  = 0.5*( ql3(IH,kte-k+kts)  + ql3(IH,kte-k+kts-1) )
            qii(k)  = 0.5*( qi3(IH,kte-k+kts)  + qi3(IH,kte-k+kts-1) )
         else
            ui(k)   = u3(IH,kte-k+kts-1)
            vi(k)   = v3(IH,kte-k+kts-1)
            thli(k) = thl3(IH,kte-k+kts-1)
            qti(k)  = qt3(IH,kte-k+kts-1)
            qvi(k)  = qv3(IH,kte-k+kts-1)
            qli(k)  = ql3(IH,kte-k+kts-1)
            qii(k)  = qi3(IH,kte-k+kts-1)
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
   rhoe(k) = rhoe3(IH,kte-k+kts-1)
   zw(k)   = zw3(IH,kte-k+kts-1)
   p(k)    = pw3(IH,kte-k+kts-1)
ENDDO

dp = p(kts-1:kte-1)-p(kts:kte)

  !
  ! compute entrainment coefficient
  !


  !
  ! get entrainment type
  !   1 ... probability of entrainment is constant
  !   2 ... probability of entrainment is a function of dTHVdz

  ! get dz/L0
    do i=1,Nup2
      do k=kts,kte
        ENTf(k,i)=((ZW(k)-ZW(k-1))/L0)
      enddo
    enddo

   ! get Poisson P(dz/L0)
  seedmf(1) = 1000000 * ( 100*thl(kte) - INT(100*thl(kte)))
  seedmf(2) = 1000000 * ( 100*thl(kte-1) - INT(100*thl(kte-1)))

  THE_SEED(1)=seedmf(1) + seedmf(2)
  THE_SEED(2)=seedmf(1) + seedmf(2)
  THE_SEED(1)=THE_SEED(1)*seedmf(1)/( seedmf(2) + 10)
  THE_SEED(2)=THE_SEED(2)*seedmf(1)/( seedmf(2) + 10)
  if(THE_SEED(1) == 0) THE_SEED(1) =  5
  if(THE_SEED(2) == 0) THE_SEED(2) = -5


if (L0 .gt. 0. ) then

   ! entrainent: Ent=Ent0/dz*P(dz/L0)
    call Poisson(1,Nup,kts,kte,ENTf,ENTi,the_seed)
    do i=1,Nup
     do k=kts,kte
       ENT(k,i) = (1.-PARAMS%STOCHFRAC) * PARAMS%Ent0/L0 &
                + PARAMS%STOCHFRAC * real(ENTi(k,i))*PARAMS%Ent0/(ZW(k)-ZW(k-1))
     enddo
    enddo
    ENT = (1.+frland(IH))*ENT  ! double entrainment over land to reduce PBLH


! increase entrainment if local minimum of THV

  do k=kts+1,kte-1
    if ( (THV(k) .lt. THV(k-1)) .and. (THV(k) .lt. THV(k+1)) ) then
           ENT(k,:)=ENT(k,:)+5.*PARAMS%ENT0/L0
!          print *,'increasing entrainment, THVs are',THV(k-1:k+1)
     endif
  enddo

!  if (entrainopt==2) ENT(kts+1:,:) = MAPL_UNDEF

else
! negative L0 means 0 entrainment
   ENT=0.
end if

! exner function
 exfh=(p/mapl_p00)**mapl_kappa
 exf=(0.5*(p(1:kte)+p(0:kte-1))/mapl_p00)**mapl_kappa


 !
 ! surface conditions
 !
   wstar=max(wstarmin,(mapl_grav/300.*wthv*pblh)**(1./3.))  ! convective velocity scale
   qstar=wqt/wstar
   thstar=wthv/wstar

   sigmaW=PARAMS%AlphaW*wstar
   sigmaQT=PARAMS%AlphaQT*qstar
   sigmaTH=PARAMS%AlphaTH*thstar

   if (PARAMS%doclasp/= 0) then
     wmin=2.*sigmaW
     wmax=2.*sigmaW

   else
     wmin=sigmaW*PARAMS%pwmin
     wmax=sigmaW*PARAMS%pwmax
   end if

       ! define surface conditions
       DO I=1,NUP2

        wlv=wmin+(wmax-wmin)/(real(NUP2))*(real(i)-1.)
        wtv=wmin+(wmax-wmin)/(real(NUP2))*real(i)

        if (PARAMS%doclasp/=0) then
          UPW(kts-1,I) = MFW(IH,I)
          UPA(kts-1,I)=MFAREA(IH,I) !0.5*(ERF(3.0/sqrt(2.))-ERF(1.0/sqrt(2.)))/real(NUP)  ! assume equal size for now
        else
          UPW(kts-1,I)=min(0.5*(wlv+wtv), 5.)  ! npa
          UPA(kts-1,I)=0.5*ERF(wtv/(sqrt(2.)*sigmaW))-0.5*ERF(wlv/(sqrt(2.)*sigmaW))
        end if

        UPU(kts-1,I)=U(kts)
        UPV(kts-1,I)=V(kts)

        if (PARAMS%doclasp/=0) then   ! if CLASP, use tile-based perturbations
          UPQT(kts-1,I)=QT(kts)+MFSRCQT(IH,I)
          UPTHV(kts-1,I)=THV(kts)+MFSRCTHL(IH,I)
        else
          UPQT(kts-1,I)=QT(kts)+0.32*UPW(kts-1,I)*sigmaQT/sigmaW
          UPTHV(kts-1,I)=THV(kts)+0.58*UPW(kts-1,I)*sigmaTH/sigmaW
        end if

       ENDDO


   !
   ! for stability make sure that the surface mass-fluxes are not more than their values computed from the surface scheme
   !

   QTsrfF=0.
   THVsrfF=0.

   DO I=1,NUP2
     QTsrfF=QTsrfF+UPW(kts-1,I)*UPA(kts-1,I)*(UPQT(kts-1,I)-QT(kts))
     THVsrfF=THVsrfF+UPW(kts-1,I)*UPA(kts-1,I)*(UPTHV(kts-1,I)-THV(kts))
   ENDDO


   if (THVsrfF .gt. wthv .and. THVsrfF .gt. 0.1) then
   ! change surface THV so that the fluxes from the mass flux equal prescribed values
        UPTHV(kts-1,:)=(UPTHV(kts-1,:)-THV(kts))*wthv/THVsrfF+THV(kts)
         print *,'adjusting surface THV perturbation by a factor',wthv/THVsrfF
   endif

   IF ( (QTsrfF .gt. wqt) .and. (wqt .gt. 0.) )  then
   ! change surface QT so that the fluxes from the mass flux equal prescribed values
   ! - we do not need to worry about the negative values as they should not exist -
        UPQT(kts-1,:)=(UPQT(kts-1,:)-QT(kts))*wqt/QTsrfF+QT(kts)
        print *,'adjusting surface QT perturbation by a factor',wqt/QTsrfF
   ENDIF



      DO I=1,NUP2
       ! compute condensation and THL,QL,QI
        call condensation_edmfA(UPTHV(kts-1,i),UPQT(kts-1,I),P(kts-1), &
        UPTHL(kts-1,I),UPQL(kts-1,i),UPQI(kts-1,i),params%ice_ramp)
     ENDDO

  !
  ! integrate updrafts
  !

         DO I=1,NUP2  ! loop over updrafts
         ! loop over vertical
         vertint:   DO k=KTS,KTE


               EntExp=exp(-ENT(K,I)*(ZW(k)-ZW(k-1)))
               EntExpU=exp(-ENT(K,I)*(ZW(k)-ZW(k-1))*PARAMS%EntWFac)

               ! thermo-dynamic variables in updraft
               QTn=QT(K)*(1-EntExp)+UPQT(K-1,I)*EntExp
               THLn=THL(K)*(1-EntExp)+UPTHL(K-1,I)*EntExp
               Un=U(K)*(1-EntExpU)+UPU(K-1,I)*EntExpU
               Vn=V(K)*(1-EntExpU)+UPV(K-1,I)*EntExpU

              ! condensation
               call condensation_edmf(QTn,THLn,P(K),THVn,QCn,wf,params%ice_ramp)

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
                 UPW(K,I)=min( sqrt(Wn2), 10. ) ! npa
                 UPTHV(K,I)=THVn
                 UPTHL(K,I)=THLn
                 UPQT(K,I)=QTn
                 UPQL(K,I)=QCn*wf
                 UPQI(K,I)=QCn*(1.-wf)
                 UPU(K,I)=Un
                 UPV(K,I)=Vn
                 UPA(K,I)=UPA(K-1,I)

               if (PARAMS%ENTRAIN==2 .and. L0>0.) then
                 ENT(K+1,I) = PARAMS%ENT0*max(1e-4,B)/max(0.1,UPW(K,I)**2)
               end if
              ELSE
                  EXIT vertint
              END IF
             ! loop over vertical
            ENDDO vertint
         ENDDO   ! loop over updrafts

         do k=kts,kte
           entx(ih,k) = sum(ENT(k,:))/Nup2
         end do


  ! CFL condition: Check that mass flux does not exceed layer mass at any level
  ! If it does, rescale updraft area.
  ! See discussion in Beljaars et al 2018 [ECMWF Tech Memo]

         factor = 1.0
         DO k=KTS,KTE
            mf = SUM(RHOE(K)*UPA(K,:)*UPW(K,:))
            if (mf .gt. dp(K)/(MAPL_GRAV*dt)) then
               factor = min(factor,dp(K)/(mf*MAPL_GRAV*dt) )
            end if
         ENDDO
         UPA = factor*UPA

         DO k=KTS,KTE
            edmfmf(IH,k) = rhoe(K)*SUM(upa(K,:)*upw(K,:))
         ENDDO

  !
  ! writing updraft properties for output
  ! all variables, except Areas are now multipled by the area
  ! to confirm with WRF grid setup we do not save the first and the last row
  !

      dry_a=0.
      moist_a=0.
      dry_w=0.
      moist_w=0.
      dry_qt=0.
      moist_qt=0.
      dry_thl=0.
      moist_thl=0.
      dry_u=0.
      moist_u=0.
      dry_v=0.
      moist_v=0.
      moist_qc=0.

      DO k=KTS-1,KTE  ! loop in vertical
#ifdef EDMF_DIAG
        w_plume1(IH,k) = upw(kte+kts-k-1,1)
        w_plume2(IH,k) = upw(kte+kts-k-1,2)
        w_plume3(IH,k) = upw(kte+kts-k-1,3)
        w_plume4(IH,k) = upw(kte+kts-k-1,4)
        w_plume5(IH,k) = upw(kte+kts-k-1,5)
        w_plume6(IH,k) = upw(kte+kts-k-1,6)
        w_plume7(IH,k) = upw(kte+kts-k-1,7)
        w_plume8(IH,k) = upw(kte+kts-k-1,8)
        w_plume9(IH,k) = upw(kte+kts-k-1,9)
        w_plume10(IH,k)= upw(kte+kts-k-1,10)
        qt_plume1(IH,k) = upqt(kte+kts-k-1,1)
        qt_plume2(IH,k) = upqt(kte+kts-k-1,2)
        qt_plume3(IH,k) = upqt(kte+kts-k-1,3)
        qt_plume4(IH,k) = upqt(kte+kts-k-1,4)
        qt_plume5(IH,k) = upqt(kte+kts-k-1,5)
        qt_plume6(IH,k) = upqt(kte+kts-k-1,6)
        qt_plume7(IH,k) = upqt(kte+kts-k-1,7)
        qt_plume8(IH,k) = upqt(kte+kts-k-1,8)
        qt_plume9(IH,k) = upqt(kte+kts-k-1,9)
        qt_plume10(IH,k)= upqt(kte+kts-k-1,10)
        thl_plume1(IH,k) = upthl(kte+kts-k-1,1)
        thl_plume2(IH,k) = upthl(kte+kts-k-1,2)
        thl_plume3(IH,k) = upthl(kte+kts-k-1,3)
        thl_plume4(IH,k) = upthl(kte+kts-k-1,4)
        thl_plume5(IH,k) = upthl(kte+kts-k-1,5)
        thl_plume6(IH,k) = upthl(kte+kts-k-1,6)
        thl_plume7(IH,k) = upthl(kte+kts-k-1,7)
        thl_plume8(IH,k) = upthl(kte+kts-k-1,8)
        thl_plume9(IH,k) = upthl(kte+kts-k-1,9)
        thl_plume10(IH,k)= upthl(kte+kts-k-1,10)
#endif
         DO I=1,NUP2 ! first sum over all i-updrafts
            IF ((UPQL(K,I)>0.) .OR. UPQI(K,I)>0.)  THEN
               moist_a(K)=moist_a(K)+UPA(K,I)
               moist_w(K)=moist_w(K)+UPA(K,I)*UPW(K,I)
               moist_qt(K)=moist_qt(K)+UPA(K,I)*UPQT(K,I)
               moist_thl(K)=moist_thl(K)+UPA(K,I)*UPTHL(K,I)
               moist_u(K)=moist_u(K)+UPA(K,I)*UPU(K,I)
               moist_v(K)=moist_v(K)+UPA(K,I)*UPV(K,I)
               moist_qc(K)=moist_qc(K)+UPA(K,I)*(UPQL(K,I)+UPQI(K,I))
            ELSE
               dry_a(K)=dry_a(K)+UPA(K,I)
               dry_w(K)=dry_w(K)+UPA(K,I)*UPW(K,I)
               dry_qt(K)=dry_qt(K)+UPA(K,I)*UPQT(K,I)
               dry_thl(K)=dry_thl(K)+UPA(K,I)*UPTHL(K,I)
               dry_u(K)=dry_u(K)+UPA(K,I)*UPU(K,I)
               dry_v(K)=dry_v(K)+UPA(K,I)*UPV(K,I)
            ENDIF

         ENDDO  ! first sum over all i-updrafts

         IF (dry_a(k)>0.) THEN
            dry_w(k)=dry_w(k)/dry_a(k)
            dry_qt(k)=dry_qt(k)/dry_a(k)
            dry_thl(k)=dry_thl(k)/dry_a(k)
            dry_u(k)=dry_u(k)/dry_a(k)
            dry_v(k)=dry_v(k)/dry_a(k)
         ELSE
            dry_w(k)=mapl_undef
            dry_qt(k)=mapl_undef
            dry_thl(k)=mapl_undef
            dry_u(k)=mapl_undef
            dry_v(k)=mapl_undef
         ENDIF

         IF (moist_a(k)>0.) THEN
            moist_w(k)=moist_w(k)/moist_a(k)
            moist_qt(k)=moist_qt(k)/moist_a(k)
            moist_thl(k)=moist_thl(k)/moist_a(k)
            moist_u(k)=moist_u(k)/moist_a(k)
            moist_v(k)=moist_v(k)/moist_a(k)
            moist_qc(k)=moist_qc(k)/moist_a(k)
         ELSE
            moist_w(k)=mapl_undef
            moist_qt(k)=mapl_undef
            moist_thl(k)=mapl_undef
            moist_u(k)=mapl_undef
            moist_v(k)=mapl_undef
            moist_qc(k)=mapl_undef
         ENDIF

   ENDDO     ! loop in vertical



  !
  ! computing variables needed for solver
  !

     s_aw=0.
     s_aws=0.
     s_awqv=0.
     s_awql=0.
     s_awqi=0.
     s_awu=0.
     s_awv=0.

     s_buoyf=0.
     s_aqt2=0.
     s_awqt=0.
     s_aw2=0.
     s_aw3=0.
     s_aqt3=0.
     s_ahl3=0.
     s_ahl2=0.
     s_awhl=0.
     s_ahlqt=0.

    DO I=1,NUP2

          DO k=KTS-1,KTE
          s_aw(K)=s_aw(K)+UPA(K,I)*UPW(K,I)
          s_aw2(K)=s_aw2(K)+UPA(K,I)*UPW(K,I)*UPW(K,I)
          s_aw3(K)=s_aw3(K)+UPA(K,I)*UPW(K,I)*UPW(K,I)*UPW(K,I)
          s_aqt2(K)=s_aqt2(K)+UPA(K,I)*(UPQT(K,I)-QTI(K))*(UPQT(K,I)-QTI(K))
          s_aqt3(K)=s_aqt3(K)+UPA(K,I)*(UPQT(K,I)-QTI(K))**3
          s_ahlqt(K)=s_ahlqt(K)+exfh(k)*UPA(K,I)*(UPQT(K,I)-QTI(K))*(UPTHL(K,i)-THLI(K))
          if (PARAMS%IMPLICIT == 1) then
             stmp = mapl_cp*exfh(k)*UPTHL(K,i) + mapl_grav*zw(k) + phis(IH) + mapl_alhl*UPQL(K,i) + UPQI(K,I)*mapl_alhs
          else
!             stmp = exfh(k)*mapl_cp*UPTHL(K,i) + UPQI(K,I)*mapl_alhs + UPQL(K,i)*mapl_alhl + mapl_grav*zw(k) - exf(k)*mapl_cp*THLI(K) - QII(K)*mapl_alhs - QLI(K)*mapl_alhl - mapl_grav*zlo(K)
!             stmp = exfh(k)*mapl_cp*UPTHL(K,i) + UPQI(K,I)*mapl_alhs + UPQL(K,i)*mapl_alhl - exfh(k)*mapl_cp*THLI(K) - QII(K)*mapl_alhs - QLI(K)*mapl_alhl
             stmp =   mapl_cp*exfh(k)*( UPTHL(K,i) - THLI(K) ) &
                    + mapl_alhl*( UPQL(K,i) - QLI(K) ) &
                    + mapl_alhs*( UPQI(K,I) - QII(K) )
          end if
          ltm=exfh(k)*(UPTHL(K,i)-THLI(K)) !+mapl_grav*zw(k)/mapl_cp
          s_aws(k)=s_aws(K)+UPA(K,i)*UPW(K,i)*stmp
          s_ahl2(k)=s_ahl2(K)+UPA(K,i)*ltm*ltm
          s_ahl3(k)=s_ahl3(K)+UPA(K,i)*ltm*ltm*ltm
          s_awhl(k)=s_awhl(K)+UPA(K,i)*UPW(K,I)*ltm
          if (PARAMS%IMPLICIT == 1) then
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
   DO K=KTS-1,KTE
      ! mass-flux diagnostic variables
      dry_a3(IH,K)=dry_a(KTE+KTS-K-1)
      moist_a3(IH,K)=moist_a(KTE+KTS-K-1)
      dry_w3(IH,K)=dry_w(KTE+KTS-K-1)
      moist_w3(IH,K)=moist_w(KTE+KTS-K-1)
      dry_qt3(IH,K)=dry_qt(KTE+KTS-K-1)
      moist_qt3(IH,K)=moist_qt(KTE+KTS-K-1)
      dry_thl3(IH,K)= dry_thl(KTE+KTS-K-1)
      moist_thl3(IH,K)=moist_thl(KTE+KTS-K-1)
      dry_u3(IH,K)=dry_u(KTE+KTS-K-1)
      moist_u3(IH,K)=moist_u(KTE+KTS-K-1)
      dry_v3(IH,K)=dry_v(KTE+KTS-K-1)
      moist_v3(IH,K)=moist_v(KTE+KTS-K-1)
      moist_qc3(IH,K)=moist_qc(KTE+KTS-K-1)
      ! outputs - variables needed for solver
      aw3(IH,K)=s_aw(KTE+KTS-K-1)
      aws3(IH,K)=s_aws(KTE+KTS-K-1)
      awqv3(IH,K)=s_awqv(KTE+KTS-K-1)
      awql3(IH,K)=s_awql(KTE+KTS-K-1)
      awqi3(IH,K)=s_awqi(KTE+KTS-K-1)
      awu3(IH,K)=s_awu(KTE+KTS-K-1)
      awv3(IH,K)=s_awv(KTE+KTS-K-1)
      ae3(IH,K)=(1.-dry_a(KTE+KTS-K-1)-moist_a(KTE+KTS-K-1))*PARAMS%EDfac

    ENDDO

! buoyancy is defined on full levels
  DO k=kts,kte
       buoyf(IH,K)=s_buoyf(KTE+KTS-K)

      mfw2(IH,K)=0.5*(s_aw2(KTE+KTS-K-1)+s_aw2(KTE+KTS-K))
      mfw3(IH,K)=0.5*(s_aw3(KTE+KTS-K-1)+s_aw3(KTE+KTS-K))
      mfhl2(IH,K)=0.5*(s_ahl2(KTE+KTS-K-1)+s_ahl2(KTE+KTS-K))
      mfqt2(IH,K)=0.5*(s_aqt2(KTE+KTS-K-1)+s_aqt2(KTE+KTS-K))
      mfqt3(IH,K)=0.5*(s_aqt3(KTE+KTS-K-1)+s_aqt3(KTE+KTS-K))
      mfhl3(IH,K)=0.5*(s_ahl3(KTE+KTS-K-1)+s_ahl3(KTE+KTS-K))
      mfwqt(IH,K)=0.5*(s_awqt(KTE+KTS-K-1)+s_awqt(KTE+KTS-K))
      mfhlqt(IH,K)=0.5*(s_ahlqt(KTE+KTS-K-1)+s_ahlqt(KTE+KTS-K))
      mfwhl(IH,K)=0.5*(s_awhl(KTE+KTS-K-1)+s_awhl(KTE+KTS-K))

  ENDDO

 END IF   !  IF ( wthv > 0.0 )

ENDDO ! loop over horizontal area


END SUBROUTINE run_edmf


subroutine calc_mf_depth(kts,kte,t,z,q,p,ztop)

  integer, intent(in   )                     :: kts, kte
  real,    intent(in   ), dimension(kts:kte) :: t, z, q, p
  real,    intent(  out)                     :: ztop

  real     :: tep,z1,z2,t1,t2,qp,pp,qsp,dqp,dqsp
  integer  :: k

  tep  = t(kte)+0.4 ! parcel values
  qp   = q(kte)

  t1   = t(kte)
  z1   = z(kte)
  ztop = z(kte)

  do k = kte-1 , kts+1, -1
    z2 = z(k)
    t2 = t(k)
    pp = p(k)

    tep   = tep - MAPL_GRAV*( z2-z1 )/MAPL_CP

    dqsp  = GEOS_DQSAT(tep , pp , qsat=qsp,  pascals=.true. )

    dqp   = max( qp - qsp, 0. )/(1.+(MAPL_ALHL/MAPL_CP)*dqsp )
    qp    = qp - dqp
    tep   = tep  + MAPL_ALHL * dqp/MAPL_CP

    if ( t2 .ge. tep ) then
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

T=EXN*THL

do i=1,NITER
   T=EXN*THV/(1.+MAPL_VIREPS*(QT-QC)-QC)
!   T=EXN*THV/(1.+mapl_epsilon*(QT-QC)-QC)
   QS=geos_qsat(T,P,pascals=.true.,ramp=ice_ramp)
   QCOLD=QC
   QC=max(0.5*QC+0.5*(QT-QS),0.)
if (abs(QC-QCOLD)<Diff) exit
enddo

 THL=(T-get_alhl(T,ice_ramp)/mapl_cp*QC)/EXN
 wf=water_f(T,ice_ramp)
 QL=QC*wf
 QI=QC*(1.-wf)

end subroutine condensation_edmfA


function  get_alhl3(T,IM,JM,LM,iceramp)

real,dimension(IM,JM,LM) ::  T,get_alhl3
real :: iceramp
integer :: IM,JM,LM
integer :: I,J,L

do i=1,im
  do j=1,jm
    do l=1,lm
        get_alhl3(i,j,l)=get_alhl(T(i,j,l),iceramp)
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
integer :: i,j,idum,p
integer,allocatable :: theseed(:)

call random_seed(SIZE=seed_len)
allocate(theseed(seed_len))

theseed(1:2)=seed
! Gfortran uses longer seeds, so fill the rest with zero
if (seed_len > 2) theseed(3:) = seed(2)


call random_seed(put=theseed)


do i=istart,iend
 do j=jstart,jend
    poi(i,j)=poidev(mu(i,j),idum)
enddo
 enddo

end subroutine Poisson



      FUNCTION poidev(xm,idum)
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
        t=t*ran1(idum)
        if (t.gt.g) goto 2
      else
        if (xm.ne.oldm) then
          oldm=xm
          sq=sqrt(2.*xm)
          alxm=log(xm)
          g=xm*alxm-gammln(xm+1.)
        endif
1       y=tan(PI*ran1(idum))
        em=sq*y+xm
        if (em.lt.0.) goto 1
        em=int(em)
        t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
        if (ran1(idum).gt.t) goto 1
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

      FUNCTION ran1(idum)
      real ran1
      integer idum

      call random_number(ran1)
      END FUNCTION ran1

end module edmf_mod
