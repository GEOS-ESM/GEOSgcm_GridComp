!#define edmf_diag 1
#define curanddevice
!#define curandhost
module edmf_mod

!
! mass flux updraft parameterization implemented by kay suselj (jpl). additional
! development by nathan arnold and david new (gmao).
!

use mapl_constantsmod, only: mapl_epsilon, mapl_grav, mapl_cp,  &
                             mapl_alhl, mapl_p00, mapl_vireps,  &
                             mapl_alhs, mapl_kappa, mapl_rgas

! use mapl_mod,          only: mapl_undef

! use geos_mod

use geos_utilsmod  ! note : this is added in for specific geos functions

use edmfparams
#ifdef curanddevice
use openacc_curand
#endif
#ifdef curandhost
use curand
#endif

implicit none

real, parameter ::     &
     wstarmin = 1.e-3, &
     zpblmin  = 100.,  &
     onethird = 1./3., &
     r        = 2., &
     mapl_undef = 1.0e15  ! note : this is the value pulled from mapl_mod

public run_edmf

contains

subroutine run_edmf(its,ite,kts,kte,dt,phis, &
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
#ifdef edmf_diag
             w_plume1,w_plume2,w_plume3,w_plume4,w_plume5, &
             w_plume6,w_plume7,w_plume8,w_plume9,w_plume10, &
             qt_plume1,qt_plume2,qt_plume3,qt_plume4,qt_plume5, &
             qt_plume6,qt_plume7,qt_plume8,qt_plume9,qt_plume10, &
             thl_plume1,thl_plume2,thl_plume3,thl_plume4,thl_plume5, &
             thl_plume6,thl_plume7,thl_plume8,thl_plume9,thl_plume10, &
#endif
             params)

! variables needed for solver:
! ae = sum_i (1-a_i)
! aw3 = sum (a_i w_i)
! aws3 = sum (a_i w_i*s_i); s=thl*cp
! aws3,awqt3,awu3,awv3 similar as above except for different variables
!
!
!mass flux variables - diagnostic outputs (on edges):
!   upa,upw,upqt,... kts:kte+1
!  dry_a,moist_a,dry_w,moist_w, ... kts:kte+1
!
!higher-order moments (needed for shoc, on mid-points):
!  buoyf= sum_i  a_i*w_i*(thv_i-<thv)*pi
!  mfw2=sum_i a_i w_i^2
!  mfw3=sum_i a_i w_i^3
!  mfwqt=sum_i a_i w_i qt_i
!  mfqt2=sum_i a_i qt_i^2
!  mfhl2=sum_i a_i h_i^2
!
! three dimensional outputs are on edges as well, but turned around
! dry_a3,moist_a3,dry_thl3, ... (its:ite,kts-1:kte)
! s_aw3,s_awthl3 ... (its:ite,kts-1:kte)

  type (edmfparams_type), intent(in) :: params
  integer, intent(in) :: its,ite,kts,kte,nup!,doclasp
  real,dimension(its:ite,kts:kte), intent(in) :: u3,v3,t3,thl3,qt3,thv3,qv3,ql3,qi3,zlo3
  real,dimension(its:ite,kts-1:kte), intent(in) :: zw3,pw3, rhoe3
  real,dimension(its:ite,kts:kte) :: mfsrcqt,mfsrcthl,mfw,mfarea
  real,dimension(its:ite), intent(in) :: ust2,wthl2,wqt2,pblh2,frland,phis
  real :: dt
  integer,dimension(its:ite) :: nup2

! outputs
  real,dimension(its:ite,kts-1:kte), intent(out) :: dry_a3, moist_a3,dry_w3,moist_w3, &
          dry_qt3,moist_qt3,dry_thl3,moist_thl3,dry_u3,moist_u3,dry_v3,moist_v3,moist_qc3

! #ifdef edmf_diag
!        real,dimension(its:ite,kts-1:kte), intent(out) :: w_plume1,w_plume2,w_plume3,w_plume4, &
!                                                          w_plume5,w_plume6,w_plume7, &
!                                                          w_plume8,w_plume9,w_plume10
!        real,dimension(its:ite,kts-1:kte), intent(out) :: qt_plume1,qt_plume2,qt_plume3,qt_plume4, &
!                                                          qt_plume5,qt_plume6,qt_plume7, &
!                                                          qt_plume8,qt_plume9,qt_plume10
!        real,dimension(its:ite,kts-1:kte), intent(out) :: thl_plume1,thl_plume2,thl_plume3,thl_plume4, &
!                                                          thl_plume5,thl_plume6,thl_plume7, &
!                                                          thl_plume8,thl_plume9,thl_plume10
! #endif

  ! outputs - variables needed for solver (s_aw - sum ai*wi, s_awphi - sum ai*wi*phii)
  real,dimension(its:ite,kts-1:kte), intent(out) :: ae3,aw3,aws3,awqv3,awql3,awqi3,awu3,awv3
   ! output - buoyancy flux: sum_i a_i*w_i*(thv_i-<thv>) ... for tke equation
  real,dimension(its:ite,kts:kte), intent(out) :: buoyf,mfw2,mfw3,mfqt3,mfhl3,mfqt2,mfwqt,mfhl2,&
                                                  mfhlqt,mfwhl,entx
  real, dimension(its:ite,kts-1:kte), intent(out) :: edmfmf
! updraft properties
  real,dimension(kts-1:kte,1:nup,its:ite) :: upw,upthl,upthv,upqt,upa,upu,upv,upqi,upql
 ! entrainment variables
  real,dimension(kts:kte,1:nup) :: entf
  real,dimension(kts:kte,1:nup,its:ite) :: ent
  integer,dimension(kts:kte,1:nup,its:ite) :: enti
! internal variables
  integer :: k,i,ih
  real :: wstar,qstar,thstar,sigmaw,sigmaqt,sigmath,z0, &
          wmin,wmax,wlv,wtv,wp
  real, dimension(its:ite) :: wthv
  real :: b,qtn,thln,thvn,qcn,un,vn,wn2,entexp,entexpu,entw,wf, totent

! internal flipped variables (geos5)

  real,dimension(kts:kte) :: qv,ql,qi
  real,dimension(kts:kte,its:ite) :: zlo,u,v,thl,thv,qt
  real,dimension(kts-1:kte,its:ite) :: zw,p
  real,dimension(kts-1:kte,its:ite) :: thli,qti
  real,dimension(kts-1:kte,its:ite) :: ui, vi,qvi,qli,qii

! internal surface cont
  real :: ust,wthl
  real,dimension(its:ite) :: wqt,pblh
  real,dimension(kts-1:kte,its:ite) :: moist_a,dry_a,dry_w,dry_qt,dry_thl,dry_u,dry_v,&
                                       moist_w,moist_qt,moist_thl,moist_u,moist_v,moist_qc
  real,dimension(kts-1:kte,its:ite) :: s_aw,s_aws,s_awqv,s_awql,s_awqi,s_awu,s_awv
  real,dimension(kts:kte,its:ite) ::  s_buoyf
  real,dimension(kts-1:kte,its:ite) :: s_aw2,s_aw3,s_aqt3,s_ahl3,s_aqt2,s_ahlqt,s_awqt,s_ahl2,s_awhl
! exner function
  real,dimension(kts:kte) :: pmid,t3tmp,zlo3tmp,qv3tmp
  real,dimension(kts:kte,its:ite) :: dp,exf
  real,dimension(kts-1:kte,its:ite) :: rhoe,exfh

  real :: ztop,stmp,ltm,mfsrf,qtsrff,thvsrff,mft,mfthvt,mf,factor
  real,dimension(its:ite) :: l0
  integer, dimension(2) :: seedmf,the_seed

  real :: ma_tot, mw_tot, mqt_tot, mthl_tot, mu_tot, mv_tot, mqc_tot
  real :: da_tot, dw_tot, dqt_tot, dthl_tot, du_tot, dv_tot

  real :: aw_tot, aws_tot, awqv_tot, awql_tot, awqi_tot, awu_tot, awv_tot
  real :: buoyf_tot
  real :: aw2_tot, aw3_tot, aqt3_tot, ahl3_tot, aqt2_tot, ahlqt_tot, awqt_tot, ahl2_tot, awhl_tot

! w parameters
  real,parameter :: &
        wa=1., &
        wb=1.5

! min values to avoid singularities
  real,parameter :: &
      wstarmin=1.e-3, &
      pblhmin=100.

  real :: tstart, tend, time1, time2, time3, time4, time5, time6
  integer :: t2_count, t3_count, t4_count, t5_count, t6_count

  time1 = 0.
  time2 = 0.
  time3 = 0.
  time4 = 0.
  time5 = 0.
  time6 = 0.

  t2_count = 0
  t3_count = 0
  t4_count = 0
  t5_count = 0
  t6_count = 0

!$acc data present(u3,v3,t3,thl3,qt3,thv3,qv3,ql3,qi3,zlo3,zw3,pw3,rhoe3,ust2,wthl2,wqt2,pblh2,frland,phis,&
!$acc             params,params%edfac,params%doclasp,params%et,params%l0fac,params%l0,params%discrete,&
!$acc             params%stochfrac,params%ent0,params%alphaw,params%alphaqt,params%alphath,params%pwmin,&
!$acc             params%pwmax,params%ice_ramp,params%entwfac,params%entrain,params%implicit,&
!$acc             entx,edmfmf,&
!$acc             dry_a3,moist_a3,dry_w3,&
!$acc             moist_w3,dry_qt3,moist_qt3,dry_thl3,moist_thl3,&
!$acc             dry_u3,moist_u3,dry_v3,moist_v3,moist_qc3,aw3,aws3,awqv3,awql3,awqi3,&
!$acc             awu3,awv3,ae3,buoyf,mfw2,mfw3,mfhl2,mfqt2,mfqt3,mfhl3,mfwqt,mfhlqt,mfwhl) &
!$acc      create(mfsrcqt,mfsrcthl,mfw,mfarea,nup2,&
!$acc             moist_a,dry_a,dry_w,dry_qt,dry_thl,dry_u,dry_v,moist_w,moist_qt,&
!$acc             moist_thl,moist_u,moist_v,moist_qc,s_aw,s_aws,s_awqv,s_awql,s_awqi,&
!$acc             s_awu,s_awv,s_aw2,s_aw3,s_aqt3,s_ahl3,&
!$acc             s_aqt2,s_ahlqt,s_awqt,s_ahl2,s_awhl,s_buoyf,&
!$acc             upa,upw,upthl,upthv,upqt,upu,upv,upqi,upql,ent,wthv,&
!$acc             zlo,u,v,thl,thv,qt,zw,p,thli,qti,ui,vi,qvi,qli,qii,wqt,pblh,&
!$acc             dp,exf,rhoe,exfh,l0,enti)

!!!  call cpu_time(tstart)

!$acc parallel
!$acc loop gang vector collapse(2)
  ! temporary, set 
  do k = kts, kte
    do i = its, ite
      mfsrcthl(i,k) = 0.
      mfsrcqt(i,k)  = 0.
      mfw(i,k)      = 0.
      mfarea(i,k)   = 0.

      ! outputs - variables needed for solver
      buoyf(i,k) =0.
      mfw2(i,k)  =0.
      mfw3(i,k)  =0.
      mfqt3(i,k) =0.
      mfhl3(i,k) =0.
      mfqt2(i,k) =0.
      mfwqt(i,k) =0.
      mfhl2(i,k) =0.
      mfhlqt(i,k)=0.
      mfwhl(i,k) =0.
      entx(i,k) = 0.
    enddo
  enddo

!$acc loop gang vector collapse(2) 
  do k = kts-1, kte
    do i = its, ite
      ! set updraft properties to zero/undef
      dry_a3(i,k)=0.
      moist_a3(i,k)=0.
      dry_w3(i,k)=mapl_undef
      moist_w3(i,k)=mapl_undef
      dry_qt3(i,k)=mapl_undef
      moist_qt3(i,k)=mapl_undef
      dry_thl3(i,k)=mapl_undef
      moist_thl3(i,k)=mapl_undef
      dry_u3(i,k)=mapl_undef
      moist_u3(i,k)=mapl_undef
      dry_v3(i,k)=mapl_undef
      moist_v3(i,k)=mapl_undef
      moist_qc3(i,k)=mapl_undef
    enddo
  enddo

!$acc loop gang vector collapse(2)
  do k = kts-1, kte
    do i = its, ite
      ! outputs - variables needed for solver
      aw3(i,k)   =0.
      aws3(i,k)  =0.
      awqv3(i,k) =0.
      awql3(i,k) =0.
      awqi3(i,k) =0.
      awu3(i,k)  =0.
      awv3(i,k)  =0.
      edmfmf(i,k)=0.
      ! this is the environmental area - by default 1.
      ae3(i,k) = params%edfac
    enddo
  enddo
!$acc end parallel

!!!  call cpu_time(tend)

!!!  time1 = tend-tstart

#ifdef edmf_diag
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

  ! sets behavior of geos_qsatlqu an geos_qsatice
  call geos_qsatset

!!!    call cpu_time(tstart)

!$acc parallel private(pmid,t3tmp,zlo3tmp,qv3tmp,ztop,wthl,ust,qv,ql,qi,entf,seedmf,the_seed,&
!$acc                  wstar,qstar,thstar,sigmaw,sigmaqt,sigmath,wmin,wmax,qtsrff,thvsrff,&
!$acc                  entexp,entexpu,qtn,thln,un,vn,b,wp,entw,wn2,factor)
!$acc loop gang
  do ih=its,ite ! loop over the horizontal dimension
    wthl=wthl2(ih)/mapl_cp
    wqt(ih)=wqt2(ih)
    ust=ust2(ih)
    pblh(ih)=pblh2(ih)

    pblh(ih)=max(pblh(ih),pblhmin)
    wthv(ih)=wthl+mapl_epsilon*thv3(ih,kte)*wqt(ih)

    ! if surface buoyancy is positive then mass-flux, otherwise not
    if ( (wthv(ih) > 0.0 .and. params%doclasp==0) .or. (any(mfsrcthl(ih,1:nup) >= -2.0) .and. params%doclasp/=0)) then

      if (params%doclasp/=0) then
        nup2(ih) = count(mfsrcthl(ih,1:nup)>=-2.0)
      else
        nup2(ih) = nup
      end if

!$acc loop vector collapse(2)
      do i = 1, nup
        do k = kts-1, kte
          upw(k,i,ih)=0.
          upthl(k,i,ih)=0.
          upthv(k,i,ih)=0.
          upqt(k,i,ih)=0.
          upa(k,i,ih)=0.
          upu(k,i,ih)=0.
          upv(k,i,ih)=0.
          upqi(k,i,ih)=0.
          upql(k,i,ih)=0.
        enddo
      enddo

!$acc loop vector collapse(2)
      do i = 1, nup
        do k = kts, kte
          ent(k,i,ih) = 0.
        enddo
      enddo

!PRK! Values change. baseDeviceRand
      ! estimate scale height for entrainment calculation
      if (params%et == 2 ) then
!$acc loop seq
        do k = kts,kte
          pmid(k) = 0.5 * (pw3(ih,k-1) + pw3(ih,k))
          t3tmp(k) = t3(ih,k)
          zlo3tmp(k) = zlo3(ih,k)
          qv3tmp(k) = qv3(ih,k)
        enddo

!!!        call calc_mf_depth(kts,kte,t3(ih,:),zlo3(ih,:),qv3(ih,:),pmid,ztop)
        call calc_mf_depth(kts,kte,t3tmp,zlo3tmp,qv3tmp,pmid,ztop)

        l0(ih) = max(min(ztop,3000.),500.) / params%l0fac
      else
        l0(ih) = params%l0
      end if  

!
! flipping variables (geos5)
!
!$acc loop vector
      do k=kts,kte
        zlo(k,ih)=zlo3(ih,kte-k+kts)
        u(k,ih)=u3(ih,kte-k+kts)
        v(k,ih)=v3(ih,kte-k+kts)
        thl(k,ih)=thl3(ih,kte-k+kts)
        thv(k,ih)=thv3(ih,kte-k+kts)
        qt(k,ih)=qt3(ih,kte-k+kts)
        qv(k)=qv3(ih,kte-k+kts)
        ql(k)=ql3(ih,kte-k+kts)
        qi(k)=qi3(ih,kte-k+kts)
        if (k<kte) then
          if (params%discrete == 0) then
            ui(k,ih)   = 0.5*( u3(ih,kte-k+kts)   + u3(ih,kte-k+kts-1) )
            vi(k,ih)   = 0.5*( v3(ih,kte-k+kts)   + v3(ih,kte-k+kts-1) )
            thli(k,ih) = 0.5*( thl3(ih,kte-k+kts) + thl3(ih,kte-k+kts-1) )
            qti(k,ih)  = 0.5*( qt3(ih,kte-k+kts)  + qt3(ih,kte-k+kts-1) )
            qvi(k,ih)  = 0.5*( qv3(ih,kte-k+kts)  + qv3(ih,kte-k+kts-1) )
            qli(k,ih)  = 0.5*( ql3(ih,kte-k+kts)  + ql3(ih,kte-k+kts-1) )
            qii(k,ih)  = 0.5*( qi3(ih,kte-k+kts)  + qi3(ih,kte-k+kts-1) )
          else
            ui(k,ih)   = u3(ih,kte-k+kts-1)
            vi(k,ih)   = v3(ih,kte-k+kts-1)
            thli(k,ih) = thl3(ih,kte-k+kts-1)
            qti(k,ih)  = qt3(ih,kte-k+kts-1)
            qvi(k,ih)  = qv3(ih,kte-k+kts-1)
            qli(k,ih)  = ql3(ih,kte-k+kts-1)
            qii(k,ih)  = qi3(ih,kte-k+kts-1)
          end if
        end if
      enddo
      ui(kte,ih)     = u(kte,ih)
      vi(kte,ih)     = v(kte,ih)
      thli(kte,ih)   = thl(kte,ih)
      qti(kte,ih)    = qt(kte,ih)
      qvi(kte,ih)    = qv(kte)
      qli(kte,ih)    = ql(kte)
      qii(kte,ih)    = qi(kte)
      ui(kts-1,ih)   = u(kts,ih)
      vi(kts-1,ih)   = v(kts,ih)
      thli(kts-1,ih) = thl(kts,ih)  ! approximate
      qti(kts-1,ih)  = qt(kts,ih)
      qvi(kts-1,ih)  = qv(kts)
      qli(kts-1,ih)  = ql(kts)
      qii(kts-1,ih)  = qi(kts)

!$acc loop vector
      do k=kts-1,kte
        rhoe(k,ih) = rhoe3(ih,kte-k+kts-1)
        zw(k,ih)   = zw3(ih,kte-k+kts-1)
        p(k,ih)    = pw3(ih,kte-k+kts-1)
      enddo

!$acc loop vector
      do k = kts, kte
        dp(k,ih) = p(k-1,ih)-p(k,ih)
      enddo

!!!      call cpu_time(tend)

!!!      time2 = time2 + tend - tstart
!!!      t2_count = t2_count + 1

!!!      call cpu_time(tstart)

      !
      ! compute entrainment coefficient
      !


      !
      ! get entrainment type
      !   1 ... probability of entrainment is constant
      !   2 ... probability of entrainment is a function of dthvdz

      ! get dz/l0

!PRK! Values change. baseDeviceRand3
! note : i'm not sure why this nested loop section cannot be included with the previous 
!        parallel section.  for some reason, entf isn't computed properly.
!$acc loop vector
      do k=kts,kte
        do i=1,nup2(ih)
          entf(k,i)=((zw(k,ih)-zw(k-1,ih))/l0(ih))
        enddo
      enddo

      ! get poisson p(dz/l0)
      seedmf(1) = 1000000 * ( 100*thl(kte,ih) - int(100*thl(kte,ih)))
      seedmf(2) = 1000000 * ( 100*thl(kte-1,ih) - int(100*thl(kte-1,ih)))

      the_seed(1)=seedmf(1) + seedmf(2)
      the_seed(2)=seedmf(1) + seedmf(2)
      the_seed(1)=the_seed(1)*seedmf(1)/( seedmf(2) + 10)
      the_seed(2)=the_seed(2)*seedmf(1)/( seedmf(2) + 10)
      if(the_seed(1) == 0) the_seed(1) =  5
      if(the_seed(2) == 0) the_seed(2) = -5

      ! print*, 'l0 = ', l0
      ! print*, 'ztop = ', ztop
      ! call exit(1)

      if (l0(ih) .gt. 0. ) then

!PRK! Values change. baseDeviceRand2
        ! entrainent: ent=ent0/dz*p(dz/l0)
        call poisson(1,nup,kts,kte,entf,enti(:,:,ih),the_seed)

!PRK! Values change. baseDeviceRand2
!$acc loop vector collapse(2)
        do i=1,nup
          do k=kts,kte
            ent(k,i,ih) = ((1.-params%stochfrac) * params%ent0/l0(ih) &
                      + params%stochfrac * real(enti(k,i,ih))*params%ent0/(zw(k,ih)-zw(k-1,ih))) &
                      * (1.+frland(ih))
            ! double entrainment over land to reduce pblh
          enddo
        enddo

        ! increase entrainment if local minimum of thv
        
!$acc loop vector collapse(2)
        do i = 1,nup
          do k=kts+1,kte-1
            if ( (thv(k,ih) .lt. thv(k-1,ih)) .and. (thv(k,ih) .lt. thv(k+1,ih)) ) then
              ent(k,i,ih)=ent(k,i,ih)+5.*params%ent0/l0(ih)
      !          print *,'increasing entrainment, thvs are',thv(k-1:k+1)
            endif 
          enddo
        enddo
        !  if (entrainopt==2) ent(kts+1:,:) = mapl_undef

      else
      ! negative l0 means 0 entrainment
!$acc loop vector collapse(2)
        do i=1,nup
          do k=kts,kte
            ent(k,i,ih) = 0.
          enddo
        enddo
      end if


!!!      call cpu_time(tend)

!!!      time3 = time3 + tend - tstart
!!!      t3_count = t3_count + 1

!!!      call cpu_time(tstart)

      !
      ! exner function
!PRK! Values change. baseDeviceRand1
!$acc loop vector
      do k = kts-1,kte
        exfh(k,ih)=(p(k,ih)/mapl_p00)**mapl_kappa
      enddo

!$acc loop vector
      do k = kts, kte
        exf(k,ih)=(0.5*(p(k,ih)+p(k-1,ih))/mapl_p00)**mapl_kappa
      enddo

      ! surface conditions
      !
!PRK! Values change. baseDeviceRand5
      wstar=max(wstarmin,(mapl_grav/300.*wthv(ih)*pblh(ih))**(1./3.))  ! convective velocity scale
      qstar=wqt(ih)/wstar
      thstar=wthv(ih)/wstar

      sigmaw=params%alphaw*wstar
      sigmaqt=params%alphaqt*qstar
      sigmath=params%alphath*thstar

      if (params%doclasp/= 0) then
        wmin=2.*sigmaw
        wmax=2.*sigmaw

      else
        wmin=sigmaw*params%pwmin
        wmax=sigmaw*params%pwmax
      end if

!PRK! Values changes. baseDeviceRand4
      ! define surface conditions
!$acc loop vector private(wlv,wtv)
      do i=1,nup2(ih)

        wlv=wmin+(wmax-wmin)/(real(nup2(ih)))*(real(i)-1.)
        wtv=wmin+(wmax-wmin)/(real(nup2(ih)))*real(i)

        if (params%doclasp/=0) then
          upw(kts-1,i,ih) = mfw(ih,i)
          upa(kts-1,i,ih)=mfarea(ih,i) !0.5*(erf(3.0/sqrt(2.))-erf(1.0/sqrt(2.)))/real(nup)  ! assume equal size for now
        else
          upw(kts-1,i,ih)=min(0.5*(wlv+wtv), 5.)  ! npa
          upa(kts-1,i,ih)=0.5*erf(wtv/(sqrt(2.)*sigmaw))-0.5*erf(wlv/(sqrt(2.)*sigmaw))
        end if

        upu(kts-1,i,ih)=u(kts,ih)
        upv(kts-1,i,ih)=v(kts,ih)

        if (params%doclasp/=0) then   ! if clasp, use tile-based perturbations
          upqt(kts-1,i,ih)=qt(kts,ih)+mfsrcqt(ih,i)
          upthv(kts-1,i,ih)=thv(kts,ih)+mfsrcthl(ih,i)
        else
          upqt(kts-1,i,ih)=qt(kts,ih)+0.32*upw(kts-1,i,ih)*sigmaqt/sigmaw
          upthv(kts-1,i,ih)=thv(kts,ih)+0.58*upw(kts-1,i,ih)*sigmath/sigmaw
        end if

      enddo

      !
      ! for stability make sure that the surface mass-fluxes are not more than their values computed from the surface scheme
      !
      qtsrff=0.
      thvsrff=0.

!$acc loop seq
      do i=1,nup2(ih)
        qtsrff=qtsrff+upw(kts-1,i,ih)*upa(kts-1,i,ih)*(upqt(kts-1,i,ih)-qt(kts,ih))
        thvsrff=thvsrff+upw(kts-1,i,ih)*upa(kts-1,i,ih)*(upthv(kts-1,i,ih)-thv(kts,ih))
      enddo

      if (thvsrff .gt. wthv(ih) .and. thvsrff .gt. 0.1) then
      ! change surface thv so that the fluxes from the mass flux equal prescribed values
!$acc loop vector
      do i = 1,nup
        upthv(kts-1,i,ih)=(upthv(kts-1,i,ih)-thv(kts,ih))*wthv(ih)/thvsrff+thv(kts,ih)
      enddo
        ! print *,'adjusting surface thv perturbation by a factor',wthv/thvsrff
      endif

      if ( (qtsrff .gt. wqt(ih)) .and. (wqt(ih) .gt. 0.) )  then
      ! change surface qt so that the fluxes from the mass flux equal prescribed values
      ! - we do not need to worry about the negative values as they should not exist -
!$acc loop vector
        do i = 1,nup
          upqt(kts-1,i,ih)=(upqt(kts-1,i,ih)-qt(kts,ih))*wqt(ih)/qtsrff+qt(kts,ih)
        enddo
        ! print *,'adjusting surface qt perturbation by a factor',wqt/qtsrff
      endif

!$acc loop vector
      do i=1,nup2(ih)
        ! compute condensation and thl,ql,qi
        call condensation_edmfa(upthv(kts-1,i,ih),upqt(kts-1,i,ih),p(kts-1,ih), &
        upthl(kts-1,i,ih),upql(kts-1,i,ih),upqi(kts-1,i,ih),params%ice_ramp)
        !print*, 'thl, ql, qi = ', upthl(kts-1,i), upql(kts-1,i), upqi(kts-1,i)
      enddo

!!!      call cpu_time(tend)

!!!      time4 = time4 + tend - tstart
!!!      t4_count = t4_count + 1

!!!      call cpu_time(tstart)
!PRK! Values change. baseDeviceRand8
      !
      ! integrate updrafts
      !
!$acc loop seq
      do i=1,nup2(ih)  ! loop over updrafts
      ! loop over vertical
!$acc loop seq
        vertint: do k=kts,kte
          entexp=exp(-ent(k,i,ih)*(zw(k,ih)-zw(k-1,ih)))
          entexpu=exp(-ent(k,i,ih)*(zw(k,ih)-zw(k-1,ih))*params%entwfac)

          ! thermo-dynamic variables in updraft
          qtn=qt(k,ih)*(1-entexp)+upqt(k-1,i,ih)*entexp
          thln=thl(k,ih)*(1-entexp)+upthl(k-1,i,ih)*entexp
          un=u(k,ih)*(1-entexpu)+upu(k-1,i,ih)*entexpu
          vn=v(k,ih)*(1-entexpu)+upv(k-1,i,ih)*entexpu

          ! note : i'm not sure why condensation_edmf cannot be incorporated to connect the parallel
          !        sections above and below it without errors. i found it's due to the value of qc
          !        not being computed properly when qt > qs
          ! condensation

!!!      write(*,*) "for ih: ",ih,", i: ",i," and k: ",k," the output is thvn: ",thvn,", qcn: ",qcn," and wf: ",wf
          call condensation_edmf(qtn,thln,p(k,ih),thvn,qcn,wf,params%ice_ramp)

          ! vertical velocity
          b=mapl_grav*(0.5*(thvn+upthv(k-1,i,ih))/thv(k,ih)-1.)
          wp=wb*ent(k,i,ih)
          if (wp==0.) then
            wn2=upw(k-1,i,ih)**2+2.*wa*b*(zw(k,ih)-zw(k-1,ih))
          else
            entw=exp(-2.*wp*(zw(k,ih)-zw(k-1,ih)))
            wn2=entw*upw(k-1,i,ih)**2+wa*b/wp*(1.-entw)
          end if

          if (wn2>0.) then
            upw(k,i,ih)=min( sqrt(wn2), 10. ) ! npa
            upthv(k,i,ih)=thvn
            upthl(k,i,ih)=thln
            upqt(k,i,ih)=qtn
            upql(k,i,ih)=qcn*wf
            upqi(k,i,ih)=qcn*(1.-wf)
            upu(k,i,ih)=un
            upv(k,i,ih)=vn
            upa(k,i,ih)=upa(k-1,i,ih)

            if (params%entrain==2 .and. l0(ih)>0.) then
              ent(k+1,i,ih) = params%ent0*max(1e-4,b)/max(0.1,upw(k,i,ih)**2)
            end if
          else
              exit vertint
          end if
            ! loop over vertical
        enddo vertint
      enddo   ! loop over updrafts
!!!      call cpu_time(tend)

!!!      time5 = time5 + tend - tstart
!!!      t5_count = t5_count + 1
!!!      call cpu_time(tstart)


!$acc loop vector private(totent)
      do k=kts,kte
        totent = 0.0
!$acc loop seq
          do i = 1,nup
            totent = totent + ent(k,i,ih)
          enddo
        entx(ih,k) = totent/nup2(ih)
        ! entx(ih,k) = sum(ent(k,:))/nup2
      end do

  ! cfl condition: check that mass flux does not exceed layer mass at any level
  ! if it does, rescale updraft area.
  ! see discussion in beljaars et al 2018 [ecmwf tech memo]

! note : look into this later for better parallelization
!PRK! Values change. baseDeviceRand9
      factor = 1.0
!$acc loop seq private(mf)
      do k=kts,kte
        mf = 0.0
!$acc loop seq
        do i = 1,nup
          mf = mf + rhoe(k,ih)*upa(k,i,ih)*upw(k,i,ih)
        enddo
        if (mf .gt. dp(k,ih)/(mapl_grav*dt)) then
          factor = min(factor,dp(k,ih)/(mf*mapl_grav*dt) )
        end if
      enddo

!$acc loop vector collapse(2)
      do i = 1,nup
        do k = kts-1, kte
          upa(k,i,ih) = factor*upa(k,i,ih)
        enddo
      enddo

!$acc loop vector private(totent)
      do k=kts,kte
        totent = 0.0
!$acc loop seq
        do i = 1,nup
          totent = totent + upa(k,i,ih)*upw(k,i,ih)
        enddo
        edmfmf(ih,k) = rhoe(k,ih)* totent
      enddo

  !
  ! writing updraft properties for output
  ! all variables, except areas are now multipled by the area
  ! to confirm with wrf grid setup we do not save the first and the last row
  !
!$acc loop vector private(ma_tot,mw_tot,mqt_tot,mthl_tot,mu_tot,mv_tot,mqc_tot,&
!$acc                     da_tot,dw_tot,dqt_tot,dthl_tot,du_tot,dv_tot)
      do k=kts-1,kte  ! loop in vertical
#ifdef edmf_diag
        w_plume1(ih,k) = upw(kte+kts-k-1,1,ih)
        w_plume2(ih,k) = upw(kte+kts-k-1,2,ih)
        w_plume3(ih,k) = upw(kte+kts-k-1,3,ih)
        w_plume4(ih,k) = upw(kte+kts-k-1,4,ih)
        w_plume5(ih,k) = upw(kte+kts-k-1,5,ih)
        w_plume6(ih,k) = upw(kte+kts-k-1,6,ih)
        w_plume7(ih,k) = upw(kte+kts-k-1,7,ih)
        w_plume8(ih,k) = upw(kte+kts-k-1,8,ih)
        w_plume9(ih,k) = upw(kte+kts-k-1,9,ih)
        w_plume10(ih,k)= upw(kte+kts-k-1,10,ih)
        qt_plume1(ih,k) = upqt(kte+kts-k-1,1,ih)
        qt_plume2(ih,k) = upqt(kte+kts-k-1,2,ih)
        qt_plume3(ih,k) = upqt(kte+kts-k-1,3,ih)
        qt_plume4(ih,k) = upqt(kte+kts-k-1,4,ih)
        qt_plume5(ih,k) = upqt(kte+kts-k-1,5,ih)
        qt_plume6(ih,k) = upqt(kte+kts-k-1,6,ih)
        qt_plume7(ih,k) = upqt(kte+kts-k-1,7,ih)
        qt_plume8(ih,k) = upqt(kte+kts-k-1,8,ih)
        qt_plume9(ih,k) = upqt(kte+kts-k-1,9,ih)
        qt_plume10(ih,k)= upqt(kte+kts-k-1,10,ih)
        thl_plume1(ih,k) = upthl(kte+kts-k-1,1,ih)
        thl_plume2(ih,k) = upthl(kte+kts-k-1,2,ih)
        thl_plume3(ih,k) = upthl(kte+kts-k-1,3,ih)
        thl_plume4(ih,k) = upthl(kte+kts-k-1,4,ih)
        thl_plume5(ih,k) = upthl(kte+kts-k-1,5,ih)
        thl_plume6(ih,k) = upthl(kte+kts-k-1,6,ih)
        thl_plume7(ih,k) = upthl(kte+kts-k-1,7,ih)
        thl_plume8(ih,k) = upthl(kte+kts-k-1,8,ih)
        thl_plume9(ih,k) = upthl(kte+kts-k-1,9,ih)
        thl_plume10(ih,k)= upthl(kte+kts-k-1,10,ih)
#endif
        ma_tot = 0.
        mw_tot = 0.
        mqt_tot = 0.
        mthl_tot = 0.
        mu_tot = 0.
        mv_tot = 0.
        mqc_tot = 0.

        da_tot = 0.
        dw_tot = 0.
        dqt_tot = 0.
        dthl_tot = 0.
        du_tot = 0.
        dv_tot = 0.

!$acc loop seq
        do i=1,nup2(ih) ! first sum over all i-updrafts
          if ((upql(k,i,ih)>0.) .or. upqi(k,i,ih)>0.)  then
            ma_tot=ma_tot+upa(k,i,ih)
            mw_tot=mw_tot+upa(k,i,ih)*upw(k,i,ih)
            mqt_tot=mqt_tot+upa(k,i,ih)*upqt(k,i,ih)
            mthl_tot=mthl_tot+upa(k,i,ih)*upthl(k,i,ih)
            mu_tot=mu_tot+upa(k,i,ih)*upu(k,i,ih)
            mv_tot=mv_tot+upa(k,i,ih)*upv(k,i,ih)
            mqc_tot=mqc_tot+upa(k,i,ih)*(upql(k,i,ih)+upqi(k,i,ih))
          else
            da_tot=da_tot+upa(k,i,ih)
            dw_tot=dw_tot+upa(k,i,ih)*upw(k,i,ih)
            dqt_tot=dqt_tot+upa(k,i,ih)*upqt(k,i,ih)
            dthl_tot=dthl_tot+upa(k,i,ih)*upthl(k,i,ih)
            du_tot=du_tot+upa(k,i,ih)*upu(k,i,ih)
            dv_tot=dv_tot+upa(k,i,ih)*upv(k,i,ih)
          endif
        enddo  ! first sum over all i-updrafts

        moist_a(k,ih)=ma_tot
        dry_a(k,ih)=da_tot

        if (dry_a(k,ih)>0.) then
          dry_w(k,ih)=dw_tot/dry_a(k,ih)
          dry_qt(k,ih)=dqt_tot/dry_a(k,ih)
          dry_thl(k,ih)=dthl_tot/dry_a(k,ih)
          dry_u(k,ih)=du_tot/dry_a(k,ih)
          dry_v(k,ih)=dv_tot/dry_a(k,ih)
        else
          dry_w(k,ih)=mapl_undef
          dry_qt(k,ih)=mapl_undef
          dry_thl(k,ih)=mapl_undef
          dry_u(k,ih)=mapl_undef
          dry_v(k,ih)=mapl_undef
        endif

        if (moist_a(k,ih)>0.) then
          moist_w(k,ih)=mw_tot/moist_a(k,ih)
          moist_qt(k,ih)=mqt_tot/moist_a(k,ih)
          moist_thl(k,ih)=mthl_tot/moist_a(k,ih)
          moist_u(k,ih)=mu_tot/moist_a(k,ih)
          moist_v(k,ih)=mv_tot/moist_a(k,ih)
          moist_qc(k,ih)=mqc_tot/moist_a(k,ih)
        else
          moist_w(k,ih)=mapl_undef
          moist_qt(k,ih)=mapl_undef
          moist_thl(k,ih)=mapl_undef
          moist_u(k,ih)=mapl_undef
          moist_v(k,ih)=mapl_undef
          moist_qc(k,ih)=mapl_undef
        endif
      enddo     ! loop in vertical

  !
  ! computing variables needed for solver
  !
!$acc loop vector private(aw_tot,aws_tot,awqv_tot,awql_tot,awqi_tot,awu_tot,awv_tot,buoyf_tot,&
!$acc                     aw2_tot,aw3_tot,aqt3_tot,ahl3_tot,aqt2_tot,ahlqt_tot,awqt_tot,ahl2_tot,awhl_tot,&
!$acc                     stmp,ltm,mfthvt,mft)
      do k=kts-1,kte
        aw_tot = 0.
        aws_tot = 0.
        awqv_tot = 0.
        awql_tot = 0.
        awqi_tot = 0.
        awu_tot = 0.
        awv_tot = 0.
        buoyf_tot = 0.
        aw2_tot = 0.
        aw3_tot = 0.
        aqt3_tot = 0.
        ahl3_tot = 0.
        aqt2_tot = 0.
        ahlqt_tot = 0.
        awqt_tot = 0.
        ahl2_tot = 0.
        awhl_tot = 0.
        
!$acc loop seq
        do i=1,nup2(ih)
          aw_tot=aw_tot+upa(k,i,ih)*upw(k,i,ih)
          aw2_tot=aw2_tot+upa(k,i,ih)*upw(k,i,ih)*upw(k,i,ih)
          aw3_tot=aw3_tot+upa(k,i,ih)*upw(k,i,ih)*upw(k,i,ih)*upw(k,i,ih)
          aqt2_tot=aqt2_tot+upa(k,i,ih)*(upqt(k,i,ih)-qti(k,ih))*(upqt(k,i,ih)-qti(k,ih))
          aqt3_tot=aqt3_tot+upa(k,i,ih)*(upqt(k,i,ih)-qti(k,ih))**3
          ahlqt_tot=ahlqt_tot+exfh(k,ih)*upa(k,i,ih)*(upqt(k,i,ih)-qti(k,ih))*(upthl(k,i,ih)-thli(k,ih))
          if (params%implicit == 1) then
            stmp = mapl_cp*exfh(k,ih)*upthl(k,i,ih) + mapl_grav*zw(k,ih) + phis(ih) + mapl_alhl*upql(k,i,ih) + upqi(k,i,ih)*mapl_alhs
          else
!             stmp = exfh(k)*mapl_cp*upthl(k,i) + upqi(k,i)*mapl_alhs + upql(k,i)*mapl_alhl + mapl_grav*zw(k) - exf(k)*mapl_cp*thli(k) - qii(k)*mapl_alhs - qli(k)*mapl_alhl - mapl_grav*zlo(k)
!             stmp = exfh(k)*mapl_cp*upthl(k,i) + upqi(k,i)*mapl_alhs + upql(k,i)*mapl_alhl - exfh(k)*mapl_cp*thli(k) - qii(k)*mapl_alhs - qli(k)*mapl_alhl
            stmp =   mapl_cp*exfh(k,ih)*( upthl(k,i,ih) - thli(k,ih) ) &
                    + mapl_alhl*( upql(k,i,ih) - qli(k,ih) ) &
                    + mapl_alhs*( upqi(k,i,ih) - qii(k,ih) )
          end if
          ltm=exfh(k,ih)*(upthl(k,i,ih)-thli(k,ih)) !+mapl_grav*zw(k)/mapl_cp
          aws_tot=aws_tot+upa(k,i,ih)*upw(k,i,ih)*stmp
          ahl2_tot=ahl2_tot+upa(k,i,ih)*ltm*ltm
          ahl3_tot=ahl3_tot+upa(k,i,ih)*ltm*ltm*ltm
          awhl_tot=awhl_tot+upa(k,i,ih)*upw(k,i,ih)*ltm
          if (params%implicit == 1) then
            awu_tot  = awu_tot  + upa(k,i,ih)*upw(k,i,ih)*upu(k,i,ih)
            awv_tot  = awv_tot  + upa(k,i,ih)*upw(k,i,ih)*upv(k,i,ih)
            awqv_tot = awqv_tot + upa(k,i,ih)*upw(k,i,ih)*(upqt(k,i,ih) - upqi(k,i,ih) - upql(k,i,ih))
            awql_tot = awql_tot + upa(k,i,ih)*upw(k,i,ih)*upql(k,i,ih)
            awqi_tot = awqi_tot + upa(k,i,ih)*upw(k,i,ih)*upqi(k,i,ih)
          else
            awu_tot  = awu_tot  + upa(k,i,ih)*upw(k,i,ih)*(upu(k,i,ih) - ui(k,ih))
            awv_tot  = awv_tot  + upa(k,i,ih)*upw(k,i,ih)*(upv(k,i,ih) - vi(k,ih))
            awqv_tot = awqv_tot + upa(k,i,ih)*upw(k,i,ih)*(upqt(k,i,ih) - upqi(k,i,ih) - upql(k,i,ih) - qvi(k,ih))
            awql_tot = awql_tot + upa(k,i,ih)*upw(k,i,ih)*(upql(k,i,ih) - qli(k,ih))
            awqi_tot = awqi_tot + upa(k,i,ih)*upw(k,i,ih)*(upqi(k,i,ih) - qii(k,ih))
          end if
          awqt_tot  = awqt_tot  + upa(k,i,ih)*upw(k,i,ih)*(upqt(k,i,ih) - qti(k,ih))

        ! mass-flux on half levels
        ! need to be careful to treat properly zeros in the upthv
          if (k > kts-1) then
            mfthvt=0.5*(upa(k-1,i,ih)*upw(k-1,i,ih)*upthv(k-1,i,ih)+upa(k,i,ih)*upw(k,i,ih)*upthv(k,i,ih))
            mft=0.5*(upa(k-1,i,ih)*upw(k-1,i,ih)+upa(k,i,ih)*upw(k,i,ih))
            buoyf_tot=buoyf_tot+(mfthvt-mft*thv(k,ih))*exf(k,ih)
          endif
        enddo

        s_aw(k,ih) = aw_tot
        s_aws(k,ih) = aws_tot
        s_awqv(k,ih) = awqv_tot
        s_awql(k,ih) = awql_tot
        s_awqi(k,ih) = awqi_tot
        s_awu(k,ih) = awu_tot
        s_awv(k,ih) = awv_tot
        s_aw2(k,ih) = aw2_tot
        s_aw3(k,ih) = aw3_tot
        s_aqt3(k,ih) = aqt3_tot
        s_ahl3(k,ih) = ahl3_tot
        s_aqt2(k,ih) = aqt2_tot
        s_ahlqt(k,ih) = ahlqt_tot
        s_awqt(k,ih) = awqt_tot
        s_ahl2(k,ih) = ahl2_tot
        s_awhl(k,ih) = awhl_tot

        if(k > kts-1) then
          s_buoyf(k,ih) = buoyf_tot
        endif
      enddo

!
! turn around the outputs and fill them in the 3d fields
!

!$acc loop vector 
      do k=kts-1,kte
        ! mass-flux diagnostic variables
        dry_a3(ih,k)=dry_a(kte+kts-k-1,ih)
        moist_a3(ih,k)=moist_a(kte+kts-k-1,ih)
        dry_w3(ih,k)=dry_w(kte+kts-k-1,ih)
        moist_w3(ih,k)=moist_w(kte+kts-k-1,ih)
        dry_qt3(ih,k)=dry_qt(kte+kts-k-1,ih)
        moist_qt3(ih,k)=moist_qt(kte+kts-k-1,ih)
        dry_thl3(ih,k)= dry_thl(kte+kts-k-1,ih)
        moist_thl3(ih,k)=moist_thl(kte+kts-k-1,ih)
        dry_u3(ih,k)=dry_u(kte+kts-k-1,ih)
        moist_u3(ih,k)=moist_u(kte+kts-k-1,ih)
        dry_v3(ih,k)=dry_v(kte+kts-k-1,ih)
        moist_v3(ih,k)=moist_v(kte+kts-k-1,ih)
        moist_qc3(ih,k)=moist_qc(kte+kts-k-1,ih)
        ! outputs - variables needed for solver
        aw3(ih,k)=s_aw(kte+kts-k-1,ih)
        aws3(ih,k)=s_aws(kte+kts-k-1,ih)
        awqv3(ih,k)=s_awqv(kte+kts-k-1,ih)
        awql3(ih,k)=s_awql(kte+kts-k-1,ih)
        awqi3(ih,k)=s_awqi(kte+kts-k-1,ih)
        awu3(ih,k)=s_awu(kte+kts-k-1,ih)
        awv3(ih,k)=s_awv(kte+kts-k-1,ih)
        ae3(ih,k)=(1.-dry_a(kte+kts-k-1,ih)-moist_a(kte+kts-k-1,ih))*params%edfac
      enddo

!$acc loop vector
      ! buoyancy is defined on full levels
      do k=kts,kte
          buoyf(ih,k)=s_buoyf(kte+kts-k,ih)

          mfw2(ih,k)=0.5*(s_aw2(kte+kts-k-1,ih)+s_aw2(kte+kts-k,ih))
          mfw3(ih,k)=0.5*(s_aw3(kte+kts-k-1,ih)+s_aw3(kte+kts-k,ih))
          mfhl2(ih,k)=0.5*(s_ahl2(kte+kts-k-1,ih)+s_ahl2(kte+kts-k,ih))
          mfqt2(ih,k)=0.5*(s_aqt2(kte+kts-k-1,ih)+s_aqt2(kte+kts-k,ih))
          mfqt3(ih,k)=0.5*(s_aqt3(kte+kts-k-1,ih)+s_aqt3(kte+kts-k,ih))
          mfhl3(ih,k)=0.5*(s_ahl3(kte+kts-k-1,ih)+s_ahl3(kte+kts-k,ih))
          mfwqt(ih,k)=0.5*(s_awqt(kte+kts-k-1,ih)+s_awqt(kte+kts-k,ih))
          mfhlqt(ih,k)=0.5*(s_ahlqt(kte+kts-k-1,ih)+s_ahlqt(kte+kts-k,ih))
          mfwhl(ih,k)=0.5*(s_awhl(kte+kts-k-1,ih)+s_awhl(kte+kts-k,ih))

      enddo

!!!      call cpu_time(tend)

!!!      time6 = time6 + tend - tstart
!!!      t6_count = t6_count + 1
    end if   !  if ( wthv > 0.0 )
  enddo ! loop over horizontal area
!$acc end parallel
!$acc end data

  ! *** timing data ***
!!!  print*, 'section 1 average time = ', time1
!!!  print*, 'section 2 average time = ', time2/real(t2_count)
!!!  print*, 'section 3 average time = ', time3/real(t3_count)
!!!  print*, 'section 4 average time = ', time4/real(t4_count)
!!!  print*, 'section 5 average time = ', time5/real(t5_count)
!!!  print*, 'section 6 average time = ', time6/real(t6_count)

end subroutine run_edmf


subroutine calc_mf_depth(kts,kte,t,z,q,p,ztop)
!$acc routine seq
  integer, intent(in   )                     :: kts, kte
  real,    intent(in   ), dimension(kts:kte) :: t, z, q, p
  real,    intent(  out)                     :: ztop

  real     :: tep,qp,qsp,dqp,dqsp
  integer  :: k

  tep  = t(kte)+0.4 ! parcel values
  qp   = q(kte)
  ztop = z(kte)

  do k = kte-1 , kts+1, -1
    tep   = tep - mapl_grav*( z(k)-z(k+1) )/mapl_cp
    dqsp  = geos_dqsat(tep , p(k) , qsat=qsp,  pascals=.true. )
    dqp   = max( qp - qsp, 0. )/(1.+(mapl_alhl/mapl_cp)*dqsp )
    qp    = qp - dqp
    tep   = tep  + mapl_alhl * dqp/mapl_cp

    if ( t(k) .ge. tep ) then
      ztop = 0.5*(z(k)+z(k+1))
      exit
    end if

  enddo  ! k loop

  return

end subroutine calc_mf_depth


subroutine condensation_edmf(qt,thl,p,thv,qc,wf,ice_ramp)
!$acc routine seq
  !
  ! zero or one condensation for edmf: calculates thv and qc
  !
  use geos_utilsmod, only : geos_qsat

  real,intent(in) :: qt,thl,p
  real,intent(in) :: ice_ramp
  real,intent(out):: thv,qc,wf


  integer :: niter,i
  real :: diff,exn,t,qs,qcold
  ! max number of iterations
  niter=50
  ! minimum difference
  diff=2.e-5

  exn=(p/mapl_p00)**mapl_kappa
  qc=0.

  t=exn*thl

  do i=1,niter
    t=exn*thl+get_alhl(t,ice_ramp)/mapl_cp*qc
  ! qsat, p is in pascal
    qs=geos_qsat(t,p,pascals=.true.,ramp=ice_ramp)
    qcold=qc
    qc=max(0.5*qc+0.5*(qt-qs),0.)
  if (abs(qc-qcold)<diff) then
    exit
  endif
  enddo
  t=exn*thl+get_alhl(t,ice_ramp)/mapl_cp*qc
  qs=geos_qsat(t,p,pascals=.true.,ramp=ice_ramp)
  qc=max(qt-qs,0.)
  thv=(thl+get_alhl(t,ice_ramp)/mapl_cp*qc/exn)*(1.+mapl_vireps*(qt-qc)-qc)
  !thv=(thl+get_alhl(t,ice_ramp)/mapl_cp*qc/exn)*(1.+(mapl_epsilon)*(qt-qc)-qc)
  wf=water_f(t,ice_ramp)
end subroutine condensation_edmf


subroutine condensation_edmfa(thv,qt,p,thl,ql,qi,ice_ramp)
!$acc routine seq
!
! zero or one condensation for edmf: calculates ql,qi from thv and qt
!

use geos_utilsmod, only : geos_qsat

real,intent(in) :: thv,qt,p
real,intent(in) :: ice_ramp
real,intent(out):: thl,ql,qi


integer :: niter,i
real :: diff,exn,t,qs,qcold,wf,qc



! max number of iterations
niter=50
! minimum difference
diff=2.e-5

exn=(p/mapl_p00)**mapl_kappa
qc=0.

t=exn*thl

do i=1,niter
   t=exn*thv/(1.+mapl_vireps*(qt-qc)-qc)
!   t=exn*thv/(1.+mapl_epsilon*(qt-qc)-qc)
   qs=geos_qsat(t,p,pascals=.true.,ramp=ice_ramp)
 !  print*,'qs = ', qs
   qcold=qc
   qc=max(0.5*qc+0.5*(qt-qs),0.)
if (abs(qc-qcold)<diff) then
  exit
endif
enddo

 thl=(t-get_alhl(t,ice_ramp)/mapl_cp*qc)/exn
 wf=water_f(t,ice_ramp)
 ql=qc*wf
 qi=qc*(1.-wf)

end subroutine condensation_edmfa


function  get_alhl3(t,im,jm,lm,iceramp)

real,dimension(im,jm,lm) ::  t,get_alhl3
real :: iceramp
integer :: im,jm,lm
integer :: i,j,l

do i=1,im
  do j=1,jm
    do l=1,lm
        get_alhl3(i,j,l)=get_alhl(t(i,j,l),iceramp)
    enddo
  enddo
enddo

end function get_alhl3

function get_alhl(t,iceramp)
!$acc routine seq
   real :: t,get_alhl,iceramp,wf
    wf=water_f(t,iceramp)
    get_alhl=wf*mapl_alhl+(1.-wf)*mapl_alhs
end function get_alhl


function water_f(t,iceramp)
!$acc routine seq
!
! computes water fraction
!
real ::t,iceramp,water_f,tw
real :: tmax,tmin

  tmax=0.
  tmin=tmax-abs(iceramp)
  tw=t-273.16

! water fraction
  if (tw>tmax) then
    water_f=1.
  else if (tw<tmin) then
    water_f=0.
  else
    water_f=(tw-tmin)/(tmax-tmin);
  end if

end function water_f


subroutine poisson(istart,iend,jstart,jend,mu,poi,seed)
#ifdef curanddevice
!$acc routine seq nohost
#endif

  integer, intent(in) :: istart,iend,jstart,jend
  real,dimension(istart:iend,jstart:jend),intent(in) :: mu
  integer, dimension(istart:iend,jstart:jend), intent(out) :: poi
  integer :: sed_len
  integer,dimension(2),  intent(in) :: seed

  integer :: seed_len
  integer :: i,j,idum,p
  integer,allocatable :: theseed(:)

#ifdef curanddevice

  type(curandstatexorwow) :: h
!  print*, "random numbers via gpu not implemented"
!  call exit(1)
!   ! using seed from input, seq = 0, offset = 0

   call curand_init(seed(1),0, 0, h)

   do i=istart,iend
     do j=jstart,jend
         poi(i,j)=poidev_gpu(mu(i,j),h)
     enddo
   enddo

#elif defined(curandhost)

  type(curandgenerator) :: gen
  integer :: stat

  stat = curandcreategeneratorhost(gen, curand_rng_pseudo_default)
  stat = curandsetpseudorandomgeneratorseed(gen, seed(1))

  do i=istart,iend
     do j=jstart,jend
         poi(i,j)=poidev_host(mu(i,j),gen,stat)
     enddo
  enddo

#else

  call random_seed(size=seed_len)
  allocate(theseed(seed_len))

  theseed(1:2)=seed
  ! gfortran uses longer seeds, so fill the rest with zero
  if (seed_len > 2) theseed(3:) = seed(2)

  call random_seed(put=theseed)

  do i=istart,iend
  do j=jstart,jend
      poi(i,j)=poidev(mu(i,j),idum)
  enddo
  enddo

#endif

end subroutine poisson

#ifdef curanddevice
      function poidev_gpu(xm,h)
!$acc routine seq nohost
      
        type(curandstatexorwow) :: h
        real poidev_gpu,xm,pi
        parameter (pi=3.141592654)
!cu    uses gammln,ran1
        real alxm,em,g,oldm,sq,t,y
        save alxm,g,oldm,sq
        data oldm /-1./
        if (xm.lt.12.)then
          if (xm.ne.oldm) then
            oldm=xm
            g=exp(-xm)
          endif
          em=-1
          t=1.
  2       em=em+1
          t=t*curand_uniform(h)
          if (t.gt.g) goto 2
        else
          if (xm.ne.oldm) then
            oldm=xm
            sq=sqrt(2.*xm)
            alxm=log(xm)
            g=xm*alxm-gammln(xm+1.)
          endif
  1       y=tan(pi*curand_uniform(h))
          em=sq*y+xm
          if (em.lt.0.) goto 1
          em=int(em)
          t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
          if (curand_uniform(h).gt.t) goto 1
        endif
        poidev_gpu=em
        return
      end function poidev_gpu

      function gammln_gpu(xx)
      real gammln_gpu,xx
      integer j
      double precision ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/76.18009172947146d0,-86.50532032941677d0, &
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
      gammln_gpu=tmp+log(stp*ser/x)
      return
      end function gammln_gpu

#elif defined(curandhost)

      function poidev_host(xm,gen,stat)

        type(curandgenerator) :: gen
        integer :: stat
        real poidev_host,xm,pi
        parameter (pi=3.141592654)
!cu    uses gammln,ran1
        real alxm,em,g,oldm,sq,t,y
        save alxm,g,oldm,sq
        data oldm /-1./
        real, dimension(1) :: randval

        if (xm.lt.12.)then
          if (xm.ne.oldm) then
            oldm=xm
            g=exp(-xm)
          endif
          em=-1
          t=1.
  2       em=em+1
          stat = curandgenerateuniform(gen,randval,1)
          t=t*randval(1)
          if (t.gt.g) goto 2
        else
          if (xm.ne.oldm) then
            oldm=xm
            sq=sqrt(2.*xm)
            alxm=log(xm)
            g=xm*alxm-gammln(xm+1.)
          endif
  1       stat = curandgenerateuniform(gen,randval,1)
          y=tan(pi*randval(1))
          em=sq*y+xm
          if (em.lt.0.) goto 1
          em=int(em)
          t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
          stat = curandgenerateuniform(gen,randval,1)
          if (randval(1).gt.t) goto 1
        endif
        poidev_host=em
        return
      end function poidev_host

#else

! note : nvfortran does not support random_number call from 'ran1' function
      function poidev(xm,idum)
      integer idum
      real poidev,xm,pi
      parameter (pi=3.141592654)
!cu    uses gammln,ran1
      real alxm,em,g,oldm,sq,t,y
      save alxm,g,oldm,sq
      data oldm /-1./
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
1       y=tan(pi*ran1(idum))
        em=sq*y+xm
        if (em.lt.0.) goto 1
        em=int(em)
        t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
        if (ran1(idum).gt.t) goto 1
      endif
      poidev=em
      return
      end function poidev

      function ran1(idum)
      real ran1
      integer idum

      call random_number(ran1)
      end function ran1

#endif

      function gammln(xx)
!$acc routine seq
      real gammln,xx
      integer j
      double precision ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp/76.18009172947146d0,-86.50532032941677d0, &
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
      end function gammln

end module edmf_mod

! nasa docket no. gsc-15,354-1, and identified as "geos-5 gcm modeling software”

! “copyright © 2008 united states government as represented by the administrator
! of the national aeronautics and space administration. all rights reserved.”

! licensed under the apache license, version 2.0 (the "license"); you may not use
! this file except in compliance with the license. you may obtain a copy of the
! license at

! http://www.apache.org/licenses/license-2.0

! unless required by applicable law or agreed to in writing, software distributed
! under the license is distributed on an "as is" basis, without warranties or
! conditions of any kind, either express or implied. see the license for the
! specific language governing permissions and limitations under the license.
