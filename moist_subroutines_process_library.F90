module Process_Library_standalone

    use MAPL_ConstantsMod
    use GEOS_UtilsMod

    implicit none

    private

    interface ICE_FRACTION
        module procedure ICE_FRACTION_3D
        module procedure ICE_FRACTION_2D
        module procedure ICE_FRACTION_1D
        module procedure ICE_FRACTION_SC
    end interface ICE_FRACTION

    ! In anvil/convective clouds
    real, parameter :: aT_ICE_ALL = 245.16
    real, parameter :: aT_ICE_MAX = 261.16
    real, parameter :: aICEFRPWR  = 2.0
    ! Over snow/ice SRF_TYPE = 2
    real, parameter :: iT_ICE_ALL = MAPL_TICE-40.0
    real, parameter :: iT_ICE_MAX = MAPL_TICE
    real, parameter :: iICEFRPWR  = 4.0
    ! Over Land     SRF_TYPE = 1
    real, parameter :: lT_ICE_ALL = 239.16
    real, parameter :: lT_ICE_MAX = 261.16
    real, parameter :: lICEFRPWR  = 2.0
    ! Over Oceans   SRF_TYPE = 0
    real, parameter :: oT_ICE_ALL = 238.16
    real, parameter :: oT_ICE_MAX = 263.16
    real, parameter :: oICEFRPWR  = 4.0

    ! parameters
    real, parameter :: EPSILON =  MAPL_H2OMW/MAPL_AIRMW
    real, parameter :: K_COND  =  2.4e-2    ! J m**-1 s**-1 K**-1
    real, parameter :: DIFFU   =  2.2e-5    ! m**2 s**-1
    ! LDRADIUS4
    real, parameter :: RHO_I   =  916.8     ! Density of ice crystal in kg/m^3
    real, parameter :: RHO_W   = 1000.0     ! Density of liquid water in kg/m^3
    real, parameter :: be      = 1./3. - 0.14
    real, parameter :: bx      = 100.* (3./(4.*MAPL_PI))**(1./3.) * 0.07*6.92
    ! combined constantc
    real, parameter :: cpbgrav = MAPL_CP/MAPL_GRAV
    real, parameter :: gravbcp = MAPL_GRAV/MAPL_CP
    real, parameter :: alhlbcp = MAPL_ALHL/MAPL_CP
    real, parameter :: alhfbcp = MAPL_ALHF/MAPL_CP
    real, parameter :: alhsbcp = MAPL_ALHS/MAPL_CP

    real, parameter :: R_AIR = 3.47e-3 !m3 Pa kg-1K-1

    real, parameter :: mapl_undef = 1.0e15  ! NOTE : This is the value pulled from MAPL_Mod

    public :: hystpdf, ICE_FRACTION, FILLQ2ZERO, BUOYANCY, EVAP3, SUBL3, FIX_UP_CLOUDS, RADCOUPLE

    contains

    subroutine BUOYANCY( T, Q, QS, DQS, DZ, ZLO, BUOY, CAPE, INHB)


        ! !DESCRIPTION: Computes the buoyancy $ g \frac{T_c-T_e}{T_e} $ at each level
        !  for a parcel raised from the surface. $T_c$ is the virtual temperature of
        !  the parcel and $T_e$ is the virtual temperature of the environment.
    
        real, dimension(:,:,:),   intent(in)  :: T, Q, QS, DQS, DZ, ZLO
        real, dimension(:,:,:),   intent(out) :: BUOY
        real, dimension(:,:),     intent(out) :: CAPE, INHB
    
        integer :: I, J, L, IM, JM, LM

        LM = size(T,3)
        IM = size(T,1)
        JM = size(T,2)
!!$acc kernels
!        BUOY(:,:,LM) =  T(:,:,LM) + gravbcp*ZLO(:,:,LM) + alhlbcp*Q(:,:,LM)
!!$acc end kernels

! !$acc parallel loop gang vector collapse(2)
!         do J = 1,JM
!             do I =1,IM
!                 BUOY(I,J,LM) =  T(I,J,LM) + gravbcp*ZLO(I,J,LM) + alhlbcp*Q(I,J,LM)
!             enddo
!         enddo
! !$acc end parallel

!!$acc parallel loop gang vector collapse(3)
!$omp target teams distribute parallel do
        do L=LM-1,1,-1
            do J = 1,JM
                do I = 1,IM
                    BUOY(I,J,L) = (T(I,J,LM) + gravbcp*ZLO(I,J,LM) + alhlbcp*Q(I,J,LM)) - (T(I,J,L) + gravbcp*ZLO(I,J,L) + alhlbcp*QS(I,J,L))
                    BUOY(I,J,L) = MAPL_GRAV*BUOY(I,J,L) / ( (1.+ alhlbcp*DQS(I,J,L))*T(I,J,L) )
                enddo
            enddo
        enddo
!$omp end target teams distribute parallel do
!!$acc end parallel loop

!!$acc parallel loop gang vector collapse(2)
!$omp target teams distribute parallel do
        do J = 1,JM
            do I = 1,IM
                BUOY(I,J,LM) = 0.0
            
                CAPE(I,J) = 0.
                INHB(I,J) = 0.
            enddo
        enddo
!$omp end target teams distribute parallel do
!!$acc end parallel loop

!!$acc parallel loop gang vector collapse(3)
!$omp target teams distribute parallel do
        do L=1,LM-1
            do J = 1, JM
                do I = 1,IM
                    if(BUOY(I,J,L)>0.) then
                    !!$acc atomic update
                    !$omp atomic update
                        CAPE(I,J) = CAPE(I,J) + BUOY(I,J,L)*DZ(I,J,L)

                    !!$acc end atomic
                    endif 

                   if (BUOY(I,J,L)<0.) then
                        !!$acc atomic update
                        !$omp atomic update
                        INHB(I,J) = INHB(I,J) - BUOY(I,J,L)*DZ(I,J,L)

                        !!$acc end atomic
                    endif
                enddo
            enddo
        end do
!$omp end target teams distribute parallel do
!!$acc end parallel loop

!!$acc parallel loop gang vector collapse(2)
!$omp target teams distribute parallel do
        do J = 1,JM
            do I = 1,IM
                if(CAPE(I,J) <= 0.0) then
                    CAPE(I,J) = mapl_undef
                    INHB(I,J) = mapl_undef
                endif
            enddo
        enddo
!$omp end target teams distribute parallel do
!!$acc end parallel loop
    end subroutine BUOYANCY

    subroutine pdffrac (flag,qtmean,sigmaqt1,sigmaqt2,qstar,clfrac)
        implicit none
  
        integer flag            ! flag to indicate shape of pdf
                                ! 1 for tophat, 2 for triangular, 3 for Gaussian
        real qtmean             ! Grid box value of q total
        real sigmaqt1           ! width of distribution (sigma)
        real sigmaqt2           ! width of distribution (sigma)
        real qstar              ! saturation q at grid box avg T
        real clfrac             ! cloud fraction (area under pdf from qs)
  
        real :: qtmode, qtmin, qtmax
  
        if(flag.eq.1) then
            if((qtmean+sigmaqt1).lt.qstar) then
                clfrac = 0.
            else
                if(sigmaqt1.gt.0.) then
                    clfrac = min((qtmean + sigmaqt1 - qstar),2.*sigmaqt1)/(2.*sigmaqt1)
                else
                    clfrac = 1.
                endif
            endif
        elseif(flag.eq.2) then
            qtmode =  qtmean + (sigmaqt1-sigmaqt2)/3.
            qtmin = max(qtmode-sigmaqt1,0.)
            qtmax = qtmode + sigmaqt2
            if(qtmax.le.qstar) then
                clfrac = 0.
            elseif ( (qtmode.le.qstar).and.(qstar.lt.qtmax) ) then
                clfrac = (qtmax-qstar)*(qtmax-qstar) / ( (qtmax-qtmin)*(qtmax-qtmode) )
            elseif ( (qtmin.le.qstar).and.(qstar.lt.qtmode) ) then
                clfrac = 1. - ( (qstar-qtmin)*(qstar-qtmin) / ( (qtmax-qtmin)*(qtmode-qtmin) ) )
            elseif ( qstar.le.qtmin ) then
                clfrac = 1.
            endif
        endif
  
        return
    end subroutine pdffrac

    subroutine pdfcondensate (flag,qtmean4,sigmaqt14,sigmaqt24,qstar4,condensate4)
        implicit none
  
        integer flag            ! flag to indicate shape of pdf
                                ! 1 for tophat, 2 for triangular
        real qtmean4            ! Grid box value of q total
        real sigmaqt14          ! width of distribution (to left)
        real sigmaqt24          ! width of distribution (to right)
        real qstar4             ! saturation q at grid box avg T
        real condensate4        ! condensate (area under (q*-qt)*pdf from qs)
  
        real *8 :: qtmode, qtmin, qtmax, constA, constB, cloudf
        real *8 :: term1, term2, term3
        real *8 :: qtmean, sigmaqt1, sigmaqt2, qstar, condensate
  
        qtmean = dble(qtmean4)
        sigmaqt1 = dble(sigmaqt14)
        sigmaqt2 = dble(sigmaqt24)
        qstar = dble(qstar4)
  
        if(flag.eq.1) then
            if(qtmean+sigmaqt1.lt.qstar) then
                condensate = 0.d0
            elseif(qstar.gt.qtmean-sigmaqt1)then
                if(sigmaqt1.gt.0.d0) then
                    condensate = (min(qtmean + sigmaqt1 - qstar,2.d0*sigmaqt1)**2)/ (4.d0*sigmaqt1)
                else
                    condensate = qtmean-qstar
                endif
            else
                condensate = qtmean-qstar
            endif
        elseif(flag.eq.2) then
            qtmode =  qtmean + (sigmaqt1-sigmaqt2)/3.d0
            qtmin = max(qtmode-sigmaqt1,0.d0)
            qtmax = qtmode + sigmaqt2
            if ( qtmax.le.qstar ) then
                condensate = 0.d0
            elseif ( (qtmode.le.qstar).and.(qstar.lt.qtmax) ) then
                constB = 2.d0 / ( (qtmax - qtmin)*(qtmax-qtmode) )
                cloudf = (qtmax-qstar)*(qtmax-qstar) / ( (qtmax-qtmin)*(qtmax-qtmode) )
                term1 = (qstar*qstar*qstar)/3.d0
                term2 = (qtmax*qstar*qstar)/2.d0
                term3 = (qtmax*qtmax*qtmax)/6.d0
                condensate = constB * (term1-term2+term3) - qstar*cloudf
            elseif ( (qtmin.le.qstar).and.(qstar.lt.qtmode) ) then
                constA = 2.d0 / ( (qtmax - qtmin)*(qtmode-qtmin) )
                cloudf = 1.d0 - ( (qstar-qtmin)*(qstar-qtmin) / ( (qtmax-qtmin)*(qtmode-qtmin) ) )
                term1 = (qstar*qstar*qstar)/3.d0
                term2 = (qtmin*qstar*qstar)/2.d0
                term3 = (qtmin*qtmin*qtmin)/6.d0
                condensate = qtmean - ( constA * (term1-term2+term3) ) - qstar*cloudf
            elseif ( qstar.le.qtmin ) then
                condensate = qtmean-qstar
            endif
        endif
        condensate4 = real(condensate)
  
        return
    end subroutine pdfcondensate

    subroutine partition_dblgss( dt,           &  ! IN
                                tabs,         &  ! INOUT
                                qwv,          &
                                qc,           &
                        !                              qi,           &
                                omega,        &  ! IN
                                zl,           &
                                pval,         &
                                total_water,  &
                                thl_first,    &
                                wthlsec,      &
                                wqwsec,       &
                        !                              wqtfac,       & ! inter-gaussian qt flux
                        !                              whlfac,       & ! inter-gaussian hl flux
                                thlsec,       &
                                qwsec,        &
                                qwthlsec,     &
                                w3var,        &
                                w_sec,        &
                                qt3,          &
                                hl3,          &
                                mffrc,        &
                                PDF_A,        &  ! INOUT
#ifdef PDFDIAG
                                PDF_SIGW1,    &  ! OUT - diagnostic only
                                PDF_SIGW2,    &
                                PDF_W1,       &
                                PDF_W2,       &
                                PDF_SIGTH1,   &
                                PDF_SIGTH2,   &
                                PDF_TH1,      &
                                PDF_TH2,      &
                                PDF_SIGQT1,   &
                                PDF_SIGQT2,   &
                                PDF_QT1,      &
                                PDF_QT2,      &
                                PDF_RQTTH,    &
                                PDF_RWTH,     &
                                PDF_RWQT,     &
#endif
                                wthv_sec,     &  ! OUT - needed elsewhere
                                wqls,         &
                                cld_sgs)

        use MAPL_ConstantsMod, only: ggr    => MAPL_GRAV,   &
                                    cp     => MAPL_CP,     &
                                    rgas   => MAPL_RGAS,   &
                                    rv     => MAPL_RVAP,   &
                                    lcond  => MAPL_ALHL,   &
                                    lfus   => MAPL_ALHS,   &
                                    pi     => MAPL_PI,     &
                                    MAPL_H2OMW, MAPL_AIRMW
        use MAPL_SatVaporMod,  only: MAPL_EQsat

        real, intent(in   )  :: DT          ! timestep [s]
        real, intent(in   )  :: tabs        ! absolute temperature [K]
        real, intent(in   )  :: qwv         ! specific humidity [kg kg-1]
        real, intent(  out)  :: qc          ! liquid+ice condensate [kg kg-1]
        real, intent(in   )  :: omega       ! resolved pressure velocity
        real, intent(in   )  :: zl          ! layer heights [m]
        real, intent(in   )  :: pval        ! layer pressure [Pa]
        real, intent(in   )  :: total_water ! total water [kg kg-1]
        real, intent(in   )  :: thl_first   ! liquid water potential temperature [K]
        real, intent(in   )  :: wthlsec     ! thl flux [K m s-1]
        real, intent(in   )  :: wqwsec      ! total water flux [kg kg-1 m s-1]
        !   real, intent(in   )  :: wqtfac      !
        !   real, intent(in   )  :: whlfac      !
        real, intent(in   )  :: thlsec
        real, intent(in   )  :: qwsec
        real, intent(in   )  :: qwthlsec
        real, intent(in   )  :: w3var       ! 3rd moment vertical velocity [m3 s-3]
        real, intent(in   )  :: qt3         ! 3rd moment qt from mass flux
        real, intent(in   )  :: hl3         ! 3rd moment hl from mass flux
        real, intent(in   )  :: w_sec       ! 2nd moment vertical velocity [m2 s-2]
        real, intent(in   )  :: mffrc       ! total EDMF updraft fraction
        !   real, intent(inout)  :: qi         ! ice condensate [kg kg-1]
        real, intent(  out)  :: cld_sgs     ! cloud fraction
        real, intent(inout)  ::    PDF_A           ! fractional area of 1st gaussian
#ifdef PDFDIAG
        real, intent(  out)  ::    PDF_SIGW1,    & ! std dev w of 1st gaussian [m s-1]
                                   PDF_SIGW2,    & ! std dev w of 2nd gaussian
                                   PDF_W1,       & ! mean vertical velocity of 1st gaussian [m s-1]
                                   PDF_W2,       & ! mean vertical velocity of 2nd gaussian [m s-1]
                                   PDF_SIGTH1,   & ! std dev pot temp of 1st gaussian [K]
                                   PDF_SIGTH2,   & ! std dev pot temp of 2nd gaussian [K]
                                   PDF_TH1,      & ! mean pot temp of 1st gaussian [K]
                                   PDF_TH2,      & ! mean pot temp of 2nd gaussian [K]
                                   PDF_SIGQT1,   & ! std dev total water of 1st gaussian [kg kg-1]
                                   PDF_SIGQT2,   & ! std dev total water of 2nd gaussian [kg kg-1]
                                   PDF_QT1,      & ! mean total water of 1st gaussian [kg kg-1]
                                   PDF_QT2,      & ! mean total water of 2nd gaussian [kg kg-1]
                                   PDF_RQTTH,    & ! QT-TH correlation coeff
                                   PDF_RWTH,     & ! W-TH correlation
                                   PDF_RWQT        ! W-QT correlation
#endif
        real, intent(  out)  :: wthv_sec
        real, intent(  out)  :: wqls


        ! Local variables

        integer i,j,k,ku,kd
        real wrk, wrk1, wrk2, wrk3, wrk4, bastoeps
        real gamaz, thv, rwqt, rwthl, wql1, wql2
        real pkap, diag_qn, diag_frac, diag_ql, diag_qi,w_first,                     &
            sqrtw2, sqrtthl, sqrtqt, w1_1, w1_2, w2_1, w2_2, thl1_1, thl1_2,        &
            thl2_1, thl2_2, qw1_1, qw1_2, qw2_1, qw2_2, aterm, onema, sm,           &
            km1, skew_w, skew_qw, skew_thl, cond_w, sqrtw2t,                                 &
            sqrtthl2_1, sqrtthl2_2, sqrtqw2_1, sqrtqw2_2, corrtest1, corrtest2,     &
            tsign, testvar, r_qwthl_1, Tl1_1, Tl1_2, esval1_1, esval1_2, esval2_1,  &
            esval2_2, om1, om2, lstarn1, lstarn2, qs1, qs2, beta1, beta2, cqt1,     &
            cqt2, s1, s2, cthl1, cthl2, std_s1, std_s2, qn1, qn2, C1, C2, ql1, ql2, &
            qi1, qi2, wqis, wqtntrgs, whlntrgs


        ! Set constants and parameters
        real, parameter :: sqrt2 = sqrt(2.0)
        real, parameter :: sqrtpii = 1.0/sqrt(pi+pi)
        real, parameter :: tbgmin = 233.16
        real, parameter :: tbgmax = MAPL_TICE
        real, parameter :: a_bg   = 1.0/(tbgmax-tbgmin)
        real, parameter :: thl_tol = 1.e-2
        real, parameter :: w_thresh = 0.001
        real, parameter :: rt_tol = 1.e-4
        real, parameter :: w_tol_sqd = 4.0e-04   ! Min vlaue of second moment of w
        real, parameter :: onebrvcp = 1.0/(rv*cp)
        real, parameter :: skew_facw = 1.2
        real, parameter :: skew_fact = 0.5
        real, parameter :: lsub = lcond+lfus
        real, parameter :: fac_cond = lcond/cp
        real, parameter :: fac_sub = lsub/cp
        real, parameter :: fac_fus = lfus/cp
        real, parameter :: gocp = ggr/cp
        real, parameter :: rog = rgas / ggr
        real, parameter :: kapa = rgas / cp
        real, parameter :: epsv=MAPL_H2OMW/MAPL_AIRMW

        real, parameter :: use_aterm_memory = 1.
        real, parameter :: tauskew = 3600. 

        ! define conserved variables
        gamaz = gocp * zl
        thv   = tabs * (1.0+epsv*qwv)
        thv   = thv*(100000.0/pval) ** kapa

        w_first = - rog * omega * thv / pval

        ! Initialize cloud variables to zero
        diag_qn   = 0.0
        diag_frac = 0.0
        diag_ql   = 0.0
        diag_qi   = 0.0

        pkap = (pval/100000.0) ** kapa


        ! Compute square roots of some variables so we don't have to do it again
        if (w_sec > 0.0) then
            sqrtw2   = sqrt(w_sec)
            Skew_w   = w3var / (sqrtw2*sqrtw2*sqrtw2)
        else
            sqrtw2   = w_thresh
            Skew_w   = 0.
        endif
        if (thlsec > 0.0) then
            sqrtthl  = sqrt(thlsec)
            skew_thl = hl3 / sqrtthl**3
        else
            sqrtthl  = 0.0
            skew_thl = 0.
        endif
        if (qwsec > 0.0) then
            sqrtqt   = sqrt(qwsec)
            skew_qw =  qt3/sqrtqt**3
        else
            sqrtqt   = 1e-4*total_water
            skew_qw  = 0.
        endif

        ! Find parameters of the double Gaussian PDF of vertical velocity

        !          aterm = pdf_a

        if (use_aterm_memory/=0) then   ! use memory in aterm and qt skewness
            aterm = pdf_a

            if (mffrc>=1e-3) then                ! if active updraft this timestep
                if (aterm<0.5) then                ! if distribution is skewed (recent updrafts)
                    aterm = max(mffrc,aterm*max(1.-DT/tauskew,0.0))
                else                               ! if distribution unskewed
                    aterm = mffrc
                end if
            else                                 ! if no active updraft
                if (aterm.lt.0.5 .and. aterm.gt.1e-3) then  ! but there is residual skewness
                    aterm = aterm*max(1.-DT/tauskew,0.0)
                else
                    aterm = 0.5
                end if
            end if

        else  ! don't use memory in aterm and qt skewness

            aterm = mffrc
            aterm = max(1e-3,min(0.99,aterm))
            if (mffrc.le.1e-3) aterm = 0.5
        end if

        onema = 1.0 - aterm


        ! If variance of w is too small or no skewness then
        !          IF (w_sec <= w_tol_sqd .or. mffrc.lt.0.01) THEN ! If variance of w is too small then
        IF (w_sec <= w_tol_sqd) THEN ! If variance of w is too small then
            Skew_w = 0.
            w1_1   = 0.
            w1_2   = 0.
            w2_1   = w_sec
            w2_2   = w_sec
            !            aterm  = 0.5
            !            onema  = 0.5
        ELSE

            ! Proportionality coefficients between widths of each vertical velocity
            ! gaussian and the sqrt of the second moment of w
            !           w2_1 = 0.4
            !           w2_2 = 0.4

            ! analytic double gaussian 2, variable sigma_w

            wrk2 = 0.667*abs(Skew_w)**0.333    ! m below A.24
            ! not used     wrk = (1+wrk2*wrk2)**3/((3.+wrk2*wrk2)*wrk2)**2  ! M in A.24

            w2_1 = (onema/(aterm*(1.+wrk2**2)))**0.5
            w2_2 = (aterm/(onema*(1.+wrk2**2)))**0.5

            w1_1 = wrk2*w2_1             ! w1_tilde in A.23
            w1_2 = -wrk2*w2_2

            ! Compute realtive weight of the first PDF "plume"
            ! See Eq A4 in Pete's dissertaion -  Ensure 0.01 < a < 0.99

            !            wrk = 1.0 - w2_1    ! 1-sigw2tilde = 1-0.4
            !            aterm = max(0.01,min(0.5*(1.-Skew_w*sqrt(1./(4.*wrk*wrk*wrk+Skew_w*Skew_w))),0.99))

            !            sqrtw2t = sqrt(wrk)

            ! Eq. A.5-A.6
            !            wrk  =   sqrt(onema/aterm)
            !            w1_1 =   sqrtw2t * wrk  ! w1tilde (A.5)
            !            w1_2 = - sqrtw2t / wrk  ! w2tilde (A.6)

            !            w2_1 = w2_1 * w_sec  ! sigma_w1 **2
            !            w2_2 = w2_2 * w_sec  ! sigma_w2 **2

        ENDIF


        !  Find parameters of the PDF of liquid/ice static energy

        ! inter-gaussian flux limited to 2x total flux
        !          whlntrgs = max(min(whlfac,2.*abs(wthlsec)),-2.*abs(wthlsec))

        IF (thlsec <= thl_tol*thl_tol .or. abs(w1_2-w1_1) <= w_thresh) THEN
            thl1_1     = thl_first
            thl1_2     = thl_first
            thl2_1     = thlsec
            thl2_2     = thlsec
            sqrtthl2_1 = sqrt(thlsec)
            sqrtthl2_2 = sqrtthl2_1

        ELSE

            !            corrtest1 = max(-1.0,min(1.0,whlntrgs/(sqrtw2*sqrtthl)))
            corrtest1 = max(-1.0,min(1.0,wthlsec/(sqrtw2*sqrtthl)))

            thl1_1 = -corrtest1 / w1_2       ! A.7
            thl1_2 = -corrtest1 / w1_1       ! A.8

            !            thl1_1 = -whlntrgs / (w1_2*sqrtthl)   !   normalized
            !            thl1_2 = -whlntrgs / (w1_1*sqrtthl)

            wrk1   = thl1_1 * thl1_1
            wrk2   = thl1_2 * thl1_2
            wrk3   = 1.0 - aterm*wrk1 - onema*wrk2
            wrk4   = skew_thl - aterm*wrk1*thl1_1 - onema*wrk2*thl1_2
            wrk    = 3. * (thl1_2-thl1_1)
            if (wrk /= 0.0) then
                thl2_1 = thlsec * min(100.,max(0.,( 3.*thl1_2*wrk3-wrk4)/(aterm*wrk))) ! A.9
                thl2_2 = thlsec * min(100.,max(0.,(-3.*thl1_1*wrk3+wrk4)/(onema*wrk))) ! A.10
            else
                !              thl2_1 = 0.0
                !              thl2_2 = 0.0
                thl2_1 = thlsec
                thl2_2 = thlsec
            endif

            thl1_1 = thl1_1*sqrtthl + thl_first    ! convert to physical units
            thl1_2 = thl1_2*sqrtthl + thl_first

            sqrtthl2_1 = sqrt(thl2_1)
            sqrtthl2_2 = sqrt(thl2_2)

        ENDIF

        ! implied correlation coefficient
#ifdef PDFDIAG
        PDF_RWTH = max(-1.,min(1.,( wthlsec/sqrtw2-aterm*(thl1_1-thl_first)*(w1_1-w_first) &
            -onema*(thl1_2-thl_first)*(w1_2-w_first) )               &
            / (aterm*sqrt(thl2_1*w2_1)+onema*sqrt(thl2_2*w2_2)) ))
#endif

        !  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO

        ! inter-gaussian flux, limited to 2x total flux
        !          wqtntrgs = max(min(wqtfac,2.*abs(wqwsec)),-2.*abs(wqwsec))

        IF (qwsec <= rt_tol*rt_tol .or. abs(w1_2-w1_1) <= w_thresh) THEN ! if no active updrafts

            if (aterm .lt. 1e-3 .or. aterm.gt.0.499 .or. Skew_qw.eq.0.) then ! if no residual skewness
                qw1_1     = total_water
                qw1_2     = total_water
                qw2_1     = qwsec
                qw2_2     = qwsec
                sqrtqw2_1 = sqrt(qw2_1)
                sqrtqw2_2 = sqrt(qw2_2)
            else
                !              qw1_1     = total_water
                !              qw1_2     = total_water
                !              qw2_1     = qwsec
                !              qw2_2     = qwsec
                wrk1 = min(10.,skew_qw*sqrtqt**3)   ! third moment qt
                qw1_1 = total_water + (wrk1/(2.*aterm-aterm**3/onema**2))**(1./3.)
                qw1_2 = (total_water -aterm*qw1_1)/onema
                qw2_1 = qwsec - min(0.5*qwsec,max(0.,(aterm/onema)*(qw1_1-total_water)**2))
                qw2_2 = qw2_1
                sqrtqw2_1 = sqrt(qw2_1)
                sqrtqw2_2 = sqrt(qw2_2)
            end if

        ELSE  ! active updrafts

            !            corrtest2 = max(-1.0,min(1.0,wqtntrgs/(sqrtw2*sqrtqt)))
            corrtest2 = max(-1.0,min(1.0,0.5*wqwsec/(sqrtw2*sqrtqt)))

            qw1_1 = - corrtest2 / w1_2            ! A.7
            qw1_2 = - corrtest2 / w1_1            ! A.8

            tsign = abs(qw1_2-qw1_1)

            wrk1  = qw1_1 * qw1_1
            wrk2  = qw1_2 * qw1_2
            wrk3  = 1.      - aterm*wrk1       - onema*wrk2
            wrk4  = Skew_qw - aterm*wrk1*qw1_1 - onema*wrk2*qw1_2
            wrk   = 3. * (qw1_2-qw1_1)

            if (wrk /= 0.0) then
                qw2_1 = qwsec * min(100.,max(0.,( 3.*qw1_2*wrk3-wrk4)/(aterm*wrk))) ! A.10
                qw2_2 = qwsec * min(100.,max(0.,(-3.*qw1_1*wrk3+wrk4)/(onema*wrk))) ! A.11
            else
                !              qw2_1 = 0.0
                !              qw2_2 = 0.0
                qw2_1 = qwsec
                qw2_2 = qwsec
            endif

            qw1_1 = qw1_1*sqrtqt + total_water
            qw1_2 = qw1_2*sqrtqt + total_water

            sqrtqw2_1 = sqrt(qw2_1)
            sqrtqw2_2 = sqrt(qw2_2)

        ENDIF   ! if qwsec small

        ! implied correlation coefficient
#ifdef PDFDIAG
        PDF_RWQT = max(-1.,min(1.,( wqwsec/sqrtw2-aterm*(qw1_1-total_water)*(w1_1-w_first) &
            -onema*(qw1_2-total_water)*(w1_2-w_first) )              &
            / (aterm*sqrt(qw2_1*w2_1)+onema*sqrt(qw2_2*w2_2)) ))
#endif

        !  CONVERT FROM TILDA VARIABLES TO "REAL" VARIABLES

        w1_1 = w1_1*sqrtw2 + w_first    ! using A.5 and A.6
        w1_2 = w1_2*sqrtw2 + w_first    ! note: this is already done for w2_x


        !=== Assign PDF diagnostics ===!

        pdf_a = aterm

#ifdef PDFDIAG
        pdf_th1 = thl1_1
        pdf_th2 = thl1_2
        pdf_sigth1 = sqrtthl2_1
        pdf_sigth2 = sqrtthl2_2

        pdf_qt1 = qw1_1
        pdf_qt2 = qw1_2
        pdf_sigqt1 = sqrtqw2_1
        pdf_sigqt2 = sqrtqw2_2

        pdf_w1 = w1_1
        pdf_w2 = w1_2
        if (w2_1.ne.0.) then
            pdf_sigw1 = w2_1*sqrtw2
            pdf_sigw2 = w2_2*sqrtw2
        else
            pdf_sigw1 = 0.0
            pdf_sigw2 = 0.0
        end if
#endif

        !==============================!


        !  FIND WITHIN-PLUME CORRELATIONS

        testvar = aterm*sqrtqw2_1*sqrtthl2_1 + onema*sqrtqw2_2*sqrtthl2_2

        IF (testvar == 0) THEN
            r_qwthl_1 = 0.
        ELSE
            r_qwthl_1 = max(-1.0,min(1.0,(qwthlsec-aterm*(qw1_1-total_water)*(thl1_1-thl_first)-onema*(qw1_2-total_water)*(thl1_2-thl_first))/testvar)) ! A.12
        ENDIF

#ifdef PDFDIAG
        pdf_rqtth = r_qwthl_1
#endif


        !  BEGIN TO COMPUTE CLOUD PROPERTY STATISTICS
        ! This section follows Bogenschutz thesis Appendix A, based on
        ! Sommeria and Deardorff (1977) and Lewellen and Yoh (1993).

        Tl1_1 = thl1_1 - gamaz
        Tl1_2 = thl1_2 - gamaz

        ! Now compute qs

        esval1_1 = 0.
        esval1_2 = 0.
        esval2_1 = 0.
        esval2_2 = 0.
        om1      = 1.
        om2      = 1.

        ! Partition based on temperature for the first plume

        IF (Tl1_1 >= tbgmax) THEN
            esval1_1 = MAPL_EQsat(Tl1_1)
            lstarn1  = lcond
        ELSE IF (Tl1_1 < tbgmin) THEN
            esval1_1 = MAPL_EQsat(Tl1_1,OverIce=.TRUE.)
            lstarn1  = lsub
        ELSE
            esval1_1 = MAPL_EQsat(Tl1_1)
            esval2_1 = MAPL_EQsat(Tl1_1,OverIce=.TRUE.)
            om1      = max(0.,min(1.,a_bg*(Tl1_1-tbgmin)))
            lstarn1  = lcond + (1.-om1)*lfus
        ENDIF

        ! this is qs evaluated at Tl
        qs1   =     om1  * (0.622*esval1_1/max(esval1_1,pval-0.378*esval1_1))      &
            + (1.-om1) * (0.622*esval2_1/max(esval2_1,pval-0.378*esval2_1))

        beta1 = (lstarn1*lstarn1*onebrvcp) / (Tl1_1*Tl1_1)

        ! Are the two plumes equal?  If so then set qs and beta
        ! in each column to each other to save computation
        IF (Tl1_1 == Tl1_2) THEN
            qs2   = qs1
            beta2 = beta1
        ELSE

            IF (Tl1_2 < tbgmin) THEN
                esval1_2 = MAPL_EQsat(Tl1_2,OverIce=.TRUE.)
                lstarn2  = lsub
            ELSE IF (Tl1_2 >= tbgmax) THEN
                esval1_2 = MAPL_EQsat(Tl1_2)
                lstarn2  = lcond
            ELSE
                esval1_2 = MAPL_EQsat(Tl1_2)
                esval2_2 = MAPL_EQsat(Tl1_2,OverIce=.TRUE.)
                om2      = max(0.,min(1.,a_bg*(Tl1_2-tbgmin)))
                lstarn2  = lcond + (1.-om2)*lfus
            ENDIF

            qs2   =     om2  * (0.622*esval1_2/max(esval1_2,pval-0.378*esval1_2))    &
            + (1.-om2) * (0.622*esval2_2/max(esval2_2,pval-0.378*esval2_2))

            beta2 = (lstarn2*lstarn2*onebrvcp) / (Tl1_2*Tl1_2)              ! A.18

        ENDIF


        !  Now compute cloud stuff -  compute s term

        cqt1  = 1.0 / (1.0+beta1*qs1)                                     ! A.19
        wrk   = (1.0+beta1*qw1_1) * cqt1

        s1    = qw1_1 - qs1* wrk                                          ! A.17
        cthl1 = cqt1*wrk*(cp/lcond)*beta1*qs1*pkap                        ! A.20

        wrk1   = cthl1 * cthl1
        wrk2   = cqt1  * cqt1
        std_s1 = sqrt(max(0.,wrk1*thl2_1+wrk2*qw2_1-2.*cthl1*sqrtthl2_1*cqt1*sqrtqw2_1*r_qwthl_1))

        qn1 = 0.
        C1  = 0.

        IF (std_s1 /= 0) THEN
            wrk = s1 / (std_s1*sqrt2)
            C1 = 0.5*(1.+erf(wrk))                                         ! A.15
            IF (C1 /= 0) qn1 = s1*C1 + (std_s1*sqrtpii)*exp(-wrk*wrk)      ! A.16
        ELSEIF (s1 > 0) THEN
            C1  = 1.0
            qn1 = s1
        ENDIF

        ! now compute non-precipitating cloud condensate

        ! If two plumes exactly equal, then just set many of these
        ! variables to themselves to save on computation.
        IF (qw1_1 == qw1_2 .and. thl2_1 == thl2_2 .and. qs1 == qs2) THEN
            s2     = s1
            cthl2  = cthl1
            cqt2   = cqt1
            std_s2 = std_s1
            C2     = C1
            qn2    = qn1
        ELSE

            cqt2   = 1.0 / (1.0+beta2*qs2)
            wrk    = (1.0+beta2*qw1_2) * cqt2
            s2     = qw1_2 - qs2*wrk
            cthl2  = wrk*cqt2*(cp/lcond)*beta2*qs2*pkap
            wrk1   = cthl2 * cthl2
            wrk2   = cqt2  * cqt2
            std_s2 = sqrt(max(0.,wrk1*thl2_2+wrk2*qw2_2-2.*cthl2*sqrtthl2_2*cqt2*sqrtqw2_2*r_qwthl_1))

            qn2 = 0.
            C2  = 0.

            IF (std_s2 /= 0) THEN
                wrk = s2 / (std_s2*sqrt2)
                C2  = 0.5*(1.+erf(wrk))
                IF (C2 /= 0) qn2 = s2*C2 + (std_s2*sqrtpii)*exp(-wrk*wrk)
            ELSEIF (s2 > 0) THEN
                C2  = 1.0
                qn2 = s2
            ENDIF

        ENDIF


        ! finally, compute the SGS cloud fraction
        diag_frac = aterm*C1 + onema*C2

        om1 = max(0.,min(1.,(Tl1_1-tbgmin)*a_bg))
        om2 = max(0.,min(1.,(Tl1_2-tbgmin)*a_bg))

        qn1 = min(qn1,qw1_1)
        qn2 = min(qn2,qw1_2)

        ql1 = qn1*om1
        ql2 = qn2*om2

        qi1 = qn1 - ql1
        qi2 = qn2 - ql2

        diag_qn = min(max(0.0, aterm*qn1 + onema*qn2), total_water)
        diag_ql = min(max(0.0, aterm*ql1 + onema*ql2), diag_qn)
        diag_qi = diag_qn - diag_ql

        !!! temporary
        !          if (abs(qc-diag_qn)>0.001) print *,'SHOC: t=',tabs,' s1=',s1,' qn1=',qn1,' qs1=',qs1,' qt1=',qw1_1


        ! Update temperature variable based on diagnosed cloud properties
        om1         = max(0.,min(1.,(tabs-tbgmin)*a_bg))
        lstarn1     = lcond + (1.-om1)*lfus
        !          tabs = thl_first - gamaz + fac_cond*(diag_ql) &
        !                            + fac_sub *(diag_qi) !&
        !  + tkesbdiss(i,j,k) * (dtn/cp)      ! tke dissipative heating
        ! Update moisture fields



        qc      = diag_ql + diag_qi
        !         qi      = diag_qi
        !         qwv     = total_water - diag_qn
        cld_sgs = diag_frac

        if (sqrtqt>0.0 .AND. sqrtw2>0.0) then
            rwqt = (1.-0.5)*wqwsec/(sqrtqt*sqrtw2)
            !            rwqt = (wqwsec)/(sqrtqt*sqrtw2)
            !            rwqt = max(-1.,min(1.,pdf_rwqt))
        else
            rwqt = 0.0
        end if
        if (sqrtthl>0.0 .AND. sqrtw2>0.0) then
            rwthl = wthlsec/(sqrtthl*sqrtw2)
            !            rwthl = max(-1.,min(1.,pdf_rwth))
        else
            rwthl = 0.0
        end if

        wql1 = C1*(cqt1*sqrt(w2_1)*sqrt(qw2_1)*rwqt-cthl1*sqrt(w2_1)*sqrt(thl2_1)*rwthl)
        wql2 = C2*(cqt2*sqrt(w2_2)*sqrt(qw2_2)*rwqt-cthl2*sqrt(w2_2)*sqrt(thl2_2)*rwthl)


        ! Compute the liquid water flux
        wqls = aterm * ((w1_1-w_first)*ql1+wql1) + onema * ((w1_2-w_first)*ql2+wql2)
        wqis = aterm * ((w1_1-w_first)*qi1) + onema * ((w1_2-w_first)*qi2)

        ! diagnostic buoyancy flux.  Includes effects from liquid water, ice
        ! condensate, liquid & ice precipitation
        wrk = epsv * thv

        bastoeps = (rv/rgas) * thv   ! thetav / epsilon

        wthv_sec = wthlsec + wrk*wqwsec                                     &
            + (fac_cond-bastoeps)*wqls                                 &
            + (fac_sub-bastoeps) *wqis

        !                          + ((lstarn1/cp)-thv(i,j,k))*0.5*(wqp_sec(i,j,kd)+wqp_sec(i,j,ku))

    end subroutine partition_dblgss

    subroutine Bergeron_Partition (           &
        DTIME            , &
        PL               , &
        TE               , &
        QV               , &
        QILS             , &
        QICN             , &
        QLLS             , &
        QLCN             , &     
        CF               , &
        AF               , &
        NL               , &
        NI               , & 
        DQALL            , &
        FQI              , &
        CNVFRC, SRF_TYPE , &
        needs_preexisting )

        real ,  intent(in   )    :: DTIME, PL, TE       !, RHCR
        real ,  intent(inout   )    ::  DQALL 
        real ,  intent(in)    :: QV, QLLS, QLCN, QICN, QILS
        real ,  intent(in)    :: CF, AF, NL, NI
        real, intent (out) :: FQI
        real, intent(in) :: CNVFRC, SRF_TYPE
        logical, intent (in)  :: needs_preexisting
     
        real  :: DC, TEFF,QCm,DEP, &
            QC, QS, RHCR, DQSL, DQSI, QI, TC, &
            DIFF, DENAIR, DENICE, AUX, &
            DCF, QTOT, LHCORR,  QL, DQI, DQL, &
            QVINC, QSLIQ, CFALL,  new_QI, new_QL, &
            QSICE, fQI_0, QS_0, DQS_0, FQA, NIX

        DIFF = 0.0     
        DEP=0.0 
        QI = QILS + QICN !neccesary because NI is for convective and large scale 
        QL = QLLS +QLCN
        QTOT=QI+QL
        FQA = 0.0
        if (QTOT .gt. 0.0) FQA = (QICN+QILS)/QTOT
        NIX= (1.0-FQA)*NI

        DQALL=DQALL/DTIME                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
        CFALL= min(CF+AF, 1.0)
        TC=TE-MAPL_TICE
        fQI_0 = fQI

        !Completelely glaciated cloud:
        if (TE .ge. iT_ICE_MAX) then   !liquid cloud
            FQI   = 0.0
        elseif(TE .le. iT_ICE_ALL) then !ice cloud
            FQI   = 1.0
        else !mixed phase cloud
            FQI   = 0.0
            if (QILS .le. 0.0) then 
                if (needs_preexisting) then
                ! new 0518 this line ensures that only preexisting ice can grow by deposition.
                ! Only works if explicit ice nucleation is available (2 moment muphysics and up)                        
                else
                    fQi  =   ice_fraction( TE, CNVFRC, SRF_TYPE )
                end if                      
                return 
            end if 
        
            QVINC=  QV 
            QSLIQ  = GEOS_QsatLQU( TE, PL*100.0 , DQ=DQSL )
            QSICE  = GEOS_QsatICE( TE, PL*100.0 , DQ=DQSI )
            QVINC =MIN(QVINC, QSLIQ) !limit to below water saturation 

            ! Calculate deposition onto preexisting ice 

            DIFF=(0.211*1013.25/(PL+0.1))*(((TE+0.1)/MAPL_TICE)**1.94)*1e-4  !From Seinfeld and Pandis 2006
            DENAIR=PL*100.0/MAPL_RGAS/TE
            DENICE= 1000.0*(0.9167 - 1.75e-4*TC -5.0e-7*TC*TC) !From PK 97
            LHcorr = ( 1.0 + DQSI*MAPL_ALHS/MAPL_CP) !must be ice deposition

            if  ((NIX .gt. 1.0) .and. (QILS .gt. 1.0e-10)) then 
                DC=max((QILS/(NIX*DENICE*MAPL_PI))**(0.333), 20.0e-6) !Assumme monodisperse size dsitribution 
            else
                DC = 20.0e-6
            end if

            TEFF= NIX*DENAIR*2.0*MAPL_PI*DIFF*DC/LHcorr ! 1/Dep time scale 

            DEP=0.0
            if ((TEFF .gt. 0.0) .and. (QILS .gt. 1.0e-14)) then 
                AUX =max(min(DTIME*TEFF, 20.0), 0.0)
                DEP=(QVINC-QSICE)*(1.0-EXP(-AUX))/DTIME
            end if
            DEP=MAX(DEP, -QILS/DTIME) !only existing ice can be sublimated

            DQI = 0.0
            DQL = 0.0
            FQI=0.0
            !Partition DQALL accounting for Bergeron-Findensen process
            if  (DQALL .ge. 0.0) then !net condensation. Note: do not allow bergeron with QLCN
                if (DEP .gt. 0.0) then 
                    DQI = min(DEP, DQALL + QLLS/DTIME)
                    DQL = DQALL - DQI
                else
                    DQL=DQALL ! could happen because the PDF allows condensation in subsaturated conditions
                    DQI = 0.0 
                end if
            end if
            if  (DQALL .lt. 0.0) then  !net evaporation. Water evaporates first regaardless of DEP   
                DQL = max(DQALL, -QLLS/DTIME)   
                DQI = max(DQALL - DQL, -QILS/DTIME)        
            end if
            if (DQALL .ne. 0.0)  FQI=max(min(DQI/DQALL, 1.0), 0.0)

        end if !=====  

    end subroutine Bergeron_Partition

    function ICE_FRACTION_3D (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
        real, intent(in) :: TEMP(:,:,:),CNV_FRACTION(:,:),SRF_TYPE(:,:)
        real :: ICEFRCT(size(TEMP,1),size(TEMP,2),size(TEMP,3))
        integer :: i,j,l
        do l=1,size(TEMP,3)
            do j=1,size(TEMP,2)
                do i=1,size(TEMP,1)
                    ICEFRCT(i,j,l) = ICE_FRACTION_SC(TEMP(i,j,l),CNV_FRACTION(i,j),SRF_TYPE(i,j))
                enddo
            enddo
        enddo
    end function ICE_FRACTION_3D

    function ICE_FRACTION_2D (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
        real, intent(in) :: TEMP(:,:),CNV_FRACTION(:,:),SRF_TYPE(:,:)
        real :: ICEFRCT(size(TEMP,1),size(TEMP,2))
        integer :: i,j
        do j=1,size(TEMP,2)
            do i=1,size(TEMP,1)
                ICEFRCT(i,j) = ICE_FRACTION_SC(TEMP(i,j),CNV_FRACTION(i,j),SRF_TYPE(i,j))
            enddo
        enddo
    end function ICE_FRACTION_2D

    function ICE_FRACTION_1D (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
        real, intent(in) :: TEMP(:),CNV_FRACTION(:),SRF_TYPE(:)
        real :: ICEFRCT(size(TEMP))
        integer :: i
        do i=1,size(TEMP)
            ICEFRCT(i) = ICE_FRACTION_SC(TEMP(i),CNV_FRACTION(i),SRF_TYPE(i))
        enddo
    end function ICE_FRACTION_1D

    function ICE_FRACTION_SC (TEMP,CNV_FRACTION,SRF_TYPE) RESULT(ICEFRCT)
        real, intent(in) :: TEMP,CNV_FRACTION,SRF_TYPE
        real             :: ICEFRCT
        real             :: tc, ptc
        real             :: ICEFRCT_C, ICEFRCT_M

        ! Anvil clouds
        ! Anvil-Convective sigmoidal function like figure 6(right)
        ! Sigmoidal functions Hu et al 2010, doi:10.1029/2009JD012384
        ICEFRCT_C  = 0.00
        if ( TEMP <= aT_ICE_ALL ) then
            ICEFRCT_C = 1.000
        else if ( (TEMP > aT_ICE_ALL) .AND. (TEMP <= aT_ICE_MAX) ) then
            ICEFRCT_C = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - aT_ICE_ALL ) / ( aT_ICE_MAX - aT_ICE_ALL ) ) )
        end if
        ICEFRCT_C = MIN(ICEFRCT_C,1.00)
        ICEFRCT_C = MAX(ICEFRCT_C,0.00)
        ICEFRCT_C = ICEFRCT_C**aICEFRPWR
#ifdef MODIS_ICE_POLY
        ! Use MODIS polynomial from Hu et al, DOI: (10.1029/2009JD012384) 
        tc = MAX(-46.0,MIN(TEMP-MAPL_TICE,46.0)) ! convert to celcius and limit range from -46:46 C
        ptc = 7.6725 + 1.0118*tc + 0.1422*tc**2 + 0.0106*tc**3 + 0.000339*tc**4 + 0.00000395*tc**5
        ICEFRCT_M = 1.0 - (1.0/(1.0 + exp(-1*ptc)))
#else
        ! Sigmoidal functions like figure 6b/6c of Hu et al 2010, doi:10.1029/2009JD012384
        if (SRF_TYPE == 2.0) then
            ! Over snow/ice
            ICEFRCT_M  = 0.00
            if ( TEMP <= iT_ICE_ALL ) then
                ICEFRCT_M = 1.000
            else if ( (TEMP > iT_ICE_ALL) .AND. (TEMP <= iT_ICE_MAX) ) then
                ICEFRCT_M = 1.00 -  ( TEMP - iT_ICE_ALL ) / ( iT_ICE_MAX - iT_ICE_ALL )
            end if
            ICEFRCT_M = MIN(ICEFRCT_M,1.00)
            ICEFRCT_M = MAX(ICEFRCT_M,0.00)
            ICEFRCT_M = ICEFRCT_M**iICEFRPWR
        else if (SRF_TYPE > 1.0) then
            ! Over Land
            ICEFRCT_M  = 0.00
            if ( TEMP <= lT_ICE_ALL ) then
                ICEFRCT_M = 1.000
            else if ( (TEMP > lT_ICE_ALL) .AND. (TEMP <= lT_ICE_MAX) ) then
                ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - lT_ICE_ALL ) / ( lT_ICE_MAX - lT_ICE_ALL ) ) )
            end if
            ICEFRCT_M = MIN(ICEFRCT_M,1.00)
            ICEFRCT_M = MAX(ICEFRCT_M,0.00)
            ICEFRCT_M = ICEFRCT_M**lICEFRPWR
        else
            ! Over Oceans
            ICEFRCT_M  = 0.00
            if ( TEMP <= oT_ICE_ALL ) then
                ICEFRCT_M = 1.000
            else if ( (TEMP > oT_ICE_ALL) .AND. (TEMP <= oT_ICE_MAX) ) then
                ICEFRCT_M = SIN( 0.5*MAPL_PI*( 1.00 - ( TEMP - oT_ICE_ALL ) / ( oT_ICE_MAX - oT_ICE_ALL ) ) )
            end if
            ICEFRCT_M = MIN(ICEFRCT_M,1.00)
            ICEFRCT_M = MAX(ICEFRCT_M,0.00)
            ICEFRCT_M = ICEFRCT_M**oICEFRPWR
        endif
#endif
        ! Combine the Convective and MODIS functions
        ICEFRCT  = ICEFRCT_M*(1.0-CNV_FRACTION) + ICEFRCT_C*(CNV_FRACTION)
    end function ICE_FRACTION_SC

    subroutine hystpdf( &
        DT          , &
        ALPHA       , &
        PDFSHAPE    , &
        CNVFRC      , &
        SRF_TYPE    , &
        PL          , &
        ZL          , &
        QV          , &
        QCl         , &
        QAl         , &
        QCi         , &
        QAi         , &
        TE          , &
        CF          , &
        AF          , &
        NL          , &
        NI          , &
        WHL         , &
        WQT         , &
        HL2         , &
        QT2         , &
        HLQT        , & 
        W3          , &
        W2          , &
        MFQT3       , &
        MFHL3       , &
        MF_FRC      , &
        PDF_A,      &  ! can remove these after development
        PDFITERS,   &
#ifdef PDFDIAG
        PDF_SIGW1,  &
        PDF_SIGW2,  &
        PDF_W1,     &
        PDF_W2,     &
        PDF_SIGHL1, &
        PDF_SIGHL2, &
        PDF_HL1,    &
        PDF_HL2,    &
        PDF_SIGQT1, &
        PDF_SIGQT2, &
        PDF_QT1,    &
        PDF_QT2,    &
        PDF_RHLQT,  &
        PDF_RWHL,   &
        PDF_RWQT,   &
#endif
        WTHV2,      &
        WQL,        &
        needs_preexisting, &
        USE_BERGERON )

        real, intent(in)    :: DT,ALPHA,PL,ZL
        integer, intent(in) :: PDFSHAPE
        real, intent(inout) :: TE,QV,QCl,QCi,CF,QAl,QAi,AF,PDF_A
        real, intent(in)    :: NL,NI,CNVFRC,SRF_TYPE
        real, intent(in)    :: WHL,WQT,HL2,QT2,HLQT,W3,W2,MF_FRC,MFQT3,MFHL3
#ifdef PDFDIAG
        real, intent(out)   :: PDF_SIGW1, PDF_SIGW2, PDF_W1, PDF_W2, &
                                PDF_SIGHL1, PDF_SIGHL2, PDF_HL1, PDF_HL2, &
                                PDF_SIGQT1, PDF_SIGQT2, PDF_QT1, PDF_QT2, &
                                PDF_RHLQT,  PDF_RWHL, PDF_RWQT
#endif
        real, intent(out)   :: WTHV2, WQL
        real, intent(out)   :: PDFITERS
        logical, intent(in) :: needs_preexisting, USE_BERGERON

        ! internal arrays
        real :: QCO, QVO, CFO, QAO, TAU,HL
        real :: QT, QMX, QMN, DQ, sigmaqt1, sigmaqt2

        real :: TEO,QSx,DQsx,QS,DQs

        real :: TEp, QSp, CFp, QVp, QCp
        real :: TEn, QSn, CFn, QVn, QCn

        real :: QCx, QVx, CFx, QAx, QC, QA, fQi
        real :: dQAi, dQAl, dQCi, dQCl, Nfac, NLv, NIv 

!      real :: fQip

        real :: tmpARR
        real :: ALHX, DQCALL
        ! internal scalars
        integer :: N, nmax

        QC = QCl + QCi
        QA = QAl + QAi
        QT  =  QC  + QA + QV  !Total water after microphysics
        tmpARR = 0.0
        nmax =  20
        QAx = 0.0

        if ( AF < 1.0 )  tmpARR = 1./(1.-AF)

        TEo = TE

        DQSx  = GEOS_DQSAT( TE, PL, QSAT=QSx )
        CFx = CF*tmpARR
        QCx = QC*tmpARR
        QVx = ( QV - QSx*AF )*tmpARR

        if ( AF >= 1.0 )  QVx = QSx*1.e-4 
        if ( AF >  0.0 )  QAx = QA/AF

        QT  = QCx + QVx

        TEp = TEo
        QSn = QSx
        TEn = TEo
        CFn = CFx
        QVn = QVx
        QCn = QCx
        DQS = DQSx

        fQi = ice_fraction( TE, CNVFRC,SRF_TYPE )
        ALHX = (1.0-fQi)*MAPL_ALHL + fQi*MAPL_ALHS

        HL = TEn + (mapl_grav/mapl_cp)*ZL - (ALHX/MAPL_CP)*QCn

        do n=1,nmax

            QVp = QVn
            QCp = QCn
            CFp = CFn
            TEp = TEn

            if(PDFSHAPE.lt.2) then

            sigmaqt1  = ALPHA*QSn
            sigmaqt2  = ALPHA*QSn

            elseif(PDFSHAPE.eq.2) then  ! triangular
                ! for triangular, symmetric: sigmaqt1 = sigmaqt2 = alpha*qsn (alpha is half width)
                ! for triangular, skewed r : sigmaqt1 < sigmaqt2
                ! try: skewed right below 500 mb
                sigmaqt1  = ALPHA*QSn
                sigmaqt2  = ALPHA*QSn
    !         elseif(pdffrac .eq. 3) then ! single gaussian

            elseif(PDFSHAPE .eq. 4) then !lognormal (sigma is dimmensionless)
                sigmaqt1 =  max(ALPHA/sqrt(3.0), 0.001)
            endif

            if (PDFSHAPE.lt.5) then
                call pdffrac(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,CFn)
                call pdfcondensate(PDFSHAPE,qt,sigmaqt1,sigmaqt2,qsn,QCn)
            elseif (PDFSHAPE.eq.5) then

                ! Update the liquid water static energy
                ALHX = (1.0-fQi)*MAPL_ALHL + fQi*MAPL_ALHS
                HL = TEn + (mapl_grav/mapl_cp)*ZL - (ALHX/MAPL_CP)*QCn

                call partition_dblgss(DT/nmax,           &
                                    TEn,          &
                                    QVn,          &
                                    QCn,          &
                                    0.0,          & ! assume OMEGA=0
                                    ZL,           &
                                    PL*100.,      &
                                    QT,           &
                                    HL,          &
                                    WHL,         &
                                    WQT,         &
                                    HL2,         &
                                    QT2,         &
                                    HLQT,        & 
                                    W3,           &
                                    W2,           &
                                    MFQT3,        &
                                    MFHL3,        &
                                    MF_FRC,       &
                                    PDF_A,        &
#ifdef PDFDIAG
                                    PDF_SIGW1,    &
                                    PDF_SIGW2,    &
                                    PDF_W1,       &
                                    PDF_W2,       &
                                    PDF_SIGHL1,   &
                                    PDF_SIGHL2,   &
                                    PDF_HL1,      &
                                    PDF_HL2,      &
                                    PDF_SIGQT1,   &
                                    PDF_SIGQT2,   &
                                    PDF_QT1,      &
                                    PDF_QT2,      &
                                    PDF_RHLQT,    &
                                    PDF_RWHL,     &
                                    PDF_RWQT,     &
#endif
                                    WTHV2,        &
                                    WQL,          &
                                    CFn)
            endif

            IF(USE_BERGERON) THEN
                DQCALL = QCn - QCp
                CF  = CFn * ( 1.-AF)
                Nfac = 100.*PL*R_AIR/TEp !density times conversion factor
                NLv = NL/Nfac
                NIv = NI/Nfac
                call Bergeron_Partition( &         !Microphysically-based partitions the new condensate
                        DT               , &
                        PL               , &
                        TEp              , &
                        QT               , &
                        QCi              , &
                        QAi              , &
                        QCl              , &
                        QAl              , &
                        CF               , &
                        AF               , &
                        NLv              , &
                        NIv              , &
                        DQCALL           , &
                        fQi              , & 
                        CNVFRC,SRF_TYPE  , &
                        needs_preexisting)
            ELSE
                fQi = ice_fraction( TEp, CNVFRC,SRF_TYPE )
            ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! These lines represent adjustments
        ! to anvil condensate due to the 
        ! assumption of a stationary TOTAL 
        ! water PDF subject to a varying 
        ! QSAT value during the iteration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
            if ( AF > 0. ) then
                QAo = QAx  ! + QSx - QS 
            else
                QAo = 0.
            end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ALHX = (1.0-fQi)*MAPL_ALHL + fQi*MAPL_ALHS

            if(PDFSHAPE.eq.1) then 
                QCn = QCp + ( QCn - QCp ) / ( 1. - (CFn * (ALPHA-1.) - (QCn/QSn))*DQS*ALHX/MAPL_CP)             
            elseif(PDFSHAPE.eq.2 .or. PDFSHAPE.eq.5) then
                ! This next line needs correcting - need proper d(del qc)/dT derivative for triangular
                ! for now, just use relaxation of 1/2.
                if (n.ne.nmax) QCn = QCp + ( QCn - QCp ) *0.5
            endif

            QVn = QVp - (QCn - QCp)
            TEn = TEp + (1.0-fQi)*(MAPL_ALHL/MAPL_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF ) &
                +      fQi* (MAPL_ALHS/MAPL_CP)*( (QCn - QCp)*(1.-AF) + (QAo-QAx)*AF )

            PDFITERS = n
            if (abs(Ten - Tep) .lt. 0.00001) exit 

            DQS  = GEOS_DQSAT( TEn, PL, QSAT=QSn )

        enddo ! qsat iteration

        CFo = CFn
        CF = CFn
        QCo = QCn
        QVo = QVn
        TEo = TEn

        ! Update prognostic variables.  Deal with special case of AF=1
        ! Temporary variables QCo, QAo become updated grid means.
        if ( AF < 1.0 ) then
            CF  = CFo * ( 1.-AF)
            QCo = QCo * ( 1.-AF)
            QAo = QAo *   AF  
        else

            ! Special case AF=1, i.e., box filled with anvil. 
            !   - Note: no guarantee QV_box > QS_box
            CF  = 0.          ! Remove any other cloud
            QAo = QA  + QC    ! Add any LS condensate to anvil type
            QCo = 0.          ! Remove same from LS   
            QT  = QAo + QV    ! Total water
            ! Now set anvil condensate to any excess of total water 
            ! over QSx (saturation value at top)
            QAo = MAX( QT - QSx, 0. )
        end if

        ! Now take {\em New} condensate and partition into ice and liquid
        ! taking care to keep both >=0 separately. New condensate can be
        ! less than old, so $\Delta$ can be < 0.

        dQCl = 0.0
        dQCi = 0.0
        dQAl = 0.0
        dQAi = 0.0

        !large scale   

        QCx   = QCo - QC
        if  (QCx .lt. 0.0) then  !net evaporation. Water evaporates first
            dQCl = max(QCx, -QCl)   
            dQCi = max(QCx - dQCl, -QCi)
        else
            dQCl  = (1.0-fQi)*QCx
            dQCi  =    fQi  * QCx
        end if

        !Anvil   
        QAx   = QAo - QA

        if  (QAx .lt. 0.0) then  !net evaporation. Water evaporates first
            dQAl = max(QAx, -QAl)   
            dQAi = max(QAx - dQAl, -QAi)
        else            
            dQAl  =  (1.0-fQi)*QAx
            dQAi  = QAx*fQi
        end if

        ! Clean-up cloud if fractions are too small
        if ( AF < 1.e-5 ) then
            dQAi = -QAi
            dQAl = -QAl
        end if
        if ( CF < 1.e-5 ) then
            dQCi = -QCi
            dQCl = -QCl
        end if

        QAi    = QAi + dQAi
        QAl    = QAl + dQAl
        QCi    = QCi + dQCi
        QCl    = QCl + dQCl
        QV     = QV  - ( dQAi+dQCi+dQAl+dQCl) 

        TE  = TE + (MAPL_ALHL*( dQAi+dQCi+dQAl+dQCl)+MAPL_ALHF*(dQAi+dQCi))/ MAPL_CP

        ! We need to take care of situations where QS moves past QA
        ! during QSAT iteration. This should be only when QA/AF is small
        ! to begin with. Effect is to make QAo negative. So, we 
        ! "evaporate" offending QA's
        !
        ! We get rid of anvil fraction also, although strictly
        ! speaking, PDF-wise, we should not do this.
        if ( QAo <= 0. ) then
            QV  = QV + QAi + QAl
            TE  = TE - (MAPL_ALHS/MAPL_CP)*QAi - (MAPL_ALHL/MAPL_CP)*QAl
            QAi = 0.
            QAl = 0.
            AF  = 0.  
        end if

    end subroutine hystpdf

    subroutine FILLQ2ZERO( Q, MASS, FILLQ  )

        ! New algorithm to fill the negative q values in a mass conserving way.
        ! Conservation of TPW was checked. Donifan Barahona
        ! Updated from FILLQ2ZERO, avoids the usage of scalars

        real, dimension(:,:,:),   intent(inout)  :: Q
        real, dimension(:,:,:),   intent(in)     :: MASS
        real, dimension(:,:),     intent(  out)  :: FILLQ
        real, dimension(:,:), allocatable        :: TPW1, TPW2, TPWC
        integer                                  :: IM,JM,LM, l

        IM = SIZE( Q, 1 )
        JM = SIZE( Q, 2 )
        LM = SIZE( Q, 3 )

        ALLOCATE(TPW1(IM, JM))
        ALLOCATE(TPW2(IM, JM))
        ALLOCATE(TPWC(IM, JM))

        TPW2 =0.0
        TPWC= 0.0
        TPW1 = SUM( Q*MASS, 3 )

        WHERE (Q < 0.0)
        Q=0.0
        END WHERE

        TPW2 = SUM( Q*MASS, 3 )

        WHERE (TPW2 > 0.0)
        TPWC=(TPW2-TPW1)/TPW2
        END WHERE

      !   print*,'TPWC = ', TPWC

        do l=1,LM
        Q(:, :, l)= Q(:, :, l)*(1.0-TPWC) !reduce Q proportionally to the increase in TPW
        end do

        FILLQ = TPW2-TPW1

        DEALLOCATE(TPW1)
        DEALLOCATE(TPW2)
        DEALLOCATE(TPWC)
    end subroutine FILLQ2ZERO

    subroutine EVAP3(&
        DT      , &
        A_EFF   , &
        RHCR    , &
        PL      , &
        TE      , &
        QV      , &
        QL      , &
        QI      , &
        F       , &
        NL      , &
        NI      , &
        QS        )

        real, intent(in   ) :: DT 
        real, intent(in   ) :: A_EFF
        real, intent(in   ) :: RHCR
        real, intent(in   ) :: PL
        real, intent(inout) :: TE
        real, intent(inout) :: QV
        real, intent(inout) :: QL,QI
        real, intent(inout) :: F
        real, intent(in   ) :: NL,NI
        real, intent(in   ) :: QS

        real :: ES,RADIUS,K1,K2,QCm,EVAP,RHx,QC  !,QS

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!         EVAPORATION OF CLOUD WATER.             !!
        !!                                                 !!
        !!  DelGenio et al (1996, J. Clim., 9, 270-303)    !!
        !!  formulation  (Eq.s 15-17)                      !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !   QS  = QSAT(         &
        !               TE    , &
        !               PL      )
     
        ES = 100.* PL * QS  / ( (EPSILON) + (1.0-(EPSILON))*QS )  ! (100's <-^ convert from mbar to Pa)

        RHx = MIN( QV/QS , 1.00 )

        K1 = (MAPL_ALHL**2) * RHO_W / ( K_COND*MAPL_RVAP*(TE**2))

        K2 = MAPL_RVAP * TE * RHO_W / ( DIFFU * (1000./PL) * ES )

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Here DIFFU is given for 1000 mb  !!
        !! so 1000./PR accounts for inc-    !!
        !! reased diffusivity at lower      !!
        !! pressure.                        !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

        if ( ( F > 0.) .and. ( QL > 0. ) ) then
            QCm=QL/F
        else
            QCm=0.
        end if

        RADIUS = LDRADIUS4(PL,TE,QCm,NL,NI,1)
        
        if ( (RHx < RHCR ) .and. (RADIUS > 0.0) ) then
            EVAP = A_EFF*QL*DT*(RHCR - RHx) / ((K1+K2)*RADIUS**2)  ! / (1.00 - RHx)
            EVAP = MIN( EVAP , QL  )
        else
            EVAP = 0.0
        end if

        QC=QL+QI
        if (QC > 0.) then
            F = F * ( QC - EVAP ) / QC
        end if

        QV = QV + EVAP
        QL = QL - EVAP
        TE = TE - alhlbcp*EVAP

    end subroutine EVAP3

    subroutine SUBL3( &
        DT        , &
        A_EFF     , &
        RHCR      , &
        PL        , &
        TE        , &
        QV        , &
        QL        , &
        QI        , &
        F         , &
        NL        , &
        NI        , &
        QS        )

        real, intent(in   ) :: DT
        real, intent(in   ) :: A_EFF
        real, intent(in   ) :: RHCR
        real, intent(in   ) :: PL
        real, intent(inout) :: TE
        real, intent(inout) :: QV
        real, intent(inout) :: QL,QI
        real, intent(inout) :: F
        real, intent(in   ) :: NL,NI
        real, intent(in   ) :: QS

        real :: ES,RADIUS,K1,K2,TEFF,QCm,SUBL,RHx,QC !, QS

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!         SUBLORATION OF CLOUD WATER.             !!
        !!                                                 !!
        !!  DelGenio et al (1996, J. Clim., 9, 270-303)    !!
        !!  formulation  (Eq.s 15-17)                      !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        !   QS  = QSAT(         &
        !               TE    , &
        !               PL      )

        ES = 100.* PL * QS  / ( (EPSILON) + (1.0-(EPSILON))*QS )  ! (100s <-^ convert from mbar to Pa)

        RHx = MIN( QV/QS , 1.00 )

        K1 = (MAPL_ALHL**2) * RHO_I / ( K_COND*MAPL_RVAP*(TE**2))

        K2 = MAPL_RVAP * TE * RHO_I / ( DIFFU * (1000./PL) * ES )

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Here DIFFU is given for 1000 mb  !!
        !! so 1000./PR accounts for inc-    !!
        !! reased diffusivity at lower      !!
        !! pressure.                        !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

        if ( ( F > 0.) .and. ( QI > 0. ) ) then
            QCm=QI/F
        else
            QCm=0.
        end if

        RADIUS = LDRADIUS4(PL,TE,QCm,NL,NI,2)
        
        if ( (RHx < RHCR) .and.(RADIUS > 0.0) ) then
            SUBL = A_EFF*QI*DT*(RHCR - RHx) / ((K1+K2)*RADIUS**2)  ! / (1.00 - RHx)
            SUBL = MIN( SUBL , QI  )
        else
            SUBL = 0.0
        end if

        QC=QL+QI
        if (QC > 0.) then
            F = F * ( QC - SUBL ) / QC
        end if

        QV = QV + SUBL
        QI = QI - SUBL
        TE = TE - alhsbcp*SUBL

    end subroutine SUBL3

    function LDRADIUS4(PL,TE,QC,NNL,NNI,ITYPE) RESULT(RADIUS)

        REAL   , INTENT(IN) :: TE,PL,QC
        REAL   , INTENT(IN) :: NNL,NNI ! #/m^3
        INTEGER, INTENT(IN) :: ITYPE
        REAL  :: RADIUS
        INTEGER, PARAMETER  :: LIQUID=1, ICE=2
        REAL :: NNX,RHO,BB,WC
 
        !- air density (kg/m^3)
        RHO = 100.*PL / (MAPL_RGAS*TE )
        IF(ITYPE == LIQUID) THEN
 
            !- liquid cloud effective radius ----- 
            !- [liu&daum, 2000 and 2005. liu et al 2008]
            !- liquid water content
            WC = RHO*QC*1000. !g/m3
            !- cloud drop number concentration #/m3
            !- from the aerosol model + ....
            NNX = NNL*1.e-6
            !- radius in meters
            RADIUS = MIN(60.e-6,MAX(2.5e-6, 1.e-6*bx*(WC/NNX)**be))
 
        ELSEIF(ITYPE == ICE) THEN
 
            !- ice cloud effective radius ----- 
            !- ice water content
            WC = RHO*QC*1000.  !g/m3
            !------ice cloud effective radius ----- [klaus wyser, 1998]
            if(TE>MAPL_TICE .or. QC <=0.) then
                BB = -2.
            else
                BB = -2. + log10(WC/50.)*(1.e-3*(MAPL_TICE-TE)**1.5)
            endif
            BB     = MIN((MAX(BB,-6.)),-2.)
            RADIUS = 377.4 + 203.3 * BB+ 37.91 * BB **2 + 2.3696 * BB **3
            RADIUS = RADIUS * 1.e-6 !- convert to meter
 
        ELSE
            STOP "WRONG HYDROMETEOR type: CLOUD = 1 OR ICE = 2"
        ENDIF
    end function LDRADIUS4

    subroutine  FIX_UP_CLOUDS( &
        QV, &
        TE, &
        QLC,&
        QIC,&
        CF, &
        QLA,&
        QIA,&
        AF  )

        real, intent(inout) :: TE,QV,QLC,CF,QLA,AF,QIC,QIA

        ! Fix if Anvil cloud fraction too small
        if (AF < 1.E-5) then
            QV  = QV + QLA + QIA
            TE  = TE - (alhlbcp)*QLA - (alhsbcp)*QIA
            AF  = 0.
            QLA = 0.
            QIA = 0.
        end if

        ! Fix if LS cloud fraction too small
        if ( CF < 1.E-5 ) then
            QV = QV + QLC + QIC
            TE = TE - (alhlbcp)*QLC - (alhsbcp)*QIC
            CF  = 0.
            QLC = 0.
            QIC = 0.
        end if
        
        ! LS LIQUID too small
        if ( QLC  < 1.E-8 ) then
            QV = QV + QLC 
            TE = TE - (alhlbcp)*QLC
            QLC = 0.
        end if
        ! LS ICE too small
        if ( QIC  < 1.E-8 ) then
            QV = QV + QIC 
            TE = TE - (alhsbcp)*QIC
            QIC = 0.
        end if

        ! Anvil LIQUID too small
        if ( QLA  < 1.E-8 ) then
            QV = QV + QLA 
            TE = TE - (alhlbcp)*QLA
            QLA = 0.
        end if
        ! Anvil ICE too small
        if ( QIA  < 1.E-8 ) then
            QV = QV + QIA 
            TE = TE - (alhsbcp)*QIA
            QIA = 0.
        end if

        ! Fix ALL cloud quants if Anvil cloud LIQUID+ICE too small
        if ( ( QLA + QIA ) < 1.E-8 ) then
            QV = QV + QLA + QIA
            TE = TE - (alhlbcp)*QLA - (alhsbcp)*QIA
            AF  = 0.
            QLA = 0.
            QIA = 0.
        end if
        ! Ditto if LS cloud LIQUID+ICE too small
        if ( ( QLC + QIC ) < 1.E-8 ) then
            QV = QV + QLC + QIC
            TE = TE - (alhlbcp)*QLC - (alhsbcp)*QIC
            CF  = 0.
            QLC = 0.
            QIC = 0.
        end if

    end subroutine FIX_UP_CLOUDS

    subroutine RADCOUPLE(  &
        TE,              & 
        PL,              & 
        CF,              & 
        AF,              & 
        QV,              &
        QClLS,           & 
        QCiLS,           & 
        QClAN,           & 
        QCiAN,           & 
        QRN_ALL,         & 
        QSN_ALL,         & 
        QGR_ALL,         &
        NL,              &
        NI,              &
        RAD_QV,          &
        RAD_QL,          &  
        RAD_QI,          & 
        RAD_QR,          & 
        RAD_QS,          & 
        RAD_QG,          &
        RAD_CF,          & 
        RAD_RL,          & 
        RAD_RI,          & 
        FAC_RL, MIN_RL, MAX_RL, &
        FAC_RI, MIN_RI, MAX_RI)

        real, intent(in ) :: TE
        real, intent(in ) :: PL
        real, intent(in ) :: AF,CF, QV, QClAN, QCiAN, QClLS, QCiLS
        real, intent(in ) :: QRN_ALL, QSN_ALL, QGR_ALL
        real, intent(in ) :: NL,NI
        real, intent(out) :: RAD_QV,RAD_QL,RAD_QI,RAD_QR,RAD_QS,RAD_QG,RAD_CF,RAD_RL,RAD_RI
        real, intent(in ) :: FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI

        ! Limits on Radii needed to ensure
        ! correct behavior of cloud optical
        ! properties currently calculated in 
        ! sorad and irrad (1e-6 m = micron)

        ! water vapor
        RAD_QV = QV

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Total cloud fraction
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        RAD_CF = MAX(MIN(CF+AF,1.0),0.0)
        if ( RAD_CF >= 1.e-5 ) then
        ! Total In-cloud liquid
        if ( (QClLS + QClAN) >= 1.e-8 ) then
            RAD_QL = ( QClLS + QClAN ) / RAD_CF
        else
            RAD_QL = 0.0
        end if
        ! Total In-cloud ice
        if ( (QCiLS + QCiAN) >= 1.e-8 ) then
            RAD_QI = ( QCiLS + QCiAN ) / RAD_CF
        else
            RAD_QI = 0.0
        end if
        ! Total In-cloud precipitation
        if (QRN_ALL >= 1.e-8 ) then
            RAD_QR = ( QRN_ALL ) / RAD_CF
        else
            RAD_QR = 0.0
        end if
        if (QSN_ALL >= 1.e-8 ) then
            RAD_QS = ( QSN_ALL ) / RAD_CF
        else
            RAD_QS = 0.0
        end if
        if (QGR_ALL >= 1.e-8 ) then
            RAD_QG = ( QGR_ALL ) / RAD_CF
        else
            RAD_QG = 0.0
        end if
        else
        RAD_CF = 0.0
        RAD_QL = 0.0
        RAD_QI = 0.0
        RAD_QR = 0.0
        RAD_QS = 0.0
        RAD_QG = 0.0
        end if
        ! cap the high end of condensates
        RAD_QL = MIN( RAD_QL, 0.01 )
        RAD_QI = MIN( RAD_QI, 0.01 )
        RAD_QR = MIN( RAD_QR, 0.01 )
        RAD_QS = MIN( RAD_QS, 0.01 )
        RAD_QG = MIN( RAD_QG, 0.01 )

        ! LIQUID RADII
        !-BRAMS formulation     
        RAD_RL = LDRADIUS4(PL,TE,RAD_QL,NL,NI,1)
        ! apply limits
        RAD_RL = MAX( MIN_RL, MIN(RAD_RL*FAC_RL, MAX_RL) )

    ! ICE RADII
        !-BRAMS formulation  
        RAD_RI = LDRADIUS4(PL,TE,RAD_QI,NL,NI,2)
    ! apply limits
        RAD_RI = MAX( MIN_RI, MIN(RAD_RI*FAC_RI, MAX_RI) )

    end subroutine RADCOUPLE
end module

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”
  
! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”
  
! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at
  
! http://www.apache.org/licenses/LICENSE-2.0
  
! Unless required by applicable law or agreed to in writing, software distributed
! under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
! CONDITIONS OF ANY KIND, either express or implied. See the License for the
! specific language governing permissions and limitations under the License.