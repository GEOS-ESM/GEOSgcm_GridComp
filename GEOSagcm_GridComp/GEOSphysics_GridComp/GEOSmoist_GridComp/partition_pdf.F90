module partition_pdf

! Total water partitioning based on SHOC double gaussian PDF
!
! Implemented in GEOS by N Arnold, based on SHOC assumed_pdf code
!

 use MAPL_ConstantsMod, only: ggr    => MAPL_GRAV,   &
                              cp     => MAPL_CP,     &
                              rgas   => MAPL_RGAS,   &
                              rv     => MAPL_RVAP,   &
                              lcond  => MAPL_ALHL,   &
                              lfus   => MAPL_ALHS,   &
                              pi     => MAPL_PI,     &
                              MAPL_H2OMW, MAPL_AIRMW

 use MAPL_Mod,          only: MAPL_UNDEF

 use MAPL_SatVaporMod,  only: MAPL_EQsat 

 implicit none

 private

 public partition_dblgss

 contains

!------------------------------------------------------------------------!
!                   Subroutine partition_dblgss                          !
!------------------------------------------------------------------------!
! Compute SGS cloud fraction and SGS condensation                        !
! using assumed analytic DOUBLE-GAUSSIAN PDF for SGS vertical velocity,  !
! moisture, and  liquid/ice water static energy, based on the            !
! general approach of  Larson et al 2002, JAS, 59, 3519-3539,            !
! and Golaz et al 2002, JAS, 59, 3540-3551                               !
! References in the comments in this code are given to                   !
! the Appendix A of Pete Bogenschutz's dissertation.                     !
!------------------------------------------------------------------------!
 subroutine partition_dblgss( dt,           &  ! IN
                              tabs,         &  ! INOUT
                              qwv,          &
                              qc,           &
!                              qi,           &
                              omega,        &  ! IN
                              zl,           &
                              pval,         &
!                              qpl,          &
!                              qpi,          &
                              total_water,  &
                              thl_first,    &
                              wthlsec,      &
                              wqwsec,       &
                              thlsec,       &
                              qwsec,        &
                              qwthlsec,     & 
                              w3var,        &   
                              w_sec,        & 
                              qt3,        &
                              mffrc,        &
                              skew_qw,      &  ! INOUT                      
                              PDF_A,        &
                              PDF_SIGW,     &  ! OUT
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
                              wthv_sec,     &
                              wqls,         &
                              cld_sgs)

   real, intent(   in)  :: DT          ! timestep [s]
!   real, intent(   in)  :: tke         ! turb kinetic energy [m2 s-2]
   real, intent(inout)  :: tabs        ! absolute temperature [K]
   real, intent(in   )  :: qwv         ! specific humidity [kg kg-1]
   real, intent(inout)  :: qc          ! liquid+ice condensate [kg kg-1]
   real, intent(   in)  :: omega       ! resolved pressure velocity
   real, intent(   in)  :: zl          ! layer heights [m]
   real, intent(   in)  :: pval        ! layer pressure [Pa]
!   real, intent(   in)  :: qpl         ! liquid precipitation [kg kg-1]
!   real, intent(   in)  :: qpi         ! ice precipitation  [kg kg-1]
   real, intent(   in)  :: total_water ! total water [kg kg-1]
   real, intent(   in)  :: thl_first   ! liquid water potential temperature [K]
   real, intent(   in)  :: wthlsec     ! thl flux [K m s-1]
   real, intent(   in)  :: wqwsec      ! total water flux [kg kg-1 m s-1]
   real, intent(   in)  :: thlsec     
   real, intent(   in)  :: qwsec 
   real, intent(   in)  :: qwthlsec 
   real, intent(   in)  :: w3var       ! 3rd moment vertical velocity [m3 s-3]
   real, intent(   in)  :: qt3       ! 3rd moment qt from mass flux
   real, intent(   in)  :: w_sec       ! 2nd moment vertical velocity [m2 s-2]
   real, intent(   in)  :: mffrc       ! total EDMF updraft fraction
!   real, intent(inout)  :: qi         ! ice condensate [kg kg-1]
   real, intent(  out)  :: cld_sgs     ! cloud fraction
   real, intent(inout)  ::    skew_qw,      & ! skewness of total water
                              PDF_A           ! fractional area of 1st gaussian
   real, intent(  out)  ::    PDF_SIGW,     & ! std dev w [m s-1]
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
                              PDF_RQTTH       ! QT-TH correlation coeff
   real, intent(  out)  :: wthv_sec
   real, intent(  out)  :: wqls


! Local variables

   integer i,j,k,ku,kd
   real wrk, wrk1, wrk2, wrk3, wrk4, bastoeps
   real gamaz, thv, rwqt, rwthl, wql1, wql2
   real pkap, diag_qn, diag_frac, diag_ql, diag_qi,w_first,                     &
        sqrtw2, sqrtthl, sqrtqt, w1_1, w1_2, w2_1, w2_2, thl1_1, thl1_2,        &
        thl2_1, thl2_2, qw1_1, qw1_2, qw2_1, qw2_2, aterm, onema, sm,           &
        km1, skew_w, cond_w, sqrtw2t,                                           &
        sqrtthl2_1, sqrtthl2_2, sqrtqw2_1, sqrtqw2_2, corrtest1, corrtest2,     &
        tsign, testvar, r_qwthl_1, Tl1_1, Tl1_2, esval1_1, esval1_2, esval2_1,  &
        esval2_2, om1, om2, lstarn1, lstarn2, qs1, qs2, beta1, beta2, cqt1,     &
        cqt2, s1, s2, cthl1, cthl2, std_s1, std_s2, qn1, qn2, C1, C2, ql1, ql2, &
        qi1, qi2, wqis


! Set constants and parameters
   real, parameter :: sqrt2 = sqrt(2.0)
   real, parameter :: sqrtpii = 1.0/sqrt(pi+pi)
   real, parameter :: tbgmin = 233.16
   real, parameter :: tbgmax = 273.16
   real, parameter :: a_bg   = 1.0/(tbgmax-tbgmin)
   real, parameter :: thl_tol = 1.e-2
   real, parameter :: w_thresh = 0.0
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


! define conserved variables
gamaz = gocp * zl
!hl = tabs + gamaz - fac_cond*(qc+qpl) - fac_fus*(qi+qpi)
!total_water = qwv + qc + qi
thv   = tabs * (1.0+epsv*qwv)
thv   = thv*(100000.0/pval) ** kapa

w_first = - rog * omega * thv / pval

!w_sec = 0.6667*tke + mfw2

! Initialize cloud variables to zero  
          diag_qn   = 0.0
          diag_frac = 0.0
          diag_ql   = 0.0
          diag_qi   = 0.0

          pkap = (pval/100000.0) ** kapa
             
              
! Compute square roots of some variables so we don't have to do it again
          if (w_sec > 0.0) then
            sqrtw2   = sqrt(w_sec)
          else
            sqrtw2   = 0.0
          endif
          if (thlsec > 0.0) then
            sqrtthl  = sqrt(thlsec)
          else
            sqrtthl  = 0.0
          endif
          if (qwsec > 0.0) then
            sqrtqt   = sqrt(qwsec)
          else
            sqrtqt   = 0.0
          endif
             

! Find parameters of the double Gaussian PDF of vertical velocity

! Skewness of vertical velocity

         Skew_w = w3var / (sqrtw2*sqrtw2*sqrtw2)

         if (use_aterm_memory/=0) then   ! use memory in aterm and qt skewness
          aterm = pdf_a
!          if (pval>90000.) print *,'before skew_qw=',skew_qw,'  aterm=',aterm,' mffrc=',mffrc

          if (mffrc>=0.01) then                ! if active updraft this timestep
            if (aterm<0.5) then                ! if distribution is skewed (recent updrafts) 
!              aterm = aterm*(1.-DT/900.)
!              skew_qw = 1.0 + (skew_qw-1.)*(1.-DT/900.) 
!              skew_qw = (aterm*skew_qw + mffrc*skew_facw*Skew_w)/(aterm+mffrc)
!              aterm = aterm + mffrc 
              aterm = mffrc  
              skew_qw = skew_facw*Skew_w
            else                               ! if distribution unskewed
              aterm = mffrc
              skew_qw = skew_facw*Skew_w
            end if
          else                                 ! if no active updraft
            if (aterm/=0.5) then               ! but there is residual skewness
              aterm = aterm*max(1.-DT/1200.,0.0)
              skew_qw = skew_qw*max(1.-DT/1200.,0.0)
            else
!              aterm = aterm*max(1.-DT/1200.,0.0)
              skew_qw = skew_qw*max(1.-DT/1200.,0.0)
!              skew_qw = skew_facw*Skew_w
            end if
          end if

         else  ! don't use memory in aterm and qt skewness
 
           aterm = mffrc
           aterm = max(0.01,min(0.99,aterm))

! two options for qt skewness:
           skew_qw = skew_facw*skew_w
!           skew_qw = qt3/sqrtqt**3
         end if
         onema = 1.0 - aterm

!        if (pval>90000.) print *,'after skew_qw=',skew_qw,'  aterm=',aterm

!          if (sqrtqt>1e-5) then
!            skew_qw = 1.0+qt3/sqrtqt**3
!          else
!            skew_qw = 1.0
!          end if


          IF (w_sec <= w_tol_sqd .or. aterm<=0.01 .or. aterm>0.499) THEN ! If variance of w is too small then
!          IF (w_sec <= w_tol_sqd ) THEN ! If variance of w is too small then
                                              ! PDF is a sum of two delta functions
            Skew_w = 0.
            w1_1   = 0.  ! normalized value of sqrt(w2)
            w1_2   = 0.
            w2_1   = w_sec
            w2_2   = w_sec
            aterm  = 0.5
            onema  = 0.5
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

            w2_1 = w2_1 * w_sec  ! sigma_w1 **2
            w2_2 = w2_2 * w_sec  ! sigma_w2 **2

          ENDIF


!  Find parameters of the  PDF of liquid/ice static energy

          IF (thlsec <= thl_tol*thl_tol .or. abs(w1_2-w1_1) <= w_thresh) THEN
            thl1_1     = thl_first
            thl1_2     = thl_first
            thl2_1     = thlsec
            thl2_2     = thlsec
            sqrtthl2_1 = sqrt(thlsec)
            sqrtthl2_2 = sqrtthl2_1

          ELSE

            corrtest1 = max(-1.0,min(1.0,wthlsec/(sqrtw2*sqrtthl)))

            thl1_1 = -corrtest1 / w1_2                 ! A.7
            thl1_2 = -corrtest1 / w1_1                 ! A.8
                
            wrk1   = thl1_1 * thl1_1
            wrk2   = thl1_2 * thl1_2
            wrk3   = 1.0 - aterm*wrk1        - onema*wrk2
            wrk4   = -skew_fact*Skew_w - aterm*wrk1*thl1_1 - onema*wrk2*thl1_2  ! testing - Moorthi
!           wrk4   =     - aterm*wrk1*thl1_1 - onema*wrk2*thl1_2
            wrk    = 3. * (thl1_2-thl1_1)
            if (wrk /= 0.0) then
              thl2_1 = thlsec * min(100.,max(0.,( 3.*thl1_2*wrk3-wrk4)/(aterm*wrk))) ! A.10
              thl2_2 = thlsec * min(100.,max(0.,(-3.*thl1_1*wrk3+wrk4)/(onema*wrk))) ! A.11
            else
              thl2_1 = 0.0
              thl2_2 = 0.0
            endif
!
            thl1_1 = thl1_1*sqrtthl + thl_first
            thl1_2 = thl1_2*sqrtthl + thl_first

            sqrtthl2_1 = sqrt(thl2_1)
            sqrtthl2_2 = sqrt(thl2_2)

          ENDIF



!  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO

          IF (qwsec <= rt_tol*rt_tol .or. abs(w1_2-w1_1) <= w_thresh) THEN
            qw1_1     = total_water !qw_first
            qw1_2     = total_water !qw_first
            qw2_1     = qwsec
            qw2_2     = qwsec
!            sqrtqw2_1 = 0.
!            sqrtqw2_2 = 0.
            sqrtqw2_1 = sqrt(qw2_1)
            sqrtqw2_2 = sqrt(qw2_2)

          ELSE

            corrtest2 = max(-1.0,min(1.0,wqwsec/(sqrtw2*sqrtqt)))

            qw1_1 = - corrtest2 / w1_2            ! A.7
            qw1_2 = - corrtest2 / w1_1            ! A.8

            tsign = abs(qw1_2-qw1_1)


!            skew_qw = skew_facw*Skew_w

!           IF (tsign > 0.4) THEN
!             Skew_qw = skew_facw*Skew_w
!           ELSE IF (tsign <= 0.2) THEN
!             Skew_qw = 0.
!           ELSE
!             Skew_qw = (skew_facw/0.2) * Skew_w * (tsign-0.2)
!           ENDIF

            wrk1  = qw1_1 * qw1_1
            wrk2  = qw1_2 * qw1_2
            wrk3  = 1.      - aterm*wrk1       - onema*wrk2
            wrk4  = Skew_qw - aterm*wrk1*qw1_1 - onema*wrk2*qw1_2
            wrk   = 3. * (qw1_2-qw1_1)

            if (wrk /= 0.0) then
              qw2_1 = qwsec * min(100.,max(0.,( 3.*qw1_2*wrk3-wrk4)/(aterm*wrk))) ! A.10
              qw2_2 = qwsec * min(100.,max(0.,(-3.*qw1_1*wrk3+wrk4)/(onema*wrk))) ! A.11
            else
              qw2_1 = 0.0
              qw2_2 = 0.0
            endif
!
            qw1_1 = qw1_1*sqrtqt + total_water !qw_first
            qw1_2 = qw1_2*sqrtqt + total_water !qw_first

            sqrtqw2_1 = sqrt(qw2_1)
            sqrtqw2_2 = sqrt(qw2_2)

          ENDIF   ! if qwsec small


!  CONVERT FROM TILDA VARIABLES TO "REAL" VARIABLES

          w1_1 = w1_1*sqrtw2 + w_first    ! using A.5 and A.6
          w1_2 = w1_2*sqrtw2 + w_first    ! note: this is already done for w2_x


!=== Assign PDF diagnostics ===!

          pdf_a = aterm
             
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
            pdf_sigw = sqrt( 0.4*w_sec )
          else
            pdf_sigw = 0.0
          end if

!==============================!


!  FIND WITHIN-PLUME CORRELATIONS 

          testvar = aterm*sqrtqw2_1*sqrtthl2_1 + onema*sqrtqw2_2*sqrtthl2_2

          IF (testvar == 0) THEN
            r_qwthl_1 = 0.
          ELSE
            r_qwthl_1 = max(-1.0,min(1.0,(qwthlsec-aterm*(qw1_1-total_water)*(thl1_1-thl_first)-onema*(qw1_2-total_water)*(thl1_2-thl_first))/testvar)) ! A.12
          ENDIF

!!! PDF params
          pdf_rqtth = r_qwthl_1
!!! end PDF params

!  BEGIN TO COMPUTE CLOUD PROPERTY STATISTICS
! This section follows Bogenschutz thesis Appendix A, based on
! Sommeria and Deardorff (1977) and Lewellen and Yoh (1993).

          Tl1_1 = thl1_1 - gamaz !+ fac_cond*qpl(i,j,k) + fac_sub*qpi(i,j,k)
          Tl1_2 = thl1_2 - gamaz !+ fac_cond*qpl(i,j,k) + fac_sub*qpi(i,j,k)

! Now compute qs

          esval1_1 = 0.
          esval1_2 = 0.
          esval2_1 = 0.
          esval2_2 = 0.
          om1      = 1.
          om2      = 1.
             
! Partition based on temperature for the first plume

!!! temporary
          Tl1_1 = max(min(Tl1_1,350.),170.)
          Tl1_2 = max(min(Tl1_2,350.),170.)
!!!

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
! temporary
!          if (abs(qs1*wrk-MAPL_EQsat(Tabs,pval))>0.001 .and. pval>2.*esval1_1 .and. pval>20000.) print *,'pval=',pval,' Tabs=',Tabs,' tl1_1=',Tl1_1,' qs1=',qs1*wrk

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
            rwqt = wqwsec/(sqrtqt*sqrtw2)
         else
            rwqt = 0.0
         end if
         if (sqrtthl>0.0 .AND. sqrtw2>0.0) then
            rwthl = wthlsec/(sqrtthl*sqrtw2)
         else
            rwthl = 0.0
         end if

         wql1 = C1*(cqt1*sqrt(w2_1)*sqrt(qw2_1)*rwqt-cthl1*sqrt(w2_1)*sqrt(thl2_1)*rwthl)
         wql2 = C2*(cqt2*sqrt(w2_2)*sqrt(qw2_2)*rwqt-cthl2*sqrt(w2_2)*sqrt(thl2_2)*rwthl)


! Compute the liquid water flux
          wqls = aterm * ((w1_1-w_first)*ql1+wql1) + onema * ((w1_2-w_first)*ql2+wql2)
!          wqls = aterm * ((w1_1-w_first)*ql1) + onema * ((w1_2-w_first)*ql2)
          wqis = aterm * ((w1_1-w_first)*qi1) + onema * ((w1_2-w_first)*qi2)
             
!         if (pval>95000.) print *,'aterm=',aterm,'  onema=',onema

             
! diagnostic buoyancy flux.  Includes effects from liquid water, ice
! condensate, liquid & ice precipitation
!         wrk = epsv * basetemp
          wrk = epsv * thv

          bastoeps = (rv/rgas) * thv   ! thetav / epsilon

!          if (pval>95000.) print *,'dblgss: wthlsec=',wthlsec,'  wqwsec=',wqwsec
!          if (pval>95000.) print *,'dblgss: wrk=',wrk,'  bastoeps=',bastoeps
!          if (pval>95000.) print *,'wqls=',wqls,'  wqis=',wqis

          wthv_sec = wthlsec + wrk*wqwsec                                     &
                   + (fac_cond-bastoeps)*wqls                                 &
                   + (fac_sub-bastoeps) *wqis

!                          + ((lstarn1/cp)-thv(i,j,k))*0.5*(wqp_sec(i,j,kd)+wqp_sec(i,j,ku))

!          if (pval>95000.) print *,'dblgss: wthv_sec=',wthv_sec

  end subroutine partition_dblgss


end module partition_pdf
