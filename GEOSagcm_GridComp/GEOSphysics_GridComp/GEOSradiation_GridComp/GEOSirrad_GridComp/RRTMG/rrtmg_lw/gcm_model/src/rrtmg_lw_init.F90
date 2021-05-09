module rrtmg_lw_init

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------
      use rrlw_wvn
      use rrtmg_lw_setcoef, only: lwatmref, lwavplank, lwavplankderiv

      implicit none

      contains

! **************************************************************************
      subroutine rrtmg_lw_ini
! **************************************************************************
!
!  Original version:       Michael J. Iacono; July, 1998
!  First revision for GCMs:   September, 1998
!  Second revision for RRTM_V3.0:  September, 2002
!
!  This subroutine performs calculations necessary for the initialization
!  of the longwave model.  Lookup tables are computed for use in the LW
!  radiative transfer, and input absorption coefficient data for each
!  spectral band are reduced from 256 g-point intervals to 140.
! **************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw
      use rrlw_tbl, only: ntbl, tblint, pade, bpade, tau_tbl, exp_tbl, tfn_tbl
      use rrlw_vsn, only: hvrini, hnamini

! ------- Local -------

      integer  :: itr, ibnd, igc, ig, ind, ipr 
      integer  :: igcsm, iprsm

      real  :: wtsum, wtsm(mg)        !
      real  :: tfn                    !

      real , parameter :: expeps = 1.e-20    ! Smallest value for exponential table

! ------- Definitions -------
!     Arrays for 10000-point look-up tables:
!     TAU_TBL Clear-sky optical depth (used in cloudy radiative transfer)
!     EXP_TBL Exponential lookup table for ransmittance
!     TFN_TBL Tau transition function; i.e. the transition of the Planck
!             function from that for the mean layer temperature to that for
!             the layer boundary temperature as a function of optical depth.
!             The "linear in tau" method is used to make the table.
!     PADE    Pade approximation constant (= 0.278)
!     BPADE   Inverse of the Pade approximation constant
!

      hvrini = '$Revision$'

! Initialize model data
      call lwdatinit
      call lwcmbdat               ! g-point interval reduction data
      call lwcldpr                ! cloud optical properties
      call lwatmref               ! reference MLS profile
      call lwavplank              ! Planck function 
      call lwavplankderiv         ! Planck function derivative wrt temp
      call lw_kgb01               ! molecular absorption coefficients
      call lw_kgb02
      call lw_kgb03
      call lw_kgb04
      call lw_kgb05
      call lw_kgb06
      call lw_kgb07
      call lw_kgb08
      call lw_kgb09
      call lw_kgb10
      call lw_kgb11
      call lw_kgb12
      call lw_kgb13
      call lw_kgb14
      call lw_kgb15
      call lw_kgb16

! Compute lookup tables for transmittance, tau transition function,
! and clear sky tau (for the cloudy sky radiative transfer).  Tau is 
! computed as a function of the tau transition function, transmittance 
! is calculated as a function of tau, and the tau transition function 
! is calculated using the linear in tau formulation at values of tau 
! above 0.01.  TF is approximated as tau/6 for tau < 0.01.  All tables 
! are computed at intervals of 0.001.  The inverse of the constant used
! in the Pade approximation to the tau transition function is set to b.

      tau_tbl(0) = 0.0 
      tau_tbl(ntbl) = 1.e10 
      exp_tbl(0) = 1.0 
      exp_tbl(ntbl) = expeps
      tfn_tbl(0) = 0.0 
      tfn_tbl(ntbl) = 1.0 
      bpade = 1.0  / pade
      do itr = 1, ntbl-1
         tfn = float(itr) / float(ntbl)
         tau_tbl(itr) = bpade * tfn / (1.  - tfn)
         exp_tbl(itr) = exp(-tau_tbl(itr))
         if (exp_tbl(itr) .le. expeps) exp_tbl(itr) = expeps
         if (tau_tbl(itr) .lt. 0.06 ) then
            tfn_tbl(itr) = tau_tbl(itr)/6. 
         else
            tfn_tbl(itr) = 1. -2. *((1. /tau_tbl(itr))-(exp_tbl(itr)/(1.-exp_tbl(itr))))
         endif
      enddo

! Perform g-point reduction from 16 per band (256 total points) to
! a band dependant number (140 total points) for all absorption
! coefficient input data and Planck fraction input data.
! Compute relative weighting for new g-point combinations.

      igcsm = 0
      do ibnd = 1,nbndlw
         iprsm = 0
         if (ngc(ibnd).lt.mg) then
            do igc = 1,ngc(ibnd) 
               igcsm = igcsm + 1
               wtsum = 0. 
               do ipr = 1, ngn(igcsm)
                  iprsm = iprsm + 1
                  wtsum = wtsum + wt(iprsm)
               enddo
               wtsm(igc) = wtsum
            enddo
            do ig = 1, ng(ibnd)
               ind = (ibnd-1)*mg + ig
               rwgt(ind) = wt(ig)/wtsm(ngm(ind))
            enddo
         else
            do ig = 1, ng(ibnd)
               igcsm = igcsm + 1
               ind = (ibnd-1)*mg + ig
               rwgt(ind) = 1.0 
            enddo
         endif
      enddo

! Reduce g-points for absorption coefficient data in each LW spectral band.

      call cmbgb1
      call cmbgb2
      call cmbgb3
      call cmbgb4
      call cmbgb5
      call cmbgb6
      call cmbgb7
      call cmbgb8
      call cmbgb9
      call cmbgb10
      call cmbgb11
      call cmbgb12
      call cmbgb13
      call cmbgb14
      call cmbgb15
      call cmbgb16

      end subroutine rrtmg_lw_ini

!***************************************************************************
      subroutine lwdatinit
!***************************************************************************

! --------- Modules ----------

      use parrrtm, only : maxxsec, maxinpx
      use rrlw_con, only: heatfac, grav, planck, boltz, &
                          clight, avogad, alosmt, gascon, radcn1, radcn2, &
                          sbcnst, secdy 
      use rrlw_vsn

      save 
 
! Longwave spectral band limits (wavenumbers)
      wavenum1(:) = (/ 10. , 350. , 500. , 630. , 700. , 820. , &
                      980. ,1080. ,1180. ,1390. ,1480. ,1800. , &
                     2080. ,2250. ,2380. ,2600. /)
      wavenum2(:) = (/350. , 500. , 630. , 700. , 820. , 980. , &
                     1080. ,1180. ,1390. ,1480. ,1800. ,2080. , &
                     2250. ,2380. ,2600. ,3250. /)
      delwave(:) =  (/340. , 150. , 130. ,  70. , 120. , 160. , &
                      100. , 100. , 210. ,  90. , 320. , 280. , &
                      170. , 130. , 220. , 650. /)

! Spectral band information
      ng(:) = (/16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/)
      nspa(:) = (/1,1,9,9,9,1,9,1,9,1,1,9,9,1,9,9/)
      nspb(:) = (/1,1,5,5,5,0,1,1,1,1,1,0,0,1,0,0/)

!     nxmol     - number of cross-sections input by user
!     ixindx(i) - index of cross-section molecule corresponding to Ith
!                 cross-section specified by user
!                 = 0 -- not allowed in rrtm
!                 = 1 -- ccl4
!                 = 2 -- cfc11
!                 = 3 -- cfc12
!                 = 4 -- cfc22
      nxmol = 4
      ixindx(1) = 1
      ixindx(2) = 2
      ixindx(3) = 3
      ixindx(4) = 4
      ixindx(5:maxinpx) = 0

! Fundamental physical constants from NIST 2002

      grav = 9.8066                         ! Acceleration of gravity
                                              ! (m s-2)
      planck = 6.62606876e-27               ! Planck constant
                                              ! (ergs s; g cm2 s-1)
      boltz = 1.3806503e-16                 ! Boltzmann constant
                                              ! (ergs K-1; g cm2 s-2 K-1)
      clight = 2.99792458e+10               ! Speed of light in a vacuum  
                                              ! (cm s-1)
      avogad = 6.02214199e+23               ! Avogadro constant
                                              ! (mol-1)
      alosmt = 2.6867775e+19                ! Loschmidt constant
                                              ! (cm-3)
      gascon = 8.31447200e+07               ! Molar gas constant
                                              ! (ergs mol-1 K-1)
      radcn1 = 1.191042722e-12              ! First radiation constant
                                              ! (W cm2 sr-1)
      radcn2 = 1.4387752                    ! Second radiation constant
                                              ! (cm K)
      sbcnst = 5.670400e-04                 ! Stefan-Boltzmann constant
                                              ! (W cm-2 K-4)
      secdy = 8.6400e4                      ! Number of seconds per day
                                              ! (s d-1)
!
!     units are generally cgs
!
!     The first and second radiation constants are taken from NIST.
!     They were previously obtained from the relations:
!          radcn1 = 2.*planck*clight*clight*1.e-07
!          radcn2 = planck*clight/boltz

      end subroutine lwdatinit

!***************************************************************************
      subroutine lwcmbdat
!***************************************************************************

      save
 
! ------- Definitions -------
!     Arrays for the g-point reduction from 256 to 140 for the 16 LW bands:
!     This mapping from 256 to 140 points has been carefully selected to 
!     minimize the effect on the resulting fluxes and cooling rates, and
!     caution should be used if the mapping is modified.  The full 256
!     g-point set can be restored with ngptlw=256, ngc=16*16, ngn=256*1., etc.
!     ngptlw  The total number of new g-points
!     ngc     The number of new g-points in each band
!     ngs     The cumulative sum of new g-points for each band
!     ngm     The index of each new g-point relative to the original
!             16 g-points for each band.  
!     ngn     The number of original g-points that are combined to make
!             each new g-point in each band.
!     ngb     The band index for each new g-point.
!     wt      RRTM weights for 16 g-points.

! ------- Data statements -------
      ngc(:) = (/10,12,16,14,16,8,12,8,12,6,8,8,4,2,2,2/)
      ngs(:) = (/10,22,38,52,68,76,88,96,108,114,122,130,134,136,138,140/)
      ngm(:) = (/1,2,3,3,4,4,5,5,6,6,7,7,8,8,9,10, &          ! band 1
                 1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 2
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 3
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,14,14, &    ! band 4
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 5
                 1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 6
                 1,1,2,2,3,4,5,6,7,8,9,10,11,11,12,12, &      ! band 7
                 1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 8
                 1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 9
                 1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6, &           ! band 10
                 1,2,3,3,4,4,5,5,6,6,7,7,7,8,8,8, &           ! band 11
                 1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 12
                 1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,4, &           ! band 13
                 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 14
                 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 15
                 1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2/)            ! band 16
      ngn(:) = (/1,1,2,2,2,2,2,2,1,1, &                       ! band 1
                 1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 2
                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 3
                 1,1,1,1,1,1,1,1,1,1,1,1,1,3, &               ! band 4
                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 5
                 2,2,2,2,2,2,2,2, &                           ! band 6
                 2,2,1,1,1,1,1,1,1,1,2,2, &                   ! band 7
                 2,2,2,2,2,2,2,2, &                           ! band 8
                 1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 9
                 2,2,2,2,4,4, &                               ! band 10
                 1,1,2,2,2,2,3,3, &                           ! band 11
                 1,1,1,1,2,2,4,4, &                           ! band 12
                 3,3,4,6, &                                   ! band 13
                 8,8, &                                       ! band 14
                 8,8, &                                       ! band 15
                 4,12/)                                       ! band 16
      ngb(:) = (/1,1,1,1,1,1,1,1,1,1, &                       ! band 1
                 2,2,2,2,2,2,2,2,2,2,2,2, &                   ! band 2
                 3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, &           ! band 3
                 4,4,4,4,4,4,4,4,4,4,4,4,4,4, &               ! band 4
                 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, &           ! band 5
                 6,6,6,6,6,6,6,6, &                           ! band 6
                 7,7,7,7,7,7,7,7,7,7,7,7, &                   ! band 7
                 8,8,8,8,8,8,8,8, &                           ! band 8
                 9,9,9,9,9,9,9,9,9,9,9,9, &                   ! band 9
                 10,10,10,10,10,10, &                         ! band 10
                 11,11,11,11,11,11,11,11, &                   ! band 11
                 12,12,12,12,12,12,12,12, &                   ! band 12
                 13,13,13,13, &                               ! band 13
                 14,14, &                                     ! band 14
                 15,15, &                                     ! band 15
                 16,16/)                                      ! band 16
      wt(:) = (/ 0.1527534276 , 0.1491729617 , 0.1420961469 , &
                 0.1316886544 , 0.1181945205 , 0.1019300893 , &
                 0.0832767040 , 0.0626720116 , 0.0424925000 , &
                 0.0046269894 , 0.0038279891 , 0.0030260086 , &
                 0.0022199750 , 0.0014140010 , 0.0005330000 , &
                 0.0000750000 /)

      end subroutine lwcmbdat

!***************************************************************************
      subroutine cmbgb1
!***************************************************************************
!
!  Original version:    MJIacono; July 1998
!  Revision for GCMs:   MJIacono; September 1998
!  Revision for RRTMG:  MJIacono, September 2002
!  Revision for F90 reformatting:  MJIacono, June 2006
!
!  The subroutines CMBGB1->CMBGB16 input the absorption coefficient
!  data for each band, which are defined for 16 g-points and 16 spectral
!  bands. The data are combined with appropriate weighting following the
!  g-point mapping arrays specified in RRTMINIT.  Plank fraction data
!  in arrays FRACREFA and FRACREFB are combined without weighting.  All
!  g-point reduced data are put into new arrays for use in RRTM.
!
!  band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                       (high key - h2o; high minor - n2)
!  note: previous versions of rrtm band 1: 
!        10-250 cm-1 (low - h2o; high - h2o)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng1
      use rrlw_kg01, only: fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, absa, ka, absb, kb, ka_mn2, kb_mn2, &
                           selfref, forref

! ------- Local -------
      integer  :: jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumk1, sumk2, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(1)
               sumk = 0.
               do ipr = 1, ngn(igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(1)
               sumk = 0.
               do ipr = 1, ngn(igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(1)
            sumk = 0.
            do ipr = 1, ngn(igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(1)
            sumk = 0.
            do ipr = 1, ngn(igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(1)
            sumk1 = 0.
            sumk2 = 0.
            do ipr = 1, ngn(igc)
               iprsm = iprsm + 1
               sumk1 = sumk1 + kao_mn2(jt,iprsm)*rwgt(iprsm)
               sumk2 = sumk2 + kbo_mn2(jt,iprsm)*rwgt(iprsm)
            enddo
            ka_mn2(jt,igc) = sumk1
            kb_mn2(jt,igc) = sumk2
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(1)
         sumf1 = 0.
         sumf2 = 0.
         do ipr = 1, ngn(igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb1

!***************************************************************************
      subroutine cmbgb2
!***************************************************************************
!
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
!
!     note: previous version of rrtm band 2: 
!           250 - 500 cm-1 (low - h2o; high - h2o)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng2
      use rrlw_kg02, only: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
                           fracrefa, fracrefb, absa, ka, absb, kb, selfref, forref

! ------- Local -------
      integer  :: jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(2)
               sumk = 0.
               do ipr = 1, ngn(ngs(1)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+16)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(2)
               sumk = 0.
               do ipr = 1, ngn(ngs(1)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+16)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(2)
            sumk = 0.
            do ipr = 1, ngn(ngs(1)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+16)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(2)
            sumk = 0.
            do ipr = 1, ngn(ngs(1)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+16)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(2)
         sumf1 = 0.
         sumf2 = 0.
         do ipr = 1, ngn(ngs(1)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb2

!***************************************************************************
      subroutine cmbgb3
!***************************************************************************
!
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
!                           (high key - h2o,co2; high minor - n2o)
!
! old band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng3
      use rrlw_kg03, only: fracrefao, fracrefbo, kao, kbo, kao_mn2o, kbo_mn2o, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, absa, ka, absb, kb, ka_mn2o, kb_mn2o, &
                           selfref, forref

! ------- Local -------
      integer  :: jn, jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(3)
                 sumk = 0.
                  do ipr = 1, ngn(ngs(2)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+32)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo
      do jn = 1,5
         do jt = 1,5
            do jp = 13,59
               iprsm = 0
               do igc = 1,ngc(3)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(2)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+32)
                  enddo
                  kb(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(3)
              sumk = 0.
               do ipr = 1, ngn(ngs(2)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mn2o(jn,jt,iprsm)*rwgt(iprsm+32)
               enddo
               ka_mn2o(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jn = 1,5
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(3)
              sumk = 0.
               do ipr = 1, ngn(ngs(2)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo_mn2o(jn,jt,iprsm)*rwgt(iprsm+32)
               enddo
               kb_mn2o(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(3)
            sumk = 0.
            do ipr = 1, ngn(ngs(2)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+32)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(3)
            sumk = 0.
            do ipr = 1, ngn(ngs(2)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+32)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(3)
            sumf = 0.
            do ipr = 1, ngn(ngs(2)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      do jp = 1,5
         iprsm = 0
         do igc = 1,ngc(3)
            sumf = 0.
            do ipr = 1, ngn(ngs(2)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefbo(iprsm,jp)
            enddo
            fracrefb(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb3

!***************************************************************************
      subroutine cmbgb4
!***************************************************************************
!
!     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
!
! old band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng4
      use rrlw_kg04, only: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
                           fracrefa, fracrefb, absa, ka, absb, kb, selfref, forref

! ------- Local -------
      integer  :: jn, jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(4)
                 sumk = 0.
                  do ipr = 1, ngn(ngs(3)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+48)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo
      do jn = 1,5
         do jt = 1,5
            do jp = 13,59
               iprsm = 0
               do igc = 1,ngc(4)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(3)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+48)
                  enddo
                  kb(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(4)
            sumk = 0.
            do ipr = 1, ngn(ngs(3)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+48)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(4)
            sumk = 0.
            do ipr = 1, ngn(ngs(3)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+48)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(4)
            sumf = 0.
            do ipr = 1, ngn(ngs(3)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      do jp = 1,5
         iprsm = 0
         do igc = 1,ngc(4)
            sumf = 0.
            do ipr = 1, ngn(ngs(3)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefbo(iprsm,jp)
            enddo
            fracrefb(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb4

!***************************************************************************
      subroutine cmbgb5
!***************************************************************************
!
!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!                           (high key - o3,co2)
!
! old band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng5
      use rrlw_kg05, only: fracrefao, fracrefbo, kao, kbo, kao_mo3, ccl4o, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, absa, ka, absb, kb, ka_mo3, ccl4, &
                           selfref, forref

! ------- Local -------
      integer  :: jn, jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(5)
                 sumk = 0.
                  do ipr = 1, ngn(ngs(4)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+64)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo
      do jn = 1,5
         do jt = 1,5
            do jp = 13,59
               iprsm = 0
               do igc = 1,ngc(5)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(4)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+64)
                  enddo
                  kb(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(5)
              sumk = 0.
               do ipr = 1, ngn(ngs(4)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mo3(jn,jt,iprsm)*rwgt(iprsm+64)
               enddo
               ka_mo3(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(5)
            sumk = 0.
            do ipr = 1, ngn(ngs(4)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+64)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(5)
            sumk = 0.
            do ipr = 1, ngn(ngs(4)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+64)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(5)
            sumf = 0.
            do ipr = 1, ngn(ngs(4)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      do jp = 1,5
         iprsm = 0
         do igc = 1,ngc(5)
            sumf = 0.
            do ipr = 1, ngn(ngs(4)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefbo(iprsm,jp)
            enddo
            fracrefb(igc,jp) = sumf
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(5)
         sumk = 0.
         do ipr = 1, ngn(ngs(4)+igc)
            iprsm = iprsm + 1
            sumk = sumk + ccl4o(iprsm)*rwgt(iprsm+64)
         enddo
         ccl4(igc) = sumk
      enddo

      end subroutine cmbgb5

!***************************************************************************
      subroutine cmbgb6
!***************************************************************************
!
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!                           (high key - nothing; high minor - cfc11, cfc12)
!
! old band 6:  820-980 cm-1 (low - h2o; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng6
      use rrlw_kg06, only: fracrefao, kao, kao_mco2, cfc11adjo, cfc12o, &
                           selfrefo, forrefo, &
                           fracrefa, absa, ka, ka_mco2, cfc11adj, cfc12, &
                           selfref, forref

! ------- Local -------
      integer  :: jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumf, sumk1, sumk2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(6)
               sumk = 0.
               do ipr = 1, ngn(ngs(5)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+80)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(6)
            sumk = 0.
            do ipr = 1, ngn(ngs(5)+igc)
               iprsm = iprsm + 1
               sumk = sumk + kao_mco2(jt,iprsm)*rwgt(iprsm+80)
            enddo
            ka_mco2(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(6)
            sumk = 0.
            do ipr = 1, ngn(ngs(5)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+80)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(6)
            sumk = 0.
            do ipr = 1, ngn(ngs(5)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+80)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(6)
         sumf = 0.
         sumk1= 0.
         sumk2= 0.
         do ipr = 1, ngn(ngs(5)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefao(iprsm)
            sumk1= sumk1+ cfc11adjo(iprsm)*rwgt(iprsm+80)
            sumk2= sumk2+ cfc12o(iprsm)*rwgt(iprsm+80)
         enddo
         fracrefa(igc) = sumf
         cfc11adj(igc) = sumk1
         cfc12(igc) = sumk2
      enddo

      end subroutine cmbgb6

!***************************************************************************
      subroutine cmbgb7
!***************************************************************************
!
!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
!                            (high key - o3; high minor - co2)
!
! old band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng7
      use rrlw_kg07, only: fracrefao, fracrefbo, kao, kbo, kao_mco2, kbo_mco2, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, absa, ka, absb, kb, ka_mco2, kb_mco2, &
                           selfref, forref

! ------- Local -------
      integer  :: jn, jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(7)
                 sumk = 0.
                  do ipr = 1, ngn(ngs(6)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+96)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo
      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(7)
               sumk = 0.
               do ipr = 1, ngn(ngs(6)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+96)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(7)
              sumk = 0.
               do ipr = 1, ngn(ngs(6)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mco2(jn,jt,iprsm)*rwgt(iprsm+96)
               enddo
               ka_mco2(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(7)
            sumk = 0.
            do ipr = 1, ngn(ngs(6)+igc)
               iprsm = iprsm + 1
               sumk = sumk + kbo_mco2(jt,iprsm)*rwgt(iprsm+96)
            enddo
            kb_mco2(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(7)
            sumk = 0.
            do ipr = 1, ngn(ngs(6)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+96)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(7)
            sumk = 0.
            do ipr = 1, ngn(ngs(6)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+96)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(7)
            sumf = 0.
            do ipr = 1, ngn(ngs(6)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(7)
         sumf = 0.
         do ipr = 1, ngn(ngs(6)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefbo(iprsm)
         enddo
         fracrefb(igc) = sumf
      enddo

      end subroutine cmbgb7

!***************************************************************************
      subroutine cmbgb8
!***************************************************************************
!
!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!                             (high key - o3; high minor - co2, n2o)
!
! old band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng8
      use rrlw_kg08, only: fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, &
                           kao_mo3, kbo, kbo_mco2, kbo_mn2o, selfrefo, forrefo, &
                           cfc12o, cfc22adjo, &
                           fracrefa, fracrefb, absa, ka, ka_mco2, ka_mn2o, &
                           ka_mo3, absb, kb, kb_mco2, kb_mn2o, selfref, forref, &
                           cfc12, cfc22adj

! ------- Local -------
      integer  :: jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumk1, sumk2, sumk3, sumk4, sumk5, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(8)
              sumk = 0.
               do ipr = 1, ngn(ngs(7)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+112)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo
      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(8)
               sumk = 0.
               do ipr = 1, ngn(ngs(7)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+112)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(8)
            sumk = 0.
            do ipr = 1, ngn(ngs(7)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+112)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(8)
            sumk = 0.
            do ipr = 1, ngn(ngs(7)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+112)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(8)
            sumk1 = 0.
            sumk2 = 0.
            sumk3 = 0.
            sumk4 = 0.
            sumk5 = 0.
            do ipr = 1, ngn(ngs(7)+igc)
               iprsm = iprsm + 1
               sumk1 = sumk1 + kao_mco2(jt,iprsm)*rwgt(iprsm+112)
               sumk2 = sumk2 + kbo_mco2(jt,iprsm)*rwgt(iprsm+112)
               sumk3 = sumk3 + kao_mo3(jt,iprsm)*rwgt(iprsm+112)
               sumk4 = sumk4 + kao_mn2o(jt,iprsm)*rwgt(iprsm+112)
               sumk5 = sumk5 + kbo_mn2o(jt,iprsm)*rwgt(iprsm+112)
            enddo
            ka_mco2(jt,igc) = sumk1
            kb_mco2(jt,igc) = sumk2
            ka_mo3(jt,igc) = sumk3
            ka_mn2o(jt,igc) = sumk4
            kb_mn2o(jt,igc) = sumk5
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(8)
         sumf1= 0.
         sumf2= 0.
         sumk1= 0.
         sumk2= 0.
         do ipr = 1, ngn(ngs(7)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
            sumk1= sumk1+ cfc12o(iprsm)*rwgt(iprsm+112)
            sumk2= sumk2+ cfc22adjo(iprsm)*rwgt(iprsm+112)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
         cfc12(igc) = sumk1
         cfc22adj(igc) = sumk2
      enddo

      end subroutine cmbgb8

!***************************************************************************
      subroutine cmbgb9
!***************************************************************************
!
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!                             (high key - ch4; high minor - n2o)!

! old band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng9
      use rrlw_kg09, only: fracrefao, fracrefbo, kao, kao_mn2o, &
                           kbo, kbo_mn2o, selfrefo, forrefo, &
                           fracrefa, fracrefb, absa, ka, ka_mn2o, &
                           absb, kb, kb_mn2o, selfref, forref

! ------- Local -------
      integer  :: jn, jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(9)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(8)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+128)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(9)
               sumk = 0.
               do ipr = 1, ngn(ngs(8)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+128)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(9)
              sumk = 0.
               do ipr = 1, ngn(ngs(8)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mn2o(jn,jt,iprsm)*rwgt(iprsm+128)
               enddo
               ka_mn2o(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(9)
            sumk = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumk = sumk + kbo_mn2o(jt,iprsm)*rwgt(iprsm+128)
            enddo
            kb_mn2o(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(9)
            sumk = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+128)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(9)
            sumk = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+128)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(9)
            sumf = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(9)
         sumf = 0.
         do ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefbo(iprsm)
         enddo
         fracrefb(igc) = sumf
      enddo

      end subroutine cmbgb9

!***************************************************************************
      subroutine cmbgb10
!***************************************************************************
!
!     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
!
! old band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng10
      use rrlw_kg10, only: fracrefao, fracrefbo, kao, kbo, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, absa, ka, absb, kb, &
                           selfref, forref

! ------- Local -------
      integer  :: jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(10)
               sumk = 0.
               do ipr = 1, ngn(ngs(9)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+144)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(10)
               sumk = 0.
               do ipr = 1, ngn(ngs(9)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+144)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(10)
            sumk = 0.
            do ipr = 1, ngn(ngs(9)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+144)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(10)
            sumk = 0.
            do ipr = 1, ngn(ngs(9)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+144)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(10)
         sumf1= 0.
         sumf2= 0.
         do ipr = 1, ngn(ngs(9)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb10

!***************************************************************************
      subroutine cmbgb11
!***************************************************************************
!
!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
!
! old band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng11
      use rrlw_kg11, only: fracrefao, fracrefbo, kao, kao_mo2, &
                           kbo, kbo_mo2, selfrefo, forrefo, &
                           fracrefa, fracrefb, absa, ka, ka_mo2, &
                           absb, kb, kb_mo2, selfref, forref

! ------- Local -------
      integer  :: jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumk1, sumk2, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(11)
               sumk = 0.
               do ipr = 1, ngn(ngs(10)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+160)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo
      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(11)
               sumk = 0.
               do ipr = 1, ngn(ngs(10)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+160)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(11)
            sumk1 = 0.
            sumk2 = 0.
            do ipr = 1, ngn(ngs(10)+igc)
               iprsm = iprsm + 1
               sumk1 = sumk1 + kao_mo2(jt,iprsm)*rwgt(iprsm+160)
               sumk2 = sumk2 + kbo_mo2(jt,iprsm)*rwgt(iprsm+160)
            enddo
            ka_mo2(jt,igc) = sumk1
            kb_mo2(jt,igc) = sumk2
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(11)
            sumk = 0.
            do ipr = 1, ngn(ngs(10)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+160)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(11)
            sumk = 0.
            do ipr = 1, ngn(ngs(10)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+160)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(11)
         sumf1= 0.
         sumf2= 0.
         do ipr = 1, ngn(ngs(10)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb11

!***************************************************************************
      subroutine cmbgb12
!***************************************************************************
!
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!
! old band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng12
      use rrlw_kg12, only: fracrefao, kao, selfrefo, forrefo, &
                           fracrefa, absa, ka, selfref, forref

! ------- Local -------
      integer  :: jn, jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(12)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(11)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+176)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(12)
            sumk = 0.
            do ipr = 1, ngn(ngs(11)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+176)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(12)
            sumk = 0.
            do ipr = 1, ngn(ngs(11)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+176)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(12)
            sumf = 0.
            do ipr = 1, ngn(ngs(11)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb12

!***************************************************************************
      subroutine cmbgb13
!***************************************************************************
!
!     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
!
! old band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng13
      use rrlw_kg13, only: fracrefao, fracrefbo, kao, kao_mco2, kao_mco, &
                           kbo_mo3, selfrefo, forrefo, &
                           fracrefa, fracrefb, absa, ka, ka_mco2, ka_mco, &
                           kb_mo3, selfref, forref

! ------- Local -------
      integer  :: jn, jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumk1, sumk2, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(13)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(12)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+192)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(13)
              sumk1 = 0.
              sumk2 = 0.
               do ipr = 1, ngn(ngs(12)+igc)
                  iprsm = iprsm + 1
                  sumk1 = sumk1 + kao_mco2(jn,jt,iprsm)*rwgt(iprsm+192)
                  sumk2 = sumk2 + kao_mco(jn,jt,iprsm)*rwgt(iprsm+192)
               enddo
               ka_mco2(jn,jt,igc) = sumk1
               ka_mco(jn,jt,igc) = sumk2
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(13)
            sumk = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumk = sumk + kbo_mo3(jt,iprsm)*rwgt(iprsm+192)
            enddo
            kb_mo3(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(13)
            sumk = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+192)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(13)
            sumk = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+192)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(13)
         sumf = 0.
         do ipr = 1, ngn(ngs(12)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefbo(iprsm)
         enddo
         fracrefb(igc) = sumf
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(13)
            sumf = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb13

!***************************************************************************
      subroutine cmbgb14
!***************************************************************************
!
!     band 14:  2250-2380 cm-1 (low - co2; high - co2)
!
! old band 14:  2250-2380 cm-1 (low - co2; high - co2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng14
      use rrlw_kg14, only: fracrefao, fracrefbo, kao, kbo, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, absa, ka, absb, kb, &
                           selfref, forref

! ------- Local -------
      integer  :: jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(14)
               sumk = 0.
               do ipr = 1, ngn(ngs(13)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+208)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(14)
               sumk = 0.
               do ipr = 1, ngn(ngs(13)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+208)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(14)
            sumk = 0.
            do ipr = 1, ngn(ngs(13)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+208)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(14)
            sumk = 0.
            do ipr = 1, ngn(ngs(13)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+208)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(14)
         sumf1= 0.
         sumf2= 0.
         do ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb14

!***************************************************************************
      subroutine cmbgb15
!***************************************************************************
!
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                              (high - nothing)
!
! old band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng15
      use rrlw_kg15, only: fracrefao, kao, kao_mn2, selfrefo, forrefo, &
                           fracrefa, absa, ka, ka_mn2, selfref, forref

! ------- Local -------
      integer  :: jn, jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(15)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(14)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+224)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(15)
              sumk = 0.
               do ipr = 1, ngn(ngs(14)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mn2(jn,jt,iprsm)*rwgt(iprsm+224)
               enddo
               ka_mn2(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(15)
            sumk = 0.
            do ipr = 1, ngn(ngs(14)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+224)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(15)
            sumk = 0.
            do ipr = 1, ngn(ngs(14)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+224)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(15)
            sumf = 0.
            do ipr = 1, ngn(ngs(14)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb15

!***************************************************************************
      subroutine cmbgb16
!***************************************************************************
!
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!
! old band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng16
      use rrlw_kg16, only: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
                           fracrefa, fracrefb, absa, ka, absb, kb, selfref, forref

! ------- Local -------
      integer  :: jn, jt, jp, igc, ipr, iprsm 
      real  :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(16)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(15)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+240)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(16)
               sumk = 0.
               do ipr = 1, ngn(ngs(15)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+240)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(16)
            sumk = 0.
            do ipr = 1, ngn(ngs(15)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+240)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(16)
            sumk = 0.
            do ipr = 1, ngn(ngs(15)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+240)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(16)
         sumf = 0.
         do ipr = 1, ngn(ngs(15)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefbo(iprsm)
         enddo
         fracrefb(igc) = sumf
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(16)
            sumf = 0.
            do ipr = 1, ngn(ngs(15)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb16

!***************************************************************************
      subroutine lwcldpr
!***************************************************************************

! --------- Modules ----------

      use rrlw_cld, only: absliq1, absice0, absice1, absice2, absice3, absice4

      save

! ABSICEn(J,IB) are the parameters needed to compute the liquid water 
! absorption coefficient in spectral region IB for ICEFLAG=n.  The units
! of ABSICEn(1,IB) are m2/g and ABSICEn(2,IB) has units (microns (m2/g)).
! For ICEFLAG = 0.

      absice0(:)= (/0.005 ,  1.0 /)

! For ICEFLAG = 1.
      absice1(1,:) = (/0.0036 , 0.0068 , 0.0003 , 0.0016 , 0.0020 /)
      absice1(2,:) = (/1.136  , 0.600  , 1.338  , 1.166  , 1.118  /)

! For ICEFLAG = 2.  In each band, the absorption
! coefficients are listed for a range of effective radii from 5.0
! to 131.0 microns in increments of 3.0 microns.
! Spherical Ice Particle Parameterization
! absorption units (abs coef/iwc): [(m^-1)/(g m^-3)]
      absice2(:,1) = (/ &
! band 1
       7.798999e-02 ,6.340479e-02 ,5.417973e-02 ,4.766245e-02 ,4.272663e-02 , &
       3.880939e-02 ,3.559544e-02 ,3.289241e-02 ,3.057511e-02 ,2.855800e-02 , &
       2.678022e-02 ,2.519712e-02 ,2.377505e-02 ,2.248806e-02 ,2.131578e-02 , &
       2.024194e-02 ,1.925337e-02 ,1.833926e-02 ,1.749067e-02 ,1.670007e-02 , &
       1.596113e-02 ,1.526845e-02 ,1.461739e-02 ,1.400394e-02 ,1.342462e-02 , &
       1.287639e-02 ,1.235656e-02 ,1.186279e-02 ,1.139297e-02 ,1.094524e-02 , &
       1.051794e-02 ,1.010956e-02 ,9.718755e-03 ,9.344316e-03 ,8.985139e-03 , &
       8.640223e-03 ,8.308656e-03 ,7.989606e-03 ,7.682312e-03 ,7.386076e-03 , &
       7.100255e-03 ,6.824258e-03 ,6.557540e-03 /)
      absice2(:,2) = (/ &
! band 2
       2.784879e-02 ,2.709863e-02 ,2.619165e-02 ,2.529230e-02 ,2.443225e-02 , &
       2.361575e-02 ,2.284021e-02 ,2.210150e-02 ,2.139548e-02 ,2.071840e-02 , &
       2.006702e-02 ,1.943856e-02 ,1.883064e-02 ,1.824120e-02 ,1.766849e-02 , &
       1.711099e-02 ,1.656737e-02 ,1.603647e-02 ,1.551727e-02 ,1.500886e-02 , &
       1.451045e-02 ,1.402132e-02 ,1.354084e-02 ,1.306842e-02 ,1.260355e-02 , &
       1.214575e-02 ,1.169460e-02 ,1.124971e-02 ,1.081072e-02 ,1.037731e-02 , &
       9.949167e-03 ,9.526021e-03 ,9.107615e-03 ,8.693714e-03 ,8.284096e-03 , &
       7.878558e-03 ,7.476910e-03 ,7.078974e-03 ,6.684586e-03 ,6.293589e-03 , &
       5.905839e-03 ,5.521200e-03 ,5.139543e-03 /)
      absice2(:,3) = (/ &
! band 3
       1.065397e-01 ,8.005726e-02 ,6.546428e-02 ,5.589131e-02 ,4.898681e-02 , &
       4.369932e-02 ,3.947901e-02 ,3.600676e-02 ,3.308299e-02 ,3.057561e-02 , &
       2.839325e-02 ,2.647040e-02 ,2.475872e-02 ,2.322164e-02 ,2.183091e-02 , &
       2.056430e-02 ,1.940407e-02 ,1.833586e-02 ,1.734787e-02 ,1.643034e-02 , &
       1.557512e-02 ,1.477530e-02 ,1.402501e-02 ,1.331924e-02 ,1.265364e-02 , &
       1.202445e-02 ,1.142838e-02 ,1.086257e-02 ,1.032445e-02 ,9.811791e-03 , &
       9.322587e-03 ,8.855053e-03 ,8.407591e-03 ,7.978763e-03 ,7.567273e-03 , &
       7.171949e-03 ,6.791728e-03 ,6.425642e-03 ,6.072809e-03 ,5.732424e-03 , &
       5.403748e-03 ,5.086103e-03 ,4.778865e-03 /)
      absice2(:,4) = (/ &
! band 4
       1.804566e-01 ,1.168987e-01 ,8.680442e-02 ,6.910060e-02 ,5.738174e-02 , &
       4.902332e-02 ,4.274585e-02 ,3.784923e-02 ,3.391734e-02 ,3.068690e-02 , &
       2.798301e-02 ,2.568480e-02 ,2.370600e-02 ,2.198337e-02 ,2.046940e-02 , &
       1.912777e-02 ,1.793016e-02 ,1.685420e-02 ,1.588193e-02 ,1.499882e-02 , &
       1.419293e-02 ,1.345440e-02 ,1.277496e-02 ,1.214769e-02 ,1.156669e-02 , &
       1.102694e-02 ,1.052412e-02 ,1.005451e-02 ,9.614854e-03 ,9.202335e-03 , &
       8.814470e-03 ,8.449077e-03 ,8.104223e-03 ,7.778195e-03 ,7.469466e-03 , &
       7.176671e-03 ,6.898588e-03 ,6.634117e-03 ,6.382264e-03 ,6.142134e-03 , &
       5.912913e-03 ,5.693862e-03 ,5.484308e-03 /)
      absice2(:,5) = (/ &
! band 5
       2.131806e-01 ,1.311372e-01 ,9.407171e-02 ,7.299442e-02 ,5.941273e-02 , &
       4.994043e-02 ,4.296242e-02 ,3.761113e-02 ,3.337910e-02 ,2.994978e-02 , &
       2.711556e-02 ,2.473461e-02 ,2.270681e-02 ,2.095943e-02 ,1.943839e-02 , &
       1.810267e-02 ,1.692057e-02 ,1.586719e-02 ,1.492275e-02 ,1.407132e-02 , &
       1.329989e-02 ,1.259780e-02 ,1.195618e-02 ,1.136761e-02 ,1.082583e-02 , &
       1.032552e-02 ,9.862158e-03 ,9.431827e-03 ,9.031157e-03 ,8.657217e-03 , &
       8.307449e-03 ,7.979609e-03 ,7.671724e-03 ,7.382048e-03 ,7.109032e-03 , &
       6.851298e-03 ,6.607615e-03 ,6.376881e-03 ,6.158105e-03 ,5.950394e-03 , &
       5.752942e-03 ,5.565019e-03 ,5.385963e-03 /)
      absice2(:,6) = (/ &
! band 6
       1.546177e-01 ,1.039251e-01 ,7.910347e-02 ,6.412429e-02 ,5.399997e-02 , &
       4.664937e-02 ,4.104237e-02 ,3.660781e-02 ,3.300218e-02 ,3.000586e-02 , &
       2.747148e-02 ,2.529633e-02 ,2.340647e-02 ,2.174723e-02 ,2.027731e-02 , &
       1.896487e-02 ,1.778492e-02 ,1.671761e-02 ,1.574692e-02 ,1.485978e-02 , &
       1.404543e-02 ,1.329489e-02 ,1.260066e-02 ,1.195636e-02 ,1.135657e-02 , &
       1.079664e-02 ,1.027257e-02 ,9.780871e-03 ,9.318505e-03 ,8.882815e-03 , &
       8.471458e-03 ,8.082364e-03 ,7.713696e-03 ,7.363817e-03 ,7.031264e-03 , &
       6.714725e-03 ,6.413021e-03 ,6.125086e-03 ,5.849958e-03 ,5.586764e-03 , &
       5.334707e-03 ,5.093066e-03 ,4.861179e-03 /)
      absice2(:,7) = (/ &
! band 7
       7.583404e-02 ,6.181558e-02 ,5.312027e-02 ,4.696039e-02 ,4.225986e-02 , &
       3.849735e-02 ,3.538340e-02 ,3.274182e-02 ,3.045798e-02 ,2.845343e-02 , &
       2.667231e-02 ,2.507353e-02 ,2.362606e-02 ,2.230595e-02 ,2.109435e-02 , &
       1.997617e-02 ,1.893916e-02 ,1.797328e-02 ,1.707016e-02 ,1.622279e-02 , &
       1.542523e-02 ,1.467241e-02 ,1.395997e-02 ,1.328414e-02 ,1.264164e-02 , &
       1.202958e-02 ,1.144544e-02 ,1.088697e-02 ,1.035218e-02 ,9.839297e-03 , &
       9.346733e-03 ,8.873057e-03 ,8.416980e-03 ,7.977335e-03 ,7.553066e-03 , &
       7.143210e-03 ,6.746888e-03 ,6.363297e-03 ,5.991700e-03 ,5.631422e-03 , &
       5.281840e-03 ,4.942378e-03 ,4.612505e-03 /)
      absice2(:,8) = (/ &
! band 8
       9.022185e-02 ,6.922700e-02 ,5.710674e-02 ,4.898377e-02 ,4.305946e-02 , &
       3.849553e-02 ,3.484183e-02 ,3.183220e-02 ,2.929794e-02 ,2.712627e-02 , &
       2.523856e-02 ,2.357810e-02 ,2.210286e-02 ,2.078089e-02 ,1.958747e-02 , &
       1.850310e-02 ,1.751218e-02 ,1.660205e-02 ,1.576232e-02 ,1.498440e-02 , &
       1.426107e-02 ,1.358624e-02 ,1.295474e-02 ,1.236212e-02 ,1.180456e-02 , &
       1.127874e-02 ,1.078175e-02 ,1.031106e-02 ,9.864433e-03 ,9.439878e-03 , &
       9.035637e-03 ,8.650140e-03 ,8.281981e-03 ,7.929895e-03 ,7.592746e-03 , &
       7.269505e-03 ,6.959238e-03 ,6.661100e-03 ,6.374317e-03 ,6.098185e-03 , &
       5.832059e-03 ,5.575347e-03 ,5.327504e-03 /)
      absice2(:,9) = (/ &
! band 9
       1.294087e-01 ,8.788217e-02 ,6.728288e-02 ,5.479720e-02 ,4.635049e-02 , &
       4.022253e-02 ,3.555576e-02 ,3.187259e-02 ,2.888498e-02 ,2.640843e-02 , &
       2.431904e-02 ,2.253038e-02 ,2.098024e-02 ,1.962267e-02 ,1.842293e-02 , &
       1.735426e-02 ,1.639571e-02 ,1.553060e-02 ,1.474552e-02 ,1.402953e-02 , &
       1.337363e-02 ,1.277033e-02 ,1.221336e-02 ,1.169741e-02 ,1.121797e-02 , &
       1.077117e-02 ,1.035369e-02 ,9.962643e-03 ,9.595509e-03 ,9.250088e-03 , &
       8.924447e-03 ,8.616876e-03 ,8.325862e-03 ,8.050057e-03 ,7.788258e-03 , &
       7.539388e-03 ,7.302478e-03 ,7.076656e-03 ,6.861134e-03 ,6.655197e-03 , &
       6.458197e-03 ,6.269543e-03 ,6.088697e-03 /)
      absice2(:,10) = (/ &
! band 10
       1.593628e-01 ,1.014552e-01 ,7.458955e-02 ,5.903571e-02 ,4.887582e-02 , &
       4.171159e-02 ,3.638480e-02 ,3.226692e-02 ,2.898717e-02 ,2.631256e-02 , &
       2.408925e-02 ,2.221156e-02 ,2.060448e-02 ,1.921325e-02 ,1.799699e-02 , &
       1.692456e-02 ,1.597177e-02 ,1.511961e-02 ,1.435289e-02 ,1.365933e-02 , &
       1.302890e-02 ,1.245334e-02 ,1.192576e-02 ,1.144037e-02 ,1.099230e-02 , &
       1.057739e-02 ,1.019208e-02 ,9.833302e-03 ,9.498395e-03 ,9.185047e-03 , &
       8.891237e-03 ,8.615185e-03 ,8.355325e-03 ,8.110267e-03 ,7.878778e-03 , &
       7.659759e-03 ,7.452224e-03 ,7.255291e-03 ,7.068166e-03 ,6.890130e-03 , &
       6.720536e-03 ,6.558794e-03 ,6.404371e-03 /)
      absice2(:,11) = (/ &
! band 11
       1.656227e-01 ,1.032129e-01 ,7.487359e-02 ,5.871431e-02 ,4.828355e-02 , &
       4.099989e-02 ,3.562924e-02 ,3.150755e-02 ,2.824593e-02 ,2.560156e-02 , &
       2.341503e-02 ,2.157740e-02 ,2.001169e-02 ,1.866199e-02 ,1.748669e-02 , &
       1.645421e-02 ,1.554015e-02 ,1.472535e-02 ,1.399457e-02 ,1.333553e-02 , &
       1.273821e-02 ,1.219440e-02 ,1.169725e-02 ,1.124104e-02 ,1.082096e-02 , &
       1.043290e-02 ,1.007336e-02 ,9.739338e-03 ,9.428223e-03 ,9.137756e-03 , &
       8.865964e-03 ,8.611115e-03 ,8.371686e-03 ,8.146330e-03 ,7.933852e-03 , &
       7.733187e-03 ,7.543386e-03 ,7.363597e-03 ,7.193056e-03 ,7.031072e-03 , &
       6.877024e-03 ,6.730348e-03 ,6.590531e-03 /)
      absice2(:,12) = (/ &
! band 12
       9.194591e-02 ,6.446867e-02 ,4.962034e-02 ,4.042061e-02 ,3.418456e-02 , &
       2.968856e-02 ,2.629900e-02 ,2.365572e-02 ,2.153915e-02 ,1.980791e-02 , &
       1.836689e-02 ,1.714979e-02 ,1.610900e-02 ,1.520946e-02 ,1.442476e-02 , &
       1.373468e-02 ,1.312345e-02 ,1.257858e-02 ,1.209010e-02 ,1.164990e-02 , &
       1.125136e-02 ,1.088901e-02 ,1.055827e-02 ,1.025531e-02 ,9.976896e-03 , &
       9.720255e-03 ,9.483022e-03 ,9.263160e-03 ,9.058902e-03 ,8.868710e-03 , &
       8.691240e-03 ,8.525312e-03 ,8.369886e-03 ,8.224042e-03 ,8.086961e-03 , &
       7.957917e-03 ,7.836258e-03 ,7.721400e-03 ,7.612821e-03 ,7.510045e-03 , &
       7.412648e-03 ,7.320242e-03 ,7.232476e-03 /)
      absice2(:,13) = (/ &
! band 13
       1.437021e-01 ,8.872535e-02 ,6.392420e-02 ,4.991833e-02 ,4.096790e-02 , &
       3.477881e-02 ,3.025782e-02 ,2.681909e-02 ,2.412102e-02 ,2.195132e-02 , &
       2.017124e-02 ,1.868641e-02 ,1.743044e-02 ,1.635529e-02 ,1.542540e-02 , &
       1.461388e-02 ,1.390003e-02 ,1.326766e-02 ,1.270395e-02 ,1.219860e-02 , &
       1.174326e-02 ,1.133107e-02 ,1.095637e-02 ,1.061442e-02 ,1.030126e-02 , &
       1.001352e-02 ,9.748340e-03 ,9.503256e-03 ,9.276155e-03 ,9.065205e-03 , &
       8.868808e-03 ,8.685571e-03 ,8.514268e-03 ,8.353820e-03 ,8.203272e-03 , &
       8.061776e-03 ,7.928578e-03 ,7.803001e-03 ,7.684443e-03 ,7.572358e-03 , &
       7.466258e-03 ,7.365701e-03 ,7.270286e-03 /)
      absice2(:,14) = (/ &
! band 14
       1.288870e-01 ,8.160295e-02 ,5.964745e-02 ,4.703790e-02 ,3.888637e-02 , &
       3.320115e-02 ,2.902017e-02 ,2.582259e-02 ,2.330224e-02 ,2.126754e-02 , &
       1.959258e-02 ,1.819130e-02 ,1.700289e-02 ,1.598320e-02 ,1.509942e-02 , &
       1.432666e-02 ,1.364572e-02 ,1.304156e-02 ,1.250220e-02 ,1.201803e-02 , &
       1.158123e-02 ,1.118537e-02 ,1.082513e-02 ,1.049605e-02 ,1.019440e-02 , &
       9.916989e-03 ,9.661116e-03 ,9.424457e-03 ,9.205005e-03 ,9.001022e-03 , &
       8.810992e-03 ,8.633588e-03 ,8.467646e-03 ,8.312137e-03 ,8.166151e-03 , &
       8.028878e-03 ,7.899597e-03 ,7.777663e-03 ,7.662498e-03 ,7.553581e-03 , &
       7.450444e-03 ,7.352662e-03 ,7.259851e-03 /)
      absice2(:,15) = (/ &
! band 15
       8.254229e-02 ,5.808787e-02 ,4.492166e-02 ,3.675028e-02 ,3.119623e-02 , &
       2.718045e-02 ,2.414450e-02 ,2.177073e-02 ,1.986526e-02 ,1.830306e-02 , &
       1.699991e-02 ,1.589698e-02 ,1.495199e-02 ,1.413374e-02 ,1.341870e-02 , &
       1.278883e-02 ,1.223002e-02 ,1.173114e-02 ,1.128322e-02 ,1.087900e-02 , &
       1.051254e-02 ,1.017890e-02 ,9.873991e-03 ,9.594347e-03 ,9.337044e-03 , &
       9.099589e-03 ,8.879842e-03 ,8.675960e-03 ,8.486341e-03 ,8.309594e-03 , &
       8.144500e-03 ,7.989986e-03 ,7.845109e-03 ,7.709031e-03 ,7.581007e-03 , &
       7.460376e-03 ,7.346544e-03 ,7.238978e-03 ,7.137201e-03 ,7.040780e-03 , &
       6.949325e-03 ,6.862483e-03 ,6.779931e-03 /)
      absice2(:,16) = (/ &
! band 16
       1.382062e-01 ,8.643227e-02 ,6.282935e-02 ,4.934783e-02 ,4.063891e-02 , &
       3.455591e-02 ,3.007059e-02 ,2.662897e-02 ,2.390631e-02 ,2.169972e-02 , &
       1.987596e-02 ,1.834393e-02 ,1.703924e-02 ,1.591513e-02 ,1.493679e-02 , &
       1.407780e-02 ,1.331775e-02 ,1.264061e-02 ,1.203364e-02 ,1.148655e-02 , &
       1.099099e-02 ,1.054006e-02 ,1.012807e-02 ,9.750215e-03 ,9.402477e-03 , &
       9.081428e-03 ,8.784143e-03 ,8.508107e-03 ,8.251146e-03 ,8.011373e-03 , &
       7.787140e-03 ,7.577002e-03 ,7.379687e-03 ,7.194071e-03 ,7.019158e-03 , &
       6.854061e-03 ,6.697986e-03 ,6.550224e-03 ,6.410138e-03 ,6.277153e-03 , &
       6.150751e-03 ,6.030462e-03 ,5.915860e-03 /)

! ICEFLAG = 3; Fu parameterization. Particle size 5 - 140 micron in 
! increments of 3 microns.
! units = m2/g
! Hexagonal Ice Particle Parameterization
! absorption units (abs coef/iwc): [(m^-1)/(g m^-3)]
      absice3(:,1) = (/ &
! band 1
       3.110649e-03 ,4.666352e-02 ,6.606447e-02 ,6.531678e-02 ,6.012598e-02 , &
       5.437494e-02 ,4.906411e-02 ,4.441146e-02 ,4.040585e-02 ,3.697334e-02 , &
       3.403027e-02 ,3.149979e-02 ,2.931596e-02 ,2.742365e-02 ,2.577721e-02 , &
       2.433888e-02 ,2.307732e-02 ,2.196644e-02 ,2.098437e-02 ,2.011264e-02 , &
       1.933561e-02 ,1.863992e-02 ,1.801407e-02 ,1.744812e-02 ,1.693346e-02 , &
       1.646252e-02 ,1.602866e-02 ,1.562600e-02 ,1.524933e-02 ,1.489399e-02 , &
       1.455580e-02 ,1.423098e-02 ,1.391612e-02 ,1.360812e-02 ,1.330413e-02 , &
       1.300156e-02 ,1.269801e-02 ,1.239127e-02 ,1.207928e-02 ,1.176014e-02 , &
       1.143204e-02 ,1.109334e-02 ,1.074243e-02 ,1.037786e-02 ,9.998198e-03 , &
       9.602126e-03 /)
      absice3(:,2) = (/ &
! band 2
       3.984966e-04 ,1.681097e-02 ,2.627680e-02 ,2.767465e-02 ,2.700722e-02 , &
       2.579180e-02 ,2.448677e-02 ,2.323890e-02 ,2.209096e-02 ,2.104882e-02 , &
       2.010547e-02 ,1.925003e-02 ,1.847128e-02 ,1.775883e-02 ,1.710358e-02 , &
       1.649769e-02 ,1.593449e-02 ,1.540829e-02 ,1.491429e-02 ,1.444837e-02 , &
       1.400704e-02 ,1.358729e-02 ,1.318654e-02 ,1.280258e-02 ,1.243346e-02 , &
       1.207750e-02 ,1.173325e-02 ,1.139941e-02 ,1.107487e-02 ,1.075861e-02 , &
       1.044975e-02 ,1.014753e-02 ,9.851229e-03 ,9.560240e-03 ,9.274003e-03 , &
       8.992020e-03 ,8.713845e-03 ,8.439074e-03 ,8.167346e-03 ,7.898331e-03 , &
       7.631734e-03 ,7.367286e-03 ,7.104742e-03 ,6.843882e-03 ,6.584504e-03 , &
       6.326424e-03 /)
      absice3(:,3) = (/ &
! band 3
       6.933163e-02 ,8.540475e-02 ,7.701816e-02 ,6.771158e-02 ,5.986953e-02 , &
       5.348120e-02 ,4.824962e-02 ,4.390563e-02 ,4.024411e-02 ,3.711404e-02 , &
       3.440426e-02 ,3.203200e-02 ,2.993478e-02 ,2.806474e-02 ,2.638464e-02 , &
       2.486516e-02 ,2.348288e-02 ,2.221890e-02 ,2.105780e-02 ,1.998687e-02 , &
       1.899552e-02 ,1.807490e-02 ,1.721750e-02 ,1.641693e-02 ,1.566773e-02 , &
       1.496515e-02 ,1.430509e-02 ,1.368398e-02 ,1.309865e-02 ,1.254634e-02 , &
       1.202456e-02 ,1.153114e-02 ,1.106409e-02 ,1.062166e-02 ,1.020224e-02 , &
       9.804381e-03 ,9.426771e-03 ,9.068205e-03 ,8.727578e-03 ,8.403876e-03 , &
       8.096160e-03 ,7.803564e-03 ,7.525281e-03 ,7.260560e-03 ,7.008697e-03 , &
       6.769036e-03 /)
      absice3(:,4) = (/ &
! band 4
       1.765735e-01 ,1.382700e-01 ,1.095129e-01 ,8.987475e-02 ,7.591185e-02 , &
       6.554169e-02 ,5.755500e-02 ,5.122083e-02 ,4.607610e-02 ,4.181475e-02 , &
       3.822697e-02 ,3.516432e-02 ,3.251897e-02 ,3.021073e-02 ,2.817876e-02 , &
       2.637607e-02 ,2.476582e-02 ,2.331871e-02 ,2.201113e-02 ,2.082388e-02 , &
       1.974115e-02 ,1.874983e-02 ,1.783894e-02 ,1.699922e-02 ,1.622280e-02 , &
       1.550296e-02 ,1.483390e-02 ,1.421064e-02 ,1.362880e-02 ,1.308460e-02 , &
       1.257468e-02 ,1.209611e-02 ,1.164628e-02 ,1.122287e-02 ,1.082381e-02 , &
       1.044725e-02 ,1.009154e-02 ,9.755166e-03 ,9.436783e-03 ,9.135163e-03 , &
       8.849193e-03 ,8.577856e-03 ,8.320225e-03 ,8.075451e-03 ,7.842755e-03 , &
       7.621418e-03 /)
      absice3(:,5) = (/ &
! band 5
       2.339673e-01 ,1.692124e-01 ,1.291656e-01 ,1.033837e-01 ,8.562949e-02 , &
       7.273526e-02 ,6.298262e-02 ,5.537015e-02 ,4.927787e-02 ,4.430246e-02 , &
       4.017061e-02 ,3.669072e-02 ,3.372455e-02 ,3.116995e-02 ,2.894977e-02 , &
       2.700471e-02 ,2.528842e-02 ,2.376420e-02 ,2.240256e-02 ,2.117959e-02 , &
       2.007567e-02 ,1.907456e-02 ,1.816271e-02 ,1.732874e-02 ,1.656300e-02 , &
       1.585725e-02 ,1.520445e-02 ,1.459852e-02 ,1.403419e-02 ,1.350689e-02 , &
       1.301260e-02 ,1.254781e-02 ,1.210941e-02 ,1.169468e-02 ,1.130118e-02 , &
       1.092675e-02 ,1.056945e-02 ,1.022757e-02 ,9.899560e-03 ,9.584021e-03 , &
       9.279705e-03 ,8.985479e-03 ,8.700322e-03 ,8.423306e-03 ,8.153590e-03 , &
       7.890412e-03 /)
      absice3(:,6) = (/ &
! band 6
       1.145369e-01 ,1.174566e-01 ,9.917866e-02 ,8.332990e-02 ,7.104263e-02 , &
       6.153370e-02 ,5.405472e-02 ,4.806281e-02 ,4.317918e-02 ,3.913795e-02 , &
       3.574916e-02 ,3.287437e-02 ,3.041067e-02 ,2.828017e-02 ,2.642292e-02 , &
       2.479206e-02 ,2.335051e-02 ,2.206851e-02 ,2.092195e-02 ,1.989108e-02 , &
       1.895958e-02 ,1.811385e-02 ,1.734245e-02 ,1.663573e-02 ,1.598545e-02 , &
       1.538456e-02 ,1.482700e-02 ,1.430750e-02 ,1.382150e-02 ,1.336499e-02 , &
       1.293447e-02 ,1.252685e-02 ,1.213939e-02 ,1.176968e-02 ,1.141555e-02 , &
       1.107508e-02 ,1.074655e-02 ,1.042839e-02 ,1.011923e-02 ,9.817799e-03 , &
       9.522962e-03 ,9.233688e-03 ,8.949041e-03 ,8.668171e-03 ,8.390301e-03 , &
       8.114723e-03 /)
      absice3(:,7) = (/ &
! band 7
       1.222345e-02 ,5.344230e-02 ,5.523465e-02 ,5.128759e-02 ,4.676925e-02 , &
       4.266150e-02 ,3.910561e-02 ,3.605479e-02 ,3.342843e-02 ,3.115052e-02 , &
       2.915776e-02 ,2.739935e-02 ,2.583499e-02 ,2.443266e-02 ,2.316681e-02 , &
       2.201687e-02 ,2.096619e-02 ,2.000112e-02 ,1.911044e-02 ,1.828481e-02 , &
       1.751641e-02 ,1.679866e-02 ,1.612598e-02 ,1.549360e-02 ,1.489742e-02 , &
       1.433392e-02 ,1.380002e-02 ,1.329305e-02 ,1.281068e-02 ,1.235084e-02 , &
       1.191172e-02 ,1.149171e-02 ,1.108936e-02 ,1.070341e-02 ,1.033271e-02 , &
       9.976220e-03 ,9.633021e-03 ,9.302273e-03 ,8.983216e-03 ,8.675161e-03 , &
       8.377478e-03 ,8.089595e-03 ,7.810986e-03 ,7.541170e-03 ,7.279706e-03 , &
       7.026186e-03 /)
      absice3(:,8) = (/ &
! band 8
       6.711058e-02 ,6.918198e-02 ,6.127484e-02 ,5.411944e-02 ,4.836902e-02 , &
       4.375293e-02 ,3.998077e-02 ,3.683587e-02 ,3.416508e-02 ,3.186003e-02 , &
       2.984290e-02 ,2.805671e-02 ,2.645895e-02 ,2.501733e-02 ,2.370689e-02 , &
       2.250808e-02 ,2.140532e-02 ,2.038609e-02 ,1.944018e-02 ,1.855918e-02 , &
       1.773609e-02 ,1.696504e-02 ,1.624106e-02 ,1.555990e-02 ,1.491793e-02 , &
       1.431197e-02 ,1.373928e-02 ,1.319743e-02 ,1.268430e-02 ,1.219799e-02 , &
       1.173682e-02 ,1.129925e-02 ,1.088393e-02 ,1.048961e-02 ,1.011516e-02 , &
       9.759543e-03 ,9.421813e-03 ,9.101089e-03 ,8.796559e-03 ,8.507464e-03 , &
       8.233098e-03 ,7.972798e-03 ,7.725942e-03 ,7.491940e-03 ,7.270238e-03 , &
       7.060305e-03 /)
      absice3(:,9) = (/ &
! band 9
       1.236780e-01 ,9.222386e-02 ,7.383997e-02 ,6.204072e-02 ,5.381029e-02 , &
       4.770678e-02 ,4.296928e-02 ,3.916131e-02 ,3.601540e-02 ,3.335878e-02 , &
       3.107493e-02 ,2.908247e-02 ,2.732282e-02 ,2.575276e-02 ,2.433968e-02 , &
       2.305852e-02 ,2.188966e-02 ,2.081757e-02 ,1.982974e-02 ,1.891599e-02 , &
       1.806794e-02 ,1.727865e-02 ,1.654227e-02 ,1.585387e-02 ,1.520924e-02 , &
       1.460476e-02 ,1.403730e-02 ,1.350416e-02 ,1.300293e-02 ,1.253153e-02 , &
       1.208808e-02 ,1.167094e-02 ,1.127862e-02 ,1.090979e-02 ,1.056323e-02 , &
       1.023786e-02 ,9.932665e-03 ,9.646744e-03 ,9.379250e-03 ,9.129409e-03 , &
       8.896500e-03 ,8.679856e-03 ,8.478852e-03 ,8.292904e-03 ,8.121463e-03 , &
       7.964013e-03 /)
      absice3(:,10) = (/ &
! band 10
       1.655966e-01 ,1.134205e-01 ,8.714344e-02 ,7.129241e-02 ,6.063739e-02 , &
       5.294203e-02 ,4.709309e-02 ,4.247476e-02 ,3.871892e-02 ,3.559206e-02 , &
       3.293893e-02 ,3.065226e-02 ,2.865558e-02 ,2.689288e-02 ,2.532221e-02 , &
       2.391150e-02 ,2.263582e-02 ,2.147549e-02 ,2.041476e-02 ,1.944089e-02 , &
       1.854342e-02 ,1.771371e-02 ,1.694456e-02 ,1.622989e-02 ,1.556456e-02 , &
       1.494415e-02 ,1.436491e-02 ,1.382354e-02 ,1.331719e-02 ,1.284339e-02 , &
       1.239992e-02 ,1.198486e-02 ,1.159647e-02 ,1.123323e-02 ,1.089375e-02 , &
       1.057679e-02 ,1.028124e-02 ,1.000607e-02 ,9.750376e-03 ,9.513303e-03 , &
       9.294082e-03 ,9.092003e-03 ,8.906412e-03 ,8.736702e-03 ,8.582314e-03 , &
       8.442725e-03 /)
      absice3(:,11) = (/ &
! band 11
       1.775615e-01 ,1.180046e-01 ,8.929607e-02 ,7.233500e-02 ,6.108333e-02 , &
       5.303642e-02 ,4.696927e-02 ,4.221206e-02 ,3.836768e-02 ,3.518576e-02 , &
       3.250063e-02 ,3.019825e-02 ,2.819758e-02 ,2.643943e-02 ,2.487953e-02 , &
       2.348414e-02 ,2.222705e-02 ,2.108762e-02 ,2.004936e-02 ,1.909892e-02 , &
       1.822539e-02 ,1.741975e-02 ,1.667449e-02 ,1.598330e-02 ,1.534084e-02 , &
       1.474253e-02 ,1.418446e-02 ,1.366325e-02 ,1.317597e-02 ,1.272004e-02 , &
       1.229321e-02 ,1.189350e-02 ,1.151915e-02 ,1.116859e-02 ,1.084042e-02 , &
       1.053338e-02 ,1.024636e-02 ,9.978326e-03 ,9.728357e-03 ,9.495613e-03 , &
       9.279327e-03 ,9.078798e-03 ,8.893383e-03 ,8.722488e-03 ,8.565568e-03 , &
       8.422115e-03 /)
      absice3(:,12) = (/ &
! band 12
       9.465447e-02 ,6.432047e-02 ,5.060973e-02 ,4.267283e-02 ,3.741843e-02 , &
       3.363096e-02 ,3.073531e-02 ,2.842405e-02 ,2.651789e-02 ,2.490518e-02 , &
       2.351273e-02 ,2.229056e-02 ,2.120335e-02 ,2.022541e-02 ,1.933763e-02 , &
       1.852546e-02 ,1.777763e-02 ,1.708528e-02 ,1.644134e-02 ,1.584009e-02 , &
       1.527684e-02 ,1.474774e-02 ,1.424955e-02 ,1.377957e-02 ,1.333549e-02 , &
       1.291534e-02 ,1.251743e-02 ,1.214029e-02 ,1.178265e-02 ,1.144337e-02 , &
       1.112148e-02 ,1.081609e-02 ,1.052642e-02 ,1.025178e-02 ,9.991540e-03 , &
       9.745130e-03 ,9.512038e-03 ,9.291797e-03 ,9.083980e-03 ,8.888195e-03 , &
       8.704081e-03 ,8.531306e-03 ,8.369560e-03 ,8.218558e-03 ,8.078032e-03 , &
       7.947730e-03 /)
      absice3(:,13) = (/ &
! band 13
       1.560311e-01 ,9.961097e-02 ,7.502949e-02 ,6.115022e-02 ,5.214952e-02 , &
       4.578149e-02 ,4.099731e-02 ,3.724174e-02 ,3.419343e-02 ,3.165356e-02 , &
       2.949251e-02 ,2.762222e-02 ,2.598073e-02 ,2.452322e-02 ,2.321642e-02 , &
       2.203516e-02 ,2.096002e-02 ,1.997579e-02 ,1.907036e-02 ,1.823401e-02 , &
       1.745879e-02 ,1.673819e-02 ,1.606678e-02 ,1.544003e-02 ,1.485411e-02 , &
       1.430574e-02 ,1.379215e-02 ,1.331092e-02 ,1.285996e-02 ,1.243746e-02 , &
       1.204183e-02 ,1.167164e-02 ,1.132567e-02 ,1.100281e-02 ,1.070207e-02 , &
       1.042258e-02 ,1.016352e-02 ,9.924197e-03 ,9.703953e-03 ,9.502199e-03 , &
       9.318400e-03 ,9.152066e-03 ,9.002749e-03 ,8.870038e-03 ,8.753555e-03 , &
       8.652951e-03 /)
      absice3(:,14) = (/ &
! band 14
       1.559547e-01 ,9.896700e-02 ,7.441231e-02 ,6.061469e-02 ,5.168730e-02 , &
       4.537821e-02 ,4.064106e-02 ,3.692367e-02 ,3.390714e-02 ,3.139438e-02 , &
       2.925702e-02 ,2.740783e-02 ,2.578547e-02 ,2.434552e-02 ,2.305506e-02 , &
       2.188910e-02 ,2.082842e-02 ,1.985789e-02 ,1.896553e-02 ,1.814165e-02 , &
       1.737839e-02 ,1.666927e-02 ,1.600891e-02 ,1.539279e-02 ,1.481712e-02 , &
       1.427865e-02 ,1.377463e-02 ,1.330266e-02 ,1.286068e-02 ,1.244689e-02 , &
       1.205973e-02 ,1.169780e-02 ,1.135989e-02 ,1.104492e-02 ,1.075192e-02 , &
       1.048004e-02 ,1.022850e-02 ,9.996611e-03 ,9.783753e-03 ,9.589361e-03 , &
       9.412924e-03 ,9.253977e-03 ,9.112098e-03 ,8.986903e-03 ,8.878039e-03 , &
       8.785184e-03 /)
      absice3(:,15) = (/ &
! band 15
       1.102926e-01 ,7.176622e-02 ,5.530316e-02 ,4.606056e-02 ,4.006116e-02 , &
       3.579628e-02 ,3.256909e-02 ,3.001360e-02 ,2.791920e-02 ,2.615617e-02 , &
       2.464023e-02 ,2.331426e-02 ,2.213817e-02 ,2.108301e-02 ,2.012733e-02 , &
       1.925493e-02 ,1.845331e-02 ,1.771269e-02 ,1.702531e-02 ,1.638493e-02 , &
       1.578648e-02 ,1.522579e-02 ,1.469940e-02 ,1.420442e-02 ,1.373841e-02 , &
       1.329931e-02 ,1.288535e-02 ,1.249502e-02 ,1.212700e-02 ,1.178015e-02 , &
       1.145348e-02 ,1.114612e-02 ,1.085730e-02 ,1.058633e-02 ,1.033263e-02 , &
       1.009564e-02 ,9.874895e-03 ,9.669960e-03 ,9.480449e-03 ,9.306014e-03 , &
       9.146339e-03 ,9.001138e-03 ,8.870154e-03 ,8.753148e-03 ,8.649907e-03 , &
       8.560232e-03 /)
      absice3(:,16) = (/ &
! band 16
       1.688344e-01 ,1.077072e-01 ,7.994467e-02 ,6.403862e-02 ,5.369850e-02 , &
       4.641582e-02 ,4.099331e-02 ,3.678724e-02 ,3.342069e-02 ,3.065831e-02 , &
       2.834557e-02 ,2.637680e-02 ,2.467733e-02 ,2.319286e-02 ,2.188299e-02 , &
       2.071701e-02 ,1.967121e-02 ,1.872692e-02 ,1.786931e-02 ,1.708641e-02 , &
       1.636846e-02 ,1.570743e-02 ,1.509665e-02 ,1.453052e-02 ,1.400433e-02 , &
       1.351407e-02 ,1.305631e-02 ,1.262810e-02 ,1.222688e-02 ,1.185044e-02 , &
       1.149683e-02 ,1.116436e-02 ,1.085153e-02 ,1.055701e-02 ,1.027961e-02 , &
       1.001831e-02 ,9.772141e-03 ,9.540280e-03 ,9.321966e-03 ,9.116517e-03 , &
       8.923315e-03 ,8.741803e-03 ,8.571472e-03 ,8.411860e-03 ,8.262543e-03 , &
       8.123136e-03 /)

! For LIQFLAG = 1. In each band, the absorption
! coefficients are listed for a range of effective radii from 2.5
! to 59.5 microns in increments of 1.0 micron.
      absliq1(:, 1) = (/ &
! band  1
       1.64047e-03 , 6.90533e-02 , 7.72017e-02 , 7.78054e-02 , 7.69523e-02 , &
       7.58058e-02 , 7.46400e-02 , 7.35123e-02 , 7.24162e-02 , 7.13225e-02 , &
       6.99145e-02 , 6.66409e-02 , 6.36582e-02 , 6.09425e-02 , 5.84593e-02 , &
       5.61743e-02 , 5.40571e-02 , 5.20812e-02 , 5.02245e-02 , 4.84680e-02 , &
       4.67959e-02 , 4.51944e-02 , 4.36516e-02 , 4.21570e-02 , 4.07015e-02 , &
       3.92766e-02 , 3.78747e-02 , 3.64886e-02 , 3.53632e-02 , 3.41992e-02 , &
       3.31016e-02 , 3.20643e-02 , 3.10817e-02 , 3.01490e-02 , 2.92620e-02 , &
       2.84171e-02 , 2.76108e-02 , 2.68404e-02 , 2.61031e-02 , 2.53966e-02 , &
       2.47189e-02 , 2.40678e-02 , 2.34418e-02 , 2.28392e-02 , 2.22586e-02 , &
       2.16986e-02 , 2.11580e-02 , 2.06356e-02 , 2.01305e-02 , 1.96417e-02 , &
       1.91682e-02 , 1.87094e-02 , 1.82643e-02 , 1.78324e-02 , 1.74129e-02 , &
       1.70052e-02 , 1.66088e-02 , 1.62231e-02 /)
      absliq1(:, 2) = (/ &
! band  2
       2.19486e-01 , 1.80687e-01 , 1.59150e-01 , 1.44731e-01 , 1.33703e-01 , &
       1.24355e-01 , 1.15756e-01 , 1.07318e-01 , 9.86119e-02 , 8.92739e-02 , &
       8.34911e-02 , 7.70773e-02 , 7.15240e-02 , 6.66615e-02 , 6.23641e-02 , &
       5.85359e-02 , 5.51020e-02 , 5.20032e-02 , 4.91916e-02 , 4.66283e-02 , &
       4.42813e-02 , 4.21236e-02 , 4.01330e-02 , 3.82905e-02 , 3.65797e-02 , &
       3.49869e-02 , 3.35002e-02 , 3.21090e-02 , 3.08957e-02 , 2.97601e-02 , &
       2.86966e-02 , 2.76984e-02 , 2.67599e-02 , 2.58758e-02 , 2.50416e-02 , &
       2.42532e-02 , 2.35070e-02 , 2.27997e-02 , 2.21284e-02 , 2.14904e-02 , &
       2.08834e-02 , 2.03051e-02 , 1.97536e-02 , 1.92271e-02 , 1.87239e-02 , &
       1.82425e-02 , 1.77816e-02 , 1.73399e-02 , 1.69162e-02 , 1.65094e-02 , &
       1.61187e-02 , 1.57430e-02 , 1.53815e-02 , 1.50334e-02 , 1.46981e-02 , &
       1.43748e-02 , 1.40628e-02 , 1.37617e-02 /)
      absliq1(:, 3) = (/ &
! band  3
       2.95174e-01 , 2.34765e-01 , 1.98038e-01 , 1.72114e-01 , 1.52083e-01 , &
       1.35654e-01 , 1.21613e-01 , 1.09252e-01 , 9.81263e-02 , 8.79448e-02 , &
       8.12566e-02 , 7.44563e-02 , 6.86374e-02 , 6.36042e-02 , 5.92094e-02 , &
       5.53402e-02 , 5.19087e-02 , 4.88455e-02 , 4.60951e-02 , 4.36124e-02 , &
       4.13607e-02 , 3.93096e-02 , 3.74338e-02 , 3.57119e-02 , 3.41261e-02 , &
       3.26610e-02 , 3.13036e-02 , 3.00425e-02 , 2.88497e-02 , 2.78077e-02 , &
       2.68317e-02 , 2.59158e-02 , 2.50545e-02 , 2.42430e-02 , 2.34772e-02 , &
       2.27533e-02 , 2.20679e-02 , 2.14181e-02 , 2.08011e-02 , 2.02145e-02 , &
       1.96561e-02 , 1.91239e-02 , 1.86161e-02 , 1.81311e-02 , 1.76673e-02 , &
       1.72234e-02 , 1.67981e-02 , 1.63903e-02 , 1.59989e-02 , 1.56230e-02 , &
       1.52615e-02 , 1.49138e-02 , 1.45791e-02 , 1.42565e-02 , 1.39455e-02 , &
       1.36455e-02 , 1.33559e-02 , 1.30761e-02 /)
      absliq1(:, 4) = (/ &
! band  4
       3.00925e-01 , 2.36949e-01 , 1.96947e-01 , 1.68692e-01 , 1.47190e-01 , &
       1.29986e-01 , 1.15719e-01 , 1.03568e-01 , 9.30028e-02 , 8.36658e-02 , &
       7.71075e-02 , 7.07002e-02 , 6.52284e-02 , 6.05024e-02 , 5.63801e-02 , &
       5.27534e-02 , 4.95384e-02 , 4.66690e-02 , 4.40925e-02 , 4.17664e-02 , &
       3.96559e-02 , 3.77326e-02 , 3.59727e-02 , 3.43561e-02 , 3.28662e-02 , &
       3.14885e-02 , 3.02110e-02 , 2.90231e-02 , 2.78948e-02 , 2.69109e-02 , &
       2.59884e-02 , 2.51217e-02 , 2.43058e-02 , 2.35364e-02 , 2.28096e-02 , &
       2.21218e-02 , 2.14700e-02 , 2.08515e-02 , 2.02636e-02 , 1.97041e-02 , &
       1.91711e-02 , 1.86625e-02 , 1.81769e-02 , 1.77126e-02 , 1.72683e-02 , &
       1.68426e-02 , 1.64344e-02 , 1.60427e-02 , 1.56664e-02 , 1.53046e-02 , &
       1.49565e-02 , 1.46214e-02 , 1.42985e-02 , 1.39871e-02 , 1.36866e-02 , &
       1.33965e-02 , 1.31162e-02 , 1.28453e-02 /)
      absliq1(:, 5) = (/ &
! band  5
       2.64691e-01 , 2.12018e-01 , 1.78009e-01 , 1.53539e-01 , 1.34721e-01 , &
       1.19580e-01 , 1.06996e-01 , 9.62772e-02 , 8.69710e-02 , 7.87670e-02 , &
       7.29272e-02 , 6.70920e-02 , 6.20977e-02 , 5.77732e-02 , 5.39910e-02 , &
       5.06538e-02 , 4.76866e-02 , 4.50301e-02 , 4.26374e-02 , 4.04704e-02 , &
       3.84981e-02 , 3.66948e-02 , 3.50394e-02 , 3.35141e-02 , 3.21038e-02 , &
       3.07957e-02 , 2.95788e-02 , 2.84438e-02 , 2.73790e-02 , 2.64390e-02 , &
       2.55565e-02 , 2.47263e-02 , 2.39437e-02 , 2.32047e-02 , 2.25056e-02 , &
       2.18433e-02 , 2.12149e-02 , 2.06177e-02 , 2.00495e-02 , 1.95081e-02 , &
       1.89917e-02 , 1.84984e-02 , 1.80269e-02 , 1.75755e-02 , 1.71431e-02 , &
       1.67283e-02 , 1.63303e-02 , 1.59478e-02 , 1.55801e-02 , 1.52262e-02 , &
       1.48853e-02 , 1.45568e-02 , 1.42400e-02 , 1.39342e-02 , 1.36388e-02 , &
       1.33533e-02 , 1.30773e-02 , 1.28102e-02 /)
      absliq1(:, 6) = (/ &
! band  6
       8.81182e-02 , 1.06745e-01 , 9.79753e-02 , 8.99625e-02 , 8.35200e-02 , &
       7.81899e-02 , 7.35939e-02 , 6.94696e-02 , 6.56266e-02 , 6.19148e-02 , &
       5.83355e-02 , 5.49306e-02 , 5.19642e-02 , 4.93325e-02 , 4.69659e-02 , &
       4.48148e-02 , 4.28431e-02 , 4.10231e-02 , 3.93332e-02 , 3.77563e-02 , &
       3.62785e-02 , 3.48882e-02 , 3.35758e-02 , 3.23333e-02 , 3.11536e-02 , &
       3.00310e-02 , 2.89601e-02 , 2.79365e-02 , 2.70502e-02 , 2.62618e-02 , &
       2.55025e-02 , 2.47728e-02 , 2.40726e-02 , 2.34013e-02 , 2.27583e-02 , &
       2.21422e-02 , 2.15522e-02 , 2.09869e-02 , 2.04453e-02 , 1.99260e-02 , &
       1.94280e-02 , 1.89501e-02 , 1.84913e-02 , 1.80506e-02 , 1.76270e-02 , &
       1.72196e-02 , 1.68276e-02 , 1.64500e-02 , 1.60863e-02 , 1.57357e-02 , &
       1.53975e-02 , 1.50710e-02 , 1.47558e-02 , 1.44511e-02 , 1.41566e-02 , &
       1.38717e-02 , 1.35960e-02 , 1.33290e-02 /)
      absliq1(:, 7) = (/ &
! band  7
       4.32174e-02 , 7.36078e-02 , 6.98340e-02 , 6.65231e-02 , 6.41948e-02 , &
       6.23551e-02 , 6.06638e-02 , 5.88680e-02 , 5.67124e-02 , 5.38629e-02 , &
       4.99579e-02 , 4.86289e-02 , 4.70120e-02 , 4.52854e-02 , 4.35466e-02 , &
       4.18480e-02 , 4.02169e-02 , 3.86658e-02 , 3.71992e-02 , 3.58168e-02 , &
       3.45155e-02 , 3.32912e-02 , 3.21390e-02 , 3.10538e-02 , 3.00307e-02 , &
       2.90651e-02 , 2.81524e-02 , 2.72885e-02 , 2.62821e-02 , 2.55744e-02 , &
       2.48799e-02 , 2.42029e-02 , 2.35460e-02 , 2.29108e-02 , 2.22981e-02 , &
       2.17079e-02 , 2.11402e-02 , 2.05945e-02 , 2.00701e-02 , 1.95663e-02 , &
       1.90824e-02 , 1.86174e-02 , 1.81706e-02 , 1.77411e-02 , 1.73281e-02 , &
       1.69307e-02 , 1.65483e-02 , 1.61801e-02 , 1.58254e-02 , 1.54835e-02 , &
       1.51538e-02 , 1.48358e-02 , 1.45288e-02 , 1.42322e-02 , 1.39457e-02 , &
       1.36687e-02 , 1.34008e-02 , 1.31416e-02 /)
      absliq1(:, 8) = (/ &
! band  8
       1.41881e-01 , 7.15419e-02 , 6.30335e-02 , 6.11132e-02 , 6.01931e-02 , &
       5.92420e-02 , 5.78968e-02 , 5.58876e-02 , 5.28923e-02 , 4.84462e-02 , &
       4.60839e-02 , 4.56013e-02 , 4.45410e-02 , 4.31866e-02 , 4.17026e-02 , &
       4.01850e-02 , 3.86892e-02 , 3.72461e-02 , 3.58722e-02 , 3.45749e-02 , &
       3.33564e-02 , 3.22155e-02 , 3.11494e-02 , 3.01541e-02 , 2.92253e-02 , &
       2.83584e-02 , 2.75488e-02 , 2.67925e-02 , 2.57692e-02 , 2.50704e-02 , &
       2.43918e-02 , 2.37350e-02 , 2.31005e-02 , 2.24888e-02 , 2.18996e-02 , &
       2.13325e-02 , 2.07870e-02 , 2.02623e-02 , 1.97577e-02 , 1.92724e-02 , &
       1.88056e-02 , 1.83564e-02 , 1.79241e-02 , 1.75079e-02 , 1.71070e-02 , &
       1.67207e-02 , 1.63482e-02 , 1.59890e-02 , 1.56424e-02 , 1.53077e-02 , &
       1.49845e-02 , 1.46722e-02 , 1.43702e-02 , 1.40782e-02 , 1.37955e-02 , &
       1.35219e-02 , 1.32569e-02 , 1.30000e-02 /)
      absliq1(:, 9) = (/ &
! band  9
       6.72726e-02 , 6.61013e-02 , 6.47866e-02 , 6.33780e-02 , 6.18985e-02 , &
       6.03335e-02 , 5.86136e-02 , 5.65876e-02 , 5.39839e-02 , 5.03536e-02 , &
       4.71608e-02 , 4.63630e-02 , 4.50313e-02 , 4.34526e-02 , 4.17876e-02 , &
       4.01261e-02 , 3.85171e-02 , 3.69860e-02 , 3.55442e-02 , 3.41954e-02 , &
       3.29384e-02 , 3.17693e-02 , 3.06832e-02 , 2.96745e-02 , 2.87374e-02 , &
       2.78662e-02 , 2.70557e-02 , 2.63008e-02 , 2.52450e-02 , 2.45424e-02 , &
       2.38656e-02 , 2.32144e-02 , 2.25885e-02 , 2.19873e-02 , 2.14099e-02 , &
       2.08554e-02 , 2.03230e-02 , 1.98116e-02 , 1.93203e-02 , 1.88482e-02 , &
       1.83944e-02 , 1.79578e-02 , 1.75378e-02 , 1.71335e-02 , 1.67440e-02 , &
       1.63687e-02 , 1.60069e-02 , 1.56579e-02 , 1.53210e-02 , 1.49958e-02 , &
       1.46815e-02 , 1.43778e-02 , 1.40841e-02 , 1.37999e-02 , 1.35249e-02 , &
       1.32585e-02 , 1.30004e-02 , 1.27502e-02 /)
      absliq1(:,10) = (/ &
! band 10
       7.97040e-02 , 7.63844e-02 , 7.36499e-02 , 7.13525e-02 , 6.93043e-02 , &
       6.72807e-02 , 6.50227e-02 , 6.22395e-02 , 5.86093e-02 , 5.37815e-02 , &
       5.14682e-02 , 4.97214e-02 , 4.77392e-02 , 4.56961e-02 , 4.36858e-02 , &
       4.17569e-02 , 3.99328e-02 , 3.82224e-02 , 3.66265e-02 , 3.51416e-02 , &
       3.37617e-02 , 3.24798e-02 , 3.12887e-02 , 3.01812e-02 , 2.91505e-02 , &
       2.81900e-02 , 2.72939e-02 , 2.64568e-02 , 2.54165e-02 , 2.46832e-02 , &
       2.39783e-02 , 2.33017e-02 , 2.26531e-02 , 2.20314e-02 , 2.14359e-02 , &
       2.08653e-02 , 2.03187e-02 , 1.97947e-02 , 1.92924e-02 , 1.88106e-02 , &
       1.83483e-02 , 1.79043e-02 , 1.74778e-02 , 1.70678e-02 , 1.66735e-02 , &
       1.62941e-02 , 1.59286e-02 , 1.55766e-02 , 1.52371e-02 , 1.49097e-02 , &
       1.45937e-02 , 1.42885e-02 , 1.39936e-02 , 1.37085e-02 , 1.34327e-02 , &
       1.31659e-02 , 1.29075e-02 , 1.26571e-02 /)
      absliq1(:,11) = (/ &
! band 11
       1.49438e-01 , 1.33535e-01 , 1.21542e-01 , 1.11743e-01 , 1.03263e-01 , &
       9.55774e-02 , 8.83382e-02 , 8.12943e-02 , 7.42533e-02 , 6.70609e-02 , &
       6.38761e-02 , 5.97788e-02 , 5.59841e-02 , 5.25318e-02 , 4.94132e-02 , &
       4.66014e-02 , 4.40644e-02 , 4.17706e-02 , 3.96910e-02 , 3.77998e-02 , &
       3.60742e-02 , 3.44947e-02 , 3.30442e-02 , 3.17079e-02 , 3.04730e-02 , &
       2.93283e-02 , 2.82642e-02 , 2.72720e-02 , 2.61789e-02 , 2.53277e-02 , &
       2.45237e-02 , 2.37635e-02 , 2.30438e-02 , 2.23615e-02 , 2.17140e-02 , &
       2.10987e-02 , 2.05133e-02 , 1.99557e-02 , 1.94241e-02 , 1.89166e-02 , &
       1.84317e-02 , 1.79679e-02 , 1.75238e-02 , 1.70983e-02 , 1.66901e-02 , &
       1.62983e-02 , 1.59219e-02 , 1.55599e-02 , 1.52115e-02 , 1.48761e-02 , &
       1.45528e-02 , 1.42411e-02 , 1.39402e-02 , 1.36497e-02 , 1.33690e-02 , &
       1.30976e-02 , 1.28351e-02 , 1.25810e-02 /)
      absliq1(:,12) = (/ &
! band 12
       3.71985e-02 , 3.88586e-02 , 3.99070e-02 , 4.04351e-02 , 4.04610e-02 , &
       3.99834e-02 , 3.89953e-02 , 3.74886e-02 , 3.54551e-02 , 3.28870e-02 , &
       3.32576e-02 , 3.22444e-02 , 3.12384e-02 , 3.02584e-02 , 2.93146e-02 , &
       2.84120e-02 , 2.75525e-02 , 2.67361e-02 , 2.59618e-02 , 2.52280e-02 , &
       2.45327e-02 , 2.38736e-02 , 2.32487e-02 , 2.26558e-02 , 2.20929e-02 , &
       2.15579e-02 , 2.10491e-02 , 2.05648e-02 , 1.99749e-02 , 1.95704e-02 , &
       1.91731e-02 , 1.87839e-02 , 1.84032e-02 , 1.80315e-02 , 1.76689e-02 , &
       1.73155e-02 , 1.69712e-02 , 1.66362e-02 , 1.63101e-02 , 1.59928e-02 , &
       1.56842e-02 , 1.53840e-02 , 1.50920e-02 , 1.48080e-02 , 1.45318e-02 , &
       1.42631e-02 , 1.40016e-02 , 1.37472e-02 , 1.34996e-02 , 1.32586e-02 , &
       1.30239e-02 , 1.27954e-02 , 1.25728e-02 , 1.23559e-02 , 1.21445e-02 , &
       1.19385e-02 , 1.17376e-02 , 1.15417e-02 /)
      absliq1(:,13) = (/ &
! band 13
       3.11868e-02 , 4.48357e-02 , 4.90224e-02 , 4.96406e-02 , 4.86806e-02 , &
       4.69610e-02 , 4.48630e-02 , 4.25795e-02 , 4.02138e-02 , 3.78236e-02 , &
       3.74266e-02 , 3.60384e-02 , 3.47074e-02 , 3.34434e-02 , 3.22499e-02 , &
       3.11264e-02 , 3.00704e-02 , 2.90784e-02 , 2.81463e-02 , 2.72702e-02 , &
       2.64460e-02 , 2.56698e-02 , 2.49381e-02 , 2.42475e-02 , 2.35948e-02 , &
       2.29774e-02 , 2.23925e-02 , 2.18379e-02 , 2.11793e-02 , 2.07076e-02 , &
       2.02470e-02 , 1.97981e-02 , 1.93613e-02 , 1.89367e-02 , 1.85243e-02 , &
       1.81240e-02 , 1.77356e-02 , 1.73588e-02 , 1.69935e-02 , 1.66392e-02 , &
       1.62956e-02 , 1.59624e-02 , 1.56393e-02 , 1.53259e-02 , 1.50219e-02 , &
       1.47268e-02 , 1.44404e-02 , 1.41624e-02 , 1.38925e-02 , 1.36302e-02 , &
       1.33755e-02 , 1.31278e-02 , 1.28871e-02 , 1.26530e-02 , 1.24253e-02 , &
       1.22038e-02 , 1.19881e-02 , 1.17782e-02 /)
      absliq1(:,14) = (/ &
! band 14
       1.58988e-02 , 3.50652e-02 , 4.00851e-02 , 4.07270e-02 , 3.98101e-02 , &
       3.83306e-02 , 3.66829e-02 , 3.50327e-02 , 3.34497e-02 , 3.19609e-02 , &
       3.13712e-02 , 3.03348e-02 , 2.93415e-02 , 2.83973e-02 , 2.75037e-02 , &
       2.66604e-02 , 2.58654e-02 , 2.51161e-02 , 2.44100e-02 , 2.37440e-02 , &
       2.31154e-02 , 2.25215e-02 , 2.19599e-02 , 2.14282e-02 , 2.09242e-02 , &
       2.04459e-02 , 1.99915e-02 , 1.95594e-02 , 1.90254e-02 , 1.86598e-02 , &
       1.82996e-02 , 1.79455e-02 , 1.75983e-02 , 1.72584e-02 , 1.69260e-02 , &
       1.66013e-02 , 1.62843e-02 , 1.59752e-02 , 1.56737e-02 , 1.53799e-02 , &
       1.50936e-02 , 1.48146e-02 , 1.45429e-02 , 1.42782e-02 , 1.40203e-02 , &
       1.37691e-02 , 1.35243e-02 , 1.32858e-02 , 1.30534e-02 , 1.28270e-02 , &
       1.26062e-02 , 1.23909e-02 , 1.21810e-02 , 1.19763e-02 , 1.17766e-02 , &
       1.15817e-02 , 1.13915e-02 , 1.12058e-02 /)
      absliq1(:,15) = (/ &
! band 15
       5.02079e-03 , 2.17615e-02 , 2.55449e-02 , 2.59484e-02 , 2.53650e-02 , &
       2.45281e-02 , 2.36843e-02 , 2.29159e-02 , 2.22451e-02 , 2.16716e-02 , &
       2.11451e-02 , 2.05817e-02 , 2.00454e-02 , 1.95372e-02 , 1.90567e-02 , &
       1.86028e-02 , 1.81742e-02 , 1.77693e-02 , 1.73866e-02 , 1.70244e-02 , &
       1.66815e-02 , 1.63563e-02 , 1.60477e-02 , 1.57544e-02 , 1.54755e-02 , &
       1.52097e-02 , 1.49564e-02 , 1.47146e-02 , 1.43684e-02 , 1.41728e-02 , &
       1.39762e-02 , 1.37797e-02 , 1.35838e-02 , 1.33891e-02 , 1.31961e-02 , &
       1.30051e-02 , 1.28164e-02 , 1.26302e-02 , 1.24466e-02 , 1.22659e-02 , &
       1.20881e-02 , 1.19131e-02 , 1.17412e-02 , 1.15723e-02 , 1.14063e-02 , &
       1.12434e-02 , 1.10834e-02 , 1.09264e-02 , 1.07722e-02 , 1.06210e-02 , &
       1.04725e-02 , 1.03269e-02 , 1.01839e-02 , 1.00436e-02 , 9.90593e-03 , &
       9.77080e-03 , 9.63818e-03 , 9.50800e-03 /)
      absliq1(:,16) = (/ &
! band 16
       5.64971e-02 , 9.04736e-02 , 8.11726e-02 , 7.05450e-02 , 6.20052e-02 , &
       5.54286e-02 , 5.03503e-02 , 4.63791e-02 , 4.32290e-02 , 4.06959e-02 , &
       3.74690e-02 , 3.52964e-02 , 3.33799e-02 , 3.16774e-02 , 3.01550e-02 , &
       2.87856e-02 , 2.75474e-02 , 2.64223e-02 , 2.53953e-02 , 2.44542e-02 , &
       2.35885e-02 , 2.27894e-02 , 2.20494e-02 , 2.13622e-02 , 2.07222e-02 , &
       2.01246e-02 , 1.95654e-02 , 1.90408e-02 , 1.84398e-02 , 1.80021e-02 , &
       1.75816e-02 , 1.71775e-02 , 1.67889e-02 , 1.64152e-02 , 1.60554e-02 , &
       1.57089e-02 , 1.53751e-02 , 1.50531e-02 , 1.47426e-02 , 1.44428e-02 , &
       1.41532e-02 , 1.38734e-02 , 1.36028e-02 , 1.33410e-02 , 1.30875e-02 , &
       1.28420e-02 , 1.26041e-02 , 1.23735e-02 , 1.21497e-02 , 1.19325e-02 , &
       1.17216e-02 , 1.15168e-02 , 1.13177e-02 , 1.11241e-02 , 1.09358e-02 , &
       1.07525e-02 , 1.05741e-02 , 1.04003e-02 /)


      absice4(:,  1) = (/ &
        &  0.6888440E-01 ,  0.6950670E-01 ,  0.6985297E-01 ,  0.7011411E-01 ,  0.7006932E-01 ,&
        &  0.6964392E-01 ,  0.6890915E-01 ,  0.6795748E-01 ,  0.6685709E-01 ,  0.6565582E-01 ,&
        &  0.6438961E-01 ,  0.6308378E-01 ,  0.6175803E-01 ,  0.6042715E-01 ,  0.5910178E-01 ,&
        &  0.5778937E-01 ,  0.5649547E-01 ,  0.5522512E-01 ,  0.5397782E-01 ,  0.5275840E-01 ,&
        &  0.5156704E-01 ,  0.5040565E-01 ,  0.4927261E-01 ,  0.4816831E-01 ,  0.4709510E-01 ,&
        &  0.4605071E-01 ,  0.4503559E-01 ,  0.4404905E-01 ,  0.4309146E-01 ,  0.4216224E-01 ,&
        &  0.4126133E-01 ,  0.4038681E-01 ,  0.3953866E-01 ,  0.3871694E-01 ,  0.3792092E-01 ,&
        &  0.3714872E-01 ,  0.3640130E-01 ,  0.3567646E-01 ,  0.3497424E-01 ,  0.3429416E-01 ,&
        &  0.3363520E-01 ,  0.3299676E-01 ,  0.3237776E-01 ,  0.3177813E-01 ,  0.3119717E-01 ,&
        &  0.3063310E-01 ,  0.3008657E-01 ,  0.2955659E-01 ,  0.2904246E-01 ,  0.2854360E-01 ,&
        &  0.2805947E-01 ,  0.2758961E-01 ,  0.2713342E-01 ,  0.2669037E-01 ,  0.2625994E-01 ,&
        &  0.2584212E-01 ,  0.2543570E-01 ,  0.2504040E-01 ,  0.2465650E-01 ,  0.2428280E-01 ,&
        &  0.2391933E-01 ,  0.2356580E-01 ,  0.2322152E-01 ,  0.2288620E-01 ,  0.2255971E-01 ,&
        &  0.2224175E-01 ,  0.2193226E-01 ,  0.2163004E-01 ,  0.2133575E-01 ,  0.2104884E-01 ,&
        &  0.2076901E-01 ,  0.2049589E-01 ,  0.2022937E-01 ,  0.1996956E-01 ,  0.1971553E-01 ,&
        &  0.1946782E-01 ,  0.1922561E-01 ,  0.1898925E-01 ,  0.1875803E-01 ,  0.1853229E-01 ,&
        &  0.1831155E-01 ,  0.1809566E-01 ,  0.1788444E-01 ,  0.1767801E-01 ,  0.1747603E-01 ,&
        &  0.1727850E-01 ,  0.1708492E-01 ,  0.1689560E-01 ,  0.1671009E-01 ,  0.1652847E-01 ,&
        &  0.1635060E-01 ,  0.1617640E-01 ,  0.1600567E-01 ,  0.1583831E-01 ,  0.1567423E-01 ,&
        &  0.1551341E-01 ,  0.1535575E-01 ,  0.1520111E-01 ,  0.1504947E-01 ,  0.1490059E-01 ,&
        &  0.1475464E-01 ,  0.1461139E-01 ,  0.1447078E-01 ,  0.1433284E-01 ,  0.1419729E-01 ,&
        &  0.1406420E-01 ,  0.1393352E-01 ,  0.1380513E-01 ,  0.1367904E-01 ,  0.1355515E-01 ,&
        &  0.1343340E-01 ,  0.1331372E-01 ,  0.1319610E-01 ,  0.1308044E-01 ,  0.1296682E-01 ,&
        &  0.1285500E-01 ,  0.1274509E-01 ,  0.1263696E-01 ,  0.1253059E-01 ,  0.1242593E-01 ,&
        &  0.1232298E-01 ,  0.1222168E-01 ,  0.1212195E-01 ,  0.1202379E-01 ,  0.1192715E-01 ,&
        &  0.1183206E-01 ,  0.1173839E-01 ,  0.1164616E-01 ,  0.1155532E-01 ,  0.1146587E-01 ,&
        &  0.1137773E-01 ,  0.1129092E-01 ,  0.1120538E-01 ,  0.1112109E-01 ,  0.1103801E-01 ,&
        &  0.1095617E-01 ,  0.1087546E-01 ,  0.1079591E-01 ,  0.1071751E-01 ,  0.1064021E-01 ,&
        &  0.1056396E-01 ,  0.1048879E-01 ,  0.1041466E-01 ,  0.1034153E-01 ,  0.1026945E-01 ,&
        &  0.1019829E-01 ,  0.1012810E-01 ,  0.1005882E-01 ,  0.9990477E-02 ,  0.9923048E-02 ,&
        &  0.9856497E-02 ,  0.9790796E-02 ,  0.9725967E-02 ,  0.9661886E-02 ,  0.9598750E-02 ,&
        &  0.9536333E-02 ,  0.9474728E-02 ,  0.9413908E-02 ,  0.9353837E-02 ,  0.9294510E-02 ,&
        &  0.9235926E-02 ,  0.9178066E-02 ,  0.9120892E-02 ,  0.9064428E-02 ,  0.9008617E-02 ,&
        &  0.8953504E-02 ,  0.8899035E-02 ,  0.8845209E-02 ,  0.8792018E-02 ,  0.8739441E-02 ,&
        &  0.8687489E-02 ,  0.8636141E-02 ,  0.8585371E-02 ,  0.8535195E-02 ,  0.8485589E-02 ,&
        &  0.8436535E-02 ,  0.8388070E-02 ,  0.8340090E-02 ,  0.8292706E-02 ,  0.8245819E-02 ,&
        &  0.8199411E-02 ,  0.8153542E-02 ,  0.8108193E-02 ,  0.8063350E-02 ,  0.8018982E-02 ,&
        &  0.7974982E-02 ,  0.7931642E-02 ,  0.7888632E-02 ,  0.7846153E-02 ,  0.7804049E-02 ,&
        &  0.7762399E-02 ,  0.7721156E-02 ,  0.7680319E-02 ,  0.7640009E-02 ,  0.7600155E-02 ,&
        &  0.7560589E-02 ,  0.7521371E-02 ,  0.7482661E-02 ,  0.7444382E-02 ,  0.7406328E-02 /)
      absice4(:,  2) = (/ &
        &  0.1808749E-01 ,  0.1932276E-01 ,  0.2023977E-01 ,  0.2114905E-01 ,  0.2201138E-01 ,&
        &  0.2278221E-01 ,  0.2344802E-01 ,  0.2401308E-01 ,  0.2448537E-01 ,  0.2487125E-01 ,&
        &  0.2517609E-01 ,  0.2540326E-01 ,  0.2555691E-01 ,  0.2564221E-01 ,  0.2566490E-01 ,&
        &  0.2563131E-01 ,  0.2554830E-01 ,  0.2542333E-01 ,  0.2526240E-01 ,  0.2507216E-01 ,&
        &  0.2485834E-01 ,  0.2462615E-01 ,  0.2437942E-01 ,  0.2412181E-01 ,  0.2385707E-01 ,&
        &  0.2358728E-01 ,  0.2331461E-01 ,  0.2304081E-01 ,  0.2276729E-01 ,  0.2249523E-01 ,&
        &  0.2222565E-01 ,  0.2195872E-01 ,  0.2169522E-01 ,  0.2143569E-01 ,  0.2118043E-01 ,&
        &  0.2092925E-01 ,  0.2068280E-01 ,  0.2044069E-01 ,  0.2020311E-01 ,  0.1997030E-01 ,&
        &  0.1974206E-01 ,  0.1951841E-01 ,  0.1929905E-01 ,  0.1908421E-01 ,  0.1887382E-01 ,&
        &  0.1866731E-01 ,  0.1846508E-01 ,  0.1826695E-01 ,  0.1807275E-01 ,  0.1788237E-01 ,&
        &  0.1769572E-01 ,  0.1751278E-01 ,  0.1733340E-01 ,  0.1715745E-01 ,  0.1698485E-01 ,&
        &  0.1681570E-01 ,  0.1664962E-01 ,  0.1648652E-01 ,  0.1632667E-01 ,  0.1616967E-01 ,&
        &  0.1601550E-01 ,  0.1586427E-01 ,  0.1571565E-01 ,  0.1556963E-01 ,  0.1542618E-01 ,&
        &  0.1528533E-01 ,  0.1514708E-01 ,  0.1501091E-01 ,  0.1487722E-01 ,  0.1474583E-01 ,&
        &  0.1461663E-01 ,  0.1448955E-01 ,  0.1436455E-01 ,  0.1424175E-01 ,  0.1412078E-01 ,&
        &  0.1400193E-01 ,  0.1388483E-01 ,  0.1376973E-01 ,  0.1365632E-01 ,  0.1354478E-01 ,&
        &  0.1343496E-01 ,  0.1332682E-01 ,  0.1322025E-01 ,  0.1311538E-01 ,  0.1301213E-01 ,&
        &  0.1291048E-01 ,  0.1281017E-01 ,  0.1271150E-01 ,  0.1261414E-01 ,  0.1251826E-01 ,&
        &  0.1242377E-01 ,  0.1233068E-01 ,  0.1223888E-01 ,  0.1214835E-01 ,  0.1205910E-01 ,&
        &  0.1197108E-01 ,  0.1188429E-01 ,  0.1179873E-01 ,  0.1171433E-01 ,  0.1163102E-01 ,&
        &  0.1154892E-01 ,  0.1146791E-01 ,  0.1138796E-01 ,  0.1130913E-01 ,  0.1123124E-01 ,&
        &  0.1115440E-01 ,  0.1107857E-01 ,  0.1100370E-01 ,  0.1092982E-01 ,  0.1085688E-01 ,&
        &  0.1078486E-01 ,  0.1071372E-01 ,  0.1064350E-01 ,  0.1057410E-01 ,  0.1050566E-01 ,&
        &  0.1043798E-01 ,  0.1037117E-01 ,  0.1030515E-01 ,  0.1023992E-01 ,  0.1017549E-01 ,&
        &  0.1011183E-01 ,  0.1004894E-01 ,  0.9986768E-02 ,  0.9925325E-02 ,  0.9864601E-02 ,&
        &  0.9804619E-02 ,  0.9745291E-02 ,  0.9686657E-02 ,  0.9628678E-02 ,  0.9571399E-02 ,&
        &  0.9514725E-02 ,  0.9458707E-02 ,  0.9403317E-02 ,  0.9348548E-02 ,  0.9294357E-02 ,&
        &  0.9240794E-02 ,  0.9187792E-02 ,  0.9135368E-02 ,  0.9083532E-02 ,  0.9032273E-02 ,&
        &  0.8981521E-02 ,  0.8931329E-02 ,  0.8881679E-02 ,  0.8832537E-02 ,  0.8783953E-02 ,&
        &  0.8735837E-02 ,  0.8688229E-02 ,  0.8641113E-02 ,  0.8594473E-02 ,  0.8548328E-02 ,&
        &  0.8502663E-02 ,  0.8457437E-02 ,  0.8412691E-02 ,  0.8368331E-02 ,  0.8324496E-02 ,&
        &  0.8281046E-02 ,  0.8238045E-02 ,  0.8195484E-02 ,  0.8153325E-02 ,  0.8111581E-02 ,&
        &  0.8070248E-02 ,  0.8029317E-02 ,  0.7988783E-02 ,  0.7948641E-02 ,  0.7908863E-02 ,&
        &  0.7869483E-02 ,  0.7830469E-02 ,  0.7791823E-02 ,  0.7753538E-02 ,  0.7715614E-02 ,&
        &  0.7678050E-02 ,  0.7640833E-02 ,  0.7603942E-02 ,  0.7567415E-02 ,  0.7531212E-02 ,&
        &  0.7495315E-02 ,  0.7459802E-02 ,  0.7424545E-02 ,  0.7389668E-02 ,  0.7355077E-02 ,&
        &  0.7320753E-02 ,  0.7286766E-02 ,  0.7253100E-02 ,  0.7219734E-02 ,  0.7186655E-02 ,&
        &  0.7153776E-02 ,  0.7121335E-02 ,  0.7089079E-02 ,  0.7057161E-02 ,  0.7025448E-02 ,&
        &  0.6994023E-02 ,  0.6962844E-02 ,  0.6931921E-02 ,  0.6901341E-02 ,  0.6871042E-02 ,&
        &  0.6840906E-02 ,  0.6810980E-02 ,  0.6781394E-02 ,  0.6752098E-02 ,  0.6722898E-02 /)
      absice4(:,  3) = (/ &
        &  0.6134844E-01 ,  0.6614562E-01 ,  0.6981280E-01 ,  0.7276729E-01 ,  0.7480936E-01 ,&
        &  0.7600386E-01 ,  0.7649013E-01 ,  0.7639786E-01 ,  0.7582567E-01 ,  0.7485179E-01 ,&
        &  0.7354883E-01 ,  0.7198476E-01 ,  0.7022718E-01 ,  0.6833836E-01 ,  0.6637321E-01 ,&
        &  0.6437675E-01 ,  0.6238502E-01 ,  0.6042662E-01 ,  0.5851643E-01 ,  0.5667127E-01 ,&
        &  0.5489796E-01 ,  0.5320195E-01 ,  0.5158193E-01 ,  0.5003763E-01 ,  0.4857029E-01 ,&
        &  0.4717415E-01 ,  0.4584703E-01 ,  0.4458477E-01 ,  0.4338476E-01 ,  0.4224324E-01 ,&
        &  0.4115709E-01 ,  0.4012140E-01 ,  0.3913360E-01 ,  0.3819148E-01 ,  0.3729217E-01 ,&
        &  0.3643164E-01 ,  0.3560931E-01 ,  0.3482125E-01 ,  0.3406615E-01 ,  0.3334230E-01 ,&
        &  0.3264759E-01 ,  0.3198040E-01 ,  0.3133877E-01 ,  0.3072184E-01 ,  0.3012829E-01 ,&
        &  0.2955574E-01 ,  0.2900426E-01 ,  0.2847244E-01 ,  0.2795912E-01 ,  0.2746342E-01 ,&
        &  0.2698442E-01 ,  0.2652142E-01 ,  0.2607354E-01 ,  0.2564005E-01 ,  0.2522026E-01 ,&
        &  0.2481393E-01 ,  0.2441975E-01 ,  0.2403736E-01 ,  0.2366677E-01 ,  0.2330685E-01 ,&
        &  0.2295745E-01 ,  0.2261820E-01 ,  0.2228836E-01 ,  0.2196760E-01 ,  0.2165572E-01 ,&
        &  0.2135240E-01 ,  0.2105746E-01 ,  0.2076980E-01 ,  0.2048994E-01 ,  0.2021732E-01 ,&
        &  0.1995168E-01 ,  0.1969257E-01 ,  0.1943992E-01 ,  0.1919375E-01 ,  0.1895319E-01 ,&
        &  0.1871875E-01 ,  0.1848959E-01 ,  0.1826607E-01 ,  0.1804747E-01 ,  0.1783412E-01 ,&
        &  0.1762554E-01 ,  0.1742162E-01 ,  0.1722211E-01 ,  0.1702717E-01 ,  0.1683648E-01 ,&
        &  0.1665000E-01 ,  0.1646723E-01 ,  0.1628851E-01 ,  0.1611341E-01 ,  0.1594198E-01 ,&
        &  0.1577409E-01 ,  0.1560962E-01 ,  0.1544846E-01 ,  0.1529047E-01 ,  0.1513555E-01 ,&
        &  0.1498369E-01 ,  0.1483481E-01 ,  0.1468877E-01 ,  0.1454555E-01 ,  0.1440490E-01 ,&
        &  0.1426700E-01 ,  0.1413165E-01 ,  0.1399878E-01 ,  0.1386839E-01 ,  0.1374025E-01 ,&
        &  0.1361442E-01 ,  0.1349083E-01 ,  0.1336940E-01 ,  0.1325011E-01 ,  0.1313288E-01 ,&
        &  0.1301767E-01 ,  0.1290437E-01 ,  0.1279305E-01 ,  0.1268351E-01 ,  0.1257590E-01 ,&
        &  0.1246997E-01 ,  0.1236584E-01 ,  0.1226335E-01 ,  0.1216252E-01 ,  0.1206332E-01 ,&
        &  0.1196568E-01 ,  0.1186960E-01 ,  0.1177498E-01 ,  0.1168185E-01 ,  0.1159014E-01 ,&
        &  0.1149985E-01 ,  0.1141093E-01 ,  0.1132334E-01 ,  0.1123704E-01 ,  0.1115207E-01 ,&
        &  0.1106831E-01 ,  0.1098580E-01 ,  0.1090448E-01 ,  0.1082434E-01 ,  0.1074534E-01 ,&
        &  0.1066747E-01 ,  0.1059068E-01 ,  0.1051498E-01 ,  0.1044035E-01 ,  0.1036678E-01 ,&
        &  0.1029417E-01 ,  0.1022258E-01 ,  0.1015197E-01 ,  0.1008230E-01 ,  0.1001361E-01 ,&
        &  0.9945779E-02 ,  0.9878875E-02 ,  0.9812825E-02 ,  0.9747642E-02 ,  0.9683326E-02 ,&
        &  0.9619834E-02 ,  0.9557145E-02 ,  0.9495276E-02 ,  0.9434108E-02 ,  0.9373835E-02 ,&
        &  0.9314229E-02 ,  0.9255388E-02 ,  0.9197296E-02 ,  0.9139908E-02 ,  0.9083220E-02 ,&
        &  0.9027225E-02 ,  0.8971911E-02 ,  0.8917253E-02 ,  0.8863264E-02 ,  0.8809891E-02 ,&
        &  0.8757176E-02 ,  0.8705065E-02 ,  0.8653563E-02 ,  0.8602656E-02 ,  0.8552338E-02 ,&
        &  0.8502601E-02 ,  0.8453442E-02 ,  0.8404816E-02 ,  0.8356771E-02 ,  0.8309252E-02 ,&
        &  0.8262245E-02 ,  0.8215820E-02 ,  0.8169831E-02 ,  0.8124417E-02 ,  0.8079469E-02 ,&
        &  0.8034971E-02 ,  0.7990990E-02 ,  0.7947500E-02 ,  0.7904486E-02 ,  0.7861921E-02 ,&
        &  0.7819697E-02 ,  0.7778111E-02 ,  0.7736838E-02 ,  0.7696060E-02 ,  0.7655641E-02 ,&
        &  0.7615650E-02 ,  0.7576041E-02 ,  0.7536822E-02 ,  0.7498102E-02 ,  0.7459822E-02 ,&
        &  0.7421801E-02 ,  0.7384120E-02 ,  0.7346914E-02 ,  0.7310124E-02 ,  0.7273538E-02 /)
      absice4(:,  4) = (/ &
        &  0.1326735E+00 ,  0.1417908E+00 ,  0.1468841E+00 ,  0.1485195E+00 ,  0.1474047E+00 ,&
        &  0.1443156E+00 ,  0.1399179E+00 ,  0.1347081E+00 ,  0.1290312E+00 ,  0.1231323E+00 ,&
        &  0.1171963E+00 ,  0.1113586E+00 ,  0.1057195E+00 ,  0.1003463E+00 ,  0.9527878E-01 ,&
        &  0.9053445E-01 ,  0.8611556E-01 ,  0.8201642E-01 ,  0.7821326E-01 ,  0.7469463E-01 ,&
        &  0.7143734E-01 ,  0.6842215E-01 ,  0.6562440E-01 ,  0.6302506E-01 ,  0.6061122E-01 ,&
        &  0.5836153E-01 ,  0.5626275E-01 ,  0.5430037E-01 ,  0.5246392E-01 ,  0.5074203E-01 ,&
        &  0.4912554E-01 ,  0.4760348E-01 ,  0.4616864E-01 ,  0.4481519E-01 ,  0.4353668E-01 ,&
        &  0.4232533E-01 ,  0.4117850E-01 ,  0.4008913E-01 ,  0.3905430E-01 ,  0.3807024E-01 ,&
        &  0.3713309E-01 ,  0.3623973E-01 ,  0.3538680E-01 ,  0.3457227E-01 ,  0.3379383E-01 ,&
        &  0.3304768E-01 ,  0.3233342E-01 ,  0.3164866E-01 ,  0.3099151E-01 ,  0.3036043E-01 ,&
        &  0.2975389E-01 ,  0.2917062E-01 ,  0.2860922E-01 ,  0.2806853E-01 ,  0.2754736E-01 ,&
        &  0.2704519E-01 ,  0.2656024E-01 ,  0.2609174E-01 ,  0.2563968E-01 ,  0.2520232E-01 ,&
        &  0.2477946E-01 ,  0.2437042E-01 ,  0.2397421E-01 ,  0.2359031E-01 ,  0.2321836E-01 ,&
        &  0.2285782E-01 ,  0.2250840E-01 ,  0.2216874E-01 ,  0.2183930E-01 ,  0.2151936E-01 ,&
        &  0.2120856E-01 ,  0.2090627E-01 ,  0.2061235E-01 ,  0.2032674E-01 ,  0.2004843E-01 ,&
        &  0.1977785E-01 ,  0.1951409E-01 ,  0.1925742E-01 ,  0.1900704E-01 ,  0.1876324E-01 ,&
        &  0.1852544E-01 ,  0.1829345E-01 ,  0.1806702E-01 ,  0.1784624E-01 ,  0.1763069E-01 ,&
        &  0.1742035E-01 ,  0.1721461E-01 ,  0.1701383E-01 ,  0.1681748E-01 ,  0.1662559E-01 ,&
        &  0.1643801E-01 ,  0.1625460E-01 ,  0.1607517E-01 ,  0.1589956E-01 ,  0.1572767E-01 ,&
        &  0.1555945E-01 ,  0.1539479E-01 ,  0.1523350E-01 ,  0.1507558E-01 ,  0.1492072E-01 ,&
        &  0.1476912E-01 ,  0.1462052E-01 ,  0.1447485E-01 ,  0.1433208E-01 ,  0.1419197E-01 ,&
        &  0.1405457E-01 ,  0.1391979E-01 ,  0.1378752E-01 ,  0.1365776E-01 ,  0.1353039E-01 ,&
        &  0.1340534E-01 ,  0.1328253E-01 ,  0.1316198E-01 ,  0.1304352E-01 ,  0.1292724E-01 ,&
        &  0.1281291E-01 ,  0.1270062E-01 ,  0.1259025E-01 ,  0.1248176E-01 ,  0.1237511E-01 ,&
        &  0.1227026E-01 ,  0.1216717E-01 ,  0.1206576E-01 ,  0.1196601E-01 ,  0.1186788E-01 ,&
        &  0.1177137E-01 ,  0.1167638E-01 ,  0.1158290E-01 ,  0.1149088E-01 ,  0.1140035E-01 ,&
        &  0.1131117E-01 ,  0.1122339E-01 ,  0.1113696E-01 ,  0.1105183E-01 ,  0.1096797E-01 ,&
        &  0.1088538E-01 ,  0.1080400E-01 ,  0.1072382E-01 ,  0.1064483E-01 ,  0.1056701E-01 ,&
        &  0.1049026E-01 ,  0.1041465E-01 ,  0.1034010E-01 ,  0.1026660E-01 ,  0.1019418E-01 ,&
        &  0.1012272E-01 ,  0.1005226E-01 ,  0.9982754E-02 ,  0.9914191E-02 ,  0.9846580E-02 ,&
        &  0.9779881E-02 ,  0.9714044E-02 ,  0.9649111E-02 ,  0.9584956E-02 ,  0.9521768E-02 ,&
        &  0.9459312E-02 ,  0.9397696E-02 ,  0.9336874E-02 ,  0.9276840E-02 ,  0.9217549E-02 ,&
        &  0.9159023E-02 ,  0.9101236E-02 ,  0.9044155E-02 ,  0.8987799E-02 ,  0.8932117E-02 ,&
        &  0.8877140E-02 ,  0.8822812E-02 ,  0.8769146E-02 ,  0.8716122E-02 ,  0.8663736E-02 ,&
        &  0.8611974E-02 ,  0.8560832E-02 ,  0.8510272E-02 ,  0.8460325E-02 ,  0.8410953E-02 ,&
        &  0.8362130E-02 ,  0.8313920E-02 ,  0.8266191E-02 ,  0.8219064E-02 ,  0.8172440E-02 ,&
        &  0.8126310E-02 ,  0.8080724E-02 ,  0.8035664E-02 ,  0.7991114E-02 ,  0.7947044E-02 ,&
        &  0.7903340E-02 ,  0.7860309E-02 ,  0.7817612E-02 ,  0.7775442E-02 ,  0.7733653E-02 ,&
        &  0.7692324E-02 ,  0.7651406E-02 ,  0.7610899E-02 ,  0.7570923E-02 ,  0.7531399E-02 ,&
        &  0.7492168E-02 ,  0.7453290E-02 ,  0.7414924E-02 ,  0.7376985E-02 ,  0.7339269E-02 /)
      absice4(:,  5) = (/ &
        &  0.3024995E+00 ,  0.3022310E+00 ,  0.2871740E+00 ,  0.2668229E+00 ,  0.2454907E+00 ,&
        &  0.2248826E+00 ,  0.2058085E+00 ,  0.1885392E+00 ,  0.1730567E+00 ,  0.1592368E+00 ,&
        &  0.1469269E+00 ,  0.1359729E+00 ,  0.1262249E+00 ,  0.1175487E+00 ,  0.1098143E+00 ,&
        &  0.1029054E+00 ,  0.9671755E-01 ,  0.9116329E-01 ,  0.8615229E-01 ,  0.8162510E-01 ,&
        &  0.7751813E-01 ,  0.7378272E-01 ,  0.7036919E-01 ,  0.6724078E-01 ,  0.6436986E-01 ,&
        &  0.6172230E-01 ,  0.5927608E-01 ,  0.5700852E-01 ,  0.5490293E-01 ,  0.5294263E-01 ,&
        &  0.5111478E-01 ,  0.4940378E-01 ,  0.4780013E-01 ,  0.4629565E-01 ,  0.4488094E-01 ,&
        &  0.4354687E-01 ,  0.4228956E-01 ,  0.4109976E-01 ,  0.3997432E-01 ,  0.3890744E-01 ,&
        &  0.3789524E-01 ,  0.3693316E-01 ,  0.3601765E-01 ,  0.3514581E-01 ,  0.3431485E-01 ,&
        &  0.3352055E-01 ,  0.3276202E-01 ,  0.3203660E-01 ,  0.3134199E-01 ,  0.3067640E-01 ,&
        &  0.3003801E-01 ,  0.2942533E-01 ,  0.2883671E-01 ,  0.2827086E-01 ,  0.2772636E-01 ,&
        &  0.2720262E-01 ,  0.2669759E-01 ,  0.2621048E-01 ,  0.2574111E-01 ,  0.2528767E-01 ,&
        &  0.2484986E-01 ,  0.2442688E-01 ,  0.2401769E-01 ,  0.2362166E-01 ,  0.2323842E-01 ,&
        &  0.2286735E-01 ,  0.2250810E-01 ,  0.2215925E-01 ,  0.2182121E-01 ,  0.2149325E-01 ,&
        &  0.2117495E-01 ,  0.2086561E-01 ,  0.2056512E-01 ,  0.2027336E-01 ,  0.1998927E-01 ,&
        &  0.1971332E-01 ,  0.1944449E-01 ,  0.1918309E-01 ,  0.1892825E-01 ,  0.1868030E-01 ,&
        &  0.1843861E-01 ,  0.1820295E-01 ,  0.1797310E-01 ,  0.1774910E-01 ,  0.1753055E-01 ,&
        &  0.1731737E-01 ,  0.1710900E-01 ,  0.1690573E-01 ,  0.1670705E-01 ,  0.1651299E-01 ,&
        &  0.1632337E-01 ,  0.1613805E-01 ,  0.1595683E-01 ,  0.1577956E-01 ,  0.1560611E-01 ,&
        &  0.1543643E-01 ,  0.1527041E-01 ,  0.1510785E-01 ,  0.1494876E-01 ,  0.1479280E-01 ,&
        &  0.1464017E-01 ,  0.1449062E-01 ,  0.1434407E-01 ,  0.1420049E-01 ,  0.1405962E-01 ,&
        &  0.1392153E-01 ,  0.1378610E-01 ,  0.1365324E-01 ,  0.1352294E-01 ,  0.1339508E-01 ,&
        &  0.1326958E-01 ,  0.1314636E-01 ,  0.1302542E-01 ,  0.1290664E-01 ,  0.1279006E-01 ,&
        &  0.1267547E-01 ,  0.1256297E-01 ,  0.1245238E-01 ,  0.1234372E-01 ,  0.1223693E-01 ,&
        &  0.1213196E-01 ,  0.1202877E-01 ,  0.1192729E-01 ,  0.1182750E-01 ,  0.1172935E-01 ,&
        &  0.1163282E-01 ,  0.1153784E-01 ,  0.1144440E-01 ,  0.1135242E-01 ,  0.1126195E-01 ,&
        &  0.1117285E-01 ,  0.1108516E-01 ,  0.1099883E-01 ,  0.1091382E-01 ,  0.1083009E-01 ,&
        &  0.1074766E-01 ,  0.1066642E-01 ,  0.1058641E-01 ,  0.1050760E-01 ,  0.1042995E-01 ,&
        &  0.1035341E-01 ,  0.1027799E-01 ,  0.1020366E-01 ,  0.1013037E-01 ,  0.1005817E-01 ,&
        &  0.9986949E-02 ,  0.9916727E-02 ,  0.9847458E-02 ,  0.9779151E-02 ,  0.9711792E-02 ,&
        &  0.9645349E-02 ,  0.9579788E-02 ,  0.9515127E-02 ,  0.9451241E-02 ,  0.9388327E-02 ,&
        &  0.9326161E-02 ,  0.9264823E-02 ,  0.9204305E-02 ,  0.9144551E-02 ,  0.9085567E-02 ,&
        &  0.9027331E-02 ,  0.8969852E-02 ,  0.8913075E-02 ,  0.8857019E-02 ,  0.8801642E-02 ,&
        &  0.8746976E-02 ,  0.8692966E-02 ,  0.8639613E-02 ,  0.8586907E-02 ,  0.8534838E-02 ,&
        &  0.8483395E-02 ,  0.8432568E-02 ,  0.8382322E-02 ,  0.8332702E-02 ,  0.8283645E-02 ,&
        &  0.8235149E-02 ,  0.8187264E-02 ,  0.8139856E-02 ,  0.8093054E-02 ,  0.8046762E-02 ,&
        &  0.8000951E-02 ,  0.7955688E-02 ,  0.7910957E-02 ,  0.7866731E-02 ,  0.7822983E-02 ,&
        &  0.7779608E-02 ,  0.7736897E-02 ,  0.7694529E-02 ,  0.7652683E-02 ,  0.7611223E-02 ,&
        &  0.7570220E-02 ,  0.7529625E-02 ,  0.7489445E-02 ,  0.7449784E-02 ,  0.7410586E-02 ,&
        &  0.7371680E-02 ,  0.7333126E-02 ,  0.7295078E-02 ,  0.7257461E-02 ,  0.7220066E-02 /)
      absice4(:,  6) = (/ &
        &  0.2880479E+00 ,  0.2709946E+00 ,  0.2465274E+00 ,  0.2226647E+00 ,  0.2012633E+00 ,&
        &  0.1825816E+00 ,  0.1664524E+00 ,  0.1525213E+00 ,  0.1404165E+00 ,  0.1298243E+00 ,&
        &  0.1205029E+00 ,  0.1122591E+00 ,  0.1049403E+00 ,  0.9841985E-01 ,  0.9259034E-01 ,&
        &  0.8735970E-01 ,  0.8264969E-01 ,  0.7839588E-01 ,  0.7453251E-01 ,  0.7101811E-01 ,&
        &  0.6780819E-01 ,  0.6486867E-01 ,  0.6216424E-01 ,  0.5966881E-01 ,  0.5736412E-01 ,&
        &  0.5522554E-01 ,  0.5323753E-01 ,  0.5138391E-01 ,  0.4965331E-01 ,  0.4803360E-01 ,&
        &  0.4651528E-01 ,  0.4508732E-01 ,  0.4374247E-01 ,  0.4247493E-01 ,  0.4127821E-01 ,&
        &  0.4014480E-01 ,  0.3907220E-01 ,  0.3805358E-01 ,  0.3708605E-01 ,  0.3616610E-01 ,&
        &  0.3529001E-01 ,  0.3445483E-01 ,  0.3365739E-01 ,  0.3289579E-01 ,  0.3216780E-01 ,&
        &  0.3146994E-01 ,  0.3080171E-01 ,  0.3016099E-01 ,  0.2954597E-01 ,  0.2895516E-01 ,&
        &  0.2838718E-01 ,  0.2784088E-01 ,  0.2731489E-01 ,  0.2680814E-01 ,  0.2631953E-01 ,&
        &  0.2584863E-01 ,  0.2539369E-01 ,  0.2495406E-01 ,  0.2452972E-01 ,  0.2411903E-01 ,&
        &  0.2372185E-01 ,  0.2333751E-01 ,  0.2296511E-01 ,  0.2260415E-01 ,  0.2225430E-01 ,&
        &  0.2191508E-01 ,  0.2158622E-01 ,  0.2126640E-01 ,  0.2095613E-01 ,  0.2065472E-01 ,&
        &  0.2036179E-01 ,  0.2007682E-01 ,  0.1979962E-01 ,  0.1953018E-01 ,  0.1926754E-01 ,&
        &  0.1901214E-01 ,  0.1876306E-01 ,  0.1852062E-01 ,  0.1828403E-01 ,  0.1805360E-01 ,&
        &  0.1782877E-01 ,  0.1760937E-01 ,  0.1739515E-01 ,  0.1718621E-01 ,  0.1698218E-01 ,&
        &  0.1678300E-01 ,  0.1658813E-01 ,  0.1639790E-01 ,  0.1621181E-01 ,  0.1602992E-01 ,&
        &  0.1585204E-01 ,  0.1567806E-01 ,  0.1550782E-01 ,  0.1534116E-01 ,  0.1517798E-01 ,&
        &  0.1501824E-01 ,  0.1486185E-01 ,  0.1470862E-01 ,  0.1455855E-01 ,  0.1441136E-01 ,&
        &  0.1426720E-01 ,  0.1412590E-01 ,  0.1398731E-01 ,  0.1385149E-01 ,  0.1371813E-01 ,&
        &  0.1358732E-01 ,  0.1345901E-01 ,  0.1333303E-01 ,  0.1320942E-01 ,  0.1308805E-01 ,&
        &  0.1296889E-01 ,  0.1285181E-01 ,  0.1273686E-01 ,  0.1262389E-01 ,  0.1251298E-01 ,&
        &  0.1240389E-01 ,  0.1229674E-01 ,  0.1219140E-01 ,  0.1208781E-01 ,  0.1198599E-01 ,&
        &  0.1188585E-01 ,  0.1178738E-01 ,  0.1169047E-01 ,  0.1159515E-01 ,  0.1150137E-01 ,&
        &  0.1140910E-01 ,  0.1131828E-01 ,  0.1122889E-01 ,  0.1114087E-01 ,  0.1105424E-01 ,&
        &  0.1096892E-01 ,  0.1088492E-01 ,  0.1080218E-01 ,  0.1072069E-01 ,  0.1064038E-01 ,&
        &  0.1056130E-01 ,  0.1048335E-01 ,  0.1040653E-01 ,  0.1033085E-01 ,  0.1025627E-01 ,&
        &  0.1018272E-01 ,  0.1011023E-01 ,  0.1003876E-01 ,  0.9968274E-02 ,  0.9898823E-02 ,&
        &  0.9830278E-02 ,  0.9762678E-02 ,  0.9695984E-02 ,  0.9630197E-02 ,  0.9565303E-02 ,&
        &  0.9501283E-02 ,  0.9438083E-02 ,  0.9375740E-02 ,  0.9314126E-02 ,  0.9253441E-02 ,&
        &  0.9193447E-02 ,  0.9134243E-02 ,  0.9075823E-02 ,  0.9018121E-02 ,  0.8961143E-02 ,&
        &  0.8904889E-02 ,  0.8849339E-02 ,  0.8794466E-02 ,  0.8740273E-02 ,  0.8686719E-02 ,&
        &  0.8633844E-02 ,  0.8581582E-02 ,  0.8529960E-02 ,  0.8478944E-02 ,  0.8428535E-02 ,&
        &  0.8378728E-02 ,  0.8329508E-02 ,  0.8280840E-02 ,  0.8232752E-02 ,  0.8185219E-02 ,&
        &  0.8138204E-02 ,  0.8091775E-02 ,  0.8045806E-02 ,  0.8000427E-02 ,  0.7955516E-02 ,&
        &  0.7911071E-02 ,  0.7867143E-02 ,  0.7823722E-02 ,  0.7780785E-02 ,  0.7738311E-02 ,&
        &  0.7696189E-02 ,  0.7654707E-02 ,  0.7613541E-02 ,  0.7572883E-02 ,  0.7532592E-02 ,&
        &  0.7492736E-02 ,  0.7453269E-02 ,  0.7414199E-02 ,  0.7375633E-02 ,  0.7337505E-02 ,&
        &  0.7299658E-02 ,  0.7262144E-02 ,  0.7225117E-02 ,  0.7188503E-02 ,  0.7152102E-02 /)
      absice4(:,  7) = (/ &
        &  0.5988429E-01 ,  0.6213512E-01 ,  0.6341273E-01 ,  0.6381485E-01 ,  0.6361549E-01 ,&
        &  0.6304067E-01 ,  0.6220275E-01 ,  0.6115729E-01 ,  0.5994790E-01 ,  0.5861807E-01 ,&
        &  0.5721333E-01 ,  0.5577239E-01 ,  0.5432652E-01 ,  0.5289872E-01 ,  0.5150447E-01 ,&
        &  0.5015296E-01 ,  0.4884947E-01 ,  0.4759751E-01 ,  0.4639460E-01 ,  0.4524262E-01 ,&
        &  0.4413919E-01 ,  0.4308327E-01 ,  0.4207106E-01 ,  0.4110039E-01 ,  0.4017112E-01 ,&
        &  0.3927908E-01 ,  0.3842300E-01 ,  0.3760041E-01 ,  0.3681007E-01 ,  0.3605022E-01 ,&
        &  0.3531944E-01 ,  0.3461519E-01 ,  0.3393644E-01 ,  0.3328240E-01 ,  0.3265181E-01 ,&
        &  0.3204254E-01 ,  0.3145472E-01 ,  0.3088620E-01 ,  0.3033656E-01 ,  0.2980512E-01 ,&
        &  0.2929087E-01 ,  0.2879296E-01 ,  0.2831038E-01 ,  0.2784294E-01 ,  0.2738995E-01 ,&
        &  0.2694988E-01 ,  0.2652316E-01 ,  0.2610897E-01 ,  0.2570668E-01 ,  0.2531582E-01 ,&
        &  0.2493591E-01 ,  0.2456660E-01 ,  0.2420738E-01 ,  0.2385787E-01 ,  0.2351768E-01 ,&
        &  0.2318675E-01 ,  0.2286419E-01 ,  0.2254977E-01 ,  0.2224375E-01 ,  0.2194521E-01 ,&
        &  0.2165415E-01 ,  0.2137040E-01 ,  0.2109344E-01 ,  0.2082305E-01 ,  0.2055916E-01 ,&
        &  0.2030159E-01 ,  0.2005026E-01 ,  0.1980428E-01 ,  0.1956418E-01 ,  0.1932957E-01 ,&
        &  0.1910021E-01 ,  0.1887583E-01 ,  0.1865637E-01 ,  0.1844193E-01 ,  0.1823181E-01 ,&
        &  0.1802645E-01 ,  0.1782523E-01 ,  0.1762840E-01 ,  0.1743545E-01 ,  0.1724667E-01 ,&
        &  0.1706167E-01 ,  0.1688035E-01 ,  0.1670260E-01 ,  0.1652849E-01 ,  0.1635783E-01 ,&
        &  0.1619058E-01 ,  0.1602634E-01 ,  0.1586541E-01 ,  0.1570742E-01 ,  0.1555245E-01 ,&
        &  0.1540040E-01 ,  0.1525120E-01 ,  0.1510470E-01 ,  0.1496085E-01 ,  0.1481957E-01 ,&
        &  0.1468085E-01 ,  0.1454461E-01 ,  0.1441077E-01 ,  0.1427930E-01 ,  0.1415001E-01 ,&
        &  0.1402306E-01 ,  0.1389828E-01 ,  0.1377560E-01 ,  0.1365505E-01 ,  0.1353639E-01 ,&
        &  0.1341974E-01 ,  0.1330503E-01 ,  0.1319216E-01 ,  0.1308115E-01 ,  0.1297194E-01 ,&
        &  0.1286447E-01 ,  0.1275865E-01 ,  0.1265455E-01 ,  0.1255203E-01 ,  0.1245118E-01 ,&
        &  0.1235180E-01 ,  0.1225400E-01 ,  0.1215766E-01 ,  0.1206279E-01 ,  0.1196933E-01 ,&
        &  0.1187727E-01 ,  0.1178658E-01 ,  0.1169719E-01 ,  0.1160912E-01 ,  0.1152233E-01 ,&
        &  0.1143681E-01 ,  0.1135248E-01 ,  0.1126936E-01 ,  0.1118741E-01 ,  0.1110664E-01 ,&
        &  0.1102696E-01 ,  0.1094841E-01 ,  0.1087093E-01 ,  0.1079450E-01 ,  0.1071911E-01 ,&
        &  0.1064475E-01 ,  0.1057136E-01 ,  0.1049896E-01 ,  0.1042754E-01 ,  0.1035707E-01 ,&
        &  0.1028748E-01 ,  0.1021882E-01 ,  0.1015106E-01 ,  0.1008416E-01 ,  0.1001815E-01 ,&
        &  0.9952932E-02 ,  0.9888554E-02 ,  0.9824975E-02 ,  0.9762178E-02 ,  0.9700185E-02 ,&
        &  0.9638965E-02 ,  0.9578467E-02 ,  0.9518731E-02 ,  0.9459637E-02 ,  0.9401373E-02 ,&
        &  0.9343728E-02 ,  0.9286805E-02 ,  0.9230558E-02 ,  0.9174970E-02 ,  0.9120035E-02 ,&
        &  0.9065741E-02 ,  0.9012085E-02 ,  0.8959038E-02 ,  0.8906613E-02 ,  0.8854765E-02 ,&
        &  0.8803529E-02 ,  0.8752861E-02 ,  0.8702757E-02 ,  0.8653214E-02 ,  0.8604226E-02 ,&
        &  0.8555784E-02 ,  0.8507868E-02 ,  0.8460473E-02 ,  0.8413609E-02 ,  0.8367250E-02 ,&
        &  0.8321369E-02 ,  0.8276025E-02 ,  0.8231102E-02 ,  0.8186728E-02 ,  0.8142780E-02 ,&
        &  0.8099272E-02 ,  0.8056236E-02 ,  0.8013674E-02 ,  0.7971568E-02 ,  0.7929876E-02 ,&
        &  0.7888515E-02 ,  0.7847751E-02 ,  0.7807292E-02 ,  0.7767302E-02 ,  0.7727637E-02 ,&
        &  0.7688394E-02 ,  0.7649507E-02 ,  0.7610995E-02 ,  0.7572962E-02 ,  0.7535344E-02 ,&
        &  0.7497969E-02 ,  0.7460913E-02 ,  0.7424326E-02 ,  0.7388133E-02 ,  0.7352122E-02 /)
      absice4(:,  8) = (/ &
        &  0.5316645E-01 ,  0.5669206E-01 ,  0.5899668E-01 ,  0.6020352E-01 ,  0.6067482E-01 ,&
        &  0.6062384E-01 ,  0.6017047E-01 ,  0.5940704E-01 ,  0.5841014E-01 ,  0.5724491E-01 ,&
        &  0.5596949E-01 ,  0.5463091E-01 ,  0.5326691E-01 ,  0.5190575E-01 ,  0.5056658E-01 ,&
        &  0.4926201E-01 ,  0.4799948E-01 ,  0.4678414E-01 ,  0.4561472E-01 ,  0.4449370E-01 ,&
        &  0.4341939E-01 ,  0.4239098E-01 ,  0.4140490E-01 ,  0.4045930E-01 ,  0.3955400E-01 ,&
        &  0.3868498E-01 ,  0.3785090E-01 ,  0.3704946E-01 ,  0.3627947E-01 ,  0.3553910E-01 ,&
        &  0.3482699E-01 ,  0.3414068E-01 ,  0.3347914E-01 ,  0.3284159E-01 ,  0.3222676E-01 ,&
        &  0.3163256E-01 ,  0.3105925E-01 ,  0.3050460E-01 ,  0.2996822E-01 ,  0.2944955E-01 ,&
        &  0.2894747E-01 ,  0.2846124E-01 ,  0.2798989E-01 ,  0.2753317E-01 ,  0.2709043E-01 ,&
        &  0.2666026E-01 ,  0.2624300E-01 ,  0.2583789E-01 ,  0.2544429E-01 ,  0.2506179E-01 ,&
        &  0.2468992E-01 ,  0.2432832E-01 ,  0.2397653E-01 ,  0.2363414E-01 ,  0.2330076E-01 ,&
        &  0.2297641E-01 ,  0.2266014E-01 ,  0.2235183E-01 ,  0.2205164E-01 ,  0.2175870E-01 ,&
        &  0.2147305E-01 ,  0.2119451E-01 ,  0.2092256E-01 ,  0.2065701E-01 ,  0.2039777E-01 ,&
        &  0.2014465E-01 ,  0.1989763E-01 ,  0.1965580E-01 ,  0.1941969E-01 ,  0.1918890E-01 ,&
        &  0.1896325E-01 ,  0.1874246E-01 ,  0.1852644E-01 ,  0.1831534E-01 ,  0.1810842E-01 ,&
        &  0.1790614E-01 ,  0.1770788E-01 ,  0.1751395E-01 ,  0.1732378E-01 ,  0.1713766E-01 ,&
        &  0.1695525E-01 ,  0.1677645E-01 ,  0.1660109E-01 ,  0.1642931E-01 ,  0.1626091E-01 ,&
        &  0.1609584E-01 ,  0.1593369E-01 ,  0.1577479E-01 ,  0.1561876E-01 ,  0.1546569E-01 ,&
        &  0.1531547E-01 ,  0.1516804E-01 ,  0.1502327E-01 ,  0.1488108E-01 ,  0.1474140E-01 ,&
        &  0.1460423E-01 ,  0.1446952E-01 ,  0.1433713E-01 ,  0.1420707E-01 ,  0.1407915E-01 ,&
        &  0.1395353E-01 ,  0.1383003E-01 ,  0.1370858E-01 ,  0.1358924E-01 ,  0.1347176E-01 ,&
        &  0.1335625E-01 ,  0.1324264E-01 ,  0.1313083E-01 ,  0.1302087E-01 ,  0.1291265E-01 ,&
        &  0.1280614E-01 ,  0.1270127E-01 ,  0.1259809E-01 ,  0.1249645E-01 ,  0.1239647E-01 ,&
        &  0.1229793E-01 ,  0.1220095E-01 ,  0.1210540E-01 ,  0.1201129E-01 ,  0.1191858E-01 ,&
        &  0.1182725E-01 ,  0.1173727E-01 ,  0.1164857E-01 ,  0.1156117E-01 ,  0.1147501E-01 ,&
        &  0.1139013E-01 ,  0.1130642E-01 ,  0.1122390E-01 ,  0.1114252E-01 ,  0.1106231E-01 ,&
        &  0.1098318E-01 ,  0.1090517E-01 ,  0.1082820E-01 ,  0.1075229E-01 ,  0.1067739E-01 ,&
        &  0.1060351E-01 ,  0.1053059E-01 ,  0.1045864E-01 ,  0.1038766E-01 ,  0.1031763E-01 ,&
        &  0.1024846E-01 ,  0.1018022E-01 ,  0.1011285E-01 ,  0.1004634E-01 ,  0.9980719E-02 ,&
        &  0.9915877E-02 ,  0.9851866E-02 ,  0.9788631E-02 ,  0.9726182E-02 ,  0.9664530E-02 ,&
        &  0.9603629E-02 ,  0.9543454E-02 ,  0.9484039E-02 ,  0.9425251E-02 ,  0.9367287E-02 ,&
        &  0.9309941E-02 ,  0.9253291E-02 ,  0.9197325E-02 ,  0.9142016E-02 ,  0.9087344E-02 ,&
        &  0.9033309E-02 ,  0.8979913E-02 ,  0.8927111E-02 ,  0.8874926E-02 ,  0.8823317E-02 ,&
        &  0.8772308E-02 ,  0.8721864E-02 ,  0.8671990E-02 ,  0.8622666E-02 ,  0.8573883E-02 ,&
        &  0.8525646E-02 ,  0.8477943E-02 ,  0.8430738E-02 ,  0.8384071E-02 ,  0.8337907E-02 ,&
        &  0.8292209E-02 ,  0.8247057E-02 ,  0.8202315E-02 ,  0.8158118E-02 ,  0.8114339E-02 ,&
        &  0.8070996E-02 ,  0.8028133E-02 ,  0.7985732E-02 ,  0.7943779E-02 ,  0.7902252E-02 ,&
        &  0.7861041E-02 ,  0.7820432E-02 ,  0.7780109E-02 ,  0.7740271E-02 ,  0.7700749E-02 ,&
        &  0.7661647E-02 ,  0.7622900E-02 ,  0.7584518E-02 ,  0.7546620E-02 ,  0.7509130E-02 ,&
        &  0.7471889E-02 ,  0.7434956E-02 ,  0.7398485E-02 ,  0.7362414E-02 ,  0.7326525E-02 /)
      absice4(:,  9) = (/ &
        &  0.6496758E-01 ,  0.7022101E-01 ,  0.7327212E-01 ,  0.7463108E-01 ,  0.7488012E-01 ,&
        &  0.7434294E-01 ,  0.7325599E-01 ,  0.7179786E-01 ,  0.7009087E-01 ,  0.6822387E-01 ,&
        &  0.6626864E-01 ,  0.6428009E-01 ,  0.6230181E-01 ,  0.6036458E-01 ,  0.5848880E-01 ,&
        &  0.5668619E-01 ,  0.5496276E-01 ,  0.5332169E-01 ,  0.5175852E-01 ,  0.5027399E-01 ,&
        &  0.4886366E-01 ,  0.4752460E-01 ,  0.4625057E-01 ,  0.4503773E-01 ,  0.4388452E-01 ,&
        &  0.4278488E-01 ,  0.4173614E-01 ,  0.4073441E-01 ,  0.3977751E-01 ,  0.3886252E-01 ,&
        &  0.3798715E-01 ,  0.3714772E-01 ,  0.3634258E-01 ,  0.3557024E-01 ,  0.3482886E-01 ,&
        &  0.3411543E-01 ,  0.3343000E-01 ,  0.3276956E-01 ,  0.3213346E-01 ,  0.3152065E-01 ,&
        &  0.3092959E-01 ,  0.3035924E-01 ,  0.2980824E-01 ,  0.2927616E-01 ,  0.2876202E-01 ,&
        &  0.2826403E-01 ,  0.2778243E-01 ,  0.2731622E-01 ,  0.2686462E-01 ,  0.2642692E-01 ,&
        &  0.2600251E-01 ,  0.2559095E-01 ,  0.2519150E-01 ,  0.2480375E-01 ,  0.2442708E-01 ,&
        &  0.2406147E-01 ,  0.2370578E-01 ,  0.2335980E-01 ,  0.2302366E-01 ,  0.2269637E-01 ,&
        &  0.2237785E-01 ,  0.2206788E-01 ,  0.2176580E-01 ,  0.2147140E-01 ,  0.2118454E-01 ,&
        &  0.2090497E-01 ,  0.2063261E-01 ,  0.2036643E-01 ,  0.2010699E-01 ,  0.1985381E-01 ,&
        &  0.1960664E-01 ,  0.1936518E-01 ,  0.1912933E-01 ,  0.1889915E-01 ,  0.1867389E-01 ,&
        &  0.1845400E-01 ,  0.1823877E-01 ,  0.1802851E-01 ,  0.1782261E-01 ,  0.1762138E-01 ,&
        &  0.1742439E-01 ,  0.1723153E-01 ,  0.1704264E-01 ,  0.1685783E-01 ,  0.1667684E-01 ,&
        &  0.1649965E-01 ,  0.1632578E-01 ,  0.1615561E-01 ,  0.1598867E-01 ,  0.1582506E-01 ,&
        &  0.1566468E-01 ,  0.1550743E-01 ,  0.1535318E-01 ,  0.1520180E-01 ,  0.1505327E-01 ,&
        &  0.1490753E-01 ,  0.1476451E-01 ,  0.1462411E-01 ,  0.1448631E-01 ,  0.1435088E-01 ,&
        &  0.1421798E-01 ,  0.1408745E-01 ,  0.1395920E-01 ,  0.1383326E-01 ,  0.1370939E-01 ,&
        &  0.1358767E-01 ,  0.1346806E-01 ,  0.1335045E-01 ,  0.1323484E-01 ,  0.1312116E-01 ,&
        &  0.1300934E-01 ,  0.1289934E-01 ,  0.1279116E-01 ,  0.1268469E-01 ,  0.1258001E-01 ,&
        &  0.1247691E-01 ,  0.1237550E-01 ,  0.1227566E-01 ,  0.1217736E-01 ,  0.1208059E-01 ,&
        &  0.1198533E-01 ,  0.1189153E-01 ,  0.1179910E-01 ,  0.1170809E-01 ,  0.1161842E-01 ,&
        &  0.1153011E-01 ,  0.1144309E-01 ,  0.1135733E-01 ,  0.1127282E-01 ,  0.1118955E-01 ,&
        &  0.1110745E-01 ,  0.1102654E-01 ,  0.1094676E-01 ,  0.1086811E-01 ,  0.1079053E-01 ,&
        &  0.1071407E-01 ,  0.1063862E-01 ,  0.1056421E-01 ,  0.1049083E-01 ,  0.1041847E-01 ,&
        &  0.1034702E-01 ,  0.1027656E-01 ,  0.1020704E-01 ,  0.1013842E-01 ,  0.1007074E-01 ,&
        &  0.1000390E-01 ,  0.9937946E-02 ,  0.9872812E-02 ,  0.9808511E-02 ,  0.9745046E-02 ,&
        &  0.9682384E-02 ,  0.9620491E-02 ,  0.9559391E-02 ,  0.9498971E-02 ,  0.9439414E-02 ,&
        &  0.9380497E-02 ,  0.9322335E-02 ,  0.9264885E-02 ,  0.9208121E-02 ,  0.9152032E-02 ,&
        &  0.9096626E-02 ,  0.9041875E-02 ,  0.8987754E-02 ,  0.8934287E-02 ,  0.8881420E-02 ,&
        &  0.8829186E-02 ,  0.8777541E-02 ,  0.8726498E-02 ,  0.8676030E-02 ,  0.8626126E-02 ,&
        &  0.8576797E-02 ,  0.8528024E-02 ,  0.8479780E-02 ,  0.8432093E-02 ,  0.8384922E-02 ,&
        &  0.8338246E-02 ,  0.8292145E-02 ,  0.8246464E-02 ,  0.8201354E-02 ,  0.8156692E-02 ,&
        &  0.8112466E-02 ,  0.8068747E-02 ,  0.8025517E-02 ,  0.7982743E-02 ,  0.7940413E-02 ,&
        &  0.7898415E-02 ,  0.7857045E-02 ,  0.7815976E-02 ,  0.7775393E-02 ,  0.7735157E-02 ,&
        &  0.7695348E-02 ,  0.7655912E-02 ,  0.7616862E-02 ,  0.7578297E-02 ,  0.7540161E-02 ,&
        &  0.7502281E-02 ,  0.7464732E-02 ,  0.7427657E-02 ,  0.7390989E-02 ,  0.7354508E-02 /)
      absice4(:, 10) = (/ &
        &  0.9604552E-01 ,  0.1033564E+00 ,  0.1063436E+00 ,  0.1065534E+00 ,  0.1050589E+00 ,&
        &  0.1025090E+00 ,  0.9936723E-01 ,  0.9592369E-01 ,  0.9234641E-01 ,  0.8873942E-01 ,&
        &  0.8517904E-01 ,  0.8171826E-01 ,  0.7839544E-01 ,  0.7523430E-01 ,  0.7224611E-01 ,&
        &  0.6943324E-01 ,  0.6679222E-01 ,  0.6431767E-01 ,  0.6199514E-01 ,  0.5981901E-01 ,&
        &  0.5777737E-01 ,  0.5586144E-01 ,  0.5405857E-01 ,  0.5236009E-01 ,  0.5076092E-01 ,&
        &  0.4925022E-01 ,  0.4782211E-01 ,  0.4646969E-01 ,  0.4518828E-01 ,  0.4397243E-01 ,&
        &  0.4281782E-01 ,  0.4171859E-01 ,  0.4067142E-01 ,  0.3967357E-01 ,  0.3872176E-01 ,&
        &  0.3781148E-01 ,  0.3694203E-01 ,  0.3610912E-01 ,  0.3531129E-01 ,  0.3454673E-01 ,&
        &  0.3381307E-01 ,  0.3310869E-01 ,  0.3243148E-01 ,  0.3178048E-01 ,  0.3115433E-01 ,&
        &  0.3055052E-01 ,  0.2996903E-01 ,  0.2940849E-01 ,  0.2886758E-01 ,  0.2834543E-01 ,&
        &  0.2784102E-01 ,  0.2735367E-01 ,  0.2688234E-01 ,  0.2642638E-01 ,  0.2598499E-01 ,&
        &  0.2555793E-01 ,  0.2514378E-01 ,  0.2474218E-01 ,  0.2435319E-01 ,  0.2397549E-01 ,&
        &  0.2360903E-01 ,  0.2325335E-01 ,  0.2290771E-01 ,  0.2257171E-01 ,  0.2224517E-01 ,&
        &  0.2192774E-01 ,  0.2161921E-01 ,  0.2131843E-01 ,  0.2102594E-01 ,  0.2074116E-01 ,&
        &  0.2046377E-01 ,  0.2019334E-01 ,  0.1992975E-01 ,  0.1967304E-01 ,  0.1942231E-01 ,&
        &  0.1917804E-01 ,  0.1893941E-01 ,  0.1870671E-01 ,  0.1847926E-01 ,  0.1825737E-01 ,&
        &  0.1804052E-01 ,  0.1782860E-01 ,  0.1762136E-01 ,  0.1741894E-01 ,  0.1722101E-01 ,&
        &  0.1702752E-01 ,  0.1683798E-01 ,  0.1665270E-01 ,  0.1647123E-01 ,  0.1629365E-01 ,&
        &  0.1611978E-01 ,  0.1594955E-01 ,  0.1578278E-01 ,  0.1561934E-01 ,  0.1545916E-01 ,&
        &  0.1530219E-01 ,  0.1514835E-01 ,  0.1499749E-01 ,  0.1484962E-01 ,  0.1470442E-01 ,&
        &  0.1456212E-01 ,  0.1442250E-01 ,  0.1428547E-01 ,  0.1415105E-01 ,  0.1401897E-01 ,&
        &  0.1388932E-01 ,  0.1376202E-01 ,  0.1363698E-01 ,  0.1351418E-01 ,  0.1339355E-01 ,&
        &  0.1327501E-01 ,  0.1315848E-01 ,  0.1304399E-01 ,  0.1293140E-01 ,  0.1282079E-01 ,&
        &  0.1271195E-01 ,  0.1260498E-01 ,  0.1249973E-01 ,  0.1239620E-01 ,  0.1229437E-01 ,&
        &  0.1219417E-01 ,  0.1209559E-01 ,  0.1199852E-01 ,  0.1190302E-01 ,  0.1180898E-01 ,&
        &  0.1171643E-01 ,  0.1162529E-01 ,  0.1153554E-01 ,  0.1144713E-01 ,  0.1136010E-01 ,&
        &  0.1127432E-01 ,  0.1118984E-01 ,  0.1110659E-01 ,  0.1102458E-01 ,  0.1094372E-01 ,&
        &  0.1086406E-01 ,  0.1078551E-01 ,  0.1070808E-01 ,  0.1063177E-01 ,  0.1055654E-01 ,&
        &  0.1048232E-01 ,  0.1040915E-01 ,  0.1033699E-01 ,  0.1026579E-01 ,  0.1019562E-01 ,&
        &  0.1012633E-01 ,  0.1005800E-01 ,  0.9990551E-02 ,  0.9923989E-02 ,  0.9858332E-02 ,&
        &  0.9793528E-02 ,  0.9729533E-02 ,  0.9666396E-02 ,  0.9603987E-02 ,  0.9542488E-02 ,&
        &  0.9481681E-02 ,  0.9421664E-02 ,  0.9362411E-02 ,  0.9303899E-02 ,  0.9246089E-02 ,&
        &  0.9188998E-02 ,  0.9132609E-02 ,  0.9076898E-02 ,  0.9021870E-02 ,  0.8967470E-02 ,&
        &  0.8913750E-02 ,  0.8860655E-02 ,  0.8808181E-02 ,  0.8756323E-02 ,  0.8705065E-02 ,&
        &  0.8654407E-02 ,  0.8604337E-02 ,  0.8554813E-02 ,  0.8505888E-02 ,  0.8457500E-02 ,&
        &  0.8409650E-02 ,  0.8362376E-02 ,  0.8315563E-02 ,  0.8269335E-02 ,  0.8223582E-02 ,&
        &  0.8178296E-02 ,  0.8133536E-02 ,  0.8089283E-02 ,  0.8045509E-02 ,  0.8002207E-02 ,&
        &  0.7959244E-02 ,  0.7916938E-02 ,  0.7874951E-02 ,  0.7833467E-02 ,  0.7792348E-02 ,&
        &  0.7751674E-02 ,  0.7711396E-02 ,  0.7671515E-02 ,  0.7632134E-02 ,  0.7593200E-02 ,&
        &  0.7554544E-02 ,  0.7516227E-02 ,  0.7478401E-02 ,  0.7440997E-02 ,  0.7403795E-02 /)
      absice4(:, 11) = (/ &
        &  0.1120205E+00 ,  0.1195659E+00 ,  0.1215598E+00 ,  0.1203989E+00 ,  0.1174667E+00 ,&
        &  0.1135672E+00 ,  0.1092127E+00 ,  0.1046952E+00 ,  0.1001720E+00 ,  0.9573516E-01 ,&
        &  0.9144754E-01 ,  0.8734994E-01 ,  0.8346878E-01 ,  0.7981736E-01 ,  0.7639769E-01 ,&
        &  0.7320409E-01 ,  0.7022630E-01 ,  0.6745320E-01 ,  0.6486455E-01 ,  0.6245110E-01 ,&
        &  0.6019717E-01 ,  0.5809085E-01 ,  0.5611662E-01 ,  0.5426347E-01 ,  0.5252461E-01 ,&
        &  0.5088725E-01 ,  0.4934417E-01 ,  0.4788693E-01 ,  0.4651003E-01 ,  0.4520687E-01 ,&
        &  0.4397245E-01 ,  0.4279992E-01 ,  0.4168537E-01 ,  0.4062553E-01 ,  0.3961669E-01 ,&
        &  0.3865369E-01 ,  0.3773556E-01 ,  0.3685754E-01 ,  0.3601797E-01 ,  0.3521462E-01 ,&
        &  0.3444500E-01 ,  0.3370715E-01 ,  0.3299877E-01 ,  0.3231877E-01 ,  0.3166556E-01 ,&
        &  0.3103643E-01 ,  0.3043129E-01 ,  0.2984862E-01 ,  0.2928701E-01 ,  0.2874544E-01 ,&
        &  0.2822283E-01 ,  0.2771834E-01 ,  0.2723097E-01 ,  0.2675988E-01 ,  0.2630426E-01 ,&
        &  0.2586380E-01 ,  0.2543702E-01 ,  0.2502349E-01 ,  0.2462325E-01 ,  0.2423494E-01 ,&
        &  0.2385843E-01 ,  0.2349324E-01 ,  0.2313862E-01 ,  0.2279409E-01 ,  0.2245949E-01 ,&
        &  0.2213439E-01 ,  0.2181861E-01 ,  0.2151091E-01 ,  0.2121187E-01 ,  0.2092086E-01 ,&
        &  0.2063755E-01 ,  0.2036147E-01 ,  0.2009252E-01 ,  0.1983069E-01 ,  0.1957508E-01 ,&
        &  0.1932616E-01 ,  0.1908309E-01 ,  0.1884616E-01 ,  0.1861468E-01 ,  0.1838890E-01 ,&
        &  0.1816836E-01 ,  0.1795290E-01 ,  0.1774228E-01 ,  0.1753662E-01 ,  0.1733558E-01 ,&
        &  0.1713913E-01 ,  0.1694672E-01 ,  0.1675871E-01 ,  0.1657463E-01 ,  0.1639452E-01 ,&
        &  0.1621823E-01 ,  0.1604568E-01 ,  0.1587667E-01 ,  0.1571108E-01 ,  0.1554883E-01 ,&
        &  0.1538988E-01 ,  0.1523411E-01 ,  0.1508142E-01 ,  0.1493174E-01 ,  0.1478484E-01 ,&
        &  0.1464089E-01 ,  0.1449966E-01 ,  0.1436109E-01 ,  0.1422518E-01 ,  0.1409166E-01 ,&
        &  0.1396063E-01 ,  0.1383200E-01 ,  0.1370565E-01 ,  0.1358160E-01 ,  0.1345975E-01 ,&
        &  0.1334003E-01 ,  0.1322237E-01 ,  0.1310677E-01 ,  0.1299311E-01 ,  0.1288147E-01 ,&
        &  0.1277162E-01 ,  0.1266367E-01 ,  0.1255748E-01 ,  0.1245304E-01 ,  0.1235030E-01 ,&
        &  0.1224925E-01 ,  0.1214983E-01 ,  0.1205194E-01 ,  0.1195563E-01 ,  0.1186083E-01 ,&
        &  0.1176754E-01 ,  0.1167565E-01 ,  0.1158518E-01 ,  0.1149609E-01 ,  0.1140837E-01 ,&
        &  0.1132194E-01 ,  0.1123682E-01 ,  0.1115295E-01 ,  0.1107031E-01 ,  0.1098886E-01 ,&
        &  0.1090863E-01 ,  0.1082951E-01 ,  0.1075154E-01 ,  0.1067468E-01 ,  0.1059892E-01 ,&
        &  0.1052418E-01 ,  0.1045051E-01 ,  0.1037786E-01 ,  0.1030619E-01 ,  0.1023554E-01 ,&
        &  0.1016580E-01 ,  0.1009701E-01 ,  0.1002912E-01 ,  0.9962142E-02 ,  0.9896057E-02 ,&
        &  0.9830839E-02 ,  0.9766445E-02 ,  0.9702916E-02 ,  0.9640107E-02 ,  0.9578238E-02 ,&
        &  0.9517062E-02 ,  0.9456682E-02 ,  0.9397073E-02 ,  0.9338205E-02 ,  0.9280046E-02 ,&
        &  0.9222627E-02 ,  0.9165908E-02 ,  0.9109871E-02 ,  0.9054518E-02 ,  0.8999816E-02 ,&
        &  0.8945790E-02 ,  0.8892382E-02 ,  0.8839618E-02 ,  0.8787465E-02 ,  0.8735923E-02 ,&
        &  0.8684986E-02 ,  0.8634641E-02 ,  0.8584852E-02 ,  0.8535657E-02 ,  0.8487010E-02 ,&
        &  0.8438895E-02 ,  0.8391378E-02 ,  0.8344316E-02 ,  0.8297840E-02 ,  0.8251843E-02 ,&
        &  0.8206324E-02 ,  0.8161333E-02 ,  0.8116846E-02 ,  0.8072853E-02 ,  0.8029319E-02 ,&
        &  0.7986145E-02 ,  0.7943621E-02 ,  0.7901418E-02 ,  0.7859729E-02 ,  0.7818406E-02 ,&
        &  0.7777532E-02 ,  0.7737046E-02 ,  0.7696967E-02 ,  0.7657392E-02 ,  0.7618273E-02 ,&
        &  0.7579426E-02 ,  0.7540921E-02 ,  0.7502915E-02 ,  0.7465327E-02 ,  0.7427948E-02 /)
      absice4(:, 12) = (/ &
        &  0.4118567E-01 ,  0.4606692E-01 ,  0.4892185E-01 ,  0.5047753E-01 ,  0.5110885E-01 ,&
        &  0.5113748E-01 ,  0.5078274E-01 ,  0.5017032E-01 ,  0.4937394E-01 ,  0.4844657E-01 ,&
        &  0.4743351E-01 ,  0.4637205E-01 ,  0.4529282E-01 ,  0.4421832E-01 ,  0.4316403E-01 ,&
        &  0.4213925E-01 ,  0.4114978E-01 ,  0.4019905E-01 ,  0.3928546E-01 ,  0.3841056E-01 ,&
        &  0.3757252E-01 ,  0.3677037E-01 ,  0.3600108E-01 ,  0.3526276E-01 ,  0.3455525E-01 ,&
        &  0.3387534E-01 ,  0.3322180E-01 ,  0.3259278E-01 ,  0.3198735E-01 ,  0.3140405E-01 ,&
        &  0.3084194E-01 ,  0.3029896E-01 ,  0.2977444E-01 ,  0.2926775E-01 ,  0.2877798E-01 ,&
        &  0.2830354E-01 ,  0.2784461E-01 ,  0.2739959E-01 ,  0.2696813E-01 ,  0.2654990E-01 ,&
        &  0.2614401E-01 ,  0.2575001E-01 ,  0.2536709E-01 ,  0.2499515E-01 ,  0.2463372E-01 ,&
        &  0.2428169E-01 ,  0.2393935E-01 ,  0.2360624E-01 ,  0.2328183E-01 ,  0.2296579E-01 ,&
        &  0.2265779E-01 ,  0.2235767E-01 ,  0.2206498E-01 ,  0.2177947E-01 ,  0.2150086E-01 ,&
        &  0.2122923E-01 ,  0.2096377E-01 ,  0.2070442E-01 ,  0.2045138E-01 ,  0.2020396E-01 ,&
        &  0.1996214E-01 ,  0.1972592E-01 ,  0.1949481E-01 ,  0.1926867E-01 ,  0.1904748E-01 ,&
        &  0.1883111E-01 ,  0.1861956E-01 ,  0.1841206E-01 ,  0.1820911E-01 ,  0.1801038E-01 ,&
        &  0.1781570E-01 ,  0.1762488E-01 ,  0.1743788E-01 ,  0.1725481E-01 ,  0.1707507E-01 ,&
        &  0.1689909E-01 ,  0.1672629E-01 ,  0.1655700E-01 ,  0.1639072E-01 ,  0.1622774E-01 ,&
        &  0.1606775E-01 ,  0.1591071E-01 ,  0.1575642E-01 ,  0.1560509E-01 ,  0.1545651E-01 ,&
        &  0.1531065E-01 ,  0.1516716E-01 ,  0.1502639E-01 ,  0.1488792E-01 ,  0.1475191E-01 ,&
        &  0.1461826E-01 ,  0.1448690E-01 ,  0.1435776E-01 ,  0.1423075E-01 ,  0.1410582E-01 ,&
        &  0.1398299E-01 ,  0.1386217E-01 ,  0.1374334E-01 ,  0.1362645E-01 ,  0.1351132E-01 ,&
        &  0.1339814E-01 ,  0.1328674E-01 ,  0.1317708E-01 ,  0.1306918E-01 ,  0.1296286E-01 ,&
        &  0.1285818E-01 ,  0.1275513E-01 ,  0.1265361E-01 ,  0.1255364E-01 ,  0.1245518E-01 ,&
        &  0.1235816E-01 ,  0.1226254E-01 ,  0.1216836E-01 ,  0.1207550E-01 ,  0.1198406E-01 ,&
        &  0.1189386E-01 ,  0.1180500E-01 ,  0.1171736E-01 ,  0.1163096E-01 ,  0.1154577E-01 ,&
        &  0.1146177E-01 ,  0.1137894E-01 ,  0.1129722E-01 ,  0.1121663E-01 ,  0.1113711E-01 ,&
        &  0.1105870E-01 ,  0.1098132E-01 ,  0.1090496E-01 ,  0.1082959E-01 ,  0.1075527E-01 ,&
        &  0.1068187E-01 ,  0.1060945E-01 ,  0.1053796E-01 ,  0.1046739E-01 ,  0.1039769E-01 ,&
        &  0.1032891E-01 ,  0.1026096E-01 ,  0.1019388E-01 ,  0.1012765E-01 ,  0.1006225E-01 ,&
        &  0.9997615E-02 ,  0.9933805E-02 ,  0.9870761E-02 ,  0.9808474E-02 ,  0.9746985E-02 ,&
        &  0.9686189E-02 ,  0.9626136E-02 ,  0.9566763E-02 ,  0.9508105E-02 ,  0.9450142E-02 ,&
        &  0.9392861E-02 ,  0.9336227E-02 ,  0.9280265E-02 ,  0.9224868E-02 ,  0.9170210E-02 ,&
        &  0.9116109E-02 ,  0.9062633E-02 ,  0.9009778E-02 ,  0.8957504E-02 ,  0.8905799E-02 ,&
        &  0.8854680E-02 ,  0.8804129E-02 ,  0.8754121E-02 ,  0.8704669E-02 ,  0.8655739E-02 ,&
        &  0.8607359E-02 ,  0.8559475E-02 ,  0.8512112E-02 ,  0.8465252E-02 ,  0.8418886E-02 ,&
        &  0.8373009E-02 ,  0.8327622E-02 ,  0.8282688E-02 ,  0.8238247E-02 ,  0.8194258E-02 ,&
        &  0.8150694E-02 ,  0.8107632E-02 ,  0.8064942E-02 ,  0.8022754E-02 ,  0.7980954E-02 ,&
        &  0.7939544E-02 ,  0.7898574E-02 ,  0.7858037E-02 ,  0.7817902E-02 ,  0.7778158E-02 ,&
        &  0.7738703E-02 ,  0.7699814E-02 ,  0.7661183E-02 ,  0.7622996E-02 ,  0.7585100E-02 ,&
        &  0.7547590E-02 ,  0.7510403E-02 ,  0.7473565E-02 ,  0.7437170E-02 ,  0.7401152E-02 ,&
        &  0.7365360E-02 ,  0.7329856E-02 ,  0.7294781E-02 ,  0.7260076E-02 ,  0.7225542E-02 /)
      absice4(:, 13) = (/ &
        &  0.7048990E-01 ,  0.7855639E-01 ,  0.8212429E-01 ,  0.8289183E-01 ,  0.8209244E-01 ,&
        &  0.8049434E-01 ,  0.7848005E-01 ,  0.7624222E-01 ,  0.7388593E-01 ,  0.7147966E-01 ,&
        &  0.6907652E-01 ,  0.6671560E-01 ,  0.6442613E-01 ,  0.6222782E-01 ,  0.6013142E-01 ,&
        &  0.5814109E-01 ,  0.5625683E-01 ,  0.5447724E-01 ,  0.5279362E-01 ,  0.5120403E-01 ,&
        &  0.4970152E-01 ,  0.4828129E-01 ,  0.4693533E-01 ,  0.4565854E-01 ,  0.4444838E-01 ,&
        &  0.4329771E-01 ,  0.4220313E-01 ,  0.4116017E-01 ,  0.4016606E-01 ,  0.3921735E-01 ,&
        &  0.3831147E-01 ,  0.3744432E-01 ,  0.3661390E-01 ,  0.3581859E-01 ,  0.3505625E-01 ,&
        &  0.3432364E-01 ,  0.3362067E-01 ,  0.3294416E-01 ,  0.3229332E-01 ,  0.3166698E-01 ,&
        &  0.3106350E-01 ,  0.3048175E-01 ,  0.2992023E-01 ,  0.2937843E-01 ,  0.2885540E-01 ,&
        &  0.2834916E-01 ,  0.2785997E-01 ,  0.2738677E-01 ,  0.2692864E-01 ,  0.2648498E-01 ,&
        &  0.2605503E-01 ,  0.2563833E-01 ,  0.2523417E-01 ,  0.2484199E-01 ,  0.2446127E-01 ,&
        &  0.2409190E-01 ,  0.2373272E-01 ,  0.2338351E-01 ,  0.2304438E-01 ,  0.2271429E-01 ,&
        &  0.2239317E-01 ,  0.2208081E-01 ,  0.2177652E-01 ,  0.2148005E-01 ,  0.2119127E-01 ,&
        &  0.2090993E-01 ,  0.2063591E-01 ,  0.2036820E-01 ,  0.2010734E-01 ,  0.1985286E-01 ,&
        &  0.1960447E-01 ,  0.1936188E-01 ,  0.1912498E-01 ,  0.1889383E-01 ,  0.1866765E-01 ,&
        &  0.1844694E-01 ,  0.1823093E-01 ,  0.1801995E-01 ,  0.1781340E-01 ,  0.1761155E-01 ,&
        &  0.1741400E-01 ,  0.1722062E-01 ,  0.1703127E-01 ,  0.1684600E-01 ,  0.1666461E-01 ,&
        &  0.1648704E-01 ,  0.1631284E-01 ,  0.1614235E-01 ,  0.1597511E-01 ,  0.1581126E-01 ,&
        &  0.1565064E-01 ,  0.1549318E-01 ,  0.1533871E-01 ,  0.1518718E-01 ,  0.1503847E-01 ,&
        &  0.1489258E-01 ,  0.1474944E-01 ,  0.1460891E-01 ,  0.1447101E-01 ,  0.1433550E-01 ,&
        &  0.1420252E-01 ,  0.1407192E-01 ,  0.1394362E-01 ,  0.1381762E-01 ,  0.1369371E-01 ,&
        &  0.1357198E-01 ,  0.1345234E-01 ,  0.1333470E-01 ,  0.1321908E-01 ,  0.1310538E-01 ,&
        &  0.1299357E-01 ,  0.1288357E-01 ,  0.1277541E-01 ,  0.1266894E-01 ,  0.1256429E-01 ,&
        &  0.1246121E-01 ,  0.1235983E-01 ,  0.1226000E-01 ,  0.1216174E-01 ,  0.1206502E-01 ,&
        &  0.1196978E-01 ,  0.1187600E-01 ,  0.1178362E-01 ,  0.1169265E-01 ,  0.1160303E-01 ,&
        &  0.1151476E-01 ,  0.1142779E-01 ,  0.1134208E-01 ,  0.1125761E-01 ,  0.1117440E-01 ,&
        &  0.1109234E-01 ,  0.1101148E-01 ,  0.1093175E-01 ,  0.1085316E-01 ,  0.1077564E-01 ,&
        &  0.1069922E-01 ,  0.1062382E-01 ,  0.1054947E-01 ,  0.1047615E-01 ,  0.1040384E-01 ,&
        &  0.1033245E-01 ,  0.1026204E-01 ,  0.1019257E-01 ,  0.1012400E-01 ,  0.1005638E-01 ,&
        &  0.9989594E-02 ,  0.9923684E-02 ,  0.9858606E-02 ,  0.9794370E-02 ,  0.9730957E-02 ,&
        &  0.9668347E-02 ,  0.9606508E-02 ,  0.9545457E-02 ,  0.9485086E-02 ,  0.9425581E-02 ,&
        &  0.9366722E-02 ,  0.9308609E-02 ,  0.9251215E-02 ,  0.9194500E-02 ,  0.9138465E-02 ,&
        &  0.9083097E-02 ,  0.9028401E-02 ,  0.8974335E-02 ,  0.8920911E-02 ,  0.8868096E-02 ,&
        &  0.8815915E-02 ,  0.8764313E-02 ,  0.8713313E-02 ,  0.8662886E-02 ,  0.8613040E-02 ,&
        &  0.8563751E-02 ,  0.8515026E-02 ,  0.8466822E-02 ,  0.8419185E-02 ,  0.8372060E-02 ,&
        &  0.8325432E-02 ,  0.8279367E-02 ,  0.8233733E-02 ,  0.8188660E-02 ,  0.8144042E-02 ,&
        &  0.8099861E-02 ,  0.8056186E-02 ,  0.8012990E-02 ,  0.7970259E-02 ,  0.7927964E-02 ,&
        &  0.7886007E-02 ,  0.7844678E-02 ,  0.7803643E-02 ,  0.7763108E-02 ,  0.7722905E-02 ,&
        &  0.7683137E-02 ,  0.7643731E-02 ,  0.7604721E-02 ,  0.7566193E-02 ,  0.7528087E-02 ,&
        &  0.7490253E-02 ,  0.7452734E-02 ,  0.7415696E-02 ,  0.7379058E-02 ,  0.7342621E-02 /)
      absice4(:, 14) = (/ &
        &  0.7435944E-01 ,  0.8300151E-01 ,  0.8688792E-01 ,  0.8781894E-01 ,  0.8698697E-01 ,&
        &  0.8522374E-01 ,  0.8298247E-01 ,  0.8049620E-01 ,  0.7788935E-01 ,  0.7523922E-01 ,&
        &  0.7260300E-01 ,  0.7002171E-01 ,  0.6752570E-01 ,  0.6513479E-01 ,  0.6285974E-01 ,&
        &  0.6070410E-01 ,  0.5866707E-01 ,  0.5674654E-01 ,  0.5493267E-01 ,  0.5322292E-01 ,&
        &  0.5160941E-01 ,  0.5008643E-01 ,  0.4864534E-01 ,  0.4728018E-01 ,  0.4598802E-01 ,&
        &  0.4476102E-01 ,  0.4359530E-01 ,  0.4248599E-01 ,  0.4142991E-01 ,  0.4042322E-01 ,&
        &  0.3946301E-01 ,  0.3854491E-01 ,  0.3766656E-01 ,  0.3682620E-01 ,  0.3602141E-01 ,&
        &  0.3524880E-01 ,  0.3450807E-01 ,  0.3379583E-01 ,  0.3311123E-01 ,  0.3245289E-01 ,&
        &  0.3181909E-01 ,  0.3120854E-01 ,  0.3061969E-01 ,  0.3005193E-01 ,  0.2950415E-01 ,&
        &  0.2897440E-01 ,  0.2846277E-01 ,  0.2796819E-01 ,  0.2748969E-01 ,  0.2702655E-01 ,&
        &  0.2657799E-01 ,  0.2614354E-01 ,  0.2572235E-01 ,  0.2531393E-01 ,  0.2491762E-01 ,&
        &  0.2453336E-01 ,  0.2415990E-01 ,  0.2379698E-01 ,  0.2344472E-01 ,  0.2310201E-01 ,&
        &  0.2276881E-01 ,  0.2244484E-01 ,  0.2212939E-01 ,  0.2182219E-01 ,  0.2152313E-01 ,&
        &  0.2123188E-01 ,  0.2094833E-01 ,  0.2067145E-01 ,  0.2040177E-01 ,  0.2013878E-01 ,&
        &  0.1988223E-01 ,  0.1963175E-01 ,  0.1938725E-01 ,  0.1914880E-01 ,  0.1891558E-01 ,&
        &  0.1868806E-01 ,  0.1846547E-01 ,  0.1824818E-01 ,  0.1803550E-01 ,  0.1782776E-01 ,&
        &  0.1762451E-01 ,  0.1742563E-01 ,  0.1723094E-01 ,  0.1704056E-01 ,  0.1685419E-01 ,&
        &  0.1667182E-01 ,  0.1649297E-01 ,  0.1631799E-01 ,  0.1614641E-01 ,  0.1597834E-01 ,&
        &  0.1581366E-01 ,  0.1565224E-01 ,  0.1549397E-01 ,  0.1533873E-01 ,  0.1518645E-01 ,&
        &  0.1503708E-01 ,  0.1489058E-01 ,  0.1474681E-01 ,  0.1460574E-01 ,  0.1446714E-01 ,&
        &  0.1433120E-01 ,  0.1419772E-01 ,  0.1406660E-01 ,  0.1393791E-01 ,  0.1381135E-01 ,&
        &  0.1368705E-01 ,  0.1356492E-01 ,  0.1344487E-01 ,  0.1332691E-01 ,  0.1321094E-01 ,&
        &  0.1309691E-01 ,  0.1298474E-01 ,  0.1287448E-01 ,  0.1276598E-01 ,  0.1265935E-01 ,&
        &  0.1255434E-01 ,  0.1245109E-01 ,  0.1234944E-01 ,  0.1224941E-01 ,  0.1215095E-01 ,&
        &  0.1205403E-01 ,  0.1195862E-01 ,  0.1186464E-01 ,  0.1177212E-01 ,  0.1168099E-01 ,&
        &  0.1159125E-01 ,  0.1150284E-01 ,  0.1141573E-01 ,  0.1132990E-01 ,  0.1124536E-01 ,&
        &  0.1116201E-01 ,  0.1107988E-01 ,  0.1099892E-01 ,  0.1091912E-01 ,  0.1084043E-01 ,&
        &  0.1076287E-01 ,  0.1068636E-01 ,  0.1061092E-01 ,  0.1053653E-01 ,  0.1046318E-01 ,&
        &  0.1039078E-01 ,  0.1031938E-01 ,  0.1024895E-01 ,  0.1017943E-01 ,  0.1011089E-01 ,&
        &  0.1004321E-01 ,  0.9976418E-02 ,  0.9910482E-02 ,  0.9845391E-02 ,  0.9781163E-02 ,&
        &  0.9717750E-02 ,  0.9655127E-02 ,  0.9593316E-02 ,  0.9532193E-02 ,  0.9471956E-02 ,&
        &  0.9412373E-02 ,  0.9353555E-02 ,  0.9295478E-02 ,  0.9238095E-02 ,  0.9181391E-02 ,&
        &  0.9125392E-02 ,  0.9070064E-02 ,  0.9015372E-02 ,  0.8961351E-02 ,  0.8907943E-02 ,&
        &  0.8855177E-02 ,  0.8803016E-02 ,  0.8751454E-02 ,  0.8700491E-02 ,  0.8650105E-02 ,&
        &  0.8600302E-02 ,  0.8551070E-02 ,  0.8502361E-02 ,  0.8454226E-02 ,  0.8406628E-02 ,&
        &  0.8359523E-02 ,  0.8312996E-02 ,  0.8266911E-02 ,  0.8221395E-02 ,  0.8176337E-02 ,&
        &  0.8131738E-02 ,  0.8087633E-02 ,  0.8044029E-02 ,  0.8000894E-02 ,  0.7958217E-02 ,&
        &  0.7915865E-02 ,  0.7874154E-02 ,  0.7832748E-02 ,  0.7791847E-02 ,  0.7751289E-02 ,&
        &  0.7711162E-02 ,  0.7671416E-02 ,  0.7632062E-02 ,  0.7593195E-02 ,  0.7554769E-02 ,&
        &  0.7516610E-02 ,  0.7478775E-02 ,  0.7441428E-02 ,  0.7404482E-02 ,  0.7367747E-02 /)
      absice4(:, 15) = (/ &
        &  0.4860886E-01 ,  0.5537690E-01 ,  0.5939421E-01 ,  0.6117365E-01 ,  0.6141138E-01 ,&
        &  0.6083780E-01 ,  0.5983812E-01 ,  0.5859718E-01 ,  0.5721138E-01 ,  0.5574136E-01 ,&
        &  0.5423342E-01 ,  0.5272206E-01 ,  0.5123409E-01 ,  0.4978795E-01 ,  0.4839501E-01 ,&
        &  0.4706126E-01 ,  0.4578897E-01 ,  0.4457895E-01 ,  0.4342648E-01 ,  0.4233154E-01 ,&
        &  0.4129016E-01 ,  0.4029983E-01 ,  0.3935565E-01 ,  0.3845455E-01 ,  0.3759540E-01 ,&
        &  0.3677373E-01 ,  0.3598761E-01 ,  0.3523422E-01 ,  0.3451207E-01 ,  0.3381907E-01 ,&
        &  0.3315369E-01 ,  0.3251331E-01 ,  0.3189681E-01 ,  0.3130325E-01 ,  0.3073138E-01 ,&
        &  0.3017913E-01 ,  0.2964653E-01 ,  0.2913160E-01 ,  0.2863381E-01 ,  0.2815254E-01 ,&
        &  0.2768677E-01 ,  0.2723581E-01 ,  0.2679860E-01 ,  0.2637500E-01 ,  0.2596439E-01 ,&
        &  0.2556534E-01 ,  0.2517821E-01 ,  0.2480227E-01 ,  0.2443695E-01 ,  0.2408183E-01 ,&
        &  0.2373645E-01 ,  0.2340052E-01 ,  0.2307358E-01 ,  0.2275526E-01 ,  0.2244517E-01 ,&
        &  0.2214337E-01 ,  0.2184898E-01 ,  0.2156183E-01 ,  0.2128212E-01 ,  0.2100906E-01 ,&
        &  0.2074261E-01 ,  0.2048270E-01 ,  0.2022879E-01 ,  0.1998072E-01 ,  0.1973845E-01 ,&
        &  0.1950176E-01 ,  0.1927064E-01 ,  0.1904427E-01 ,  0.1882311E-01 ,  0.1860685E-01 ,&
        &  0.1839526E-01 ,  0.1818811E-01 ,  0.1798534E-01 ,  0.1778708E-01 ,  0.1759264E-01 ,&
        &  0.1740245E-01 ,  0.1721595E-01 ,  0.1703339E-01 ,  0.1685430E-01 ,  0.1667893E-01 ,&
        &  0.1650695E-01 ,  0.1633828E-01 ,  0.1617277E-01 ,  0.1601057E-01 ,  0.1585145E-01 ,&
        &  0.1569541E-01 ,  0.1554205E-01 ,  0.1539170E-01 ,  0.1524395E-01 ,  0.1509895E-01 ,&
        &  0.1495658E-01 ,  0.1481679E-01 ,  0.1467943E-01 ,  0.1454446E-01 ,  0.1441181E-01 ,&
        &  0.1428149E-01 ,  0.1415341E-01 ,  0.1402749E-01 ,  0.1390375E-01 ,  0.1378197E-01 ,&
        &  0.1366232E-01 ,  0.1354463E-01 ,  0.1342886E-01 ,  0.1331502E-01 ,  0.1320292E-01 ,&
        &  0.1309264E-01 ,  0.1298412E-01 ,  0.1287728E-01 ,  0.1277217E-01 ,  0.1266867E-01 ,&
        &  0.1256676E-01 ,  0.1246638E-01 ,  0.1236757E-01 ,  0.1227020E-01 ,  0.1217438E-01 ,&
        &  0.1207990E-01 ,  0.1198687E-01 ,  0.1189518E-01 ,  0.1180483E-01 ,  0.1171579E-01 ,&
        &  0.1162804E-01 ,  0.1154157E-01 ,  0.1145626E-01 ,  0.1137219E-01 ,  0.1128930E-01 ,&
        &  0.1120758E-01 ,  0.1112698E-01 ,  0.1104748E-01 ,  0.1096906E-01 ,  0.1089174E-01 ,&
        &  0.1081542E-01 ,  0.1074015E-01 ,  0.1066588E-01 ,  0.1059260E-01 ,  0.1052026E-01 ,&
        &  0.1044890E-01 ,  0.1037843E-01 ,  0.1030889E-01 ,  0.1024025E-01 ,  0.1017250E-01 ,&
        &  0.1010558E-01 ,  0.1003951E-01 ,  0.9974295E-02 ,  0.9909872E-02 ,  0.9846289E-02 ,&
        &  0.9783447E-02 ,  0.9721389E-02 ,  0.9660077E-02 ,  0.9599508E-02 ,  0.9539680E-02 ,&
        &  0.9480577E-02 ,  0.9422157E-02 ,  0.9364453E-02 ,  0.9307345E-02 ,  0.9251029E-02 ,&
        &  0.9195290E-02 ,  0.9140226E-02 ,  0.9085805E-02 ,  0.9031995E-02 ,  0.8978802E-02 ,&
        &  0.8926218E-02 ,  0.8874246E-02 ,  0.8822831E-02 ,  0.8772007E-02 ,  0.8721734E-02 ,&
        &  0.8672044E-02 ,  0.8622875E-02 ,  0.8574254E-02 ,  0.8526157E-02 ,  0.8478578E-02 ,&
        &  0.8431531E-02 ,  0.8384986E-02 ,  0.8338917E-02 ,  0.8293363E-02 ,  0.8248281E-02 ,&
        &  0.8203656E-02 ,  0.8159554E-02 ,  0.8115833E-02 ,  0.8072644E-02 ,  0.8029863E-02 ,&
        &  0.7987492E-02 ,  0.7945571E-02 ,  0.7904110E-02 ,  0.7863070E-02 ,  0.7822440E-02 ,&
        &  0.7782107E-02 ,  0.7742366E-02 ,  0.7702901E-02 ,  0.7663889E-02 ,  0.7625191E-02 ,&
        &  0.7586881E-02 ,  0.7548926E-02 ,  0.7511328E-02 ,  0.7474182E-02 ,  0.7437435E-02 ,&
        &  0.7400930E-02 ,  0.7364721E-02 ,  0.7328955E-02 ,  0.7293573E-02 ,  0.7258368E-02 /)
      absice4(:, 16) = (/ &
        &  0.1371110E+00 ,  0.1263437E+00 ,  0.1138219E+00 ,  0.1031942E+00 ,  0.9464203E-01 ,&
        &  0.8767045E-01 ,  0.8183866E-01 ,  0.7683717E-01 ,  0.7245801E-01 ,  0.6856512E-01 ,&
        &  0.6507240E-01 ,  0.6191926E-01 ,  0.5906262E-01 ,  0.5646807E-01 ,  0.5410580E-01 ,&
        &  0.5194927E-01 ,  0.4997518E-01 ,  0.4816421E-01 ,  0.4649395E-01 ,  0.4495197E-01 ,&
        &  0.4352304E-01 ,  0.4219576E-01 ,  0.4095744E-01 ,  0.3979897E-01 ,  0.3871448E-01 ,&
        &  0.3769469E-01 ,  0.3673410E-01 ,  0.3582698E-01 ,  0.3496906E-01 ,  0.3415615E-01 ,&
        &  0.3338472E-01 ,  0.3265037E-01 ,  0.3195064E-01 ,  0.3128339E-01 ,  0.3064628E-01 ,&
        &  0.3003615E-01 ,  0.2945239E-01 ,  0.2889207E-01 ,  0.2835424E-01 ,  0.2783764E-01 ,&
        &  0.2734070E-01 ,  0.2686233E-01 ,  0.2640112E-01 ,  0.2595654E-01 ,  0.2552762E-01 ,&
        &  0.2511272E-01 ,  0.2471195E-01 ,  0.2432435E-01 ,  0.2394912E-01 ,  0.2358571E-01 ,&
        &  0.2323346E-01 ,  0.2289199E-01 ,  0.2256062E-01 ,  0.2223899E-01 ,  0.2192652E-01 ,&
        &  0.2162318E-01 ,  0.2132805E-01 ,  0.2104083E-01 ,  0.2076170E-01 ,  0.2048975E-01 ,&
        &  0.2022498E-01 ,  0.1996714E-01 ,  0.1971574E-01 ,  0.1947053E-01 ,  0.1923143E-01 ,&
        &  0.1899821E-01 ,  0.1877081E-01 ,  0.1854839E-01 ,  0.1833139E-01 ,  0.1811944E-01 ,&
        &  0.1791232E-01 ,  0.1770979E-01 ,  0.1751175E-01 ,  0.1731830E-01 ,  0.1712875E-01 ,&
        &  0.1694354E-01 ,  0.1676205E-01 ,  0.1658457E-01 ,  0.1641057E-01 ,  0.1624031E-01 ,&
        &  0.1607347E-01 ,  0.1590993E-01 ,  0.1574957E-01 ,  0.1559250E-01 ,  0.1543850E-01 ,&
        &  0.1528754E-01 ,  0.1513927E-01 ,  0.1499397E-01 ,  0.1485124E-01 ,  0.1471122E-01 ,&
        &  0.1457379E-01 ,  0.1443890E-01 ,  0.1430642E-01 ,  0.1417625E-01 ,  0.1404839E-01 ,&
        &  0.1392277E-01 ,  0.1379936E-01 ,  0.1367808E-01 ,  0.1355890E-01 ,  0.1344163E-01 ,&
        &  0.1332644E-01 ,  0.1321315E-01 ,  0.1310174E-01 ,  0.1299219E-01 ,  0.1288432E-01 ,&
        &  0.1277823E-01 ,  0.1267383E-01 ,  0.1257107E-01 ,  0.1246995E-01 ,  0.1237042E-01 ,&
        &  0.1227242E-01 ,  0.1217589E-01 ,  0.1208087E-01 ,  0.1198723E-01 ,  0.1189508E-01 ,&
        &  0.1180423E-01 ,  0.1171478E-01 ,  0.1162661E-01 ,  0.1153972E-01 ,  0.1145410E-01 ,&
        &  0.1136972E-01 ,  0.1128655E-01 ,  0.1120451E-01 ,  0.1112365E-01 ,  0.1104392E-01 ,&
        &  0.1096532E-01 ,  0.1088778E-01 ,  0.1081130E-01 ,  0.1073585E-01 ,  0.1066146E-01 ,&
        &  0.1058802E-01 ,  0.1051558E-01 ,  0.1044411E-01 ,  0.1037357E-01 ,  0.1030394E-01 ,&
        &  0.1023524E-01 ,  0.1016739E-01 ,  0.1010042E-01 ,  0.1003431E-01 ,  0.9969062E-02 ,&
        &  0.9904598E-02 ,  0.9840964E-02 ,  0.9778117E-02 ,  0.9716044E-02 ,  0.9654765E-02 ,&
        &  0.9594198E-02 ,  0.9534383E-02 ,  0.9475265E-02 ,  0.9416858E-02 ,  0.9359164E-02 ,&
        &  0.9302156E-02 ,  0.9245805E-02 ,  0.9190131E-02 ,  0.9135040E-02 ,  0.9080692E-02 ,&
        &  0.9026886E-02 ,  0.8973733E-02 ,  0.8921186E-02 ,  0.8869237E-02 ,  0.8817863E-02 ,&
        &  0.8767076E-02 ,  0.8716858E-02 ,  0.8667196E-02 ,  0.8618085E-02 ,  0.8569493E-02 ,&
        &  0.8521454E-02 ,  0.8473919E-02 ,  0.8426909E-02 ,  0.8380398E-02 ,  0.8334377E-02 ,&
        &  0.8288868E-02 ,  0.8243834E-02 ,  0.8199251E-02 ,  0.8155166E-02 ,  0.8111528E-02 ,&
        &  0.8068328E-02 ,  0.8025627E-02 ,  0.7983285E-02 ,  0.7941457E-02 ,  0.7900015E-02 ,&
        &  0.7858957E-02 ,  0.7818337E-02 ,  0.7778154E-02 ,  0.7738378E-02 ,  0.7698989E-02 ,&
        &  0.7659886E-02 ,  0.7621341E-02 ,  0.7583064E-02 ,  0.7545216E-02 ,  0.7507672E-02 ,&
        &  0.7470503E-02 ,  0.7433670E-02 ,  0.7397174E-02 ,  0.7361118E-02 ,  0.7325431E-02 ,&
        &  0.7289980E-02 ,  0.7254812E-02 ,  0.7220077E-02 ,  0.7185708E-02 ,  0.7151500E-02 /)

      end subroutine lwcldpr

end module rrtmg_lw_init
