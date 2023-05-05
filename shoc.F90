module shoc

! Implementation of the Simplified High Order Closure (SHOC) scheme 
! of Bogenschutz and Krueger (2013), J. Adv. Model. Earth Syst, 5, 195-211,
! doi: 10.1002/jame.200118. (further referred to as BK13)
! in a single column form suitable for use in a GCM physics package. 
! Alex Belochitski, heavily based on the code of Peter Bogenschutz.
! S Moorthi - optimization, cleanup, improve and customize for gsm
!
! N Arnold - implemented in GEOS
!

 use MAPL_ConstantsMod, only: ggr    => MAPL_GRAV,   &
                              cp     => MAPL_CP,     &
                              rgas   => MAPL_RGAS,   &
                              rv     => MAPL_RVAP,   &
                              lcond  => MAPL_ALHL,   &
                              lfus   => MAPL_ALHS,   &
                              pi     => MAPL_PI,     &
                              MAPL_H2OMW, MAPL_AIRMW

!  use MAPL,              only: MAPL_UNDEF

! use MAPL_SatVaporMod,  only: MAPL_EQsat 

! use shocparams

 implicit none

 private

 public update_moments

 contains

 subroutine update_moments( IM, JM, LM, & ! in
                             DT,      &  ! in
                             SH,      &  ! in
                             EVAP,    &  ! in
                             ZL,      &  ! in
                             KH,      &  ! in
                             TKE,     &  ! in
                             ISOTROPY, &  ! in
                             QT,      &  ! in
                             HL,      &  ! in
                             MFFRC,   &  ! in
                             MFQT2,   &  ! in
                             MFQT3,   &  ! in
                             MFHL2,   &  ! in
                             MFHL3,   &  ! in
                             MFW2,    &  ! in
                             MFW3,    &  ! in
                             MFWQT,   &  ! in
                             MFWHL,   &  ! in
                             MFHLQT,  &  ! in
                             qt2,     &  ! inout
                             qt3,     &  ! inout
                             hl2,     &  ! out
                             hl3,     &  ! out
                             w2,      &  ! out
                             w3,      &  ! out
                             wqt,     &  ! out
                             whl,     &  ! out
                             hlqt,    &  ! out
                           hl2tune,   &  ! tuning parameters
                           qt2tune,   &
                           hlqt2tune, &
                           qt2scale,  &
                           qt3_tscale )


    integer, intent(in   ) :: IM, JM, LM       ! dimensions
    real,    intent(in   ) :: DT               ! timestep [s]
    real,    intent(in   ) :: SH   (IM,JM)     ! surface sensible heat flux
    real,    intent(in   ) :: EVAP (IM,JM)     ! surface evaporation
    real,    intent(in   ) :: ZL   (IM,JM,LM)  ! heights [m]
    real,    intent(in   ) :: KH   (IM,JM,0:LM)  ! diffusivity
    real,    intent(in   ) :: TKE  (IM,JM,LM)  ! turbulent kinetic energy
    real,    intent(in   ) :: ISOTROPY(IM,JM,LM)  ! isotropy timescale
    real,    intent(in   ) :: QT   (IM,JM,LM)  ! total water
    real,    intent(in   ) :: HL   (IM,JM,LM)  ! liquid water static energy
    real,    intent(in   ) :: MFFRC(IM,JM,LM)  ! mass flux area fraction 
    real,    intent(in   ) :: MFQT2(IM,JM,LM)  ! 
    real,    intent(in   ) :: MFQT3(IM,JM,LM)  ! 
    real,    intent(in   ) :: MFHL2(IM,JM,LM)  ! 
    real,    intent(in   ) :: MFHL3(IM,JM,LM)  ! 
    real,    intent(in   ) :: MFW2 (IM,JM,LM)  ! 
    real,    intent(in   ) :: MFW3 (IM,JM,LM)  ! 
    real,    intent(in   ) :: MFWQT(IM,JM,LM)  ! 
    real,    intent(in   ) :: MFWHL(IM,JM,LM)  ! 
    real,    intent(in   ) :: MFHLQT(IM,JM,LM) ! 
    real,    intent(inout) :: qt2  (IM,JM,LM)  ! total water variance
    real,    intent(inout) :: qt3  (IM,JM,LM)  ! third moment of total water
    real,    intent(  out) :: hl2  (IM,JM,LM)  ! liquid water static energy variance
    real,    intent(  out) :: hl3  (IM,JM,LM)  ! third moment static energy
    real,    intent(  out) :: w2   (IM,JM,LM)  ! vertical velocity variance
    real,    intent(  out) :: w3   (IM,JM,LM)  ! third moment vertical velocity
    real,    intent(  out) :: wqt  (IM,JM,LM)  ! vertical flux of total water 
    real,    intent(  out) :: whl  (IM,JM,LM)  ! vertical flux of liquid water static energy
    real,    intent(  out) :: hlqt (IM,JM,LM)  ! total water, static energy covariance
    real,    intent(in   ) :: HL2TUNE,     &   ! tuning parameters
                              HLQT2TUNE,   &
                              QT2SCALE,    &
                              QT2TUNE,     &
                              QT3_TSCALE   
   
    real, parameter :: HL2MIN = 0.0025
    real, parameter :: HL2MAX = 2.0

    ! Local variables
    integer :: i, j, k, kd, ku
    real :: wrk1, wrk2, wrk3, sm, onemmf
    real, dimension(IM,JM,0:LM) :: qt2_edge, &
                                   hl2_edge, &
                                   wqt_edge, &
                                   whl_edge, &
                                   hlqt_edge,&
                                   qtgrad

!$acc data present(SH, EVAP, ZL, KH, TKE, ISOTROPY, QT, HL, MFFRC, MFQT2, MFQT3) &
!$acc      present(MFHL2, MFHL3, MFW2, MFW3, MFWQT, MFWHL, MFHLQT) &
!$acc      present(qt2, qt3) &
!$acc      present(hl2, hl3, w2, w3, wqt, whl, hlqt) &
!$acc      create(qt2_edge, hl2_edge, wqt_edge, whl_edge, hlqt_edge, qtgrad)

!$acc parallel
!$acc loop gang vector collapse(3) private(wrk1, wrk2, wrk3, sm)
    ! define resolved gradients on edges
    do k=1,LM-1
        do j = 1,JM
            do i = 1,IM
                wrk1 = 1.0 / (ZL(i,j,k)-ZL(i,j,k+1))
                wrk3 = KH(i,j,k) * wrk1

                sm   = 0.5*(ISOTROPY(i,j,k)+ISOTROPY(i,j,k+1))*wrk1*wrk3 !Tau*Kh/dz^2

                ! SGS vertical flux liquid/ice water static energy. Eq 1 in BK13                                                        
                wrk1            = HL(i,j,k) - HL(i,j,k+1)
                whl_edge(i,j,k) = - wrk3 * wrk1

                ! SGS vertical flux of total water. Eq 2 in BK13                                                                        
                wrk2            = QT(i,j,k) - QT(i,j,k+1)
                wqt_edge(i,j,k) = - wrk3 * wrk2

                ! Second moment of liquid/ice water static energy. Eq 4 in BK13 
                hl2_edge(i,j,k) = HL2TUNE * sm * wrk1 * wrk1

                ! Second moment of total water mixing ratio.  Eq 3 in BK13
                qtgrad(i,j,k)   = wrk2 / (ZL(i,j,k)-ZL(i,j,k+1))
                qt2_edge(i,j,k) = KH(i,j,k)*qtgrad(i,j,k)**2

                ! Covariance of total water mixing ratio and liquid/ice water static energy.  Eq 5 in BK13
                hlqt_edge(i,j,k) = HLQT2TUNE * sm * wrk1 * wrk2
            enddo
        enddo
    end do
!$acc end parallel

    ! set lower boundary conditions
!$acc parallel
!$acc loop gang vector collapse(2)
    do j = 1,JM
        do i = 1,IM
            whl_edge(i,j,LM)  = SH(i,j)/cp
            wqt_edge(i,j,LM)  = EVAP(i,j)
            hl2_edge(i,j,LM)  = hl2_edge(i,j,LM-1)
            qt2_edge(i,j,LM)  = qt2_edge(i,j,LM-1)
            hlqt_edge(i,j,LM) = hlqt_edge(i,j,LM-1)
            qtgrad(i,j,LM)    = qtgrad(i,j,LM-1)
            qtgrad(i,j,0)     = qtgrad(i,j,1)
        enddo
    enddo

!$acc loop gang vector collapse(3) private(kd, ku, onemmf, wrk1, wrk2, wrk3)
    do k=1,LM
        do j = 1,JM
            do i = 1,IM
                kd = k-1
                ku = k
                if (k==1) kd = k

                onemmf = 1.0 - MFFRC(i,j,k)

                w2(i,j,k) = onemmf*0.667*TKE(i,j,k) + MFW2(i,j,k)

                hl2(i,j,k) = onemmf*0.5*( hl2_edge(i,j,kd) + hl2_edge(i,j,ku) ) + MFHL2(i,j,k)

                wrk1 = 0.5*(qt2_edge(i,j,kd)+qt2_edge(i,j,ku))              ! averaging ED gradient production term
                wrk2 = 0.5*MFWQT(i,j,k)*0.5*(qtgrad(i,j,kd)+qtgrad(i,j,ku)) ! MF gradient production term
                qt2(i,j,k) = qt2(i,j,k) + DT*(wrk1-wrk2) 

                wrk3 = QT2TUNE*sqrt(0.01+TKE(i,j,k))/(QT2SCALE*0.4*ZL(i,j,k)/(0.4*ZL(i,j,k)+QT2SCALE))
                qt2(i,j,k) = qt2(i,j,k) / (1. + DT*wrk3)

                hlqt(i,j,k) = onemmf*0.5*( hlqt_edge(i,j,kd) + hlqt_edge(i,j,ku) ) + MFHLQT(i,j,k)

                wqt(i,j,k)  = onemmf*0.5*( wqt_edge(i,j,kd) + wqt_edge(i,j,ku) ) + MFWQT(i,j,k)

                whl(i,j,k)  = onemmf*0.5*( whl_edge(i,j,kd) + whl_edge(i,j,ku) ) + MFWHL(i,j,k)

                ! Restrict QT variance, 1-25% of total water.
                qt2(i,j,k) = max(min(qt2(i,j,k),(0.25*QT(i,j,k))**2),(0.01*QT(i,j,k))**2)
                hl2(i,j,k) = max(min(hl2(i,j,k),HL2MAX),HL2MIN)

                ! Ensure realizibility
        !        hl2 = max(hl2,whl*whl/max(w2,0.1))
        !        qt2 = max(qt2,wqt*wqt/max(w2,0.1))
                hlqt(i,j,k) = sign( min( abs(hlqt(i,j,k)), sqrt(hl2(i,j,k)*qt2(i,j,k)) ), hlqt(i,j,k) )

                qt3(i,j,k) = max( MFQT3(i,j,k), qt3(i,j,k)*max(1.-DT/QT3_TSCALE,0.0) )
                hl3(i,j,k) = MFHL3(i,j,k)
                w3(i,j,k)  = MFW3(i,j,k)
            enddo
        enddo
    end do

!$acc end parallel

!$acc end data
 end subroutine update_moments

end module shoc

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