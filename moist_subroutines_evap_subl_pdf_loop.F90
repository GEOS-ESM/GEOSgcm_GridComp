module evap_subl_pdf_loop

    use Process_Library_standalone
    use MAPL_ConstantsMod
    use GEOS_UtilsMod

    implicit none

    public evap_subl_pdf_loop_standalone

    private

    contains

    subroutine evap_subl_pdf_loop_standalone(IM, JM, LM, do_qa, USE_BERGERON, PDFSHAPE, DT_MOIST, CCW_EVAP_EFF, &
        CCI_EVAP_EFF, maxrhcrit2D, turnrhcrit2D, minrhcrit2D, TROPP, CNV_FRC, SRF_TYPE, &
        PLmb, PLEmb, ZL0, NACTL, NACTI, WHL, WQT, HL2, &
        QT2, HLQT, W3, W2, QT3, HL, HL3, EDMF_FRC, QST3, PDFITERS, WTHV2, WQL, &
        Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN, PDF_A, SUBLC, RHX, EVAPC)

        implicit none

        real, parameter :: mapl_undef = 1.0e15  ! NOTE : This is the value pulled from MAPL_Mod

        integer, intent(in) :: IM, JM, LM

        logical, intent(IN) :: do_qa, USE_BERGERON

        integer, intent(IN) :: PDFSHAPE
        real, intent(IN)    :: DT_MOIST, CCW_EVAP_EFF, CCI_EVAP_EFF

        real, dimension(:,:),   intent(IN) :: maxrhcrit2D, turnrhcrit2D, minrhcrit2D, TROPP, CNV_FRC, SRF_TYPE
        real, dimension(:,:,:), intent(IN) :: PLmb, ZL0, NACTL, NACTI, WHL, WQT, HL2, &
            QT2, HLQT, W3, W2, QT3, HL, HL3, EDMF_FRC, &
            QST3
        real, dimension(1:IM,1:JM,0:LM), intent(IN) :: PLEmb
        real, dimension(:,:,:), intent(INOUT) ::  PDFITERS, WTHV2, WQL, Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN, PDF_A
        real, dimension(:,:,:), intent(OUT)   ::  SUBLC, RHX, EVAPC

        ! Local variables
        integer :: I, J, L
        real    :: turnrhcrit_up
        real    :: ALPHAl, ALPHAu, ALPHA, RHCRIT
        real    :: EVAPC_C, SUBLC_C

        call update_ESTBLX_to_GPU

        ! evap/subl/pdf
        !$acc parallel loop gang vector collapse(3) private(ALPHA, ALPHAl, ALPHAu, turnrhcrit_up, RHCRIT, &
        !$acc                                               EVAPC_C, SUBLC_C)
        do L=1,LM
            do J=1,JM
                do I=1,IM
                    if (.not. do_qa) then ! if not doing the evap/subl/pdf inside of GFDL-MP 
                        ! Send the condensates through the pdf after convection
                        !  Use Slingo-Ritter (1985) formulation for critical relative humidity
                        ALPHA = maxrhcrit2D(I,J)
                        ! lower turn from maxrhcrit
                        if (PLmb(i,j,l) .le. turnrhcrit2D(I,J)) then
                            ALPHAl = minrhcrit2D(I,J)
                        else
                            if (L.eq.LM) then
                                ALPHAl = maxrhcrit2D(I,J)
                            else
                                ALPHAl = minrhcrit2D(I,J) + (maxrhcrit2D(I,J)-minrhcrit2D(I,J))/(19.) * &
                                        ((atan( (2.*(PLmb(i,j,l)-turnrhcrit2D(I,J))/min(100., PLEmb(i,j,LM)-turnrhcrit2D(I,J))-1.) * &
                                        tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
                            endif
                        endif
                        ! upper turn back to maxrhcrit
                        turnrhcrit_up = TROPP(i,j)/100.0
                        IF (turnrhcrit_up == MAPL_UNDEF) turnrhcrit_up = 100.
                        if (PLmb(i,j,l) .le. turnrhcrit_up) then
                            ALPHAu = maxrhcrit2D(I,J)
                        else
                            ALPHAu = maxrhcrit2D(I,J) - (maxrhcrit2D(I,J)-minrhcrit2D(I,J))/(19.) * &
                                    ((atan( (2.*(PLmb(i,j,l)-turnrhcrit_up)/( turnrhcrit2D(I,J)-turnrhcrit_up)-1.) * &
                                    tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
                        endif
                        ! combine and limit
                        ALPHA = min( 0.25, 1.0 - min(max(ALPHAl,ALPHAu),1.) ) ! restrict RHcrit to > 75% 
                        ! Put condensates in touch with the PDF
                        call hystpdf( &
                                    DT_MOIST       , &
                                    ALPHA          , &
                                    PDFSHAPE       , &
                                    CNV_FRC(I,J)   , &
                                    SRF_TYPE(I,J)  , &
                                    PLmb(I,J,L)    , &
                                    ZL0(I,J,L)     , &
                                    Q(I,J,L)       , &
                                    QLLS(I,J,L)    , &
                                    QLCN(I,J,L)    , &
                                    QILS(I,J,L)    , &
                                    QICN(I,J,L)    , &
                                    T(I,J,L)       , &
                                    CLLS(I,J,L)    , &
                                    CLCN(I,J,L)    , &
                                    NACTL(I,J,L)   , &
                                    NACTI(I,J,L)   , &
                                    WHL(I,J,L)     , &
                                    WQT(I,J,L)     , &
                                    HL2(I,J,L)     , &
                                    QT2(I,J,L)     , &
                                    HLQT(I,J,L)    , &
                                    W3(I,J,L)      , &
                                    W2(I,J,L)      , &
                                    QT3(I,J,L)     , &
                                    HL3(I,J,L)     , &
                                    EDMF_FRC(I,J,L), &
                                    PDF_A(I,J,L)   , &
                                    PDFITERS(I,J,L), &
                                    WTHV2(I,J,L)   , &
                                    WQL(I,J,L)     , &
                                    .false.        , & 
                                    USE_BERGERON)

                        ! evaporation for CN/LS
                        RHCRIT = 1.0
                        ! EVAPC(I,J,L) = Q(I,J,L)

                        call EVAP3 (         &
                                DT_MOIST      , &
                                CCW_EVAP_EFF  , &
                                RHCRIT        , &
                                PLmb(I,J,L)  , &
                                T(I,J,L)  , &
                                Q(I,J,L)  , &
                                QLCN(I,J,L)  , &
                                QICN(I,J,L)  , &
                                CLCN(I,J,L)  , &
                                NACTL(I,J,L)  , &
                                NACTI(I,J,L)  , &
                                QST3(I,J,L), &
                                EVAPC_C  )

                        EVAPC(I,J,L) = EVAPC_C / DT_MOIST
                        ! EVAPC(I,J,L) = (Q(I,J,L) - EVAPC(I,J,L)) / DT_MOIST
                        ! sublimation for CN/LS
                        RHCRIT = 1.0 - ALPHA
                        ! SUBLC(I,J,L) =   Q(I,J,L)

                        call SUBL3 (        &
                                DT_MOIST      , &
                                CCI_EVAP_EFF  , &
                                RHCRIT        , &
                                PLmb(I,J,L)  , &
                                    T(I,J,L)  , &
                                    Q(I,J,L)  , &
                                QLCN(I,J,L)  , &
                                QICN(I,J,L)  , &
                                CLCN(I,J,L)  , &
                                NACTL(I,J,L)  , &
                                NACTI(I,J,L)  , &
                                QST3(I,J,L),  &
                                SUBLC_C)

                        SUBLC(I,J,L) = SUBLC_C / DT_MOIST
                        ! SUBLC(I,J,L) = (Q(I,J,L) - SUBLC(I,J,L)) / DT_MOIST
                        ! cleanup clouds
                        call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L) )
                        RHX(I,J,L) = Q(I,J,L)/GEOS_QSAT( T(I,J,L), PLmb(I,J,L) )
                    endif
                end do ! IM loop
            end do ! JM loop
        end do ! LM loop
        !$acc end parallel loop
    end subroutine

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
