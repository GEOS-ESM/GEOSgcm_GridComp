module evap_subl_pdf_loop

    use Process_Library_standalone
    use MAPL_ConstantsMod
    use GEOS_UtilsMod

    implicit none

    public evap_subl_pdf_loop_standalone

    private

    contains

    subroutine evap_subl_pdf_loop_standalone(IM, JM, LM, do_qa, USE_BERGERON, PDFSHAPE, DT_MOIST, CCW_EVAP_EFF, &
        CCI_EVAP_EFF, dw_ocean, dw_land, turnrhcrit, CNV_FRC, SRF_TYPE, EIS, KLCL, AREA, &
        PLmb, PLEmb, ZL0, NACTL, NACTI, WHL, WQT, HL2, &
        QT2, HLQT, W3, W2, QT3, HL3, EDMF_FRC, QST3, PDFITERS, WTHV2, WQL, &
        Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN, PDF_A, SUBLC, RHX, EVAPC, RHCRIT3D)

        implicit none

        real, parameter :: mapl_undef = 1.0e15  ! NOTE : This is the value pulled from MAPL_Mod

        integer, intent(in) :: IM, JM, LM

        logical, intent(IN) :: do_qa, USE_BERGERON

        integer, intent(IN) :: PDFSHAPE
        real, intent(IN)    :: DT_MOIST, CCW_EVAP_EFF, CCI_EVAP_EFF, dw_ocean, dw_land
        real, intent(INOUT) :: turnrhcrit

        integer, dimension(:,:),   intent(IN) :: KLCL
        real, dimension(:,:),   intent(IN) :: CNV_FRC, SRF_TYPE, EIS, AREA
        real, dimension(:,:,:), intent(IN) :: PLmb, ZL0, NACTL, NACTI, WHL, WQT, HL2, &
            QT2, HLQT, W3, W2, QT3, HL3, EDMF_FRC, &
            QST3
        real, dimension(:,:,:), pointer, intent(INOUT) :: RHCRIT3D
        real, dimension(1:IM,1:JM,0:LM), intent(IN) :: PLEmb
        real, dimension(:,:,:), intent(INOUT) ::  Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN, PDF_A
        real, dimension(:,:,:), intent(OUT)   ::  WTHV2, WQL, PDFITERS, SUBLC, RHX, EVAPC

        ! Local variables
        integer :: I, J, L, II, JJ
        real    :: facEIS, minrhcrit
        real    :: ALPHA, RHCRIT, EVAPC_C, SUBLC_C

        call update_ESTBLX_to_GPU

        ! evap/subl/pdf
        !$acc parallel loop gang vector collapse(3) private(facEIS, minrhcrit, RHCRIT, &
        !$acc                                               EVAPC_C, SUBLC_C, ALPHA, I, J, L, II, JJ)
        do L=1,LM
            do J=1,JM
                do I=1,IM
                    ! Send the condensates through the pdf after convection
                    facEIS = MAX(0.0,MIN(1.0,EIS(I,J)/10.0))**2
                    ! determine combined minrhcrit in stable/unstable regimes
                    minrhcrit  = (1.0-dw_ocean)*(1.0-facEIS) + (1.0-dw_land)*facEIS
                    if (turnrhcrit <= 0.0) then
                        check : do JJ = 1, J
                            do II = 1, I
                                turnrhcrit  = PLmb(II, JJ, KLCL(II,JJ)) - 250.0 ! 250mb above the LCL
                                if (turnrhcrit > 0.0) exit check 
                            enddo
                        enddo check
                    endif
                    ! Use Slingo-Ritter (1985) formulation for critical relative humidity
                    RHCRIT = 1.0
                    ! lower turn from maxrhcrit=1.0
                    if (PLmb(i,j,l) .le. turnrhcrit) then
                        RHCRIT = minrhcrit
                    else
                        if (L.eq.LM) then
                            RHCRIT = 1.0
                        else
                            RHCRIT = minrhcrit + (1.0-minrhcrit)/(19.) * &
                                    ((atan( (2.*(PLmb(i,j,l)-turnrhcrit)/(PLEmb(i,j,LM)-turnrhcrit)-1.) * &
                                    tan(20.*MAPL_PI/21.-0.5*MAPL_PI) ) + 0.5*MAPL_PI) * 21./MAPL_PI - 1.)
                        endif
                    endif
                    ! include grid cell area scaling and limit RHcrit to > 70% 
                    ALPHA = max(0.0,min(0.30, (1.0-RHCRIT)*SQRT(SQRT(AREA(I,J)/1.e10)) ) )
                    ! fill RHCRIT export
                    if (associated(RHCRIT3D)) RHCRIT3D(I,J,L) = 1.0-ALPHA
                    ! Put condensates in touch with the PDF
                    if (.not. do_qa) then ! if not doing cloud pdf inside of GFDL-MP 
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
                        RHX(I,J,L) = Q(I,J,L)/GEOS_QSAT( T(I,J,L), PLmb(I,J,L) )
                        ! meltfrz new condensates
                        call MELTFRZ ( DT_MOIST     , &
                                       CNV_FRC(I,J) , &
                                       SRF_TYPE(I,J), &
                                       T(I,J,L)     , &
                                       QLCN(I,J,L)  , &
                                       QICN(I,J,L) )
                        call MELTFRZ ( DT_MOIST     , &
                                       CNV_FRC(I,J) , &
                                       SRF_TYPE(I,J), &
                                       T(I,J,L)     , &
                                       QLLS(I,J,L)  , &
                                       QILS(I,J,L) )
                    endif
                    ! evaporation for CN
                    if (CCW_EVAP_EFF > 0.0) then ! else evap done inside GFDL
                        RHCRIT = 1.0
                        EVAPC(I,J,L) = Q(I,J,L)
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
                        ! EVAPC(I,J,L) = ( Q(I,J,L) - EVAPC(I,J,L) ) / DT_MOIST
                        ! Note : I'm not sure why the above calculation always results in 0.
                        !        I'm explicitly passing the difference between the Q and EVAPC from EVAP3
                        !        into EVAPC_C
                        EVAPC(I,J,L) = EVAPC_C / DT_MOIST
                    endif
                    ! sublimation for CN
                    if (CCI_EVAP_EFF > 0.0) then ! else subl done inside GFDL
                        RHCRIT = 1.0 - ALPHA
                        SUBLC(I,J,L) = Q(I,J,L)
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
                                QST3(I,J,L), &
                                SUBLC_C  )
                        ! SUBLC(I,J,L) = ( Q(I,J,L) - SUBLC(I,J,L) ) / DT_MOIST
                                ! Note : I'm not sure why the above calculation always results in 0.
                        !        I'm explicitly passing the difference between the Q and SUBLC from SUBL3
                        !        into SUBLC_C
                        SUBLC(I,J,L) = SUBLC_C / DT_MOIST
                    endif
                    ! cleanup clouds
                    call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L) )
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
