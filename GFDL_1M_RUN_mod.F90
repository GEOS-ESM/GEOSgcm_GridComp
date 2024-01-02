module GFDL_1M_RUN_module

    implicit none

    private

    ! Internals
    real, pointer, dimension(:,:,:), public :: Q, QLLS, QLCN, CLLS, CLCN, QILS, QICN, QRAIN, QSNOW, QGRAUPEL
    real, pointer, dimension(:,:,:), public :: NACTL, NACTI
    ! Imports
    real, pointer, dimension(:,:,:), public :: ZLE, PLE, T, U, V, W, KH
    real, pointer, dimension(:,:), public   :: AREA, FRLAND, TS, DTSX, SH, EVAP, KPBLSC
    real, pointer, dimension(:,:,:), public :: HL2, HL3, QT2, QT3, W2, W3, HLQT, WQT, WQL, WHL, EDMF_FRC
    real, pointer, dimension(:,:,:), public :: WTHV2
    real, pointer, dimension(:,:,:), public :: OMEGA

    real, pointer, dimension(:,:  ), public :: PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL
    real, pointer, dimension(:,:  ), public :: LS_PRCP, LS_SNR, ICE, FRZR, CNV_FRC, SRF_TYPE
    real, pointer, dimension(:,:,:), public :: DQVDT_macro, DQIDT_macro, DQLDT_macro, DQADT_macro, DQRDT_macro, DQSDT_macro, DQGDT_macro
    real, pointer, dimension(:,:,:), public ::  DUDT_macro,  DVDT_macro,  DTDT_macro
    real, pointer, dimension(:,:,:), public :: DQVDT_micro, DQIDT_micro, DQLDT_micro, DQADT_micro, DQRDT_micro, DQSDT_micro, DQGDT_micro
    real, pointer, dimension(:,:,:), public ::  DUDT_micro,  DVDT_micro,  DTDT_micro
    real, pointer, dimension(:,:,:), public :: RAD_CF, RAD_QV, RAD_QL, RAD_QI, RAD_QR, RAD_QS, RAD_QG
    real, pointer, dimension(:,:,:), public :: CLDREFFL, CLDREFFI
    real, pointer, dimension(:,:,:), public :: EVAPC, SUBLC
    real, pointer, dimension(:,:,:), public :: RHX, REV_LS, RSU_LS
    real, pointer, dimension(:,:,:), public :: PFL_LS, PFL_AN
    real, pointer, dimension(:,:,:), public :: PFI_LS, PFI_AN
    real, pointer, dimension(:,:,:), public :: PDF_A, PDFITERS
    real, pointer, dimension(:,:,:), public :: RHCRIT3D
    real, pointer, dimension(:,:), public   :: EIS, LTS
    real, pointer, dimension(:,:), public   :: DBZ_MAX, DBZ_1KM, DBZ_TOP, DBZ_M10C
    real, pointer, dimension(:,:,:), public :: PTR3D
    real, pointer, dimension(:,:  ), public :: PTR2D

    real*8, public :: DT_R8
    
    ! Local resource variables
    real, public    :: TURNRHCRIT
    real, public    :: CCW_EVAP_EFF
    real, public    :: CCI_EVAP_EFF
    integer, public :: PDFSHAPE
    real, public    :: ANV_ICEFALL 
    real, public    :: LS_ICEFALL
    real, public   :: FAC_RL
    real, public    :: MIN_RL
    real, public    :: MAX_RL
    real, public    :: FAC_RI
    real, public    :: MIN_RI
    real, public    :: MAX_RI
    logical, public :: LHYDROSTATIC
    logical, public :: LPHYS_HYDROSTATIC
    logical, public :: do_qa
 
    public :: data_setup

    contains

    subroutine data_setup(IM, JM, LM, dirName, rank_str)

        implicit none

        integer, intent(IN) :: IM, JM, LM

        character*100, intent(IN) :: dirName, rank_str

        integer :: fileID

        ! Imports and Internal array allocations
        allocate(Q(IM, JM, LM))
        allocate(QRAIN(IM, JM, LM))
        allocate(QSNOW(IM, JM, LM))
        allocate(QGRAUPEL(IM, JM, LM))
        allocate(QLLS(IM, JM, LM))
        allocate(QLCN(IM, JM, LM))
        allocate(CLCN(IM, JM, LM))
        allocate(CLLS(IM, JM, LM))
        allocate(QILS(IM, JM, LM))
        allocate(QICN(IM, JM, LM))
        allocate(NACTL(IM, JM, LM))
        allocate(NACTI(IM, JM, LM))
        allocate(AREA(IM, JM))
        allocate(ZLE(IM, JM, 0:LM))
        allocate(PLE(IM, JM, 0:LM))
        allocate(T(IM, JM, LM))
        allocate(U(IM, JM, LM))
        allocate(V(IM, JM, LM))
        allocate(W(IM, JM, LM))
        allocate(FRLAND(IM, JM))
        allocate(KH(IM, JM, 0:LM))
        allocate(EDMF_FRC(IM, JM, LM))
        allocate(W2(IM, JM, LM))
        allocate(W3(IM, JM, LM))
        allocate(WQT(IM, JM, LM))
        allocate(WHL(IM, JM, LM))
        allocate(HL2(IM, JM, LM))
        allocate(HL3(IM, JM, LM))
        allocate(QT2(IM, JM, LM))
        allocate(QT3(IM, JM, LM))
        allocate(HLQT(IM, JM, LM))
        allocate(TS(IM, JM))
        allocate(KPBLSC(IM, JM))
        allocate(SH(IM, JM))
        allocate(EVAP(IM, JM))
        allocate(OMEGA(IM, JM, LM))

        ! Export array allocation
        allocate(RAD_CF(IM, JM, LM))
        allocate(RAD_QV(IM, JM, LM))
        allocate(RAD_QL(IM, JM, LM))
        allocate(RAD_QI(IM, JM, LM))
        allocate(RAD_QR(IM, JM, LM))
        allocate(RAD_QS(IM, JM, LM))
        allocate(RAD_QG(IM, JM, LM))
        allocate(CLDREFFL(IM, JM, LM))
        allocate(CLDREFFI(IM, JM, LM))
        allocate(CNV_FRC(IM, JM))
        allocate(SRF_TYPE(IM, JM))
        allocate(EVAPC(IM, JM, LM))
        allocate(SUBLC(IM, JM, LM))
        allocate(PRCP_RAIN(IM, JM))
        allocate(PRCP_SNOW(IM, JM))
        allocate(PRCP_ICE(IM, JM))
        allocate(PRCP_GRAUPEL(IM, JM))
        allocate(LS_PRCP(IM, JM))
        allocate(LS_SNR(IM, JM))
        allocate(ICE(IM, JM))
        allocate(FRZR(IM, JM))
        allocate(RHX(IM, JM, LM))
        allocate(REV_LS(IM, JM, LM))
        allocate(RSU_LS(IM, JM, LM))
        allocate(PFL_AN(IM, JM, 0:LM))
        allocate(PFL_LS(IM, JM, 0:LM))
        allocate(PFI_AN(IM, JM, 0:LM))
        allocate(PFI_LS(IM, JM, 0:LM))
        allocate(PDF_A(IM, JM, LM))
        allocate(WTHV2(IM, JM, LM))
        allocate(WQL(IM, JM, LM))
        allocate(PDFITERS(IM, JM, LM))
        allocate(LTS(IM, JM))
        allocate(EIS(IM, JM))
        allocate(DQVDT_macro(IM, JM, LM))
        allocate(DQIDT_macro(IM, JM, LM))
        allocate(DQLDT_macro(IM, JM, LM))
        allocate(DQADT_macro(IM, JM, LM))
        allocate(DQRDT_macro(IM, JM, LM))
        allocate(DQSDT_macro(IM, JM, LM))
        allocate(DQGDT_macro(IM, JM, LM))
        allocate(DUDT_macro(IM, JM, LM))
        allocate(DVDT_macro(IM, JM, LM))
        allocate(DTDT_macro(IM, JM, LM))
        allocate(RHCRIT3D(IM, JM, LM))
        allocate(DQVDT_micro(IM, JM, LM))
        allocate(DQIDT_micro(IM, JM, LM))
        allocate(DQLDT_micro(IM, JM, LM))
        allocate(DQADT_micro(IM, JM, LM))
        allocate(DQRDT_micro(IM, JM, LM))
        allocate(DQSDT_micro(IM, JM, LM))
        allocate(DQGDT_micro(IM, JM, LM))
        allocate(DUDT_micro(IM, JM, LM))
        allocate(DVDT_micro(IM, JM, LM))
        allocate(DTDT_micro(IM, JM, LM))

        open(newunit=fileID, file=trim(dirName) // '/Q_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) Q
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QRAIN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QRAIN
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QSNOW_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QSNOW
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QGRAUPEL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QGRAUPEL
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QLLS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QLLS
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QLCN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QLCN
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CLCN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CLCN
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CLLS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CLLS
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QILS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QILS
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QICN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QICN
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/NACTL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) NACTL
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/NACTI_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) NACTI
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/AREA_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) AREA
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/ZLE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ZLE
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PLE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) PLE
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/T_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) T
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/U_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) U
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/V_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) V
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/W_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) W
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/FRLAND_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) FRLAND
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/KH_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) KH
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/EDMF_FRC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) EDMF_FRC
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/W2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) W2
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/W3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) W3
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/WQT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) WQT
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/WHL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) WHL
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/HL2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) HL2
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/HL3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) HL3
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QT2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QT2
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QT3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QT3
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/HLQT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) HLQT
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/TS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) TS
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/KPBLSC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) KPBLSC
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/SH_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) SH
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/EVAP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) EVAP
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/OMEGA_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) OMEGA
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DT_R8_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DT_R8
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/TURNRHCRIT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) TURNRHCRIT
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CCW_EVAP_EFF_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CCW_EVAP_EFF
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CCI_EVAP_EFF_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CCI_EVAP_EFF
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PDFSHAPE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) PDFSHAPE
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/ANV_ICEFALL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ANV_ICEFALL
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/LS_ICEFALL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) LS_ICEFALL
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/FAC_RL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) FAC_RL
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/MIN_RL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) MIN_RL
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/MAX_RL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) MAX_RL
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/FAC_RI_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) FAC_RI
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/MIN_RI_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) MIN_RI
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/MAX_RI_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) MAX_RI
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/LHYDROSTATIC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) LHYDROSTATIC
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/LPHYS_HYDROSTATIC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) LPHYS_HYDROSTATIC
        close(fileID)

        ! Note : do_qa is a module variable in gfdl_cloud_microphys.F90
        do_qa = .false.

        !$acc data copyin(DT_MOIST, PDFSHAPE, CNV_FRC, SRF_TYPE, PLmb, ZL0, &
        !$acc             NACTL, NACTI, WHL, WQT, HL2, QT2, &
        !$acc             HLQT, W3, W2, QT3, HL3, EDMF_FRC, USE_BERGERON, &
        !$acc             QST3, CCI_EVAP_EFF, CCW_EVAP_EFF, &
        !$acc             PLEmb) &
        !$acc      copy(Q, QLLS, QLCN, QILS, QICN, T, CLLS, &
        !$acc           CLCN, PDF_A, WTHV2, WQL) &
        !$acc      copyout(EVAPC,RHX, SUBLC, PDFITERS)


        !$acc end data


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