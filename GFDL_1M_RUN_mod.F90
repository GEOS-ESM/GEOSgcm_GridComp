module GFDL_1M_RUN_module

    use GEOS_UtilsMod
    use MAPL_PhysicalConstantsMod
    use GEOSmoist_Process_Library
    use gfdl2_cloud_microphys_mod

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

    ! Exports
    real, pointer, dimension(:,:  ), public :: PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL
    real, pointer, dimension(:,:  ), public :: LS_PRCP, LS_SNR, ICE, FRZR, CNV_FRC, SRF_TYPE
    real, pointer, dimension(:,:,:), public :: DQVDT_macro, DQIDT_macro, DQLDT_macro, DQADT_macro, DQRDT_macro, DQSDT_macro, DQGDT_macro
    real, pointer, dimension(:,:,:), public :: DUDT_macro,  DVDT_macro,  DTDT_macro
    real, pointer, dimension(:,:,:), public :: DQVDT_micro, DQIDT_micro, DQLDT_micro, DQADT_micro, DQRDT_micro, DQSDT_micro, DQGDT_micro
    real, pointer, dimension(:,:,:), public :: DUDT_micro,  DVDT_micro,  DTDT_micro
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

    ! Module variable from aer_actv_single_moment.F90
    logical, public :: USE_BERGERON
 
    public :: data_setup, GFDL_1M_RUN_driver, compare_results

    contains

    subroutine data_setup(IM, JM, LM, dirName, rank_str)

        implicit none

        integer, intent(IN) :: IM, JM, LM

        character*100, intent(IN) :: dirName, rank_str

        integer :: fileID

        logical :: file_exists

        print*, 'Running data_setup...'

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

        inquire(file=trim(dirName) // '/DBZ_MAX_' // trim(rank_str) // '.out', exist=file_exists)
        if(file_exists) then
            allocate(DBZ_MAX(IM, JM))
        endif

        inquire(file=trim(dirName) // '/DBZ_1KM_' // trim(rank_str) // '.out', exist=file_exists)
        if(file_exists) then
            allocate(DBZ_1KM(IM, JM))
        endif

        inquire(file=trim(dirName) // '/DBZ_TOP_' // trim(rank_str) // '.out', exist=file_exists)
        if(file_exists) then
            allocate(DBZ_TOP(IM, JM))
        endif

        inquire(file=trim(dirName) // '/DBZ_M10C_' // trim(rank_str) // '.out', exist=file_exists)
        if(file_exists) then
            allocate(DBZ_M10C(IM, JM))
        endif

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

        open(newunit=fileID, file=trim(dirName) // '/test/mp_time_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) mp_time
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/t_min_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) t_min
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/t_sub_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) t_sub
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_r2g_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_r2g
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_smlt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_smlt
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_g2r_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_g2r
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/dw_land_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dw_land
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/dw_ocean_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dw_ocean
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/vi_fac_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vi_fac
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/vr_fac_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vr_fac
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/vs_fac_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vs_fac
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/vg_fac_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vg_fac
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/ql_mlt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ql_mlt
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/do_qa_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_qa
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/fix_negative_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) fix_negative
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/vi_max_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vi_max
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/vs_max_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vs_max
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/vg_max_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vg_max
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/vr_max_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) vr_max
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/qs_mlt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qs_mlt
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/qs0_crt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qs0_crt
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/qi_gen_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qi_gen
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/ql0_max_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ql0_max
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/qi0_max_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qi0_max
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/qi0_crt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qi0_crt
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/qr0_crt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qr0_crt
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/fast_sat_adj_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) fast_sat_adj
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/rh_inc_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) rh_inc
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/rh_ins_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) rh_ins
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/rh_inr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) rh_inr
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/const_vi_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) const_vi
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/const_vs_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) const_vs
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/const_vg_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) const_vg
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/const_vr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) const_vr
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/use_ccn_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) use_ccn
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/rthreshu_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) rthreshu
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/rthreshs_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) rthreshs
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/ccn_l_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ccn_l
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/ccn_o_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ccn_o
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/qc_crt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qc_crt
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_g2v_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_g2v
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_v2g_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_v2g
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_s2v_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_s2v
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_v2s_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_v2s
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_revp_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_revp
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_frz_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_frz
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/do_bigg_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_bigg
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/do_evap_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_evap
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/do_subl_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_subl
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/sat_adj0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) sat_adj0
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/c_piacr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_piacr
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_imlt_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_imlt
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_v2l_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_v2l
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_l2v_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_l2v
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_i2v_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_i2v
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_i2s_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_i2s
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/tau_l2r_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) tau_l2r
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/qi_lim_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) qi_lim
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/ql_gen_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ql_gen
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/c_paut_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_paut
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/c_psaci_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_psaci
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/c_pgacs_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_pgacs
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/c_pgaci_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_pgaci
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/z_slope_liq_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) z_slope_liq
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/z_slope_ice_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) z_slope_ice
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/prog_ccn_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) prog_ccn
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/c_cracw_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) c_cracw
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/alin_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) alin
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/clin_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) clin
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/preciprad_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) preciprad
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/cld_min_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) cld_min
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/use_ppm_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) use_ppm
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/mono_prof_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) mono_prof
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/do_sedi_heat_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_sedi_heat
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/sedi_transport_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) sedi_transport
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/do_sedi_w_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) do_sedi_w
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/dt_fr_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dt_fr
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/de_ice_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) de_ice
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/icloud_f_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) icloud_f
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/irain_f_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) irain_f
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/test/mp_print_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) mp_print
        close(fileID)

        ! Note : do_qa is a module variable in gfdl_cloud_microphys.F90
        ! do_qa = .false.

    end subroutine

    subroutine GFDL_1M_RUN_driver(IM, JM, LM, dirName, rank_str)

        integer, intent(in) :: IM,JM,LM

        character*32, intent(IN) :: dirName, rank_str

        ! Local
        real, allocatable, dimension(:,:,:)  :: U0, V0
        real, allocatable, dimension(:,:,:)  :: PLEmb, ZLE0
        real, allocatable, dimension(:,:,:)  :: PLmb,  ZL0
        real, allocatable, dimension(:,:,:)  :: DZ, DZET, DP, MASS, iMASS
        real, allocatable, dimension(:,:,:)  :: DQST3, QST3
        real, allocatable, dimension(:,:,:)  :: DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
                                                DQSDTmic, DQGDTmic, DQADTmic, &
                                                DUDTmic,  DVDTmic,  DTDTmic
        real, allocatable, dimension(:,:,:)  :: TMP3D
        real, allocatable, dimension(:,:)    :: frland2D
        real, allocatable, dimension(:,:)    :: TMP2D
        integer, allocatable, dimension(:,:) :: KLCL

        real                                 :: DT_MOIST

        ! Local variables
        real    :: facEIS
        real    :: minrhcrit, ALPHA, RHCRIT
        ! integer :: IM,JM,LM
        integer :: I, J, L

        integer :: fileID

        logical :: file_exists

        real, allocatable, dimension(:,:,:) :: PLEmb_ref, PLmb_ref, ZLE0_ref, ZL0_ref, &
                                               DZET_ref, DQST3_ref, DP_ref, MASS_ref, U0_ref, V0_ref, &
                                               QST3_ref, TMP3D_ref

        real, allocatable, dimension(:,:,:) :: T_ref, Q_ref, QLLS_ref, QILS_ref, CLLS_ref, QLCN_ref, QICN_ref, CLCN_ref, PDF_A_ref, &
                                               WTHV2_ref, WQL_ref, PDFITERS_ref, EVAPC_ref, SUBLC_ref
        

        real, allocatable, dimension(:,:) :: EIS_ref, LTS_ref
        
        integer, allocatable, dimension(:,:) :: KLCL_ref

        print*, 'Running GFDL_1M_RUN_driver...'

        DT_MOIST = DT_R8

        ! Allocatables
        ! Edge variables 
        ALLOCATE ( ZLE0 (IM,JM,0:LM) )
        ALLOCATE ( PLEmb(IM,JM,0:LM) )
        ! Layer variables
        ALLOCATE ( U0   (IM,JM,LM  ) )
        ALLOCATE ( V0   (IM,JM,LM  ) )
        ALLOCATE ( ZL0  (IM,JM,LM  ) )
        ALLOCATE ( PLmb (IM,JM,LM  ) )
        ALLOCATE ( DZET (IM,JM,LM  ) )
        ALLOCATE ( DZ   (IM,JM,LM  ) )
        ALLOCATE ( DP   (IM,JM,LM  ) )
        ALLOCATE ( MASS (IM,JM,LM  ) )
        ALLOCATE ( iMASS(IM,JM,LM  ) )
        ALLOCATE ( DQST3(IM,JM,LM  ) )
        ALLOCATE (  QST3(IM,JM,LM  ) )
        ALLOCATE ( TMP3D(IM,JM,LM  ) )
        ! Local tendencies
        ALLOCATE ( DQVDTmic(IM,JM,LM  ) )
        ALLOCATE ( DQLDTmic(IM,JM,LM  ) )
        ALLOCATE ( DQIDTmic(IM,JM,LM  ) )
        ALLOCATE ( DQRDTmic(IM,JM,LM  ) )
        ALLOCATE ( DQSDTmic(IM,JM,LM  ) )
        ALLOCATE ( DQGDTmic(IM,JM,LM  ) )
        ALLOCATE ( DQADTmic(IM,JM,LM  ) )
        ALLOCATE (  DUDTmic(IM,JM,LM  ) )
        ALLOCATE (  DVDTmic(IM,JM,LM  ) )
        ALLOCATE (  DTDTmic(IM,JM,LM  ) )
        ! 2D Variables
        ALLOCATE ( frland2D     (IM,JM) ) 
        ALLOCATE ( KLCL         (IM,JM) )
        ALLOCATE ( TMP2D        (IM,JM) )


        ALLOCATE ( PLEmb_ref(IM,JM,0:LM) )
        ALLOCATE ( PLmb_ref (IM,JM,LM  ) )
        ALLOCATE ( ZLE0_ref (IM,JM,0:LM) )
        ALLOCATE ( ZL0_ref  (IM,JM,LM  ) )
        ALLOCATE ( DZET_ref (IM,JM,LM  ) )
        ALLOCATE ( DQST3_ref(IM,JM,LM  ) )
        ALLOCATE ( DP_ref   (IM,JM,LM  ) )
        ALLOCATE ( MASS_ref (IM,JM,LM  ) )
        ALLOCATE ( U0_ref   (IM,JM,LM  ) )
        ALLOCATE ( V0_ref   (IM,JM,LM  ) )
        ALLOCATE (  QST3_ref(IM,JM,LM  ) )

        ALLOCATE ( KLCL_ref         (IM,JM) )
        ALLOCATE ( TMP3D_ref (IM,JM,LM  ) )
        allocate(LTS_ref(IM, JM))
        allocate(EIS_ref(IM, JM))

        allocate(T_ref(IM, JM, LM))
        allocate(Q_ref(IM, JM, LM))
        allocate(QLLS_ref(IM, JM, LM))
        allocate(QILS_ref(IM, JM, LM))
        allocate(CLLS_ref(IM, JM, LM))
        allocate(QLCN_ref(IM, JM, LM))
        allocate(QICN_ref(IM, JM, LM))
        allocate(CLCN_ref(IM, JM, LM))
        allocate(PDF_A_ref(IM, JM, LM))
        allocate(WTHV2_ref(IM, JM, LM))
        allocate(WQL_ref(IM, JM, LM))
        allocate(PDFITERS_ref(IM, JM, LM))
        allocate(EVAPC_ref(IM, JM, LM))
        allocate(SUBLC_ref(IM, JM, LM))

        ! Derived States
        PLEmb    =  PLE*.01
        PLmb     = 0.5*(PLEmb(:,:,0:LM-1) + PLEmb(:,:,1:LM))
        DO L=0,LM
            ZLE0(:,:,L)= ZLE(:,:,L) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
        END DO
        ZL0      = 0.5*(ZLE0(:,:,0:LM-1) + ZLE0(:,:,1:LM) ) ! Layer Height (m) above the surface
        DZET     =     (ZLE0(:,:,0:LM-1) - ZLE0(:,:,1:LM) ) ! Layer thickness (m)
        ! Note : DQST3 value isn't used anywhere.  QST3 is used below
        DQST3    = GEOS_DQSAT(T, PLmb, QSAT=QST3)
        DP       = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )
        MASS     = DP/MAPL_GRAV
        ! Note : iMASS value isn't used anywhere
        iMASS    = 1.0/MASS
        U0       = U
        V0       = V

        ! ** Check Variables
            
        open(newunit=fileID, file=trim(dirName) // '/PLEmb_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) PLEmb_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PLmb_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) PLmb_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/ZLE0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ZLE0_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/ZL0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ZL0_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DZET_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DZET_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQST3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DQST3_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DP_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/MASS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) MASS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/U0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) U0_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/V0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) V0_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QST3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QST3_ref
        close(fileID)

        print*,'Compare sum(diff(PLEmb)) = ',sum(PLEmb_ref - PLEmb)
        print*,'Compare sum(PLEmb) = ',sum(PLEmb)
        print*,'Compare sum(PLEmb_ref) = ',sum(PLEmb_ref)
        print*,'***'
        print*,'Compare sum(diff(PLmb)) = ',sum(PLmb_ref - PLmb)
        print*,'Compare sum(PLmb) = ',sum(PLmb)
        print*,'Compare sum(PLmb_ref) = ',sum(PLmb_ref)
        print*,'***'
        print*,'Compare sum(diff(ZLE0)) = ',sum(ZLE0_ref - ZLE0)
        print*,'Compare sum(ZLE0) = ',sum(ZLE0)
        print*,'Compare sum(ZLE0_ref) = ',sum(ZLE0_ref)
        print*,'***'
        print*,'Compare sum(diff(ZL0)) = ',sum(ZL0_ref - ZL0)
        print*,'Compare sum(ZL0) = ',sum(ZL0)
        print*,'Compare sum(ZL0_ref) = ',sum(ZL0_ref)
        print*,'***'
        print*,'Compare sum(diff(DZET)) = ',sum(DZET_ref - DZET)
        print*,'Compare sum(DZET) = ',sum(DZET)
        print*,'Compare sum(DZET_ref) = ',sum(DZET_ref)
        print*,'***'
        ! print*,'Compare sum(diff(DQST3)) = ',sum(DQST3_ref - DQST3)
        ! print*,'Compare sum(DQST3) = ',sum(DQST3)
        ! print*,'Compare sum(DQST3_ref) = ',sum(DQST3_ref)
        ! print*,'***'
        print*,'Compare sum(diff(DP)) = ',sum(DP_ref - DP)
        print*,'Compare sum(DP) = ',sum(DP)
        print*,'Compare sum(DP_ref) = ',sum(DP_ref)
        print*,'***'
        print*,'Compare sum(diff(MASS)) = ',sum(MASS_ref - MASS)
        print*,'Compare sum(MASS) = ',sum(MASS)
        print*,'Compare sum(MASS_ref) = ',sum(MASS_ref)
        print*,'***'
        print*,'Compare sum(diff(U0)) = ',sum(U0_ref - U0)
        print*,'Compare sum(U0) = ',sum(U0)
        print*,'Compare sum(U0_ref) = ',sum(U0_ref)
        print*,'***'
        print*,'Compare sum(diff(V0)) = ',sum(V0_ref - V0)
        print*,'Compare sum(V0) = ',sum(V0)
        print*,'Compare sum(V0_ref) = ',sum(V0_ref)
        print*,'***'
        print*,'Compare sum(diff(QST3)) = ',sum(abs(QST3_ref - QST3))
        print*,'Compare sum(QST3) = ',sum(QST3)
        print*,'Compare sum(QST3_ref) = ',sum(QST3_ref)
        print*,'***'

        call exit(1)

        ! Note : QST3 has a slight variation in standalone calcuation when compared to Discover calcuations
        ! ***

        KLCL = FIND_KLCL( T, Q, PLmb, IM, JM, LM ) 

        if (associated(PTR2D)) then
            do J=1,JM
               do I=1,IM
                    PTR2D(I,J) = ZL0(I,J,KLCL(I,J))
               end do
            end do
        endif
        TMP3D = (100.0*PLmb/MAPL_P00)**(MAPL_KAPPA)
        call FIND_EIS(T/TMP3D, QST3, T, ZL0, PLEmb, KLCL, IM, JM, LM, LTS, EIS)

        open(newunit=fileID, file=trim(dirName) // '/KLCL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) KLCL_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/TMP3D_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) TMP3D_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/LTS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) LTS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/EIS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) EIS_ref
        close(fileID)

        print*,'Compare sum(diff(KLCL)) = ',sum(KLCL_ref - KLCL)
        print*,'Compare sum(KLCL) = ',sum(KLCL)
        print*,'Compare sum(KLCL_ref) = ',sum(KLCL_ref)
        print*,'***'
        print*,'Compare sum(diff(TMP3D)) = ',sum(TMP3D_ref - TMP3D)
        print*,'Compare sum(TMP3D) = ',sum(TMP3D)
        print*,'Compare sum(TMP3D_ref) = ',sum(TMP3D_ref)
        print*,'***'
        print*,'Compare sum(diff(LTS)) = ',sum(LTS_ref - LTS)
        print*,'Compare sum(LTS) = ',sum(LTS)
        print*,'Compare sum(LTS_ref) = ',sum(LTS_ref)
        print*,'***'
        print*,'Compare sum(diff(EIS)) = ',sum(EIS_ref - EIS)
        print*,'Compare sum(EIS) = ',sum(EIS)
        print*,'Compare sum(EIS_ref) = ',sum(EIS_ref)
        print*,'***'


        DUDT_macro=U
        DVDT_macro=V
        DTDT_macro=T
        DQVDT_macro=Q
        DQLDT_macro=QLCN+QLLS
        DQIDT_macro=QICN+QILS
        DQADT_macro=CLCN+CLLS
        DQRDT_macro=QRAIN
        DQSDT_macro=QSNOW
        DQGDT_macro=QGRAUPEL

         ! Include shallow precip condensates if present
        ! call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_PRC3', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) then
          QRAIN = QRAIN + PTR3D*DT_MOIST
        endif
        ! call MAPL_GetPointer(EXPORT, PTR3D,  'SHLW_SNO3', RC=STATUS); VERIFY_(STATUS)
        if (associated(PTR3D)) then 
          QSNOW = QSNOW + PTR3D*DT_MOIST
        endif

        do L=1,LM
            do J=1,JM
                do I=1,IM
                    ! Send the condensates through the pdf after convection
                    facEIS = MAX(0.0,MIN(1.0,EIS(I,J)/10.0))**2
                    ! determine combined minrhcrit in stable/unstable regimes
                    minrhcrit  = (1.0-dw_ocean)*(1.0-facEIS) + (1.0-dw_land)*facEIS
                    if (turnrhcrit <= 0.0) then
                        ! determine the turn pressure using the LCL
                        turnrhcrit  = PLmb(I, J, KLCL(I,J)) - 250.0 ! 250mb above the LCL
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
                            QST3(I,J,L)  )
                        EVAPC(I,J,L) = ( Q(I,J,L) - EVAPC(I,J,L) ) / DT_MOIST
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
                            QST3(I,J,L)  )
                        SUBLC(I,J,L) = ( Q(I,J,L) - SUBLC(I,J,L) ) / DT_MOIST
                    endif
                    ! cleanup clouds
                    call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L) )
                end do ! IM loop
            end do ! JM loop
        end do ! LM loop

        open(newunit=fileID, file=trim(dirName) // '/T_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) T_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/Q_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) Q_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QLLS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) QLLS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QILS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) QILS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CLLS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) CLLS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QLCN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) QLCN_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/QICN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) QICN_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CLCN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) CLCN_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PDF_A_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PDF_A_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/WTHV2_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) WTHV2_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/WQL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) WQL_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PDFITERS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PDFITERS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/EVAPC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) EVAPC_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/SUBLC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) SUBLC_ref
        close(fileID)

        print*,'Compare sum(diff(T)) = ',sum(T_ref - T)
        print*,'Compare sum(T) = ',sum(T)
        print*,'Compare sum(T_ref) = ',sum(T_ref)
        print*,'***'
        print*,'Compare sum(diff(Q)) = ',sum(Q_ref - Q)
        print*,'Compare sum(Q) = ',sum(Q)
        print*,'Compare sum(Q_ref) = ',sum(Q_ref)
        print*,'***'
        print*,'Compare sum(diff(QLLS)) = ',sum(QLLS_ref - QLLS)
        print*,'Compare sum(QLLS) = ',sum(QLLS)
        print*,'Compare sum(QLLS_ref) = ',sum(QLLS_ref)
        print*,'***'
        print*,'Compare sum(diff(QILS)) = ',sum(QILS_ref - QILS)
        print*,'Compare sum(QILS) = ',sum(QILS)
        print*,'Compare sum(QILS_ref) = ',sum(QILS_ref)
        print*,'***'
        print*,'Compare sum(diff(CLLS)) = ',sum(CLLS_ref - CLLS)
        print*,'Compare sum(CLLS) = ',sum(CLLS)
        print*,'Compare sum(CLLS_ref) = ',sum(CLLS_ref)
        print*,'***'
        print*,'Compare sum(diff(QLCN)) = ',sum(QLCN_ref - QLCN)
        print*,'Compare sum(QLCN) = ',sum(QLCN)
        print*,'Compare sum(QLCN_ref) = ',sum(QLCN_ref)
        print*,'***'
        print*,'Compare sum(diff(QICN)) = ',sum(QICN_ref - QICN)
        print*,'Compare sum(QICN) = ',sum(QICN)
        print*,'Compare sum(QICN_ref) = ',sum(QICN_ref)
        print*,'***'
        print*,'Compare sum(diff(CLCN)) = ',sum(CLCN_ref - CLCN)
        print*,'Compare sum(CLCN) = ',sum(CLCN)
        print*,'Compare sum(CLCN_ref) = ',sum(CLCN_ref)
        print*,'***'
        print*,'Compare sum(diff(PDF_A)) = ',sum(PDF_A_ref - PDF_A)
        print*,'Compare sum(PDF_A) = ',sum(PDF_A)
        print*,'Compare sum(PDF_A_ref) = ',sum(PDF_A_ref)
        print*,'***'
        print*,'Compare sum(diff(WTHV2)) = ',sum(WTHV2_ref - WTHV2)
        print*,'Compare sum(WTHV2) = ',sum(WTHV2)
        print*,'Compare sum(WTHV2_ref) = ',sum(WTHV2_ref)
        print*,'***'
        print*,'Compare sum(diff(WQL)) = ',sum(WQL_ref - WQL)
        print*,'Compare sum(WQL) = ',sum(WQL)
        print*,'Compare sum(WQL_ref) = ',sum(WQL_ref)
        print*,'***'
        print*,'Compare sum(diff(PDFITERS)) = ',sum(PDFITERS_ref - PDFITERS)
        print*,'Compare sum(PDFITERS) = ',sum(PDFITERS)
        print*,'Compare sum(PDFITERS_ref) = ',sum(PDFITERS_ref)
        print*,'***'
        print*,'Compare sum(diff(EVAPC)) = ',sum(EVAPC_ref - EVAPC)
        print*,'Compare sum(EVAPC) = ',sum(EVAPC)
        print*,'Compare sum(EVAPC_ref) = ',sum(EVAPC_ref)
        print*,'***'
        print*,'Compare sum(diff(SUBLC)) = ',sum(SUBLC_ref - SUBLC)
        print*,'Compare sum(SUBLC) = ',sum(SUBLC)
        print*,'Compare sum(SUBLC_ref) = ',sum(SUBLC_ref)
        print*,'***'
        
        call exit(1)

        ! Update macrophysics tendencies
        DUDT_macro=( U         - DUDT_macro)/DT_MOIST
        DVDT_macro=( V         - DVDT_macro)/DT_MOIST
        DTDT_macro=( T         - DTDT_macro)/DT_MOIST
        DQVDT_macro=( Q         -DQVDT_macro)/DT_MOIST
        DQLDT_macro=((QLCN+QLLS)-DQLDT_macro)/DT_MOIST
        DQIDT_macro=((QICN+QILS)-DQIDT_macro)/DT_MOIST
        DQADT_macro=((CLCN+CLLS)-DQADT_macro)/DT_MOIST
        DQRDT_macro=( QRAIN     -DQRDT_macro)/DT_MOIST
        DQSDT_macro=( QSNOW     -DQSDT_macro)/DT_MOIST
        DQGDT_macro=( QGRAUPEL  -DQGDT_macro)/DT_MOIST

        ! Zero-out microphysics tendencies
        DQVDT_micro = Q
        DQLDT_micro = QLLS + QLCN
        DQIDT_micro = QILS + QICN
        DQRDT_micro = QRAIN
        DQSDT_micro = QSNOW
        DQGDT_micro = QGRAUPEL
        DQADT_micro = CLLS + CLCN
        DUDT_micro = U
        DVDT_micro = V
        DTDT_micro = T

        ! Delta-Z layer thickness (gfdl expects this to be negative)
        DZ = -1.0*DZET
        ! Zero-out local microphysics tendencies
        DQVDTmic = 0.
        DQLDTmic = 0.
        DQRDTmic = 0.
        DQIDTmic = 0.
        DQSDTmic = 0.
        DQGDTmic = 0.
        DQADTmic = 0.
        DUDTmic = 0.
        DVDTmic = 0.
        DTDTmic = 0.
       ! Zero-out 3D Precipitation Fluxes 
        ! Ice
        PFI_LS = 0.
        ! Liquid
        PFL_LS = 0.
        ! Cloud
        RAD_CF = MIN(CLCN+CLLS,1.0)
        ! Liquid
        RAD_QL = QLCN+QLLS
        ! Ice
        RAD_QI = QICN+QILS
        ! VAPOR
        RAD_QV = Q
        ! RAIN
        RAD_QR = QRAIN
        ! SNOW
        RAD_QS = QSNOW
        ! GRAUPEL
        RAD_QG = QGRAUPEL

        ! Subroutines calls that set up Cloud Microphysics module variables
        call setup_con
        call setupm
        
        call gfdl_cloud_microphys_driver( &
            ! Input water/cloud species and liquid+ice CCN [NACTL+NACTI (#/m^3)]
            RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, (NACTL+NACTI), &
            ! Output tendencies
            DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
            DQSDTmic, DQGDTmic, DQADTmic, DTDTmic, &
            ! Input fields
            T, W, U, V, DUDTmic, DVDTmic, DZ, DP, &
            ! constant inputs
            AREA, DT_MOIST, frland2D, CNV_FRC, SRF_TYPE, EIS, &
            RHCRIT3D, ANV_ICEFALL, LS_ICEFALL, &
            ! Output rain re-evaporation and sublimation
            REV_LS, RSU_LS, & 
            ! Output precipitates
            PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
            ! Output mass flux during sedimentation (Pa kg/kg)
            PFL_LS(:,:,1:LM), PFI_LS(:,:,1:LM), &
            ! constant grid/time information
            LHYDROSTATIC, LPHYS_HYDROSTATIC, &
            1,IM, 1,JM, 1,LM, 1, LM)

        ! Apply tendencies
        T = T + DTDTmic * DT_MOIST
        U = U + DUDTmic * DT_MOIST
        V = V + DVDTmic * DT_MOIST
        
        ! Apply moist/cloud species tendencies
        RAD_QV = RAD_QV + DQVDTmic * DT_MOIST
        RAD_QL = RAD_QL + DQLDTmic * DT_MOIST
        RAD_QR = RAD_QR + DQRDTmic * DT_MOIST
        RAD_QI = RAD_QI + DQIDTmic * DT_MOIST
        RAD_QS = RAD_QS + DQSDTmic * DT_MOIST
        RAD_QG = RAD_QG + DQGDTmic * DT_MOIST
        RAD_CF = MIN(1.0,MAX(0.0,RAD_CF + DQADTmic * DT_MOIST))
        
        ! Redistribute CN/LS CF/QL/QI
        call REDISTRIBUTE_CLOUDS(RAD_CF, RAD_QL, RAD_QI, CLCN, CLLS, QLCN, QLLS, QICN, QILS, RAD_QV, T)

        ! Convert precip diagnostics from mm/day to kg m-2 s-1
        PRCP_RAIN    = MAX(PRCP_RAIN    / 86400.0, 0.0)
        PRCP_SNOW    = MAX(PRCP_SNOW    / 86400.0, 0.0)
        PRCP_ICE     = MAX(PRCP_ICE     / 86400.0, 0.0)
        PRCP_GRAUPEL = MAX(PRCP_GRAUPEL / 86400.0, 0.0)
        
        ! Fill GEOS precip diagnostics
        LS_PRCP = PRCP_RAIN
        LS_SNR  = PRCP_SNOW
        ICE     = PRCP_ICE + PRCP_GRAUPEL
        FRZR    = 0.0
        
        ! Convert precipitation fluxes from (Pa kg/kg) to (kg m-2 s-1)
        PFL_LS = PFL_LS/(MAPL_GRAV*DT_MOIST)
        PFI_LS = PFI_LS/(MAPL_GRAV*DT_MOIST)
        
        ! Redistribute precipitation fluxes for chemistry
        TMP3D =  MIN(1.0,MAX(QLCN/MAX(RAD_QL,1.E-8),0.0))
        PFL_AN(:,:,1:LM) = PFL_LS(:,:,1:LM) * TMP3D
        PFL_LS(:,:,1:LM) = PFL_LS(:,:,1:LM) - PFL_AN(:,:,1:LM)
        TMP3D =  MIN(1.0,MAX(QICN/MAX(RAD_QI,1.E-8),0.0))
        PFI_AN(:,:,1:LM) = PFI_LS(:,:,1:LM) * TMP3D
        PFI_LS(:,:,1:LM) = PFI_LS(:,:,1:LM) - PFI_AN(:,:,1:LM)
        
        ! cleanup suspended precipitation condensates
        call FIX_NEGATIVE_PRECIP(RAD_QR, RAD_QS, RAD_QG)
        
        ! Fill vapor/rain/snow/graupel state
        Q        = RAD_QV
        QRAIN    = RAD_QR
        QSNOW    = RAD_QS
        QGRAUPEL = RAD_QG

        ! Radiation Coupling
        do L = 1, LM
            do J = 1, JM
                do I = 1, IM
                    ! cleanup clouds
                    call FIX_UP_CLOUDS( Q(I,J,L), T(I,J,L), QLLS(I,J,L), QILS(I,J,L), CLLS(I,J,L), QLCN(I,J,L), QICN(I,J,L), CLCN(I,J,L) )
                    ! get radiative properties
                    call RADCOUPLE ( T(I,J,L), PLmb(I,J,L), CLLS(I,J,L), CLCN(I,J,L), &
                        Q(I,J,L), QLLS(I,J,L), QILS(I,J,L), QLCN(I,J,L), QICN(I,J,L), QRAIN(I,J,L), QSNOW(I,J,L), QGRAUPEL(I,J,L), NACTL(I,J,L), NACTI(I,J,L), &
                        RAD_QV(I,J,L), RAD_QL(I,J,L), RAD_QI(I,J,L), RAD_QR(I,J,L), RAD_QS(I,J,L), RAD_QG(I,J,L), RAD_CF(I,J,L), &
                        CLDREFFL(I,J,L), CLDREFFI(I,J,L), &
                        FAC_RL, MIN_RL, MAX_RL, FAC_RI, MIN_RI, MAX_RI)
                    if (do_qa) RHX(I,J,L) = Q(I,J,L)/GEOS_QSAT( T(I,J,L), PLmb(I,J,L) )
                enddo
            enddo
        enddo

        call FILLQ2ZERO(RAD_QV, MASS, TMP2D)
        call FILLQ2ZERO(RAD_QL, MASS, TMP2D)
        call FILLQ2ZERO(RAD_QI, MASS, TMP2D)
        call FILLQ2ZERO(RAD_QR, MASS, TMP2D)
        call FILLQ2ZERO(RAD_QS, MASS, TMP2D)
        call FILLQ2ZERO(RAD_QG, MASS, TMP2D)
        call FILLQ2ZERO(RAD_CF, MASS, TMP2D)
        where (RAD_QI .le. 0.0)
        CLDREFFI = MAPL_UNDEF
        end where
        where (RAD_QL .le. 0.0)
        CLDREFFL = MAPL_UNDEF
        end where

        ! Update microphysics tendencies
        DQVDT_micro = ( Q          - DQVDT_micro) / DT_MOIST
        DQLDT_micro = ((QLLS+QLCN) - DQLDT_micro) / DT_MOIST
        DQIDT_micro = ((QILS+QICN) - DQIDT_micro) / DT_MOIST
        DQADT_micro = ((CLLS+CLCN) - DQADT_micro) / DT_MOIST
        DQRDT_micro = ( QRAIN      - DQRDT_micro) / DT_MOIST
        DQSDT_micro = ( QSNOW      - DQSDT_micro) / DT_MOIST
        DQGDT_micro = ( QGRAUPEL   - DQGDT_micro) / DT_MOIST
         DUDT_micro = ( U          -  DUDT_micro) / DT_MOIST
         DVDT_micro = ( V          -  DVDT_micro) / DT_MOIST
         DTDT_micro = ( T          -  DTDT_micro) / DT_MOIST

        if(associated(PTR3D)) PTR3D = DQRDT_macro + DQRDT_micro

        if(associated(PTR3D)) then
            call dissipative_ke_heating(IM,JM,LM, MASS,U0,V0, &
                                        DUDT_macro+DUDT_micro,&
                                        DVDT_macro+DVDT_micro,PTR3D)
        endif

        if (associated(PTR3D) .OR. &
            associated(DBZ_MAX) .OR. associated(DBZ_1KM) .OR. associated(DBZ_TOP) .OR. associated(DBZ_M10C)) then

            call CALCDBZ(TMP3D,100*PLmb,T,Q,QRAIN,QSNOW,QGRAUPEL,IM,JM,LM,1,0,1)
            if (associated(PTR3D)) PTR3D = TMP3D

            if (associated(DBZ_MAX)) then
                DBZ_MAX=-9999.0
                DO L=1,LM ; DO J=1,JM ; DO I=1,IM
                    DBZ_MAX(I,J) = MAX(DBZ_MAX(I,J),TMP3D(I,J,L))
                END DO ; END DO ; END DO
            endif

            if (associated(DBZ_1KM)) then  
                call cs_interpolator(1, IM, 1, JM, LM, TMP3D, 1000., ZLE0, DBZ_1KM, -20.)
            endif
                   
            if (associated(DBZ_TOP)) then
                DBZ_TOP=MAPL_UNDEF  
                DO J=1,JM ; DO I=1,IM
                    DO L=LM,1,-1
                        if (ZLE0(i,j,l) >= 25000.) continue
                        if (TMP3D(i,j,l) >= 18.5 ) then
                            DBZ_TOP(I,J) = ZLE0(I,J,L)
                            exit
                        endif
                    END DO
                END DO ; END DO
            endif 
                   
            if (associated(DBZ_M10C)) then
                DBZ_M10C=MAPL_UNDEF  
                DO J=1,JM ; DO I=1,IM
                    DO L=LM,1,-1
                        if (ZLE0(i,j,l) >= 25000.) continue
                        if (T(i,j,l) <= MAPL_TICE-10.0) then
                            DBZ_M10C(I,J) = TMP3D(I,J,L)
                            exit
                        endif
                    END DO
                END DO ; END DO
            endif        

        endif

    end subroutine

    subroutine compare_results(IM, JM, LM, dirName, rank_str)

        integer, intent(in)       :: IM, JM, LM
        character*100, intent(IN) :: dirName, rank_str
        integer                   :: fileID

        logical                   :: file_exists

        real, dimension(:,:  ), allocatable :: PRCP_RAIN_ref, PRCP_SNOW_ref, PRCP_ICE_ref, PRCP_GRAUPEL_ref
        real, dimension(:,:  ), allocatable :: LS_PRCP_ref, LS_SNR_ref, ICE_ref, FRZR_ref, CNV_FRC_ref, SRF_TYPE_ref
        real, dimension(:,:,:), allocatable :: DQVDT_macro_ref, DQIDT_macro_ref, DQLDT_macro_ref, DQADT_macro_ref, DQRDT_macro_ref, DQSDT_macro_ref, DQGDT_macro_ref
        real, dimension(:,:,:), allocatable :: DUDT_macro_ref,  DVDT_macro_ref,  DTDT_macro_ref
        real, dimension(:,:,:), allocatable :: DQVDT_micro_ref, DQIDT_micro_ref, DQLDT_micro_ref, DQADT_micro_ref, DQRDT_micro_ref, DQSDT_micro_ref, DQGDT_micro_ref
        real, dimension(:,:,:), allocatable ::  DUDT_micro_ref,  DVDT_micro_ref,  DTDT_micro_ref
        real, dimension(:,:,:), allocatable :: RAD_CF_ref, RAD_QV_ref, RAD_QL_ref, RAD_QI_ref, RAD_QR_ref, RAD_QS_ref, RAD_QG_ref
        real, dimension(:,:,:), allocatable :: CLDREFFL_ref, CLDREFFI_ref
        real, dimension(:,:,:), allocatable :: EVAPC_ref, SUBLC_ref
        real, dimension(:,:,:), allocatable :: RHX_ref, REV_LS_ref, RSU_LS_ref
        real, dimension(:,:,:), allocatable :: PFL_LS_ref, PFL_AN_ref
        real, dimension(:,:,:), allocatable :: PFI_LS_ref, PFI_AN_ref
        real, dimension(:,:,:), allocatable :: PDF_A_ref, PDFITERS_ref
        real, dimension(:,:,:), allocatable :: RHCRIT3D_ref, WTHV2_ref, WQL_ref
        real, dimension(:,:), allocatable   :: EIS_ref, LTS_ref
        real, dimension(:,:), allocatable   :: DBZ_MAX_ref, DBZ_1KM_ref, DBZ_TOP_ref, DBZ_M10C_ref
        real, dimension(:,:,:), allocatable :: PTR3D_ref
        real, dimension(:,:  ), allocatable :: PTR2D_ref

        print*, 'Running compare_results...'

        allocate(RAD_CF_ref(IM, JM, LM))
        allocate(RAD_QV_ref(IM, JM, LM))
        allocate(RAD_QL_ref(IM, JM, LM))
        allocate(RAD_QI_ref(IM, JM, LM))
        allocate(RAD_QR_ref(IM, JM, LM))
        allocate(RAD_QS_ref(IM, JM, LM))
        allocate(RAD_QG_ref(IM, JM, LM))
        allocate(CLDREFFL_ref(IM, JM, LM))
        allocate(CLDREFFI_ref(IM, JM, LM))
        allocate(CNV_FRC_ref(IM, JM))
        allocate(SRF_TYPE_ref(IM, JM))
        allocate(EVAPC_ref(IM, JM, LM))
        allocate(SUBLC_ref(IM, JM, LM))
        allocate(PRCP_RAIN_ref(IM, JM))
        allocate(PRCP_SNOW_ref(IM, JM))
        allocate(PRCP_ICE_ref(IM, JM))
        allocate(PRCP_GRAUPEL_ref(IM, JM))
        allocate(LS_PRCP_ref(IM, JM))
        allocate(LS_SNR_ref(IM, JM))
        allocate(ICE_ref(IM, JM))
        allocate(FRZR_ref(IM, JM))
        allocate(RHX_ref(IM, JM, LM))
        allocate(REV_LS_ref(IM, JM, LM))
        allocate(RSU_LS_ref(IM, JM, LM))
        allocate(PFL_AN_ref(IM, JM, 0:LM))
        allocate(PFL_LS_ref(IM, JM, 0:LM))
        allocate(PFI_AN_ref(IM, JM, 0:LM))
        allocate(PFI_LS_ref(IM, JM, 0:LM))
        allocate(PDF_A_ref(IM, JM, LM))
        allocate(WTHV2_ref(IM, JM, LM))
        allocate(WQL_ref(IM, JM, LM))
        allocate(PDFITERS_ref(IM, JM, LM))
        allocate(LTS_ref(IM, JM))
        allocate(EIS_ref(IM, JM))
        allocate(DQVDT_macro_ref(IM, JM, LM))
        allocate(DQIDT_macro_ref(IM, JM, LM))
        allocate(DQLDT_macro_ref(IM, JM, LM))
        allocate(DQADT_macro_ref(IM, JM, LM))
        allocate(DQRDT_macro_ref(IM, JM, LM))
        allocate(DQSDT_macro_ref(IM, JM, LM))
        allocate(DQGDT_macro_ref(IM, JM, LM))
        allocate(DUDT_macro_ref(IM, JM, LM))
        allocate(DVDT_macro_ref(IM, JM, LM))
        allocate(DTDT_macro_ref(IM, JM, LM))
        allocate(RHCRIT3D_ref(IM, JM, LM))
        allocate(DQVDT_micro_ref(IM, JM, LM))
        allocate(DQIDT_micro_ref(IM, JM, LM))
        allocate(DQLDT_micro_ref(IM, JM, LM))
        allocate(DQADT_micro_ref(IM, JM, LM))
        allocate(DQRDT_micro_ref(IM, JM, LM))
        allocate(DQSDT_micro_ref(IM, JM, LM))
        allocate(DQGDT_micro_ref(IM, JM, LM))
        allocate(DUDT_micro_ref(IM, JM, LM))
        allocate(DVDT_micro_ref(IM, JM, LM))
        allocate(DTDT_micro_ref(IM, JM, LM))

        allocate(DBZ_M10C_ref(IM, JM))

        inquire(file=trim(dirName) // '/DBZ_MAX_' // trim(rank_str) // '.out', exist=file_exists)
        if(file_exists) then
            allocate(DBZ_MAX_ref(IM, JM))
            open(newunit=fileID, file=trim(dirName) // '/DBZ_MAX_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
            read(fileID) DBZ_MAX_ref
            close(fileID)
        endif

        inquire(file=trim(dirName) // '/DBZ_1KM_' // trim(rank_str) // '.out', exist=file_exists)
        if(file_exists) then
            allocate(DBZ_1KM_ref(IM, JM))
            open(newunit=fileID, file=trim(dirName) // '/DBZ_1KM_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
            read(fileID) DBZ_1KM_ref
            close(fileID)
        endif

        inquire(file=trim(dirName) // '/DBZ_TOP_' // trim(rank_str) // '.out', exist=file_exists)
        if(file_exists) then
            allocate(DBZ_TOP_ref(IM, JM))
            open(newunit=fileID, file=trim(dirName) // '/DBZ_TOP_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
            read(fileID) DBZ_TOP_ref
            close(fileID)
        endif

        inquire(file=trim(dirName) // '/DBZ_M10C_' // trim(rank_str) // '.out', exist=file_exists)
        if(file_exists) then
            allocate(DBZ_M10C_ref(IM, JM))
            open(newunit=fileID, file=trim(dirName) // '/DBZ_M10C_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
            read(fileID) DBZ_M10C_ref
            close(fileID)
        endif

        open(newunit=fileID, file=trim(dirName) // '/RAD_CF_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_CF_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_QV_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QV_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_QL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QL_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_QI_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QI_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_QR_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QR_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_QS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RAD_QG_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QG_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CLDREFFL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) CLDREFFL_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CLDREFFI_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) CLDREFFI_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/CNV_FRC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) CNV_FRC_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/SRF_TYPE_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) SRF_TYPE_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/EVAPC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) EVAPC_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/SUBLC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) SUBLC_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PRCP_RAIN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PRCP_RAIN_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PRCP_SNOW_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PRCP_SNOW_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PRCP_ICE_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PRCP_ICE_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PRCP_GRAUPEL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PRCP_GRAUPEL_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/LS_PRCP_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) LS_PRCP_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/LS_SNR_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) LS_SNR_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/ICE_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) ICE_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/FRZR_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) FRZR_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RHX_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RHX_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/REV_LS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) REV_LS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RSU_LS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RSU_LS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PFL_AN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PFL_AN_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PFL_LS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PFL_LS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PFI_AN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PFI_AN_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PFI_LS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PFI_LS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PDF_A_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PDF_A_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/WTHV2_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) WTHV2_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/WQL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) WQL_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/PDFITERS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PDFITERS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/LTS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) LTS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/EIS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) EIS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQVDT_macro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQVDT_macro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQIDT_macro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQIDT_macro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQLDT_macro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQLDT_macro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQADT_macro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQADT_macro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQRDT_macro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQRDT_macro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQSDT_macro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQSDT_macro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQGDT_macro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQGDT_macro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DUDT_macro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DUDT_macro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DVDT_macro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DVDT_macro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DTDT_macro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DTDT_macro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RHCRIT3D_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RHCRIT3D_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQVDT_micro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQVDT_micro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQIDT_micro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQIDT_micro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQLDT_micro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQLDT_micro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQADT_micro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQADT_micro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQRDT_micro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQRDT_micro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQSDT_micro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQSDT_micro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DQGDT_micro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQGDT_micro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DUDT_micro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DUDT_micro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DVDT_micro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DVDT_micro_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/DTDT_micro_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DTDT_micro_ref
        close(fileID)

        print*,'Compare sum(diff(RAD_CF)) = ',sum(RAD_CF_ref - RAD_CF)
        print*,'Compare sum(RAD_CF) = ',sum(RAD_CF)
        print*,'Compare sum(RAD_CF_ref) = ',sum(RAD_CF_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QV)) = ',sum(RAD_QV_ref - RAD_QV)
        print*,'Compare sum(RAD_QV) = ',sum(RAD_QV)
        print*,'Compare sum(RAD_QV_ref) = ',sum(RAD_QV_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QL)) = ',sum(RAD_QL_ref - RAD_QL)
        print*,'Compare sum(RAD_QL) = ',sum(RAD_QL)
        print*,'Compare sum(RAD_QL_ref) = ',sum(RAD_QL_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QI)) = ',sum(RAD_QI_ref - RAD_QI)
        print*,'Compare sum(RAD_QI) = ',sum(RAD_QI)
        print*,'Compare sum(RAD_QI_ref) = ',sum(RAD_QI_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QR)) = ',sum(RAD_QR_ref - RAD_QR)
        print*,'Compare sum(RAD_QR) = ',sum(RAD_QR)
        print*,'Compare sum(RAD_QR_ref) = ',sum(RAD_QR_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QS)) = ',sum(RAD_QS_ref - RAD_QS)
        print*,'Compare sum(RAD_QS) = ',sum(RAD_QS)
        print*,'Compare sum(RAD_QS_ref) = ',sum(RAD_QS_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QG)) = ',sum(RAD_QG_ref - RAD_QG)
        print*,'Compare sum(RAD_QG) = ',sum(RAD_QG)
        print*,'Compare sum(RAD_QG_ref) = ',sum(RAD_QG_ref)
        print*,'***'
        print*,'Compare sum(diff(CLDREFFL)) = ',sum(CLDREFFL_ref - CLDREFFL)
        print*,'Compare sum(CLDREFFL) = ',sum(CLDREFFL)
        print*,'Compare sum(CLDREFFL_ref) = ',sum(CLDREFFL_ref)
        print*,'***'
        print*,'Compare sum(diff(CLDREFFI)) = ',sum(CLDREFFI_ref - CLDREFFI)
        print*,'Compare sum(CLDREFFI) = ',sum(CLDREFFI)
        print*,'Compare sum(CLDREFFI_ref) = ',sum(CLDREFFI_ref)
        print*,'***'
        print*,'Compare sum(diff(CNV_FRC)) = ',sum(CNV_FRC_ref - CNV_FRC)
        print*,'Compare sum(CNV_FRC) = ',sum(CNV_FRC)
        print*,'Compare sum(CNV_FRC_ref) = ',sum(CNV_FRC_ref)
        print*,'***'
        print*,'Compare sum(diff(SRF_TYPE)) = ',sum(SRF_TYPE_ref - SRF_TYPE)
        print*,'Compare sum(SRF_TYPE) = ',sum(SRF_TYPE)
        print*,'Compare sum(SRF_TYPE_ref) = ',sum(SRF_TYPE_ref)
        print*,'***'
        print*,'Compare sum(diff(EVAPC)) = ',sum(EVAPC_ref - EVAPC)
        print*,'Compare sum(EVAPC) = ',sum(EVAPC)
        print*,'Compare sum(EVAPC_ref) = ',sum(EVAPC_ref)
        print*,'***'
        print*,'Compare sum(diff(SUBLC)) = ',sum(SUBLC_ref - SUBLC)
        print*,'Compare sum(SUBLC) = ',sum(SUBLC)
        print*,'Compare sum(SUBLC_ref) = ',sum(SUBLC_ref)
        print*,'***'
        print*,'Compare sum(diff(PRCP_RAIN)) = ',sum(PRCP_RAIN_ref - PRCP_RAIN)
        print*,'Compare sum(PRCP_RAIN) = ',sum(PRCP_RAIN)
        print*,'Compare sum(PRCP_RAIN_ref) = ',sum(PRCP_RAIN_ref)
        print*,'***'
        print*,'Compare sum(diff(PRCP_SNOW)) = ',sum(PRCP_SNOW_ref - PRCP_SNOW)
        print*,'Compare sum(PRCP_SNOW) = ',sum(PRCP_SNOW)
        print*,'Compare sum(PRCP_SNOW_ref) = ',sum(PRCP_SNOW_ref)
        print*,'***'
        print*,'Compare sum(diff(PRCP_ICE)) = ',sum(PRCP_ICE_ref - PRCP_ICE)
        print*,'Compare sum(PRCP_ICE) = ',sum(PRCP_ICE)
        print*,'Compare sum(PRCP_ICE_ref) = ',sum(PRCP_ICE_ref)
        print*,'***'
        print*,'Compare sum(diff(PRCP_GRAUPEL)) = ',sum(PRCP_GRAUPEL_ref - PRCP_GRAUPEL)
        print*,'Compare sum(PRCP_GRAUPEL) = ',sum(PRCP_GRAUPEL)
        print*,'Compare sum(PRCP_GRAUPEL_ref) = ',sum(PRCP_GRAUPEL_ref)
        print*,'***'
        print*,'Compare sum(diff(LS_PRCP)) = ',sum(LS_PRCP_ref - LS_PRCP)
        print*,'Compare sum(LS_PRCP) = ',sum(LS_PRCP)
        print*,'Compare sum(LS_PRCP_ref) = ',sum(LS_PRCP_ref)
        print*,'***'
        print*,'Compare sum(diff(LS_SNR)) = ',sum(LS_SNR_ref - LS_SNR)
        print*,'Compare sum(LS_SNR) = ',sum(LS_SNR)
        print*,'Compare sum(LS_SNR_ref) = ',sum(LS_SNR_ref)
        print*,'***'
        print*,'Compare sum(diff(ICE)) = ',sum(ICE_ref - ICE)
        print*,'Compare sum(ICE) = ',sum(ICE)
        print*,'Compare sum(ICE_ref) = ',sum(ICE_ref)
        print*,'***'
        print*,'Compare sum(diff(FRZR)) = ',sum(FRZR_ref - FRZR)
        print*,'Compare sum(FRZR) = ',sum(FRZR)
        print*,'Compare sum(FRZR_ref) = ',sum(FRZR_ref)
        print*,'***'
        print*,'Compare sum(diff(RHX)) = ',sum(RHX_ref - RHX)
        print*,'Compare sum(RHX) = ',sum(RHX)
        print*,'Compare sum(RHX_ref) = ',sum(RHX_ref)
        print*,'***'
        print*,'Compare sum(diff(REV_LS)) = ',sum(REV_LS_ref - REV_LS)
        print*,'Compare sum(REV_LS) = ',sum(REV_LS)
        print*,'Compare sum(REV_LS_ref) = ',sum(REV_LS_ref)
        print*,'***'
        print*,'Compare sum(diff(RSU_LS)) = ',sum(RSU_LS_ref - RSU_LS)
        print*,'Compare sum(RSU_LS) = ',sum(RSU_LS)
        print*,'Compare sum(RSU_LS_ref) = ',sum(RSU_LS_ref)
        print*,'***'
        print*,'Compare sum(diff(PFL_AN)) = ',sum(PFL_AN_ref - PFL_AN)
        print*,'Compare sum(PFL_AN) = ',sum(PFL_AN)
        print*,'Compare sum(PFL_AN_ref) = ',sum(PFL_AN_ref)
        print*,'***'
        print*,'Compare sum(diff(PFL_LS)) = ',sum(PFL_LS_ref - PFL_LS)
        print*,'Compare sum(PFL_LS) = ',sum(PFL_LS)
        print*,'Compare sum(PFL_LS_ref) = ',sum(PFL_LS_ref)
        print*,'***'
        print*,'Compare sum(diff(PFI_AN)) = ',sum(PFI_AN_ref - PFI_AN)
        print*,'Compare sum(PFI_AN) = ',sum(PFI_AN)
        print*,'Compare sum(PFI_AN_ref) = ',sum(PFI_AN_ref)
        print*,'***'
        print*,'Compare sum(diff(PFI_LS)) = ',sum(PFI_LS_ref - PFI_LS)
        print*,'Compare sum(PFI_LS) = ',sum(PFI_LS)
        print*,'Compare sum(PFI_LS_ref) = ',sum(PFI_LS_ref)
        print*,'***'
        print*,'Compare sum(diff(PDF_A)) = ',sum(PDF_A_ref - PDF_A)
        print*,'Compare sum(PDF_A) = ',sum(PDF_A)
        print*,'Compare sum(PDF_A_ref) = ',sum(PDF_A_ref)
        print*,'***'
        print*,'Compare sum(diff(WTHV2)) = ',sum(WTHV2_ref - WTHV2)
        print*,'Compare sum(WTHV2) = ',sum(WTHV2)
        print*,'Compare sum(WTHV2_ref) = ',sum(WTHV2_ref)
        print*,'***'
        print*,'Compare sum(diff(WQL)) = ',sum(WQL_ref - WQL)
        print*,'Compare sum(WQL) = ',sum(WQL)
        print*,'Compare sum(WQL_ref) = ',sum(WQL_ref)
        print*,'***'
        print*,'Compare sum(diff(PDFITERS)) = ',sum(PDFITERS_ref - PDFITERS)
        print*,'Compare sum(PDFITERS) = ',sum(PDFITERS)
        print*,'Compare sum(PDFITERS_ref) = ',sum(PDFITERS_ref)
        print*,'***'
        print*,'Compare sum(diff(LTS)) = ',sum(LTS_ref - LTS)
        print*,'Compare sum(LTS) = ',sum(LTS)
        print*,'Compare sum(LTS_ref) = ',sum(LTS_ref)
        print*,'***'
        print*,'Compare sum(diff(EIS)) = ',sum(EIS_ref - EIS)
        print*,'Compare sum(EIS) = ',sum(EIS)
        print*,'Compare sum(EIS_ref) = ',sum(EIS_ref)
        print*,'***'
        print*,'Compare sum(diff(DQVDT_macro)) = ',sum(DQVDT_macro_ref - DQVDT_macro)
        print*,'Compare sum(DQVDT_macro) = ',sum(DQVDT_macro)
        print*,'Compare sum(DQVDT_macro_ref) = ',sum(DQVDT_macro_ref)
        print*,'***'
        print*,'Compare sum(diff(DQIDT_macro)) = ',sum(DQIDT_macro_ref - DQIDT_macro)
        print*,'Compare sum(DQIDT_macro) = ',sum(DQIDT_macro)
        print*,'Compare sum(DQIDT_macro_ref) = ',sum(DQIDT_macro_ref)
        print*,'***'
        print*,'Compare sum(diff(DQLDT_macro)) = ',sum(DQLDT_macro_ref - DQLDT_macro)
        print*,'Compare sum(DQLDT_macro) = ',sum(DQLDT_macro)
        print*,'Compare sum(DQLDT_macro_ref) = ',sum(DQLDT_macro_ref)
        print*,'***'
        print*,'Compare sum(diff(DQADT_macro)) = ',sum(DQADT_macro_ref - DQADT_macro)
        print*,'Compare sum(DQADT_macro) = ',sum(DQADT_macro)
        print*,'Compare sum(DQADT_macro_ref) = ',sum(DQADT_macro_ref)
        print*,'***'
        print*,'Compare sum(diff(DQRDT_macro)) = ',sum(DQRDT_macro_ref - DQRDT_macro)
        print*,'Compare sum(DQRDT_macro) = ',sum(DQRDT_macro)
        print*,'Compare sum(DQRDT_macro_ref) = ',sum(DQRDT_macro_ref)
        print*,'***'
        print*,'Compare sum(diff(DQSDT_macro)) = ',sum(DQSDT_macro_ref - DQSDT_macro)
        print*,'Compare sum(DQSDT_macro) = ',sum(DQSDT_macro)
        print*,'Compare sum(DQSDT_macro_ref) = ',sum(DQSDT_macro_ref)
        print*,'***'
        print*,'Compare sum(diff(DQGDT_macro)) = ',sum(DQGDT_macro_ref - DQGDT_macro)
        print*,'Compare sum(DQGDT_macro) = ',sum(DQGDT_macro)
        print*,'Compare sum(DQGDT_macro_ref) = ',sum(DQGDT_macro_ref)
        print*,'***'
        print*,'Compare sum(diff(DUDT_macro)) = ',sum(DUDT_macro_ref - DUDT_macro)
        print*,'Compare sum(DUDT_macro) = ',sum(DUDT_macro)
        print*,'Compare sum(DUDT_macro_ref) = ',sum(DUDT_macro_ref)
        print*,'***'
        print*,'Compare sum(diff(DVDT_macro)) = ',sum(DVDT_macro_ref - DVDT_macro)
        print*,'Compare sum(DVDT_macro) = ',sum(DVDT_macro)
        print*,'Compare sum(DVDT_macro_ref) = ',sum(DVDT_macro_ref)
        print*,'***'
        print*,'Compare sum(diff(DTDT_macro)) = ',sum(DTDT_macro_ref - DTDT_macro)
        print*,'Compare sum(DTDT_macro) = ',sum(DTDT_macro)
        print*,'Compare sum(DTDT_macro_ref) = ',sum(DTDT_macro_ref)
        print*,'***'
        print*,'Compare sum(diff(RHCRIT3D)) = ',sum(RHCRIT3D_ref - RHCRIT3D)
        print*,'Compare sum(RHCRIT3D) = ',sum(RHCRIT3D)
        print*,'Compare sum(RHCRIT3D_ref) = ',sum(RHCRIT3D_ref)
        print*,'***'
        print*,'Compare sum(diff(DQVDT_micro)) = ',sum(DQVDT_micro_ref - DQVDT_micro)
        print*,'Compare sum(DQVDT_micro) = ',sum(DQVDT_micro)
        print*,'Compare sum(DQVDT_micro_ref) = ',sum(DQVDT_micro_ref)
        print*,'***'
        print*,'Compare sum(diff(DQIDT_micro)) = ',sum(DQIDT_micro_ref - DQIDT_micro)
        print*,'Compare sum(DQIDT_micro) = ',sum(DQIDT_micro)
        print*,'Compare sum(DQIDT_micro_ref) = ',sum(DQIDT_micro_ref)
        print*,'***'
        print*,'Compare sum(diff(DQLDT_micro)) = ',sum(DQLDT_micro_ref - DQLDT_micro)
        print*,'Compare sum(DQLDT_micro) = ',sum(DQLDT_micro)
        print*,'Compare sum(DQLDT_micro_ref) = ',sum(DQLDT_micro_ref)
        print*,'***'
        print*,'Compare sum(diff(DQADT_micro)) = ',sum(DQADT_micro_ref - DQADT_micro)
        print*,'Compare sum(DQADT_micro) = ',sum(DQADT_micro)
        print*,'Compare sum(DQADT_micro_ref) = ',sum(DQADT_micro_ref)
        print*,'***'
        print*,'Compare sum(diff(DQRDT_micro)) = ',sum(DQRDT_micro_ref - DQRDT_micro)
        print*,'Compare sum(DQRDT_micro) = ',sum(DQRDT_micro)
        print*,'Compare sum(DQRDT_micro_ref) = ',sum(DQRDT_micro_ref)
        print*,'***'
        print*,'Compare sum(diff(DQSDT_micro)) = ',sum(DQSDT_micro_ref - DQSDT_micro)
        print*,'Compare sum(DQSDT_micro) = ',sum(DQSDT_micro)
        print*,'Compare sum(DQSDT_micro_ref) = ',sum(DQSDT_micro_ref)
        print*,'***'
        print*,'Compare sum(diff(DQGDT_micro)) = ',sum(DQGDT_micro_ref - DQGDT_micro)
        print*,'Compare sum(DQGDT_micro) = ',sum(DQGDT_micro)
        print*,'Compare sum(DQGDT_micro_ref) = ',sum(DQGDT_micro_ref)
        print*,'***'
        print*,'Compare sum(diff(DUDT_micro)) = ',sum(DUDT_micro_ref - DUDT_micro)
        print*,'Compare sum(DUDT_micro) = ',sum(DUDT_micro)
        print*,'Compare sum(DUDT_micro_ref) = ',sum(DUDT_micro_ref)
        print*,'***'
        print*,'Compare sum(diff(DVDT_micro)) = ',sum(DVDT_micro_ref - DVDT_micro)
        print*,'Compare sum(DVDT_micro) = ',sum(DVDT_micro)
        print*,'Compare sum(DVDT_micro_ref) = ',sum(DVDT_micro_ref)
        print*,'***'
        print*,'Compare sum(diff(DTDT_micro)) = ',sum(DTDT_micro_ref - DTDT_micro)
        print*,'Compare sum(DTDT_micro) = ',sum(DTDT_micro)
        print*,'Compare sum(DTDT_micro_ref) = ',sum(DTDT_micro_ref)
        print*,'***'

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