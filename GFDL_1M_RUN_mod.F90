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

    real, pointer, dimension(:,:  ), public :: PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL
    real, pointer, dimension(:,:  ), public :: LS_PRCP, LS_SNR, ICE, FRZR, CNV_FRC, SRF_TYPE
    real, pointer, dimension(:,:,:), public :: DQVDT_macro, DQIDT_macro, DQLDT_macro, DQADT_macro, DQRDT_macro, DQSDT_macro, DQGDT_macro
    real, pointer, dimension(:,:,:), public :: DUDT_macro,  DVDT_macro,  DTDT_macro
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

    ! Module variable from aer_actv_single_moment.F90
    logical, public :: USE_BERGERON
 
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

    subroutine GFDL_1M_RUN_driver(IM, JM, LM)

        integer, intent(in) :: IM,JM,LM

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

        ! Derived States
        PLEmb    =  PLE*.01
        PLmb     = 0.5*(PLEmb(:,:,0:LM-1) + PLEmb(:,:,1:LM))
        DO L=0,LM
        ZLE0(:,:,L)= ZLE(:,:,L) - ZLE(:,:,LM)   ! Edge Height (m) above the surface
        END DO
        ZL0      = 0.5*(ZLE0(:,:,0:LM-1) + ZLE0(:,:,1:LM) ) ! Layer Height (m) above the surface
        DZET     =     (ZLE0(:,:,0:LM-1) - ZLE0(:,:,1:LM) ) ! Layer thickness (m)
        DQST3    = GEOS_DQSAT(T, PLmb, QSAT=QST3)
        DP       = ( PLE(:,:,1:LM)-PLE(:,:,0:LM-1) )
        MASS     = DP/MAPL_GRAV
        iMASS    = 1.0/MASS
        U0       = U
        V0       = V

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