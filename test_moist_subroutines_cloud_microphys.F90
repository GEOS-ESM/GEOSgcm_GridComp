module test_cloud_microphys_subroutines

    use moist_subroutines_cloud_microphys
    use timing_module

    implicit none

    public test_gfdl_cloud_microphys_driver

    private

    real, dimension(:,:,:), allocatable :: RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QI_ref, RAD_QS, RAD_QS_ref, RAD_QG, RAD_CF, &
        DQVDTmic, DQVDTmic_ref, DQLDTmic, DQLDTmic_ref, DQRDTmic, DQRDTmic_ref, DQIDTmic, DQIDTmic_ref, &
        DQSDTmic, DQSDTmic_ref, DQGDTmic, DQGDTmic_ref, DQADTmic, DQADTmic_ref, DTDTmic, DTDTmic_ref, &
        W, W_ref, U, V, DUDTmic, DUDTmic_ref, DVDTmic, DVDTmic_ref, DZ, DP, NACTL, NACTI, &
        REV_LS, REV_LS_ref, RSU_LS, RSU_LS_ref, PFL_LS, PFL_LS_ref, PFI_LS, PFI_LS_ref, T

    real, dimension(:,:), allocatable :: AREA, FRLAND, CNV_FRC, SRF_TYPE, PRCP_RAIN, PRCP_RAIN_ref, PRCP_SNOW, PRCP_SNOW_ref, &
    PRCP_ICE, PRCP_ICE_ref, PRCP_GRAUPEL, PRCP_GRAUPEL_ref

    real :: DT_MOIST, ANV_ICEFALL, LS_ICEFALL

    logical :: LHYDROSTATIC, LPHYS_HYDROSTATIC

    contains

    subroutine test_gfdl_cloud_microphys_driver(IM, JM, LM, dirName, rank_str)
        integer :: IM, JM, LM, fileID, itf, ktf, its, ite, kts, kte, mtp, mxp, mzp, ii, nmp, plume
        character*100 :: dirName, rank_str

        print*,'Testing gfdl_cloud_microphys_driver'

        allocate(RAD_QV(IM, JM, LM))
        allocate(RAD_QL(IM, JM, LM))
        allocate(RAD_QR(IM, JM, LM))
        allocate(RAD_QI(IM, JM, LM))
        allocate(RAD_QS(IM, JM, LM))
        allocate(RAD_QG(IM, JM, LM))
        allocate(RAD_CF(IM, JM, LM))
        allocate(DQVDTmic(IM, JM, LM))
        allocate(DQVDTmic_ref(IM, JM, LM))
        allocate(DQLDTmic(IM, JM, LM))
        allocate(DQLDTmic_ref(IM, JM, LM))
        allocate(DQRDTmic(IM, JM, LM))
        allocate(DQRDTmic_ref(IM, JM, LM))
        allocate(DQIDTmic(IM, JM, LM))
        allocate(DQIDTmic_ref(IM, JM, LM))
        allocate(DQSDTmic(IM, JM, LM))
        allocate(DQSDTmic_ref(IM, JM, LM))
        allocate(DQGDTmic(IM, JM, LM))
        allocate(DQGDTmic_ref(IM, JM, LM))
        allocate(DQADTmic(IM, JM, LM))
        allocate(DQADTmic_ref(IM, JM, LM))
        allocate(DTDTmic(IM, JM, LM))
        allocate(T(IM, JM, LM))
        allocate(W(IM, JM, LM))
        allocate(W_ref(IM, JM, LM))
        allocate(U(IM, JM, LM))
        allocate(V(IM, JM, LM))
        allocate(DUDTmic(IM, JM, LM))
        allocate(DUDTmic_ref(IM, JM, LM))
        allocate(DVDTmic(IM, JM, LM))
        allocate(DVDTmic_ref(IM, JM, LM))
        allocate(DZ(IM, JM, LM))
        allocate(DP(IM, JM, LM))
        allocate(AREA(IM, JM))
        allocate(FRLAND(IM, JM))
        allocate(CNV_FRC(IM, JM))
        allocate(SRF_TYPE(IM, JM))
        allocate(NACTL(IM, JM, LM))
        allocate(NACTI(IM, JM, LM))
        allocate(RAD_QI_ref(IM, JM, LM))
        allocate(RAD_QS_ref(IM, JM, LM))
        allocate(DTDTmic_ref(IM, JM, LM))
        allocate(REV_LS(IM, JM, LM))
        allocate(REV_LS_ref(IM, JM, LM))
        allocate(RSU_LS(IM, JM, LM))
        allocate(RSU_LS_ref(IM, JM, LM))
        allocate(PRCP_RAIN(IM, JM))
        allocate(PRCP_RAIN_ref(IM, JM))
        allocate(PRCP_SNOW(IM, JM))
        allocate(PRCP_SNOW_ref(IM, JM))
        allocate(PRCP_ICE(IM, JM))
        allocate(PRCP_ICE_ref(IM, JM))
        allocate(PRCP_GRAUPEL(IM, JM))
        allocate(PRCP_GRAUPEL_ref(IM, JM))
        allocate(PFL_LS(IM, JM, LM))
        allocate(PFL_LS_ref(IM, JM, LM))
        allocate(PFI_LS(IM, JM, LM))
        allocate(PFI_LS_ref(IM, JM, LM))

        open(newunit=fileID, file=trim(dirName) // '/RAD_QV_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) RAD_QV
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In RAD_QV = ', RAD_QV
        
        open(newunit=fileID, file=trim(dirName) // '/RAD_QL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) RAD_QL
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In RAD_QL = ', RAD_QL
        
        open(newunit=fileID, file=trim(dirName) // '/RAD_QR_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) RAD_QR
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In RAD_QR = ', RAD_QR
        
        open(newunit=fileID, file=trim(dirName) // '/RAD_QI_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) RAD_QI
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In RAD_QI = ', RAD_QI
        
        open(newunit=fileID, file=trim(dirName) // '/RAD_QS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) RAD_QS
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In RAD_QS = ', RAD_QS
        
        open(newunit=fileID, file=trim(dirName) // '/RAD_QG_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) RAD_QG
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In RAD_QG = ', RAD_QG
        
        open(newunit=fileID, file=trim(dirName) // '/RAD_CF_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) RAD_CF
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In RAD_CF = ', RAD_CF
        
        open(newunit=fileID, file=trim(dirName) // '/NACTL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) NACTL
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In NACTL = ', NACTL
        
        open(newunit=fileID, file=trim(dirName) // '/NACTI_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) NACTI
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In NACTI = ', NACTI
        
        open(newunit=fileID, file=trim(dirName) // '/DQVDTmic_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DQVDTmic
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DQVDTmic = ', DQVDTmic
        
        open(newunit=fileID, file=trim(dirName) // '/DQLDTmic_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DQLDTmic
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DQLDTmic = ', DQLDTmic
        
        open(newunit=fileID, file=trim(dirName) // '/DQRDTmic_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DQRDTmic
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DQRDTmic = ', DQRDTmic
        
        open(newunit=fileID, file=trim(dirName) // '/DQIDTmic_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DQIDTmic
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DQIDTmic = ', DQIDTmic
        
        open(newunit=fileID, file=trim(dirName) // '/DQSDTmic_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DQSDTmic
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DQSDTmic = ', DQSDTmic
        
        open(newunit=fileID, file=trim(dirName) // '/DQGDTmic_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DQGDTmic
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DQGDTmic = ', DQGDTmic
        
        open(newunit=fileID, file=trim(dirName) // '/DQADTmic_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DQADTmic
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DQADTmic = ', DQADTmic
        
        open(newunit=fileID, file=trim(dirName) // '/DTDTmic_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DTDTmic
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DTDTmic = ', DTDTmic
        
        open(newunit=fileID, file=trim(dirName) // '/T_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) T
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In T = ', T
        
        open(newunit=fileID, file=trim(dirName) // '/W_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) W
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In W = ', W
        
        open(newunit=fileID, file=trim(dirName) // '/U_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) U
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In U = ', U
        
        open(newunit=fileID, file=trim(dirName) // '/V_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) V
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In V = ', V
        
        open(newunit=fileID, file=trim(dirName) // '/DUDTmic_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DUDTmic
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DUDTmic = ', DUDTmic
        
        open(newunit=fileID, file=trim(dirName) // '/DVDTmic_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DVDTmic
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DVDTmic = ', sum(DVDTmic)
        
        open(newunit=fileID, file=trim(dirName) // '/DZ_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DZ
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DZ = ', DZ
        
        open(newunit=fileID, file=trim(dirName) // '/DP_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DP
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DP = ', DP
        
        open(newunit=fileID, file=trim(dirName) // '/AREA_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) AREA
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In AREA = ', AREA
        
        open(newunit=fileID, file=trim(dirName) // '/DT_MOIST_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DT_MOIST
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In DT_MOIST = ', DT_MOIST
        
        open(newunit=fileID, file=trim(dirName) // '/FRLAND_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) FRLAND
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In FRLAND = ', FRLAND
        
        open(newunit=fileID, file=trim(dirName) // '/CNV_FRC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CNV_FRC
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In CNV_FRC = ', CNV_FRC
        
        open(newunit=fileID, file=trim(dirName) // '/SRF_TYPE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) SRF_TYPE
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In SRF_TYPE = ', SRF_TYPE
        
        open(newunit=fileID, file=trim(dirName) // '/ANV_ICEFALL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ANV_ICEFALL
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In ANV_ICEFALL = ', ANV_ICEFALL
        
        open(newunit=fileID, file=trim(dirName) // '/LS_ICEFALL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) LS_ICEFALL
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In LS_ICEFALL = ', LS_ICEFALL
        
        open(newunit=fileID, file=trim(dirName) // '/LHYDROSTATIC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) LHYDROSTATIC
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In LHYDROSTATIC = ', LHYDROSTATIC
        
        open(newunit=fileID, file=trim(dirName) // '/LPHYS_HYDROSTATIC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) LPHYS_HYDROSTATIC
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In LPHYS_HYDROSTATIC = ', LPHYS_HYDROSTATIC

        call update_microphys_constants(dirName, rank_str)

        call start_timing()

        call gfdl_cloud_microphys_driver( &
            ! Input water/cloud species and liquid+ice CCN [NACTL+NACTI]
            RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, (NACTL+NACTI)/1.e6, &
            ! Output tendencies
            DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
            DQSDTmic, DQGDTmic, DQADTmic, DTDTmic, &
            ! Input fields
            T, W, U, V, DUDTmic, DVDTmic, DZ, DP, &
            ! constant inputs
            AREA, DT_MOIST, FRLAND, CNV_FRC, SRF_TYPE, &
            ANV_ICEFALL, LS_ICEFALL, &
            ! Output rain re-evaporation and sublimation
            REV_LS, RSU_LS, & 
            ! Output precipitates
            PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
            ! Output mass flux during sedimentation (Pa kg/kg)
            PFL_LS(:,:,1:LM), PFI_LS(:,:,1:LM), &
            ! constant grid/time information
            LHYDROSTATIC, LPHYS_HYDROSTATIC, &
            1,IM, 1,JM, 1,LM, 1, LM)

        call end_timing()

        call print_timing()

        open(newunit=fileID, file=trim(dirName) // '/RAD_QI_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QI_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out RAD_QI_ref = ', RAD_QI_ref
        
        open(newunit=fileID, file=trim(dirName) // '/RAD_QS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RAD_QS_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out RAD_QS_ref = ', RAD_QS_ref
        
        open(newunit=fileID, file=trim(dirName) // '/DTDTmic_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DTDTmic_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out DTDTmic_ref = ', DTDTmic_ref
        
        open(newunit=fileID, file=trim(dirName) // '/DQADTmic_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQADTmic_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out DQADTmic_ref = ', DQADTmic_ref
        
        open(newunit=fileID, file=trim(dirName) // '/DUDTmic_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DUDTmic_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out DUDTmic_ref = ', DUDTmic_ref
        
        open(newunit=fileID, file=trim(dirName) // '/DVDTmic_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DVDTmic_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out DVDTmic_ref = ', DVDTmic_ref
        
        open(newunit=fileID, file=trim(dirName) // '/W_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) W_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out W_ref = ', W_ref
        
        open(newunit=fileID, file=trim(dirName) // '/DQVDTmic_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQVDTmic_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out DQVDTmic_ref = ', DQVDTmic_ref
        
        open(newunit=fileID, file=trim(dirName) // '/DQLDTmic_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQLDTmic_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out DQLDTmic_ref = ', DQLDTmic_ref
        
        open(newunit=fileID, file=trim(dirName) // '/DQRDTmic_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQRDTmic_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out DQRDTmic_ref = ', DQRDTmic_ref
        
        open(newunit=fileID, file=trim(dirName) // '/DQIDTmic_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQIDTmic_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out DQIDTmic_ref = ', DQIDTmic_ref
        
        open(newunit=fileID, file=trim(dirName) // '/DQSDTmic_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQSDTmic_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out DQSDTmic_ref = ', DQSDTmic_ref
        
        open(newunit=fileID, file=trim(dirName) // '/DQGDTmic_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) DQGDTmic_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out DQGDTmic_ref = ', DQGDTmic_ref
        
        open(newunit=fileID, file=trim(dirName) // '/PRCP_RAIN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PRCP_RAIN_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out PRCP_RAIN_ref = ', PRCP_RAIN_ref
        
        open(newunit=fileID, file=trim(dirName) // '/PRCP_SNOW_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PRCP_SNOW_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out PRCP_SNOW_ref = ', PRCP_SNOW_ref
        
        open(newunit=fileID, file=trim(dirName) // '/PRCP_ICE_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PRCP_ICE_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out PRCP_ICE_ref = ', PRCP_ICE_ref
        
        open(newunit=fileID, file=trim(dirName) // '/PRCP_GRAUPEL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PRCP_GRAUPEL_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out PRCP_GRAUPEL_ref = ', PRCP_GRAUPEL_ref
        
        open(newunit=fileID, file=trim(dirName) // '/PFL_LS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PFL_LS_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out PFL_LS_ref = ', PFL_LS_ref
        
        open(newunit=fileID, file=trim(dirName) // '/PFI_LS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PFI_LS_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out PFI_LS_ref = ', PFI_LS_ref
        
        open(newunit=fileID, file=trim(dirName) // '/REV_LS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) REV_LS_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out REV_LS_ref = ', REV_LS_ref
        
        open(newunit=fileID, file=trim(dirName) // '/RSU_LS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RSU_LS_ref
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': Out RSU_LS_ref = ', RSU_LS_ref

        print*,'Compare sum(diff(RAD_QI)) = ',sum(RAD_QI_ref - RAD_QI)
        print*,'Compare sum(RAD_QI) = ',sum(RAD_QI)
        print*,'Compare sum(RAD_QI_ref) = ',sum(RAD_QI_ref)
        print*,'***'
        print*,'Compare sum(diff(RAD_QS)) = ',sum(RAD_QS_ref - RAD_QS)
        print*,'Compare sum(RAD_QS) = ',sum(RAD_QS)
        print*,'Compare sum(RAD_QS_ref) = ',sum(RAD_QS_ref)
        print*,'***'
        print*,'Compare sum(diff(DTDTmic)) = ',sum(DTDTmic_ref - DTDTmic)
        print*,'Compare sum(DTDTmic) = ',sum(DTDTmic)
        print*,'Compare sum(DTDTmic_ref) = ',sum(DTDTmic_ref)
        print*,'***'
        print*,'Compare sum(diff(DQADTmic)) = ',sum(DQADTmic_ref - DQADTmic)
        print*,'Compare sum(DQADTmic) = ',sum(DQADTmic)
        print*,'Compare sum(DQADTmic_ref) = ',sum(DQADTmic_ref)
        print*,'***'
        print*,'Compare sum(diff(DUDTmic)) = ',sum(DUDTmic_ref - DUDTmic)
        print*,'Compare sum(DUDTmic) = ',sum(DUDTmic)
        print*,'Compare sum(DUDTmic_ref) = ',sum(DUDTmic_ref)
        print*,'***'
        print*,'Compare sum(diff(DVDTmic)) = ',sum(DVDTmic_ref - DVDTmic)
        print*,'Compare sum(DVDTmic) = ',sum(DVDTmic)
        print*,'Compare sum(DVDTmic_ref) = ',sum(DVDTmic_ref)
        print*,'***'
        print*,'Compare sum(diff(W)) = ',sum(W_ref - W)
        print*,'Compare sum(W) = ',sum(W)
        print*,'Compare sum(W_ref) = ',sum(W_ref)
        print*,'***'
        print*,'Compare sum(diff(DQVDTmic)) = ',sum(DQVDTmic_ref - DQVDTmic)
        print*,'Compare sum(DQVDTmic) = ',sum(DQVDTmic)
        print*,'Compare sum(DQVDTmic_ref) = ',sum(DQVDTmic_ref)
        print*,'***'
        print*,'Compare sum(diff(DQLDTmic)) = ',sum(DQLDTmic_ref - DQLDTmic)
        print*,'Compare sum(DQLDTmic) = ',sum(DQLDTmic)
        print*,'Compare sum(DQLDTmic_ref) = ',sum(DQLDTmic_ref)
        print*,'***'
        print*,'Compare sum(diff(DQRDTmic)) = ',sum(DQRDTmic_ref - DQRDTmic)
        print*,'Compare sum(DQRDTmic) = ',sum(DQRDTmic)
        print*,'Compare sum(DQRDTmic_ref) = ',sum(DQRDTmic_ref)
        print*,'***'
        print*,'Compare sum(diff(DQIDTmic)) = ',sum(DQIDTmic_ref - DQIDTmic)
        print*,'Compare sum(DQIDTmic) = ',sum(DQIDTmic)
        print*,'Compare sum(DQIDTmic_ref) = ',sum(DQIDTmic_ref)
        print*,'***'
        print*,'Compare sum(diff(DQSDTmic)) = ',sum(DQSDTmic_ref - DQSDTmic)
        print*,'Compare sum(DQSDTmic) = ',sum(DQSDTmic)
        print*,'Compare sum(DQSDTmic_ref) = ',sum(DQSDTmic_ref)
        print*,'***'
        print*,'Compare sum(diff(DQGDTmic)) = ',sum(DQGDTmic_ref - DQGDTmic)
        print*,'Compare sum(DQGDTmic) = ',sum(DQGDTmic)
        print*,'Compare sum(DQGDTmic_ref) = ',sum(DQGDTmic_ref)
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
        print*,'Compare sum(diff(PFL_LS)) = ',sum(PFL_LS_ref - PFL_LS)
        print*,'Compare sum(PFL_LS) = ',sum(PFL_LS)
        print*,'Compare sum(PFL_LS_ref) = ',sum(PFL_LS_ref)
        print*,'***'
        print*,'Compare sum(diff(PFI_LS)) = ',sum(PFI_LS_ref - PFI_LS)
        print*,'Compare sum(PFI_LS) = ',sum(PFI_LS)
        print*,'Compare sum(PFI_LS_ref) = ',sum(PFI_LS_ref)
        print*,'***'
        print*,'Compare sum(diff(REV_LS)) = ',sum(REV_LS_ref - REV_LS)
        print*,'Compare sum(REV_LS) = ',sum(REV_LS)
        print*,'Compare sum(REV_LS_ref) = ',sum(REV_LS_ref)
        print*,'***'
        print*,'Compare sum(diff(RSU_LS)) = ',sum(RSU_LS_ref - RSU_LS)
        print*,'Compare sum(RSU_LS) = ',sum(RSU_LS)
        print*,'Compare sum(RSU_LS_ref) = ',sum(RSU_LS_ref)
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