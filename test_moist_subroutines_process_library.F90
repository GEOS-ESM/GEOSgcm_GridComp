module test_process_library_subroutines

    use Process_Library_standalone
    use timing_module

    implicit none

    public test_FILLQ2ZERO, test_BUOYANCY, test_hystpdf

    private

    real, dimension(:,:,:), allocatable :: MASS, Q, Q_comp
    real, dimension(:,:,:), allocatable :: T, QST3, DQST3, DZET, ZL0, BYNCY, BYNCY_comp
    real, dimension(:,:),   allocatable :: TMP2D, TMP2D_comp
    real, dimension(:,:),   allocatable :: CAPE, CAPE_comp, INHB, INHB_comp

    real :: DT_MOIST, ALPHA, CNV_FRC_IJ, SRF_TYPE_IJ, PLmb_IJL, ZL0_IJL, Q_IJL, Q_IJL_ref, &
            QLLS_IJL, QLLS_IJL_ref, QLCN_IJL, QLCN_IJL_ref, QILS_IJL, QILS_IJL_ref, &
            QICN_IJL, QICN_IJL_ref, T_IJL, T_IJL_ref, CLLS_IJL, CLLS_IJL_ref, CLCN_IJL, CLCN_IJL_ref, &
            NACTL_IJL, NACTI_IJL, WHL_IJL, WQT_IJL, HL2_IJL, QT2_IJL, HLQT_IJL, W3_IJL, W2_IJL, &
            QT3_IJL, HL3_IJL, EDMF_FRC_IJL, PDF_A_IJL, PDF_A_IJL_ref, PDFITERS_IJL, PDFITERS_IJL_ref, &
            WTHV2_IJL, WTHV2_IJL_ref, WQL_IJL, WQL_IJL_ref

    integer :: PDFSHAPE

    logical :: USE_BERGERON
    contains

    subroutine test_FILLQ2ZERO(IM, JM, LM, dirName, rank_str)

        integer :: IM, JM, LM, fileID
        character*100 :: dirName, rank_str

        allocate(Q(IM, JM, LM))
        allocate(Q_comp(IM, JM, LM))
        allocate(MASS(IM, JM, LM))
        allocate(TMP2D(IM, JM))
        allocate(TMP2D_comp(IM, JM))

        open(newunit=fileID, file=trim(dirName) // "/Q_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) Q
        close(fileID)
        ! write(*,*) 'sum(abs(Q)) = ', sum(abs(Q))

        open(newunit=fileID, file=trim(dirName) // "/MASS_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) MASS
        close(fileID)
        ! write(*,*) 'sum(abs(MASS)) = ', sum(abs(MASS))

        print*,'Testing Fillq2zero Subroutine'

        call start_timing()
        call FILLQ2ZERO(Q, MASS, TMP2D)
        call end_timing()
        call print_timing()

        open(newunit=fileID, file=trim(dirName) // "/Q_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) Q_comp
        close(fileID)
        

        open(newunit=fileID, file=trim(dirName) // "/TMP2D_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID)  TMP2D_comp
        close(fileID)

        write(*,*) 'sum(Q-Q_comp) = ', sum(Q-Q_comp)
        write(*,*) 'sum(TMP2D-TMP2D_comp) = ', sum(TMP2D-TMP2D_comp)

    end subroutine

    subroutine test_BUOYANCY(IM, JM, LM, dirName, rank_str)

        integer :: IM, JM, LM, fileID
        character*100 :: dirName, rank_str

        allocate(T(IM, JM, LM))
        allocate(Q(IM, JM, LM))
        allocate(QST3(IM, JM, LM))
        allocate(DQST3(IM, JM, LM))
        allocate(DZET(IM, JM, LM))
        allocate(ZL0(IM, JM, LM))
        allocate(BYNCY(IM, JM, LM))
        allocate(BYNCY_comp(IM, JM, LM))
        allocate(CAPE(IM, JM))
        allocate(CAPE_comp(IM, JM))
        allocate(INHB(IM, JM))
        allocate(INHB_comp(IM, JM))

        open(newunit=fileID, file=trim(dirName) // "/T_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) T
        close(fileID)
        ! write(*,*) 'sum(abs(T)) = ', sum(abs(T))

        open(newunit=fileID, file=trim(dirName) // "/Q_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) Q
        close(fileID)
        ! write(*,*) 'sum(abs(Q)) = ', sum(abs(Q))

        open(newunit=fileID, file=trim(dirName) // "/QST3_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) QST3
        close(fileID)
        ! write(*,*) 'sum(abs(QST3)) = ', sum(abs(QST3))

        open(newunit=fileID, file=trim(dirName) // "/DQST3_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) DQST3
        close(fileID)
        ! write(*,*) 'sum(abs(DQST3)) = ', sum(abs(DQST3))

        open(newunit=fileID, file=trim(dirName) // "/DZET_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) DZET
        close(fileID)
        ! write(*,*) 'sum(abs(DZET)) = ', sum(abs(DZET))

        open(newunit=fileID, file=trim(dirName) // "/ZL0_" // trim(rank_str) // ".in", &
            status='old', form="unformatted", action="read")
        read(fileID) ZL0
        close(fileID)
        ! write(*,*) 'sum(abs(ZL0)) = ', sum(abs(ZL0))

!!$acc data copyin(T, Q, QST3, DQST3, DZET, ZL0) &
!!acc      copyout(BYNCY, CAPE, INHB)

!$omp target data map(to:T, Q, QST3, DQST3, DZET, ZL0) &
!$omp             map(from:BYNCY, CAPE, INHB)
        print*,'Testing Buoyancy Subroutine'

        call start_timing()
        call BUOYANCY( T, Q, QST3, DQST3, DZET, ZL0, BYNCY, CAPE, INHB)
        call end_timing()
        call print_timing()

!$omp end target data
!!$acc end data
        open(newunit=fileID, file=trim(dirName) // "/BYNCY_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) BYNCY_comp
        close(fileID)
        ! write(*,*) 'sum(abs(BYNCY)) = ', sum(abs(BYNCY_comp))

        open(newunit=fileID, file=trim(dirName) // "/CAPE_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) CAPE_comp
        close(fileID)
        ! write(*,*) 'sum(abs(CAPE)) = ', sum(abs(CAPE_comp))

        open(newunit=fileID, file=trim(dirName) // "/INHB_" // trim(rank_str) // ".out", &
            status='old', form="unformatted", action="read")
        read(fileID) INHB_comp
        close(fileID)
        ! write(*,*) 'sum(abs(INHB)) = ', sum(abs(INHB_comp))

        print*,'Compare sum(diff(INHB)) = ',sum(INHB - INHB_comp)
        print*,'Compare sum(INHB) = ',sum(INHB)
        print*,'Compare sum(INHB_comp) = ',sum(INHB_comp)
        print*,'***'
        print*,'Compare sum(diff(CAPE)) = ',sum(CAPE - CAPE_comp)
        print*,'Compare sum(CAPE) = ',sum(CAPE)
        print*,'Compare sum(CAPE_comp) = ',sum(CAPE_comp)
        print*,'***'
        print*,'Compare sum(diff(BYNCY)) = ',sum(BYNCY - BYNCY_comp)
        print*,'Compare sum(BYNCY) = ',sum(BYNCY)
        print*,'Compare sum(BYNCY_comp) = ',sum(BYNCY_comp)
        print*,'***'
        
        ! Encoding results for CI
        print*, '#CI#VAR|INHB#DIFF|',sum(INHB - INHB_comp)
        print*, '#CI#VAR|INHB#NEW|',sum(INHB)
        print*, '#CI#VAR|INHB#REF|',sum(INHB_comp)
        print*, '#CI#VAR|INHB#THRSH|',sum(INHB_comp)*1.0e-9
        
        print*, '#CI#VAR|CAPE#DIFF|',sum(CAPE - CAPE_comp)
        print*, '#CI#VAR|CAPE#NEW|',sum(CAPE)
        print*, '#CI#VAR|CAPE#REF|',sum(CAPE_comp)
        print*, '#CI#VAR|CAPE#THRSH|',sum(CAPE_comp)*1.0e-9
        
        print*, '#CI#VAR|BYNCY#DIFF|',sum(BYNCY - BYNCY_comp)
        print*, '#CI#VAR|BYNCY#NEW|',sum(BYNCY)
        print*, '#CI#VAR|BYNCY#REF|',sum(BYNCY_comp)
        print*, '#CI#VAR|BYNCY#THRSH|',sum(BYNCY_comp)*1.0e-9

    end subroutine

    subroutine test_hystpdf(IM, JM, LM, dirName, rank_str)

        integer :: IM, JM, LM, fileID
        character*100 :: dirName, rank_str

        print*, 'Testing hystpdf subroutine'

        open(newunit=fileID, file=trim(dirName) // '/DT_MOIST_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DT_MOIST
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In DT_MOIST = ', DT_MOIST

        open(newunit=fileID, file=trim(dirName) // '/ALPHA_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ALPHA
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In ALPHA = ', ALPHA

        open(newunit=fileID, file=trim(dirName) // '/PDFSHAPE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) PDFSHAPE
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In PDFSHAPE = ', PDFSHAPE

        open(newunit=fileID, file=trim(dirName) // '/CNV_FRC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CNV_FRC_IJ
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In CNV_FRC_IJ = ', CNV_FRC_IJ

        open(newunit=fileID, file=trim(dirName) // '/SRF_TYPE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) SRF_TYPE_IJ
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In SRF_TYPE_IJ = ', SRF_TYPE_IJ

        open(newunit=fileID, file=trim(dirName) // '/PLmb_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) PLmb_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In PLmb_IJL = ', PLmb_IJL

        open(newunit=fileID, file=trim(dirName) // '/ZL0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ZL0_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In ZL0_IJL = ', ZL0_IJL

        open(newunit=fileID, file=trim(dirName) // '/Q_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) Q_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In Q_IJL = ', Q_IJL

        open(newunit=fileID, file=trim(dirName) // '/QLLS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QLLS_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In QLLS_IJL = ', QLLS_IJL

        open(newunit=fileID, file=trim(dirName) // '/QLCN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QLCN_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In QLCN_IJL = ', QLCN_IJL

        open(newunit=fileID, file=trim(dirName) // '/QILS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QILS_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In QILS_IJL = ', QILS_IJL

        open(newunit=fileID, file=trim(dirName) // '/QICN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QICN_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In QICN_IJL = ', QICN_IJL

        open(newunit=fileID, file=trim(dirName) // '/T_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) T_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In T_IJL = ', T_IJL

        open(newunit=fileID, file=trim(dirName) // '/CLLS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CLLS_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In CLLS_IJL = ', CLLS_IJL

        open(newunit=fileID, file=trim(dirName) // '/CLCN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CLCN_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In CLCN_IJL = ', CLCN_IJL

        open(newunit=fileID, file=trim(dirName) // '/NACTL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) NACTL_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In NACTL_IJL = ', NACTL_IJL

        open(newunit=fileID, file=trim(dirName) // '/NACTI_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) NACTI_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In NACTI_IJL = ', NACTI_IJL

        open(newunit=fileID, file=trim(dirName) // '/WHL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) WHL_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In WHL_IJL = ', WHL_IJL

        open(newunit=fileID, file=trim(dirName) // '/WQT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) WQT_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In WQT_IJL = ', WQT_IJL

        open(newunit=fileID, file=trim(dirName) // '/HL2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) HL2_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In HL2_IJL = ', HL2_IJL

        open(newunit=fileID, file=trim(dirName) // '/QT2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QT2_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In QT2_IJL = ', QT2_IJL

        open(newunit=fileID, file=trim(dirName) // '/HLQT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) HLQT_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In HLQT_IJL = ', HLQT_IJL

        open(newunit=fileID, file=trim(dirName) // '/W3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) W3_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In W3_IJL = ', W3_IJL

        open(newunit=fileID, file=trim(dirName) // '/W2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) W2_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In W2_IJL = ', W2_IJL

        open(newunit=fileID, file=trim(dirName) // '/QT3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QT3_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In QT3_IJL = ', QT3_IJL

        open(newunit=fileID, file=trim(dirName) // '/HL3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) HL3_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In HL3_IJL = ', HL3_IJL

        open(newunit=fileID, file=trim(dirName) // '/EDMF_FRC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) EDMF_FRC_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In EDMF_FRC_IJL = ', EDMF_FRC_IJL

        open(newunit=fileID, file=trim(dirName) // '/PDF_A_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) PDF_A_IJL
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In PDF_A_IJL = ', PDF_A_IJL

        open(newunit=fileID, file=trim(dirName) // '/USE_BERGERON_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) USE_BERGERON
        close(fileID)
        ! print*, 'Rank ', trim(rank_str),': In USE_BERGERON = ', USE_BERGERON

        call start_timing()
        call hystpdf( &
                      DT_MOIST       , &
                      ALPHA          , &
                      PDFSHAPE       , &
                      CNV_FRC_IJ   , &
                      SRF_TYPE_IJ  , &
                      PLmb_IJL    , &
                      ZL0_IJL     , &
                      Q_IJL       , &
                      QLLS_IJL    , &
                      QLCN_IJL    , &
                      QILS_IJL    , &
                      QICN_IJL    , &
                      T_IJL       , &
                      CLLS_IJL    , &
                      CLCN_IJL    , &
                      NACTL_IJL   , &
                      NACTI_IJL   , &
                      WHL_IJL     , &
                      WQT_IJL     , &
                      HL2_IJL     , &
                      QT2_IJL     , &
                      HLQT_IJL    , &
                      W3_IJL      , &
                      W2_IJL      , &
                      QT3_IJL     , &
                      HL3_IJL     , &
                      EDMF_FRC_IJL, &
                      PDF_A_IJL   , &
                      PDFITERS_IJL, &
                      WTHV2_IJL   , &
                      WQL_IJL     , &
                      .false.        , & 
                      USE_BERGERON)
        call end_timing()
        call print_timing()

        open(newunit=fileID, file=trim(dirName) // '/Q_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) Q_IJL_ref
        close(fileID)
        print*, 'Rank ', trim(rank_str),': Out Q_IJL_ref = ', Q_IJL_ref
        
        open(newunit=fileID, file=trim(dirName) // '/QLLS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) QLLS_IJL_ref
        close(fileID)
        print*, 'Rank ', trim(rank_str),': Out QLLS_IJL_ref = ', QLLS_IJL_ref
        
        open(newunit=fileID, file=trim(dirName) // '/QLCN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) QLCN_IJL_ref
        close(fileID)
        print*, 'Rank ', trim(rank_str),': Out QLCN_IJL_ref = ', QLCN_IJL_ref
        
        open(newunit=fileID, file=trim(dirName) // '/QILS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) QILS_IJL_ref
        close(fileID)
        print*, 'Rank ', trim(rank_str),': Out QILS_IJL_ref = ', QILS_IJL_ref
        
        open(newunit=fileID, file=trim(dirName) // '/QICN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) QICN_IJL_ref
        close(fileID)
        print*, 'Rank ', trim(rank_str),': Out QICN_IJL_ref = ', QICN_IJL_ref
        
        open(newunit=fileID, file=trim(dirName) // '/T_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) T_IJL_ref
        close(fileID)
        print*, 'Rank ', trim(rank_str),': Out T_IJL_ref = ', T_IJL_ref
        
        open(newunit=fileID, file=trim(dirName) // '/CLLS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) CLLS_IJL_ref
        close(fileID)
        print*, 'Rank ', trim(rank_str),': Out CLLS_IJL_ref = ', CLLS_IJL_ref
        
        open(newunit=fileID, file=trim(dirName) // '/CLCN_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) CLCN_IJL_ref
        close(fileID)
        print*, 'Rank ', trim(rank_str),': Out CLCN_IJL_ref = ', CLCN_IJL_ref
        
        open(newunit=fileID, file=trim(dirName) // '/PDF_A_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PDF_A_IJL_ref
        close(fileID)
        print*, 'Rank ', trim(rank_str),': Out PDF_A_IJL_ref = ', PDF_A_IJL_ref
        
        open(newunit=fileID, file=trim(dirName) // '/PDFITERS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PDFITERS_IJL_ref
        close(fileID)
        print*, 'Rank ', trim(rank_str),': Out PDFITERS_IJL_ref = ', PDFITERS_IJL_ref
        
        open(newunit=fileID, file=trim(dirName) // '/WTHV2_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) WTHV2_IJL_ref
        close(fileID)
        print*, 'Rank ', trim(rank_str),': Out WTHV2_IJL_ref = ', WTHV2_IJL_ref
        
        open(newunit=fileID, file=trim(dirName) // '/WQL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) WQL_IJL_ref
        close(fileID)
        print*, 'Rank ', trim(rank_str),': Out WQL_IJL_ref = ', WQL_IJL_ref

        print*,'Compare Q_IJL : ',Q_IJL_ref-Q_IJL
        print*,'Compare QLLS_IJL : ',QLLS_IJL_ref-QLLS_IJL
        print*,'Compare QLCN_IJL : ',QLCN_IJL_ref-QLCN_IJL
        print*,'Compare QILS_IJL : ',QILS_IJL_ref-QILS_IJL
        print*,'Compare QICN_IJL : ',QICN_IJL_ref-QICN_IJL
        print*,'Compare T_IJL : ',T_IJL_ref-T_IJL
        print*,'Compare CLLS_IJL : ',CLLS_IJL_ref-CLLS_IJL
        print*,'Compare CLCN_IJL : ',CLCN_IJL_ref-CLCN_IJL
        print*,'Compare PDF_A_IJL : ',PDF_A_IJL_ref-PDF_A_IJL
        print*,'Compare PDFITERS_IJL : ',PDFITERS_IJL_ref-PDFITERS_IJL
        print*,'Compare WTHV2_IJL : ',WTHV2_IJL_ref-WTHV2_IJL
        print*,'Compare WQL_IJL : ',WQL_IJL_ref-WQL_IJL
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