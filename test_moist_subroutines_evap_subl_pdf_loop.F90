module test_evap_subl_pdf_subroutines

    use evap_subl_pdf_loop
    use timing_module

    implicit none

    public test_evap_subl_pdf_loop

    private

    contains

    subroutine test_evap_subl_pdf_loop(IM, JM, LM, dirName, rank_str)

        implicit none

        integer, intent(IN) :: IM, JM, LM

        character*100, intent(IN) :: dirName, rank_str

        real :: DT_MOIST, CCW_EVAP_EFF, CCI_EVAP_EFF, dw_ocean, dw_land, turnrhcrit

        integer :: PDFSHAPE, fileID

        logical :: do_qa, USE_BERGERON

        integer, dimension(:,:), allocatable :: KLCL
        real, dimension(:,:), allocatable :: CNV_FRC, SRF_TYPE, EIS, AREA
        real, dimension(:,:,:), allocatable :: PLmb, PLEmb, ZL0, NACTL, NACTI, WHL, WQT, HL2, &
            QT2, HLQT, W3, W2, QT3, HL3, EDMF_FRC, QST3
        real, dimension(:,:,:), allocatable, target     :: RHCRIT3D
        real, dimension(:,:,:), allocatable ::  PDFITERS, WTHV2, WQL, Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN, PDF_A
        real, dimension(:,:,:), allocatable ::  PDFITERS_ref, WTHV2_ref, WQL_ref, Q_ref, T_ref, QLLS_ref, QILS_ref, &
            CLLS_ref, QLCN_ref, QICN_ref, CLCN_ref, PDF_A_ref, SUBLC, SUBLC_ref, RHX, RHX_ref, EVAPC, EVAPC_ref

        real, dimension(:,:,:), pointer :: RHCRIT3D_ptr

        allocate(CNV_FRC     (IM,JM))
        allocate(SRF_TYPE    (IM,JM))
        allocate(EIS(IM, JM))
        allocate(KLCL(IM, JM))
        allocate(AREA(IM, JM))

        allocate(PLmb  (IM,JM,LM))
        ALLOCATE(PLEmb (IM,JM,0:LM) )
        allocate(ZL0   (IM, JM, LM))
        allocate(NACTL(IM, JM, LM))
        allocate(NACTI(IM, JM, LM))
        allocate(WHL  (IM, JM, LM))
        allocate(WQT  (IM, JM, LM))
        allocate(HL2  (IM, JM, LM))
        allocate(QT2  (IM, JM, LM))
        allocate(HLQT (IM, JM, LM))
        allocate(W3   (IM, JM, LM))
        allocate(W2   (IM, JM, LM))
        allocate(QT3  (IM, JM, LM))
        allocate(HL3  (IM, JM, LM))
        allocate(EDMF_FRC (IM, JM, LM))
        allocate(QST3 (IM, JM, LM))
        allocate(PDFITERS   (IM, JM, LM))
        allocate(WTHV2   (IM, JM, LM))
        allocate(WQL (IM, JM, LM))
        allocate(Q    (IM, JM, LM))
        allocate(T    (IM, JM, LM))
        allocate(QLLS   (IM, JM, LM))
        allocate(QILS   (IM, JM, LM))
        allocate(CLLS   (IM, JM, LM))
        allocate(QLCN   (IM, JM, LM))
        allocate(QICN   (IM, JM, LM))
        allocate(CLCN   (IM, JM, LM))
        allocate(PDF_A   (IM, JM, LM))
        allocate(RHCRIT3D (IM, JM, LM))

        allocate(PDFITERS_ref   (IM, JM, LM))
        allocate(WTHV2_ref   (IM, JM, LM))
        allocate(WQL_ref (IM, JM, LM))
        allocate(Q_ref    (IM, JM, LM))
        allocate(T_ref    (IM, JM, LM))
        allocate(QLLS_ref   (IM, JM, LM))
        allocate(QILS_ref   (IM, JM, LM))
        allocate(CLLS_ref   (IM, JM, LM))
        allocate(QLCN_ref   (IM, JM, LM))
        allocate(QICN_ref   (IM, JM, LM))
        allocate(CLCN_ref   (IM, JM, LM))
        allocate(PDF_A_ref   (IM, JM, LM))
        allocate(SUBLC(IM, JM, LM))
        allocate(RHX(IM, JM, LM))
        allocate(EVAPC(IM, JM, LM))
        allocate(SUBLC_ref(IM, JM, LM))
        allocate(RHX_ref(IM, JM, LM))
        allocate(EVAPC_ref(IM, JM, LM))

        RHCRIT3D_ptr => RHCRIT3D

        open(newunit=fileID, file=trim(dirName) // '/DT_MOIST_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) DT_MOIST
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(DT_MOIST) = ', DT_MOIST

        open(newunit=fileID, file=trim(dirName) // '/PDFSHAPE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) PDFSHAPE
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': PDFSHAPE = ', PDFSHAPE

        open(newunit=fileID, file=trim(dirName) // '/CNV_FRC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CNV_FRC
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(CNV_FRC) = ', sum(CNV_FRC)

        open(newunit=fileID, file=trim(dirName) // '/SRF_TYPE_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) SRF_TYPE
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(SRF_TYPE) = ', sum(SRF_TYPE)

        open(newunit=fileID, file=trim(dirName) // '/PLmb_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) PLmb
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(PLmb) = ', sum(PLmb)

        open(newunit=fileID, file=trim(dirName) // '/ZL0_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) ZL0
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(ZL0) = ', sum(ZL0)

        open(newunit=fileID, file=trim(dirName) // '/Q_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) Q
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(Q) = ', sum(Q)

        open(newunit=fileID, file=trim(dirName) // '/QLLS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QLLS
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(QLLS) = ', sum(QLLS)

        open(newunit=fileID, file=trim(dirName) // '/QLCN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QLCN
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(QLCN) = ', sum(QLCN)

        open(newunit=fileID, file=trim(dirName) // '/QILS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QILS
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(QILS) = ', sum(QILS)

        open(newunit=fileID, file=trim(dirName) // '/QICN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QICN
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(QICN) = ', sum(QICN)

        open(newunit=fileID, file=trim(dirName) // '/T_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) T
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(T) = ', sum(T)

        open(newunit=fileID, file=trim(dirName) // '/CLLS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CLLS
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(CLLS) = ', sum(CLLS)

        open(newunit=fileID, file=trim(dirName) // '/CLCN_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CLCN
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(CLCN) = ', sum(CLCN)

        open(newunit=fileID, file=trim(dirName) // '/NACTL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) NACTL
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(NACTL) = ', sum(NACTL)

        open(newunit=fileID, file=trim(dirName) // '/NACTI_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) NACTI
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(NACTI) = ', sum(NACTI)

        open(newunit=fileID, file=trim(dirName) // '/WHL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) WHL
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(WHL) = ', sum(WHL)

        open(newunit=fileID, file=trim(dirName) // '/WQT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) WQT
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(WQT) = ', sum(WQT)

        open(newunit=fileID, file=trim(dirName) // '/HL2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) HL2
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(HL2) = ', sum(HL2)

        open(newunit=fileID, file=trim(dirName) // '/QT2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QT2
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(QT2) = ', sum(QT2)

        open(newunit=fileID, file=trim(dirName) // '/HLQT_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) HLQT
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(HLQT) = ', sum(HLQT)

        open(newunit=fileID, file=trim(dirName) // '/W3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) W3
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(W3) = ', sum(W3)

        open(newunit=fileID, file=trim(dirName) // '/W2_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) W2
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(W2) = ', sum(W2)

        open(newunit=fileID, file=trim(dirName) // '/QT3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QT3
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(QT3) = ', sum(QT3)

        open(newunit=fileID, file=trim(dirName) // '/HL3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) HL3
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(HL3) = ', sum(HL3)

        open(newunit=fileID, file=trim(dirName) // '/EDMF_FRC_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) EDMF_FRC
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(EDMF_FRC) = ', sum(EDMF_FRC)

        open(newunit=fileID, file=trim(dirName) // '/PDF_A_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) PDF_A
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(PDF_A) = ', sum(PDF_A)

        open(newunit=fileID, file=trim(dirName) // '/USE_BERGERON_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) USE_BERGERON
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In USE_BERGERON = ', USE_BERGERON

        open(newunit=fileID, file=trim(dirName) // '/CCW_EVAP_EFF_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CCW_EVAP_EFF
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(CCW_EVAP_EFF) = ', CCW_EVAP_EFF

        open(newunit=fileID, file=trim(dirName) // '/QST3_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) QST3
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(QST3) = ', sum(QST3)

        open(newunit=fileID, file=trim(dirName) // '/CCI_EVAP_EFF_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) CCI_EVAP_EFF
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(CCI_EVAP_EFF) = ', CCI_EVAP_EFF

        open(newunit=fileID, file=trim(dirName) // '/PLEmb_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) PLEmb
        close(fileID)
        ! write(*,*) 'Rank ', trim(rank_str),': In sum(PLEmb) = ', sum(PLEmb)

        open(newunit=fileID, file=trim(dirName) // '/EIS_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) EIS
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/dw_ocean_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dw_ocean
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/dw_land_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) dw_land
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/KLCL_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) KLCL
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RHCRIT3D_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) RHCRIT3D
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/AREA_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) AREA
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/turnrhcrit_' // trim(rank_str) // '.in', status='old', form='unformatted', action='read')
        read(fileID) turnrhcrit
        close(fileID)

        ! Note : do_qa is manually set and may not represent a particular problem setup
        do_qa = .false.

        !$acc data copyin(DT_MOIST, PDFSHAPE, CNV_FRC, SRF_TYPE, PLmb, ZL0, &
        !$acc             NACTL, NACTI, WHL, WQT, HL2, QT2, &
        !$acc             HLQT, W3, W2, QT3, HL3, EDMF_FRC, USE_BERGERON, &
        !$acc             QST3, CCI_EVAP_EFF, CCW_EVAP_EFF, &
        !$acc             PLEmb) &
        !$acc      copy(Q, QLLS, QLCN, QILS, QICN, T, CLLS, &
        !$acc           CLCN, PDF_A, WTHV2, WQL) &
        !$acc      copyout(EVAPC,RHX, SUBLC, PDFITERS)

        call start_timing()

        call evap_subl_pdf_loop_standalone(IM, JM, LM, do_qa, USE_BERGERON, PDFSHAPE, DT_MOIST, CCW_EVAP_EFF, &
            CCI_EVAP_EFF, dw_ocean, dw_land, turnrhcrit, CNV_FRC, SRF_TYPE, EIS, KLCL, AREA, &
            PLmb, PLEmb, ZL0, NACTL, NACTI, WHL, WQT, HL2, &
            QT2, HLQT, W3, W2, QT3, HL3, EDMF_FRC, QST3, PDFITERS, WTHV2, WQL, &
            Q, T, QLLS, QILS, CLLS, QLCN, QICN, CLCN, PDF_A, SUBLC, RHX, EVAPC, RHCRIT3D_ptr)

        call end_timing()

        !$acc end data

        call print_timing()

        open(newunit=fileID, file=trim(dirName) // '/PDFITERS_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) PDFITERS_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/WTHV2_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) WTHV2_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/WQL_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) WQL_ref
        close(fileID)
        
        open(newunit=fileID, file=trim(dirName) // '/Q_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) Q_ref
        close(fileID)
        
        open(newunit=fileID, file=trim(dirName) // '/T_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) T_ref
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

        open(newunit=fileID, file=trim(dirName) // '/SUBLC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) SUBLC_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/RHX_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) RHX_ref
        close(fileID)

        open(newunit=fileID, file=trim(dirName) // '/EVAPC_' // trim(rank_str) // '.out', status='old', form='unformatted', action='read')
        read(fileID) EVAPC_ref
        close(fileID)

        print*,'Compare sum(diff(PDFITERS)) = ',sum(PDFITERS_ref - PDFITERS)
        print*,'Compare sum(PDFITERS) = ',sum(PDFITERS)
        print*,'Compare sum(PDFITERS_ref) = ',sum(PDFITERS_ref)
        print*,'***'
        print*,'Compare sum(diff(WTHV2)) = ',sum(WTHV2_ref - WTHV2)
        print*,'Compare sum(WTHV2) = ',sum(WTHV2)
        print*,'Compare sum(WTHV2_ref) = ',sum(WTHV2_ref)
        print*,'***'
        print*,'Compare sum(diff(WQL)) = ',sum(WQL_ref - WQL)
        print*,'Compare sum(WQL) = ',sum(WQL)
        print*,'Compare sum(WQL_ref) = ',sum(WQL_ref)
        print*,'***'
        print*,'Compare sum(diff(Q)) = ',sum(Q_ref - Q)
        print*,'Compare sum(Q) = ',sum(Q)
        print*,'Compare sum(Q_ref) = ',sum(Q_ref)
        print*,'***'
        print*,'Compare sum(diff(T)) = ',sum(T_ref - T)
        print*,'Compare sum(T) = ',sum(T)
        print*,'Compare sum(T_ref) = ',sum(T_ref)
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
        print*,'Compare sum(diff(SUBLC)) = ',sum(SUBLC_ref - SUBLC)
        print*,'Compare sum(SUBLC) = ',sum(SUBLC)
        print*,'Compare sum(SUBLC_ref) = ',sum(SUBLC_ref)
        print*,'***'
        print*,'Compare sum(diff(RHX)) = ',sum(RHX_ref - RHX)
        print*,'Compare sum(RHX) = ',sum(RHX)
        print*,'Compare sum(RHX_ref) = ',sum(RHX_ref)
        print*,'***'
        print*,'Compare sum(diff(EVAPC)) = ',sum(EVAPC_ref - EVAPC)
        print*,'Compare sum(EVAPC) = ',sum(EVAPC)
        print*,'Compare sum(EVAPC_ref) = ',sum(EVAPC_ref)
        print*,'***'

        ! Encoding results for CI
        print*, '#CI#VAR|PDFITERS#DIFF|',sum(PDFITERS_ref - PDFITERS)
        print*, '#CI#VAR|PDFITERS#NEW|',sum(PDFITERS)
        print*, '#CI#VAR|PDFITERS#REF|',sum(PDFITERS_ref)

        print*, '#CI#VAR|WTHV2#DIFF|',sum(WTHV2_ref - WTHV2)
        print*, '#CI#VAR|WTHV2#NEW|',sum(WTHV2)
        print*, '#CI#VAR|WTHV2#REF|',sum(WTHV2_ref)

        print*, '#CI#VAR|WQL#DIFF|',sum(WQL_ref - WQL)
        print*, '#CI#VAR|WQL#NEW|',sum(WQL)
        print*, '#CI#VAR|WQL#REF|',sum(WQL_ref)

        print*, '#CI#VAR|Q#DIFF|',sum(Q_ref - Q)
        print*, '#CI#VAR|Q#NEW|',sum(Q)
        print*, '#CI#VAR|Q#REF|',sum(Q_ref)

        print*, '#CI#VAR|T#DIFF|',sum(T_ref - T)
        print*, '#CI#VAR|T#NEW|',sum(T)
        print*, '#CI#VAR|T#REF|',sum(T_ref)

        print*, '#CI#VAR|QLLS#DIFF|',sum(QLLS_ref - QLLS)
        print*, '#CI#VAR|QLLS#NEW|',sum(QLLS)
        print*, '#CI#VAR|QLLS#REF|',sum(QLLS_ref)

        print*, '#CI#VAR|QILS#DIFF|',sum(QILS_ref - QILS)
        print*, '#CI#VAR|QILS#NEW|',sum(QILS)
        print*, '#CI#VAR|QILS#REF|',sum(QILS_ref)

        print*, '#CI#VAR|CLLS#DIFF|',sum(CLLS_ref - CLLS)
        print*, '#CI#VAR|CLLS#NEW|',sum(CLLS)
        print*, '#CI#VAR|CLLS#REF|',sum(CLLS_ref)

        print*, '#CI#VAR|QLCN#DIFF|',sum(QLCN_ref - QLCN)
        print*, '#CI#VAR|QLCN#NEW|',sum(QLCN)
        print*, '#CI#VAR|QLCN#REF|',sum(QLCN_ref)

        print*, '#CI#VAR|QICN#DIFF|',sum(QICN_ref - QICN)
        print*, '#CI#VAR|QICN#NEW|',sum(QICN)
        print*, '#CI#VAR|QICN#REF|',sum(QICN_ref)

        print*, '#CI#VAR|CLCN#DIFF|',sum(CLCN_ref - CLCN)
        print*, '#CI#VAR|CLCN#NEW|',sum(CLCN)
        print*, '#CI#VAR|CLCN#REF|',sum(CLCN_ref)

        print*, '#CI#VAR|PDF_A#DIFF|',sum(PDF_A_ref - PDF_A)
        print*, '#CI#VAR|PDF_A#NEW|',sum(PDF_A)
        print*, '#CI#VAR|PDF_A#REF|',sum(PDF_A_ref)

        print*, '#CI#VAR|SUBLC#DIFF|',sum(SUBLC_ref - SUBLC)
        print*, '#CI#VAR|SUBLC#NEW|',sum(SUBLC)
        print*, '#CI#VAR|SUBLC#REF|',sum(SUBLC_ref)

        print*, '#CI#VAR|RHX#DIFF|',sum(RHX_ref - RHX)
        print*, '#CI#VAR|RHX#NEW|',sum(RHX)
        print*, '#CI#VAR|RHX#REF|',sum(RHX_ref)

        print*, '#CI#VAR|EVAPC#DIFF|',sum(EVAPC_ref - EVAPC)
        print*, '#CI#VAR|EVAPC#NEW|',sum(EVAPC)
        print*, '#CI#VAR|EVAPC#REF|',sum(EVAPC_ref)
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