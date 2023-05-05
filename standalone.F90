program run_edmf_standalone

    use edmfparams
    use edmf_mod

    implicit none

    integer :: IM, JM, LM, NUMUP, fileID
    real    :: DT
    character*100 :: dirName, rank_str

    real, dimension(:,:),   allocatable :: PHIS, USTAR, SH, EVAP, frland, zpbl
    real, dimension(:,:,:), allocatable :: Z, ZLE, PLE, RHOE, U, V, T, THL, THV
    real, dimension(:,:,:), allocatable :: QT, Q, QL, QI
    real, dimension(:,:,:), allocatable :: ae3, aw3, aws3, awqv3, awql3, awqi3, awu3, awv3
    real, dimension(:,:,:), allocatable :: mfw2, mfw3, mfqt3, mfhl3, mfwqt, mfqt2, mfhl2, mfhlqt, mfwhl
    real, dimension(:,:,:), allocatable :: edmfdrya, edmfmoista, edmfdryw, edmfmoistw, edmfdryqt, edmfmoistqt
    real, dimension(:,:,:), allocatable :: edmfdrythl, edmfmoistthl, edmfdryu, edmfmoistu, edmfdryv, edmfmoistv
    real, dimension(:,:,:), allocatable :: edmfmoistqc, buoyf, edmf_ent, edmf_mf

    real, dimension(:,:,:), allocatable :: ae3_ref, aw3_ref, aws3_ref, awqv3_ref
    real, dimension(:,:,:), allocatable :: awql3_ref, awqi3_ref, awu3_ref, awv3_ref
    real, dimension(:,:,:), allocatable :: mfw2_ref, mfw3_ref, mfqt3_ref, mfhl3_ref, mfwqt_ref
    real, dimension(:,:,:), allocatable :: mfqt2_ref, mfhl2_ref, mfhlqt_ref, mfwhl_ref
    real, dimension(:,:,:), allocatable :: edmfdrya_ref, edmfmoista_ref, edmfdryw_ref, edmfmoistw_ref
    real, dimension(:,:,:), allocatable :: edmfdryqt_ref, edmfmoistqt_ref
    real, dimension(:,:,:), allocatable :: edmfdrythl_ref, edmfmoistthl_ref, edmfdryu_ref, edmfmoistu_ref
    real, dimension(:,:,:), allocatable :: edmfdryv_ref, edmfmoistv_ref
    real, dimension(:,:,:), allocatable :: edmfmoistqc_ref, buoyf_ref, edmf_ent_ref, edmf_mf_ref
    real :: tStart, tEnd, time

    type(EDMFPARAMS_TYPE) :: EDMFPARAMS_var

    dirName = './c90-90-90-72'
    rank_str = '1'

    ! C12 Data Set
    if (dirName == '/workdir/Turbulence/datasets/run_edmf/c12-6-6-72') then
        IM = 6
        JM = 6
        LM = 72
        DT = 900.0
    else if (dirName == './c90-90-90-72') then
        IM = 90
        JM = 90
        LM = 72
        DT = 450.0
    else if (dirName == '/workdir/Turbulence/datasets/run_edmf/c180-180-180-72') then
        IM = 180
        JM = 180
        LM = 72
        DT = 450.0
    endif

    
    NUMUP = 10

    EDMFPARAMS_var%DISCRETE = 0
    EDMFPARAMS_var%IMPLICIT = 1
    EDMFPARAMS_var%ENTRAIN  = 1
    EDMFPARAMS_var%DOCLASP  = 0
    EDMFPARAMS_var%THERMAL_PLUME = 0
    EDMFPARAMS_var%TEST = 0
    EDMFPARAMS_var%DEBUG = 0
    EDMFPARAMS_var%ET = 2
    EDMFPARAMS_var%L0 = 100.0
    EDMFPARAMS_var%L0fac = 10.0
    EDMFPARAMS_var%STOCHFRAC = 0.50
    EDMFPARAMS_var%ENTWFAC = 0.3333
    EDMFPARAMS_var%EDFAC = 1.0
    EDMFPARAMS_var%ENT0 = 1.0
    EDMFPARAMS_var%ALPHATH = 2.89
    EDMFPARAMS_var%ALPHAQT = 2.89
    EDMFPARAMS_var%ALPHAW = 0.572
    EDMFPARAMS_var%PWMAX = 3.0
    EDMFPARAMS_var%PWMIN = 1.0
    EDMFPARAMS_var%WA = 0.0
    EDMFPARAMS_var%WB = 0.0
    EDMFPARAMS_var%AU0 = 0.0
    EDMFPARAMS_var%CTH1 = 0.0
    EDMFPARAMS_var%CTH2 = 0.0
    EDMFPARAMS_var%RH0_QB = 0.0
    EDMFPARAMS_var%C_KH_MF = 0.0
    EDMFPARAMS_var%ICE_RAMP = -40.0

!   Inputs arrays
    allocate(PHIS(IM, JM))
    allocate(Z   (IM, JM, LM))
    allocate(ZLE (IM, JM, 0:LM))
    allocate(PLE (IM, JM, 0:LM))
    allocate(RHOE(IM, JM, 0:LM))
    allocate(U   (IM, JM, LM))
    allocate(V   (IM, JM, LM))
    allocate(T   (IM, JM, LM))
    allocate(THL (IM, JM, LM))
    allocate(THV (IM, JM, LM))
    allocate(QT  (IM, JM, LM))
    allocate(Q   (IM, JM, LM))
    allocate(QL  (IM, JM, LM))
    allocate(QI  (IM, JM, LM))
    allocate(USTAR(IM, JM))
    allocate(SH   (IM, JM))
    allocate(EVAP (IM, JM))
    allocate(frland(IM, JM))
    allocate(zpbl  (IM, JM))

!   Output arrays
    allocate(ae3  (IM, JM, 0:LM))
    allocate(aw3  (IM, JM, 0:LM))
    allocate(aws3 (IM, JM, 0:LM))
    allocate(awqv3(IM, JM, 0:LM))
    allocate(awql3(IM, JM, 0:LM))
    allocate(awqi3(IM, JM, 0:LM))
    allocate(awu3 (IM, JM, 0:LM))
    allocate(awv3 (IM, JM, 0:LM))

    allocate(mfw2 (IM, JM, LM))
    allocate(mfw3 (IM, JM, LM))
    allocate(mfqt3(IM, JM, LM))
    allocate(mfhl3(IM, JM, LM))
    allocate(mfwqt(IM, JM, LM))
    allocate(mfqt2(IM, JM, LM))
    allocate(mfhl2(IM, JM, LM))
    allocate(mfhlqt(IM, JM, LM))
    allocate(mfwhl (IM, JM, LM))

    allocate(edmfdrya    (IM, JM, 0:LM))
    allocate(edmfmoista  (IM, JM, 0:LM))
    allocate(edmfdryw    (IM, JM, 0:LM))
    allocate(edmfmoistw  (IM, JM, 0:LM))
    allocate(edmfdryqt   (IM, JM, 0:LM))
    allocate(edmfmoistqt (IM, JM, 0:LM))
    allocate(edmfdrythl  (IM, JM, 0:LM))
    allocate(edmfmoistthl(IM, JM, 0:LM))
    allocate(edmfdryu    (IM, JM, 0:LM))
    allocate(edmfmoistu  (IM, JM, 0:LM))
    allocate(edmfdryv    (IM, JM, 0:LM))
    allocate(edmfmoistv  (IM, JM, 0:LM))
    allocate(edmfmoistqc (IM, JM, 0:LM))
    allocate(buoyf       (IM, JM, LM))
    allocate(edmf_ent    (IM, JM, LM))
    allocate(edmf_mf     (IM, JM, 0:LM))

!   Output Reference Arrays
    allocate(ae3_ref  (IM, JM, 0:LM))
    allocate(aw3_ref  (IM, JM, 0:LM))
    allocate(aws3_ref (IM, JM, 0:LM))
    allocate(awqv3_ref(IM, JM, 0:LM))
    allocate(awql3_ref(IM, JM, 0:LM))
    allocate(awqi3_ref(IM, JM, 0:LM))
    allocate(awu3_ref (IM, JM, 0:LM))
    allocate(awv3_ref (IM, JM, 0:LM))

    allocate(mfw2_ref (IM, JM, LM))
    allocate(mfw3_ref (IM, JM, LM))
    allocate(mfqt3_ref(IM, JM, LM))
    allocate(mfhl3_ref(IM, JM, LM))
    allocate(mfwqt_ref(IM, JM, LM))
    allocate(mfqt2_ref(IM, JM, LM))
    allocate(mfhl2_ref(IM, JM, LM))
    allocate(mfhlqt_ref(IM, JM, LM))
    allocate(mfwhl_ref (IM, JM, LM))

    allocate(edmfdrya_ref    (IM, JM, 0:LM))
    allocate(edmfmoista_ref  (IM, JM, 0:LM))
    allocate(edmfdryw_ref    (IM, JM, 0:LM))
    allocate(edmfmoistw_ref  (IM, JM, 0:LM))
    allocate(edmfdryqt_ref   (IM, JM, 0:LM))
    allocate(edmfmoistqt_ref (IM, JM, 0:LM))
    allocate(edmfdrythl_ref  (IM, JM, 0:LM))
    allocate(edmfmoistthl_ref(IM, JM, 0:LM))
    allocate(edmfdryu_ref    (IM, JM, 0:LM))
    allocate(edmfmoistu_ref  (IM, JM, 0:LM))
    allocate(edmfdryv_ref    (IM, JM, 0:LM))
    allocate(edmfmoistv_ref  (IM, JM, 0:LM))
    allocate(edmfmoistqc_ref (IM, JM, 0:LM))
    allocate(buoyf_ref       (IM, JM, LM))
    allocate(edmf_ent_ref    (IM, JM, LM))
    allocate(edmf_mf_ref     (IM, JM, 0:LM))

!   Reading Input Arrays

    open(newunit=fileID, file=trim(dirName) // "/PHIS_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) PHIS
    close(fileID)
    !write(*,*) 'sum(PHIS) = ', sum(PHIS)

    open(newunit=fileID, file=trim(dirName) // "/Z_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) Z
    close(fileID)
    !write(*,*) 'sum(Z) = ', sum(Z)

    open(newunit=fileID, file=trim(dirName) // "/ZLE_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) ZLE
    close(fileID)
    !write(*,*) 'sum(ZLE) = ', sum(ZLE)

    open(newunit=fileID, file=trim(dirName) // "/PLE_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) PLE
    close(fileID)
    !write(*,*) 'sum(PLE) = ', sum(PLE)

    open(newunit=fileID, file=trim(dirName) // "/RHOE_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) RHOE
    close(fileID)
    !write(*,*) 'sum(RHOE) = ', sum(RHOE)

    open(newunit=fileID, file=trim(dirName) // "/U_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) U
    close(fileID)
    !write(*,*) 'sum(U) = ', sum(U)

    open(newunit=fileID, file=trim(dirName) // "/V_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) V
    close(fileID)
    !write(*,*) 'sum(V) = ', sum(V)

    open(newunit=fileID, file=trim(dirName) // "/T_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) T
    close(fileID)
    !write(*,*) 'sum(T) = ', sum(T)

    open(newunit=fileID, file=trim(dirName) // "/THL_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) THL
    close(fileID)
    !write(*,*) 'sum(THL) = ', sum(THL)

    open(newunit=fileID, file=trim(dirName) // "/THV_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) THV
    close(fileID)
    !write(*,*) 'sum(THV) = ', sum(THV)

    open(newunit=fileID, file=trim(dirName) // "/QT_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QT
    close(fileID)
    !write(*,*) 'sum(QT) = ', sum(QT)

    open(newunit=fileID, file=trim(dirName) // "/Q_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) Q
    close(fileID)
    !write(*,*) 'sum(Q) = ', sum(Q)

    open(newunit=fileID, file=trim(dirName) // "/QL_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QL
    close(fileID)
    !write(*,*) 'sum(QL) = ', sum(QL)

    open(newunit=fileID, file=trim(dirName) // "/QI_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QI
    close(fileID)
    !write(*,*) 'sum(QI) = ', sum(QI)

    open(newunit=fileID, file=trim(dirName) // "/USTAR_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) USTAR
    close(fileID)
    !write(*,*) 'sum(USTAR) = ', sum(USTAR)

    open(newunit=fileID, file=trim(dirName) // "/SH_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) SH
    close(fileID)
    !write(*,*) 'sum(SH) = ', sum(SH)

    open(newunit=fileID, file=trim(dirName) // "/EVAP_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) EVAP
    close(fileID)
    !write(*,*) 'sum(EVAP) = ', sum(EVAP)

    open(newunit=fileID, file=trim(dirName) // "/FRLAND_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) frland
    close(fileID)
    !write(*,*) 'sum(frland) = ', sum(frland)

    open(newunit=fileID, file=trim(dirName) // "/ZPBL_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) zpbl
    close(fileID)
    !write(*,*) 'sum(zpbl) = ', sum(zpbl)

!   Reading output reference solution arrays

    open(newunit=fileID, file=trim(dirName) // "/ae3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) ae3_ref
    close(fileID)
    !write(*,*) 'sum(ae3) = ', sum(ae3_ref)

    open(newunit=fileID, file=trim(dirName) // "/aw3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) aw3_ref
    close(fileID)
    !write(*,*) 'sum(aw3) = ', sum(aw3_ref)

    open(newunit=fileID, file=trim(dirName) // "/aws3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) aws3_ref
    close(fileID)
    !write(*,*) 'sum(aws3) = ', sum(aws3_ref)

    open(newunit=fileID, file=trim(dirName) // "/awqv3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) awqv3_ref
    close(fileID)
    !write(*,*) 'sum(awqv3) = ', sum(awqv3_ref)

    open(newunit=fileID, file=trim(dirName) // "/awql3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) awql3_ref
    close(fileID)
    !write(*,*) 'sum(awql3) = ', sum(awql3_ref)

    open(newunit=fileID, file=trim(dirName) // "/awqi3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) awqi3_ref
    close(fileID)
    !write(*,*) 'sum(awqi3) = ', sum(awqi3_ref)

    open(newunit=fileID, file=trim(dirName) // "/awu3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) awu3_ref
    close(fileID)
    !write(*,*) 'sum(awu3) = ', sum(awu3_ref)

    open(newunit=fileID, file=trim(dirName) // "/awv3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) awv3_ref
    close(fileID)
    !write(*,*) 'sum(awv3) = ', sum(awv3_ref)

    open(newunit=fileID, file=trim(dirName) // "/mfw2_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) mfw2_ref
    close(fileID)
    !write(*,*) 'sum(mfw2) = ', sum(mfw2_ref)

    open(newunit=fileID, file=trim(dirName) // "/mfw3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) mfw3_ref
    close(fileID)
    !write(*,*) 'sum(mfw3) = ', sum(mfw3_ref)

    open(newunit=fileID, file=trim(dirName) // "/mfqt3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) mfqt3_ref
    close(fileID)
    !write(*,*) 'sum(mfqt3) = ', sum(mfqt3_ref)

    open(newunit=fileID, file=trim(dirName) // "/mfhl3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) mfhl3_ref
    close(fileID)
    !write(*,*) 'sum(mfhl3) = ', sum(mfhl3_ref)

    open(newunit=fileID, file=trim(dirName) // "/mfwqt_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) mfwqt_ref
    close(fileID)
    !write(*,*) 'sum(mfwqt) = ', sum(mfwqt_ref)

    open(newunit=fileID, file=trim(dirName) // "/mfqt2_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) mfqt2_ref
    close(fileID)
    !write(*,*) 'sum(mfqt2) = ', sum(mfqt2_ref)

    open(newunit=fileID, file=trim(dirName) // "/mfhl2_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) mfhl2_ref
    close(fileID)
    !write(*,*) 'sum(mfhl2) = ', sum(mfhl2_ref)

    open(newunit=fileID, file=trim(dirName) // "/mfhlqt_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) mfhlqt_ref
    close(fileID)
    !write(*,*) 'sum(mfhlqt) = ', sum(mfhlqt_ref)

    open(newunit=fileID, file=trim(dirName) // "/mfwhl_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) mfwhl_ref
    close(fileID)
    !write(*,*) 'sum(mfwhl) = ', sum(mfwhl_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfdrya_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfdrya_ref
    close(fileID)
    !write(*,*) 'sum(edmfdrya) = ', sum(edmfdrya_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfmoista_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfmoista_ref
    close(fileID)
    !write(*,*) 'sum(edmfmoista) = ', sum(edmfmoista_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfdryw_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfdryw_ref
    close(fileID)
    !write(*,*) 'sum(edmfdryw) = ', sum(edmfdryw_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfmoistw_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfmoistw_ref
    close(fileID)
    !write(*,*) 'sum(edmfmoistw) = ', sum(edmfmoistw_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfdryqt_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfdryqt_ref
    close(fileID)
    !write(*,*) 'sum(edmfdryqt) = ', sum(edmfdryqt_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfmoistqt_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfmoistqt_ref
    close(fileID)
    !write(*,*) 'sum(edmfmoistqt) = ', sum(edmfmoistqt_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfdrythl_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfdrythl_ref
    close(fileID)
    !write(*,*) 'sum(edmfdrythl) = ', sum(edmfdrythl_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfmoistthl_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfmoistthl_ref
    close(fileID)
    !write(*,*) 'sum(edmfmoistthl) = ', sum(edmfmoistthl_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfdryu_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfdryu_ref
    close(fileID)
    !write(*,*) 'sum(edmfdryu) = ', sum(edmfdryu_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfmoistu_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfmoistu_ref
    close(fileID)
    !write(*,*) 'sum(edmfmoistu) = ', sum(edmfmoistu_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfdryv_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfdryv_ref
    close(fileID)
    !write(*,*) 'sum(edmfdryv) = ', sum(edmfdryv_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfmoistv_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfmoistv_ref
    close(fileID)
    !write(*,*) 'sum(edmfmoistv) = ', sum(edmfmoistv_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmfmoistqc_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmfmoistqc_ref
    close(fileID)
    !write(*,*) 'sum(edmfmoistqc) = ', sum(edmfmoistqc_ref)

    open(newunit=fileID, file=trim(dirName) // "/buoyf_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) buoyf_ref
    close(fileID)
    !write(*,*) 'sum(buoyf) = ', sum(buoyf_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmf_ent_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmf_ent_ref
    close(fileID)
    !write(*,*) 'sum(edmf_ent) = ', sum(edmf_ent_ref)

    open(newunit=fileID, file=trim(dirName) // "/edmf_mf_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) edmf_mf_ref
    close(fileID)
    !write(*,*) 'sum(edmf_mf) = ', sum(edmf_mf_ref)

!$acc data copyin(phis,z,zle,ple,rhoe,u,v,t,thl,thv,qt,q,ql,qi,ustar,sh,evap,frland,zpbl,&
!$acc             EDMFPARAMS_var,EDMFPARAMS_var%edfac,EDMFPARAMS_var%doclasp,EDMFPARAMS_var%et,&
!$acc             EDMFPARAMS_var%l0fac,EDMFPARAMS_var%l0,EDMFPARAMS_var%discrete,EDMFPARAMS_var%stochfrac,&
!$acc             EDMFPARAMS_var%ent0,EDMFPARAMS_var%alphaw,EDMFPARAMS_var%alphaqt,EDMFPARAMS_var%alphath,&
!$acc             EDMFPARAMS_var%pwmin,EDMFPARAMS_var%pwmax,EDMFPARAMS_var%ice_ramp,EDMFPARAMS_var%entwfac,&
!$acc             EDMFPARAMS_var%entrain,EDMFPARAMS_var%implicit) &
!$acc      copyout(ae3,aw3,aws3,awqv3,awql3,awqi3,awu3,awv3,mfw2,mfw3,mfqt3,mfhl3,mfwqt,&
!$acc              mfqt2,mfhl2,mfhlqt,mfwhl,edmfdrya,edmfmoista,edmfdryw,edmfmoistw,&
!$acc              edmfdryqt,edmfmoistqt,edmfdrythl,edmfmoistthl,edmfdryu,edmfmoistu,&
!$acc              edmfdryv,edmfmoistv,edmfmoistqc,buoyf,edmf_ent,edmf_mf)

    call cpu_time(tStart)

    call RUN_EDMF(1, IM*JM, 1, LM, DT,      & ! in
               PHIS, Z, ZLE, PLE, RHOE,       & ! in
               NUMUP, U, V, T, THL, THV, QT,  & ! in
               Q, QL, QI, USTAR,              & ! in
               SH, EVAP, frland, zpbl,        & ! in
               ae3, aw3, aws3, awqv3,         & ! for trisolver 
               awql3, awqi3, awu3, awv3,      & ! for trisolver
               mfw2,mfw3,mfqt3,mfhl3,mfwqt,   & ! for ADG PDF
               mfqt2, mfhl2, mfhlqt, mfwhl,   & ! for ADG PDF
               edmfdrya, edmfmoista,          & ! diag
               edmfdryw, edmfmoistw,          & ! diag
               edmfdryqt, edmfmoistqt,        & ! diag
               edmfdrythl, edmfmoistthl,      & ! diag
               edmfdryu, edmfmoistu,          & ! diag
               edmfdryv, edmfmoistv,          & ! diag
               edmfmoistqc,                   & ! diag
               buoyf, edmf_ent, edmf_mf,      & ! diag  
               EDMFPARAMS_var )

    call cpu_time(tEnd)
    time = tEnd - tStart
    print*, 'Total Time = ', time

!$acc end data

    write(*,*) 'Sum diff abs ae3 = ', sum(abs(ae3) - abs(ae3_ref)), sum(abs(ae3)), sum(abs(ae3_ref))
    write(*,*) 'Sum diff abs aw3 = ', sum(abs(aw3) - abs(aw3_ref)), sum(abs(aw3)), sum(abs(aw3_ref))
    write(*,*) 'Sum diff abs aws3 = ', sum(abs(aws3) - abs(aws3_ref)), sum(abs(aws3)),  sum(abs(aws3_ref))
    write(*,*) 'Sum diff abs awqv3 = ', sum(abs(awqv3) - abs(awqv3_ref)), sum(abs(awqv3)), sum(abs(awqv3_ref))
    write(*,*) 'Sum diff abs awql3 = ', sum(abs(awql3) - abs(awql3_ref)), sum(abs(awql3)), sum(abs(awql3_ref))
    write(*,*) 'Sum diff abs awqi3 = ', sum(abs(awqi3) - abs(awqi3_ref)), sum(abs(awqi3)), sum(abs(awqi3_ref))
    write(*,*) 'Sum diff abs awu3 = ', sum(abs(awu3) - abs(awu3_ref)), sum(abs(awu3)), sum(abs(awu3_ref))
    write(*,*) 'Sum diff abs awv3 = ', sum(abs(awv3) - abs(awv3_ref)), sum(abs(awv3)), sum(abs(awv3_ref))

    write(*,*) 'Sum diff abs mfw2 = ', sum(abs(mfw2) - abs(mfw2_ref)), sum(abs(mfw2)), sum(abs(mfw2_ref))
    write(*,*) 'Sum diff abs mfw3 = ', sum(abs(mfw3) - abs(mfw3_ref)), sum(abs(mfw3)), sum(abs(mfw3_ref))
    write(*,*) 'Sum diff abs mfqt3 = ', sum(abs(mfqt3) - abs(mfqt3_ref)), sum(abs(mfqt3)), sum((mfqt3_ref))
    write(*,*) 'Sum diff abs mfhl3 = ', sum(abs(mfhl3) - abs(mfhl3_ref)), sum(abs(mfhl3)), sum((mfhl3_ref))
    write(*,*) 'Sum diff abs mfwqt = ', sum(abs(mfwqt) - abs(mfwqt_ref)), sum(abs(mfwqt)), sum(abs(mfwqt_ref))
    write(*,*) 'Sum diff abs mfqt2 = ', sum(abs(mfqt2) - abs(mfqt2_ref)), sum(abs(mfqt2)), sum(abs(mfqt2_ref))
    write(*,*) 'Sum diff abs mfhl2 = ', sum(abs(mfhl2) - abs(mfhl2_ref)), sum(abs(mfhl2)), sum(abs(mfhl2_ref))
    write(*,*) 'Sum diff abs mfhlqt = ', sum(abs(mfhlqt) - abs(mfhlqt_ref)), sum(abs(mfhlqt)), sum(abs(mfhlqt_ref))
    write(*,*) 'Sum diff abs mfwhl = ', sum(abs(mfwhl) - abs(mfwhl_ref)), sum(abs(mfwhl)), sum(abs(mfwhl_ref))

    write(*,*) 'Sum diff abs edmfdrya = ', sum(abs(edmfdrya) - abs(edmfdrya_ref)), sum(abs(edmfdrya)), sum(abs(edmfdrya_ref))
    write(*,*) 'Sum diff abs edmfmoista = ', sum(abs(edmfmoista) - abs(edmfmoista_ref)), sum(abs(edmfmoista)), sum(abs(edmfmoista_ref))
    write(*,*) 'Sum diff abs edmfdryw= ', sum(abs(edmfdryw) - abs(edmfdryw_ref)), sum(abs(edmfdryw)), sum(abs(edmfdryw_ref))
    write(*,*) 'Sum diff abs edmfmoistw = ', sum(abs(edmfmoistw) - abs(edmfmoistw_ref)), sum(abs(edmfmoistw)), sum(abs(edmfmoistw_ref))
    write(*,*) 'Sum diff abs edmfdryqt = ', sum(abs(edmfdryqt) - abs(edmfdryqt_ref)), sum(abs(edmfdryqt)), sum(abs(edmfdryqt_ref))
    write(*,*) 'Sum diff abs edmfmoistqt = ', sum(abs(edmfmoistqt) - abs(edmfmoistqt_ref)), sum(abs(edmfmoistqt)), sum(abs(edmfmoistqt_ref))
    write(*,*) 'Sum diff abs edmfdrythl = ', sum(abs(edmfdrythl) - abs(edmfdrythl_ref)), sum(abs(edmfdrythl)), sum(abs(edmfdrythl_ref))
    write(*,*) 'Sum diff abs edmfmoistthl = ', sum(abs(edmfmoistthl) - abs(edmfmoistthl_ref)), sum(abs(edmfmoistthl)), sum(abs(edmfmoistthl_ref))
    write(*,*) 'Sum diff abs edmfdryu = ', sum(abs(edmfdryu) - abs(edmfdryu_ref)), sum(abs(edmfdryu)), sum(abs(edmfdryu_ref))
    write(*,*) 'Sum diff abs edmfmoistu = ', sum(abs(edmfmoistu) - abs(edmfmoistu_ref)), sum(abs(edmfmoistu)), sum(abs(edmfmoistu_ref))
    write(*,*) 'Sum diff abs edmfdryv = ', sum(abs(edmfdryv) - abs(edmfdryv_ref)), sum(abs(edmfdryv)), sum(abs(edmfdryv_ref))
    write(*,*) 'Sum diff abs edmfmoistv = ', sum(abs(edmfmoistv) - abs(edmfmoistv_ref)), sum(abs(edmfmoistv)), sum(abs(edmfmoistv_ref))
    write(*,*) 'Sum diff abs edmfmoistqc = ', sum(abs(edmfmoistqc) - abs(edmfmoistqc_ref)), sum(abs(edmfmoistqc)), sum(abs(edmfmoistqc_ref))
    write(*,*) 'Sum diff abs buoyf = ', sum(abs(buoyf) - abs(buoyf_ref)), sum(abs(buoyf)), sum(abs(buoyf_ref))
    write(*,*) 'Sum diff abs edmf_ent = ', sum(abs(edmf_ent) - abs(edmf_ent_ref)), sum(abs(edmf_ent)), sum(abs(edmf_ent_ref))
    write(*,*) 'Sum diff abs edmf_mf = ', sum(abs(edmf_mf) - abs(edmf_mf_ref)), sum(abs(edmf_mf)), sum(abs(edmf_mf_ref))

end program

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
