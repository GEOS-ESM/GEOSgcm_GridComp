program postlock_standalone

    use MAPL_ConstantsMod
    use EDMFPARAMS

    implicit none

    integer :: IM, JM, LM, I, J, L, KPBLMIN, PBLHT_OPTION, DO_SHOC, locmax, fileID
    real    :: PDFSHAPE, dthvdz, maxdthvdz, DT, maxkh
    real    :: minlval, LAMBDADISS
    logical :: CALC_TCZPBL, CALC_ZPBL2, CALC_ZPBL10p, file_exists

    character*100 :: dirName, rank_str

    real, parameter :: tcri_crit = 0.25
    real, parameter :: ri_crit = 0.00
    real, parameter :: ri_crit2 = 0.20
    real, parameter :: mapl_undef = 1.0e15  ! NOTE : This is the value pulled from MAPL_Mod

    ! *** Arrays for postlock **

    ! *** Input Arrays ***
    real, dimension(:),     allocatable :: PREF

    real, dimension(:,:),   allocatable :: FRLAND, CT, CQ, CU

    real, dimension(:,:,:), allocatable :: TKESHOC, KH, THV, Z, KM, U, V
    real, dimension(:,:,:), allocatable :: QLLS, T, Q, TH, KHSFC, ZLE, RI, PLO, RDZ
    real, dimension(:,:,:), allocatable :: DMI, PLE, TV, AE3, RHOE, AWS3, AWQV3, AWQL3, AWQI3
    real, dimension(:,:,:), allocatable :: AWU3, AWV3, HL, AW3

    ! *** Output Arrays ***

    real, dimension(:),     allocatable :: temparray, thetav
    real, dimension(:),     allocatable :: temparray_ref, thetav_ref

    real, dimension(:,:),   allocatable :: ZPBL,  thetavs, thetavh, uv2h, KPBL, KPBLTC
    real, dimension(:,:),   allocatable :: ZPBL2, KPBL2, KPBL_SC, ZPBL10p, KPBL10p, TCZPBL

    real, dimension(:,:),   allocatable :: ZPBL_ref,  thetavs_ref, thetavh_ref, uv2h_ref, KPBL_ref, KPBLTC_ref
    real, dimension(:,:),   allocatable :: ZPBL2_ref, KPBL2_ref, KPBL_SC_ref, ZPBL10p_ref, KPBL10p_ref, TCZPBL_ref

    real, dimension(:,:,:), allocatable :: ISOTROPY, tcrib
    real, dimension(:,:,:), allocatable :: CKS, AKS, TKH, AKQ, CKQ, EKV, AKV
    real, dimension(:,:,:), allocatable :: CKV, BKS, BKQ, BKV, AKSS, AKUU, RHOAW3
    real, dimension(:,:,:), allocatable :: AKQQ, CKSS, CKQQ, CKUU, BKSS, BKQQ, BKUU, YS
    real, dimension(:,:,:), allocatable :: YQV, YQL, YQI, YU, YV

    real, dimension(:,:,:), allocatable :: ISOTROPY_ref, tcrib_ref
    real, dimension(:,:,:), allocatable :: CKS_ref, AKS_ref, TKH_ref, AKQ_ref, CKQ_ref, EKV_ref, AKV_ref
    real, dimension(:,:,:), allocatable :: CKV_ref, BKS_ref, BKQ_ref, BKV_ref, AKSS_ref, AKUU_ref, RHOAW3_ref
    real, dimension(:,:,:), allocatable :: AKQQ_ref, CKSS_ref, CKQQ_ref, CKUU_ref, BKSS_ref, BKQQ_ref, BKUU_ref, YS_ref
    real, dimension(:,:,:), allocatable :: YQV_ref, YQL_ref, YQI_ref, YU_ref, YV_ref
    

    real, dimension(:,:),   allocatable :: PPBL_ref, ZPBLRI_ref, ZPBLRI2_ref, ZPBLTHV_ref, ZPBLHTKE_ref
    real, dimension(:,:),   pointer :: PPBL, ZPBLHTKE, ZPBLRI, ZPBLRI2, ZPBLTHV

    real, dimension(:,:,:), allocatable :: TKE_ref, AKSODT_ref, CKSODT_ref, AKQODT_ref, CKQODT_ref, AKVODT_ref, CKVODT_ref
    real, dimension(:,:,:), pointer :: TKE, AKSODT, CKSODT, AKQODT, CKQODT, AKVODT, CKVODT

    type (EDMFPARAMS_TYPE) :: EDMFPARAMS_

    real :: t1, t2

    dirName = './c90-90-90-72'
    rank_str = '0'

    ! C12 Data Set
    if (dirName == './c12-6-6-72') then
        IM = 6
        JM = 6
        LM = 72
        DT = 900.0
    else if (dirName == './c90-90-90-72') then
        IM = 90
        JM = 90
        LM = 72
        DT = 450.0
    else if(dirName == './c180-180-180-72') then
        IM = 180
        JM = 180
        LM = 72
        DT = 450.0
    endif

    DO_SHOC = 1
    LAMBDADISS = 50.0
    PDFSHAPE = 5.0
    CALC_TCZPBL = .true.
    CALC_ZPBL2 = .false.
    CALC_ZPBL10p = .true.
    PBLHT_OPTION = 4
    EDMFPARAMS_%IMPLICIT = 1
    EDMFPARAMS_%DISCRETE = 0

    PPBL     => null()
    ZPBLHTKE => null()
    ZPBLRI   => null()
    ZPBLRI2  => null()
    ZPBLTHV  => null()
    TKE      => null()
    AKSODT   => null()
    CKSODT   => null()
    AKQODT   => null()
    CKQODT   => null()
    AKVODT   => null()
    CKVODT   => null()

    ! *** Input Arrays ***
    allocate(TKESHOC(IM, JM, LM))
    allocate(KH     (IM, JM, 0:LM))
    allocate(THV    (IM, JM, LM))
    allocate(Z      (IM, JM, LM))
    allocate(KM     (IM, JM, 0:LM))
    allocate(U      (IM, JM, LM))
    allocate(V      (IM, JM, LM))
    allocate(QLLS   (IM, JM, LM))
    allocate(PREF   (0:LM))
    allocate(T      (IM, JM, LM))
    allocate(Q      (IM, JM, LM))
    allocate(TH     (IM, JM, LM))
    allocate(KHSFC  (IM, JM, 0:LM))
    allocate(ZLE    (IM, JM, 0:LM))
    allocate(RI     (IM, JM, 0:LM))
    allocate(FRLAND (IM, JM))
    allocate(PLO    (IM, JM, LM))
    allocate(RDZ    (IM, JM, LM-1))
    allocate(DMI    (IM, JM, LM))
    allocate(CT     (IM, JM))
    allocate(PLE    (IM, JM, 0:LM))
    allocate(TV     (IM, JM, LM))
    allocate(CQ     (IM, JM))
    allocate(CU     (IM, JM))
    allocate(AE3    (IM, JM, 0:LM))
    allocate(RHOE   (IM, JM, 0:LM))
    allocate(AWS3   (IM, JM, 0:LM))
    allocate(AWQV3  (IM, JM, 0:LM))
    allocate(AWQL3  (IM, JM, 0:LM))
    allocate(AWQI3  (IM, JM, 0:LM))
    allocate(AWU3   (IM, JM, 0:LM))
    allocate(AWV3   (IM, JM, 0:LM))
    allocate(AW3    (IM, JM, 0:LM))
    allocate(TKH(IM, JM, 0:LM))
    allocate(tcrib(IM, JM, LM))

    INQUIRE(FILE=trim(dirName) // "/TKESHOC_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/TKESHOC_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) TKESHOC
    close(fileID)
    !write(*,*) 'sum(TKESHOC) = ', sum((TKESHOC)), size(TKESHOC,1), size(TKESHOC,2), size(TKESHOC,3)

    INQUIRE(FILE=trim(dirName) // "/KH_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/KH_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) KH
    close(fileID)
    ! write(*,*) 'sum(KH) = ', sum((KH)), size(KH,1), size(KH,2), size(KH,3)

    INQUIRE(FILE=trim(dirName) // "/THV_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/THV_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) THV
    close(fileID)
    !write(*,*) 'sum(THV) = ', sum((THV)), size(THV,1), size(THV,2), size(THV,3)

    INQUIRE(FILE=trim(dirName) // "/Z_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/Z_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) Z
    close(fileID)
    ! write(*,*) 'sum(Z) = ', sum((Z)), size(Z,1), size(Z,2), size(Z,3)

    INQUIRE(FILE=trim(dirName) // "/KM_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/KM_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) KM
    close(fileID)
    !write(*,*) 'sum(KM) = ', sum((KM)), size(KM,1), size(KM,2), size(KM,3)

    INQUIRE(FILE=trim(dirName) // "/U_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/U_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) U
    close(fileID)
    ! write(*,*) 'sum(U) = ', sum((U)), size(U,1), size(U,2), size(U,3)
     
    INQUIRE(FILE=trim(dirName) // "/V_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/V_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) V
    close(fileID)
    ! write(*,*) 'sum(V) = ', sum((V)), size(V,1), size(V,2), size(V,3)

    INQUIRE(FILE=trim(dirName) // "/QLLS_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/QLLS_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QLLS
    close(fileID)
    !write(*,*) 'sum(QLLS) = ', sum((QLLS)), size(QLLS,1), size(QLLS,2), size(QLLS,3)

    INQUIRE(FILE=trim(dirName) // "/PREF_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/PREF_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) PREF
    close(fileID)
    !write(*,*) 'sum(PREF) = ', sum((PREF)), size(PREF,1)
    
    INQUIRE(FILE=trim(dirName) // "/T_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/T_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) T
    close(fileID)
    !write(*,*) 'sum(T) = ', sum((T)), size(T,1), size(T,2), size(T,3)

    INQUIRE(FILE=trim(dirName) // "/Q_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/Q_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) Q
    close(fileID)
    !write(*,*) 'sum(Q) = ', sum((Q)), size(Q,1), size(Q,2), size(Q,3)

    INQUIRE(FILE=trim(dirName) // "/TH_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/TH_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) TH
    close(fileID)
    !write(*,*) 'sum(TH) = ', sum((TH)), size(TH,1), size(TH,2), size(TH,3)

    INQUIRE(FILE=trim(dirName) // "/KHSFC_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/KHSFC_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) KHSFC
    close(fileID)
    !write(*,*) 'sum(KHSFC) = ', sum((KHSFC)), size(KHSFC,1), size(KHSFC,2), size(KHSFC,3)

    INQUIRE(FILE=trim(dirName) // "/ZLE_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/ZLE_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) ZLE
    close(fileID)
    !write(*,*) 'sum(ZLE) = ', sum((ZLE)), size(ZLE,1), size(ZLE,2), size(ZLE,3)

    INQUIRE(FILE=trim(dirName) // "/RI_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/RI_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) RI
    close(fileID)
    !write(*,*) 'sum(RI) = ', sum((RI)), size(RI,1), size(RI,2), size(RI,3)

    INQUIRE(FILE=trim(dirName) // "/FRLAND_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/FRLAND_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) FRLAND
    close(fileID)
    !write(*,*) 'sum(FRLAND) = ', sum((FRLAND)), size(FRLAND,1), size(FRLAND,2)

    INQUIRE(FILE=trim(dirName) // "/PLO_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/PLO_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) PLO
    close(fileID)
    !write(*,*) 'sum(PLO) = ', sum((PLO)), size(PLO,1), size(PLO,2), size(PLO,3)

    INQUIRE(FILE=trim(dirName) // "/RDZ_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/RDZ_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) RDZ
    close(fileID)
    !write(*,*) 'sum(RDZ) = ', sum((RDZ)), size(RDZ,1), size(RDZ,2), size(RDZ,3)

    INQUIRE(FILE=trim(dirName) // "/DMI_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/DMI_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) DMI
    close(fileID)
    !write(*,*) 'sum(DMI)) = ', sum((DMI)), size(DMI,1), size(DMI,2), size(DMI,3)

    INQUIRE(FILE=trim(dirName) // "/CT_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/CT_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) CT
    close(fileID)
    !write(*,*) 'sum(CT) = ', sum((CT)), size(CT,1), size(CT,2)

    INQUIRE(FILE=trim(dirName) // "/PLE_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/PLE_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) PLE
    close(fileID)
    !write(*,*) 'sum(PLE) = ', sum((PLE)), size(PLE,1), size(PLE,2), size(PLE,3)

    INQUIRE(FILE=trim(dirName) // "/TV_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/TV_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) TV
    close(fileID)
    !write(*,*) 'sum(TV) = ', sum((TV)), size(TV,1), size(TV,2), size(TV,3)

    INQUIRE(FILE=trim(dirName) // "/CQ_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/CQ_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) CQ
    close(fileID)
    !write(*,*) 'sum(CQ) = ', sum((CQ)), size(CQ,1), size(CQ,2)

    INQUIRE(FILE=trim(dirName) // "/CU_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/CU_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) CU
    close(fileID)
    !write(*,*) 'sum(CU) = ', sum((CU)), size(CU,1), size(CU,2)

    INQUIRE(FILE=trim(dirName) // "/AE3_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/AE3_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) AE3
    close(fileID)
    !write(*,*) 'sum(AE3) = ', sum((AE3)), size(AE3,1), size(AE3,2), size(AE3,3)

    INQUIRE(FILE=trim(dirName) // "/RHOE_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/RHOE_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) RHOE
    close(fileID)
    !write(*,*) 'sum(RHOE) = ', sum((RHOE)), size(RHOE,1), size(RHOE,2), size(RHOE,3)

    INQUIRE(FILE=trim(dirName) // "/AWS3_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/AWS3_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) AWS3
    close(fileID)
    !write(*,*) 'sum(AWS3) = ', sum((AWS3)), size(AWS3,1), size(AWS3,2), size(AWS3,3)

    INQUIRE(FILE=trim(dirName) // "/AWQV3_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/AWQV3_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) AWQV3
    close(fileID)
    !write(*,*) 'sum(AWQV3) = ', sum((AWQV3)), size(AWQV3,1), size(AWQV3,2), size(AWQV3,3)

    INQUIRE(FILE=trim(dirName) // "/AWQL3_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/AWQL3_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) AWQL3
    close(fileID)
    !write(*,*) 'sum(AWQL3) = ', sum((AWQL3)), size(AWQL3,1), size(AWQL3,2), size(AWQL3,3)

    INQUIRE(FILE=trim(dirName) // "/AWQI3_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/AWQI3_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) AWQI3
    close(fileID)
    !write(*,*) 'sum(AWQI3) = ', sum((AWQI3)), size(AWQI3,1), size(AWQI3,2), size(AWQI3,3)

    INQUIRE(FILE=trim(dirName) // "/AWU3_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/AWU3_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) AWU3
    close(fileID)
    !write(*,*) 'sum(AWU3) = ', sum((AWU3)), size(AWU3,1), size(AWU3,2), size(AWU3,3)

    INQUIRE(FILE=trim(dirName) // "/AWV3_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/AWV3_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) AWV3
    close(fileID)
    !write(*,*) 'sum(AWV3) = ', sum((AWV3)), size(AWV3,1), size(AWV3,2), size(AWV3,3)

    INQUIRE(FILE=trim(dirName) // "/AW3_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/AW3_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) AW3
    close(fileID)
    !write(*,*) 'sum(AW3) = ', sum((AW3)), size(AW3,1), size(AW3,2), size(AW3,3)

    INQUIRE(FILE=trim(dirName) // "/TKH_" // trim(rank_str) // ".in", EXIST=file_exists)
    open(newunit=fileID, FILE=trim(dirName) // "/TKH_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) TKH
    close(fileID)
    ! write(*,*) 'sum(TKH) = ', sum((TKH)), size(TKH,1), size(TKH,2), size(TKH,3)

    ! *** Output Arrays ***
    INQUIRE(FILE=trim(dirName) // "/TKE_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/TKE_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(TKE_ref(IM, JM, 0:LM))
        allocate(TKE(IM, JM, 0:LM))
        read(fileID) TKE_ref
        close(fileID)
        !write(*,*) 'sum(TKE) = ', sum(abs(TKE_ref)), size(TKE_ref,1), size(TKE_ref,2), size(TKE_ref,3), associated(TKE)
     endif

     INQUIRE(FILE=trim(dirName) // "/ISOTROPY_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/ISOTROPY_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(ISOTROPY(IM, JM, LM))
        allocate(ISOTROPY_ref(IM, JM, LM))
        read(fileID) ISOTROPY_ref
        close(fileID)
        !write(*,*) 'sum(ISOTROPY) = ', sum(abs(ISOTROPY_ref)), size(ISOTROPY,1), size(ISOTROPY,2), size(ISOTROPY,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/ZPBL_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/ZPBL_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(ZPBL(IM, JM))
        allocate(ZPBL_ref(IM, JM))
        read(fileID) ZPBL_ref
        close(fileID)
        !write(*,*) 'sum(ZPBL) = ', sum(abs(ZPBL_ref)), size(ZPBL,1), size(ZPBL,2)
     endif

     INQUIRE(FILE=trim(dirName) // "/PPBL_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/PPBL_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(PPBL_ref(IM, JM))
        allocate(PPBL(IM, JM))
        read(fileID) PPBL_ref
        close(fileID)
        !write(*,*) 'sum(PPBL) = ', sum(abs(PPBL_ref)), size(PPBL_ref,1), size(PPBL_ref,2), associated(PPBL)
     endif

     INQUIRE(FILE=trim(dirName) // "/TCZPBL_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/TCZPBL_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(TCZPBL(IM, JM))
        allocate(TCZPBL_ref(IM, JM))
        read(fileID) TCZPBL_ref
        close(fileID)
        !write(*,*) 'sum(TCZPBL) = ', sum(abs(TCZPBL_ref)), size(TCZPBL,1), size(TCZPBL,2)
     endif

     INQUIRE(FILE=trim(dirName) // "/thetavs_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/thetavs_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(thetavs(IM, JM))
        allocate(thetavs_ref(IM, JM))
        read(fileID) thetavs_ref
        close(fileID)
        !write(*,*) 'sum(thetavs) = ', sum(abs(thetavs_ref)), size(thetavs,1), size(thetavs,2)
     endif

     INQUIRE(FILE=trim(dirName) // "/tcrib_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/tcrib_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        ! allocate(tcrib(IM, JM, LM))
        allocate(tcrib_ref(IM, JM, LM))
        read(fileID) tcrib_ref
        close(fileID)
        !write(*,*) 'sum(tcrib) = ', sum(abs(tcrib_ref)), size(tcrib,1), size(tcrib,2), size(tcrib,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/thetavh_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/thetavh_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(thetavh(IM, JM))
        allocate(thetavh_ref(IM, JM))
        read(fileID) thetavh_ref
        close(fileID)
        !write(*,*) 'sum(thetavh) = ', sum(abs(thetavh_ref)), size(thetavh,1), size(thetavh,2)
     endif

     INQUIRE(FILE=trim(dirName) // "/uv2h_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/uv2h_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(uv2h(IM, JM))
        allocate(uv2h_ref(IM, JM))
        read(fileID) uv2h_ref
        close(fileID)
        !write(*,*) 'sum(uv2h) = ', sum(abs(uv2h_ref)), size(uv2h,1), size(uv2h,2)
     endif

     INQUIRE(FILE=trim(dirName) // "/KPBLTC_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/KPBLTC_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(KPBLTC(IM, JM))
        allocate(KPBLTC_ref(IM, JM))
        read(fileID) KPBLTC_ref
        close(fileID)
        !write(*,*) 'sum(KPBLTC) = ', sum(abs(KPBLTC_ref)), size(KPBLTC,1), size(KPBLTC,2)
     endif

     INQUIRE(FILE=trim(dirName) // "/ZPBL2_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/ZPBL2_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(ZPBL2(IM, JM))
        allocate(ZPBL2_ref(IM, JM))
        read(fileID) ZPBL2_ref
        close(fileID)
        !write(*,*) 'sum(ZPBL2) = ', sum(abs(ZPBL2_ref)), size(ZPBL2,1), size(ZPBL2,2)
     endif

     INQUIRE(FILE=trim(dirName) // "/KPBL2_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/KPBL2_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(KPBL2(IM, JM))
        allocate(KPBL2_ref(IM, JM))
        read(fileID) KPBL2_ref
        close(fileID)
        !write(*,*) 'sum(KPBL2) = ', sum(abs(KPBL2_ref)), size(KPBL2,1), size(KPBL2,2)
     endif

     INQUIRE(FILE=trim(dirName) // "/KPBL_SC_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/KPBL_SC_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(KPBL_SC(IM, JM))
        allocate(KPBL_SC_ref(IM, JM))
        read(fileID) KPBL_SC_ref
        close(fileID)
        !write(*,*) 'sum(KPBL_SC) = ', sum(abs(KPBL_SC_ref)), size(KPBL_SC,1), size(KPBL_SC,2)
     endif

     INQUIRE(FILE=trim(dirName) // "/temparray_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/temparray_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(temparray(1:LM+1))
        allocate(temparray_ref(1:LM+1))
        read(fileID) temparray_ref
        close(fileID)
        !write(*,*) 'sum(temparray) = ', sum(abs(temparray_ref)), size(temparray,1)
     endif

     INQUIRE(FILE=trim(dirName) // "/ZPBL10p_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/ZPBL10p_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(ZPBL10p(IM, JM))
        allocate(ZPBL10p_ref(IM, JM))
        read(fileID) ZPBL10p_ref
        close(fileID)
        !write(*,*) 'sum(ZPBL10p) = ', sum(abs(ZPBL10p_ref)), size(ZPBL10p,1), size(ZPBL10p,2)
     endif

     INQUIRE(FILE=trim(dirName) // "/KPBL10p_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/KPBL10p_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(KPBL10p(IM, JM))
        allocate(KPBL10p_ref(IM, JM))
        read(fileID) KPBL10p_ref
        close(fileID)
        !write(*,*) 'sum(KPBL10p) = ', sum(abs(KPBL10p_ref)), size(KPBL10p,1), size(KPBL10p,2)
     endif

      INQUIRE(FILE=trim(dirName) // "/ZPBLHTKE_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/ZPBLHTKE_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(ZPBLHTKE_ref(IM, JM))
        allocate(ZPBLHTKE(IM, JM))
        read(fileID) ZPBLHTKE_ref
        close(fileID)
        !write(*,*) 'sum(ZPBLHTKE) = ', sum(abs(ZPBLHTKE)), size(ZPBLHTKE,1), size(ZPBLHTKE,2), associated(ZPBLHTKE)
     endif

      INQUIRE(FILE=trim(dirName) // "/ZPBLRI_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/ZPBLRI_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(ZPBLRI_ref(IM, JM))
        allocate(ZPBLRI(IM, JM))
        read(fileID) ZPBLRI_ref
        close(fileID)
        !write(*,*) 'sum(ZPBLRI) = ', sum(abs(ZPBLRI_ref)), size(ZPBLRI_ref,1), size(ZPBLRI_ref,2), associated(ZPBLRI)
     endif

     INQUIRE(FILE=trim(dirName) // "/ZPBLRI2_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/ZPBLRI2_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(ZPBLRI2_ref(IM, JM))
        allocate(ZPBLRI2(IM, JM))
        read(fileID) ZPBLRI2_ref
        close(fileID)
        !write(*,*) 'sum(ZPBLRI2) = ', sum(abs(ZPBLRI2_ref)), size(ZPBLRI2,1), size(ZPBLRI2,2), associated(ZPBLRI2)
     endif

     INQUIRE(FILE=trim(dirName) // "/ZPBLTHV_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/ZPBLTHV_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(ZPBLTHV_ref(IM, JM))
        allocate(ZPBLTHV(IM, JM))
        read(fileID) ZPBLTHV_ref
        close(fileID)
        !write(*,*) 'sum(ZPBLTHV) = ', sum(abs(ZPBLTHV_ref)), size(ZPBLTHV_ref,1), size(ZPBLTHV_ref,2), associated(ZPBLTHV)
     endif

     INQUIRE(FILE=trim(dirName) // "/thetav_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/thetav_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(thetav(0:LM))
        allocate(thetav_ref(0:LM))
        read(fileID) thetav_ref
        close(fileID)
        !write(*,*) 'sum(thetav) = ', sum(abs(thetav_ref)), size(thetav,1)
     endif

     INQUIRE(FILE=trim(dirName) // "/KPBL_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/KPBL_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(KPBL(IM, JM))
        allocate(KPBL_ref(IM, JM))
        read(fileID) KPBL_ref
        close(fileID)
        !write(*,*) 'sum(KPBL) = ', sum(abs(KPBL_ref)), size(KPBL,1), size(KPBL,2)
     endif

     INQUIRE(FILE=trim(dirName) // "/CKS_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/CKS_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(CKS(IM, JM, LM))
        allocate(CKS_ref(IM, JM, LM))
        read(fileID) CKS_ref
        close(fileID)
        !write(*,*) 'sum(CKS) = ', sum(abs(CKS_ref)), size(CKS,1), size(CKS,2), size(CKS,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/AKS_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/AKS_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(AKS(IM, JM, LM))
        allocate(AKS_ref(IM, JM, LM))
        read(fileID) AKS_ref
        close(fileID)
        !write(*,*) 'sum(AKS) = ', sum(abs(AKS_ref)), size(AKS,1), size(AKS,2), size(AKS,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/TKH_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/TKH_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(TKH_ref(IM, JM, 0:LM))
        read(fileID) TKH_ref
        close(fileID)
        !write(*,*) 'sum(TKH) = ', sum(abs(TKH_ref)), size(TKH,1), size(TKH,2), size(TKH,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/AKQ_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/AKQ_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(AKQ(IM, JM, LM))
        allocate(AKQ_ref(IM, JM, LM))
        read(fileID) AKQ_ref
        close(fileID)
        !write(*,*) 'sum(AKQ) = ', sum(abs(AKQ_ref)), size(AKQ,1), size(AKQ,2), size(AKQ,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/CKQ_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/CKQ_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(CKQ(IM, JM, LM))
        allocate(CKQ_ref(IM, JM, LM))
        read(fileID) CKQ_ref
        close(fileID)
        !write(*,*) 'sum(CKQ) = ', sum(abs(CKQ_ref)), size(CKQ,1), size(CKQ,2), size(CKQ,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/EKV_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/EKV_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(EKV(IM, JM, LM))
        allocate(EKV_ref(IM, JM, LM))
        read(fileID) EKV_ref
        close(fileID)
        !write(*,*) 'sum(EKV) = ', sum(abs(EKV_ref)), size(EKV,1), size(EKV,2), size(EKV,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/AKV_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/AKV_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(AKV(IM, JM, LM))
        allocate(AKV_ref(IM, JM, LM))
        read(fileID) AKV_ref
        close(fileID)
        !write(*,*) 'sum(AKV) = ', sum(abs(AKV_ref)), size(AKV,1), size(AKV,2), size(AKV,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/CKV_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/CKV_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(CKV(IM, JM, LM))
        allocate(CKV_ref(IM, JM, LM))
        read(fileID) CKV_ref
        close(fileID)
        !write(*,*) 'sum(CKV) = ', sum(abs(CKV_ref)), size(CKV,1), size(CKV,2), size(CKV,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/BKS_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/BKS_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(BKS(IM, JM, LM))
        allocate(BKS_ref(IM, JM, LM))
        read(fileID) BKS_ref
        close(fileID)
        !write(*,*) 'sum(BKS) = ', sum(abs(BKS_ref)), size(BKS,1), size(BKS,2), size(BKS,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/BKQ_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/BKQ_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(BKQ(IM, JM, LM))
        allocate(BKQ_ref(IM, JM, LM))
        read(fileID) BKQ_ref
        close(fileID)
        !write(*,*) 'sum(BKQ) = ', sum(abs(BKQ_ref)), size(BKQ,1), size(BKQ,2), size(BKQ,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/BKV_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/BKV_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(BKV(IM, JM, LM))
        allocate(BKV_ref(IM, JM, LM))
        read(fileID) BKV_ref
        close(fileID)
        !write(*,*) 'sum(BKV) = ', sum(abs(BKV_ref)), size(BKV,1), size(BKV,2), size(BKV,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/AKSS_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/AKSS_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(AKSS(IM, JM, LM))
        allocate(AKSS_ref(IM, JM, LM))
        read(fileID) AKSS_ref
        close(fileID)
        !write(*,*) 'sum(AKSS) = ', sum(abs(AKSS_ref)), size(AKSS,1), size(AKSS,2), size(AKSS,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/AKUU_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/AKUU_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(AKUU(IM, JM, LM))
        allocate(AKUU_ref(IM, JM, LM))
        read(fileID) AKUU_ref
        close(fileID)
        !write(*,*) 'sum(AKUU) = ', sum(abs(AKUU_ref)), size(AKUU,1), size(AKUU,2), size(AKUU,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/RHOAW3_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/RHOAW3_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(RHOAW3(IM, JM, 0:LM))
        allocate(RHOAW3_ref(IM, JM, 0:LM))
        read(fileID) RHOAW3_ref
        close(fileID)
        !write(*,*) 'sum(RHOAW3) = ', sum(abs(RHOAW3_ref)), size(RHOAW3,1), size(RHOAW3,2), size(RHOAW3,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/AKQQ_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/AKQQ_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(AKQQ(IM, JM, LM))
        allocate(AKQQ_ref(IM, JM, LM))
        read(fileID) AKQQ_ref
        close(fileID)
        !write(*,*) 'sum(AKQQ) = ', sum(abs(AKQQ_ref)), size(AKQQ,1), size(AKQQ,2), size(AKQQ,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/CKSS_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/CKSS_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(CKSS(IM, JM, LM))
        allocate(CKSS_ref(IM, JM, LM))
        read(fileID) CKSS_ref
        close(fileID)
        !write(*,*) 'sum(CKSS) = ', sum(abs(CKSS_ref)), size(CKSS,1), size(CKSS,2), size(CKSS,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/CKQQ_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/CKQQ_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(CKQQ(IM, JM, LM))
        allocate(CKQQ_ref(IM, JM, LM))
        read(fileID) CKQQ_ref
        close(fileID)
        !write(*,*) 'sum(CKQQ) = ', sum(abs(CKQQ_ref)), size(CKQQ,1), size(CKQQ,2), size(CKQQ,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/CKUU_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/CKUU_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(CKUU(IM, JM, LM))
        allocate(CKUU_ref(IM, JM, LM))
        read(fileID) CKUU_ref
        close(fileID)
        !write(*,*) 'sum(CKUU) = ', sum(abs(CKUU_ref)), size(CKUU,1), size(CKUU,2), size(CKUU,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/BKSS_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/BKSS_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(BKSS(IM, JM, LM))
        allocate(BKSS_ref(IM, JM, LM))
        read(fileID) BKSS_ref
        close(fileID)
        !write(*,*) 'sum(BKSS) = ', sum(abs(BKSS_ref)), size(BKSS,1), size(BKSS,2), size(BKSS,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/BKQQ_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/BKQQ_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(BKQQ(IM, JM, LM))
        allocate(BKQQ_ref(IM, JM, LM))
        read(fileID) BKQQ_ref
        close(fileID)
        !write(*,*) 'sum(BKQQ) = ', sum(abs(BKQQ_ref)), size(BKQQ,1), size(BKQQ,2), size(BKQQ,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/BKUU_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/BKUU_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(BKUU(IM, JM, LM))
        allocate(BKUU_ref(IM, JM, LM))
        read(fileID) BKUU_ref
        close(fileID)
        !write(*,*) 'sum(BKUU) = ', sum(abs(BKUU_ref)), size(BKUU,1), size(BKUU,2), size(BKUU,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/YS_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/YS_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(YS(IM, JM, LM))
        allocate(YS_ref(IM, JM, LM))
        read(fileID) YS_ref
        close(fileID)
        !write(*,*) 'sum(YS) = ', sum(abs(YS_ref)), size(YS,1), size(YS,2), size(YS,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/YQV_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/YQV_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(YQV(IM, JM, LM))
        allocate(YQV_ref(IM, JM, LM))
        read(fileID) YQV_ref
        close(fileID)
        !write(*,*) 'sum(YQV) = ', sum(abs(YQV_ref)), size(YQV,1), size(YQV,2), size(YQV,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/YQL_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/YQL_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(YQL(IM, JM, LM))
        allocate(YQL_ref(IM, JM, LM))
        read(fileID) YQL_ref
        close(fileID)
        !write(*,*) 'sum(YQL) = ', sum(abs(YQL_ref)), size(YQL,1), size(YQL,2), size(YQL,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/YQI_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/YQI_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(YQI(IM, JM, LM))
        allocate(YQI_ref(IM, JM, LM))
        read(fileID) YQI_ref
        close(fileID)
        !write(*,*) 'sum(YQI) = ', sum(abs(YQI_ref)), size(YQI,1), size(YQI,2), size(YQI,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/YU_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/YU_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(YU(IM, JM, LM))
        allocate(YU_ref(IM, JM, LM))
        read(fileID) YU_ref
        close(fileID)
        !write(*,*) 'sum(YU) = ', sum(abs(YU_ref)), size(YU,1), size(YU,2), size(YU,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/YV_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/YV_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(YV(IM, JM, LM))
        allocate(YV_ref(IM, JM, LM))
        read(fileID) YV_ref
        close(fileID)
        !write(*,*) 'sum(YV) = ', sum(abs(YV_ref)), size(YV,1), size(YV,2), size(YV,3)
     endif

     INQUIRE(FILE=trim(dirName) // "/AKSODT_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/AKSODT_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(AKSODT_ref(IM, JM, LM))
        allocate(AKSODT(IM, JM, LM))
        read(fileID) AKSODT_ref
        close(fileID)
        !write(*,*) 'sum(AKSODT) = ', sum(abs(AKSODT_ref)), size(AKSODT_ref,1), size(AKSODT_ref,2), size(AKSODT_ref,3), associated(AKSODT)
     endif

     INQUIRE(FILE=trim(dirName) // "/CKSODT_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/CKSODT_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(CKSODT_ref(IM, JM, LM))
        allocate(CKSODT(IM, JM, LM))
        read(fileID) CKSODT_ref
        close(fileID)
        !write(*,*) 'sum(CKSODT) = ', sum(abs(CKSODT_ref)), size(CKSODT_ref,1), size(CKSODT_ref,2), size(CKSODT_ref,3), associated(CKSODT)
     endif

     INQUIRE(FILE=trim(dirName) // "/AKQODT_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/AKQODT_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(AKQODT_ref(IM, JM, LM))
        allocate(AKQODT(IM, JM, LM))
        read(fileID) AKQODT_ref
        close(fileID)
        !write(*,*) 'sum(AKQODT) = ', sum(abs(AKQODT_ref)), size(AKQODT_ref,1), size(AKQODT_ref,2), size(AKQODT_ref,3), associated(AKQODT)
     endif

     INQUIRE(FILE=trim(dirName) // "/CKQODT_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/CKQODT_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(CKQODT_ref(IM, JM, LM))
        allocate(CKQODT(IM, JM, LM))
        read(fileID) CKQODT_ref
        close(fileID)
        !write(*,*) 'sum(CKQODT) = ', sum(abs(CKQODT_ref)), size(CKQODT_ref,1), size(CKQODT_ref,2), size(CKQODT_ref,3), associated(CKQODT)
     endif

     INQUIRE(FILE=trim(dirName) // "/AKVODT_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/AKVODT_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(AKVODT_ref(IM, JM, LM))
        allocate(AKVODT(IM, JM, LM))
        read(fileID) AKVODT_ref
        close(fileID)
        !write(*,*) 'sum(AKVODT) = ', sum(abs(AKVODT_ref)), size(AKVODT_ref,1), size(AKVODT_ref,2), size(AKVODT_ref,3), associated(AKVODT)
     endif

     INQUIRE(FILE=trim(dirName) // "/CKVODT_" // trim(rank_str) // ".out", EXIST=file_exists)
      if(file_exists.eqv..true.) then
        open(newunit=fileID, FILE=trim(dirName) // "/CKVODT_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
        allocate(CKVODT_ref(IM, JM, LM))
        allocate(CKVODT(IM, JM, LM))
        read(fileID) CKVODT_ref
        close(fileID)
        !write(*,*) 'sum(CKVODT) = ', sum(abs(CKVODT_ref)), size(CKVODT_ref,1), size(CKVODT_ref,2), size(CKVODT_ref,3), associated(CKVODT)
     endif

!$acc data copyin(TKESHOC, KH, THV, Z, KM, U, V, QLLS, PREF, T, Q, TH, KHSFC, ZLE, RI, FRLAND, PLO, RDZ, DMI, CT, PLE, TV, CQ, CU, AE3, &
!$acc             RHOE, AWS3, AWQV3, AWQL3, AWQI3, AWU3, AW3) &
!$acc     copyout(TKE, ISOTROPY, ZPBL, PPBL, TCZPBL, thetavs, tcrib, thetavh, uv2h, KPBLTC, ZPBL2, KPBL2, KPBL_SC, ZPBL10p, &
!$acc             KPBL10p, ZPBLHTKE, ZPBLRI, ZPBLRI2, ZPBLTHV, thetav, KPBL, CKS, AKS, AKQ, CKQ, EKV, AKV, CKV, BKS, BKQ, BKV, AKSS, AKUU, &
!$acc             RHOAW3, AKQQ, CKSS, CKQQ, CKUU, BKSS, BKQQ, BKUU, YS, YQV, YQL, YQI, YU, YV, AKSODT, CKSODT, AKQODT, CKQODT, AKVODT, CKVODT) &
!$acc     create(AWV3) copy(TKH)

    call cpu_time(t1)
    
    ! TKE 
    if (associated(TKE)) then         ! Reminder: TKE is on model edges
        if (DO_SHOC /= 0) then          !           TKESHOC is not.
!$acc kernels
          TKE(:,:,1:LM-1) = 0.5*(TKESHOC(:,:,1:LM-1)+TKESHOC(:,:,2:LM))
          TKE(:,:,0) = 1e-6
          TKE(:,:,LM) = 1e-6
!$acc end kernels
        else
!$acc kernels
          TKE = MAPL_UNDEF
          do L = 1,LM-1
            TKE(:,:,L) = ( LAMBDADISS * &
            ( -1.*(KH(:,:,L)*MAPL_GRAV/((THV(:,:,L) + THV(:,:,L+1))*0.5) *  ((THV(:,:,L) - THV(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1)))) +  &
            (KM(:,:,L)*((U(:,:,L) - U(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1)))*((U(:,:,L) - U(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1))))  +  &
            (KM(:,:,L)*((V(:,:,L) - V(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1)))*((V(:,:,L) - V(:,:,L+1))/(Z(:,:,L) - Z(:,:,L+1)))) )) ** 2
            TKE(:,:,L) = TKE(:,:,L) ** (1./3.)
          enddo

          ! If not running SHOC, estimate ISOTROPY from KH and TKE,
          ! based on Eq. 7 from Bogenschutz and Krueger (2013).
          ! This is a placeholder to allow use of the double-gaussian
          ! PDF without SHOC, but should be tested and revised!
          ISOTROPY(:,:,LM) = KH(:,:,LM-1) / max(0.01,0.1*TKE(:,:,LM-1))
          ISOTROPY(:,:,1) = KH(:,:,1) / max(0.01,0.1*TKE(:,:,1))
          do L = 2,LM-1
            ISOTROPY(:,:,L) = (KH(:,:,L)+KH(:,:,L-1)) / (0.1*(TKE(:,:,L)+TKE(:,:,L-1)))
          end do
          ISOTROPY = max(10.,min(2000.,ISOTROPY))
!$acc end kernels
        end if
    end if ! TKE

    ! Update the higher order moments required for the ADG PDF
    ! NOTE : Commented this portion out since the output values from update_moments
    !        are not used within postlock
    ! if (PDFSHAPE.eq.5) then
    !     HL = T + (mapl_grav*Z - mapl_alhl*QLLS)/mapl_cp
    !     call update_moments(IM, JM, LM, DT, &
    !                         SH,             &  ! in
    !                         EVAP,           &
    !                         Z,              &
    !                         KH,             &
    !                         TKESHOC,        &
    !                         ISOTROPY,       &
    !                         QT,             &
    !                         HL,             &
    !                         EDMF_FRC,       &
    !                         MFQT2,          &
    !                         MFQT3,          &
    !                         MFHL2,          &
    !                         MFHL3,          &
    !                         MFW2,           &
    !                         MFW3,           &
    !                         MFWQT,          &
    !                         MFWHL,          &
    !                         MFHLQT,         &
    !                         qt2,            &  ! inout
    !                         qt3,            &
    !                         hl2,            &  ! out
    !                         hl3,            &
    !                         w2,             &
    !                         w3,             &
    !                         wqt,            &
    !                         whl,            &
    !                         hlqt,           &
    !                         hl2tune,        &  ! tuning parameters
    !                         qt2tune,        &
    !                         hlqt2tune,      &
    !                         qt2scale,       &
    !                         qt3_tscale )
    ! end if
!$acc kernels
    KPBLMIN  = count(PREF < 50000.)

    ZPBL = MAPL_UNDEF
    if (associated(PPBL)) PPBL = MAPL_UNDEF
!$acc end kernels

    if (CALC_TCZPBL) then
!$acc kernels
        TCZPBL = MAPL_UNDEF
!$acc end kernels
        if (LM .eq. 72) then
!$acc kernels
            thetavs = T(:,:,LM)*(1.0+MAPL_VIREPS*Q(:,:,LM)/(1.0-Q(:,:,LM)))*(TH(:,:,LM)/T(:,:,LM))
            tcrib(:,:,LM) = 0.0
            do J = 1, JM
                do I = 1, IM
                    do L=LM-1,1,-1
                        thetavh(I,J) = T(I,J,L)*(1.0+MAPL_VIREPS*Q(I,J,L)/(1.0-Q(I,J,L)))*(TH(I,J,L)/T(I,J,L))
                        uv2h(I,J) = max(U(I,J,L)**2+V(I,J,L)**2,1.0E-8)
                        tcrib(I,J,L) = MAPL_GRAV*(thetavh(I,J)-thetavs(I,J))*Z(I,J,L)/(thetavs(I,J)*uv2h(I,J))
                        if (tcrib(I,J,L) >= tcri_crit) then
                            TCZPBL(I,J) = Z(I,J,L+1)+(tcri_crit-tcrib(I,J,L+1))/(tcrib(I,J,L)-tcrib(I,J,L+1))*(Z(I,J,L)-Z(I,J,L+1))
                            KPBLTC(I,J) = float(L)
                            exit
                        end if
                    end do
                end do
            end do
!$acc end kernels
        else
!$acc kernels
            tcrib(:,:,LM) = 0.0
            do I = 1, IM
                do J = 1, JM
                    do L=LM-1,1,-1
                        uv2h(I,J) = max(U(I,J,L)**2+V(I,J,L)**2,1.0E-8)
                        tcrib(I,J,L) = MAPL_GRAV*(THV(I,J,L)-THV(I,J,LM))*Z(I,J,L)/(THV(I,J,LM)*uv2h(I,J))
                        if (tcrib(I,J,L) >= tcri_crit) then
                            TCZPBL(I,J) = Z(I,J,L+1)+(tcri_crit-tcrib(I,J,L+1))/(tcrib(I,J,L)-tcrib(I,J,L+1))*(Z(I,J,L)-Z(I,J,L+1))
                            KPBLTC(I,J) = float(L)
                            exit
                        end if
                    end do
                end do
            end do
!$acc end kernels
        endif
!$acc kernels
        do J = 1,JM
            do I = 1,IM
                if(TCZPBL(I,J) < 0.) then
                    TCZPBL(I,J) = Z(I,J,LM)
                    KPBLTC(I,J) = float(LM)
                endif
            enddo
        enddo
!$acc end kernels
    end if ! CALC_TCZPBL

    if (CALC_ZPBL2) then
!$acc kernels
        ZPBL2 = MAPL_UNDEF

        do I = 1, IM
            do J = 1, JM
                do L=LM,2,-1
                    if ((KH(I,J,L-1) < 2.).and.(KH(I,J,L) >= 2.).and.(ZPBL2(I,J)==MAPL_UNDEF)) then
                        ZPBL2(I,J) = Z(I,J,L)
                        KPBL2(I,J) = float(L)
                    end if
                end do
            end do
        end do

        where ( ZPBL2 .eq. MAPL_UNDEF )
            ZPBL2 = Z(:,:,LM)
            KPBL2 = float(LM)
        end where
        ZPBL2 = MIN(ZPBL2,Z(:,:,KPBLMIN))
!$acc end kernels
    end if ! CALC_ZPBL2

!$acc kernels
    ! Calc KPBL using surface turbulence, for use in shallow scheme
    KPBL_SC = MAPL_UNDEF
!$acc end kernels

!$acc parallel

!$acc loop gang vector collapse(2) private(temparray)
    do I = 1, IM
        do J = 1, JM
            if (DO_SHOC==0) then
                temparray(1:LM+1) = KHSFC(I,J,0:LM)
            else
                temparray(1:LM+1) = KH(I,J,0:LM)
            end if

            maxkh = maxval(temparray)
!$acc loop seq
            do L=LM-1,2,-1
                if ( (temparray(L) < 0.1*maxkh) .and. (temparray(L+1) >= 0.1*maxkh)  &
                .and. (KPBL_SC(I,J) == MAPL_UNDEF ) ) then
                    KPBL_SC(I,J) = float(L)
                    exit
                end if
            end do
            if (  KPBL_SC(I,J) .eq. MAPL_UNDEF .or. (maxkh.lt.1.)) then
                KPBL_SC(I,J) = float(LM)
            endif
        end do
    end do

!$acc end parallel

    if (CALC_ZPBL10p) then

!$acc kernels
        ZPBL10p = MAPL_UNDEF
!$acc end kernels

!$acc parallel
!$acc loop gang vector collapse(2) private(temparray)
        do I = 1, IM
            do J = 1, JM
                temparray(1:LM+1) = KH(I,J,0:LM)
!$acc loop seq
                do L = LM,2,-1
                    locmax = maxloc(temparray,1)
                    ! Suggested fix for potential out of bounds array access
                    locmax = min(locmax,LM)
                    minlval = max(0.001,0.0001*temparray(locmax))
                    if(temparray(locmax-1)<minlval.and.temparray(locmax+1)<minlval) temparray(locmax) = minlval
                enddo
                maxkh = temparray(LM)
!$acc loop seq
                do L = LM-1,2,-1
                    if(temparray(L)>maxkh) maxkh = temparray(L)
                    if(temparray(L-1)<minlval) exit
                end do
!$acc loop seq
                do L=LM-1,2,-1
                    if ( (temparray(L) < 0.1*maxkh) .and. (temparray(L+1) >= 0.1*maxkh)  &
                    .and. (ZPBL10p(I,J) == MAPL_UNDEF ) ) then
                        ZPBL10p(I,J) = ZLE(I,J,L)+ &
                    ((ZLE(I,J,L-1)-ZLE(I,J,L))/(temparray(L)-temparray(L+1))) * (0.1*maxkh-temparray(L+1))
                        KPBL10p(I,J) = float(L)
                        exit
                    end if
                end do
                if (  ZPBL10p(I,J) .eq. MAPL_UNDEF .or. (maxkh.lt.1.)) then
                    ZPBL10p(I,J) = Z(I,J,LM)
                    KPBL10p(I,J) = float(LM)
                endif
            end do
        end do
!$acc end parallel

!$acc kernels
        ZPBL10p = MIN(ZPBL10p,Z(:,:,KPBLMIN))
!$acc end kernels
    end if ! CALC_ZPBL10p

    ! HTKE pbl height
    if (associated(ZPBLHTKE)) then
!$acc kernels
        ZPBLHTKE = MAPL_UNDEF
!$acc end kernels
    end if ! ZPBLHTKE

    ! RI local diagnostic for pbl height thresh 0.
    if (associated(ZPBLRI)) then
!$acc kernels
        ZPBLRI = MAPL_UNDEF
        where (RI(:,:,LM-1)>ri_crit) ZPBLRI = Z(:,:,LM)

        do I = 1, IM
        do J = 1, JM
            do L=LM-1,1,-1
                if( (RI(I,J,L-1)>ri_crit) .and. (ZPBLRI(I,J) == MAPL_UNDEF) ) then
                    ZPBLRI(I,J) = Z(I,J,L+1)+(ri_crit-RI(I,J,L))/(RI(I,J,L-1)-RI(I,J,L))*(Z(I,J,L)-Z(I,J,L+1))
                end if
            end do
        end do 
        end do 

        where ( ZPBLRI .eq. MAPL_UNDEF ) ZPBLRI = Z(:,:,LM)
        ZPBLRI = MIN(ZPBLRI,Z(:,:,KPBLMIN))
        where ( ZPBLRI < 0.0 ) ZPBLRI = Z(:,:,LM)
!$acc end kernels
    end if ! ZPBLRI

    ! RI local diagnostic for pbl height thresh 0.2
    if (associated(ZPBLRI2)) then
!$acc kernels
        ZPBLRI2 = MAPL_UNDEF
        where (RI(:,:,LM-1) > ri_crit2) ZPBLRI2 = Z(:,:,LM)

        do I = 1, IM
            do J = 1, JM
                do L=LM-1,1,-1
                    if( (RI(I,J,L-1)>ri_crit2) .and. (ZPBLRI2(I,J) == MAPL_UNDEF) ) then
                        ZPBLRI2(I,J) = Z(I,J,L+1)+(ri_crit2-RI(I,J,L))/(RI(I,J,L-1)-RI(I,J,L))*(Z(I,J,L)-Z(I,J,L+1))
                    end if
                end do
            end do
        end do

        where ( ZPBLRI2 .eq. MAPL_UNDEF ) ZPBLRI2 = Z(:,:,LM)
        ZPBLRI2 = MIN(ZPBLRI2,Z(:,:,KPBLMIN))
        where ( ZPBLRI2 < 0.0 ) ZPBLRI2 = Z(:,:,LM)
!$acc end kernels
    end if ! ZPBLRI2

    ! thetav gradient based pbl height diagnostic
    if (associated(ZPBLTHV)) then
!$acc kernels
        ZPBLTHV = MAPL_UNDEF
        do I = 1, IM
            do J = 1, JM
                do L=LM,1,-1
                    thetav(L) = T(I,J,L)*(1.0*MAPL_VIREPS*Q(I,J,L)/(1.0-Q(I,J,L)))*(TH(I,J,L)/T(I,J,L))
                end do

                maxdthvdz = 0

                do L=LM-1,1,-1
                    if(Z(I,J,L)<=Z(I,J,KPBLMIN)) then
                        dthvdz = (thetav(L+1)-thetav(L))/(Z(I,J,L+1)-Z(I,J,L))
                        if(dthvdz>maxdthvdz) then
                            maxdthvdz = dthvdz
                            ZPBLTHV(I,J) = 0.5*(Z(I,J,L+1)+Z(I,J,L))
                        end if
                    end if
                end do

            end do 
        end do 
!$acc end kernels
    end if ! ZPBLTHV

    SELECT CASE(PBLHT_OPTION)

        CASE( 1 )
!$acc kernels
            ZPBL = ZPBL2
            KPBL = KPBL2
!$acc end kernels
        CASE( 2 )
!$acc kernels
            ZPBL = ZPBL10p
            KPBL = KPBL10P
!$acc end kernels
        CASE( 3 )
!$acc kernels
            ZPBL = TCZPBL
            KPBL = KPBLTC
!$acc end kernels
        CASE( 4 )
!$acc kernels
            do J = 1,JM
                do I = 1,IM
                    if(FRLAND(I,J) > 0) then
                        ZPBL(I,J) = TCZPBL(I,J)
                        KPBL(I,J) = KPBLTC(I,J)

                    else
                        ZPBL(I,J) = ZPBL10p(I,J)
                        KPBL(I,J) = KPBL10p(I,J)
                    endif
                enddo
            enddo
!$acc end kernels
            ! WHERE (FRLAND(:,:)>0)
            ! ZPBL = TCZPBL
            ! KPBL = KPBLTC

            ! ELSEWHERE
            ! ZPBL = ZPBL10p
            ! KPBL = KPBL10P

            ! END WHERE

    END SELECT

!$acc kernels
    ZPBL = MIN(ZPBL,Z(:,:,KPBLMIN))
    KPBL = MAX(KPBL,float(KPBLMIN))
!$acc end kernels
    if (associated(PPBL)) then
!$acc kernels
        do J = 1, JM
            do I = 1, IM
                PPBL(I,J) = PLO(I,J,nint(KPBL(I,J)))
                PPBL(I,J) = MAX(PPBL(I,J), PLO(I,J,KPBLMIN))
            end do
        end do
        ! PPBL = MAX(PPBL,PLO(:,:,KPBLMIN))
!$acc end kernels
    end if

    ! Second difference coefficients for scalars; RDZ is RHO/DZ, DMI is (G DT)/DP
    ! ---------------------------------------------------------------------------

!$acc kernels
    do J = 1, JM
        do I = 1, IM
            AKS(I,J,1) = 0.0
            CKS(I,J,LM) = -CT(I,J) * DMI(I,J,LM)
        enddo
    enddo

    DO L = 1, LM-1
        do J = 1, JM
            do I = 1,IM
                CKS(I,J,L) = -KH(I,J,L) * RDZ(I,J,L)
                AKS(I,J,L+1) = CKS(I,J,L) * DMI(I,J,L+1)
                CKS(I,J,L) = CKS(I,J,L) * DMI(I,J,L)
            enddo
        enddo
    enddo

    ! CKS(:,:,1:LM-1) = -KH(:,:,1:LM-1) * RDZ(:,:,1:LM-1)
    ! AKS(:,:,1     ) = 0.0
    ! AKS(:,:,2:LM  ) = CKS(:,:,1:LM-1) * DMI(:,:,2:LM  )
    ! CKS(:,:,1:LM-1) = CKS(:,:,1:LM-1) * DMI(:,:,1:LM-1)
    ! CKS(:,:,  LM  ) = -CT             * DMI(:,:,  LM  )

    ! Fill KH at level LM+1 with CT * RDZ for diagnostic output
    ! ---------------------------------------------------------
    do J = 1,JM
        do I = 1,IM
            KH(I,J,LM) = CT(I,J) * (PLE(I,J,LM)/(MAPL_RGAS * TV(I,J,LM))) / Z(I,J,LM)
            TKH(I,J,LM) = KH(I,J,LM)
        enddo
    enddo

    ! KH(:,:,LM) = CT * (PLE(:,:,LM)/(MAPL_RGAS * TV(:,:,LM))) / Z(:,:,LM)
    ! TKH(:,:,LM) = KH(:,:,LM)

    ! Water vapor can differ at the surface
    !--------------------------------------

    ! AKQ(:,:,:)  = AKS(:,:,:)
    ! CKQ(:,:,:)  = CKS(:,:,:)

    do L = 1,LM
        do J = 1,JM
            do I = 1,IM
                AKQ(I,J,L) = AKS(I,J,L)
                CKQ(I,J,L) = CKS(I,J,L)
            enddo
        enddo
    enddo

    do J = 1,JM
        do I = 1,IM
            CKQ(I,J,LM) = -CQ(I,J) * DMI(I,J,LM)
            AKV(I,J,1)  = 0.0
            CKV(I,J,LM) = - CU(I,J) * DMI(I,J,LM)
            EKV(I,J,LM) = MAPL_GRAV * CU(I,J)

            ! Fill KM at level LM with CU * RDZ for diagnostic output
            ! -------------------------------------------------------

            KM(I,J,LM) = CU(I,J) * (PLE(I,J,LM)/(MAPL_RGAS * TV(I,J,LM))) / Z(I,J,LM)
        enddo
    enddo
    ! CKQ(:,:,LM) = -CQ * DMI(:,:,LM)

    ! Second difference coefficients for winds
    ! EKV is saved to use in the frictional heating calc.
    ! ---------------------------------------------------
    
    do L = 1,LM-1
        do J = 1,JM
            do I = 1,IM
                EKV(I,J,L)   = -KM(I,J,L)  * RDZ(I,J,L)
                AKV(I,J,L+1) =  EKV(I,J,L) * DMI(I,J,L+1)
                CKV(I,J,L)   =  EKV(I,J,L) * DMI(I,J,L)
                EKV(I,J,L)   = -MAPL_GRAV  * EKV(I,J,L)
            enddo
        enddo
    enddo

    ! EKV(:,:,1:LM-1) = -KM(:,:,1:LM-1) * RDZ(:,:,1:LM-1)
    ! AKV(:,:,1     ) = 0.0
    ! AKV(:,:,2:LM  ) = EKV(:,:,1:LM-1) * DMI(:,:,2:LM  )
    ! CKV(:,:,1:LM-1) = EKV(:,:,1:LM-1) * DMI(:,:,1:LM-1)
    ! EKV(:,:,1:LM-1) = -MAPL_GRAV      * EKV(:,:,1:LM-1)

    ! CKV(:,:,  LM  ) = -  CU           * DMI(:,:,  LM  )
    ! EKV(:,:,  LM  ) =  MAPL_GRAV      * CU

    ! Fill KM at level LM with CU * RDZ for diagnostic output
    ! -------------------------------------------------------

    ! KM(:,:,LM) = CU * (PLE(:,:,LM)/(MAPL_RGAS * TV(:,:,LM))) / Z(:,:,LM)

    ! Setup the tridiagonal matrix
    ! ----------------------------

    do L = 1,LM
        do J = 1,JM
            do I = 1,IM
                BKS(I,J,L) = 1.00 - (AKS(I,J,L)+CKS(I,J,L))
                BKQ(I,J,L) = 1.00 - (AKQ(I,J,L)+CKQ(I,J,L))
                BKV(I,J,L) = 1.00 - (AKV(I,J,L)+CKV(I,J,L))
            enddo
        enddo
    enddo

    ! BKS = 1.00 - (AKS+CKS)
    ! BKQ = 1.00 - (AKQ+CKQ)
    ! BKV = 1.00 - (AKV+CKV)

    !
    ! A,B,C,D-s for mass flux
    !
        
    do J = 1,JM
        do I = 1,IM
            AKSS(I,J,1)=0.0
            AKUU(I,J,1)=0.0        
        enddo
    enddo
    ! AKSS(:,:,1)=0.0
    ! !  AKQQ(:,:,1)=0.0
    ! AKUU(:,:,1)=0.0


    RHOAW3=RHOE*AW3

    if (EDMFPARAMS_%IMPLICIT == 1 .and. EDMFPARAMS_%DISCRETE == 0) then
        ! AKSS(:,:,2:LM) = - KH(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,2:LM) &
        !                 - 0.5*DMI(:,:,2:LM)*RHOAW3(:,:,1:LM-1)
        ! AKUU(:,:,2:LM) = - KM(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,2:LM) &
        !                 - 0.5*DMI(:,:,2:LM)*RHOAW3(:,:,1:LM-1)
        do L = 1,LM-1
            do J = 1,JM
                do I = 1,IM
                    AKSS(I,J,L+1) = - KH(I,J,L)*RDZ(I,J,L)*AE3(I,J,L)*DMI(I,J,L+1) &
                                    - 0.5*DMI(I,J,L+1)*RHOAW3(I,J,L)
                    AKUU(I,J,L+1) = - KM(I,J,L)*RDZ(I,J,L)*AE3(I,J,L)*DMI(I,J,L+1) &
                                    - 0.5*DMI(I,J,L+1)*RHOAW3(I,J,L)
                enddo
            enddo
        enddo
    else
        AKSS(:,:,2:LM) = - KH(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,2:LM)
        AKUU(:,:,2:LM) = - KM(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,2:LM)
    end if
    AKQQ(:,:,:) = AKSS(:,:,:)

    do J = 1,JM
        do I = 1, IM
            CKSS(I,J,LM)=-CT(I,J)*DMI(I,J,LM)
            CKQQ(I,J,LM)=-CQ(I,J)*DMI(I,J,LM)
            CKUU(I,J,LM)=-CU(I,J)*DMI(I,J,LM)
        enddo
    enddo
    ! CKSS(:,:,LM)=-CT*DMI(:,:,LM)
    ! CKQQ(:,:,LM)=-CQ*DMI(:,:,LM)
    ! CKUU(:,:,LM)=-CU*DMI(:,:,LM)
  
    if (EDMFPARAMS_%IMPLICIT == 1 .and. EDMFPARAMS_%DISCRETE == 0) then
        ! CKSS(:,:,1:LM-1) = - KH(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,1:LM-1) &
        !                     + 0.5*DMI(:,:,1:LM-1)*RHOAW3(:,:,1:LM-1)
        ! CKUU(:,:,1:LM-1) = - KM(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,1:LM-1) &
        !                     + 0.5*DMI(:,:,1:LM-1)*RHOAW3(:,:,1:LM-1)

        do L = 1,LM-1
            do J = 1,JM
                do I = 1,IM
                    CKSS(I,J,L) = - KH(I,J,L)*RDZ(I,J,L)*AE3(I,J,L)*DMI(I,J,L) &
                    + 0.5*DMI(I,J,L)*RHOAW3(I,J,L)
                    CKUU(I,J,L) = - KM(I,J,L)*RDZ(I,J,L)*AE3(I,J,L)*DMI(I,J,L) &
                    + 0.5*DMI(I,J,L)*RHOAW3(I,J,L)

                enddo
            enddo
        enddo
    else
        CKSS(:,:,1:LM-1) = - KH(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,1:LM-1)
        CKUU(:,:,1:LM-1) = - KM(:,:,1:LM-1)*RDZ(:,:,1:LM-1)*AE3(:,:,1:LM-1)*DMI(:,:,1:LM-1)
    end if
    CKQQ(:,:,1:LM-1) = CKSS(:,:,1:LM-1)  
 
    do L = 1,LM
        do J = 1,JM
            do I = 1,IM
                BKSS(I,J,L) = 1.0 - (CKSS(I,J,L)+AKSS(I,J,L))
                BKQQ(I,J,L) = 1.0 - (CKQQ(I,J,L)+AKQQ(I,J,L))
                BKUU(I,J,L) = 1.0 - (CKUU(I,J,L)+AKUU(I,J,L))
            enddo
        enddo
    enddo
    ! BKSS = 1.0 - (CKSS+AKSS)
    ! BKQQ = 1.0 - (CKQQ+AKQQ)
    ! BKUU = 1.0 - (CKUU+AKUU)

    ! Add mass flux contribution
    
    if (EDMFPARAMS_%IMPLICIT == 1) then
        if (EDMFPARAMS_%DISCRETE == 0) then
            ! BKSS(:,:,LM) = BKSS(:,:,LM) - DMI(:,:,LM)*RHOAW3(:,:,LM-1)
            ! BKQQ(:,:,LM) = BKQQ(:,:,LM) - DMI(:,:,LM)*RHOAW3(:,:,LM-1)
            ! BKUU(:,:,LM) = BKUU(:,:,LM) - DMI(:,:,LM)*RHOAW3(:,:,LM-1)

            do J = 1,JM
                do I = 1,IM
                    BKSS(I,J,LM) = BKSS(I,J,LM) - DMI(I,J,LM)*RHOAW3(I,J,LM-1)
                    BKQQ(I,J,LM) = BKQQ(I,J,LM) - DMI(I,J,LM)*RHOAW3(I,J,LM-1)
                    BKUU(I,J,LM) = BKUU(I,J,LM) - DMI(I,J,LM)*RHOAW3(I,J,LM-1)
                enddo
            enddo   

            ! BKSS(:,:,1:LM-1) = BKSS(:,:,1:LM-1) + DMI(:,:,1:LM-1)*( RHOAW3(:,:,1:LM-1) - RHOAW3(:,:,0:LM-2) )
            ! BKQQ(:,:,1:LM-1) = BKQQ(:,:,1:LM-1) + DMI(:,:,1:LM-1)*( RHOAW3(:,:,1:LM-1) - RHOAW3(:,:,0:LM-2) )
            ! BKUU(:,:,1:LM-1) = BKUU(:,:,1:LM-1) + DMI(:,:,1:LM-1)*( RHOAW3(:,:,1:LM-1) - RHOAW3(:,:,0:LM-2) ) 

            do L = 1,LM-1
                do J = 1,JM
                    do I = 1,IM
                        BKSS(I,J,L) = BKSS(I,J,L) + DMI(I,J,L)*( RHOAW3(I,J,L) - RHOAW3(I,J,L-1) )
                        BKQQ(I,J,L) = BKQQ(I,J,L) + DMI(I,J,L)*( RHOAW3(I,J,L) - RHOAW3(I,J,L-1) )
                        BKUU(I,J,L) = BKUU(I,J,L) + DMI(I,J,L)*( RHOAW3(I,J,L) - RHOAW3(I,J,L-1) ) 
                    enddo
                enddo
            enddo
        else if (EDMFPARAMS_%DISCRETE == 1) then
            AKSS(:,:,2:LM) = AKSS(:,:,2:LM) - DMI(:,:,2:LM)*RHOAW3(:,:,1:LM-1)
            AKQQ(:,:,2:LM) = AKQQ(:,:,2:LM) - DMI(:,:,2:LM)*RHOAW3(:,:,1:LM-1)
            AKUU(:,:,2:LM) = AKUU(:,:,2:LM) - DMI(:,:,2:LM)*RHOAW3(:,:,1:LM-1)

            BKSS(:,:,2:LM-1) = BKSS(:,:,2:LM-1) + DMI(:,:,2:LM-1)*RHOAW3(:,:,2:LM-1)
            BKQQ(:,:,2:LM-1) = BKQQ(:,:,2:LM-1) + DMI(:,:,2:LM-1)*RHOAW3(:,:,2:LM-1)
            BKUU(:,:,2:LM-1) = BKUU(:,:,2:LM-1) + DMI(:,:,2:LM-1)*RHOAW3(:,:,2:LM-1)
        end if
    end if

    ! Y-s ... these are rhs - mean value - surface flux 
    ! (these are added in the diffuse and vrtisolve)


    !
    ! 2:LM -> 1:LM-1, 1:LM-1 -> 0:LM-2
    !

    do J = 1, JM
        do I = 1, IM
            YS(I,J,LM)  = -DMI(I,J,LM)*RHOE(I,J,LM-1)*AWS3(I,J,LM-1)
            YQV(I,J,LM) = -DMI(I,J,LM)*RHOE(I,J,LM-1)*AWQV3(I,J,LM-1)
            YQL(I,J,LM) = -DMI(I,J,LM)*RHOE(I,J,LM-1)*AWQL3(I,J,LM-1)

            YQI(I,J,LM) = -DMI(I,J,LM)*AWQI3(I,J,LM-1)*RHOE(I,J,LM-1)
            YU(I,J,LM)  = -DMI(I,J,LM)*AWU3(I,J,LM-1)*RHOE(I,J,LM-1)
            YV(I,J,LM)  = -DMI(I,J,LM)*AWV3(I,J,LM-1)*RHOE(I,J,LM-1)
        enddo
    enddo

    do L = 1, LM-1
        do J = 1, JM
            do I = 1,IM
                YS(I,J,L)  = DMI(I,J,L)*(  RHOE(I,J,L)*AWS3(I,J,L)  - RHOE(I,J,L-1)*AWS3(I,J,L-1) )
                YQV(I,J,L) = DMI(I,J,L)*(  RHOE(I,J,L)*AWQV3(I,J,L) - RHOE(I,J,L-1)*AWQV3(I,J,L-1) )
                YQL(I,J,L) = DMI(I,J,L)*(  RHOE(I,J,L)*AWQL3(I,J,L) - RHOE(I,J,L-1)*AWQL3(I,J,L-1) )

                YQI(I,J,L) = DMI(I,J,L)*( AWQI3(I,J,L)*RHOE(I,J,L) - RHOE(I,J,L-1)*AWQI3(I,J,L-1) )
                YU(I,J,L)  = DMI(I,J,L)*( AWU3(I,J,L)*RHOE(I,J,L)  - RHOE(I,J,L-1)*AWU3(I,J,L-1) )
                YV(I,J,L)  = DMI(I,J,L)*( AWV3(I,J,L)*RHOE(I,J,L)  - RHOE(I,J,L-1)*AWV3(I,J,L-1) )
            enddo
        enddo
    enddo

    ! YS(:,:,LM)  = -DMI(:,:,LM)*RHOE(:,:,LM-1)*AWS3(:,:,LM-1)
    ! YQV(:,:,LM) = -DMI(:,:,LM)*RHOE(:,:,LM-1)*AWQV3(:,:,LM-1)
    ! YQL(:,:,LM) = -DMI(:,:,LM)*RHOE(:,:,LM-1)*AWQL3(:,:,LM-1)

    ! YS(:,:,1:LM-1)  = DMI(:,:,1:LM-1)*(  RHOE(:,:,1:LM-1)*AWS3(:,:,1:LM-1)  - RHOE(:,:,0:LM-2)*AWS3(:,:,0:LM-2) )
    ! YQV(:,:,1:LM-1) = DMI(:,:,1:LM-1)*(  RHOE(:,:,1:LM-1)*AWQV3(:,:,1:LM-1) - RHOE(:,:,0:LM-2)*AWQV3(:,:,0:LM-2) )
    ! YQL(:,:,1:LM-1) = DMI(:,:,1:LM-1)*(  RHOE(:,:,1:LM-1)*AWQL3(:,:,1:LM-1) - RHOE(:,:,0:LM-2)*AWQL3(:,:,0:LM-2) )

    ! YQI(:,:,LM) = -DMI(:,:,LM)*AWQI3(:,:,LM-1)*RHOE(:,:,LM-1)
    ! YU(:,:,LM)  = -DMI(:,:,LM)*AWU3(:,:,LM-1)*RHOE(:,:,LM-1)
    ! YV(:,:,LM)  = -DMI(:,:,LM)*AWV3(:,:,LM-1)*RHOE(:,:,LM-1)

    ! YQI(:,:,1:LM-1) = DMI(:,:,1:LM-1)*( AWQI3(:,:,1:LM-1)*RHOE(:,:,1:LM-1) - RHOE(:,:,0:LM-2)*AWQI3(:,:,0:LM-2) )
    ! YU(:,:,1:LM-1)  = DMI(:,:,1:LM-1)*( AWU3(:,:,1:LM-1)*RHOE(:,:,1:LM-1)  - RHOE(:,:,0:LM-2)*AWU3(:,:,0:LM-2) )
    ! YV(:,:,1:LM-1)  = DMI(:,:,1:LM-1)*( AWV3(:,:,1:LM-1)*RHOE(:,:,1:LM-1)  - RHOE(:,:,0:LM-2)*AWV3(:,:,0:LM-2) )

    ! Add prescribed surface fluxes
! Note : Arrays within USE_SCM_SURF, if they are not used or called above, are not allocated for the OpenACC implementation
#ifdef USE_SCM_SURF
    if ( SCM_SL /= 0 .and. SCM_SL_FLUX == 1 ) then
        YS(:,:,LM)  = YS(:,:,LM)  + DMI(:,:,LM)*SH(:,:)/RHOE(:,:,LM)
        YQV(:,:,LM) = YQV(:,:,LM) + DMI(:,:,LM)*EVAP(:,:)/RHOE(:,:,LM)
    end if
#endif

    ! Add the topographic roughness term
    ! ----------------------------------

    if (associated(AKSODT)) then
        AKSODT = -AKS/DT
        AKSODT(:,:,1) = 0.0
    end if

    if (associated(CKSODT)) then
        CKSODT = -CKS/DT
        CKSODT(:,:,LM) = 0.0
    end if

    if (associated(AKQODT)) then
        AKQODT = -AKQ/DT
        AKQODT(:,:,1) = 0.0
    end if

    if (associated(CKQODT)) then
        CKQODT = -CKQ/DT
        CKQODT(:,:,LM) = 0.0
    end if

    if (associated(AKVODT)) AKVODT = -AKV/DT
    if (associated(CKVODT)) CKVODT = -CKV/DT

!$acc end kernels

    call cpu_time(t2)

!$acc end data

    print*,'POSTLOCK Runtime = ', t2-t1

    write(*,*) 'SUM DIFF ZPBL = ',    sum(ZPBL - ZPBL_ref),       sum(ZPBL),    sum(ZPBL_ref)
    write(*,*) 'SUM DIFF thetavs = ', sum(thetavs - thetavs_ref), sum(thetavs), sum(thetavs_ref)
    write(*,*) 'SUM DIFF thetavh = ', sum(thetavh - thetavh_ref), sum(thetavh), sum(thetavh_ref)
    write(*,*) 'SUM DIFF uv2h = ', sum(uv2h - uv2h_ref), sum(uv2h), sum(uv2h_ref)
    write(*,*) 'SUM DIFF KPBL = ', sum(KPBL - KPBL_ref), sum(KPBL), sum(KPBL_ref)
    write(*,*) 'SUM DIFF KPBLTC = ', sum(KPBLTC - KPBLTC_ref), sum(KPBLTC), sum(KPBLTC_ref)
    if(CALC_ZPBL2) then
        write(*,*) 'SUM DIFF ZPBL2 = ', sum(ZPBL2 - ZPBL2_ref), sum(ZPBL2), sum(ZPBL2_ref)
        write(*,*) 'SUM DIFF KPBL2 = ', sum(KPBL2 - KPBL2_ref), sum(KPBL2), sum(KPBL2_ref)
    endif
    write(*,*) 'SUM DIFF KPBL_SC = ', sum(KPBL_SC - KPBL_SC_ref), sum(KPBL_SC), sum(KPBL_SC_ref)
    write(*,*) 'SUM DIFF ZPBL10p = ', sum(ZPBL10p - ZPBL10p_ref), sum(ZPBL10p), sum(ZPBL10p_ref)
    write(*,*) 'SUM DIFF KPBL10p = ', sum(KPBL10p - KPBL10p_ref), sum(KPBL10p), sum(KPBL10p_ref)
    write(*,*) 'SUM DIFF TCZPBL = ', sum(TCZPBL - TCZPBL_ref), sum(TCZPBL), sum(TCZPBL_ref)
    if(DO_SHOC.eq.0) write(*,*) 'SUM DIFF ISOTROPY = ', sum(ISOTROPY - ISOTROPY_ref), sum(ISOTROPY), sum(ISOTROPY_ref)
    ! write(*,*) 'SUM DIFF tcrib = ', sum(tcrib - tcrib_ref), sum(tcrib), sum(tcrib_ref)
    write(*,*) 'SUM DIFF CKS = ', sum(CKS - CKS_ref), sum(CKS), sum(CKS_ref)
    write(*,*) 'SUM DIFF AKS = ', sum(AKS - AKS_ref), sum(AKS), sum(AKS_ref)
    write(*,*) 'SUM DIFF AKQ = ', sum(AKQ - AKQ_ref), sum(AKQ), sum(AKQ_ref)
    write(*,*) 'SUM DIFF CKQ = ', sum(CKQ - CKQ_ref), sum(CKQ), sum(CKQ_ref)
    write(*,*) 'SUM DIFF EKV = ', sum(EKV - EKV_ref), sum(EKV), sum(EKV_ref)
    write(*,*) 'SUM DIFF AKV = ', sum(AKV - AKV_ref), sum(AKV), sum(AKV_ref)
    write(*,*) 'SUM DIFF CKV = ', sum(CKV - CKV_ref), sum(CKV), sum(CKV_ref)
    write(*,*) 'SUM DIFF BKS = ', sum(BKS - BKS_ref), sum(BKS), sum(BKS_ref)
    write(*,*) 'SUM DIFF BKQ = ', sum(BKQ - BKQ_ref), sum(BKQ), sum(BKQ_ref)
    write(*,*) 'SUM DIFF BKV = ', sum(BKV - BKV_ref), sum(BKV), sum(BKV_ref)
    write(*,*) 'SUM DIFF AKSS = ', sum(AKSS - AKSS_ref), sum(AKSS), sum(AKSS_ref)
    write(*,*) 'SUM DIFF AKQQ = ', sum(AKQQ - AKQQ_ref), sum(AKQQ), sum(AKQQ_ref)
    write(*,*) 'SUM DIFF CKSS = ', sum(CKSS - CKSS_ref), sum(CKSS), sum(CKSS_ref)
    write(*,*) 'SUM DIFF CKQQ = ', sum(CKQQ - CKQQ_ref), sum(CKQQ), sum(CKQQ_ref)
    write(*,*) 'SUM DIFF CKUU = ', sum(CKUU - CKUU_ref), sum(CKUU), sum(CKUU_ref)
    write(*,*) 'SUM DIFF BKSS = ', sum(BKSS - BKSS_ref), sum(BKSS), sum(BKSS_ref)
    write(*,*) 'SUM DIFF BKQQ = ', sum(BKQQ - BKQQ_ref), sum(BKQQ), sum(BKQQ_ref)
    write(*,*) 'SUM DIFF BKUU = ', sum(BKUU - BKUU_ref), sum(BKUU), sum(BKUU_ref)
    write(*,*) 'SUM DIFF YS = ', sum(YS - YS_ref), sum(YS), sum(YS_ref)
    write(*,*) 'SUM DIFF TKH = ', sum(TKH - TKH_ref), sum(TKH), sum(TKH_ref)

    if (associated(PPBL)) write(*,*) 'SUM DIFF PPBL = ', sum(PPBL - PPBL_ref), sum(PPBL), sum(PPBL_ref)
    if (associated(ZPBLRI)) write(*,*) 'SUM DIFF ZPBLRI = ', sum(ZPBLRI - ZPBLRI_ref), sum(ZPBLRI), sum(ZPBLRI_ref)
    if (associated(ZPBLRI2)) write(*,*) 'SUM DIFF ZPBLRI2 = ', sum(ZPBLRI2 - ZPBLRI2_ref), sum(ZPBLRI2), sum(ZPBLRI2_ref)
    if (associated(ZPBLTHV)) write(*,*) 'SUM DIFF ZPBLTHV = ', sum(ZPBLTHV - ZPBLTHV_ref), sum(ZPBLTHV), sum(ZPBLTHV_ref)
    if (associated(ZPBLHTKE)) write(*,*) 'SUM DIFF ZPBLHTKE = ', sum(ZPBLHTKE - ZPBLHTKE_ref), sum(ZPBLHTKE), sum(ZPBLHTKE_ref)
    if (associated(TKE)) write(*,*) 'SUM DIFF TKE = ', sum(TKE - TKE_ref), sum(TKE), sum(TKE_ref)
    if (associated(AKSODT)) write(*,*) 'SUM DIFF AKSODT = ', sum(AKSODT - AKSODT_ref), sum(AKSODT), sum(AKSODT_ref)
    if (associated(CKSODT)) write(*,*) 'SUM DIFF CKSODT = ', sum(CKSODT - CKSODT_ref), sum(CKSODT), sum(CKSODT_ref)
    if (associated(AKQODT)) write(*,*) 'SUM DIFF AKQODT = ', sum(AKQODT - AKQODT_ref), sum(AKQODT), sum(AKQODT_ref)
    if (associated(CKQODT)) write(*,*) 'SUM DIFF CKQODT = ', sum(CKQODT - CKQODT_ref), sum(CKQODT), sum(CKQODT_ref)
    if (associated(AKVODT)) write(*,*) 'SUM DIFF AKVODT = ', sum(AKVODT - AKVODT_ref), sum(AKVODT), sum(AKVODT_ref)
    if (associated(CKVODT)) write(*,*) 'SUM DIFF CKVODT = ', sum(CKVODT - CKVODT_ref), sum(CKVODT), sum(CKVODT_ref)

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