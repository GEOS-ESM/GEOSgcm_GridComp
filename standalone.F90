program run_shoc_standalone

    use SHOCPARAMS
    use shoc

    implicit none

    integer :: IM, JM, LM, fileID
    real    :: DT
    character*100 :: dirName, rank_str

    real, dimension(:,:,:), allocatable :: DMI, DMI_I, PLO, ZLE, Z, U, V, OMEGA, T, Q
    real, dimension(:,:,:), allocatable :: QI, QL, QPI, QPL, QA, WTHV2, PRANDTLSHOC
    real, dimension(:,:,:), allocatable :: TKESHOC, TKH, ISOTROPY
    real, dimension(:,:,:), allocatable :: TKESHOC_ref, TKH_ref, ISOTROPY_ref

    real, dimension(:,:,:), pointer :: TKEDISS, TKEBUOY, TKESHEAR, TKETRANS
    real, dimension(:,:,:), pointer :: LSHOC, LSHOC1, LSHOC2, LSHOC3, LSHOC_CLD, LSHOC_CLR
    real, dimension(:,:,:), pointer :: BRUNTSHOC, SHEARSHOC

    real :: ts, te

    type(SHOCPARAMS_TYPE) :: SHOCPARAMS_var

    character(len=12), dimension(:), allocatable :: args
    integer :: num_args

    dirName = './c180-180-180-72'

    num_args = command_argument_count()
    if (num_args.eq.0) then
        print*, 'Need to include rank argument'
        call exit(1)
    endif
    allocate(args(num_args)) 

    call get_command_argument(1,args(1))
    rank_str = trim(args(1))

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

    TKEDISS => null()
    TKEBUOY => null()
    TKESHEAR => null()
    TKETRANS => null()
    LSHOC => null()
    LSHOC1 => null()
    LSHOC2 => null()
    LSHOC3 => null()
    LSHOC_CLD => null()
    LSHOC_CLR => null()
    BRUNTSHOC => null()
    SHEARSHOC => null()

    SHOCPARAMS_var%LAMBDA   = 4.0e-02
    SHOCPARAMS_var%TSCALE   = 400.0
    SHOCPARAMS_var%VONK     = 0.40
    SHOCPARAMS_var%CKVAL    = 0.10
    SHOCPARAMS_var%CEFAC    = 1.0
    SHOCPARAMS_var%CESFAC   = 4.0
    SHOCPARAMS_var%CLDLEN   = 0
    SHOCPARAMS_var%SUS12LEN = 1
    SHOCPARAMS_var%BUOYOPT  = 2

!   Inputs arrays
    allocate(DMI(IM, JM, LM))
    allocate(DMI_I(IM, JM, LM))
    allocate(PLO(IM, JM, LM))
    allocate(ZLE(IM, JM, 0:LM))
    allocate(Z(IM, JM, LM))
    allocate(U(IM, JM, LM))
    allocate(V(IM, JM, LM))
    allocate(OMEGA(IM, JM, LM))
    allocate(T(IM, JM, LM))
    allocate(Q(IM, JM, LM))
    allocate(QI(IM, JM, LM))
    allocate(QL(IM, JM, LM))
    allocate(QPI(IM, JM, LM))
    allocate(QPL(IM, JM, LM))
    allocate(QA(IM, JM, LM))
    allocate(WTHV2(IM, JM, LM))
    allocate(PRANDTLSHOC(IM, JM, LM))

!   Output arrays
    allocate(TKESHOC(IM, JM, LM))
    allocate(TKH(IM, JM, 0:LM))
    allocate(ISOTROPY(IM, JM, LM))

    allocate(TKESHOC_ref(IM, JM, LM))
    allocate(TKH_ref(IM, JM, 0:LM))
    allocate(ISOTROPY_ref(IM, JM, LM))

!   Reading Input Arrays

    open(newunit=fileID, file=trim(dirName) // "/DMI_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) DMI
    close(fileID)
    !write(*,*) 'sum(DMI) = ', sum(DMI)

    open(newunit=fileID, file=trim(dirName) // "/PLO_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) PLO
    close(fileID)
    !write(*,*) 'sum(PLO) = ', sum(PLO)

    open(newunit=fileID, file=trim(dirName) // "/ZLE_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) ZLE
    close(fileID)
    !write(*,*) 'sum(ZLE) = ', sum(ZLE)

    open(newunit=fileID, file=trim(dirName) // "/Z_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) Z
    close(fileID)
    !write(*,*) 'sum(Z) = ', sum(Z)

    open(newunit=fileID, file=trim(dirName) // "/U_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) U
    close(fileID)
    !write(*,*) 'sum(U) = ', sum(U)

    open(newunit=fileID, file=trim(dirName) // "/V_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) V
    close(fileID)
    !write(*,*) 'sum(V) = ', sum(V)

    open(newunit=fileID, file=trim(dirName) // "/OMEGA_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) OMEGA
    close(fileID)
    !write(*,*) 'sum(OMEGA) = ', sum(OMEGA)

    open(newunit=fileID, file=trim(dirName) // "/T_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) T
    close(fileID)
    !write(*,*) 'sum(T) = ', sum(T)

    open(newunit=fileID, file=trim(dirName) // "/Q_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) Q
    close(fileID)
    !write(*,*) 'sum(Q) = ', sum(Q)

    open(newunit=fileID, file=trim(dirName) // "/QI_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QI
    close(fileID)
    !write(*,*) 'sum(QI) = ', sum(QI)

    open(newunit=fileID, file=trim(dirName) // "/QL_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QL
    close(fileID)
    !write(*,*) 'sum(QL) = ', sum(QL)

    open(newunit=fileID, file=trim(dirName) // "/QPI_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QPI
    close(fileID)
    !write(*,*) 'sum(QPI) = ', sum(QPI)

    open(newunit=fileID, file=trim(dirName) // "/QPL_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QPL
    close(fileID)
    !write(*,*) 'sum(QPL) = ', sum(QPL)

    open(newunit=fileID, file=trim(dirName) // "/QA_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QA
    close(fileID)
    !write(*,*) 'sum(QA) = ', sum(QA)

    open(newunit=fileID, file=trim(dirName) // "/WTHV2_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) WTHV2
    close(fileID)
    !write(*,*) 'sum(WTHV2) = ', sum(WTHV2)

    open(newunit=fileID, file=trim(dirName) // "/PRANDTLSHOC_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) PRANDTLSHOC
    close(fileID)
    !write(*,*) 'sum(PRANDTLSCHO) = ', sum(PRANDTLSHOC)

    open(newunit=fileID, file=trim(dirName) // "/TKESHOC_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) TKESHOC
    close(fileID)
    !write(*,*) 'sum(TKESHOC) = ', sum(TKESHOC)

    open(newunit=fileID, file=trim(dirName) // "/TKH_" // trim(rank_str) // ".in", status='old', form="unformatted", action="read")
    read(fileID) TKH
    close(fileID)
    !write(*,*) 'sum(TKH) = ', sum(TKH)

    ! Reading output reference solution arrays

    open(newunit=fileID, file=trim(dirName) // "/TKESHOC_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) TKESHOC_ref
    close(fileID)
    !write(*,*) 'sum(TKESHOC_ref) = ', sum(TKESHOC_ref)

    open(newunit=fileID, file=trim(dirName) // "/TKH_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) TKH_ref
    close(fileID)
    !write(*,*) 'sum(TKH_ref) = ', sum(TKH_ref)

    open(newunit=fileID, file=trim(dirName) // "/ISOTROPY_" // trim(rank_str) // ".out", status='old', form="unformatted", action="read")
    read(fileID) ISOTROPY_ref
    close(fileID)
    !write(*,*) 'sum(ISOTROPY_ref) = ', sum(ISOTROPY_ref)

   DMI_I(:,:,1:LM) = DT/DMI(:,:,1:LM)

!$acc data copyin(DMI_I,PLO,ZLE,Z,U,V,OMEGA,T,Q,QI,QL,QPI,QPL,QA,WTHV2) &
!$acc      copyout(ISOTROPY) &
!$acc      copy(PRANDTLSHOC,TKESHOC,TKH,TKEDISS,TKEBUOY,TKESHEAR,TKETRANS,LSHOC,LSHOC_CLR,&
!$acc           LSHOC_CLD,LSHOC1,LSHOC2,LSHOC3,BRUNTSHOC,SHEARSHOC,SHOCPARAMS_var,&
!$acc           SHOCPARAMS_var%LAMBDA,SHOCPARAMS_var%TSCALE,SHOCPARAMS_var%VONK,&
!$acc           SHOCPARAMS_var%CKVAL,SHOCPARAMS_var%CEFAC,SHOCPARAMS_var%CESFAC,&
!$acc           SHOCPARAMS_var%CLDLEN,SHOCPARAMS_var%SUS12LEN,SHOCPARAMS_var%BUOYOPT)

    call cpu_time(ts)

    call RUN_SHOC( IM, JM, LM, LM+1, DT, &  
                       !== Inputs ==
                       DMI_I(:,:,1:LM),      &
                       PLO(:,:,1:LM),         &
                       ZLE(:,:,0:LM),         &
                       Z(:,:,1:LM),           &
                       U(:,:,1:LM),           &
                       V(:,:,1:LM),           &
                       OMEGA(:,:,1:LM),       &
                       T(:,:,1:LM),           &
                       Q(:,:,1:LM),           &
                       QI(:,:,1:LM),          &
                       QL(:,:,1:LM),          &
                       QPI(:,:,1:LM),         &
                       QPL(:,:,1:LM),         &
                       QA(:,:,1:LM),          &
                       WTHV2(:,:,1:LM),       &
                       PRANDTLSHOC(:,:,1:LM), &
                       !== Input-Outputs ==
                       TKESHOC(:,:,1:LM),     &
                       TKH(:,:,1:LM),         &
                       !== Outputs ==
                       ISOTROPY(:,:,1:LM),    &
                       !== Diagnostics ==  ! not used elsewhere
                       TKEDISS,               &
                       TKEBUOY,               &
                       TKESHEAR,              &
                       TKETRANS,              &
                       LSHOC,                 &
                       LSHOC_CLR,             &
                       LSHOC_CLD,             &
                       LSHOC1,                &
                       LSHOC2,                &
                       LSHOC3,                &
                       BRUNTSHOC,             &
                       SHEARSHOC,             &
                       !== Tuning params ==
                       SHOCPARAMS_var )

    call cpu_time(te)
!$acc end data


    write(*,*) 'Sum diff abs TKESHOC = ', sum(TKESHOC - TKESHOC_ref), sum(abs(TKESHOC)), sum(abs(TKESHOC_ref))
    write(*,*) 'Sum diff abs TKEH = ', sum(TKH - TKH_ref), sum(abs(TKH)), sum(abs(TKH_ref))
    write(*,*) 'Sum diff abs ISOTROPY = ', sum(ISOTROPY - ISOTROPY_ref), sum(abs(ISOTROPY)), sum(abs(ISOTROPY_ref))
    print*,'****'
    print*,'RUN_SHOC runtime = ', te-ts

end program

! NASA Docket No. GSC-15,354-1, and identified as "GEOS-5 GCM Modeling Software”

! “Copyright © 2008 United States Government as represented by the Administrator
! of the National Aeronautics and Space Administration. All Rights Reserved.”

! Licensed under the Apache License, Version 2.0 (the "License"); you may not use
! this file except in compliance with the License. You may obtain a copy of the
! License at

! http://www.apache.org/licenses/LICENSE-2.0