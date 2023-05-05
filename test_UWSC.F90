program test_UWSC

    use uwshcu
    use SHLWPARAMS

    implicit none

    character (len=100) :: dirName

    integer   :: idim, k0, ncnst, IM, JM, fileID, dotransport, ii, jj
    real      :: dt, thlsrc_pert
    real :: start, end

    real, dimension(:,:,:), allocatable :: pifc0_inv     !  Environmental pressure at the interfaces [ Pa ]
    real, dimension(:,:,:), allocatable :: zifc0_inv     !  Environmental height at the interfaces   [ m ]
    real, dimension(:,:,:), allocatable :: exnifc0_inv   !  Exner function at the interfaces
    real, dimension(:,:,:), allocatable :: pmid0_inv     !  Environmental pressure at the layer mid-point [ Pa ]
    real, dimension(:,:,:), allocatable :: zmid0_inv     !  Environmental height at the layer mid-point [ m ]
    real, dimension(:,:,:), allocatable :: exnmid0_inv   !  Exner function at the layer mid-point
    real, dimension(:,:,:), allocatable :: dp0_inv       !  Environmental layer pressure thickness [ Pa ] > 0.

    real, dimension(:,:,:), allocatable :: u0_inv          !  Environmental zonal wind [ m/s ]
    real, dimension(:,:,:), allocatable :: v0_inv          !  Environmental meridional wind [ m/s ]
    real, dimension(:,:,:), allocatable :: qv0_inv         !  Environmental water vapor specific humidity [ kg/kg ]
    real, dimension(:,:,:), allocatable :: ql0_inv         !  Environmental liquid water specific humidity [ kg/kg ]
    real, dimension(:,:,:), allocatable :: qi0_inv         !  Environmental ice specific humidity [ kg/kg ]
    real, dimension(:,:,:), allocatable :: th0_inv         !  Environmental temperature [ K ]
    real, dimension(:,:,:), allocatable :: tke_inv         !  Turbulent kinetic energy at the interfaces [ m2/s2 ]
                                                         !  at the previous time step [ fraction ]

    real, dimension(:,:,:,:), allocatable :: tr0_inv   !  Environmental tracers [ #, kg/kg ]

    real, dimension(:,:,:), allocatable :: umf_inv     !  Updraft mass flux at interfaces [kg/m2/s]
    real, dimension(:,:,:), allocatable :: qvten_inv   !  Tendency of water vapor specific humidity [ kg/kg/s ]
    real, dimension(:,:,:), allocatable :: qlten_inv   !  Tendency of liquid water specific humidity [ kg/kg/s ]
    real, dimension(:,:,:), allocatable :: qiten_inv   !  Tendency of ice specific humidity [ kg/kg/s ]
    real, dimension(:,:,:), allocatable :: thten_inv   !  Tendency of potential temperature [ K/s ]
    real, dimension(:,:,:), allocatable :: uten_inv    !  Tendency of zonal wind [ m/s2 ]
    real, dimension(:,:,:), allocatable :: vten_inv    !  Tendency of meridional wind [ m/s2 ]

    real, dimension(:,:,:), allocatable :: qrten_inv   !  Tendency of rain water specific humidity [ kg/kg/s ]
    real, dimension(:,:,:), allocatable :: qsten_inv   !  Tendency of snow specific humidity [ kg/kg/s ]
    real, dimension(:,:,:), allocatable :: cufrc_inv   !  Shallow cumulus cloud fraction at the layer mid-point [ fraction ]
    real, dimension(:,:,:), allocatable :: qcu_inv     !  Liquid+ice specific humidity within cumulus updraft [ kg/kg ]
    real, dimension(:,:,:), allocatable :: qlu_inv     !  Liquid water specific humidity within cumulus updraft [ kg/kg ]
    real, dimension(:,:,:), allocatable :: qiu_inv     !  Ice specific humidity within cumulus updraft [ kg/kg ]
    real, dimension(:,:,:), allocatable :: qc_inv      !  Tendency of cumulus condensate detrained into the environment [ kg/kg/s ]
    real, dimension(:,:,:), allocatable :: wu_inv
    real, dimension(:,:,:), allocatable :: qtu_inv
    real, dimension(:,:,:), allocatable :: thlu_inv
    real, dimension(:,:,:), allocatable :: thvu_inv
    real, dimension(:,:,:), allocatable :: uu_inv
    real, dimension(:,:,:), allocatable :: vu_inv
    real, dimension(:,:,:), allocatable :: fer_inv
    real, dimension(:,:,:), allocatable :: fdr_inv
    real, dimension(:,:,:), allocatable :: xc_inv
    real, dimension(:,:,:), allocatable :: qldet_inv
    real, dimension(:,:,:), allocatable :: qidet_inv
    real, dimension(:,:,:), allocatable :: qlsub_inv
    real, dimension(:,:,:), allocatable :: qisub_inv
    real, dimension(:,:,:), allocatable :: ndrop_inv
    real, dimension(:,:,:), allocatable :: nice_inv

    real, dimension(:,:), allocatable :: cbmf             !  Cumulus base mass flux [ kg/m2/s ]
    real, dimension(:,:), allocatable :: cnt_inv          !  Cumulus top  interface index, cnt = kpen [ no ]
    real, dimension(:,:), allocatable :: cnb_inv          !  Cumulus base interface index, cnb = krel - 1 [ no ]
    real, dimension(:,:), allocatable :: cin
    real, dimension(:,:), allocatable :: plcl
    real, dimension(:,:), allocatable :: plfc
    real, dimension(:,:), allocatable :: pinv
    real, dimension(:,:), allocatable :: prel
    real, dimension(:,:), allocatable :: pbup
    real, dimension(:,:), allocatable :: wlcl
    real, dimension(:,:), allocatable :: qtsrc
    real, dimension(:,:), allocatable :: thlsrc
    real, dimension(:,:), allocatable :: thvlsrc
    real, dimension(:,:), allocatable :: tkeavg
    real, dimension(:,:), allocatable :: cldtop
    real, dimension(:,:), allocatable :: kpbl_inv           !  Height of PBL [ m ]
    real, dimension(:,:), allocatable :: cush               !  Convective scale height [m]
    
    real, dimension(:,:), allocatable :: cush_test
    real, dimension(:,:,:,:), allocatable :: tr0_inv_test
    real, dimension(:,:,:), allocatable :: umf_inv_test

    type(shlwparam_type) :: shlwparams_obj

    character *100 :: BUFFER

    ! read(*,*) IM, JM, k0, ncnst, and directory of initialized and comparison data
    call getarg(1,BUFFER)
    read(BUFFER,*) IM
    call getarg(2,BUFFER)
    read(BUFFER,*) JM
    call getarg(3,BUFFER)
    read(BUFFER,*) k0
    call getarg(4,BUFFER)
    read(BUFFER,*) ncnst
    call getarg(5,BUFFER)
    read(BUFFER,'(A)') dirName
!!! 48 48 72 43 ../Datasets/Shallow/C48-L72/
    ! IM = 48
    ! JM = 48
    ! k0 = 72
    ! ncnst = 43
    ! dirName = "../Datasets/Shallow/C48-L72/"

    write(*,*) 'IM = ', IM
    write(*,*) 'JM = ', JM
    write(*,*) 'k0 = ', k0
    write(*,*) 'ncnst = ', ncnst
    write(*,*) 'dirName = ', trim(dirName)

    allocate(pifc0_inv(IM, JM, k0+1))
    allocate(zifc0_inv(IM, JM, k0+1))
    allocate(exnifc0_inv(IM, JM, k0+1))
    allocate(pmid0_inv(IM, JM, k0))
    allocate(zmid0_inv(IM, JM, k0))
    allocate(exnmid0_inv(IM, JM, k0))
    allocate(dp0_inv(IM, JM, k0))

    allocate(u0_inv(IM, JM, k0))
    allocate(v0_inv(IM, JM, k0))
    allocate(qv0_inv(IM, JM, k0))
    allocate(ql0_inv(IM, JM, k0))
    allocate(qi0_inv(IM, JM, k0))
    allocate(th0_inv(IM, JM, k0))
    allocate(tke_inv(IM, JM, k0+1))

    allocate(tr0_inv(IM, JM, k0, ncnst))
    allocate(tr0_inv_test(IM, JM, k0, ncnst))

    allocate(umf_inv(IM, JM, k0+1))
    allocate(umf_inv_test(Im, JM, k0+1))

    allocate(qvten_inv(IM, JM, k0))
    allocate(qlten_inv(IM, JM, k0))
    allocate(qiten_inv(IM, JM, k0))
    allocate(thten_inv(IM, JM, k0))
    allocate(uten_inv(IM, JM, k0))
    allocate(vten_inv(IM, JM, k0))

    allocate(qrten_inv(IM, JM, k0))
    allocate(qsten_inv(IM, JM, k0))
    allocate(cufrc_inv(IM, JM, k0))
    allocate(qcu_inv(IM, JM, k0))
    allocate(qlu_inv(IM, JM, k0))
    allocate(qiu_inv(IM, JM, k0))
    allocate(qc_inv(IM, JM, k0))
    allocate(wu_inv(IM, JM, k0+1))
    allocate(qtu_inv(IM, JM, k0+1))
    allocate(thlu_inv(IM, JM, k0+1))
    allocate(thvu_inv(IM, JM, k0+1))
    allocate(uu_inv(IM, JM, k0+1))
    allocate(vu_inv(IM, JM, k0+1))
    allocate(fer_inv(IM, JM, k0))
    allocate(fdr_inv(IM, JM, k0))
    allocate(xc_inv(IM, JM, k0))
    allocate(qldet_inv(IM, JM, k0))
    allocate(qidet_inv(IM, JM, k0))
    allocate(qlsub_inv(IM, JM, k0))
    allocate(qisub_inv(IM, JM, k0))
    allocate(ndrop_inv(IM, JM, k0))
    allocate(nice_inv(IM, JM, k0))

    allocate(cbmf(IM, JM))
    allocate(cnt_inv(IM, JM))
    allocate(cnb_inv(IM, JM))
    allocate(cin(IM, JM))
    allocate(plcl(IM, JM))
    allocate(plfc(IM, JM))
    allocate(pinv(IM, JM))
    allocate(prel(IM, JM))
    allocate(pbup(IM, JM))
    allocate(wlcl(IM, JM))
    allocate(qtsrc(IM, JM))
    allocate(thlsrc(IM, JM))
    allocate(thvlsrc(IM, JM))
    allocate(tkeavg(IM, JM))
    allocate(cldtop(IM, JM))
    allocate(kpbl_inv(IM, JM))
    allocate(cush(IM, JM))
    allocate(cush_test(IM, JM))

    open(newunit=fileID, file=trim(dirName) // "/plo.dat", status='old', form="unformatted", action="read")
    read(fileID) pmid0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/zlo.dat", status='old', form="unformatted", action="read")
    read(fileID) zmid0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/pk.dat", status='old', form="unformatted", action="read")
    read(fileID) exnmid0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/ple.dat", status='old', form="unformatted", action="read")
    read(fileID) pifc0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/zle.dat", status='old', form="unformatted", action="read")
    read(fileID) zifc0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/pke.dat", status='old', form="unformatted", action="read")
    read(fileID) exnifc0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/dp.dat", status='old', form="unformatted", action="read")
    read(fileID) dp0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/u1.dat", status='old', form="unformatted", action="read")
    read(fileID) u0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/v1.dat", status='old', form="unformatted", action="read")
    read(fileID) v0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/q1.dat", status='old', form="unformatted", action="read")
    read(fileID) qv0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/qlls.dat", status='old', form="unformatted", action="read")
    read(fileID) ql0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/qils.dat", status='old', form="unformatted", action="read")
    read(fileID) qi0_inv
    close(fileID)
    
    open(newunit=fileID, file=trim(dirName) // "/th1.dat", status='old', form="unformatted", action="read")
    read(fileID) th0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/tke.dat", status='old', form="unformatted", action="read")
    read(fileID) tke_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/kpblsc.dat", status='old', form="unformatted", action="read")
    read(fileID) kpbl_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/cush_in.dat", status='old', form="unformatted", action="read")
    read(fileID) cush
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/xho_in.dat", status='old', form="unformatted", action="read")
    read(fileID) tr0_inv
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/Outputs/cush_out.dat", status='old', form="unformatted", action="read")
    read(fileID) cush_test
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/Outputs/xho_out.dat", status='old', form="unformatted", action="read")
    read(fileID) tr0_inv_test
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/Outputs/umf_sc.dat", status='old', form="unformatted", action="read")
    read(fileID) umf_inv_test
    close(fileID)

    umf_inv = -1.0
    qvten_inv = -1.0
    qlten_inv = -1.0
    qiten_inv = -1.0
    thten_inv = -1.0
    uten_inv = -1.0
    vten_inv = -1.0

    qrten_inv = -1.0
    qsten_inv = -1.0
    cufrc_inv = -1.0
    qcu_inv = -1.0
    qlu_inv = -1.0
    qiu_inv = -1.0
    qc_inv = -1.0
    wu_inv = -1.0
    qtu_inv = -1.0
    thlu_inv = -1.0
    thvu_inv = -1.0
    uu_inv = -1.0
    vu_inv = -1.0
    fer_inv = -1.0
    fdr_inv = -1.0
    xc_inv = -1.0
    qldet_inv = -1.0
    qidet_inv = -1.0
    qlsub_inv = -1.0
    qisub_inv = -1.0
    ndrop_inv = -1.0
    nice_inv = -1.0

    cbmf = -1.0
    cnt_inv = -1.0
    cnb_inv = -1.0
    cin = -1.0
    plcl = -1.0
    plfc = -1.0
    pinv = -1.0
    prel = -1.0
    pbup = -1.0
    wlcl = -1.0
    qtsrc = -1.0
    thlsrc = -1.0
    thvlsrc = -1.0
    tkeavg = -1.0
    cldtop = -1.0

    ! Set SHLWPARAMS_OBJ values to "defaults" from GEOS_MoistGridComp.F90
    shlwparams_obj%niter_xc = 2
    shlwparams_obj%iter_cin = 2
    shlwparams_obj%use_CINcin = 1
    shlwparams_obj%use_self_detrain = 0
    shlwparams_obj%use_momenflx = 1
    shlwparams_obj%use_cumpenent = 1
    shlwparams_obj%scverbose = 0
    shlwparams_obj%rpen = 3.0
    shlwparams_obj%rle = 0.1
    shlwparams_obj%rkm = 12.0
    shlwparams_obj%rkfre = 1.0
    shlwparams_obj%rmaxfrac = 0.1
    shlwparams_obj%mumin1 = 0.906
    shlwparams_obj%rbuoy = 1.0
    shlwparams_obj%rdrag = 1.0
    shlwparams_obj%epsvarw = 5.0e-4
    shlwparams_obj%pgfc = 0.7
    shlwparams_obj%criqc = 1.0e-3
    shlwparams_obj%kevp = 2.0e-6
    shlwparams_obj%rdrop = 8.0e-6
    shlwparams_obj%frc_rasn = 1.0

    idim = IM * JM

    dotransport = 1
    thlsrc_pert = 0.0
    dt = 450.0
    pmid0_inv = pmid0_inv*100
    write(*,*) 'Calling compute_uwshcu_inv...'
!$acc data copyin(pmid0_inv,zmid0_inv,exnmid0_inv,pifc0_inv,zifc0_inv,&
!$acc      exnifc0_inv,dp0_inv,u0_inv,v0_inv,qv0_inv,ql0_inv,qi0_inv,&
!$acc      th0_inv,tke_inv,kpbl_inv)&
!$acc      copy(cush,tr0_inv) &
!$acc      copyout(umf_inv,qvten_inv,qlten_inv,qiten_inv,thten_inv,uten_inv,&
!$acc      vten_inv,qrten_inv,qsten_inv,cufrc_inv,fer_inv,fdr_inv, &
!$acc      qldet_inv,qidet_inv,qlsub_inv,qisub_inv,ndrop_inv,nice_inv)
    call cpu_time(start)
    call compute_uwshcu_inv(idim, k0, ncnst, dt, &
    pmid0_inv, zmid0_inv, exnmid0_inv, pifc0_inv, zifc0_inv, exnifc0_inv, dp0_inv, &
    u0_inv, v0_inv, qv0_inv, ql0_inv, qi0_inv, th0_inv, tke_inv, kpbl_inv, &
    thlsrc_pert, &
    cush, tr0_inv, &
    umf_inv, qvten_inv, qlten_inv, qiten_inv, &
    thten_inv, uten_inv, vten_inv, qrten_inv, &
    qsten_inv, cufrc_inv, fer_inv, fdr_inv, &
    qldet_inv, qidet_inv, qlsub_inv, qisub_inv, &
    ndrop_inv, nice_inv, &
#ifdef UWDIAG
         qcu_inv, qlu_inv, qiu_inv, cbmf, qc_inv,                   & ! DIAGNOSTIC ONLY
         cnt_inv, cnb_inv, cin, plcl, plfc, pinv, prel, pbup,       &
         wlcl, qtsrc, thlsrc, thvlsrc, tkeavg, cldtop, wu_inv,      &
         qtu_inv, thlu_inv, thvu_inv, uu_inv, vu_inv, xc_inv,       &
#endif
    dotransport, shlwparams_obj)
    call cpu_time(end)
!$acc end data
    write(*, *) 'done.'
    write(*, '(a, f10.4, a)') 'Time taken: ', end-start, 's.'

    !write(100, *) cush_test
    !write(200, *) cush
    write(*,*) "cush Difference Sum = ", sum(cush - cush_test)
    print*,'sum(cush) = ', sum(cush)
    print*,'sum(cust_test) = ', sum(cush_test)
    write(*,*) "tr0_inv Difference Sum = ", sum(tr0_inv - tr0_inv_test)
    write(*,*) "umf_inv Difference Sum = ", sum(umf_inv - umf_inv_test)

end program
