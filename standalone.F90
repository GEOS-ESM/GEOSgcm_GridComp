program ncar_gwd_standalone

    use gw_drag_ncar, only : gw_intr_ncar
    use gw_convect,   only : BeresSourceDesc, gw_beres_init
    use gw_common,    only : GWBand, gw_common_init, gw_newtonian_set
    use gw_oro,       only : gw_oro_init
    use gw_rdg,       only : gw_rdg_init
    use MAPL_ConstantsMod

    implicit none

    integer               :: IM, JM, LM, fileID, status, PGWV, RC, nrdg
    integer               :: NCAR_BKG_PGWV, NCAR_ORO_PGWV, NCAR_NRDG
    real                  :: DT, bgstressmax, effgwbkg, effgworo
    real                  :: t_start, t_end
    real                  :: NCAR_PRNDL, NCAR_QBO_HDEPTH_SCALING
    real                  :: NCAR_HR_CF, NCAR_BKG_GW_DC, NCAR_BKG_FCRIT2
    real                  :: NCAR_BKG_WAVELENGTH, NCAR_DC_BERES_SRC_LEVEL
    real                  :: NCAR_SC_BERES_SRC_LEVEL
    real                  :: NCAR_ET_TAUBGND, NCAR_BKG_TNDMAX
    real                  :: NCAR_ORO_GW_DC, NCAR_ORO_FCRIT2, NCAR_ORO_WAVELENGTH
    real                  :: NCAR_ORO_SOUTH_FAC, NCAR_ORO_TNDMAX, NCAR_EFFGWORO
    character*100         :: Errstring, dirName, rank
    logical               :: NCAR_TAU_TOP_ZERO, NCAR_DC_BERES, NCAR_SC_BERES

    type(BeresSourceDesc) :: beres_dc_desc, beres_sc_desc
    type(GWBand)          :: beres_band, oro_band, rdg_band

    real, dimension(:,:,:), allocatable :: PLE, T, Q, U, V, HT_dpc, ZI
    real, dimension(:,:,:), allocatable :: PMID, PDEL, RPDEL, PILN, ZM, PMLN
    real, dimension(:,:,:), allocatable :: DUDT_GWD, DVDT_GWD, DTDT_GWD, DUDT_ORG, DVDT_ORG, DTDT_ORG
    real, dimension(:,:,:), allocatable :: TAUXO_3D, TAUYO_3D, TAUXB_3D, TAUYB_3D, FEO_3D, FEB_3D, FEPO_3D, FEPB_3D
    real, dimension(:,:,:), allocatable :: DUBKGSRC, DVBKGSRC, DTBKGSRC
    real, dimension(:,:,:), allocatable :: ANGLL, HWDTH, KWVRDG, CLNGT, EFFRDG
    real, dimension(:,:),   allocatable :: LATS, SGH, TAUXO_TMP, TAUYO_TMP, TAUXB_TMP, TAUYB_TMP
    real, dimension(:,:),   allocatable :: GBXAR_TMP, GBXAR
    real, dimension(:),     allocatable :: PREF, alpha

    real, dimension(:,:,:), allocatable :: DUDT_GWD_ref, DVDT_GWD_ref, DTDT_GWD_ref, DUDT_ORG_ref, DVDT_ORG_ref, DTDT_ORG_ref
    real, dimension(:,:,:), allocatable :: TAUXO_3D_ref, TAUYO_3D_ref, TAUXB_3D_ref, TAUYB_3D_ref
    real, dimension(:,:,:), allocatable :: FEO_3D_ref, FEB_3D_ref, FEPO_3D_ref, FEPB_3D_ref
    real, dimension(:,:,:), allocatable :: DUBKGSRC_ref, DVBKGSRC_ref, DTBKGSRC_ref
    real, dimension(:,:),   allocatable :: TAUXO_TMP_ref, TAUYO_TMP_ref, TAUXB_TMP_ref, TAUYB_TMP_ref
    real, dimension(:),     allocatable :: alpha_ref


    real, dimension(:,:,:), allocatable :: PDEL_ref, PILN_ref, RPDEL_ref, PMID_ref, PMLN_ref, ZI_ref, ZM_ref
    real, dimension(:,:,:), allocatable :: ANGLL_ref, KWVRDG_ref, EFFRDG_ref
    real, dimension(:,:),   allocatable :: GBXAR_TMP_ref
    dirName = './new_c180_data/ncar_gwd'
    
    if(command_argument_count().ne.1) then
        print*, 'Missing arguments : <executable> <rank>'
        call exit(1)
    else
        call get_command_argument(1, rank)
    endif
    
    if (dirName == './data-c12-6-6-72') then
    ! C12 Array Sizes
        IM = 6
        JM = 6
        LM = 72
        DT = 900.0
    else if (dirName == './data-c90-90-90-72') then
        IM = 90
        JM = 90
        LM = 72
        DT = 450.0
    else if (dirName == './data-c180-180-180-72') then
        IM = 180
        JM = 180
        LM = 72
        DT = 450.0
    else if (dirName == './new_c180_data/ncar_gwd') then
        IM = 180
        JM = 180
        LM = 72
        DT = 450.0
    endif

    PGWV = 4

    bgstressmax = 0.900
    effgwbkg    = 0.125
    effgworo    = 0.250

    NCAR_TAU_TOP_ZERO = .TRUE.

    allocate(LATS  (IM, JM))

    open(newunit=fileID, file=trim(dirName) // "/NCAR_PRNDL_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_PRNDL
    ! print*,'NCAR_PRNDL = ', NCAR_PRNDL
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_QBO_HDEPTH_SCALING_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_QBO_HDEPTH_SCALING
    ! print*,'NCAR_QBO_HDEPTH_SCALING = ',NCAR_QBO_HDEPTH_SCALING
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_HR_CF_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_HR_CF
    ! print*,'NCAR_HR_CF = ', NCAR_HR_CF
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_BKG_PGWV_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_BKG_PGWV
    ! print*,'NCAR_BKG_PGWV = ', NCAR_BKG_PGWV
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_BKG_GW_DC_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_BKG_GW_DC
    ! print*,'NCAR_BKG_GW_DC = ', NCAR_BKG_GW_DC
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_BKG_FCRIT2_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_BKG_FCRIT2
    ! print*,'NCAR_BKG_FCRIT2 = ', NCAR_BKG_FCRIT2
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_BKG_WAVELENGTH_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_BKG_WAVELENGTH
    ! print*,'NCAR_BKG_WAVELENGTH = ', NCAR_BKG_WAVELENGTH
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_DC_BERES_SRC_LEVEL_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_DC_BERES_SRC_LEVEL
    ! print*,'NCAR_DC_BERES_SRC_LEVEL = ', NCAR_DC_BERES_SRC_LEVEL
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_ET_TAUBGND_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_ET_TAUBGND
    ! print*,'NCAR_ET_TAUBGND = ', NCAR_ET_TAUBGND
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_BKG_TNDMAX_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_BKG_TNDMAX
    ! print*,'NCAR_BKG_TNDMAX = ', NCAR_BKG_TNDMAX
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_DC_BERES_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_DC_BERES
    ! print*,'NCAR_DC_BERES = ', NCAR_DC_BERES
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/LATS_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) LATS
    ! print*,'sum(LATS) = ', sum(abs(LATS))
    close(fileID)

    call gw_common_init( NCAR_TAU_TOP_ZERO , 1 , &
                        MAPL_GRAV , &
                        MAPL_RGAS , &
                        MAPL_CP , &
                        NCAR_PRNDL, NCAR_QBO_HDEPTH_SCALING, NCAR_HR_CF, ERRstring )

    call gw_beres_init( trim(dirName) // "/newmfspectra40_dc25.nc" ,  &
        beres_band, &
        beres_dc_desc, &
        NCAR_BKG_PGWV, NCAR_BKG_GW_DC, NCAR_BKG_FCRIT2, &
        NCAR_BKG_WAVELENGTH, NCAR_DC_BERES_SRC_LEVEL, &
        1000.0, .TRUE., NCAR_ET_TAUBGND, NCAR_BKG_TNDMAX, NCAR_DC_BERES, &
        IM*JM, LATS)

        ! print*,'Rank = ', trim(rank), ': beres_band%ngwv = ', beres_band%ngwv
        ! print*,'Rank = ', trim(rank), ': beres_band%dc = ', beres_band%dc
        ! print*,'Rank = ', trim(rank), ': beres_band%sum(abs(cref)  = ', sum(abs(beres_band%cref))
        ! print*,'Rank = ', trim(rank), ': beres_band%fcrit2 = ', beres_band%fcrit2
        ! print*,'Rank = ', trim(rank), ': beres_band%kwv = ', beres_band%kwv
        ! print*,'Rank = ', trim(rank), ': beres_band%effkwv = ', beres_band%effkwv

        ! print*,'Rank = ', trim(rank), ': beres_dc_dsrc%active = ', beres_dc_desc%active
        ! print*,'Rank = ', trim(rank), ': beres_dc_dsrc%storm_shift = ', beres_dc_desc%storm_shift
        ! print*,'Rank = ', trim(rank), ': beres_dc_dsrc%min_hdepth = ', beres_dc_desc%min_hdepth
        ! print*,'Rank = ', trim(rank), ': beres_dc_dsrc%spectrum_source = ', beres_dc_desc%spectrum_source
        ! print*,'Rank = ', trim(rank), ': beres_dc_dsrc%k = ', sum(abs(beres_dc_desc%k))
        ! print*,'Rank = ', trim(rank), ': beres_dc_dsrc%tndmax = ', beres_dc_desc%tndmax
        ! print*,'Rank = ', trim(rank), ': beres_dc_dsrc%maxh = ', beres_dc_desc%maxuh
        ! print*,'Rank = ', trim(rank), ': beres_dc_dsrc%hd = ', sum(abs(beres_dc_desc%hd))
        ! print*,'Rank = ', trim(rank), ': beres_dc_dsrc%mfcc = ', sum(abs(beres_dc_desc%mfcc))
        ! print*,'Rank = ', trim(rank), ': beres_dc_dsrc%taubck = ', sum(abs(beres_dc_desc%taubck))

    open(newunit=fileID, file=trim(dirName) // "/NCAR_SC_BERES_SRC_LEVEL_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_SC_BERES_SRC_LEVEL
    ! print*,'NCAR_BKG_TNDMAX = ', NCAR_BKG_TNDMAX
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_SC_BERES_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_SC_BERES
    ! print*,'NCAR_BKG_TNDMAX = ', NCAR_BKG_TNDMAX
    close(fileID)

    call gw_beres_init( trim(dirName) // "/newmfspectra40_dc25.nc",  &
            beres_band,  &
            beres_sc_desc,  &
            NCAR_BKG_PGWV, NCAR_BKG_GW_DC, NCAR_BKG_FCRIT2,  &
            NCAR_BKG_WAVELENGTH, NCAR_SC_BERES_SRC_LEVEL, &
            0.0, .FALSE., NCAR_ET_TAUBGND, NCAR_BKG_TNDMAX, NCAR_SC_BERES, &
            IM*JM, LATS)

    ! print*,'Rank = ', trim(rank), ': beres_band%ngwv = ', beres_band%ngwv
    ! print*,'Rank = ', trim(rank), ': beres_band%dc = ', beres_band%dc
    ! print*,'Rank = ', trim(rank), ': beres_band%sum(abs(cref)  = ', sum(abs(beres_band%cref))
    ! print*,'Rank = ', trim(rank), ': beres_band%fcrit2 = ', beres_band%fcrit2
    ! print*,'Rank = ', trim(rank), ': beres_band%kwv = ', beres_band%kwv
    ! print*,'Rank = ', trim(rank), ': beres_band%effkwv = ', beres_band%effkwv

    ! print*,'Rank = ', trim(rank), ': beres_sc_desc%active = ', beres_sc_desc%active
    ! print*,'Rank = ', trim(rank), ': beres_sc_desc%storm_shift = ', beres_sc_desc%storm_shift
    ! print*,'Rank = ', trim(rank), ': beres_sc_desc%min_hdepth = ', beres_sc_desc%min_hdepth
    ! print*,'Rank = ', trim(rank), ': beres_sc_desc%spectrum_source = ', beres_sc_desc%spectrum_source
    ! print*,'Rank = ', trim(rank), ': beres_sc_desc%k = ', sum(abs(beres_sc_desc%k))
    ! print*,'Rank = ', trim(rank), ': beres_sc_desc%tndmax = ', beres_sc_desc%tndmax
    ! print*,'Rank = ', trim(rank), ': beres_sc_desc%maxh = ', beres_sc_desc%maxuh
    ! print*,'Rank = ', trim(rank), ': beres_sc_desc%hd = ', sum(abs(beres_sc_desc%hd))
    ! print*,'Rank = ', trim(rank), ': beres_sc_desc%mfcc = ', sum(abs(beres_sc_desc%mfcc))
    ! print*,'Rank = ', trim(rank), ': beres_sc_desc%taubck = ', sum(abs(beres_sc_desc%taubck))
!     call gw_oro_init(oro_band)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_ORO_GW_DC_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_ORO_GW_DC
    ! print*,'NCAR_ORO_GW_DC = ', NCAR_ORO_GW_DC
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_ORO_FCRIT2_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_ORO_FCRIT2
    ! print*,'NCAR_ORO_FCRIT2 = ', NCAR_ORO_FCRIT2
    close(fileID)


    open(newunit=fileID, file=trim(dirName) // "/NCAR_ORO_FCRIT2_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_ORO_FCRIT2
    ! print*,'NCAR_ORO_FCRIT2 = ', NCAR_ORO_FCRIT2
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_ORO_WAVELENGTH_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_ORO_WAVELENGTH
    ! print*,'NCAR_ORO_WAVELENGTH = ', NCAR_ORO_WAVELENGTH
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_ORO_PGWV_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_ORO_PGWV
    ! print*,'NCAR_ORO_PGWV = ', NCAR_ORO_PGWV
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_ORO_SOUTH_FAC_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_ORO_SOUTH_FAC
    ! print*,'NCAR_ORO_SOUTH_FAC = ', NCAR_ORO_SOUTH_FAC
    close(fileID)


    call gw_oro_init ( oro_band, NCAR_ORO_GW_DC, &
                        NCAR_ORO_FCRIT2, NCAR_ORO_WAVELENGTH, NCAR_ORO_PGWV, &
                        NCAR_ORO_SOUTH_FAC )

    ! print*,'Rank = ', trim(rank), ': oro_band%ngwv = ', oro_band%ngwv
    ! print*,'Rank = ', trim(rank), ': oro_band%dc = ', oro_band%dc
    ! print*,'Rank = ', trim(rank), ': oro_band%sum(abs(cref)  = ', sum(abs(oro_band%cref))
    ! print*,'Rank = ', trim(rank), ': oro_band%fcrit2 = ', oro_band%fcrit2
    ! print*,'Rank = ', trim(rank), ': oro_band%kwv = ', oro_band%kwv
    ! print*,'Rank = ', trim(rank), ': oro_band%effkwv = ', oro_band%effkwv

    open(newunit=fileID, file=trim(dirName) // "/NCAR_NRDG_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_NRDG
    ! print*,'NCAR_NRDG = ', NCAR_NRDG
    close(fileID)

    ! Ridge Scheme
    if(NCAR_NRDG > 0) then
        print*, "Running Ridge Scheme"
        open(newunit=fileID, file=trim(dirName) // "/NCAR_ORO_TNDMAX_v2_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
        read(fileID) NCAR_ORO_TNDMAX
        ! print*,'NCAR_ORO_TNDMAX_v2 = ', NCAR_ORO_TNDMAX
        close(fileID)

        NCAR_ORO_TNDMAX = NCAR_ORO_TNDMAX/86400.0
        call gw_rdg_init ( rdg_band, NCAR_ORO_GW_DC, NCAR_ORO_FCRIT2, NCAR_ORO_WAVELENGTH, NCAR_ORO_TNDMAX, NCAR_ORO_PGWV )
    endif

    allocate(alpha(LM+1))
    allocate(alpha_ref(LM+1))
    allocate(PREF (LM+1))

    open(newunit=fileID, file=trim(dirName) // "/PREF_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) PREF
    ! print*,'sum(abs(PREF)) = ', sum(abs(PREF))
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/ALPHA_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) alpha
    ! print*,'sum(abs(alpha)) = ', sum(abs(alpha))
    close(fileID)

    call gw_newtonian_set(LM, PREF, alpha)

    open(newunit=fileID, file=trim(dirName) // "/ALPHA_" // trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) alpha_ref
    ! print*,'sum(abs(alpha)) = ', sum(abs(alpha))
    close(fileID)

    ! print*,'sum(alpha - alpha_ref) = ', sum(alpha - alpha_ref)
    ! print*,'sum(alpha_ref) = ', sum(alpha_ref)

!     ! Input arrays
    allocate(PLE   (IM, JM, LM+1))
    allocate(T     (IM, JM, LM))
    allocate(Q     (IM, JM, LM))
!     allocate(U     (IM, JM, LM))
!     allocate(V     (IM, JM, LM))
!     allocate(HT_dpc(IM, JM, LM))
!     allocate(SGH   (IM, JM))
!     allocate(PREF  (        LM+1))
    allocate(PMID  (IM, JM, LM))
    allocate(PMID_ref  (IM, JM, LM))
    allocate(PDEL  (IM, JM, LM))
    allocate(PDEL_ref  (IM, JM, LM))
    allocate(RPDEL (IM, JM, LM))
    allocate(RPDEL_ref (IM, JM, LM))
    allocate(PILN  (IM, JM, LM+1))
    allocate(PILN_ref  (IM, JM, LM+1))
    allocate(PMLN  (IM, JM, LM))
    allocate(PMLN_ref  (IM, JM, LM))
    allocate(ZI( IM, JM, LM+1))
    allocate(ZI_ref( IM, JM, LM+1))
    allocate(ZM    (IM, JM, LM))
    allocate(ZM_ref    (IM, JM, LM))
!     ! Output arrays
!     allocate(DUDT_GWD (IM, JM, LM))
!     allocate(DVDT_GWD (IM, JM, LM))
!     allocate(DTDT_GWD (IM, JM, LM))
!     allocate(DUDT_ORG (IM, JM, LM))
!     allocate(DVDT_ORG (IM, JM, LM))
!     allocate(DTDT_ORG (IM, JM, LM))
!     allocate(TAUXO_TMP(IM, JM))
!     allocate(TAUYO_TMP(IM, JM))
!     allocate(TAUXO_3D (IM, JM, LM+1))
!     allocate(TAUYO_3D (IM, JM, LM+1))
!     allocate(FEO_3D   (IM, JM, LM+1))
!     allocate(TAUXB_TMP(IM, JM))
!     allocate(TAUYB_TMP(IM, JM))
!     allocate(TAUXB_3D (IM, JM, LM+1))
!     allocate(TAUYB_3D (IM, JM, LM+1))
!     allocate(FEB_3D   (IM, JM, LM+1))
!     allocate(FEPO_3D  (IM, JM, LM+1))
!     allocate(FEPB_3D  (IM, JM, LM+1))
!     allocate(DUBKGSRC (IM, JM, LM))
!     allocate(DVBKGSRC (IM, JM, LM))
!     allocate(DTBKGSRC (IM, JM, LM))

!     ! Output Array References
!     allocate(DUDT_GWD_ref (IM, JM, LM))
!     allocate(DVDT_GWD_ref (IM, JM, LM))
!     allocate(DTDT_GWD_ref (IM, JM, LM))
!     allocate(DUDT_ORG_ref (IM, JM, LM))
!     allocate(DVDT_ORG_ref (IM, JM, LM))
!     allocate(DTDT_ORG_ref (IM, JM, LM))
!     allocate(TAUXO_TMP_ref(IM, JM))
!     allocate(TAUYO_TMP_ref(IM, JM))
!     allocate(TAUXO_3D_ref (IM, JM, LM+1))
!     allocate(TAUYO_3D_ref (IM, JM, LM+1))
!     allocate(FEO_3D_ref   (IM, JM, LM+1))
!     allocate(TAUXB_TMP_ref(IM, JM))
!     allocate(TAUYB_TMP_ref(IM, JM))
!     allocate(TAUXB_3D_ref (IM, JM, LM+1))
!     allocate(TAUYB_3D_ref (IM, JM, LM+1))
!     allocate(FEB_3D_ref   (IM, JM, LM+1))
!     allocate(FEPO_3D_ref  (IM, JM, LM+1))
!     allocate(FEPB_3D_ref  (IM, JM, LM+1))
!     allocate(DUBKGSRC_ref (IM, JM, LM))
!     allocate(DVBKGSRC_ref (IM, JM, LM))
!     allocate(DTBKGSRC_ref (IM, JM, LM))

    open(newunit=fileID, file=trim(dirName) // "/PLE_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) PLE
    close(fileID)

!     !write(*,*) 'sum(PLE) = ', sum(PLE)

    open(newunit=fileID, file=trim(dirName) // "/T_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) T
    close(fileID)

!     !write(*,*) 'sum(T) = ', sum(T)

    open(newunit=fileID, file=trim(dirName) // "/Q_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) Q
    close(fileID)

!     open(newunit=fileID, file=trim(dirName) // "/U.in", status='old', form="unformatted", action="read")
!     read(fileID) U
!     close(fileID)

!     !write(*,*) 'sum(U) = ', sum(U)

!     open(newunit=fileID, file=trim(dirName) // "/V.in", status='old', form="unformatted", action="read")
!     read(fileID) V
!     close(fileID)

!     !write(*,*) 'sum(V) = ', sum(V)

!     open(newunit=fileID, file=trim(dirName) // "/ht_dpc.in", status='old', form="unformatted", action="read")
!     read(fileID) HT_dpc
!     close(fileID)

!     !write(*,*) 'sum(HT_dpc) = ', sum(HT_dpc)

!     open(newunit=fileID, file=trim(dirName) // "/sgh.in", status='old', form="unformatted", action="read")
!     read(fileID) SGH
!     close(fileID)

!     !write(*,*) 'sum(SGH) = ', sum(SGH)

!     open(newunit=fileID, file=trim(dirName) // "/pref.in", status='old', form="unformatted", action="read")
!     read(fileID) PREF
!     close(fileID)

!     !write(*,*) 'sum(PREF) = ', sum(PREF)

    ! open(newunit=fileID, file=trim(dirName) // "/PMID_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    ! read(fileID) PMID
    ! close(fileID)

!     !write(*,*) 'sum(PMID) = ', sum(PMID)

    ! open(newunit=fileID, file=trim(dirName) // "/PDEL_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    ! read(fileID) PDEL
    ! close(fileID)

!     !write(*,*) 'sum(PDEL) = ', sum(PDEL)

    ! open(newunit=fileID, file=trim(dirName) // "/PRDEL_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    ! read(fileID) RPDEL
    ! close(fileID)

!     !write(*,*) 'sum(RPDEL) = ', sum(RPDEL)

    ! open(newunit=fileID, file=trim(dirName) // "/PILN_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    ! read(fileID) PILN
    ! close(fileID)

!     !write(*,*) 'sum(PILN) = ', sum(PILN)

    ! open(newunit=fileID, file=trim(dirName) // "/PMLN_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    ! read(fileID) PMLN
    ! close(fileID)

!     open(newunit=fileID, file=trim(dirName) // "/zm.in", status='old', form="unformatted", action="read")
!     read(fileID) ZM
!     close(fileID)

!     !write(*,*) 'sum(ZM) = ', sum(ZM)

    ! open(newunit=fileID, file=trim(dirName) // "/LATS_v2_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    ! read(fileID) LATS
    ! close(fileID)

!     !write(*,*) 'sum(LATS) = ', sum(LATS)

    CALL PREGEO(IM*JM,   LM,   &
        PLE, LATS,   PMID,  PDEL, RPDEL,     PILN,     PMLN)

    ! Compute ZM
    !-------------

    call GEOPOTENTIAL( IM*JM, LM,                  &
        PILN,  PMLN,  PLE,   PMID, PDEL, RPDEL,   &
        T,     Q,     ZI,    ZM                   )

    ! open(newunit=fileID, file=trim(dirName) // "/ZI_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    ! read(fileID) ZI_ref
    ! close(fileID)

    ! open(newunit=fileID, file=trim(dirName) // "/ZM_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    ! read(fileID) ZM_ref
    ! close(fileID)

    ! print*,'sum(ZI_ref-ZI) = ', sum(ZI_ref-ZI)
    ! print*,'sum(ZM_ref-ZM) = ', sum(ZM_ref-ZM)

    allocate(GBXAR_TMP(IM, JM))
    allocate(GBXAR_TMP_ref(IM, JM))
    allocate(GBXAR    (IM, JM))
    allocate(ANGLL    (IM, JM, 16))
    allocate(ANGLL_ref(IM, JM, 16))
    allocate(HWDTH    (IM, JM, 16))
    allocate(CLNGT    (IM, JM, 16))
    allocate(KWVRDG(IM, JM, 16))
    allocate(KWVRDG_ref(IM, JM, 16))
    allocate(EFFRDG(IM, JM, 16))
    allocate(EFFRDG_ref(IM, JM, 16))

    open(newunit=fileID, file=trim(dirName) // "/GBXAR_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) GBXAR
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/ANGLL_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) ANGLL
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/HWDTH_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) HWDTH
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/CLNGT_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) CLNGT
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/KWVRDG_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) KWVRDG
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/EFFRDG_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) EFFRDG
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_EFFGWORO_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_EFFGWORO
    close(fileID)


    GBXAR_TMP = GBXAR * (MAPL_RADIUS/1000.)**2 ! transform to km^2
    WHERE (ANGLL < -180)
        ANGLL = 0.0
    END WHERE

    do nrdg = 1, NCAR_NRDG
        KWVRDG(:,:,nrdg) = 0.001/(HWDTH(:,:,nrdg)+0.001)
        EFFRDG(:,:,nrdg) = NCAR_EFFGWORO*(HWDTH(:,:,nrdg)*CLNGT(:,:,nrdg))/GBXAR_TMP
    enddo

    ! open(newunit=fileID, file=trim(dirName) // "/ANGLL_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    ! read(fileID) ANGLL_ref
    ! close(fileID)

    ! open(newunit=fileID, file=trim(dirName) // "/KWVRDG_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    ! read(fileID) KWVRDG_ref
    ! close(fileID)

    ! open(newunit=fileID, file=trim(dirName) // "/EFFRDG_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    ! read(fileID) EFFRDG_ref
    ! close(fileID)

    ! open(newunit=fileID, file=trim(dirName) // "/GBXAR_TMP_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    ! read(fileID) GBXAR_TMP_ref
    ! close(fileID)

    ! print*,'sum(ANGLL_ref-ANGLL) = ', sum(ANGLL_ref-ANGLL)
    ! print*,'sum(KWVRDG_ref-KWVRDG) = ', sum(KWVRDG_ref-KWVRDG)
    ! print*,'sum(EFFRDG_ref-EFFRDG) = ', sum(EFFRDG_ref-EFFRDG)
    ! print*,'sum(GBXAR_TMP_ref-GBXAR_TMP) = ', sum(GBXAR_TMP_ref-GBXAR_TMP)
    ! print*,'NCAR_NRDG = ', NCAR_NRDG

! !$acc  data copyin(oro_band, beres_band, beres_desc) &
! !!$acc       copyin(beres_band%cref, beres_band%ngwv, beres_band%dc) &
! !!$acc       copyin(oro_band%cref, oro_band%ngwv) &
! !!$acc       copyin(beres_desc%hd, beres_desc%mfcc, beres_desc%k) &
! !$acc       copyin(PLE, T, U, V, HT_dpc) &
! !$acc       copyin(SGH, PREF, PMID, PDEL, RPDEL, PILN, ZM, LATS) &
! !$acc       copyout(DUDT_GWD, DVDT_GWD, DTDT_GWD) &
! !$acc       copyout(DUDT_ORG, DVDT_ORG, DTDT_ORG) &
! !$acc       copyout(TAUXO_TMP, TAUYO_TMP, TAUXO_3D, TAUYO_3D, FEO_3D) &
! !$acc       copyout(TAUXB_TMP, TAUYB_TMP, TAUXB_3D, TAUYB_3D, FEB_3D) &
! !$acc       copyout(FEPO_3D, FEPB_3D) &
! !$acc       create(DUBKGSRC, DVBKGSRC, DTBKGSRC)
! !$acc  data copyin(beres_band%cref, beres_band%ngwv, beres_band%dc) &
! !$acc       copyin(oro_band%cref, oro_band%ngwv) &
! !$acc       copyin(beres_desc%hd, beres_desc%mfcc, beres_desc%k)

!     call cpu_time(t_start)

!     call gw_intr_ncar(IM*JM,    LM,         DT,                  &
!             PGWV,      beres_desc, beres_band, oro_band,            &
!             PLE,       T,          U,          V,      HT_dpc,      &
!             SGH,       PREF,                                        &
!             PMID,      PDEL,       RPDEL,      PILN,   ZM,    LATS, &
!             DUDT_GWD,  DVDT_GWD,   DTDT_GWD,                        &
!             DUDT_ORG,  DVDT_ORG,   DTDT_ORG,                        &
!             TAUXO_TMP, TAUYO_TMP,  TAUXO_3D,   TAUYO_3D,  FEO_3D,   &
!             TAUXB_TMP, TAUYB_TMP,  TAUXB_3D,   TAUYB_3D,  FEB_3D,   &
!             FEPO_3D,   FEPB_3D,    DUBKGSRC,   DVBKGSRC,  DTBKGSRC, &
!             BGSTRESSMAX, effgworo, effgwbkg,   RC            )

!     call cpu_time(t_end)
! !$acc end data
! !$acc end data

!     open(newunit=fileID, file=trim(dirName) // "/dudt_gwd.out", status='old', form="unformatted", action="read")
!     read(fileID) DUDT_GWD_ref
!     close(fileID)

!     !write(*,*) 'sum(DUDT_GWD_ref) = ', sum(DUDT_GWD_ref)

!     open(newunit=fileID, file=trim(dirName) // "/dvdt_gwd.out", status='old', form="unformatted", action="read")
!     read(fileID) DVDT_GWD_ref
!     close(fileID)

!     !write(*,*) 'sum(DVDT_GWD_ref) = ', sum(DVDT_GWD_ref)

!     open(newunit=fileID, file=trim(dirName) // "/dtdt_gwd.out", status='old', form="unformatted", action="read")
!     read(fileID) DTDT_GWD_ref
!     close(fileID)

!     !write(*,*) 'sum(DTDT_GWD_ref) = ', sum(DTDT_GWD_ref)

!     open(newunit=fileID, file=trim(dirName) // "/dudt_org.out", status='old', form="unformatted", action="read")
!     read(fileID) DUDT_ORG_ref
!     close(fileID)

!     !write(*,*) 'sum(DUDT_ORG_ref) = ', sum(DUDT_ORG_ref)

!     open(newunit=fileID, file=trim(dirName) // "/dvdt_org.out", status='old', form="unformatted", action="read")
!     read(fileID) DVDT_ORG_ref
!     close(fileID)

!     !write(*,*) 'sum(DVDT_ORG_ref) = ', sum(DVDT_ORG_ref)

!     open(newunit=fileID, file=trim(dirName) // "/dtdt_org.out", status='old', form="unformatted", action="read")
!     read(fileID) DTDT_ORG_ref
!     close(fileID)

!     !write(*,*) 'sum(DTDT_ORG_ref) = ', sum(DTDT_ORG_ref)

!     open(newunit=fileID, file=trim(dirName) // "/tauxo_tmp.out", status='old', form="unformatted", action="read")
!     read(fileID) TAUXO_TMP_ref
!     close(fileID)

!     !write(*,*) 'sum(TAUXO_TMP_ref) = ', sum(TAUXO_TMP_ref)

!     open(newunit=fileID, file=trim(dirName) // "/tauyo_tmp.out", status='old', form="unformatted", action="read")
!     read(fileID) TAUYO_TMP_ref
!     close(fileID)

!     !write(*,*) 'sum(TAUYO_TMP_ref) = ', sum(TAUYO_TMP_ref)

!     open(newunit=fileID, file=trim(dirName) // "/tauxo_3d.out", status='old', form="unformatted", action="read")
!     read(fileID) TAUXO_3D_ref
!     close(fileID)

!     !write(*,*) 'sum(TAUXO_3D_ref) = ', sum(TAUXO_3D_ref)

!     open(newunit=fileID, file=trim(dirName) // "/tauyo_3d.out", status='old', form="unformatted", action="read")
!     read(fileID) TAUYO_3D_ref
!     close(fileID)

!     !write(*,*) 'sum(TAUYO_3D_ref) = ', sum(TAUYO_3D_ref)

!     open(newunit=fileID, file=trim(dirName) // "/feo_3d.out", status='old', form="unformatted", action="read")
!     read(fileID) FEO_3D_ref
!     close(fileID)

!     !write(*,*) 'sum(FEO_3D_ref) = ', sum(FEO_3D_ref)

!     open(newunit=fileID, file=trim(dirName) // "/tauxb_tmp.out", status='old', form="unformatted", action="read")
!     read(fileID) TAUXB_TMP_ref
!     close(fileID)

!     !write(*,*) 'sum(TAUXB_TMP_ref) = ', sum(TAUXB_TMP_ref)

!     open(newunit=fileID, file=trim(dirName) // "/tauyb_tmp.out", status='old', form="unformatted", action="read")
!     read(fileID) TAUYB_TMP_ref
!     close(fileID)

!     !write(*,*) 'sum(TAUYB_TMP_ref) = ', sum(TAUYB_TMP_ref)

!     open(newunit=fileID, file=trim(dirName) // "/tauxb_3d.out", status='old', form="unformatted", action="read")
!     read(fileID) TAUXB_3D_ref
!     close(fileID)

!     !write(*,*) 'sum(TAUXB_3D_ref) = ', sum(TAUXB_3D_ref)

!     open(newunit=fileID, file=trim(dirName) // "/tauyb_3d.out", status='old', form="unformatted", action="read")
!     read(fileID) TAUYB_3D_ref
!     close(fileID)

!     !write(*,*) 'sum(TAUYB_3D_ref) = ', sum(TAUYB_3D_ref)

!     open(newunit=fileID, file=trim(dirName) // "/feb_3d.out", status='old', form="unformatted", action="read")
!     read(fileID) FEB_3D_ref
!     close(fileID)

!     !write(*,*) 'sum(FEB_3D_ref) = ', sum(FEB_3D_ref)

!     open(newunit=fileID, file=trim(dirName) // "/fepo_3d.out", status='old', form="unformatted", action="read")
!     read(fileID) FEPO_3D_ref
!     close(fileID)

!     !write(*,*) 'sum(FEPO_3D_ref) = ', sum(FEPO_3D_ref)

!     open(newunit=fileID, file=trim(dirName) // "/fepb_3d.out", status='old', form="unformatted", action="read")
!     read(fileID) FEPB_3D_ref
!     close(fileID)

!     !write(*,*) 'sum(FEPB_3D_ref) = ', sum(FEPB_3D_ref)

!     open(newunit=fileID, file=trim(dirName) // "/dubkgsrc.out", status='old', form="unformatted", action="read")
!     read(fileID) DUBKGSRC_ref
!     close(fileID)

!     !write(*,*) 'sum(DUBKGSRC_ref) = ', sum(DUBKGSRC_ref)

!     open(newunit=fileID, file=trim(dirName) // "/dvbkgsrc.out", status='old', form="unformatted", action="read")
!     read(fileID) DVBKGSRC_ref
!     close(fileID)

!     !write(*,*) 'sum(DVBKGSRC_ref) = ', sum(DVBKGSRC_ref)

!     open(newunit=fileID, file=trim(dirName) // "/dtbkgsrc.out", status='old', form="unformatted", action="read")
!     read(fileID) DTBKGSRC_ref
!     close(fileID)

!     !write(*,*) 'sum(DTBKGSRC_ref) = ', sum(DTBKGSRC_ref)

!     write(*,*) 'Sum Abs Diff DUDT_GWD = ', sum(abs(DUDT_GWD_ref)) - sum(abs(DUDT_GWD)), sum(abs(DUDT_GWD_ref)), sum(abs(DUDT_GWD))
!     write(*,*) 'Sum Abs Diff DVDT_GWD = ', sum(abs(DVDT_GWD_ref)) - sum(abs(DVDT_GWD)), sum(abs(DVDT_GWD_ref)), sum(abs(DVDT_GWD))
!     write(*,*) 'Sum Abs Diff DTDT_GWD = ', sum(abs(DTDT_GWD_ref)) - sum(abs(DTDT_GWD)), sum(abs(DTDT_GWD_ref)), sum(abs(DTDT_GWD))
!     write(*,*) 'Sum Abs Diff DUDT_ORG = ', sum(abs(DUDT_ORG_ref)) - sum(abs(DUDT_ORG)), sum(abs(DUDT_ORG_ref)), sum(abs(DUDT_ORG))
!     write(*,*) 'Sum Abs Diff DVDT_ORG = ', sum(abs(DVDT_ORG_ref)) - sum(abs(DVDT_ORG)), sum(abs(DVDT_ORG_ref)), sum(abs(DVDT_ORG))
!     write(*,*) 'Sum Abs Diff DTDT_ORG = ', sum(abs(DTDT_ORG_ref)) - sum(abs(DTDT_ORG)), sum(abs(DTDT_ORG_ref)), sum(abs(DTDT_ORG))

!     write(*,*) 'Execution time of gw_intr_ncar = ', t_end - t_start
    
    contains


  subroutine pregeo(pcols,pver,&
    ple,lats,pmid,pdel,rpdel,piln,pmln)

    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!

    integer, intent(in) :: pcols    ! Number of longitudes
    integer, intent(in) :: pver     ! Number of vertical layers

    real,    intent(in) :: ple (pcols,pver+1)    ! Interface pressures
    real,    intent(in) :: lats(pcols)           ! latitude in radian

! Output arguments

    real,    intent(out) :: pmid  (pcols,pver)   ! Midpoint pressures
    real,    intent(out) :: pdel  (pcols,pver)   ! layer thickness
    real,    intent(out) :: rpdel (pcols,pver)   ! inverse of layer thickness
    real,    intent(out) :: piln  (pcols,pver+1) ! Log interface pressures
    real,    intent(out) :: pmln  (pcols,pver)   ! Log midpoint pressures

!
!---------------------------Local variables-----------------------------
!
    integer :: i,k

    real    :: hvsd  ! Efficiency factor

    real, parameter :: PI_GWD  = 4.0*atan(1.0)  ! This is *not* MAPL_PI

!
!-----------------------------------------------------------------------
!

! Form pressure factors
!----------------------

    I_LOOP: DO I = 1, PCOLS

       DO K = 1, PVER
           PMID(I,K) = 0.5*(  PLE(I,K  ) + PLE(I,K+1) )
           PDEL(I,K) =        PLE(I,K+1) - PLE(I,K  )
          RPDEL(I,K) = 1.0 / PDEL(I,K)
           PILN(I,K) = log(   PLE(I,K) )
           PMLN(I,K) = log(  PMID(I,K) ) !
       END DO
       PILN(I,PVER+1)  = log( PLE(I,PVER+1)  )
    END DO I_LOOP

  end subroutine pregeo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine geopotential(pcols  , pver   ,                   &
    piln   , pmln   , pint  , pmid   , pdel   , rpdel  , &
    t      , q      , zi     , zm     )

!-----------------------------------------------------------------------
!
! Purpose:
! Compute the geopotential height (above the surface) at the midpoints and
! interfaces using the input temperatures and pressures.
! Author: B.Boville, Feb 2001 from earlier code by Boville and S.J. Lin
!
!-----------------------------------------------------------------------

implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
integer, intent(in) :: pcols                ! Number of longitudes
integer, intent(in) :: pver                 ! Number of vertical layers

real,    intent(in) :: piln (pcols,pver+1)  ! Log interface pressures
real,    intent(in) :: pmln (pcols,pver)    ! Log midpoint pressures
real,    intent(in) :: pint (pcols,pver+1)  ! Interface pressures
real,    intent(in) :: pmid (pcols,pver)    ! Midpoint pressures
real,    intent(in) :: pdel (pcols,pver)    ! layer thickness
real,    intent(in) :: rpdel(pcols,pver)    ! inverse of layer thickness
real,    intent(in) :: t    (pcols,pver)    ! temperature
real,    intent(in) :: q    (pcols,pver)    ! specific humidity

! Output arguments

real,    intent(out) :: zi(pcols,pver+1)    ! Height above surface at interfaces
real,    intent(out) :: zm(pcols,pver)      ! Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
!
logical  :: fvdyn              ! finite volume dynamics
integer  :: i,k                ! Lon, level indices
real     :: hkk                ! diagonal element of hydrostatic matrix
real     :: hkl                ! off-diagonal element
real     :: tv                 ! virtual temperature
real     :: tvfac              ! Tv/T

real, parameter :: ROG     = MAPL_RGAS/MAPL_GRAV
!
!-----------------------------------------------------------------------
!

! Set dynamics flag

fvdyn = .true.

! The surface height is zero by definition.

I_LOOP: do i = 1, pcols

  zi(i,pver+1) = 0.0

! Compute zi, zm from bottom up.
! Note, zi(i,k) is the interface above zm(i,k)

  do k = pver, 1, -1

! First set hydrostatic elements consistent with dynamics

     if (fvdyn) then
        hkl = piln(i,k+1) - piln(i,k)
        hkk = piln(i,k+1) - pmln(i,k)
     else
        hkl = pdel(i,k) / pmid(i,k)
        hkk = 0.5 * hkl
     end if

! Now compute tv, zm, zi

     tvfac   = 1. + MAPL_VIREPS * q(i,k)
     tv      = t(i,k) * tvfac

     zm(i,k) = zi(i,k+1) + ROG * tv * hkk
     zi(i,k) = zi(i,k+1) + ROG * tv * hkl
  end do
end do I_LOOP

return
end subroutine geopotential
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
