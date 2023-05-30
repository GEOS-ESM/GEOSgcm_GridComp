program gwd_standalone

    use gw_drag_ncar, only : gw_intr_ncar
    use gw_convect,   only : BeresSourceDesc, gw_beres_init
    use gw_common,    only : GWBand, gw_common_init, gw_newtonian_set
    use gw_oro,       only : gw_oro_init
    use gw_rdg,       only : gw_rdg_init
    use MAPL_ConstantsMod

    implicit none

    integer               :: I, J, L, IM, JM, LM, fileID, status, PGWV, RC, nrdg
    integer               :: NCAR_BKG_PGWV, NCAR_ORO_PGWV, NCAR_NRDG, ikpbl
    real, parameter       :: dxmax_ss = 12000.0, dxmin_ss =  3000.0
    real                  :: DT, bgstressmax, effgwbkg, effgworo,  H0, HH, Z1, TAU1
    real                  :: t_start, t_end
    real                  :: NCAR_PRNDL, NCAR_QBO_HDEPTH_SCALING
    real                  :: NCAR_HR_CF, NCAR_BKG_GW_DC, NCAR_BKG_FCRIT2
    real                  :: NCAR_BKG_WAVELENGTH, NCAR_DC_BERES_SRC_LEVEL
    real                  :: NCAR_SC_BERES_SRC_LEVEL
    real                  :: NCAR_ET_TAUBGND, NCAR_BKG_TNDMAX
    real                  :: NCAR_ORO_GW_DC, NCAR_ORO_FCRIT2, NCAR_ORO_WAVELENGTH
    real                  :: NCAR_ORO_SOUTH_FAC, NCAR_ORO_TNDMAX, NCAR_EFFGWORO
    real                  :: NCAR_EFFGWBKG, effbeljaars, limbeljaars, tcrib, var_temp, wsp
    character*100         :: Errstring, dirName, rank
    logical               :: NCAR_TAU_TOP_ZERO, NCAR_DC_BERES, NCAR_SC_BERES

    type(BeresSourceDesc) :: beres_dc_desc, beres_sc_desc
    type(GWBand)          :: beres_band, oro_band, rdg_band

    real, dimension(:,:,:), allocatable :: PLE, T, Q, U, V, HT_dpc, ZI
    real, dimension(:,:,:), allocatable :: PMID, PDEL, RPDEL, PILN, ZM, PMLN
    real, dimension(:,:,:), allocatable :: DUDT_GWD_NCAR, DVDT_GWD_NCAR, DTDT_GWD_NCAR, DUDT_ORG_NCAR, DVDT_ORG_NCAR, DTDT_ORG_NCAR
    real, dimension(:,:,:), allocatable :: DUDT_GWD_GEOS, DVDT_GWD_GEOS, DTDT_GWD_GEOS, DUDT_ORG_GEOS, DVDT_ORG_GEOS, DTDT_ORG_GEOS
    real, dimension(:,:,:), allocatable :: DUDT_GWD, DVDT_GWD, DTDT_GWD, DUDT_ORG, DVDT_ORG, DTDT_ORG
    real, dimension(:,:,:), allocatable :: TAUXO_3D, TAUYO_3D, TAUXB_3D, TAUYB_3D, FEO_3D, FEB_3D, FEPO_3D, FEPB_3D
    real, dimension(:,:,:), allocatable :: DUBKGSRC, DVBKGSRC, DTBKGSRC
    real, dimension(:,:,:), allocatable :: ANGLL, HWDTH, KWVRDG, CLNGT, EFFRDG
    real, dimension(:,:,:), allocatable :: HT_dc, HT_sc, QLDT_mst, QIDT_mst, MXDIS, ANIXY, THV, DUDT_TOFD, DVDT_TOFD
    real, dimension(:,:,:), allocatable :: DUDT_TOT, DVDT_TOT, DTDT_TOT, DUDT_RAH, DVDT_RAH, DTDT_RAH
    real, dimension(:,:),   allocatable :: LATS, SGH, TAUXO_TMP_NCAR, TAUYO_TMP_NCAR, TAUXB_TMP_NCAR, TAUYB_TMP_NCAR
    real, dimension(:,:),   allocatable :: TAUXO_TMP_GEOS, TAUYO_TMP_GEOS, TAUXB_TMP_GEOS, TAUYB_TMP_GEOS
    real, dimension(:,:),   allocatable :: TAUXO_TMP, TAUYO_TMP, TAUXB_TMP, TAUYB_TMP
    real, dimension(:,:),   allocatable :: GBXAR_TMP, GBXAR, PHIS, a2, VARFLT, AREA, HeFold
    real, dimension(:,:),   allocatable :: PEGWD_X, PEORO_X, PERAY_X, PEBKG_X, BKGERR_X
    real, dimension(:,:),   allocatable :: KEGWD_X, KEORO_X,  KERAY_X,  KEBKG_X, KERES_X
    real, dimension(:),     allocatable :: PREF, alpha

    real, dimension(:,:,:), allocatable :: DUDT_GWD_NCAR_ref, DVDT_GWD_NCAR_ref, DTDT_GWD_NCAR_ref, DUDT_ORG_NCAR_ref, DVDT_ORG_NCAR_ref, DTDT_ORG_NCAR_ref
    real, dimension(:,:,:), allocatable :: TAUXO_3D_ref, TAUYO_3D_ref, TAUXB_3D_ref, TAUYB_3D_ref
    real, dimension(:,:,:), allocatable :: FEO_3D_ref, FEB_3D_ref, FEPO_3D_ref, FEPB_3D_ref
    real, dimension(:,:,:), allocatable :: DUBKGSRC_ref, DVBKGSRC_ref, DTBKGSRC_ref, DUDT_TOFD_ref, DVDT_TOFD_ref
    real, dimension(:,:),   allocatable :: TAUXO_TMP_NCAR_ref, TAUYO_TMP_NCAR_ref, TAUXB_TMP_NCAR_ref, TAUYB_TMP_NCAR_ref
    real, dimension(:),     allocatable :: alpha_ref


    real, dimension(:,:,:), allocatable :: PDEL_ref, PILN_ref, RPDEL_ref, PMID_ref, PMLN_ref, ZI_ref, ZM_ref
    real, dimension(:,:,:), allocatable :: ANGLL_ref, KWVRDG_ref, EFFRDG_ref
    real, dimension(:,:),   allocatable :: GBXAR_TMP_ref
    real, dimension(:,:),   allocatable :: PEGWD_X_ref, PEORO_X_ref, PERAY_X_ref, PEBKG_X_ref, BKGERR_X_ref
    real, dimension(:,:),   allocatable :: KEGWD_X_ref, KEORO_X_ref,  KERAY_X_ref,  KEBKG_X_ref, KERES_X_ref
    real, dimension(:,:,:), allocatable :: DUDT_TOT_ref, DVDT_TOT_ref, DTDT_TOT_ref, DUDT_RAH_ref, DVDT_RAH_ref, DTDT_RAH_ref

    dirName = './new_c180_data/ncar_gwd'
    
    if(command_argument_count().ne.1) then
        print*, 'Missing arguments : <executable> <rank>'
        call exit(1)
    else
        call get_command_argument(1, rank)
    endif
    

    ! Note : The setup for IM, JM, LM, and DT is specific for C180.  Other resolutions may not verify
    IM = 180
    JM = 180
    LM = 72
    DT = 600.0

    PGWV = 4

    bgstressmax = 0.900
    ! Note : These values below can be changed
    effgwbkg    = 0.125
    effgworo    = 0.250
    effbeljaars = 0.1
    limbeljaars = 400.0/86400.0

    NCAR_TAU_TOP_ZERO = .TRUE.

    allocate(LATS  (IM, JM))
    allocate(PLE   (IM, JM, LM+1))
    allocate(T     (IM, JM, LM))
    allocate(Q     (IM, JM, LM))
    allocate(U     (IM, JM, LM))
    allocate(V     (IM, JM, LM))
    allocate(HT_dc (IM, JM, LM))
    allocate(HT_sc (IM, JM, LM))
    allocate(QLDT_mst (IM, JM, LM))
    allocate(QIDT_mst (IM, JM, LM))
    allocate(MXDIS(IM, JM, 16))
    allocate(ANIXY(IM, JM, 16))
    allocate(PHIS(IM, JM))
    allocate(SGH   (IM, JM))
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

    allocate(DUDT_GWD_NCAR (IM, JM, LM))
    allocate(DVDT_GWD_NCAR (IM, JM, LM))
    allocate(DTDT_GWD_NCAR (IM, JM, LM))
    allocate(DUDT_ORG_NCAR (IM, JM, LM))
    allocate(DVDT_ORG_NCAR (IM, JM, LM))
    allocate(DTDT_ORG_NCAR (IM, JM, LM))
    allocate(TAUXO_TMP_NCAR(IM, JM))
    allocate(TAUYO_TMP_NCAR(IM, JM))
    allocate(TAUXB_TMP_NCAR(IM, JM))
    allocate(TAUYB_TMP_NCAR(IM, JM))

    allocate(DUDT_GWD_GEOS (IM, JM, LM))
    allocate(DVDT_GWD_GEOS (IM, JM, LM))
    allocate(DTDT_GWD_GEOS (IM, JM, LM))
    allocate(DUDT_ORG_GEOS (IM, JM, LM))
    allocate(DVDT_ORG_GEOS (IM, JM, LM))
    allocate(DTDT_ORG_GEOS (IM, JM, LM))
    allocate(TAUXO_TMP_GEOS(IM, JM))
    allocate(TAUYO_TMP_GEOS(IM, JM))
    allocate(TAUXB_TMP_GEOS(IM, JM))
    allocate(TAUYB_TMP_GEOS(IM, JM))

    allocate(DUDT_GWD (IM, JM, LM))
    allocate(DVDT_GWD (IM, JM, LM))
    allocate(DTDT_GWD (IM, JM, LM))
    allocate(DUDT_ORG (IM, JM, LM))
    allocate(DVDT_ORG (IM, JM, LM))
    allocate(DTDT_ORG (IM, JM, LM))
    allocate(TAUXO_TMP(IM, JM))
    allocate(TAUYO_TMP(IM, JM))
    allocate(TAUXB_TMP(IM, JM))
    allocate(TAUYB_TMP(IM, JM))

    allocate(DUDT_GWD_NCAR_ref (IM, JM, LM))
    allocate(DVDT_GWD_NCAR_ref (IM, JM, LM))
    allocate(DTDT_GWD_NCAR_ref (IM, JM, LM))
    allocate(DUDT_ORG_NCAR_ref (IM, JM, LM))
    allocate(DVDT_ORG_NCAR_ref (IM, JM, LM))
    allocate(DTDT_ORG_NCAR_ref (IM, JM, LM))
    allocate(TAUXO_TMP_NCAR_ref(IM, JM))
    allocate(TAUYO_TMP_NCAR_ref(IM, JM))
    allocate(TAUXB_TMP_NCAR_ref(IM, JM))
    allocate(TAUYB_TMP_NCAR_ref(IM, JM))

    allocate(THV(IM, JM, LM))
    allocate(a2(IM, JM))
    allocate(VARFLT(IM, JM))
    allocate(HeFold(IM, JM))
    allocate(DUDT_TOFD (IM, JM, LM))
    allocate(DVDT_TOFD (IM, JM, LM))
    allocate(DUDT_TOFD_ref (IM, JM, LM))
    allocate(DVDT_TOFD_ref (IM, JM, LM))
    allocate(AREA(IM, JM))
    allocate(DUDT_TOT (IM, JM, LM))
    allocate(DVDT_TOT (IM, JM, LM))
    allocate(DTDT_TOT (IM, JM, LM))
    allocate(DUDT_RAH (IM, JM, LM))
    allocate(DVDT_RAH (IM, JM, LM))
    allocate(DTDT_RAH (IM, JM, LM))
    allocate(PEGWD_X(IM, JM))
    allocate(PEORO_X(IM, JM))
    allocate(PERAY_X(IM, JM))
    allocate(PEBKG_X(IM, JM))
    allocate(BKGERR_X(IM, JM))
    allocate(KEGWD_X(IM, JM))
    allocate(KEORO_X(IM, JM))
    allocate(KERAY_X(IM, JM))
    allocate(KEBKG_X(IM, JM))
    allocate(KERES_X(IM, JM))
    allocate(alpha(LM+1))
    allocate(PREF (LM+1))

    allocate(DUDT_TOT_REF (IM, JM, LM))
    allocate(DVDT_TOT_REF (IM, JM, LM))
    allocate(DTDT_TOT_REF (IM, JM, LM))
    allocate(DUDT_RAH_REF (IM, JM, LM))
    allocate(DVDT_RAH_REF (IM, JM, LM))
    allocate(DTDT_RAH_REF (IM, JM, LM))
    allocate(PEGWD_X_REF(IM, JM))
    allocate(PEORO_X_REF(IM, JM))
    allocate(PERAY_X_REF(IM, JM))
    allocate(PEBKG_X_REF(IM, JM))
    allocate(BKGERR_X_REF(IM, JM))
    allocate(KEGWD_X_REF(IM, JM))
    allocate(KEORO_X_REF(IM, JM))
    allocate(KERAY_X_REF(IM, JM))
    allocate(KEBKG_X_REF(IM, JM))
    allocate(KERES_X_REF(IM, JM))
    allocate(alpha_ref(LM+1))
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

    open(newunit=fileID, file=trim(dirName) // "/NCAR_SC_BERES_SRC_LEVEL_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_SC_BERES_SRC_LEVEL
    ! print*,'NCAR_BKG_TNDMAX = ', NCAR_BKG_TNDMAX
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_SC_BERES_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_SC_BERES
    ! print*,'NCAR_BKG_TNDMAX = ', NCAR_BKG_TNDMAX
    close(fileID)

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

    open(newunit=fileID, file=trim(dirName) // "/NCAR_NRDG_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_NRDG
    ! print*,'NCAR_NRDG = ', NCAR_NRDG
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/PREF_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) PREF
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/ALPHA_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) alpha
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/PLE_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) PLE
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/T_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) T
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/Q_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) Q
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/U_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) U
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/V_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) V
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/HT_dc_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) HT_dc
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/HT_sc_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) HT_sc
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/QLDT_mst_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QLDT_mst
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/QIDT_mst_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) QIDT_mst
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/MXDIS_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) MXDIS
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/ANIXY_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) ANIXY
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/PHIS_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) PHIS
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_EFFGWBKG_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_EFFGWBKG
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/SGH_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) SGH
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/NCAR_ORO_TNDMAX_v2_" // trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) NCAR_ORO_TNDMAX
    close(fileID)

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

    open(newunit=fileID, file=trim(dirName) // "/AREA_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) AREA
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/VARFLT_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) VARFLT
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/H0_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) H0
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/HH_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) HH
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/Z1_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) Z1
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/TAU1_"// trim(rank) // ".in", status='old', form="unformatted", action="read")
    read(fileID) TAU1
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

    call gw_beres_init( trim(dirName) // "/newmfspectra40_dc25.nc",  &
            beres_band,  &
            beres_sc_desc,  &
            NCAR_BKG_PGWV, NCAR_BKG_GW_DC, NCAR_BKG_FCRIT2,  &
            NCAR_BKG_WAVELENGTH, NCAR_SC_BERES_SRC_LEVEL, &
            0.0, .FALSE., NCAR_ET_TAUBGND, NCAR_BKG_TNDMAX, NCAR_SC_BERES, &
            IM*JM, LATS)

    call gw_oro_init ( oro_band, NCAR_ORO_GW_DC, &
                        NCAR_ORO_FCRIT2, NCAR_ORO_WAVELENGTH, NCAR_ORO_PGWV, &
                        NCAR_ORO_SOUTH_FAC )

    ! Ridge Scheme
    if(NCAR_NRDG > 0) then
        print*, "Running Ridge Scheme"

        NCAR_ORO_TNDMAX = NCAR_ORO_TNDMAX/86400.0
        call gw_rdg_init ( rdg_band, NCAR_ORO_GW_DC, NCAR_ORO_FCRIT2, NCAR_ORO_WAVELENGTH, NCAR_ORO_TNDMAX, NCAR_ORO_PGWV )
    endif

    call gw_newtonian_set(LM, PREF, alpha)

    open(newunit=fileID, file=trim(dirName) // "/ALPHA_" // trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) alpha_ref
    close(fileID)

    print*,'sum(alpha - alpha_ref) = ', sum(alpha - alpha_ref)

    CALL PREGEO(IM*JM,   LM,   &
        PLE, LATS,   PMID,  PDEL, RPDEL,     PILN,     PMLN)

    ! Compute ZM
    !-------------

    call GEOPOTENTIAL( IM*JM, LM,                  &
        PILN,  PMLN,  PLE,   PMID, PDEL, RPDEL,   &
        T,     Q,     ZI,    ZM                   )

    open(newunit=fileID, file=trim(dirName) // "/ZI_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) ZI_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/ZM_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) ZM_ref
    close(fileID)

    print*,'sum(ZI_ref-ZI) = ', sum(ZI_ref-ZI)
    print*,'sum(ZM_ref-ZM) = ', sum(ZM_ref-ZM)



    GBXAR_TMP = GBXAR * (MAPL_RADIUS/1000.)**2 ! transform to km^2
    WHERE (ANGLL < -180)
        ANGLL = 0.0
    END WHERE

    do nrdg = 1, NCAR_NRDG
        KWVRDG(:,:,nrdg) = 0.001/(HWDTH(:,:,nrdg)+0.001)
        EFFRDG(:,:,nrdg) = NCAR_EFFGWORO*(HWDTH(:,:,nrdg)*CLNGT(:,:,nrdg))/GBXAR_TMP
    enddo

    open(newunit=fileID, file=trim(dirName) // "/ANGLL_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) ANGLL_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/KWVRDG_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) KWVRDG_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/EFFRDG_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) EFFRDG_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/GBXAR_TMP_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) GBXAR_TMP_ref
    close(fileID)

    print*,'sum(ANGLL_ref-ANGLL) = ', sum(ANGLL_ref-ANGLL)
    print*,'sum(KWVRDG_ref-KWVRDG) = ', sum(KWVRDG_ref-KWVRDG)
    print*,'sum(EFFRDG_ref-EFFRDG) = ', sum(EFFRDG_ref-EFFRDG)
    print*,'sum(GBXAR_TMP_ref-GBXAR_TMP) = ', sum(GBXAR_TMP_ref-GBXAR_TMP)

    DUDT_GWD_NCAR = 0.0
    DVDT_GWD_NCAR = 0.0
    DTDT_GWD_NCAR = 0.0
    TAUXB_TMP_NCAR = 0.0
    TAUYB_TMP_NCAR = 0.0
    DUDT_ORG_NCAR = 0.0
    DVDT_ORG_NCAR = 0.0
    DTDT_ORG_NCAR = 0.0
    TAUXO_TMP_NCAR = 0.0
    TAUYO_TMP_NCAR = 0.0
    
    if ( (NCAR_EFFGWORO /= 0.0) .OR. (NCAR_EFFGWBKG /= 0.0) ) then
        call gw_intr_ncar(IM*JM,    LM,         DT,     NCAR_NRDG,   &
                    beres_dc_desc, beres_sc_desc, &
                    beres_band, oro_band, rdg_band, &
                    PLE,       T,          U,          V,                   &
                    HT_dc,     HT_sc,      QLDT_mst+QIDT_mst,               &
                    SGH,       MXDIS,      HWDTH,      CLNGT,  ANGLL,       &
                    ANIXY,     GBXAR_TMP,  KWVRDG,     EFFRDG, PREF,        &
                    PMID,      PDEL,       RPDEL,      PILN,   ZM,    LATS, &
                    PHIS,                                                   &
                    DUDT_GWD_NCAR,  DVDT_GWD_NCAR,   DTDT_GWD_NCAR,         &
                    DUDT_ORG_NCAR,  DVDT_ORG_NCAR,   DTDT_ORG_NCAR,         &
                    TAUXO_TMP_NCAR, TAUYO_TMP_NCAR,  &
                    TAUXB_TMP_NCAR, TAUYB_TMP_NCAR,  &
                    NCAR_EFFGWORO, &
                    NCAR_EFFGWBKG, alpha, &
                    RC)
    endif

    open(newunit=fileID, file=trim(dirName) // "/DUDT_GWD_NCAR_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DUDT_GWD_NCAR_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/DVDT_GWD_NCAR_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DVDT_GWD_NCAR_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/DTDT_GWD_NCAR_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DTDT_GWD_NCAR_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/DUDT_ORG_NCAR_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DUDT_ORG_NCAR_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/DVDT_ORG_NCAR_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DVDT_ORG_NCAR_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/DTDT_ORG_NCAR_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DTDT_ORG_NCAR_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/TAUXO_TMP_NCAR_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) TAUXO_TMP_NCAR_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/TAUYO_TMP_NCAR_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) TAUYO_TMP_NCAR_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/TAUXB_TMP_NCAR_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) TAUXB_TMP_NCAR_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/TAUYB_TMP_NCAR_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) TAUYB_TMP_NCAR_ref
    close(fileID)
    
    print*,'sum(DUDT_GWD_NCAR_ref - DUDT_GWD_NCAR) = ', sum(DUDT_GWD_NCAR_ref - DUDT_GWD_NCAR)
    print*,'sum(DVDT_GWD_NCAR_ref - DVDT_GWD_NCAR) = ', sum(DVDT_GWD_NCAR_ref - DVDT_GWD_NCAR)
    print*,'sum(DTDT_GWD_NCAR_ref - DTDT_GWD_NCAR) = ', sum(DTDT_GWD_NCAR_ref - DTDT_GWD_NCAR)
    print*,'sum(DUDT_ORG_NCAR_ref - DUDT_ORG_NCAR) = ', sum(DUDT_ORG_NCAR_ref - DUDT_ORG_NCAR)
    print*,'sum(DVDT_ORG_NCAR_ref - DVDT_ORG_NCAR) = ', sum(DVDT_ORG_NCAR_ref - DVDT_ORG_NCAR)
    print*,'sum(DTDT_ORG_NCAR_ref - DTDT_ORG_NCAR) = ', sum(DTDT_ORG_NCAR_ref - DTDT_ORG_NCAR)
    print*,'sum(TAUXO_TMP_NCAR_ref - TAUXO_TMP_NCAR) = ', sum(TAUXO_TMP_NCAR_ref - TAUXO_TMP_NCAR)
    print*,'sum(TAUYO_TMP_NCAR_ref - TAUYO_TMP_NCAR) = ', sum(TAUYO_TMP_NCAR_ref - TAUYO_TMP_NCAR)
    print*,'sum(TAUXB_TMP_NCAR_ref - TAUXB_TMP_NCAR) = ', sum(TAUXB_TMP_NCAR_ref - TAUXB_TMP_NCAR)
    print*,'sum(TAUYB_TMP_NCAR_ref - TAUYB_TMP_NCAR) = ', sum(TAUYB_TMP_NCAR_ref - TAUYB_TMP_NCAR)

    !$acc data copyin(PREF, PDEL, U, V, DUDT_GWD, DVDT_GWD, DTDT_GWD, DUDT_ORG, DVDT_ORG, DTDT_ORG, &
    !$acc             T, Q, PMID, ZM, VARFLT, AREA, &
    !$acc             DUDT_GWD_NCAR, DVDT_GWD_NCAR, DTDT_GWD_NCAR, TAUXB_TMP_NCAR, TAUYB_TMP_NCAR, &
    !$acc             DUDT_ORG_NCAR, DVDT_ORG_NCAR, DTDT_ORG_NCAR, TAUXO_TMP_NCAR, TAUYO_TMP_NCAR) &
    !$acc      copyout(DUDT_TOT, DVDT_TOT, DTDT_TOT, DUDT_RAH, DVDT_RAH, DTDT_RAH, &
    !$acc              PEGWD_X, PEORO_X, PERAY_X, PEBKG_X, KEGWD_X, KEORO_X,  &
    !$acc              KERAY_X, KEBKG_X, KERES_X, BKGERR_X, &
    !$acc              DUDT_TOFD, DVDT_TOFD) &
    !$acc      create(THV, Hefold, a2, DUDT_GWD_GEOS, DVDT_GWD_GEOS, DTDT_GWD_GEOS, &
    !$acc             TAUXB_TMP_GEOS, TAUYB_TMP_GEOS, DUDT_ORG_GEOS, DVDT_ORG_GEOS, &
    !$acc             DTDT_ORG_GEOS, TAUXO_TMP_GEOS, TAUYO_TMP_GEOS)

    ! Use GEOS GWD only for Extratropical background sources...
    !$acc kernels
    DUDT_GWD_GEOS = 0.0
    DVDT_GWD_GEOS = 0.0
    DTDT_GWD_GEOS = 0.0
    TAUXB_TMP_GEOS = 0.0
    TAUYB_TMP_GEOS = 0.0
    DUDT_ORG_GEOS = 0.0
    DVDT_ORG_GEOS = 0.0
    DTDT_ORG_GEOS = 0.0
    TAUXO_TMP_GEOS = 0.0
    TAUYO_TMP_GEOS = 0.0

    ! *** Bill Putnam mentioned that gw_intr won't be used as much, so I'm
    ! *** not porting it at the moment.
    ! if ( (GEOS_EFFGWORO /= 0.0) .OR. (GEOS_EFFGWBKG /= 0.0) ) then
    !     call gw_intr   (IM*JM,      LM,         DT,                  &
    !          GEOS_PGWV,                                              &
    !          PLE,       T,          U,          V,      SGH,   PREF, &
    !          PMID,      PDEL,       RPDEL,      PILN,   ZM,    LATS, &
    !          DUDT_GWD_GEOS,  DVDT_GWD_GEOS,   DTDT_GWD_GEOS,         &
    !          DUDT_ORG_GEOS,  DVDT_ORG_GEOS,   DTDT_ORG_GEOS,         &
    !          TAUXO_TMP_GEOS, TAUYO_TMP_GEOS,  TAUXO_3D,   TAUYO_3D,  FEO_3D,   &
    !          TAUXB_TMP_GEOS, TAUYB_TMP_GEOS,  TAUXB_3D,   TAUYB_3D,  FEB_3D,   &
    !          FEPO_3D,   FEPB_3D,    DUBKGSRC,   DVBKGSRC,  DTBKGSRC, &
    !          GEOS_BGSTRESS, &
    !          GEOS_EFFGWORO, &
    !          GEOS_EFFGWBKG, &
    !          RC)
    !    endif

    ! ! Total
    DUDT_GWD=DUDT_GWD_GEOS+DUDT_GWD_NCAR
    DVDT_GWD=DVDT_GWD_GEOS+DVDT_GWD_NCAR
    DTDT_GWD=DTDT_GWD_GEOS+DTDT_GWD_NCAR
    ! Background
    TAUXB_TMP=TAUXB_TMP_GEOS+TAUXB_TMP_NCAR
    TAUYB_TMP=TAUYB_TMP_GEOS+TAUYB_TMP_NCAR
    ! Orographic
    DUDT_ORG=DUDT_ORG_GEOS+DUDT_ORG_NCAR
    DVDT_ORG=DVDT_ORG_GEOS+DVDT_ORG_NCAR
    DTDT_ORG=DTDT_ORG_GEOS+DTDT_ORG_NCAR
    TAUXO_TMP=TAUXO_TMP_GEOS+TAUXO_TMP_NCAR
    TAUYO_TMP=TAUYO_TMP_GEOS+TAUYO_TMP_NCAR
    !$acc end kernels

    if (effbeljaars > 0.0) then
        !$acc kernels
        THV = T * (1.0 + MAPL_VIREPS * Q) / ( (PMID/MAPL_P00)**MAPL_KAPPA )
        !$acc end kernels

        !$acc parallel loop gang vector collapse(2) &
        !$acc private(ikpbl, tcrib)
        DO J=1,JM
            DO I=1,IM    
                ! Find the PBL height
                ikpbl = LM
                do L=LM-1,1,-1
                    tcrib = MAPL_GRAV*(THV(I,J,L)-THV(I,J,LM))*ZM(I,J,L)/ &
                            (THV(I,J,LM)*MAX(U(I,J,L)**2+V(I,J,L)**2,1.0E-8))
                    if (tcrib >= 0.25) then
                        ikpbl = L
                        exit
                    end if
                end do
                ! determine the efolding height
                a2(i,j)=effbeljaars * 1.08371722e-7 * VARFLT(i,j) * &
                        MAX(0.0,MIN(1.0,dxmax_ss*(1.-dxmin_ss/SQRT(AREA(i,j))/(dxmax_ss-dxmin_ss))))
                ! Revise e-folding height based on PBL height and topographic std. dev.
                Hefold(i,j) = 1500.0 !MIN(MAX(2*SQRT(VARFLT(i,j)),ZM(i,j,ikpbl)),1500.)
            END DO
        END DO
        !$acc end parallel

        !$acc parallel loop gang vector collapse(3) &
        !$acc private(var_temp, wsp)
        DO L=1, LM
            DO J=1,JM
                DO I=1,IM
                    var_temp = 0.0
                    if (a2(i,j) > 0.0 .AND. ZM(I,J,L) < 4.0*Hefold(i,j)) then
                        wsp      = SQRT(U(i,j,l)**2 + V(i,j,l)**2)
                        wsp      = SQRT(MIN(wsp/25.0,1.0))*MAX(25.0,wsp) ! enhance winds below 25 m/s
                        var_temp = ZM(I,J,L)/Hefold(i,j)
                        var_temp = exp(-var_temp*sqrt(var_temp))*(var_temp**(-1.2))
                        var_temp = wsp*a2(i,j)*(var_temp/Hefold(i,j))
                        !  Note:  This is a semi-implicit treatment of the time differencing
                        !  per Beljaars et al. (2004, QJRMS) doi: 10.1256/qj.03.73
                        DUDT_TOFD(i,j,l) = - var_temp*U(i,j,l)/(1. + var_temp*DT)
                        DVDT_TOFD(i,j,l) = - var_temp*V(i,j,l)/(1. + var_temp*DT)
                        ! Apply Tendency Limiter
                        if (abs(DUDT_TOFD(i,j,l)) > limbeljaars) then
                            DUDT_TOFD(i,j,l) = (limbeljaars/abs(DUDT_TOFD(i,j,l))) * DUDT_TOFD(i,j,l)
                        end if
                        if (abs(DVDT_TOFD(i,j,l)) > limbeljaars) then
                            DVDT_TOFD(i,j,l) = (limbeljaars/abs(DVDT_TOFD(i,j,l))) * DVDT_TOFD(i,j,l)
                        end if
                    else
                        DUDT_TOFD(i,j,l) = 0.0
                        DVDT_TOFD(i,j,l) = 0.0
                    end if
                END DO
            END DO
        END DO
        !$acc end parallel

        !$acc kernels
        DUDT_GWD=DUDT_GWD+DUDT_TOFD
        DVDT_GWD=DVDT_GWD+DVDT_TOFD
        !$acc end kernels
        !deallocate( THV )
    else
        !$acc kernels
        DUDT_TOFD=0.0
        DVDT_TOFD=0.0
        !$acc end kernels
    endif

    open(newunit=fileID, file=trim(dirName) // "/DUDT_TOFD_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DUDT_TOFD_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/DVDT_TOFD_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DVDT_TOFD_ref
    close(fileID)

    !$acc update host(DUDT_TOFD, DVDT_TOFD)

    print*,'sum(DUDT_TOFD_ref - DUDT_TOFD) = ', sum(DUDT_TOFD_ref - DUDT_TOFD)
    print*,'sum(DVDT_TOFD_ref - DVDT_TOFD) = ', sum(DVDT_TOFD_ref - DVDT_TOFD)

    CALL POSTINTR(IM*JM, LM, DT, H0, HH, Z1, TAU1, &
        PREF,     &
        PDEL,     &
        U,        &
        V,        &
        DUDT_GWD, &
        DVDT_GWD, &
        DTDT_GWD, &
        DUDT_ORG, &
        DVDT_ORG, &
        DTDT_ORG, &
        DUDT_TOT, &
        DVDT_TOT, &
        DTDT_TOT, &
        DUDT_RAH, &
        DVDT_RAH, &
        DTDT_RAH, &
        PEGWD_X,  &
        PEORO_X,  &
        PERAY_X,  &
        PEBKG_X,  &
        KEGWD_X,  &
        KEORO_X,  &
        KERAY_X,  &
        KEBKG_X,  &
        KERES_X,  &
        BKGERR_X  )

!$acc end data

    open(newunit=fileID, file=trim(dirName) // "/DUDT_TOT_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DUDT_TOT_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/DVDT_TOT_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DVDT_TOT_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/DTDT_TOT_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DTDT_TOT_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/DUDT_RAH_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DUDT_RAH_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/DVDT_RAH_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DVDT_RAH_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/DTDT_RAH_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) DTDT_RAH_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/PEGWD_X_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) PEGWD_X_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/PEORO_X_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) PEORO_X_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/PERAY_X_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) PERAY_X_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/PEBKG_X_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) PEBKG_X_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/KEGWD_X_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) KEGWD_X_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/KEORO_X_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) KEORO_X_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/KERAY_X_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) KERAY_X_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/KEBKG_X_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) KEBKG_X_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/KERES_X_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) KERES_X_ref
    close(fileID)

    open(newunit=fileID, file=trim(dirName) // "/BKGERR_X_"// trim(rank) // ".out", status='old', form="unformatted", action="read")
    read(fileID) BKGERR_X_ref
    close(fileID)

    print*,'sum(DUDT_TOT_ref - DUDT_TOT) = ', sum(DUDT_TOT_ref - DUDT_TOT)
    print*,'sum(DVDT_TOT_ref - DVDT_TOT) = ', sum(DVDT_TOT_ref - DVDT_TOT)
    print*,'sum(DTDT_TOT_ref - DTDT_TOT) = ', sum(DTDT_TOT_ref - DTDT_TOT)
    print*,'sum(DUDT_RAH_ref - DUDT_RAH) = ', sum(DUDT_RAH_ref - DUDT_RAH)
    print*,'sum(DVDT_RAH_ref - DVDT_RAH) = ', sum(DVDT_RAH_ref - DVDT_RAH)
    print*,'sum(DTDT_RAH_ref - DTDT_RAH) = ', sum(DTDT_RAH_ref - DTDT_RAH)
    print*,'sum(PEGWD_X_ref - PEGWD_X) = ', sum(PEGWD_X_ref - PEGWD_X)
    print*,'sum(PEORO_X_ref - PEORO_X) = ', sum(PEORO_X_ref - PEORO_X)
    print*,'sum(PERAY_X_ref - PERAY_X) = ', sum(PERAY_X_ref - PERAY_X)
    print*,'sum(PEBKG_X_ref - PEBKG_X) = ', sum(PEBKG_X_ref - PEBKG_X)
    print*,'sum(KEGWD_X_ref - KEGWD_X) = ', sum(KEGWD_X_ref - KEGWD_X)
    print*,'sum(KEORO_X_ref - KEORO_X) = ', sum(KEORO_X_ref - KEORO_X)
    print*,'sum(KERAY_X_ref - KERAY_X) = ', sum(KERAY_X_ref - KERAY_X)
    print*,'sum(KEBKG_X_ref - KEBKG_X) = ', sum(KEBKG_X_ref - KEBKG_X)
    print*,'sum(KERES_X_ref - KERES_X) = ', sum(KERES_X_ref - KERES_X)
    print*,'sum(BKGERR_X_ref - BKGERR_X) = ', sum(BKGERR_X_ref - BKGERR_X)
    
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

subroutine postintr(pcols,pver,dt, h0, hh, z1, tau1, &
        pref, &
        pdel, &
        u, &
        v, &
        dudt_gwd, &
        dvdt_gwd, &
        dtdt_gwd, &
        dudt_org, &
        dvdt_org, &
        dtdt_org, &

        ! Outputs
        dudt_tot, &
        dvdt_tot, &
        dtdt_tot, &
        dudt_rah, &
        dvdt_rah, &
        dtdt_rah, &
        pegwd, &
        peoro, &
        peray, &
        pebkg, &
        kegwd, &
        keoro, &
        keray, &
        kebkg, &
        keres, &
        bkgerr )

    implicit none

    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !

    integer, intent(in) :: PCOLS ! Number of longitudes
    integer, intent(in) :: PVER  ! Number of vertical layers
    real,    intent(in) :: DT    ! Time step
    real,    intent(in) :: H0, HH, Z1, TAU1 ! Rayleigh friction parameters

    real,    intent(in) :: PREF(PVER+1)
    real,    intent(in) :: PDEL(PCOLS,PVER)
    real,    intent(in) :: U(PCOLS,PVER)
    real,    intent(in) :: V(PCOLS,PVER)

    real,    intent(in) :: DUDT_GWD(PCOLS,PVER)
    real,    intent(in) :: DVDT_GWD(PCOLS,PVER)
    real,    intent(in) :: DTDT_GWD(PCOLS,PVER)
    real,    intent(in) :: DUDT_ORG(PCOLS,PVER)
    real,    intent(in) :: DVDT_ORG(PCOLS,PVER)
    real,    intent(in) :: DTDT_ORG(PCOLS,PVER)

    real,    intent(out) :: DUDT_TOT(PCOLS,PVER)
    real,    intent(out) :: DVDT_TOT(PCOLS,PVER)
    real,    intent(out) :: DTDT_TOT(PCOLS,PVER)
    real,    intent(out) :: DUDT_RAH(PCOLS,PVER)
    real,    intent(out) :: DVDT_RAH(PCOLS,PVER)
    real,    intent(out) :: DTDT_RAH(PCOLS,PVER)
    real,    intent(out) :: PEGWD(PCOLS)
    real,    intent(out) :: PEORO(PCOLS)
    real,    intent(out) :: PERAY(PCOLS)
    real,    intent(out) :: PEBKG(PCOLS)
    real,    intent(out) :: KEGWD(PCOLS)
    real,    intent(out) :: KEORO(PCOLS)
    real,    intent(out) :: KERAY(PCOLS)
    real,    intent(out) :: KEBKG(PCOLS)
    real,    intent(out) :: KERES(PCOLS)
    real,    intent(out) :: BKGERR(PCOLS)

    !
    !---------------------------Local variables-----------------------------
    !
    integer :: i,k
    real :: zref, kray, ts, te
    !
    !-----------------------------------------------------------------------
    !

    call cpu_time(ts)

!$acc parallel loop gang copyin(PCOLS)
    DO I = 1, PCOLS
        PEGWD(I)  = 0.0
        PEORO(I)  = 0.0
        PERAY(I)  = 0.0
        PEBKG(I)  = 0.0
        KEGWD(I)  = 0.0
        KEORO(I)  = 0.0
        KERAY(I)  = 0.0
        KEBKG(I)  = 0.0
        KERES(I)  = 0.0
        BKGERR(I) = 0.0
    enddo
!$acc end parallel

!$acc parallel loop gang vector collapse(2) &
!$acc private(ZREF, KRAY)
    DO K = 1, PVER
        DO I = 1, PCOLS

            ! Rayleigh friction
            !------------------
            if (TAU1 > 0.0) then
            ZREF     = H0 * LOG(MAPL_P00/(0.5*(PREF(K)+PREF(K+1))))
            KRAY     = (1.0/TAU1)*( 1.0 - TANH( (Z1-ZREF)/HH ) )
            KRAY     = KRAY/(1+DT*KRAY)
            DUDT_RAH(I,K) = -U(I,K)*KRAY
            DVDT_RAH(I,K) = -V(I,K)*KRAY
            DTDT_RAH(I,K) = - ((U(I,K) + (0.5*DT)*DUDT_RAH(I,K))*DUDT_RAH(I,K) + &
                                (V(I,K) + (0.5*DT)*DVDT_RAH(I,K))*DVDT_RAH(I,K)   ) * (1.0/MAPL_CP)
            else
            DUDT_RAH(I,K) = 0.0
            DVDT_RAH(I,K) = 0.0
            DTDT_RAH(I,K) = 0.0
            endif

            DUDT_TOT(I,K) = DUDT_RAH(I,K) + DUDT_GWD(I,K)
            DVDT_TOT(I,K) = DVDT_RAH(I,K) + DVDT_GWD(I,K)
            DTDT_TOT(I,K) = DTDT_RAH(I,K) + DTDT_GWD(I,K)

            ! KE diagnostics
            !----------------
            !$acc atomic update
            PEGWD(I) = PEGWD(I) +  DTDT_TOT(I,K)               *PDEL(I,K)*(MAPL_CP/MAPL_GRAV)
            !$acc end atomic

            !$acc atomic update
            PEORO(I) = PEORO(I) +  DTDT_ORG(I,K)               *PDEL(I,K)*(MAPL_CP/MAPL_GRAV)
            !$acc end atomic
            
            !$acc atomic update
            PERAY(I) = PERAY(I) +  DTDT_RAH(I,K)               *PDEL(I,K)*(MAPL_CP/MAPL_GRAV)
            !$acc end atomic

            !$acc atomic update
            PEBKG(I) = PEBKG(I) + (DTDT_GWD(I,K)-DTDT_ORG(I,K))*PDEL(I,K)*(MAPL_CP/MAPL_GRAV)
            !$acc end atomic

            !$acc atomic update
            KEGWD(I) = KEGWD(I) + ((U(I,K)+(0.5*DT)*DUDT_TOT(I,K))*DUDT_TOT(I,K) +   &
                                    (V(I,K)+(0.5*DT)*DVDT_TOT(I,K))*DVDT_TOT(I,K) ) * PDEL(I,K)*(1.0/MAPL_GRAV)
            !$acc end atomic

            !$acc atomic update
            KEORO(I) = KEORO(I) + ((U(I,K)+(0.5*DT)*DUDT_ORG(I,K))*DUDT_ORG(I,K) +   &
                                    (V(I,K)+(0.5*DT)*DVDT_ORG(I,K))*DVDT_ORG(I,K) ) * PDEL(I,K)*(1.0/MAPL_GRAV)
            !$acc end atomic

            !$acc atomic update
            KERAY(I) = KERAY(I) + ((U(I,K)+(0.5*DT)*DUDT_RAH(I,K))*DUDT_RAH(I,K) +   &
                                    (V(I,K)+(0.5*DT)*DVDT_RAH(I,K))*DVDT_RAH(I,K) ) * PDEL(I,K)*(1.0/MAPL_GRAV)
            !$acc end atomic
            
            !$acc atomic update
            KEBKG(I) = KEBKG(I) + ((U(I,K)+(0.5*DT)*(DUDT_GWD(I,K) - DUDT_ORG(I,K)))*(DUDT_GWD(I,K) - DUDT_ORG(I,K)) +     &
                                    (V(I,K)+(0.5*DT)*(DVDT_GWD(I,K) - DVDT_ORG(I,K)))*(DVDT_GWD(I,K) - DVDT_ORG(I,K))   ) * &
                                    PDEL(I,K)*(1.0/MAPL_GRAV)
            !$acc end atomic
        END DO
    END DO
!$acc end parallel

!$acc parallel loop gang
    DO I = 1, PCOLS
        BKGERR(I) = -( PEBKG(I) + KEBKG(I) )
        KERES(I)  =    PEGWD(I) + KEGWD(I) + BKGERR(I)
    enddo
!$acc end parallel

    call cpu_time(te)

    print*,'postintr time = ', te-ts
end subroutine postintr
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
