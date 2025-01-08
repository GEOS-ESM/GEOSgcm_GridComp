module pymoist_interface_mod

   use iso_c_binding, only: c_int, c_float, c_double, c_bool

   implicit none

   private
   public :: pymoist_interface_f_init
   public :: pymoist_interface_f_run_AerActivation
   public :: pymoist_interface_f_run_GFDL1M
   public :: pymoist_interface_f_run_gfdl1m_driver
   public :: pymoist_interface_f_finalize
   public :: make_moist_flags_C_interop
   public :: moist_flags_interface_type

   !-----------------------------------------------------------------------
   ! Shadow C interoperable config struct for FV. See `fv_arrays.f90` for
   ! the original structure, docs and default values
   !-----------------------------------------------------------------------
   type, bind(c) :: moist_flags_interface_type
      ! Grid information
      integer(kind=c_int) :: npx
      integer(kind=c_int) :: npy
      integer(kind=c_int) :: npz
      integer(kind=c_int) :: layout_x
      integer(kind=c_int) :: layout_y
      integer(kind=c_int) :: n_tiles
      ! Aer Activation
      integer(kind=c_int) :: n_modes
      ! GFDL_1M driver
      logical(kind=c_bool) :: phys_hydrostatic
      logical(kind=c_bool) :: hydrostatic
      logical(kind=c_bool) :: do_qa
      logical(kind=c_bool) :: fix_negative
      logical(kind=c_bool) :: fast_sat_adj
      logical(kind=c_bool) :: const_vi
      logical(kind=c_bool) :: const_vs
      logical(kind=c_bool) :: const_vg
      logical(kind=c_bool) :: const_vr
      logical(kind=c_bool) :: use_ccn
      logical(kind=c_bool) :: do_bigg
      logical(kind=c_bool) :: do_evap
      logical(kind=c_bool) :: do_subl
      logical(kind=c_bool) :: z_slope_liq
      logical(kind=c_bool) :: z_slope_ice
      logical(kind=c_bool) :: prog_ccn
      logical(kind=c_bool) :: preciprad
      logical(kind=c_bool) :: use_ppm
      logical(kind=c_bool) :: mono_prof
      logical(kind=c_bool) :: do_sedi_heat
      logical(kind=c_bool) :: sedi_transport
      logical(kind=c_bool) :: do_sedi_w
      logical(kind=c_bool) :: de_ice
      logical(kind=c_bool) :: mp_print
      real(kind=c_float) :: dt_moist
      real(kind=c_float) :: mp_time
      real(kind=c_float) :: t_min
      real(kind=c_float) :: t_sub
      real(kind=c_float) :: tau_r2g
      real(kind=c_float) :: tau_smlt
      real(kind=c_float) :: tau_g2r
      real(kind=c_float) :: dw_land
      real(kind=c_float) :: dw_ocean
      real(kind=c_float) :: vi_fac
      real(kind=c_float) :: vr_fac
      real(kind=c_float) :: vs_fac
      real(kind=c_float) :: vg_fac
      real(kind=c_float) :: ql_mlt
      real(kind=c_float) :: vi_max
      real(kind=c_float) :: vs_max
      real(kind=c_float) :: vg_max
      real(kind=c_float) :: vr_max
      real(kind=c_float) :: qs_mlt
      real(kind=c_float) :: qs0_crt
      real(kind=c_float) :: qi_gen
      real(kind=c_float) :: ql0_max
      real(kind=c_float) :: qi0_max
      real(kind=c_float) :: qi0_crt
      real(kind=c_float) :: qr0_crt
      real(kind=c_float) :: rh_inc
      real(kind=c_float) :: rh_ins
      real(kind=c_float) :: rh_inr
      real(kind=c_float) :: rthreshu
      real(kind=c_float) :: rthreshs
      real(kind=c_float) :: ccn_l
      real(kind=c_float) :: ccn_o
      real(kind=c_float) :: qc_crt
      real(kind=c_float) :: tau_g2v
      real(kind=c_float) :: tau_v2g
      real(kind=c_float) :: tau_s2v
      real(kind=c_float) :: tau_v2s
      real(kind=c_float) :: tau_revp
      real(kind=c_float) :: tau_frz
      real(kind=c_float) :: sat_adj0
      real(kind=c_float) :: c_piacr
      real(kind=c_float) :: tau_imlt
      real(kind=c_float) :: tau_v2l
      real(kind=c_float) :: tau_l2v
      real(kind=c_float) :: tau_i2v
      real(kind=c_float) :: tau_i2s
      real(kind=c_float) :: tau_l2r
      real(kind=c_float) :: qi_lim
      real(kind=c_float) :: ql_gen
      real(kind=c_float) :: c_paut
      real(kind=c_float) :: c_psaci
      real(kind=c_float) :: c_pgacs
      real(kind=c_float) :: c_pgaci
      real(kind=c_float) :: c_cracw
      real(kind=c_float) :: alin
      real(kind=c_float) :: clin
      real(kind=c_float) :: cld_min
      real(kind=c_float) :: icloud_f
      real(kind=c_float) :: irain_f
      ! Magic number
      integer(kind=c_int) :: make_flags_C_interop = 123456789
   end type


   interface

      subroutine pymoist_interface_f_init( moist_flags) bind(c, name='pymoist_interface_c_init')

         import c_int, c_float, c_double, c_bool, moist_flags_interface_type

         implicit none
         type(moist_flags_interface_type), intent(in) :: moist_flags

      end subroutine pymoist_interface_f_init

      subroutine pymoist_interface_f_run_AerActivation( &
         aero_dgn, aero_num, aero_hygroscopicity, aero_sigma, &
         frland, nn_ocean, nn_land, &
         t, plo, &
         qicn, qils, qlcn, qlls, &
         vvel, tke, &
         nacti, nwfa, nact ) bind(c, name='pymoist_interface_c_run_AerActivation')

         import c_int, c_float, c_double

         implicit none

         ! Input
         real(kind=c_float), dimension(*), intent(in) :: aero_dgn
         real(kind=c_float), dimension(*), intent(in) :: aero_num
         real(kind=c_float), dimension(*), intent(in) :: aero_hygroscopicity
         real(kind=c_float), dimension(*), intent(in) :: aero_sigma

         real(kind=c_float), dimension(*), intent(in) :: frland
         real(kind=c_float), value, intent(in) :: nn_ocean
         real(kind=c_float), value, intent(in) :: nn_land

         real(kind=c_float), dimension(*), intent(in) :: t
         real(kind=c_float), dimension(*), intent(in) :: plo

         real(kind=c_float), dimension(*), intent(in) :: qicn
         real(kind=c_float), dimension(*), intent(in) :: qils
         real(kind=c_float), dimension(*), intent(in) :: qlcn
         real(kind=c_float), dimension(*), intent(in) :: qlls

         real(kind=c_float), dimension(*), intent(in) :: vvel
         real(kind=c_float), dimension(*), intent(in) :: tke

         ! Output
         real(kind=c_float), dimension(*), intent(in) :: nacti
         real(kind=c_float), dimension(*), intent(in) :: nwfa
         real(kind=c_float), dimension(*), intent(in) :: nact

      end subroutine pymoist_interface_f_run_AerActivation

      subroutine pymoist_interface_f_run_GFDL1M( &
         dw_land, dw_ocean, PDFSHAPE, TURNRHCRIT_PARAM, &
         DT_MOIST, CCW_EVAP_EFF, CCI_EVAP_EFF, &
         LMELTFRZ, &
         AREA, CNV_FRC, SRF_TYPE, &
         KLCL, &
         EIS, PLmb, PLEmb, NACTL, NACTI, QST,&
         T, Q, QLCN, QICN, QLLS, QILS, CLLS, CLCN, &
         SUBLC, EVAPC, RHX )  bind(c, name='pymoist_interface_c_run_GFDL1M')

         import c_int, c_float, c_bool

         implicit none

         ! Input

         ! All parameters should be in `init`
         real(kind=c_float), value, intent(in) :: dw_land, dw_ocean ! Namelist flags
         real(kind=c_float), value, intent(in) :: TURNRHCRIT_PARAM, DT_MOIST, CCW_EVAP_EFF, CCI_EVAP_EFF
         integer(kind=c_int), value, intent(in) :: PDFSHAPE, LMELTFRZ


         real(kind=c_float), dimension(*), intent(in) :: AREA, CNV_FRC, SRF_TYPE, EIS, PLmb, PLEmb, NACTL, NACTI, QST
         integer(kind=c_int), dimension(*), intent(in) :: KLCL

         ! InOut
         real(kind=c_float), dimension(*), intent(inout) :: T, Q, QLCN, QICN, QLLS, QILS, CLCN, CLLS

         ! Output
         real(kind=c_float), dimension(*), intent(out) :: SUBLC, EVAPC, RHX

      end subroutine pymoist_interface_f_run_GFDL1M

      subroutine pymoist_interface_f_run_GFDL1M_driver( &
         RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, NACTAll, &
         DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
         DQSDTmic, DQGDTmic, DQADTmic, DTDTmic, &
         T, W, U, V, DUDTmic, DVDTmic, DZ, DP, &
         AREA, DT_MOIST, FRLAND, CNV_FRC, SRF_TYPE, EIS, &
         RHCRIT3D, ANV_ICEFALL, LS_ICEFALL, &
         REV_LS, RSU_LS, &
         PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
         PFL_LS, PFI_LS, &
         LHYDROSTATIC, LPHYS_HYDROSTATIC ) bind(c, name='pymoist_interface_c_run_GFDL1M_driver')

         import c_int, c_float, c_bool

         implicit none

         ! Input
         real(kind=c_float), dimension(*), intent(in) :: RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, NACTAll
         real(kind=c_float), dimension(*), intent(in) :: T, W, U, V, DUDTmic, DVDTmic, DZ, DP

         real(kind=c_float), dimension(*), intent(in) :: AREA, FRLAND, CNV_FRC, SRF_TYPE, EIS, RHCRIT3D
         real(kind=c_float), value, intent(in) :: DT_MOIST, LS_ICEFALL, ANV_ICEFALL
         integer(kind=c_int), value, intent(in) :: LHYDROSTATIC, LPHYS_HYDROSTATIC ! Actually bool but LOGICAL(4) and (1) are fighting

         ! Output
         real(kind=c_float), dimension(*), intent(inout) :: DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, DQSDTmic, DQGDTmic, DQADTmic, DTDTmic
         real(kind=c_float), dimension(*), intent(inout) :: REV_LS, RSU_LS
         real(kind=c_float), dimension(*), intent(inout) :: PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL
         real(kind=c_float), dimension(*), intent(inout) :: PFL_LS, PFI_LS

      end subroutine pymoist_interface_f_run_GFDL1M_driver


      subroutine pymoist_interface_f_finalize() bind(c, name='pymoist_interface_c_finalize')
      end subroutine pymoist_interface_f_finalize

   end interface

contains

   subroutine make_moist_flags_C_interop(npx, npy, npz , nx, ny, n_tiles, n_modes, &
      phys_hydrostatic, hydrostatic, do_qa, fix_negative, fast_sat_adj, &
      const_vi, const_vs, const_vg, const_vr, use_ccn, do_bigg, do_evap, &
      do_subl, z_slope_liq, z_slope_ice, prog_ccn, preciprad, use_ppm, &
      mono_prof, do_sedi_heat, sedi_transport, do_sedi_w, de_ice, mp_print, &
      dt_moist, mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, &
      dw_land, dw_ocean, vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, vi_max, &
      vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max, &
      qi0_crt, qr0_crt, rh_inc, rh_ins, rh_inr, rthreshu, rthreshs, ccn_l, &
      ccn_o, qc_crt, tau_g2v, tau_v2g, tau_s2v, tau_v2s, tau_revp, tau_frz, &
      sat_adj0, c_piacr, tau_imlt, tau_v2l, tau_l2v, tau_i2v, tau_i2s, &
      tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs, c_pgaci, c_cracw, &
      alin, clin, cld_min, icloud_f, irain_f, moist_flags)

      integer, intent(in) :: nx, ny, npx, npy, npz, n_tiles
      ! Aer Activation
      integer, intent(in) :: n_modes
      ! GFDL_1M driver
      logical, intent(in) :: phys_hydrostatic, hydrostatic, do_qa, fix_negative, fast_sat_adj
      logical, intent(in) :: const_vi, const_vs, const_vg, const_vr, use_ccn, do_bigg, do_evap
      logical, intent(in) :: do_subl, z_slope_liq, z_slope_ice, prog_ccn, preciprad, use_ppm
      logical, intent(in) :: mono_prof, do_sedi_heat, sedi_transport, do_sedi_w, de_ice, mp_print
      real, intent(in) :: dt_moist, mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r
      real, intent(in) :: dw_land, dw_ocean, vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, vi_max
      real, intent(in) :: vs_max, vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max
      real, intent(in) :: qi0_crt, qr0_crt, rh_inc, rh_ins, rh_inr, rthreshu, rthreshs, ccn_l
      real, intent(in) :: ccn_o, qc_crt, tau_g2v, tau_v2g, tau_s2v, tau_v2s, tau_revp, tau_frz
      real, intent(in) :: sat_adj0, c_piacr, tau_imlt, tau_v2l, tau_l2v, tau_i2v, tau_i2s
      real, intent(in) :: tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, c_pgacs, c_pgaci, c_cracw
      real, intent(in) :: alin, clin, cld_min, icloud_f, irain_f


      type(moist_flags_interface_type), intent(out) :: moist_flags

      moist_flags%npx = npx
      moist_flags%npy = npy
      moist_flags%npz = npz
      moist_flags%layout_x = nx
      moist_flags%layout_y = ny
      moist_flags%n_tiles = n_tiles
      moist_flags%n_modes = n_modes
      moist_flags%phys_hydrostatic = phys_hydrostatic
      moist_flags%hydrostatic = hydrostatic
      moist_flags%do_qa = do_qa
      moist_flags%fix_negative = fix_negative
      moist_flags%fast_sat_adj = fast_sat_adj
      moist_flags%const_vi = const_vi
      moist_flags%const_vs = const_vs
      moist_flags%const_vg = const_vg
      moist_flags%const_vr = const_vr
      moist_flags%use_ccn = use_ccn
      moist_flags%do_bigg = do_bigg
      moist_flags%do_evap = do_evap
      moist_flags%do_subl = do_subl
      moist_flags%z_slope_liq = z_slope_liq
      moist_flags%z_slope_ice = z_slope_ice
      moist_flags%prog_ccn = prog_ccn
      moist_flags%preciprad = preciprad
      moist_flags%use_ppm = use_ppm
      moist_flags%mono_prof = mono_prof
      moist_flags%do_sedi_heat = do_sedi_heat
      moist_flags%sedi_transport = sedi_transport
      moist_flags%do_sedi_w = do_sedi_w
      moist_flags%de_ice = de_ice
      moist_flags%mp_print = mp_print
      moist_flags%dt_moist = dt_moist
      moist_flags%mp_time = mp_time
      moist_flags%t_min = t_min
      moist_flags%t_sub = t_sub
      moist_flags%tau_r2g = tau_r2g
      moist_flags%tau_smlt = tau_smlt
      moist_flags%tau_g2r = tau_g2r
      moist_flags%dw_land = dw_land
      moist_flags%dw_ocean = dw_ocean
      moist_flags%vi_fac = vi_fac
      moist_flags%vr_fac = vr_fac
      moist_flags%vs_fac = vs_fac
      moist_flags%vg_fac = vg_fac
      moist_flags%ql_mlt = ql_mlt
      moist_flags%vi_max = vi_max
      moist_flags%vs_max = vs_max
      moist_flags%vg_max = vg_max
      moist_flags%vr_max = vr_max
      moist_flags%qs_mlt = qs_mlt
      moist_flags%qs0_crt = qs0_crt
      moist_flags%qi_gen = qi_gen
      moist_flags%ql0_max = ql0_max
      moist_flags%qi0_max = qi0_max
      moist_flags%qi0_crt = qi0_crt
      moist_flags%qr0_crt = qr0_crt
      moist_flags%rh_inc = rh_inc
      moist_flags%rh_ins = rh_ins
      moist_flags%rh_inr = rh_inr
      moist_flags%rthreshu = rthreshu
      moist_flags%rthreshs = rthreshs
      moist_flags%ccn_l = ccn_l
      moist_flags%ccn_o = ccn_o
      moist_flags%qc_crt = qc_crt
      moist_flags%tau_g2v = tau_g2v
      moist_flags%tau_v2g = tau_v2g
      moist_flags%tau_s2v = tau_s2v
      moist_flags%tau_v2s = tau_v2s
      moist_flags%tau_revp = tau_revp
      moist_flags%tau_frz = tau_frz
      moist_flags%sat_adj0 = sat_adj0
      moist_flags%c_piacr = c_piacr
      moist_flags%tau_imlt = tau_imlt
      moist_flags%tau_v2l = tau_v2l
      moist_flags%tau_l2v = tau_l2v
      moist_flags%tau_i2v = tau_i2v
      moist_flags%tau_i2s = tau_i2s
      moist_flags%tau_l2r = tau_l2r
      moist_flags%qi_lim = qi_lim
      moist_flags%ql_gen = ql_gen
      moist_flags%c_paut = c_paut
      moist_flags%c_psaci = c_psaci
      moist_flags%c_pgacs = c_pgacs
      moist_flags%c_pgaci = c_pgaci
      moist_flags%c_cracw = c_cracw
      moist_flags%alin = alin
      moist_flags%clin = clin
      moist_flags%cld_min = cld_min
      moist_flags%icloud_f = icloud_f
      moist_flags%irain_f = irain_f

   end subroutine make_moist_flags_C_interop

end module pymoist_interface_mod
