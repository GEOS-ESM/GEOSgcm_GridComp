module pymoist_interface_mod

   use iso_c_binding, only: c_int, c_float, c_double, c_bool, c_ptr

   implicit none

   private
   public :: pymoist_interface_f_init
   public :: gfdl_1m_interface_f_init
   public :: pymoist_interface_f_run_AerActivation
   public :: pymoist_interface_f_run_GFDL_1M
   public :: pymoist_interface_f_finalize
   public :: make_moist_flags_C_interop
   public :: make_gfdl_1m_flags_C_interop
   public :: moist_flags_interface_type
   public :: gfdl_1m_flags_interface_type

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
      ! Magic number
      integer(kind=c_int) :: make_flags_C_interop = 123456789
   end type

   type, bind(c) :: gfdl_1m_flags_interface_type
      ! GFDL_1M driver
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
      logical(kind=c_bool) :: do_qa
      logical(kind=c_bool) :: fix_negative
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
      logical(kind=c_bool) :: fast_sat_adj
      real(kind=c_float) :: rh_inc
      real(kind=c_float) :: rh_ins
      real(kind=c_float) :: rh_inr
      logical(kind=c_bool) :: const_vi
      logical(kind=c_bool) :: const_vs
      logical(kind=c_bool) :: const_vg
      logical(kind=c_bool) :: const_vr
      logical(kind=c_bool) :: use_ccn
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
      logical(kind=c_bool) :: do_bigg
      logical(kind=c_bool) :: do_evap
      logical(kind=c_bool) :: do_subl
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
      logical(kind=c_bool) :: z_slope_liq
      logical(kind=c_bool) :: z_slope_ice
      logical(kind=c_bool) :: prog_ccn
      real(kind=c_float) :: c_cracw
      real(kind=c_float) :: alin
      real(kind=c_float) :: clin
      logical(kind=c_bool) :: preciprad
      real(kind=c_float) :: cld_min
      logical(kind=c_bool) :: use_ppm
      logical(kind=c_bool) :: mono_prof
      logical(kind=c_bool) :: do_sedi_heat
      logical(kind=c_bool) :: sedi_transport
      logical(kind=c_bool) :: do_sedi_w
      logical(kind=c_bool) :: de_ice
      integer(kind=c_int) :: icloud_f
      integer(kind=c_int) :: irain_f
      logical(kind=c_bool) :: mp_print
      logical(kind=c_bool) :: use_bergeron
      ! Magic number
      integer(kind=c_int) :: make_flags_C_interop = 123456789
   end type

   interface

      subroutine pymoist_interface_f_init(IMPORT, EXPORT, INTERNAL, MAPL_COMP, moist_flags) bind(c, name='pymoist_interface_c_init')

         import c_ptr, c_int, c_float, c_double, c_bool, moist_flags_interface_type

         implicit none
         type(moist_flags_interface_type), intent(in) :: moist_flags
         type(c_ptr), intent(in) :: IMPORT, EXPORT, INTERNAL, MAPL_COMP

      end subroutine pymoist_interface_f_init

      subroutine gfdl_1m_interface_f_init(gfdl_1m_flags) bind(c, name='gfdl_1m_interface_c_init')

         import c_int, c_float, c_double, c_bool, gfdl_1m_flags_interface_type

         implicit none
         type(gfdl_1m_flags_interface_type), intent(in) :: gfdl_1m_flags

      end subroutine gfdl_1m_interface_f_init

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

      subroutine pymoist_interface_f_run_GFDL_1M() bind(c, name='pymoist_interface_c_run_GFDL_1M')
      end subroutine pymoist_interface_f_run_GFDL_1M

      subroutine pymoist_interface_f_finalize() bind(c, name='pymoist_interface_c_finalize')
      end subroutine pymoist_interface_f_finalize

   end interface

contains

   subroutine make_moist_flags_C_interop(npx, npy, npz , nx, ny, n_tiles, n_modes, moist_flags)

      integer, intent(in) :: nx, ny, npx, npy, npz, n_tiles
      ! Aer Activation
      integer, intent(in) :: n_modes


      type(moist_flags_interface_type), intent(out) :: moist_flags

      moist_flags%npx = npx
      moist_flags%npy = npy
      moist_flags%npz = npz
      moist_flags%layout_x = nx
      moist_flags%layout_y = ny
      moist_flags%n_tiles = n_tiles
      moist_flags%n_modes = n_modes

   end subroutine make_moist_flags_C_interop

   subroutine make_gfdl_1m_flags_C_interop(dt_moist, mp_time, t_min, t_sub, tau_r2g, tau_smlt, tau_g2r, &
      dw_land, dw_ocean, vi_fac, vr_fac, vs_fac, vg_fac, ql_mlt, do_qa, fix_negative, vi_max, vs_max, &
      vg_max, vr_max, qs_mlt, qs0_crt, qi_gen, ql0_max, qi0_max, qi0_crt, qr0_crt, fast_sat_adj, rh_inc, &
      rh_ins, rh_inr, const_vi, const_vs, const_vg, const_vr, use_ccn, rthreshu, rthreshs, ccn_l, ccn_o, &
      qc_crt, tau_g2v, tau_v2g, tau_s2v, tau_v2s, tau_revp, tau_frz, do_bigg, do_evap, do_subl, sat_adj0, &
      c_piacr, tau_imlt, tau_v2l, tau_l2v, tau_i2v, tau_i2s, tau_l2r, qi_lim, ql_gen, c_paut, c_psaci, &
      c_pgacs, c_pgaci, z_slope_liq, z_slope_ice, prog_ccn, c_cracw, alin, clin, preciprad, cld_min, &
      use_ppm, mono_prof, do_sedi_heat, sedi_transport, do_sedi_w, de_ice, icloud_f, irain_f, mp_print, &
      use_bergeron, gfdl_1m_flags)

      ! GFDL_1M driver
      logical, intent(in) :: use_bergeron, do_qa, fix_negative, fast_sat_adj
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
      real, intent(in) :: alin, clin, cld_min
      integer, intent(in) :: icloud_f, irain_f

      type(gfdl_1m_flags_interface_type), intent(out) :: gfdl_1m_flags

      gfdl_1m_flags%dt_moist = dt_moist
      gfdl_1m_flags%mp_time = mp_time
      gfdl_1m_flags%t_min = t_min
      gfdl_1m_flags%t_sub = t_sub
      gfdl_1m_flags%tau_r2g = tau_r2g
      gfdl_1m_flags%tau_smlt = tau_smlt
      gfdl_1m_flags%tau_g2r = tau_g2r
      gfdl_1m_flags%dw_land = dw_land
      gfdl_1m_flags%dw_ocean = dw_ocean
      gfdl_1m_flags%vi_fac = vi_fac
      gfdl_1m_flags%vr_fac = vr_fac
      gfdl_1m_flags%vs_fac = vs_fac
      gfdl_1m_flags%vg_fac = vg_fac
      gfdl_1m_flags%ql_mlt = ql_mlt
      gfdl_1m_flags%do_qa = do_qa
      gfdl_1m_flags%fix_negative = fix_negative
      gfdl_1m_flags%vi_max = vi_max
      gfdl_1m_flags%vs_max = vs_max
      gfdl_1m_flags%vg_max = vg_max
      gfdl_1m_flags%vr_max = vr_max
      gfdl_1m_flags%qs_mlt = qs_mlt
      gfdl_1m_flags%qs0_crt = qs0_crt
      gfdl_1m_flags%qi_gen = qi_gen
      gfdl_1m_flags%ql0_max = ql0_max
      gfdl_1m_flags%qi0_max = qi0_max
      gfdl_1m_flags%qi0_crt = qi0_crt
      gfdl_1m_flags%qr0_crt = qr0_crt
      gfdl_1m_flags%fast_sat_adj = fast_sat_adj
      gfdl_1m_flags%rh_inc = rh_inc
      gfdl_1m_flags%rh_ins = rh_ins
      gfdl_1m_flags%rh_inr = rh_inr
      gfdl_1m_flags%const_vi = const_vi
      gfdl_1m_flags%const_vs = const_vs
      gfdl_1m_flags%const_vg = const_vg
      gfdl_1m_flags%const_vr = const_vr
      gfdl_1m_flags%use_ccn = use_ccn
      gfdl_1m_flags%rthreshu = rthreshu
      gfdl_1m_flags%rthreshs = rthreshs
      gfdl_1m_flags%ccn_l = ccn_l
      gfdl_1m_flags%ccn_o = ccn_o
      gfdl_1m_flags%qc_crt = qc_crt
      gfdl_1m_flags%tau_g2v = tau_g2v
      gfdl_1m_flags%tau_v2g = tau_v2g
      gfdl_1m_flags%tau_s2v = tau_s2v
      gfdl_1m_flags%tau_v2s = tau_v2s
      gfdl_1m_flags%tau_revp = tau_revp
      gfdl_1m_flags%tau_frz = tau_frz
      gfdl_1m_flags%do_bigg = do_bigg
      gfdl_1m_flags%do_evap = do_evap
      gfdl_1m_flags%do_subl = do_subl
      gfdl_1m_flags%sat_adj0 = sat_adj0
      gfdl_1m_flags%c_piacr = c_piacr
      gfdl_1m_flags%tau_imlt = tau_imlt
      gfdl_1m_flags%tau_v2l = tau_v2l
      gfdl_1m_flags%tau_l2v = tau_l2v
      gfdl_1m_flags%tau_i2v = tau_i2v
      gfdl_1m_flags%tau_i2s = tau_i2s
      gfdl_1m_flags%tau_l2r = tau_l2r
      gfdl_1m_flags%qi_lim = qi_lim
      gfdl_1m_flags%ql_gen = ql_gen
      gfdl_1m_flags%c_paut = c_paut
      gfdl_1m_flags%c_psaci = c_psaci
      gfdl_1m_flags%c_pgacs = c_pgacs
      gfdl_1m_flags%c_pgaci = c_pgaci
      gfdl_1m_flags%z_slope_liq = z_slope_liq
      gfdl_1m_flags%z_slope_ice = z_slope_ice
      gfdl_1m_flags%prog_ccn = prog_ccn
      gfdl_1m_flags%c_cracw = c_cracw
      gfdl_1m_flags%alin = alin
      gfdl_1m_flags%clin = clin
      gfdl_1m_flags%preciprad = preciprad
      gfdl_1m_flags%cld_min = cld_min
      gfdl_1m_flags%use_ppm = use_ppm
      gfdl_1m_flags%mono_prof = mono_prof
      gfdl_1m_flags%do_sedi_heat = do_sedi_heat
      gfdl_1m_flags%sedi_transport = sedi_transport
      gfdl_1m_flags%do_sedi_w = do_sedi_w
      gfdl_1m_flags%de_ice = de_ice
      gfdl_1m_flags%icloud_f = icloud_f
      gfdl_1m_flags%irain_f = irain_f
      gfdl_1m_flags%mp_print = mp_print
      gfdl_1m_flags%use_bergeron = use_bergeron

   end subroutine make_gfdl_1m_flags_C_interop

end module pymoist_interface_mod
