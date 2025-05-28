module pymoist_interface_mod

   use iso_c_binding, only: c_int, c_float, c_double, c_bool

   implicit none

   !> Single precision real numbers, 6 digits, range 10⁻³⁷ to 10³⁷-1; 32 bits
   integer, parameter :: sp = selected_real_kind(6, 37)
   !> Double precision real numbers, 15 digits, range 10⁻³⁰⁷ to 10³⁰⁷-1; 64 bits
   integer, parameter :: dp = selected_real_kind(15, 307)

   private
   public :: pymoist_interface_f_init
   public :: pymoist_interface_f_finalize
   public :: gfdl_1m_interface_f_init
   public :: pymoist_interface_f_run_GFDL1M
   public :: pymoist_interface_f_run_GFDL_1M_driver
   public :: pymoist_interface_f_run_AerActivation
   public :: compute_uwshcu_f_init
   public :: compute_uwshcu_f_run_compute_uwshcu
   public :: make_moist_flags_C_interop
   public :: make_gfdl_1m_flags_C_interop
   public :: moist_flags_interface_type
   public :: gfdl_1m_flags_interface_type

   !!-----------------------------------------------------------------------
   !! Type structs                                                        !! 
   !!-----------------------------------------------------------------------
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
      integer(kind=c_int) :: icloud_f
      integer(kind=c_int) :: irain_f
      ! Magic number
      integer(kind=c_int) :: make_flags_C_interop = 123456789
   end type

   interface

      !!-----------------------------------------------------------------------
      !! Interface: pyMoist                                                  !! 
      !!-----------------------------------------------------------------------
      subroutine pymoist_interface_f_init(moist_flags) bind(c, name='pymoist_interface_c_init')

         import c_int, c_float, c_double, c_bool, moist_flags_interface_type

         implicit none
         type(moist_flags_interface_type), intent(in) :: moist_flags

      end subroutine pymoist_interface_f_init

      subroutine pymoist_interface_f_finalize() bind(c, name='pymoist_interface_c_finalize')
      end subroutine pymoist_interface_f_finalize

      !!-----------------------------------------------------------------------
      !! Interface: GFDL One moment microphysics                             !! 
      !!-----------------------------------------------------------------------

      subroutine gfdl_1m_interface_f_init(gfdl_1m_flags) bind(c, name='gfdl_1m_interface_c_init')

         import c_int, c_float, c_double, c_bool, gfdl_1m_flags_interface_type

         implicit none
         type(gfdl_1m_flags_interface_type), intent(in) :: gfdl_1m_flags

      end subroutine gfdl_1m_interface_f_init

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

      subroutine pymoist_interface_f_run_GFDL_1M_driver( &
         RAD_QV, RAD_QL, RAD_QR, RAD_QI, RAD_QS, RAD_QG, RAD_CF, NACTAll, &
         DQVDTmic, DQLDTmic, DQRDTmic, DQIDTmic, &
         DQSDTmic, DQGDTmic, DQADTmic, DTDTmic, &
         T, W, U, V, DUDTmic, DVDTmic, DZ, DP, &
         AREA, DT_MOIST, FRLAND, CNV_FRC, SRF_TYPE, EIS, &
         RHCRIT3D, ANV_ICEFALL, LS_ICEFALL, &
         REV_LS, RSU_LS, &
         PRCP_RAIN, PRCP_SNOW, PRCP_ICE, PRCP_GRAUPEL, &
         PFL_LS, PFI_LS, &
         LHYDROSTATIC, LPHYS_HYDROSTATIC ) bind(c, name='pymoist_interface_c_run_GFDL_1M_driver')

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

      end subroutine pymoist_interface_f_run_GFDL_1M_driver

      !!-----------------------------------------------------------------------
      !! Interface: Aerosol Activation                                       !! 
      !!-----------------------------------------------------------------------

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

      !!-----------------------------------------------------------------------
      !! Interface: UW Shallow Convection                                    !! 
      !!-----------------------------------------------------------------------

      subroutine compute_uwshcu_f_init( &
         !inputs
         ncnst, &
         k0, &
         windsrcavg &
         !inouts
         !outputs
         ) bind(c, name='compute_uwshcu_c_init')
         import c_int, c_float, c_double
         implicit none

         integer(kind=c_int), value, intent(in) :: ncnst

         integer(kind=c_int), value, intent(in) :: k0

         integer(kind=c_int), value, intent(in) :: windsrcavg

      end subroutine compute_uwshcu_f_init

      subroutine compute_uwshcu_f_run_compute_uwshcu( &
         !inputs
         dotransport, &
         k0, &
         windsrcavg, &
         qtsrchgt, &
         qtsrc_fac, &
         thlsrc_fac, &
         frc_rasn, &
         rbuoy, &
         epsvarw, &
         use_CINcin, &
         mumin1, &
         rmaxfrac, &
         PGFc, &
         dt, &
         niter_xc, &
         criqc, &
         rle, &
         cridist_opt, &
         mixscale, &
         rdrag, &
         rkm, &
         use_self_detrain, &
         detrhgt, &
         use_cumpenent, &
         rpen, &
         use_momenflx, &
         rdrop, &
         iter_cin, &
         pifc0_inv, pifc0_inv_dim_sizes, pifc0_inv_rank, &
         zifc0_inv, zifc0_inv_dim_sizes, zifc0_inv_rank, &
         pmid0_inv, pmid0_inv_dim_sizes, pmid0_inv_rank, &
         zmid0_inv, zmid0_inv_dim_sizes, zmid0_inv_rank, &
         kpbl_inv, kpbl_inv_dim_sizes, kpbl_inv_rank, &
         exnmid0_inv, exnmid0_inv_dim_sizes, exnmid0_inv_rank, &
         exnifc0_inv, exnifc0_inv_dim_sizes, exnifc0_inv_rank, &
         dp0_inv, dp0_inv_dim_sizes, dp0_inv_rank, &
         u0_inv, u0_inv_dim_sizes, u0_inv_rank, &
         v0_inv, v0_inv_dim_sizes, v0_inv_rank, &
         qv0_inv, qv0_inv_dim_sizes, qv0_inv_rank, &
         ql0_inv, ql0_inv_dim_sizes, ql0_inv_rank, &
         qi0_inv, qi0_inv_dim_sizes, qi0_inv_rank, &
         t0_inv, t0_inv_dim_sizes, t0_inv_rank, &
         frland_in, frland_in_dim_sizes, frland_in_rank, &
         tke_inv, tke_inv_dim_sizes, tke_inv_rank, &
         rkfre, rkfre_dim_sizes, rkfre_rank, &
         cush, cush_dim_sizes, cush_rank, &
         shfx, shfx_dim_sizes, shfx_rank, &
         evap, evap_dim_sizes, evap_rank, &
         cnvtr, cnvtr_dim_sizes, cnvtr_rank, &
         CNV_Tracers, CNV_Tracers_dim_sizes, CNV_Tracers_rank, &
         !inouts
         !outputs
         umf_inv, umf_inv_dim_sizes, umf_inv_rank, &
         dcm_inv, dcm_inv_dim_sizes, dcm_inv_rank, &
         qtflx_inv, qtflx_inv_dim_sizes, qtflx_inv_rank, &
         slflx_inv, slflx_inv_dim_sizes, slflx_inv_rank, &
         uflx_inv, uflx_inv_dim_sizes, uflx_inv_rank, &
         vflx_inv, vflx_inv_dim_sizes, vflx_inv_rank, &
         qvten_inv, qvten_inv_dim_sizes, qvten_inv_rank, &
         qlten_inv, qlten_inv_dim_sizes, qlten_inv_rank, &
         qiten_inv, qiten_inv_dim_sizes, qiten_inv_rank, &
         tten_inv, tten_inv_dim_sizes, tten_inv_rank, &
         uten_inv, uten_inv_dim_sizes, uten_inv_rank, &
         vten_inv, vten_inv_dim_sizes, vten_inv_rank, &
         qrten_inv, qrten_inv_dim_sizes, qrten_inv_rank, &
         qsten_inv, qsten_inv_dim_sizes, qsten_inv_rank, &
         cufrc_inv, cufrc_inv_dim_sizes, cufrc_inv_rank, &
         fer_inv, fer_inv_dim_sizes, fer_inv_rank, &
         fdr_inv, fdr_inv_dim_sizes, fdr_inv_rank, &
         ndrop_inv, ndrop_inv_dim_sizes, ndrop_inv_rank, &
         nice_inv, nice_inv_dim_sizes, nice_inv_rank, &
         qldet_inv, qldet_inv_dim_sizes, qldet_inv_rank, &
         qlsub_inv, qlsub_inv_dim_sizes, qlsub_inv_rank, &
         qidet_inv, qidet_inv_dim_sizes, qidet_inv_rank, &
         qisub_inv, qisub_inv_dim_sizes, qisub_inv_rank, &
         tpert_out, tpert_out_dim_sizes, tpert_out_rank, &
         qpert_out, qpert_out_dim_sizes, qpert_out_rank &
         ) bind(c, name='compute_uwshcu_c_run_compute_uwshcu')
         import c_int, c_float, c_double
         implicit none

         integer(kind=c_int), value, intent(in) :: dotransport

         integer(kind=c_int), value, intent(in) :: k0

         integer(kind=c_int), value, intent(in) :: windsrcavg

         real(kind=c_float), value, intent(in) :: qtsrchgt

         real(kind=c_float), value, intent(in) :: qtsrc_fac

         real(kind=c_float), value, intent(in) :: thlsrc_fac

         real(kind=c_float), value, intent(in) :: frc_rasn

         real(kind=c_float), value, intent(in) :: rbuoy

         real(kind=c_float), value, intent(in) :: epsvarw

         integer(kind=c_int), value, intent(in) :: use_CINcin

         real(kind=c_float), value, intent(in) :: mumin1

         real(kind=c_float), value, intent(in) :: rmaxfrac

         real(kind=c_float), value, intent(in) :: PGFc

         real(kind=c_float), value, intent(in) :: dt

         integer(kind=c_int), value, intent(in) :: niter_xc

         real(kind=c_float), value, intent(in) :: criqc

         real(kind=c_float), value, intent(in) :: rle

         integer(kind=c_int), value, intent(in) :: cridist_opt

         real(kind=c_float), value, intent(in) :: mixscale

         real(kind=c_float), value, intent(in) :: rdrag

         real(kind=c_float), value, intent(in) :: rkm

         integer(kind=c_int), value, intent(in) :: use_self_detrain

         real(kind=c_float), value, intent(in) :: detrhgt

         integer(kind=c_int), value, intent(in) :: use_cumpenent

         real(kind=c_float), value, intent(in) :: rpen

         integer(kind=c_int), value, intent(in) :: use_momenflx

         real(kind=c_float), value, intent(in) :: rdrop

         integer(kind=c_int), value, intent(in) :: iter_cin

         real(kind=c_float), dimension(*), intent(in) :: pifc0_inv
         integer(kind=c_int), dimension(*), intent(in) :: pifc0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: pifc0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: zifc0_inv
         integer(kind=c_int), dimension(*), intent(in) :: zifc0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: zifc0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: pmid0_inv
         integer(kind=c_int), dimension(*), intent(in) :: pmid0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: pmid0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: zmid0_inv
         integer(kind=c_int), dimension(*), intent(in) :: zmid0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: zmid0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: kpbl_inv
         integer(kind=c_int), dimension(*), intent(in) :: kpbl_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: kpbl_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: exnmid0_inv
         integer(kind=c_int), dimension(*), intent(in) :: exnmid0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: exnmid0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: exnifc0_inv
         integer(kind=c_int), dimension(*), intent(in) :: exnifc0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: exnifc0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: dp0_inv
         integer(kind=c_int), dimension(*), intent(in) :: dp0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: dp0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: u0_inv
         integer(kind=c_int), dimension(*), intent(in) :: u0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: u0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: v0_inv
         integer(kind=c_int), dimension(*), intent(in) :: v0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: v0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: qv0_inv
         integer(kind=c_int), dimension(*), intent(in) :: qv0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: qv0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: ql0_inv
         integer(kind=c_int), dimension(*), intent(in) :: ql0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: ql0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: qi0_inv
         integer(kind=c_int), dimension(*), intent(in) :: qi0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: qi0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: t0_inv
         integer(kind=c_int), dimension(*), intent(in) :: t0_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: t0_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: frland_in
         integer(kind=c_int), dimension(*), intent(in) :: frland_in_dim_sizes
         integer(kind=c_int), value, intent(in) :: frland_in_rank

         real(kind=c_float), dimension(*), intent(in) :: tke_inv
         integer(kind=c_int), dimension(*), intent(in) :: tke_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: tke_inv_rank

         real(kind=c_float), dimension(*), intent(in) :: rkfre
         integer(kind=c_int), dimension(*), intent(in) :: rkfre_dim_sizes
         integer(kind=c_int), value, intent(in) :: rkfre_rank

         real(kind=c_float), dimension(*), intent(in) :: cush
         integer(kind=c_int), dimension(*), intent(in) :: cush_dim_sizes
         integer(kind=c_int), value, intent(in) :: cush_rank

         real(kind=c_float), dimension(*), intent(in) :: shfx
         integer(kind=c_int), dimension(*), intent(in) :: shfx_dim_sizes
         integer(kind=c_int), value, intent(in) :: shfx_rank

         real(kind=c_float), dimension(*), intent(in) :: evap
         integer(kind=c_int), dimension(*), intent(in) :: evap_dim_sizes
         integer(kind=c_int), value, intent(in) :: evap_rank

         real(kind=c_float), dimension(*), intent(in) :: cnvtr
         integer(kind=c_int), dimension(*), intent(in) :: cnvtr_dim_sizes
         integer(kind=c_int), value, intent(in) :: cnvtr_rank

         real(kind=c_float), dimension(*), intent(in) :: CNV_Tracers
         integer(kind=c_int), dimension(*), intent(in) :: CNV_Tracers_dim_sizes
         integer(kind=c_int), value, intent(in) :: CNV_Tracers_rank

         real(kind=c_float), dimension(*), intent(out) :: umf_inv
         integer(kind=c_int), dimension(*), intent(in) :: umf_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: umf_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: dcm_inv
         integer(kind=c_int), dimension(*), intent(in) :: dcm_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: dcm_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: qtflx_inv
         integer(kind=c_int), dimension(*), intent(in) :: qtflx_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: qtflx_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: slflx_inv
         integer(kind=c_int), dimension(*), intent(in) :: slflx_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: slflx_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: uflx_inv
         integer(kind=c_int), dimension(*), intent(in) :: uflx_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: uflx_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: vflx_inv
         integer(kind=c_int), dimension(*), intent(in) :: vflx_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: vflx_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: qvten_inv
         integer(kind=c_int), dimension(*), intent(in) :: qvten_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: qvten_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: qlten_inv
         integer(kind=c_int), dimension(*), intent(in) :: qlten_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: qlten_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: qiten_inv
         integer(kind=c_int), dimension(*), intent(in) :: qiten_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: qiten_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: tten_inv
         integer(kind=c_int), dimension(*), intent(in) :: tten_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: tten_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: uten_inv
         integer(kind=c_int), dimension(*), intent(in) :: uten_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: uten_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: vten_inv
         integer(kind=c_int), dimension(*), intent(in) :: vten_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: vten_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: qrten_inv
         integer(kind=c_int), dimension(*), intent(in) :: qrten_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: qrten_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: qsten_inv
         integer(kind=c_int), dimension(*), intent(in) :: qsten_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: qsten_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: cufrc_inv
         integer(kind=c_int), dimension(*), intent(in) :: cufrc_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: cufrc_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: fer_inv
         integer(kind=c_int), dimension(*), intent(in) :: fer_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: fer_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: fdr_inv
         integer(kind=c_int), dimension(*), intent(in) :: fdr_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: fdr_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: ndrop_inv
         integer(kind=c_int), dimension(*), intent(in) :: ndrop_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: ndrop_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: nice_inv
         integer(kind=c_int), dimension(*), intent(in) :: nice_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: nice_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: qldet_inv
         integer(kind=c_int), dimension(*), intent(in) :: qldet_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: qldet_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: qlsub_inv
         integer(kind=c_int), dimension(*), intent(in) :: qlsub_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: qlsub_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: qidet_inv
         integer(kind=c_int), dimension(*), intent(in) :: qidet_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: qidet_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: qisub_inv
         integer(kind=c_int), dimension(*), intent(in) :: qisub_inv_dim_sizes
         integer(kind=c_int), value, intent(in) :: qisub_inv_rank

         real(kind=c_float), dimension(*), intent(out) :: tpert_out
         integer(kind=c_int), dimension(*), intent(in) :: tpert_out_dim_sizes
         integer(kind=c_int), value, intent(in) :: tpert_out_rank

         real(kind=c_float), dimension(*), intent(out) :: qpert_out
         integer(kind=c_int), dimension(*), intent(in) :: qpert_out_dim_sizes
         integer(kind=c_int), value, intent(in) :: qpert_out_rank

      end subroutine compute_uwshcu_f_run_compute_uwshcu

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

   subroutine make_gfdl_1m_flags_C_interop(phys_hydrostatic, hydrostatic, do_qa, fix_negative, fast_sat_adj, &
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
      alin, clin, cld_min, icloud_f, irain_f, gfdl_1m_flags)

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
      real, intent(in) :: alin, clin, cld_min
      integer, intent(in) :: icloud_f, irain_f

      type(gfdl_1m_flags_interface_type), intent(out) :: gfdl_1m_flags

      gfdl_1m_flags%phys_hydrostatic = phys_hydrostatic
      gfdl_1m_flags%hydrostatic = hydrostatic
      gfdl_1m_flags%do_qa = do_qa
      gfdl_1m_flags%fix_negative = fix_negative
      gfdl_1m_flags%fast_sat_adj = fast_sat_adj
      gfdl_1m_flags%const_vi = const_vi
      gfdl_1m_flags%const_vs = const_vs
      gfdl_1m_flags%const_vg = const_vg
      gfdl_1m_flags%const_vr = const_vr
      gfdl_1m_flags%use_ccn = use_ccn
      gfdl_1m_flags%do_bigg = do_bigg
      gfdl_1m_flags%do_evap = do_evap
      gfdl_1m_flags%do_subl = do_subl
      gfdl_1m_flags%z_slope_liq = z_slope_liq
      gfdl_1m_flags%z_slope_ice = z_slope_ice
      gfdl_1m_flags%prog_ccn = prog_ccn
      gfdl_1m_flags%preciprad = preciprad
      gfdl_1m_flags%use_ppm = use_ppm
      gfdl_1m_flags%mono_prof = mono_prof
      gfdl_1m_flags%do_sedi_heat = do_sedi_heat
      gfdl_1m_flags%sedi_transport = sedi_transport
      gfdl_1m_flags%do_sedi_w = do_sedi_w
      gfdl_1m_flags%de_ice = de_ice
      gfdl_1m_flags%mp_print = mp_print
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
      gfdl_1m_flags%rh_inc = rh_inc
      gfdl_1m_flags%rh_ins = rh_ins
      gfdl_1m_flags%rh_inr = rh_inr
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
      gfdl_1m_flags%c_cracw = c_cracw
      gfdl_1m_flags%alin = alin
      gfdl_1m_flags%clin = clin
      gfdl_1m_flags%cld_min = cld_min
      gfdl_1m_flags%icloud_f = icloud_f
      gfdl_1m_flags%irain_f = irain_f

   end subroutine make_gfdl_1m_flags_C_interop

end module pymoist_interface_mod
