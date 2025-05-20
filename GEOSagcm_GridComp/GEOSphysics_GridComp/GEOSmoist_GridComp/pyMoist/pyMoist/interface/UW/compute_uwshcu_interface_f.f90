module compute_uwshcu_interface_mod

   use iso_c_binding, only: c_int, c_float, c_double, c_bool

   implicit none

   !> Single precision real numbers, 6 digits, range 10⁻³⁷ to 10³⁷-1; 32 bits
   integer, parameter :: sp = selected_real_kind(6, 37)
   !> Double precision real numbers, 15 digits, range 10⁻³⁰⁷ to 10³⁰⁷-1; 64 bits
   integer, parameter :: dp = selected_real_kind(15, 307)

   private
   public :: compute_uwshcu_f_init
   public :: compute_uwshcu_f_run_compute_uwshcu

   interface

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

         integer(kind=c_int), dimension(*), intent(in) :: kpbl_inv
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

end module compute_uwshcu_interface_mod
