MODULE ConvParGFCore
  USE PlumeStateModule
USE PlumeDiagnosticsModule
  IMPLICIT NONE

CONTAINS

  ! Main convection driver
  SUBROUTINE GF2020_DRV(mxp,myp,mzp,mtp,nmp         &
       ,ims,ime, jms,jme, kms,kme              &
       ,its,ite, jts,jte, kts,kte              &
       ,flip                                   &
       ,mynum                 &
       ,dt                    &
       ,dx2d                  &
       ,stochastic_sig        &
       ,zm                    &
       ,zt                    &
       ,dm                    &
       ,lons                  &
       ,lats                  &
       ,aot500                &
       ,temp2m                &
       ,sflux_r               &
       ,sflux_t               &
       ,topt                  &
       ,xland                 &
       ,sfc_press             &
       ,kpbl                  &
       ,cnvfrc                &
       ,srftype               &
       ,col_sat               &
       ,u                     &
       ,v                     &
       ,om                    &
       ,temp                  &
       ,press                 &
       ,rvap                  &
       ,mp_ice                &
       ,mp_liq                &
       ,mp_cf                 &
       ,ccn_in                &
       ,curr_rvap             &
       !---- forcings---
       ,buoy_exc              &
       ,rthften               &! gsf_t
       ,rqvften               &! gsf_q
       ,rth_advten            &!advf_t
       ,rthblten              &!sgsf_t
       ,rqvblten              &!sgsf_q
       !---- output ----
       ,conprr                &
       ,lightn_dens           &
       ,rthcuten              &
       ,rqvcuten              &
       ,rqccuten              &
       ,rucuten               &
       ,rvcuten               &
       ,sub_mpqi              &
       ,sub_mpql              &
       ,sub_mpcf              &
       ,rbuoycuten            &
       ,rchemcuten            &
       ,revsu_gf              &
       ,prfil_gf              &
       ,do_this_column        &
       ,plume_states          &
       ,AA0,AA1,AA2,AA3,AA1_BL,AA1_CIN,TAU_BL,TAU_DP,TAU_MD)
    ! ... subroutine body with all code from ConvPar_GF2020.F90 (as previously mapped) ...
  END SUBROUTINE GF2020_DRV

  ! Convection core, implementation of CUP_gf
  SUBROUTINE CUP_gf(...)
    ! ... All CUP_gf core logic migrated here.
  END SUBROUTINE CUP_gf

  ! Any internal helpers required by GF2020_DRV and CUP_gf
  ! (gf2020_drv_init_slice_arrays, gf2020_drv_load_surface_data, etc.)
  ! Place here as CONTAINS

END MODULE ConvParGFCore
