
module catch_iau

  ! module for "incremental analysis update" of Catchment in tile-space
  
  ! reichle+csdraper, 3 Apr 2012

  use catch_constants, ONLY:                 &
       N_SNOW        => CATCH_N_SNOW,        &
       N_GT          => CATCH_N_GT     

  use lsm_routines, ONLY: catch_calc_soil_moist
  
  implicit none
  
  private
  
  public :: apply_catch_iau, check_catch_progn
  
contains
  
  ! ***********************************************************************
  
  subroutine apply_catch_iau( NTILES,                                      & 
       VEG, DZSF, VGWMAX, CDCR1, CDCR2, PSIS, BEE, POROS, WPWET,           & 
       ARS1, ARS2, ARS3, ARA1, ARA2, ARA3, ARA4, ARW1, ARW2, ARW3, ARW4,   & 
       TC1_INC, TC2_INC, TC4_INC, QC1_INC, QC2_INC, QC4_INC,               & 
       CAPAC_INC, CATDEF_INC, RZEXC_INC, SRFEXC_INC,                       & 
       GHTCNT_INC, WESNN_INC, HTSNNN_INC, SNDZN_INC,                       &
       TC1, TC2, TC4, QC1, QC2, QC4,                                       & 
       CAPAC, CATDEF, RZEXC, SRFEXC, 	                                   & 
       GHTCNT, WESNN, HTSNNN, SNDZN  )
    
    implicit none
    
    integer,                           intent(in)    :: NTILES
    
    ! CATCHMENT MODEL PARAMETERS 
    
    integer, dimension(       NTILES), intent(in)    :: VEG
    real,    dimension(       NTILES), intent(in)    :: DZSF, VGWMAX, CDCR1, CDCR2
    real,    dimension(       NTILES), intent(in)    :: PSIS, BEE, POROS, WPWET
    real,    dimension(       NTILES), intent(in)    :: ARS1, ARS2, ARS3
    real,    dimension(       NTILES), intent(in)    :: ARA1, ARA2, ARA3, ARA4
    real,    dimension(       NTILES), intent(in)    :: ARW1, ARW2, ARW3, ARW4

    ! CATCHMENT MODEL PROGNOSTIC INCREMENTS
    
    real,    dimension(       NTILES), intent(in)    :: TC1_INC, TC2_INC, TC4_INC
    real,    dimension(       NTILES), intent(in)    :: QC1_INC, QC2_INC, QC4_INC
    real,    dimension(       NTILES), intent(in)    :: CAPAC_INC, CATDEF_INC
    real,    dimension(       NTILES), intent(in)    :: RZEXC_INC, SRFEXC_INC
    real,    dimension(N_GT,  NTILES), intent(in)    :: GHTCNT_INC
    real,    dimension(N_SNOW,NTILES), intent(in)    :: WESNN_INC, HTSNNN_INC, SNDZN_INC
    
    ! CATCHMENT MODEL PROGNOSTICS
    
    real,    dimension(       NTILES), intent(inout) :: TC1, TC2, TC4, QC1, QC2, QC4
    real,    dimension(       NTILES), intent(inout) :: CAPAC, CATDEF, RZEXC, SRFEXC
    real,    dimension(N_GT,  NTILES), intent(inout) :: GHTCNT
    real,    dimension(N_SNOW,NTILES), intent(inout) :: WESNN, HTSNNN, SNDZN

    ! ---------------------------------------------------------------------------

    ! update prognostic variables 
    
    TC1	    = TC1    + TC1_INC
    TC2	    = TC2    + TC2_INC
    TC4	    = TC4    + TC4_INC

    QC1	    = QC1    + QC1_INC
    QC2	    = QC2    + QC2_INC
    QC4	    = QC4    + QC4_INC

    CAPAC   = CAPAC  + CAPAC_INC

    CATDEF  = CATDEF + CATDEF_INC
    RZEXC   = RZEXC  + RZEXC_INC
    SRFEXC  = SRFEXC + SRFEXC_INC
    
    GHTCNT  = GHTCNT + GHTCNT_INC
    
    WESNN   = WESNN  + WESNN_INC
    HTSNNN  = HTSNNN + HTSNNN_INC
    SNDZN   = SNDZN  + SNDZN_INC
    
    ! make sure that updated prognostics are OK
    
    call check_catch_progn( NTILES,                                          &
         VEG, DZSF, VGWMAX, CDCR1, CDCR2, PSIS, BEE, POROS, WPWET,           & 
         ARS1, ARS2, ARS3, ARA1, ARA2, ARA3, ARA4, ARW1, ARW2, ARW3, ARW4,   & 
         TC1, TC2, TC4, QC1, QC2, QC4,                                       & 
         CAPAC, CATDEF, RZEXC, SRFEXC, 	                                     &  
         GHTCNT, WESNN, HTSNNN, SNDZN  )
    
  end subroutine apply_catch_iau
  
  ! ***********************************************************************
  
  subroutine check_catch_progn( NTILES,                                    &
       VEG, DZSF, VGWMAX, CDCR1, CDCR2, PSIS, BEE, POROS, WPWET,           & 
       ARS1, ARS2, ARS3, ARA1, ARA2, ARA3, ARA4, ARW1, ARW2, ARW3, ARW4,   & 
       TC1, TC2, TC4, QC1, QC2, QC4,                                       & 
       CAPAC, CATDEF, RZEXC, SRFEXC,                                       & 
       GHTCNT, WESNN, HTSNNN, SNDZN  )
    
    ! check Catchment prognostic variables for physical constraints, re-set
    !  if constraints are violated
    !
    ! reichle,  2 Aug 2005
    ! reichle,  5 Feb 2008 - moved from clsm_ensdrv_pert_routines.F90 and 
    !                        added workaround "catdef>1."
    ! reichle,  3 Apr 2012 - moved to here from file "clsm_ensdrv_drv_routines.F90" 
    !                         of "lana" directory (offline LDAS) and revised for use
    !                         without off-line "catch_types" structures
    !                        NOTE: in LDASsa use wrapper check_cat_progn() in file
    !                              clsm_ensdrv_drv_routines.F90
    !
    ! -------------------------------------------------------------------
    
    implicit none

    integer,                           intent(in)    :: NTILES
    
    ! CATCHMENT MODEL PARAMETERS 
    
    integer, dimension(       NTILES), intent(in)    :: VEG
    real,    dimension(       NTILES), intent(in)    :: DZSF, VGWMAX, CDCR1, CDCR2
    real,    dimension(       NTILES), intent(in)    :: PSIS, BEE, POROS, WPWET
    real,    dimension(       NTILES), intent(in)    :: ARS1, ARS2, ARS3
    real,    dimension(       NTILES), intent(in)    :: ARA1, ARA2, ARA3, ARA4
    real,    dimension(       NTILES), intent(in)    :: ARW1, ARW2, ARW3, ARW4

    ! CATCHMENT MODEL PROGNOSTICS
    
    real,    dimension(       NTILES), intent(inout) :: TC1, TC2, TC4, QC1, QC2, QC4
    real,    dimension(       NTILES), intent(inout) :: CAPAC, CATDEF, RZEXC, SRFEXC
    real,    dimension(N_GT,  NTILES), intent(inout) :: GHTCNT
    real,    dimension(N_SNOW,NTILES), intent(inout) :: WESNN, HTSNNN, SNDZN

    ! ----------------------------------------------------------------
    
    ! local variables
    
    integer :: i, k 
    
    real, dimension(NTILES) :: ar1, ar2, ar4

    ! ----------------------------------------------------------------
    
    ! check for violations of physical constraints and correct accordingly
    
    do i=1,NTILES
       
       ! tc1,tc2,tc4 - no checks implemented
       
       ! enforce qc>=0,      maybe qc <= some number ? 
       
       qc1(i)    = max( qc1(i),   0.)  
       qc2(i)    = max( qc2(i),   0.)
       qc4(i)    = max( qc4(i),   0.)
       
       ! enforce capac>=0,   maybe capac <= satcap ?
       
       capac(i)  = max( capac(i), 0.)
       
       ! checks on soil moisture states see below!! (call to calc_soil_moist())
       
       ! no checks on ground heat content implemented
       !
       ! ghtcnt(1:N_gt,i)
       
       do k=1,N_snow
          
          ! snow water equivalent >= 0
          
          wesnn(k,i) = max(wesnn(k,i), 0.)
          
          ! snow heat content <= 0 ???
          
          !! htsnn(k,i) = min(htsnn(k,i), 0.)
          
          ! snow depth >= 0
          
          sndzn(k,i)  = max(sndzn(k,i), 0.)
          
       end do
       
    end do

    ! check soil moisture states (done as part of calculation of
    ! soil moisture content)
    ! reichle, 6 Feb 2004

    ! NOTE: calc_soil_moist() was moved into catchment.F90 (when GEOS5
    !       went from catchment.f to catchment.F90).  The constraint
    !       in calc_soil_moist() was originally catdef>0., but this
    !       proved insufficient when the code was compiled on discover
    !       with "-openmp" because of unprotected divisions by zero
    !       in partition().  (See comment dated 26 March 2007 in the old
    !       catchment.f)
    !       Here, preface the call to calc_soil_moist() with the appropriate
    !       lower bound so that the off-line driver can be used with older 
    !       versions of catchment.F90 (at least version 1.37 and earlier).
    !       IMPORTANT: This *will* mess up the optional diagnostic
    !       "werror", which is not used here but may be in the future.
    !       reichle - 5 Feb 2008
    
    ! call to revised subroutine catch_calc_soil_moist() -- which includes the
    ! lower bound on catdef, - reichle, 3 Apr 2012
        
    call catch_calc_soil_moist( &
         NTILES,veg,dzsf,vgwmax,cdcr1,cdcr2,psis,bee,poros,wpwet, &
         ars1,ars2,ars3,ara1,ara2, &
         ara3,ara4,arw1,arw2,arw3,arw4, &
         srfexc,rzexc,catdef, &
         ar1, ar2, ar4 )
    
  end subroutine check_catch_progn

end module catch_iau


! ==================== EOF ================================================
