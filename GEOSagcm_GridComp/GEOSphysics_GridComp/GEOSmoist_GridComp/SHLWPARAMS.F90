      module SHLWPARAMS

      implicit none

! Define SHLWPARAM
! ---------------
      integer use_tracer_transp_uw   ! transport tracers in UW

      type SHLWPARAM_TYPE
           integer  :: niter_xc           ! Number xc iterations
           integer  :: iter_cin           ! Number iterations for implicit CIN
           integer  :: use_CINcin         ! if true, calc CIN thru L..
           integer  :: use_self_detrain   ! 
           integer  :: use_momenflx       ! Perform momentum transport
           integer  :: use_cumpenent      ! Cumulative penetrative entrainment
           integer  :: scverbose          ! activate print statements
           integer  :: windsrcavg         ! Source air uses PBL mean momentum
           real     :: rpen               ! Penentrative entrainment factor
           real     :: rle
           real     :: rkm                ! Factor controlling lateral mixing rate
           real     :: mixscale           ! Controls vertical structure of mixing
           real     :: detrhgt            ! Mixing rate increases above this height
           real     :: rkfre              ! Vertical velocity fraction of tke
           real     :: rmaxfrac           ! Maximum core updraft fraction
           real     :: mumin1             ! 
           real     :: rbuoy              ! Non-hydro pressure effect on updraft
           real     :: rdrag              ! Drag coefficient
           real     :: epsvarw            ! Variance of PBL w by mesoscale
           real     :: PGFc               ! Pressure gradient force
           real     :: criqc              ! Updraft maximum condensate 
           real     :: frc_rasn           ! Precip fraction of expelled condensate
           real     :: kevp               ! Evaporative efficiency
           real     :: rdrop              ! liquid drop radius
           real     :: thlsrc_fac         ! Scaling factor for thlsrc perturbation
           real     :: qtsrc_fac          ! Scaling factor for qtsrc perturbation
      endtype SHLWPARAM_TYPE

  end module SHLWPARAMS
