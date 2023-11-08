program MicrophysicsDriver

  use gfdl2_cloud_microphys_mod, only: gfdl_cloud_microphys_init, gfdl_cloud_microphys_driver

  implicit none
  
  integer, parameter :: irank = 4

  character(len=256) :: file_name
  integer :: file_handle
  
  ! microphysics variabls
  ! -scalars-
  logical :: hydrostatic, phys_hydrostatic
  integer :: iis, iie, jjs, jje !< physics window
  integer :: kks, kke !< vertical dimension
  integer :: ktop, kbot !< vertical compute domain
  real :: dt_in !< physics time step
  real :: anv_icefall, lsc_icefall
  ! -arrays-
  real, allocatable, dimension (:, :) :: area, land, cnv_fraction, srf_type, eis
  real, allocatable, dimension (:, :, :) :: rhcrit
  real, allocatable, dimension (:, :, :) :: delp, dz, uin, vin, pt
  real, allocatable, dimension (:, :, :) :: qv, ql, qr, qg, qa, qn
  real, allocatable, dimension (:, :, :) :: qi, qs
  real, allocatable, dimension (:, :, :) :: pt_dt, qa_dt, udt, vdt, w
  real, allocatable, dimension (:, :, :) :: qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt
  real, allocatable, dimension (:, :) :: rain, snow, ice, graupel
  real, allocatable, dimension (:, :, :) :: m2_rain, m2_sol ! Rain and Ice fluxes (Pa kg/kg)
  real, allocatable, dimension (:, :, :) :: revap ! Rain evaporation
  real, allocatable, dimension (:, :, :) :: isubl ! Ice sublimation

    ! Read input data
  write(file_name, '(a26, i1, a4)') 'input-data/microphys_data.', irank, '.bin'
  open(newunit = file_handle, file = file_name, form = 'unformatted', status = 'old')
  ! -scalars-
  read(file_handle) &
       iis, iie, jjs, jje, kks, kke, ktop, kbot, dt_in, &
       anv_icefall, lsc_icefall, hydrostatic, phys_hydrostatic
  print *, 'iis/iie/jjs/jje/kks/kke/ktop/kbot: ', iis, iie, jjs, jje, kks, kke, ktop, kbot
  ! -arrays-
  ! --first, allocate memory--
  allocate(area(iis:iie, jjs:jje))
  allocate(land, cnv_fraction, srf_type, eis, mold=area)
  allocate(rhcrit(iis:iie, jjs:jje, kks:kke))
  allocate(&
       delp, dz, uin, vin, pt, &
       qv, ql, qr, qg, qa, qn, &
       qi, qs, &
       pt_dt, qa_dt, udt, vdt, w, &
       qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, &
       mold=rhcrit)
  allocate(rain, snow, ice, graupel, mold=area)
  allocate(m2_rain, m2_sol, revap, isubl, mold=rhcrit)
  ! --now, read--
  read(file_handle) &
       ! intent(in)
       area, land, cnv_fraction, srf_type, eis, rhcrit, &
       delp, dz, uin, vin, pt, &
       qv, ql, qr, qg, qa, qn, &
       ! intent(inout)
       qi, qs, &
       pt_dt, qa_dt, udt, vdt, w, &
       qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt
  close(file_handle)

  ! print some info about input data
  print *, 'dt_in: ', dt_in
  print *, 'anv/lsc_icefall: ', anv_icefall, lsc_icefall
  print *, 'hydrostatic/phys_hydrostatic: ', hydrostatic, phys_hydrostatic
  print *, 'area (shape/min/max/sum): ', shape(area), minval(area), maxval(area), sum(area)
  print *, 'delp (shape/min/max/sum): ', shape(delp), minval(delp), maxval(delp), sum(delp)
  print *, 'qv (shape/min/max/sum): ', shape(qv), minval(qv), maxval(qv), sum(qv)
  print *, 'qi (shape/min/max/sum): ', shape(qi), minval(qi), maxval(qi), sum(qi)
  print *, 'pt_dt (shape/min/max/sum): ', shape(pt_dt), minval(pt_dt), maxval(pt_dt), sum(pt_dt)
  print *, 'qv_dt (shape/min/max/sum): ', shape(qv_dt), minval(qv_dt), maxval(qv_dt), sum(qv_dt)

  call gfdl_cloud_microphys_init()

  call gfdl_cloud_microphys_driver ( &
       ! intent (in)
       qv, ql, qr, &
       ! intent (inout)
       qi, qs, &
       ! intent (in)
       qg, qa, qn, &
       ! intent (inout)
       qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt, qa_dt, pt_dt, &
       ! intent (in)
       pt, &
       ! intent (inout)
       w, &
       ! intent (in)
       uin, vin, &
       ! intent (inout)
       udt, vdt, &
       ! intent (in)
       dz, delp, &
       area, dt_in, land, cnv_fraction, srf_type, eis, rhcrit, anv_icefall, lsc_icefall, &
       ! intent (out)
       revap, isubl, rain, snow, ice, graupel, m2_rain, m2_sol, &
       ! intent (in)
       hydrostatic, phys_hydrostatic, &
       iis, iie, jjs, jje, kks, kke, ktop, kbot)

end program MicrophysicsDriver
