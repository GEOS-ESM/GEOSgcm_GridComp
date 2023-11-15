program MicrophysicsDriver

  use gfdl2_cloud_microphys_mod, only: gfdl_cloud_microphys_init, gfdl_cloud_microphys_driver
  use gfdl2_cloud_microphys_orig_mod, only: &
       gfdl_cloud_microphys_orig_init => gfdl_cloud_microphys_init, &
       gfdl_cloud_microphys_orig_driver => gfdl_cloud_microphys_driver

  implicit none

  integer, parameter :: irank = 4
  character(len=*), parameter :: fmt = '(1x, a10, 1x, a1, 1x, a5, 1x, a1, 1x, e15.9)'

  character(len=256) :: file_name
  integer :: file_handle
  real :: start, finish

  ! microphysics variabls
  ! -scalars-
  logical :: hydrostatic, phys_hydrostatic
  integer :: iis, iie, jjs, jje !< physics window
  integer :: kks, kke !< vertical dimension
  integer :: ktop, kbot !< vertical compute domain
  real :: dt_in !< physics time step
  real :: anv_icefall, lsc_icefall
  ! -arrays-
  ! --for-original--
  real, allocatable, dimension (:, :) :: area1, land1, cnv_fraction1, srf_type1, eis1
  real, allocatable, dimension (:, :, :) :: rhcrit1
  real, allocatable, dimension (:, :, :) :: delp1, dz1, uin1, vin1, pt1
  real, allocatable, dimension (:, :, :) :: qv1, ql1, qr1, qg1, qa1, qn1
  real, allocatable, dimension (:, :, :) :: qi1, qs1
  real, allocatable, dimension (:, :, :) :: pt_dt1, qa_dt1, udt1, vdt1, w1
  real, allocatable, dimension (:, :, :) :: qv_dt1, ql_dt1, qr_dt1, qi_dt1, qs_dt1, qg_dt1
  real, allocatable, dimension (:, :) :: rain1, snow1, ice1, graupel1
  real, allocatable, dimension (:, :, :) :: m2_rain1, m2_sol1 ! Rain and Ice fluxes (Pa kg/kg)
  real, allocatable, dimension (:, :, :) :: revap1 ! Rain evaporation
  real, allocatable, dimension (:, :, :) :: isubl1 ! Ice sublimation
  ! --for-new--
  real, allocatable, dimension (:, :) :: area2, land2, cnv_fraction2, srf_type2, eis2
  real, allocatable, dimension (:, :, :) :: rhcrit2
  real, allocatable, dimension (:, :, :) :: delp2, dz2, uin2, vin2, pt2
  real, allocatable, dimension (:, :, :) :: qv2, ql2, qr2, qg2, qa2, qn2
  real, allocatable, dimension (:, :, :) :: qi2, qs2
  real, allocatable, dimension (:, :, :) :: pt_dt2, qa_dt2, udt2, vdt2, w2
  real, allocatable, dimension (:, :, :) :: qv_dt2, ql_dt2, qr_dt2, qi_dt2, qs_dt2, qg_dt2
  real, allocatable, dimension (:, :) :: rain2, snow2, ice2, graupel2
  real, allocatable, dimension (:, :, :) :: m2_rain2, m2_sol2 ! Rain and Ice fluxes (Pa kg/kg)
  real, allocatable, dimension (:, :, :) :: revap2 ! Rain evaporation
  real, allocatable, dimension (:, :, :) :: isubl2 ! Ice sublimation

  ! Read data and call original version
  ! -Read-input-data-
  write(file_name, '(a26, i1, a4)') 'input-data/microphys_data.', irank, '.bin'
  open(newunit = file_handle, file = file_name, form = 'unformatted', status = 'old')
  ! -scalars-
  read(file_handle) &
       iis, iie, jjs, jje, kks, kke, ktop, kbot, dt_in, &
       anv_icefall, lsc_icefall, hydrostatic, phys_hydrostatic
  ! print *, 'iis/iie/jjs/jje/kks/kke/ktop/kbot: ', iis, iie, jjs, jje, kks, kke, ktop, kbot
  ! -arrays-
  ! --first, allocate memory--
  allocate(area1(iis:iie, jjs:jje))
  allocate(land1, cnv_fraction1, srf_type1, eis1, mold=area1)
  allocate(rhcrit1(iis:iie, jjs:jje, kks:kke))
  allocate(&
       delp1, dz1, uin1, vin1, pt1, &
       qv1, ql1, qr1, qg1, qa1, qn1, &
       qi1, qs1, &
       pt_dt1, qa_dt1, udt1, vdt1, w1, &
       qv_dt1, ql_dt1, qr_dt1, qi_dt1, qs_dt1, qg_dt1, &
       mold=rhcrit1)
  allocate(rain1, snow1, ice1, graupel1, mold=area1)
  allocate(m2_rain1, m2_sol1, revap1, isubl1, mold=rhcrit1)
  ! --now, read--
  read(file_handle) &
       ! intent(in)
       area1, land1, cnv_fraction1, srf_type1, eis1, rhcrit1, &
       delp1, dz1, uin1, vin1, pt1, &
       qv1, ql1, qr1, qg1, qa1, qn1, &
       ! intent(inout)
       qi1, qs1, &
       pt_dt1, qa_dt1, udt1, vdt1, w1, &
       qv_dt1, ql_dt1, qr_dt1, qi_dt1, qs_dt1, qg_dt1
  close(file_handle)

  ! ! -print-some-info-about-input-data-
  ! print *, 'dt_in: ', dt_in
  ! print *, 'anv/lsc_icefall: ', anv_icefall, lsc_icefall
  ! print *, 'hydrostatic/phys_hydrostatic: ', hydrostatic, phys_hydrostatic
  ! print *, 'area (shape/min/max/sum): ', shape(area1), minval(area1), maxval(area1), sum(area1)
  ! print *, 'delp (shape/min/max/sum): ', shape(delp1), minval(delp1), maxval(delp1), sum(delp1)
  ! print *, 'qv (shape/min/max/sum): ', shape(qv1), minval(qv1), maxval(qv1), sum(qv1)
  ! print *, 'qi (shape/min/max/sum): ', shape(qi1), minval(qi1), maxval(qi1), sum(qi1)
  ! print *, 'pt_dt (shape/min/max/sum): ', shape(pt_dt1), minval(pt_dt1), maxval(pt_dt1), sum(pt_dt1)
  ! print *, 'qv_dt (shape/min/max/sum): ', shape(qv_dt1), minval(qv_dt1), maxval(qv_dt1), sum(qv_dt1)

  call gfdl_cloud_microphys_init()

  call cpu_time(start)
  call gfdl_cloud_microphys_orig_driver ( &
       ! intent (in)
       qv1, ql1, qr1, &
       ! intent (inout)
       qi1, qs1, &
       ! intent (in)
       qg1, qa1, qn1, &
       ! intent (inout)
       qv_dt1, ql_dt1, qr_dt1, qi_dt1, qs_dt1, qg_dt1, qa_dt1, pt_dt1, &
       ! intent (in)
       pt1, &
       ! intent (inout)
       w1, &
       ! intent (in)
       uin1, vin1, &
       ! intent (inout)
       udt1, vdt1, &
       ! intent (in)
       dz1, delp1, &
       area1, dt_in, land1, cnv_fraction1, srf_type1, eis1, rhcrit1, anv_icefall, lsc_icefall, &
       ! intent (out)
       revap1, isubl1, rain1, snow1, ice1, graupel1, m2_rain1, m2_sol1, &
       ! intent (in)
       hydrostatic, phys_hydrostatic, &
       iis, iie, jjs, jje, kks, kke, ktop, kbot)
  call cpu_time(finish)
  print *, 'Time taken (cpu): ', finish - start, 's'

  ! Read data and call new (gpu) version
  ! -Read-input-data-
  write(file_name, '(a26, i1, a4)') 'input-data/microphys_data.', irank, '.bin'
  open(newunit = file_handle, file = file_name, form = 'unformatted', status = 'old')
  ! -scalars-
  read(file_handle) &
       iis, iie, jjs, jje, kks, kke, ktop, kbot, dt_in, &
       anv_icefall, lsc_icefall, hydrostatic, phys_hydrostatic
  ! print *, 'iis/iie/jjs/jje/kks/kke/ktop/kbot: ', iis, iie, jjs, jje, kks, kke, ktop, kbot
  ! -arrays-
  ! --first, allocate memory--
  allocate(area2(iis:iie, jjs:jje))
  allocate(land2, cnv_fraction2, srf_type2, eis2, mold=area2)
  allocate(rhcrit2(iis:iie, jjs:jje, kks:kke))
  allocate(&
       delp2, dz2, uin2, vin2, pt2, &
       qv2, ql2, qr2, qg2, qa2, qn2, &
       qi2, qs2, &
       pt_dt2, qa_dt2, udt2, vdt2, w2, &
       qv_dt2, ql_dt2, qr_dt2, qi_dt2, qs_dt2, qg_dt2, &
       mold=rhcrit2)
  allocate(rain2, snow2, ice2, graupel2, mold=area2)
  allocate(m2_rain2, m2_sol2, revap2, isubl2, mold=rhcrit2)
  ! --now, read--
  read(file_handle) &
       ! intent(in)
       area2, land2, cnv_fraction2, srf_type2, eis2, rhcrit2, &
       delp2, dz2, uin2, vin2, pt2, &
       qv2, ql2, qr2, qg2, qa2, qn2, &
       ! intent(inout)
       qi2, qs2, &
       pt_dt2, qa_dt2, udt2, vdt2, w2, &
       qv_dt2, ql_dt2, qr_dt2, qi_dt2, qs_dt2, qg_dt2
  close(file_handle)

  ! ! -print-some-info-about-input-data-
  ! print *, 'dt_in: ', dt_in
  ! print *, 'anv/lsc_icefall: ', anv_icefall, lsc_icefall
  ! print *, 'hydrostatic/phys_hydrostatic: ', hydrostatic, phys_hydrostatic
  ! print *, 'area (shape/min/max/sum): ', shape(area2), minval(area2), maxval(area2), sum(area2)
  ! print *, 'delp (shape/min/max/sum): ', shape(delp2), minval(delp2), maxval(delp2), sum(delp2)
  ! print *, 'qv (shape/min/max/sum): ', shape(qv2), minval(qv2), maxval(qv2), sum(qv2)
  ! print *, 'qi (shape/min/max/sum): ', shape(qi2), minval(qi2), maxval(qi2), sum(qi2)
  ! print *, 'pt_dt (shape/min/max/sum): ', shape(pt_dt2), minval(pt_dt2), maxval(pt_dt2), sum(pt_dt2)
  ! print *, 'qv_dt (shape/min/max/sum): ', shape(qv_dt2), minval(qv_dt2), maxval(qv_dt2), sum(qv_dt2)

  call gfdl_cloud_microphys_init()

  call cpu_time(start)
  call gfdl_cloud_microphys_driver ( &
       ! intent (in)
       qv2, ql2, qr2, &
       ! intent (inout)
       qi2, qs2, &
       ! intent (in)
       qg2, qa2, qn2, &
       ! intent (inout)
       qv_dt2, ql_dt2, qr_dt2, qi_dt2, qs_dt2, qg_dt2, qa_dt2, pt_dt2, &
       ! intent (in)
       pt2, &
       ! intent (inout)
       w2, &
       ! intent (in)
       uin2, vin2, &
       ! intent (inout)
       udt2, vdt2, &
       ! intent (in)
       dz2, delp2, &
       area2, dt_in, land2, cnv_fraction2, srf_type2, eis2, rhcrit2, anv_icefall, lsc_icefall, &
       ! intent (out)
       revap2, isubl2, rain2, snow2, ice2, graupel2, m2_rain2, m2_sol2, &
       ! intent (in)
       hydrostatic, phys_hydrostatic, &
       iis, iie, jjs, jje, kks, kke, ktop, kbot)
  call cpu_time(finish)
  print *, 'Time taken (gpu): ', finish - start, 's'

  print *, ''
  print *, '-----------|-------|----------------'
  print *, '       var |  type |     error'
  print *, '-----------|-------|----------------'
  
  write(*, fmt) 'qv', '|','in', '|', norm2(qv2-qv1) ! /norm2(qi1)
  print *, '-----------|-------|----------------'

  write(*, fmt) 'qi', '|','inout', '|', norm2(qi2-qi1) ! /norm2(qi1)
  write(*, fmt) 'qs', '|','inout', '|', norm2(qs2-qs1) ! /norm2(qs1)
  write(*, fmt) 'qv_dt', '|','inout', '|', norm2(qv_dt2-qv_dt1) ! /norm2(qv_dt1)
  write(*, fmt) 'ql_dt', '|','inout', '|', norm2(ql_dt2-ql_dt1) ! /norm2(ql_dt1)
  write(*, fmt) 'qr_dt', '|','inout', '|', norm2(qr_dt2-qr_dt1) ! /norm2(qr_dt1)
  write(*, fmt) 'qi_dt', '|','inout', '|', norm2(qi_dt2-qi_dt1) ! /norm2(qi_dt1)
  write(*, fmt) 'qs_dt', '|','inout', '|', norm2(qs_dt2-qs_dt1) ! /norm2(qs_dt1)
  write(*, fmt) 'qg_dt', '|','inout', '|', norm2(qg_dt2-qg_dt1) ! /norm2(qg_dt1)
  write(*, fmt) 'qa_dt', '|','inout', '|', norm2(qa_dt2-qa_dt1) ! /norm2(qa_dt1)
  write(*, fmt) 'w', '|','inout', '|', norm2(w2-w1) ! /norm2(w1)
  write(*, fmt) 'udt', '|','inout', '|', norm2(udt2-udt1) ! /norm2(udt1)
  write(*, fmt) 'vdt', '|','inout', '|', norm2(vdt2-vdt1) ! /norm2(vdt1)
  print *, '-----------|-------|----------------'

  write(*, fmt) 'revap', '|', 'out', '|', norm2(revap2-revap1) ! /norm2(revap1)
  write(*, fmt) 'isubl', '|', 'out', '|', norm2(isubl2-isubl1) ! /norm2(isubl1)
  write(*, fmt) 'rain', '|', 'out', '|', norm2(rain2-rain1) ! /norm2(rain1)
  write(*, fmt) 'snow', '|', 'out', '|', norm2(snow2-snow1) ! /norm2(snow1)
  write(*, fmt) 'graupel', '|', 'out', '|', norm2(graupel2-graupel1) ! /norm2(graupel1)
  write(*, fmt) 'm2_rain', '|', 'out', '|', norm2(m2_rain2-m2_rain1) ! /norm2(m2_rain1)
  write(*, fmt) 'm2_sol', '|', 'out', '|', norm2(m2_sol2-m2_sol1) ! /norm2(m2_sol1)
  print *, '-----------|-------|----------------'

end program MicrophysicsDriver
