program MicrophysicsDriver

  use input_mod, only: InputScalars_T, InputArrays_T, get_data_from_file
  use gfdl2_cloud_microphys_mod, only: gfdl_cloud_microphys_init, gfdl_cloud_microphys_driver
  use gfdl2_cloud_microphys_orig_mod, only: gfdl_cloud_microphys_cpu_init => gfdl_cloud_microphys_init
  use gfdl2_cloud_microphys_orig_mod, only: gfdl_cloud_microphys_cpu_driver => gfdl_cloud_microphys_driver

  implicit none

  integer, parameter :: irank = 4
  character(len=*), parameter :: fmt = '(1x, a10, 1x, a1, 1x, a5, 1x, a1, 1x, e15.9)'

  character(len=256) :: file_name
  integer :: file_handle
  real :: start, finish

  ! microphysics variabls
  type(InputScalars_T) :: sclr1, sclr2
  type(InputArrays_T) :: arr1, arr2

  ! Input file
  write(file_name, '(a26, i1, a4)') 'input-data/microphys_data.', irank, '.bin'
  
  ! Read data and call original (cpu) version
  call get_data_from_file(file_name, sclr1, arr1)
  print *, 'CPU:'
  call sclr1%write_scalars()
  call arr1%write_arrays()

  call gfdl_cloud_microphys_cpu_init()

  call cpu_time(start)
  call gfdl_cloud_microphys_cpu_driver ( &
       ! intent (in)
       arr1%qv, arr1%ql, arr1%qr, &
       ! intent (inout)
       arr1%qi, arr1%qs, &
       ! intent (in)
       arr1%qg, arr1%qa, arr1%qn, &
       ! intent (inout)
       arr1%qv_dt, arr1%ql_dt, arr1%qr_dt, arr1%qi_dt, arr1%qs_dt, arr1%qg_dt, arr1%qa_dt, arr1%pt_dt, &
       ! intent (in)
       arr1%pt, &
       ! intent (inout)
       arr1%w, &
       ! intent (in)
       arr1%uin, arr1%vin, &
       ! intent (inout)
       arr1%udt, arr1%vdt, &
       ! intent (in)
       arr1%dz, arr1%delp, &
       arr1%area, sclr1%dt_in, arr1%land, arr1%cnv_fraction, arr1%srf_type, arr1%eis, arr1%rhcrit, &
       sclr1%anv_icefall, sclr1%lsc_icefall, &
       ! intent (out)
       arr1%revap, arr1%isubl, arr1%rain, arr1%snow, arr1%ice, arr1%graupel, arr1%m2_rain, arr1%m2_sol, &
       ! intent (in)
       sclr1%hydrostatic, sclr1%phys_hydrostatic, &
       sclr1%iis, sclr1%iie, sclr1%jjs, sclr1%jje, sclr1%kks, sclr1%kke, sclr1%ktop, sclr1%kbot)
  call cpu_time(finish)
  print *, 'Time taken (cpu): ', finish - start, 's'

  ! Read data and call new (gpu) version
  call get_data_from_file(file_name, sclr2, arr2)
  print *, ''
  print *, 'GPU:'
  call sclr2%write_scalars()
  call arr2%write_arrays()

  call gfdl_cloud_microphys_init()

  call cpu_time(start)
  call gfdl_cloud_microphys_driver ( &
       ! intent (in)
       arr2%qv, arr2%ql, arr2%qr, &
       ! intent (inout)
       arr2%qi, arr2%qs, &
       ! intent (in)
       arr2%qg, arr2%qa, arr2%qn, &
       ! intent (inout)
       arr2%qv_dt, arr2%ql_dt, arr2%qr_dt, arr2%qi_dt, arr2%qs_dt, arr2%qg_dt, arr2%qa_dt, arr2%pt_dt, &
       ! intent (in)
       arr2%pt, &
       ! intent (inout)
       arr2%w, &
       ! intent (in)
       arr2%uin, arr2%vin, &
       ! intent (inout)
       arr2%udt, arr2%vdt, &
       ! intent (in)
       arr2%dz, arr2%delp, &
       arr2%area, sclr2%dt_in, arr2%land, arr2%cnv_fraction, arr2%srf_type, arr2%eis, arr2%rhcrit, &
       sclr2%anv_icefall, sclr2%lsc_icefall, &
       ! intent (out)
       arr2%revap, arr2%isubl, arr2%rain, arr2%snow, arr2%ice, arr2%graupel, arr2%m2_rain, arr2%m2_sol, &
       ! intent (in)
       sclr2%hydrostatic, sclr2%phys_hydrostatic, &
       sclr2%iis, sclr2%iie, sclr2%jjs, sclr2%jje, sclr2%kks, sclr2%kke, sclr2%ktop, sclr2%kbot)
  call cpu_time(finish)
  print *, 'Time taken (gpu): ', finish - start, 's'

  print *, ''
  print *, '-----------|-------|----------------'
  print *, '       var |  type |     error'
  print *, '-----------|-------|----------------'
  
  write(*, fmt) 'qv', '|','in', '|', norm2(arr2%qv-arr1%qv) ! /norm2(arr1%qv)
  print *, '-----------|-------|----------------'

  write(*, fmt) 'qi', '|','inout', '|', norm2(arr2%qi-arr1%qi) ! /norm2(arr1%qi)
  write(*, fmt) 'qs', '|','inout', '|', norm2(arr2%qs-arr1%qs) ! /norm2(arr1%qs)
  write(*, fmt) 'qv_dt', '|','inout', '|', norm2(arr2%qv_dt-arr1%qv_dt) ! /norm2(arr1%qv_dt)
  write(*, fmt) 'ql_dt', '|','inout', '|', norm2(arr2%ql_dt-arr1%ql_dt) ! /norm2(arr1%ql_dt)
  write(*, fmt) 'qr_dt', '|','inout', '|', norm2(arr2%qr_dt-arr1%qr_dt) ! /norm2(arr1%qr_dt)
  write(*, fmt) 'qi_dt', '|','inout', '|', norm2(arr2%qi_dt-arr1%qi_dt) ! /norm2(arr1%qi_dt)
  write(*, fmt) 'qs_dt', '|','inout', '|', norm2(arr2%qs_dt-arr1%qs_dt) ! /norm2(arr1%qs_dt)
  write(*, fmt) 'qg_dt', '|','inout', '|', norm2(arr2%qg_dt-arr1%qg_dt) ! /norm2(arr1%qg_dt)
  write(*, fmt) 'qa_dt', '|','inout', '|', norm2(arr2%qa_dt-arr1%qa_dt) ! /norm2(arr1%qa_dt)
  write(*, fmt) 'w', '|','inout', '|', norm2(arr2%w-arr1%w) ! /norm2(arr1%w)
  write(*, fmt) 'udt', '|','inout', '|', norm2(arr2%udt-arr1%udt) ! /norm2(arr1%udt)
  write(*, fmt) 'vdt', '|','inout', '|', norm2(arr2%vdt-arr1%vdt) ! /norm2(arr1%vdt)
  print *, '-----------|-------|----------------'

  write(*, fmt) 'revap', '|', 'out', '|', norm2(arr2%revap-arr1%revap) ! /norm2(arr1%revap)
  write(*, fmt) 'isubl', '|', 'out', '|', norm2(arr2%isubl-arr1%isubl) ! /norm2(arr1%isubl)
  write(*, fmt) 'rain', '|', 'out', '|', norm2(arr2%rain-arr1%rain) ! /norm2(arr1%rain)
  write(*, fmt) 'snow', '|', 'out', '|', norm2(arr2%snow-arr1%snow) ! /norm2(arr1%snow)
  write(*, fmt) 'graupel', '|', 'out', '|', norm2(arr2%graupel-arr1%graupel) ! /norm2(arr1%graupel)
  write(*, fmt) 'm2_rain', '|', 'out', '|', norm2(arr2%m2_rain-arr1%m2_rain) ! /norm2(arr1%m2_rain)
  write(*, fmt) 'm2_sol', '|', 'out', '|', norm2(arr2%m2_sol-arr1%m2_sol) ! /norm2(arr1%m2_sol)
  print *, '-----------|-------|----------------'

end program MicrophysicsDriver
