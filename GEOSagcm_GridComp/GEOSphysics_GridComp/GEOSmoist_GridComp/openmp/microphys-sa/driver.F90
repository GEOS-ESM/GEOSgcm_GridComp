program MicrophysicsDriver

  use input_mod, only: InputScalars_T, InputArrays_T, get_data_from_file
  use output_mod, only: OutputArrays_T, write_difference
  use gfdl2_cloud_microphys_mod, only: gfdl_cloud_microphys_init, gfdl_cloud_microphys_driver
  use gfdl2_cloud_microphys_cpu_mod, only: gfdl_cloud_microphys_cpu_init => gfdl_cloud_microphys_init
  use gfdl2_cloud_microphys_cpu_mod, only: gfdl_cloud_microphys_cpu_driver => gfdl_cloud_microphys_driver

  implicit none

  integer, parameter :: irank = 4
  character(len=*), parameter :: fmt = '(1x, a10, 1x, a1, 1x, a5, 1x, a1, 1x, e15.9)'

  character(len=256) :: file_name
  integer :: file_handle
  real :: start, finish

  ! microphysics variabls
  type(InputScalars_T) :: sclr1, sclr2
  type(InputArrays_T) :: inarr1, inarr2
  type(OutputArrays_T) :: outarr1, outarr2

  ! Input file
  write(file_name, '(a26, i1, a4)') 'input-data/microphys_data.', irank, '.bin'
  
  ! Read data and call original (cpu) version
  call get_data_from_file(file_name, sclr1, inarr1)
  print *, 'CPU:'
  call sclr1%write_scalars()
  call inarr1%write_arrays()
  outarr1 = OutputArrays_T(sclr1%iis, sclr1%iie, sclr1%jjs, sclr1%jje, sclr1%kks, sclr1%kke)

  call gfdl_cloud_microphys_cpu_init()

  call cpu_time(start)
  call gfdl_cloud_microphys_cpu_driver ( &
       ! intent (in)
       inarr1%qv, inarr1%ql, inarr1%qr, &
       ! intent (inout)
       inarr1%qi, inarr1%qs, &
       ! intent (in)
       inarr1%qg, inarr1%qa, inarr1%qn, &
       ! intent (inout)
       inarr1%qv_dt, inarr1%ql_dt, inarr1%qr_dt, inarr1%qi_dt, &
       inarr1%qs_dt, inarr1%qg_dt, inarr1%qa_dt, inarr1%pt_dt, &
       ! intent (in)
       inarr1%pt, &
       ! intent (inout)
       inarr1%w, &
       ! intent (in)
       inarr1%uin, inarr1%vin, &
       ! intent (inout)
       inarr1%udt, inarr1%vdt, &
       ! intent (in)
       inarr1%dz, inarr1%delp, inarr1%area, sclr1%dt_in, inarr1%land, inarr1%cnv_fraction, &
       inarr1%srf_type, inarr1%eis, inarr1%rhcrit, sclr1%anv_icefall, sclr1%lsc_icefall, &
       ! intent (out)
       outarr1%revap, outarr1%isubl, &
       outarr1%rain, outarr1%snow, outarr1%ice, outarr1%graupel, &
       outarr1%m2_rain, outarr1%m2_sol, &
       ! intent (in)
       sclr1%hydrostatic, sclr1%phys_hydrostatic, &
       sclr1%iis, sclr1%iie, sclr1%jjs, sclr1%jje, sclr1%kks, sclr1%kke, sclr1%ktop, sclr1%kbot)
  call cpu_time(finish)
  print *, 'Time taken (cpu): ', finish - start, 's'

  ! Read data and call new (gpu) version
  call get_data_from_file(file_name, sclr2, inarr2)
  print *, ''
  print *, 'GPU:'
  call sclr2%write_scalars()
  call inarr2%write_arrays()
  outarr2 = OutputArrays_T(sclr2%iis, sclr2%iie, sclr2%jjs, sclr2%jje, sclr2%kks, sclr2%kke)

  call gfdl_cloud_microphys_init()

  call cpu_time(start)
  call gfdl_cloud_microphys_driver ( &
       ! intent (in)
       inarr2%qv, inarr2%ql, inarr2%qr, &
       ! intent (inout)
       inarr2%qi, inarr2%qs, &
       ! intent (in)
       inarr2%qg, inarr2%qa, inarr2%qn, &
       ! intent (inout)
       inarr2%qv_dt, inarr2%ql_dt, inarr2%qr_dt, inarr2%qi_dt, &
       inarr2%qs_dt, inarr2%qg_dt, inarr2%qa_dt, inarr2%pt_dt, &
       ! intent (in)
       inarr2%pt, &
       ! intent (inout)
       inarr2%w, &
       ! intent (in)
       inarr2%uin, inarr2%vin, &
       ! intent (inout)
       inarr2%udt, inarr2%vdt, &
       ! intent (in)
       inarr2%dz, inarr2%delp, inarr2%area, sclr2%dt_in, inarr2%land, inarr2%cnv_fraction, &
       inarr2%srf_type, inarr2%eis, inarr2%rhcrit, sclr2%anv_icefall, sclr2%lsc_icefall, &
       ! intent (out)
       outarr2%revap, outarr2%isubl, &
       outarr2%rain, outarr2%snow, outarr2%ice, outarr2%graupel, &
       outarr2%m2_rain, outarr2%m2_sol, &
       ! intent (in)
       sclr2%hydrostatic, sclr2%phys_hydrostatic, &
       sclr2%iis, sclr2%iie, sclr2%jjs, sclr2%jje, sclr2%kks, sclr2%kke, sclr2%ktop, sclr2%kbot)
  call cpu_time(finish)
  print *, 'Time taken (gpu): ', finish - start, 's'

  print *, ''
  print *, '-----------|-------|----------------'
  print *, '       var |  type |     error'
  print *, '-----------|-------|----------------'
  
  write(*, fmt) 'qv', '|','in', '|', norm2(inarr2%qv-inarr1%qv) ! /norm2(inarr1%qv)
  print *, '-----------|-------|----------------'

  write(*, fmt) 'qi', '|','inout', '|', norm2(inarr2%qi-inarr1%qi) ! /norm2(inarr1%qi)
  write(*, fmt) 'qs', '|','inout', '|', norm2(inarr2%qs-inarr1%qs) ! /norm2(inarr1%qs)
  write(*, fmt) 'qv_dt', '|','inout', '|', norm2(inarr2%qv_dt-inarr1%qv_dt) ! /norm2(inarr1%qv_dt)
  write(*, fmt) 'ql_dt', '|','inout', '|', norm2(inarr2%ql_dt-inarr1%ql_dt) ! /norm2(inarr1%ql_dt)
  write(*, fmt) 'qr_dt', '|','inout', '|', norm2(inarr2%qr_dt-inarr1%qr_dt) ! /norm2(inarr1%qr_dt)
  write(*, fmt) 'qi_dt', '|','inout', '|', norm2(inarr2%qi_dt-inarr1%qi_dt) ! /norm2(inarr1%qi_dt)
  write(*, fmt) 'qs_dt', '|','inout', '|', norm2(inarr2%qs_dt-inarr1%qs_dt) ! /norm2(inarr1%qs_dt)
  write(*, fmt) 'qg_dt', '|','inout', '|', norm2(inarr2%qg_dt-inarr1%qg_dt) ! /norm2(inarr1%qg_dt)
  write(*, fmt) 'qa_dt', '|','inout', '|', norm2(inarr2%qa_dt-inarr1%qa_dt) ! /norm2(inarr1%qa_dt)
  write(*, fmt) 'w', '|','inout', '|', norm2(inarr2%w-inarr1%w) ! /norm2(inarr1%w)
  write(*, fmt) 'udt', '|','inout', '|', norm2(inarr2%udt-inarr1%udt) ! /norm2(inarr1%udt)
  write(*, fmt) 'vdt', '|','inout', '|', norm2(inarr2%vdt-inarr1%vdt) ! /norm2(inarr1%vdt)
  print *, '-----------|-------|----------------'

  ! write(*, fmt) 'revap', '|', 'out', '|', norm2(inarr2%revap-inarr1%revap) ! /norm2(inarr1%revap)
  ! write(*, fmt) 'isubl', '|', 'out', '|', norm2(inarr2%isubl-inarr1%isubl) ! /norm2(inarr1%isubl)
  ! write(*, fmt) 'rain', '|', 'out', '|', norm2(inarr2%rain-inarr1%rain) ! /norm2(inarr1%rain)
  ! write(*, fmt) 'snow', '|', 'out', '|', norm2(inarr2%snow-inarr1%snow) ! /norm2(inarr1%snow)
  ! write(*, fmt) 'graupel', '|', 'out', '|', norm2(inarr2%graupel-inarr1%graupel) ! /norm2(inarr1%graupel)
  ! write(*, fmt) 'm2_rain', '|', 'out', '|', norm2(inarr2%m2_rain-inarr1%m2_rain) ! /norm2(inarr1%m2_rain)
  ! write(*, fmt) 'm2_sol', '|', 'out', '|', norm2(inarr2%m2_sol-inarr1%m2_sol) ! /norm2(inarr1%m2_sol)
  ! print *, '-----------|-------|----------------'
  call write_difference(outarr1, outarr2)
  
end program MicrophysicsDriver
