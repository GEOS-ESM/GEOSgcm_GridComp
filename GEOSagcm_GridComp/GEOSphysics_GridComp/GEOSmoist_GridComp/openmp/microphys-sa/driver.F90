program MicrophysicsDriver

  use mpi

  implicit none

  integer :: irank, mpi_err

  call MPI_Init(mpi_err)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, mpi_err)

  call serial_driver(irank)

  call MPI_Finalize(mpi_err)

contains

  subroutine serial_driver(irank)

    use input_mod, only: InputScalars_T, InputArrays_T, get_data_from_file
    use input_mod, only: write_inout_difference => write_difference
    use output_mod, only: OutputArrays_T, write_output_difference => write_difference
    use gfdl2_cloud_microphys_mod, only: gfdl_cloud_microphys_init, gfdl_cloud_microphys_driver
    use gfdl2_cloud_microphys_cpu_mod, only: gfdl_cloud_microphys_cpu_init => gfdl_cloud_microphys_init
    use gfdl2_cloud_microphys_cpu_mod, only: gfdl_cloud_microphys_cpu_driver => gfdl_cloud_microphys_driver

    implicit none

    ! Arguments
    integer, intent(in):: irank

    ! Locals
    character(len=*), parameter :: fmt = '(1x, a1, i2, a1, 1x, a, f10.7, 1x, a1)'
    character(len=256) :: file_name
    integer :: file_handle, nranks, i
    real :: start, finish, cpu_time_, gpu_time_
    type(InputScalars_T) :: sclr1, sclr2
    type(InputArrays_T) :: inarr1, inarr2
    type(OutputArrays_T) :: outarr1, outarr2

    ! Input file
    write(file_name, '(a26, i2.2, a4)') 'input-data/microphys_data.', irank, '.bin'

    ! Read data and call original (cpu) version
    call get_data_from_file(file_name, sclr1, inarr1)
    ! print *, 'CPU:'
    ! call sclr1%write_scalars()
    ! call inarr1%write_arrays()
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
    cpu_time_ = finish - start
    ! call outarr1%write_arrays()

    ! Read data and call new (gpu) version
    call get_data_from_file(file_name, sclr2, inarr2)
    ! print *, ''
    ! print *, 'GPU:'
    ! call sclr2%write_scalars()
    ! call inarr2%write_arrays()
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
    gpu_time_ = finish - start
    ! call outarr2%write_arrays()

    call MPI_Comm_size(MPI_COMM_WORLD, nranks, mpi_err)
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
    do i = 0, nranks-1
       if (i == irank) then
          write(*, *)
          write(*, fmt) '[', i, ']', 'Time taken (cpu):', cpu_time_, 's'
          write(*, fmt) '[', i, ']', 'Time taken (gpu):', gpu_time_, 's'
          ! Write errors to stdout
          call write_inout_difference(inarr1, inarr2)
          call write_output_difference(outarr1, outarr2)
       end if
       call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
    end do

  end subroutine serial_driver

end program MicrophysicsDriver
