program audit_step4_vs_statsgo
  ! to compile: ifort -O2 -g audit_step4_vs_statsgo.f90 -o audit_step4_vs_statsgo.x
  ! to run : ./audit_step4_vs_statsgo.x
  ! summaries of how step4 compares to statsgo used to fill
  implicit none
  integer, parameter :: nc = 43200
  integer, parameter :: nr = 21600

  character(len=*), parameter :: dir_out = &
  '/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/HWSDv2-NGDC-STATSGO/'

  character(len=*), parameter :: dir_stat = &
  '/discover/nobackup/projects/gmao/bcs_shared/preprocessing_bcs_inputs/land/soil/soil_properties/v2/STATSGO2/statsgo_'

  character(len=*), parameter :: fmask = 'data_sources.msk'

  call audit_var('clay_top', 'clay_top_30sec.dat')
  call audit_var('sand_top', 'sand_top_30sec.dat')
  call audit_var('clay_sub', 'clay_sub_30sec.dat')
  call audit_var('sand_sub', 'sand_sub_30sec.dat')
  call audit_var('oc_top',   'oc_top_30sec.dat')
  call audit_var('oc_sub',   'oc_sub_30sec.dat')

contains

  logical function is_statsgo_used(maskval)
    implicit none
    integer(4), intent(in) :: maskval
    integer(4) :: m
    m = mod(maskval, 100_4)
    is_statsgo_used = (m == 22_4 .or. m == 82_4 .or. m == 28_4)
  end function is_statsgo_used

  subroutine audit_var(label, varfile)
    implicit none
    character(len=*), intent(in) :: label, varfile

    integer :: i, j
    integer :: u_mask, u_out, u_stat
    integer(4), allocatable :: mask_row(:)
    real(4),    allocatable :: out_row(:), stat_row(:)
    real(4) :: max_abs, absd
    integer :: nstatrow
    integer(8) :: n_used, n_used_stat_valid
    integer(8) :: n_mis_all, n_mis_stat_valid
    real(8) :: hole_frac

    allocate(mask_row(nc))
    allocate(out_row(nc))
    allocate(stat_row(nc))

    u_mask = 10
    u_out  = 11
    u_stat = 12

    open(u_mask, file=trim(dir_out)//trim(fmask),    form='unformatted', status='old', action='read')
    open(u_out,  file=trim(dir_out)//trim(varfile),  form='unformatted', status='old', action='read')

    ! STATSGO: F77 sequential with 4-byte markers per row; likely endian mismatch
    open(u_stat, file=trim(dir_stat)//trim(varfile), form='unformatted', status='old', action='read', convert='little_endian')

    max_abs   = 0.0_4
    n_used            = 0_8
    n_used_stat_valid = 0_8
    n_mis_all         = 0_8
    n_mis_stat_valid  = 0_8    

    do j = 1, nr
      read(u_mask) (mask_row(i), i=1,nc)
      read(u_out)  (out_row(i),  i=1,nc)
      read(u_stat) (stat_row(i), i=1,nc)

!      if (j == 1 .or. j == nr/2 .or. j == nr) then
      if (j == 15600 .or. j == 15000 .or. j == 16200) then
        print *, '--- ', trim(label), ' row j=', j
        print *, 'STAT(r4) min/max=', minval(stat_row), maxval(stat_row), &
                 ' neg=', count(stat_row < 0.0_4), ' gt100=', count(stat_row > 100.0_4)
        print *, 'OUT (r4) min/max=', minval(out_row),  maxval(out_row),  &
                 ' neg=', count(out_row  < 0.0_4), ' gt100=', count(out_row  > 100.0_4)
        print *, 'STAT samples:', stat_row(1), stat_row(1000), stat_row(20000)
        print *, 'OUT  samples:', out_row(1),  out_row(1000),  out_row(20000)
      end if
      nstatrow = 0
      do i = 1, nc
        if (is_statsgo_used(mask_row(i))) nstatrow = nstatrow + 1
      end do

      if (nstatrow > 0 .and. (j == 1 .or. mod(j,2000)==0)) then
        print *, 'ROW j=', j, ' statsgo_used_pixels_in_row=', nstatrow, &
                 ' STAT min/max=', minval(stat_row), maxval(stat_row), &
                 ' OUT min/max=', minval(out_row), maxval(out_row)
      end if

      do i = 1, nc
        if (is_statsgo_used(mask_row(i))) then
          n_used = n_used + 1_8
      
          absd = abs(out_row(i) - stat_row(i))
          if (absd > 0.0_4) n_mis_all = n_mis_all + 1_8
      
          if (stat_row(i) >= 0.0_4) then
            n_used_stat_valid = n_used_stat_valid + 1_8
            if (absd > 0.0_4) n_mis_stat_valid = n_mis_stat_valid + 1_8
          end if
        end if
      end do

    end do

    close(u_mask); close(u_out); close(u_stat)

    hole_frac = 0.0d0
    if (n_used > 0_8) hole_frac = dble(n_used - n_used_stat_valid) / dble(n_used)

    print *, 'RESULT ', trim(label)
    print *, '  STATSGO-used pixels           =', n_used
    print *, '  STATSGO-used & STAT valid     =', n_used_stat_valid
    print *, '  Mismatch (all STATSGO-used)   =', n_mis_all
    print *, '  Mismatch (STAT valid only)    =', n_mis_stat_valid
    print *, '  Hole fraction (STAT missing within used mask)=', hole_frac


    deallocate(mask_row, out_row, stat_row)
  end subroutine audit_var

end program audit_step4_vs_statsgo

