program main
!main purpose: Reads the Pfafstetter code dataset and generates files for the connectivity of catchments in the routing network.

  use constant, only: nc, nupmax

  implicit none

  ! Declare allocatable arrays for routing and Pfafstetter information:
  integer, allocatable, dimension(:)      :: downid, finalid
  real*8, allocatable, dimension(:)       :: pfaf         ! Pfafstetter number for each catchment
  integer, allocatable, dimension(:,:)    :: pfaf_digit, upstream
  integer*8, allocatable, dimension(:)    :: res          ! Temporary storage for digit extraction
  integer, allocatable, dimension(:)      :: pfaf_last, pfaf_msk, code, behind
  integer, allocatable, dimension(:)      :: first, last, nup, nts, nts_old
  real, allocatable, dimension(:)         :: pfaf_area, pfaf_acar, pfaf_acar_old
  
  ! Declare loop and temporary variables:
  integer                                 :: i, j, jj, k, p, down, cur, idx, num, ok, samed, did, nmax
  integer                                 :: fulli(12), fullj(12)
  real                                    :: val(5)
  
  ! Define file path for input routing data:
  character(len=900)                      :: file_path !"input/Pfafcatch-routing.dat"   

  if (command_argument_count() /= 1) then
      print *, "no <file_path> found"
      stop
  endif
  call get_command_argument(1, file_path)

  !---------------------------------------------------------------------------
  ! Read routing data from the input file:
  open(77, file=trim(file_path), form="formatted", status="old")
  read(77, *) num
  
  ! Allocate arrays based on the total number of catchments (nc):
  allocate(downid(nc), finalid(nc), pfaf(nc), pfaf_digit(nc,12), res(nc), &
           pfaf_last(nc), pfaf_msk(nc), pfaf_area(nc))
  allocate(first(nc), last(nc))
  
  do i = 1, nc
     read(77, *) idx, pfaf(i), val(1:5), pfaf_area(i)
  end do
  
  !---------------------------------------------------------------------------
  ! Separate the Pfafstetter number into its 12 individual digits.
  res = int8(pfaf)  ! Convert Pfafstetter numbers to 64-bit integers
  pfaf_digit(:,1) = res / (int8(10) ** int8(11))
  do i = 2, 12
     res = res - int8(10) ** int8(13-i) * int8(pfaf_digit(:, i-1))
     pfaf_digit(:, i) = res / (int8(10) ** int8(12-i))
  end do

  !---------------------------------------------------------------------------
  ! Determine the positions of the last nonzero digit (pfaf_last)
  ! and the position of the last digit that is not 1 (stored in 'last').
  first = 2   ! Initialize 'first' to 2 by default
  last = 2    ! Initialize 'last' to 2 by default
  do i = 1, nc
     do j = 12, 1, -1
        if (pfaf_digit(i, j) /= 0) then
           pfaf_last(i) = j
           do k = 0, j-1
              if (pfaf_digit(i, j-k) /= 1) then
                 last(i) = j - k
                 exit
              endif
           end do
           exit
        endif
     end do
  end do
  do i = 1, nc
     if (last(i) <= 1) last(i) = 2
  end do

  !---------------------------------------------------------------------------
  ! Determine the position of the final zero that has nonzero digits after it.
  do i = 1, nc
     do j = last(i), 2, -1
        if (pfaf_digit(i, j) == 0) then
           first(i) = j
           exit
        endif
     end do
  end do
  
  !---------------------------------------------------------------------------
  ! Determine the immediate downstream catchment for each catchment.
  do i = 1, nc
     
     if (first(i) > last(i) - 1) then
        ! No valid downstream digit exists; mark as terminal (sink)
        downid(i) = -1
     else
        
        allocate(code(1 : last(i) - first(i)))
        code = pfaf_digit(i, first(i) : last(i)-1)
        if (any(code == 2) .or. any(code == 4) .or. any(code == 6) .or. any(code == 8)) then 
           ! If any digit in the extracted part is even, then the catchment is non-coastal.
           fulli = pfaf_digit(i, :)
           do j = i-1, 1, -1  ! Loop backward to find a catchment just downstream of catchment i
              ok = 1
              fullj = pfaf_digit(j, :)
              samed = 0
              do k = 1, min(pfaf_last(i), pfaf_last(j))
                 if (fulli(k) == fullj(k)) then
                    samed = samed + 1
                 else
                    exit
                 endif
              end do  ! End of k loop: number of matching leading digits stored in 'samed'
              if (samed + 1 <= pfaf_last(j)) then
                 ! Check that none of catchment j's remaining digits (after the common part)
                 ! are even, which would indicate a branching downstream.
                 allocate(behind(1 : pfaf_last(j) - samed))
                 behind = fullj(samed+1 : pfaf_last(j))
                 if (any(mod(behind, 2) == 0)) ok = 0
                 deallocate(behind)
              else
                 ok = 0
              endif
              if (ok == 1) then
                 downid(i) = j  ! Found the immediate downstream catchment for catchment i
                 exit
              endif
           end do  ! End of j loop
        else
           downid(i) = -1  ! If extracted digits are not even, mark as sink (or coastal)
        endif
        deallocate(code)
        
     endif  ! End if for determining downstream catchment for catchment i
     
  end do

  !---------------------------------------------------------------------------
  ! Write the downstream catchment IDs to an output file:
  open(88, file="output/downstream_1D_new_noadj.txt")
  do i = 1, nc
     write(88, *) downid(i)
  end do

  ! Write catchment areas to an output file:
  open(88, file="output/Pfaf_area.txt")
  do i = 1, nc
     write(88, *) pfaf_area(i)
  end do

  !---------------------------------------------------------------------------
  ! Build an upstream connectivity matrix:
  allocate(upstream(nupmax, nc), nup(nc))
  nup = 0
  upstream = -1
  do i = 1, nc
    did = downid(i)
    if (did >= 1) then
      nup(did) = nup(did) + 1
      upstream(nup(did), did) = i
    end if
  end do
  open(88, file="output/upstream_1D.txt")
  do i = 1, nc
    write(88, '(34(I8))') upstream(:, i)
  end do
  open(88, file="output/Pfaf_upnum.txt")
  do i = 1, nc
    write(88, *) nup(i)
  end do

  !---------------------------------------------------------------------------
  ! Calculate the number of steps (nts) from each catchment to the sink:
  allocate(nts(nc), pfaf_acar(nc))
  nts = -9999
  do i = 1, nc
    k = 0
    cur = i
    do while (downid(cur) /= -1)
      k = k + 1
      cur = downid(cur)
    end do
    nts(i) = k
  end do
  open(88, file="output/Pfaf_tosink.txt")
  do i = 1, nc
    write(88, *) nts(i)
  end do
  
  !---------------------------------------------------------------------------
  ! Aggregate catchment areas along the flow network:
  nmax = maxval(nts)
  pfaf_acar = pfaf_area
  do j = nmax, 1, -1
    do i = 1, nc  
      if (nts(i) == j) then
        did = downid(i)
        pfaf_acar(did) = pfaf_acar(did) + pfaf_acar(i)
      endif
    end do
  end do
  open(88, file="temp/Pfaf_acar.txt")
  do i = 1, nc
    write(88, *) pfaf_acar(i)
  end do

end program main
