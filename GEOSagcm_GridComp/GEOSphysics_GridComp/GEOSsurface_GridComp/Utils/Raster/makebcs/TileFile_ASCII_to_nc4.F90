!
! Utility program that converts ASCII-formatted *.til file and catchment.def file into a single nc4 file
!
! Usage TileFile_ASCII-to-nc4.x tile_file catchmentdef_file
!  LIS: TileFile_ASCII_to_nc4.x tile_file --lis native_im native_jm bounds_file
!
! wjiang, rreichle, 29 Nov 2024

program TileFile_ASCII_to_nc4
  use, intrinsic :: iso_fortran_env, only: REAL64 
  use MAPL,                          only: MAPL_WriteTilingNC4, MAPL_ease_extent 
  use netcdf
  use LogRectRasterizeMod,           only: MAPL_UNDEF_R8
  
  implicit none
  
  character(512)                 :: arg
  integer                        :: i, unit, unit2
  
  character(:),      allocatable :: tile_file
  character(:),      allocatable :: catchmentdef_file
  real(REAL64),      allocatable :: rTable(:,:)
  integer,           allocatable :: iTable(:,:)
  character(128)                 :: gName1, gName2
  character(len=512)             :: tmpline
  
  character(:),      allocatable :: array(:)
  character(len=:),  allocatable :: filenameNC4
  
  real                           :: cell_area
  
  integer                        :: n_tile, n_grid, n_lon1, n_lat1, n_cat, tmp_in1, tmp_in2
  integer                        :: n_lon2, n_lat2, nx, ny, num, ll, maxcat 
  logical                        :: file_exists, lis_mode
  integer                        :: lis_im, lis_jm
  character(512)                 :: bounds_file


  ! ----------------------------------------------------------------------
  !
  ! process command-line arguments

  CALL get_command_argument(1, arg)
  tile_file  = trim(arg)
  CALL get_command_argument(2, arg)
  catchmentdef_file = trim(arg)
  lis_mode = trim(arg) == '--lis'
  
  if (lis_mode) then
     call read_lis_arguments(lis_im, lis_jm, bounds_file)
  end if  
  ! ----------------------------------------------------------------------
  !
  ! open and read *.til ASCII file

  open (newunit=unit, file=trim(tile_file), form='formatted', action='read')

  read (unit,"(A)") tmpline            ! header line 1: N_tile [maxcat]  nx  ny             (see below)
  read (unit,*) N_grid             ! header line 2: N_grid [=1 for EASE, =2 otherwise]
  read (unit,*) gName1             ! header line 3: name  of atm grid
  read (unit,*) n_lon1             ! header line 4: N_lon of atm grid
  read (unit,*) n_lat1             ! header line 5: N_lat of atm grid
  
  ! special treatment needed for header line 1 because maxcat is not included in legacy bcs
  
  call split(tmpline, array, " ")  
  read(array(1), *) n_tile
  num = size(array)
  ll = 0
  if (num == 4) then
     ll = 1
     read(array(2), *) maxcat      ! number of Pfafstetter catchments 
  else
     maxcat = -1                   ! maxcat not available in legacy bcs
  endif
  
  read(array(2+ll), *) nx          ! N_lon of raster grid
  read(array(3+ll), *) ny          ! N_lat of raster grid
  
  if (N_grid == 1) then

     ! EASE grid tile space

     ! in some legacy bcs, dummy ocean grid info is included in header (despite N_grid=1);
     ! read next line and decide if it is dummy header or info for first tile
     
     read (unit,"(A)") tmpline
     if (index(tmpline,'OCEAN')/=0) then
        read (unit,*) 
        read (unit,*) 
        read (unit,"(A)") tmpline
     endif

  else

     ! lat/lon or cube-sphere tile space

     read (unit,*) gName2
     read (unit,*) n_lon2
     read (unit,*) n_lat2
     read (unit,"(A)") tmpline     ! read info for first tile (to accommodate legacy EASE grid issues above)

  endif

  if (.not. lis_mode .and. N_grid == 2) then
     if (index(gName2, 'LIS') /= 0) then
        error stop 'LIS geometry requires explicit --lis conversion'
     end if
  end if
  
  if (lis_mode) then
  
     if (N_grid /= 2 .or. index(gName2, 'LIS') == 0) then
        error stop 'Expected two-grid LIS geometry'
     end if
  
     if (maxcat /= lis_im * lis_jm) then
        error stop 'ASCII native land-ID bound does not match LIS dimensions'
     end if
  
  end if  

  allocate(iTable(N_tile,0:7))
  if (lis_mode) iTable = 0
  allocate(rTable(N_tile,10))
  
  rTable = MAPL_UNDEF_r8  

  ! read ASCII tile file (NOTE: Info for first tile is already in tmpline!)
  
  if ( index(gName1, 'EASE') /=0 ) then  ! EASE grid tile space
     
     read (tmpline,*) iTable(1,0), iTable(1,4), rTable(1,1), rTable(1,2),   &
                      iTable(1,2), iTable(1,3), rTable(1,4)

     do i = 2, N_tile
        read (unit,*) iTable(i,0), iTable(i,4), rTable(i,1), rTable(i,2),   &
                      iTable(i,2), iTable(i,3), rTable(i,4)
     enddo

     ! rTable(:,4) is tile area fraction within grid cell (fr), convert to area;
     ! get fr back in WriteTilingNC4

     call MAPL_ease_extent(gName1, tmp_in1, tmp_in2, cell_area=cell_area)  ! get EASE grid cell area
     
     rTable(:,3) = rTable(:,4)*cell_area
     rTable(:,4) = cell_area 
     
  else  ! lat/lon or cube-sphere tile space
     
     read (tmpline,*) iTable(1,0), rTable(1,3), rTable(1,1), rTable(1,2),   &
                      iTable(1,2), iTable(1,3), rTable(1,4), iTable(1,6),   &
                      iTable(1,4), iTable(1,5), rTable(1,5), iTable(1,7)
     
     do i = 2, N_tile
        read (unit,*) iTable(i,0), rTable(i,3), rTable(i,1), rTable(i,2),   &
                      iTable(i,2), iTable(i,3), rTable(i,4), iTable(i,6),   &
                      iTable(i,4), iTable(i,5), rTable(i,5), iTable(i,7)
     enddo

     ! re-define rTable(:,4) and rTable(:,5).
     ! fr will be re-created in WriteTilingNC4

     where (rTable(:,4) /=0.0)
       rTable(:,4) = rTable(:,3)/rTable(:,4)
     endwhere

     where (rTable(:,5) /=0.0)
       rTable(:,5) = rTable(:,3)/rTable(:,5)
     endwhere

  endif

  close(unit)
  
  ! ----------------------------------------------------------------------
  !
  ! open and read catchment.def ASCII file

  file_exists = .false.
  if (.not. lis_mode) inquire(file=trim(catchmentdef_file),exist=file_exists)

  if (file_exists) then

     open (newunit=unit, file=trim(catchmentdef_file), form='formatted', action='read')

     read(unit, *) n_cat    ! number of *land* tiles

     do i = 1, n_cat
        read(unit, *)   &
          tmp_in1,      & 
          tmp_in2,      &
          rTable(i, 6), &
          rTable(i, 7), &
          rTable(i, 8), &
          rTable(i, 9), &
          rTable(i,10)
     enddo

     close(unit)

  endif
  if (lis_mode) then
     call read_lis_bounds(bounds_file, rTable, iTable, maxcat)
  end if

  ! assemble name of nc4 file

  ll = index(tile_file, '.til')
  filenameNC4 = tile_file(1:ll)//'nc4'

  ! write nc4 file
  
  if (N_grid == 1) then
     call MAPL_WriteTilingNC4(filenameNc4, [gName1        ], [n_lon1        ], [n_lat1        ], nx, ny, iTable, rTable, N_PfafCat=maxcat) 
  else
     call MAPL_WriteTilingNC4(filenameNc4, [gName1, gName2], [n_lon1, n_lon2], [n_lat1, n_lat2], nx, ny, iTable, rTable, N_PfafCat=maxcat) 
  endif
  if (lis_mode) call label_lis_geometry(filenameNC4,iTable(:,0),iTable(:,4),lis_im,lis_jm)

contains

  subroutine nc_check(status, context)
  
    integer, intent(in)        :: status
    character(*), intent(in)   :: context
  
    if (status /= nf90_noerr) then
       write(*, '(A)') trim(context)//': '//trim(nf90_strerror(status))
       error stop 'Cannot finish LIS geometry NetCDF'
    end if
  
  end subroutine nc_check  

  subroutine label_lis_geometry(filename, typ, source_id, im, jm)

    character(*), intent(in) :: filename
    integer, intent(in) :: typ(:), source_id(:), im, jm

    integer :: ncid, legacy_var, ndims, tile_count
    integer :: dimids(nf90_max_var_dims)
    integer :: lis_index_var, lis_i_var, lis_j_var
    integer, allocatable :: lis_index(:), lis_i(:), lis_j(:)

    ! Keep MAPL's existing geometry schema for compatibility and add
    ! explicit LIS indices. Legacy pfaf_index and N_PfafCat are
    ! compatibility slots only, not Pfafstetter hydrological identifiers.

    call nc_check( &
         nf90_open(trim(filename), nf90_write, ncid), 'open converted geometry')

    call nc_check( &
         nf90_inq_varid(ncid, 'pfaf_index', legacy_var), 'MAPL geometry pfaf_index slot')

    call nc_check( &
         nf90_inquire_variable(ncid, legacy_var, ndims=ndims, dimids=dimids), 'tile dimension')

    if (ndims /= 1) error stop 'Unexpected MAPL geometry index rank'

    call nc_check( &
         nf90_inquire_dimension(ncid, dimids(1), len=tile_count), 'tile count')

    if (tile_count /= size(typ)) error stop 'MAPL geometry tile count mismatch'

    call nc_check(nf90_redef(ncid), 'metadata definition')

    call nc_check( &
         nf90_put_att(ncid, nf90_global, 'land_discretization', 'LIS'), 'land_discretization')

    call nc_check( &
         nf90_put_att(ncid, nf90_global, 'geometry_only', 1), 'geometry_only')

    call nc_check( &
         nf90_put_att(ncid, nf90_global, 'N_LIS_IDs', im * jm), 'N_LIS_IDs')

    call nc_check( &
         nf90_put_att(ncid, nf90_global, 'LIS_im', im), 'LIS_im')

    call nc_check( &
         nf90_put_att(ncid, nf90_global, 'LIS_jm', jm), 'LIS_jm')

    call nc_check( &
         nf90_put_att(ncid, nf90_global, 'land_id_convention', &
         'LIS_ID=(row-1)*LIS_im+col; one-based; rows south-to-north, columns west-to-east'), &
         'ID convention')

    call nc_check( &
         nf90_put_att(ncid, nf90_global, 'routing_warning', &
         'Not Pfafstetter geometry. Do not use with GEOS Catchment/ROUTE without a separate mapping.'), &
         'routing warning')

    call nc_check( &
         nf90_put_att(ncid, nf90_global, 'legacy_Pfaf_metadata', &
         'For land only, pfaf_index and N_PfafCat are compatibility slots for LIS IDs/count; use lis_index.'), &
         'legacy note')

    call nc_check( &
         nf90_put_att(ncid, nf90_global, 'elevation_note', &
         'Tile elevation is undefined; LIS supplies land parameters independently.'), &
         'elevation note')

    call nc_check( &
         nf90_put_att(ncid, legacy_var, 'long_name', &
         'Legacy source-index slot: native LIS cell ID on land; NOT a Pfafstetter index'), &
         'legacy variable label')

    call nc_check( &
         nf90_def_var(ncid, 'lis_index', nf90_int, dimids(1:1), lis_index_var), 'lis_index')

    call nc_check( &
         nf90_def_var(ncid, 'lis_i', nf90_int, dimids(1:1), lis_i_var), 'lis_i')

    call nc_check( &
         nf90_def_var(ncid, 'lis_j', nf90_int, dimids(1:1), lis_j_var), 'lis_j')

    call nc_check( &
         nf90_put_att(ncid, lis_index_var, 'long_name', &
         'Native LIS cell ID; zero for nonland tiles'), 'lis_index name')

    call nc_check( &
         nf90_put_att(ncid, lis_i_var, 'long_name', &
         'One-based LIS longitude column; zero for nonland'), 'lis_i name')

    call nc_check( &
         nf90_put_att(ncid, lis_j_var, 'long_name', &
         'One-based LIS latitude row; zero for nonland'), 'lis_j name')

    call nc_check(nf90_enddef(ncid), 'end definitions')

    allocate(lis_index(size(typ)), lis_i(size(typ)), lis_j(size(typ)))

    lis_index = 0
    lis_i = 0
    lis_j = 0

    where (typ == 100)
       lis_index = source_id
       lis_i = modulo(source_id - 1, im) + 1
       lis_j = (source_id - 1) / im + 1
    end where

    call nc_check(nf90_put_var(ncid, lis_index_var, lis_index), 'write lis_index')
    call nc_check(nf90_put_var(ncid, lis_i_var, lis_i), 'write lis_i')
    call nc_check(nf90_put_var(ncid, lis_j_var, lis_j), 'write lis_j')

    call nc_check(nf90_close(ncid), 'close LIS geometry')

  end subroutine label_lis_geometry
  
  subroutine split(input_line,array,delimiters,order,nulls)
     
    character(len=*),intent(in)              :: input_line 
    character(len=*),optional,intent(in)     :: delimiters 
    character(len=*),optional,intent(in)     :: order      
    character(len=*),optional,intent(in)     :: nulls      
    character(len=:),allocatable,intent(out) :: array(:)   
    
    integer                       :: n                     
    integer,allocatable           :: ibegin(:)             
    integer,allocatable           :: iterm(:)              
    character(len=:),allocatable  :: dlim                  
    character(len=:),allocatable  :: ordr                  
    character(len=:),allocatable  :: nlls                  
    integer                       :: ii,iiii               
    integer                       :: icount                
    integer                       :: ilen                  
    integer                       :: i10,i20,i30           
    integer                       :: icol                  
    integer                       :: idlim                 
    integer                       :: ifound                
    integer                       :: inotnull              
    integer                       :: ireturn               
    integer                       :: imax                      
    
    ! decide on value for optional DELIMITERS parameter
    if (present(delimiters)) then                                     ! optional delimiter list was present
       if(delimiters/='')then                                         ! if DELIMITERS was specified and not null use it
          dlim=delimiters
       else                                                           ! DELIMITERS was specified on call as empty string
          dlim=' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0) ! use default delimiter when not specified
       endif
    else                                                              ! no delimiter value was specified
       dlim=' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0)    ! use default delimiter when not specified
    endif
    idlim=len(dlim)                                                   ! dlim a lot of blanks on some machines if dlim is a big string
    
    if(present(order))then; ordr=adjustl(order); else; ordr='sequential'; endif ! decide on value for optional ORDER parameter
       
    if(present(nulls))then; nlls=adjustl(nulls); else; nlls='ignore'    ; endif ! optional parameter
           
    n=len(input_line)+1                        ! max number of strings INPUT_LINE could split into if all delimiter
    allocate(ibegin(n))                        ! allocate enough space to hold starting location of tokens if string all tokens
    allocate(iterm(n))                         ! allocate enough space to hold ending location of tokens if string all tokens
    ibegin(:)=1
    iterm(:)=1
    
    ilen=len(input_line)                                           ! ILEN is the column position of the last non-blank character
    icount=0                                                       ! how many tokens found
    inotnull=0                                                     ! how many tokens found not composed of delimiters
    imax=0                                                         ! length of longest token found
    
    select case (ilen)
       
    case (0)                                                       ! command was totally blank
       
    case default                                                   ! there is at least one non-delimiter in INPUT_LINE if get here
       icol=1                                                      ! initialize pointer into input line
       INFINITE: do i30=1,ilen,1                                   ! store into each array element
          ibegin(i30)=icol                                         ! assume start new token on the character
          if(index(dlim(1:idlim),input_line(icol:icol))==0)then    ! if current character is not a delimiter
             iterm(i30)=ilen                                       ! initially assume no more tokens
             do i10=1,idlim                                        ! search for next delimiter
                ifound=index(input_line(ibegin(i30):ilen),dlim(i10:i10))
                IF(ifound>0)then
                   iterm(i30)=min(iterm(i30),ifound+ibegin(i30)-2)
                endif
             enddo
             icol=iterm(i30)+2                                     ! next place to look as found end of this token
             inotnull=inotnull+1                                   ! increment count of number of tokens not composed of delimiters
          else                                                     ! character is a delimiter for a null string
             iterm(i30)=icol-1                                     ! record assumed end of string. Will be less than beginning
             icol=icol+1                                           ! advance pointer into input string
          endif
          imax=max(imax,iterm(i30)-ibegin(i30)+1)
          icount=i30                                               ! increment count of number of tokens found
          if(icol>ilen)then                                        ! no text left
             exit INFINITE
          endif
       enddo INFINITE
       
    end select
    
    select case (trim(adjustl(nlls)))
    case ('ignore','','ignoreend')
       ireturn=inotnull
    case default
       ireturn=icount
    end select
    allocate(character(len=imax) :: array(ireturn))                ! allocate the array to return
    !allocate(array(ireturn))                                      ! allocate the array to turn
    
    select case (trim(adjustl(ordr)))                              ! decide which order to store tokens
    case ('reverse','right') ; ii=ireturn ; iiii=-1                ! last to first
    case default             ; ii=1       ; iiii=1                 ! first to last
    end select
    
    do i20=1,icount                                                ! fill the array with the tokens that were found
       if(iterm(i20)<ibegin(i20))then
          select case (trim(adjustl(nlls)))
          case ('ignore','','ignoreend')
          case default
             array(ii)=' '
             ii=ii+iiii
          end select
       else
          array(ii)=input_line(ibegin(i20):iterm(i20))
          ii=ii+iiii
       endif
    enddo
  end subroutine split

  subroutine read_lis_arguments(im, jm, bounds_file)
  
    integer, intent(out)          :: im
    integer, intent(out)          :: jm
    character(*), intent(out)     :: bounds_file
  
    integer                       :: status
    character(512)                :: argument
  
    if (command_argument_count() /= 5) then
       error stop 'LIS conversion requires native_im native_jm bounds_file'
    end if
  
    call get_command_argument(3, argument)
    read(argument, *, iostat=status) im
  
    if (status /= 0) then
       error stop 'Invalid LIS longitude count'
    end if
  
    call get_command_argument(4, argument)
    read(argument, *, iostat=status) jm
  
    if (status /= 0) then
       error stop 'Invalid LIS latitude count'
    end if
  
    if (im <= 0 .or. jm <= 0) then
       error stop 'Invalid LIS dimensions'
    end if
  
    call get_command_argument(5, bounds_file)
  
  end subroutine read_lis_arguments  

  subroutine read_lis_bounds(filename, rTable, iTable, max_land_id)
  
    character(*), intent(in)       :: filename
    real(REAL64), intent(inout)    :: rTable(:,:)
    integer, intent(in)            :: iTable(:,0:)
    integer, intent(in)            :: max_land_id
  
    integer                        :: unit
    integer                        :: n_tile_bounds
    integer                        :: i
  
    open(newunit=unit, file=trim(filename), status='old', action='read')
  
    read(unit, *) n_tile_bounds
  
    if (n_tile_bounds /= size(iTable, 1)) then
       error stop 'Tile bounds count mismatch'
    end if
  
    do i = 1, size(iTable, 1)
       read(unit, *) rTable(i,6:9)
    end do
  
    close(unit)
  
    do i = 1, size(iTable, 1)
  
       if (iTable(i,0) == 100) then
          if (iTable(i,4) < 1 .or. iTable(i,4) > max_land_id) then
             error stop 'Invalid LIS ID in ASCII tiling'
          end if
       end if
  
    end do
  
  end subroutine read_lis_bounds  
end program

! ======================= EOF ====================================================
