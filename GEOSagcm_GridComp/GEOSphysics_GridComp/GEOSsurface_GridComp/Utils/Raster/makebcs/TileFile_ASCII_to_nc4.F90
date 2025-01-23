!
! Utility program that converts ASCII-formatted *.til file and catchment.def file into a single nc4 file
!
! Usage TileFile_ASCII-to-nc4.x tile_file catchmentdef_file
!
! wjiang, rreichle, 29 Nov 2024

program TileFile_ASCII_to_nc4
  
  use MAPL
  use LogRectRasterizeMod, only: WriteTilingNC4
  use EASE_conv,           only: ease_extent
  
  implicit none
  
  character(512)                 :: arg
  integer                        :: i, unit, unit2
  
  character(:),      allocatable :: tile_file
  character(:),      allocatable :: catchmentdef_file
  real(kind=8),      allocatable :: rTable(:,:)
  integer,           allocatable :: iTable(:,:)
  character(128)                 :: gName1, gName2
  character(len=512)             :: tmpline
  
  character(:),      allocatable :: array(:)
  character(len=:),  allocatable :: filenameNC4
  
  real                           :: cell_area
  
  integer                        :: n_tile, n_grid, n_lon1, n_lat1, n_cat, tmp_in1, tmp_in2
  integer                        :: n_lon2, n_lat2, nx, ny, num, ll, maxcat
  logical                        :: file_exists
  
  ! ----------------------------------------------------------------------
  !
  ! process command-line arguments

  CALL get_command_argument(1, arg)
  tile_file  = trim(arg)
  CALL get_command_argument(2, arg)
  catchmentdef_file = trim(arg)
  
  ! ----------------------------------------------------------------------
  !
  ! open and read *.til ASCII file

  open (newunit=unit, file=trim(tile_file), form='formatted', action='read')

  read (unit,*) tmpline            ! header line 1: N_tile [maxcat]  nx  ny             (see below)
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
     
     read (unit,*) tmpline
     if (index(tmpline,'OCEAN')/=0) then
        read (unit,*) 
        read (unit,*) 
        read (unit,*) tmpline
     endif

  else

     ! lat/lon or cube-sphere tile space

     read (unit,*) gName2
     read (unit,*) n_lon2
     read (unit,*) n_lat2
     read (unit,*) tmpline     ! read info for first tile (to accommodate legacy EASE grid issues above)

  endif

  allocate(iTable(N_tile,0:7))
  allocate(rTable(N_tile,10))
  
  rTable = MAPL_UNDEF  

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

     call ease_extent(gName1, tmp_in1, tmp_in2, cell_area=cell_area)  ! get EASE grid cell area
     
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

  inquire( file= trim(catchmentdef_file), exist=file_exists)

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

  ! assemble name of nc4 file

  ll = index(tile_file, '.til')
  filenameNC4 = tile_file(1:ll)//'nc4'

  ! write nc4 file
  
  if (N_grid == 1) then
     call WriteTilingNC4(filenameNc4, [gName1        ], [n_lon1        ], [n_lat1        ], nx, ny, iTable, rTable, N_PfafCat=maxcat) 
  else
     call WriteTilingNC4(filenameNc4, [gName1, gName2], [n_lon1, n_lon2], [n_lat1, n_lat2], nx, ny, iTable, rTable, N_PfafCat=maxcat) 
  endif
  
contains
  
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
  
end program

! ======================= EOF ====================================================
