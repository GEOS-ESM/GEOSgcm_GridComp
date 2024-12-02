!
! Usage TileToNC4.x tile_file catch_file
!

program TileToNC4
   use MAPL
   use LogRectRasterizeMod
   use EASE_conv, only: ease_extent
   implicit none
   character(512) :: arg
   integer :: i, unit, unit2

   character(:), allocatable :: tile_file
   character(:), allocatable :: catch_file
   real(kind=8), allocatable   :: rTable(:,:)
   integer, allocatable        :: iTable(:,:)
   character(128) :: gName1, gName2
   character(len=512) :: tmpline

   character(:), allocatable :: array(:)
   character(len=:), allocatable :: filenameNC4
   real :: cell_area
   integer :: n_tile, n_grid, n_lon1, n_lat1, n_cat, tmp_in1, tmp_in2
   integer :: n_lon2, n_lat2, nx, ny, num, ll, maxcat
   logical :: file_exists
 
   CALL get_command_argument(1, arg)
   tile_file  = trim(arg)
   CALL get_command_argument(2, arg)
   catch_file = trim(arg)

   open (newunit=unit, file=trim(tile_file), form='formatted', action='read')
   read (unit,*) tmpline
   read (unit,*) N_grid
   read (unit,*) gName1
   read (unit,*) n_lon1
   read (unit,*) n_lat1

   call split(tmpline, array, " ")  
   read(array(1), *) n_tile
   num = size(array)
   ll = 0
   if (num == 4) then
      ll = 1
      read(array(2), *) maxcat
   else
      maxcat = -1
   endif

   read(array(2+ll), *) nx
   read(array(3+ll), *) ny

   if (N_grid == 1) then
     read (unit,*) tmpline
     if (index(tmpline,'OCEAN')/=0) then
       read (unit,*) 
       read (unit,*) 
       read (unit,*) tmpline
     endif
   else
     read (unit,*) gName2
     read (unit,*) n_lon2
     read (unit,*) n_lat2
     read (unit,*) tmpline
   endif
   if ( index(gName1, 'EASE') /=0) then
     call ease_extent(gName1, tmp_in1, tmp_in2, cell_area=cell_area)
   endif
   ! At this point, the first line is already in tmpline  
   allocate(iTable(N_tile,0:5))
   allocate(rTable(N_tile,10))
   rTable = MAPL_UNDEF
   if ( index(gName1, 'EASE') /=0 ) then
      i = 1
      read(tmpline,*)   iTable(i,0), iTable(i,4), rTable(i,1), rTable(i,2), &
                        iTable(i,2), iTable(i,3), rTable(i,3)
      do i = 2, N_tile
        read (unit,*) tmpline
        read(tmpline,*) iTable(i,0), iTable(i,4), rTable(i,1), rTable(i,2), &
                        iTable(i,2), iTable(i,3), rTable(i,4)
      enddo
      rTable(:,3) = rTable(:,4)*cell_area
      ! rTable(:,4) is fr, now set to area
      ! The fr will be got back in WriteTilingNC4
      rTable(:,4) = cell_area 
   else
      i = 1
      read(tmpline,*)   iTable(i,0), rTable(i,3), rTable(i,1), rTable(i,2), &
                        iTable(i,2), iTable(i,3), rTable(i,4), iTable(i,6), &
                        iTable(i,4), iTable(i,5), rTable(i,5), iTable(i,7)
 
      do i = 2, N_tile
        read (unit,*)   tmpline
        read(tmpline,*) iTable(i,0), rTable(i,3), rTable(i,1), rTable(i,2), &
                        iTable(i,2), iTable(i,3), rTable(i,4), iTable(i,6), &
                        iTable(i,4), iTable(i,5), rTable(i,5), iTable(i,7)
      enddo
   endif
   close(unit)

   inquire(file= trim(catch_file), exist=file_exists)
   if (file_exists) then
      open (newunit=unit, file=trim(catch_file), form='formatted', action='read')
      read(unit, *) n_cat
      do i = 1, n_cat
         read(unit, *)  tmp_in1,      & 
                     tmp_in2,      &
                     rTable(i, 6), &
                     rTable(i, 7), &
                     rTable(i, 8), &
                     rTable(i, 9), &
                     rTable(i, 10)
      enddo
      close(unit)   
  endif

  ll = index(tile_file, '.til')
  filenameNC4 = tile_file(1:ll)//'nc4'
  if (N_grid == 1) then
     call WriteTilingNC4(filenameNc4, [gName1], [n_lon1],[n_lat1], nx, ny, iTable, rTable, maxcat = maxcat) 
  else
     call WriteTilingNC4(filenameNc4, [gName1, gName2], [n_lon1, n_lon2],[n_lat1, n_lat2], nx, ny, iTable, rTable, maxcat=maxcat) 
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
           if(delimiters/='')then                                       ! if DELIMITERS was specified and not null use it
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
   
       case (0)                                                      ! command was totally blank
   
       case default                                                   ! there is at least one non-delimiter in INPUT_LINE if get here
           icol=1                                                      ! initialize pointer into input line
           INFINITE: do i30=1,ilen,1                                   ! store into each array element
               ibegin(i30)=icol                                         ! assume start new token on the character
               if(index(dlim(1:idlim),input_line(icol:icol))==0)then  ! if current character is not a delimiter
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
               if(icol>ilen)then                                     ! no text left
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
       !allocate(array(ireturn))                                       ! allocate the array to turn
   
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
