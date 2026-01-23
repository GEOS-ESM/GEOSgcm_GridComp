program main

  use routing_constant,only : nc
  implicit none
  
  integer,allocatable,dimension(:)   :: downid,finalid
  real*8,allocatable,dimension(:)    :: pfaf
  integer,allocatable,dimension(:,:) :: pfaf_digit
  integer*8,allocatable,dimension(:) :: res
  integer,allocatable,dimension(:)   :: pfaf_last,pfaf_msk,code,behind
  integer,allocatable,dimension(:)   :: first,last
  
  integer :: i,j,jj,k,p,down,cur,idx,num,ok,samed
  integer :: fulli(12),fullj(12)
  real    :: val(9)
  
  character(len=100) :: file_path !/discover/nobackup/projects/gmao/bcs_shared/make_bcs_inputs/land/topo/v1/SRTM-TopoData/Pfafcatch-routing.dat

! Get downstream catchment and final destination ID for each catchment, and determine whether it directs to an ocean or inland lake.
! downid=Pfafstetter index of catchment just downstream
! finalid=Pfafstetter index of catchment at outlet point
! pfaf= Pfafstetter number for catchment
! pfaf_digit= The 12 digits in a Pfafstetter number, separated
! pfaf_last= The index of the last nonzero digit in a Pfafstetter number (counting from the left)
! pfaf_msk =1 for non-sink catchments, 2 for sink catchments with endpoints in ocean, =3 for sink catchments with endpoints in interior lake
! last= The index of the last digit in a Pfafstetter number after removing any 11..000 tail.
! first= The index of the last zero (but not the zero at the very end).  However, if there are no zeroes until the end, first =2 (the second index, since the first index indicates the continent).

  if (command_argument_count() /= 1) then
      print *, "no <file_path> found"
      stop
  endif
  call get_command_argument(1, file_path)

  open(77,file=file_path, form="formatted", status="old")
  read(77,*)num
  
  allocate(downid(nc),finalid(nc),pfaf(nc),pfaf_digit(nc,12),res(nc),pfaf_last(nc),pfaf_msk(nc))
  allocate(first(nc),last(nc))
  
  do i=1,nc
     read(77,*)idx,pfaf(i)
  enddo
  
! Separate Pfafstetter number into individual digits  
  res=int8(pfaf)
  pfaf_digit(:,1)=res/(int8(10)**int8(11))
  do i=2,12
     res=res-int8(10)**int8(13-i)*int8(pfaf_digit(:,i-1))
     pfaf_digit(:,i)=res/(int8(10)**int8(12-i))
  enddo

! Determine positions of last nonzero digit (pfaf_last) and the last digit that’s neither 0 nor 1 (at the end)  
  first=2
  last=2
  do i=1,nc
     do j=12,1,-1
        if(pfaf_digit(i,j)/=0)then
           pfaf_last(i)=j
           do k=0,j-1
              if(pfaf_digit(i,j-k)/=1)then
                 last(i)=j-k
                 exit
              endif
           enddo
           exit
        endif
     enddo
  enddo
  do i=1,nc
     if(last(i)<=1) last(i)=2
  enddo

! Determine position of final zero that has some nonzero digits after it  
  do i=1,nc
     do j=last(i),2,-1
        if(pfaf_digit(i,j)==0)then
           first(i)=j
           exit
        endif
     enddo
  enddo
  
  do i=1,nc
     
     if(first(i)>last(i)-1)then
        downid(i)=-1
     else
        
        allocate(code(1:last(i)-first(i)))
        code=pfaf_digit(i,first(i):last(i)-1)
        if(any(code==2).or.any(code==4).or.any(code==6).or.any(code==8))then 
      ! If all digits (after the first) are odd, the Pfafstetter logic implies that the catchment will be on the coast.
           fulli=pfaf_digit(i,:)
           do j=i-1,1,-1  ! Test each catchment to see if it lies just downstream of catchment i 
              ok=1
              fullj=pfaf_digit(j,:)
              samed=0
              do k=1,min(pfaf_last(i),pfaf_last(j)) ! Determine the index (samed) up to which the Pfaf numbers of catchment I and j match  
                 if(fulli(k)==fullj(k))then
                    samed=samed+1
                 else
                    exit
                 endif
              enddo ! end k loop
              if(samed+1<=pfaf_last(j))then
              ! Check that none of catchment j’s indices (after samed) are even, which would imply a downstream branching off from the river on which catchment i lies.  
                 allocate(behind(1:pfaf_last(j)-samed))
                 behind=fullj(samed+1:pfaf_last(j))
                 if(any(mod(behind,2)==0)) ok=0
                 deallocate(behind)
              else
                 ok=0
              endif
              if(ok==1)then
                 downid(i)=j
                 exit
              endif
           enddo ! end j loop
        else
           downid(i)=-1  
        endif
        deallocate(code)
        
     endif ! end i loop
     
     
  enddo


  open(88,file="outputs/Pfaf_downid.txt")
  do i=1,nc
     write(88,*)downid(i)
  enddo

! Keep “moving downstream” until you find the catchment with no downstream catchment:  
  do i=1,nc
     cur=i	
     down=downid(i)
     do while(down/=-1)
  	cur=down
        down=downid(cur)
     enddo
     finalid(i)=cur
  enddo
  
  open(88,file="outputs/Pfaf_finalID.txt")
  do i=1,nc
     write(88,*)finalid(i)
  enddo
  
  
  ! Set masks: 1 = has downstream catchment; 2 = drains to ocean; 3 = drains to inland lake  
  do i=1,nc
     if(downid(i)/=-1)then
        pfaf_msk(i)=1
     else
  	allocate(code(1:pfaf_last(i)))
        code=pfaf_digit(i,1:pfaf_last(i))
        if(any(code==0))then
           pfaf_msk(i)=3
        else
           pfaf_msk(i)=2
        endif
        deallocate(code)
     end if
  enddo
  
  open(88,file="outputs/Pfaf_msk.txt")
  do i=1,nc
     write(88,*)pfaf_msk(i)
  enddo
  
  
  
end program main
