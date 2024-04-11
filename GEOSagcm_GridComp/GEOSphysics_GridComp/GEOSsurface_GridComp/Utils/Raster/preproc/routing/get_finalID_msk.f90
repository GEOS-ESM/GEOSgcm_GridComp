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
  
  res=int8(pfaf)
  pfaf_digit(:,1)=res/(int8(10)**int8(11))
  do i=2,12
     res=res-int8(10)**int8(13-i)*int8(pfaf_digit(:,i-1))
     pfaf_digit(:,i)=res/(int8(10)**int8(12-i))
  enddo
  
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
           fulli=pfaf_digit(i,:)
           do j=i-1,1,-1
              ok=1
              fullj=pfaf_digit(j,:)
              samed=0
              do k=1,min(pfaf_last(i),pfaf_last(j))
                 if(fulli(k)==fullj(k))then
                    samed=samed+1
                 else
                    exit
                 endif
              enddo
              if(samed+1<=pfaf_last(j))then
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
           enddo
        else
           downid(i)=-1  
        endif
        deallocate(code)
        
     endif
     
     
  enddo


  open(88,file="outputs/Pfaf_downid.txt")
  do i=1,nc
     write(88,*)downid(i)
  enddo
  
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
