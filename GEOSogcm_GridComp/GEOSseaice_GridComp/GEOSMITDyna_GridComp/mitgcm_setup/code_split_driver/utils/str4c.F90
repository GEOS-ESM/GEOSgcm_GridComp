       module str4c_mod

       contains

!      Function that puts a string into a 1-byte integer for passing to C
!      Allocates memory for iarr that must be freed.
       subroutine str4c( iarr, str )
!      Args
       integer*1, pointer :: iarr(:)
       character*(*) str

       allocate(iarr(len(str)+1))
       iarr = transfer(str,iarr)
       iarr(len(str)+1)=0

       end subroutine str4c

       end module str4c_mod
