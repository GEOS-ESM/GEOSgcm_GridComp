
!#define _memdiag

module memory
#ifdef _CUDA


use iso_c_binding
use cudafor
type adr
    integer*8 :: loc
    integer*8 :: size 
    integer*8 :: gap
    integer :: cindex = 0
    integer :: cnum = 0
    integer :: oindex = 0
    integer :: agn = 0
    type(c_ptr) :: locp
end type

type adrd
    type(c_devptr) :: loc
    integer*8 :: size
    real, device, allocatable :: ar(:)
end type


type(adr) :: plist(500)
type(adr) :: clist(100)
type(adrd) :: dlist(100)
integer :: np = 0
integer :: nc = 0
integer :: acgap = 4
type(c_devptr) :: cpointer

integer :: ddnp = 0
real, device, allocatable :: ddar(:)
real, device :: ddtemp(1)
integer :: ddsizec = 0
integer :: ddindex = 0
integer :: ddflush = 0



interface dbal
    module procedure dbalr, dbalr2, dbalr3, dbali, dbali2, dbali3
end interface 

interface dbcp
    module procedure dbcpi1, dbcpi2, dbcpi3, dbcpr1, dbcpr2, dbcpr3
end interface 

interface ddbxeg
    module procedure ddbxegi, ddbxegr
end interface

contains

subroutine ddbxegi( a, x, y , pt)
    integer, allocatable, device :: a(:,:)
    integer :: x,y
    type(c_devptr), intent(out) :: pt
    

    if (ddflush == 0) then
        
        ddsizec = ddsizec + (x*y)
        !pt = c_devloc( ddtemp(1) )

    else
        
        pt = c_devloc( ddar( ddindex ) )
        ddindex = ddindex + (x*y)
       
    end if
end subroutine



subroutine ddbxegr( a, x, y , pt)
    real, allocatable, device :: a(:,:)
    integer :: x,y
    type(c_devptr), intent(out) :: pt
    

    if (ddflush == 0) then
        
        ddsizec = ddsizec + (x*y)
        pt = c_devloc( ddtemp(1) )

    else
        
        pt = c_devloc( ddar( ddindex ) )
        ddindex = ddindex + (x*y)
       
    end if
end subroutine

subroutine dflush()
    allocate( ddar( ddsizec + 1 ) )
    
    ddflush = 1
    ddindex = 1
end subroutine

subroutine dclean()
    deallocate( ddar )
    ddindex = 0
    ddsizec = 0
    ddflush = 0
end subroutine

    
subroutine dbgenr( p, s )
    real, intent(in) :: p(*)
    integer, intent(in) :: s
    np = np + 1
    plist(np)%loc = loc(p(1))
    plist(np)%locp = c_loc(p(1))
    plist(np)%size = s
    plist(np)%gap = 0
    plist(np)%oindex = np
#ifdef _memdiag
    print *, "index ", np
    print *, "real allocation ", np, " loc: ", plist(np)%loc, " size: ", plist(np)%size
#endif
end subroutine

subroutine dbgeni( p, s )
    integer, intent(in) :: p(*)
    integer, intent(in) :: s
    np = np + 1
    plist(np)%loc = loc(p(1))
    plist(np)%locp = c_loc(p(1))
    plist(np)%size = s
    plist(np)%gap = 0
    plist(np)%oindex = np
#ifdef _memdiag
    print *, "index ", np   
    print *, "integer allocation ", np, " loc: ", plist(np)%loc, " size: ", plist(np)%size
#endif
end subroutine

subroutine dbalr( p )
    real, intent(in) :: p(:)
    call dbgenr( p, size(p) * 4)
end subroutine

subroutine dbalr2( p)
    real, intent(in) :: p(:,:)
    call dbgenr( p, size(p) * 4)
end subroutine

subroutine dbalr3( p)
    real, intent(in) :: p(:,:,:)
    call dbgenr( p, size(p) * 4)
end subroutine

subroutine dbali( p )
    integer, intent(in) :: p(:)
    call dbgeni( p, size(p) * 4)
end subroutine

subroutine dbali2( p )
    integer, intent(in) :: p(:,:)
    call dbgeni( p, size(p) * 4)
end subroutine

subroutine dbali3( p )
    integer, intent(in) :: p(:,:,:)
    call dbgeni( p, size(p) * 4)
end subroutine


subroutine dbflushrg()
    integer :: i,j
    integer*8 :: loc, size, oin
    type(c_ptr) :: locp, cpt
    integer :: cpti
#ifdef _memdiag
    print *, "analyzing memory"
    print *, "sorting entries"
#endif
    do j = 1, np
        do i = 1, np-1

            if (plist(i)%loc > plist(i+1)%loc) then
                loc = plist(i)%loc
                locp = plist(i)%locp
                size = plist(i)%size
                oin = plist(i)%oindex

                plist(i)%loc = plist(i+1)%loc
                plist(i)%locp = plist(i+1)%locp
                plist(i)%size = plist(i+1)%size
                plist(i)%oindex = plist(i+1)%oindex
                plist(i+1)%loc = loc
                plist(i+1)%locp = locp
                plist(i+1)%size = size
                plist(i+1)%oindex = oin
            end if

        end do
    end do

    do i = 1, np - 1
        plist(i)%gap = plist(i+1)%loc - (plist(i)%loc + plist(i)%size)
    end do
    plist(np)%gap = 9999999
#ifdef _memdiag
    print *, "sorted elements"
#endif    
    do i = 1, np 
#ifdef _memdiag
        print *, plist(i)%loc, plist(i)%size, plist(i)%gap
#endif
        if (plist(i)%gap < 0) then
            print *, "ERROR! Memory overlap found at index ", plist(i)%oindex
            stop
        end if
    end do
#ifdef _memdiag
    print *, "analyzing contiguous regions"
#endif
    nc = 1
    clist(1)%loc = plist(1)%loc
    clist(1)%cindex = 1
    do i = 1, np
        plist(i)%cnum = nc
        plist(i)%cindex = clist(nc)%size/4 

        if (plist(i)%gap > acgap) then
            clist(nc)%size = clist(nc)%size + plist(i)%size
            if (i < np) then
                clist(nc+1)%loc = plist(i+1)%loc
                clist(nc+1)%cindex = i+1
            end if
            nc = nc + 1
        else
            clist(nc)%size = clist(nc)%size + plist(i)%size + plist(i)%gap
        end if        

    end do
    nc = nc - 1

#ifdef _memdiag
    print *, "contiguous regions", nc
    print *, "number alloc/copy reduced to ", 100.0 * real(nc)/real(np), "%"

    do i = 1, nc 
        print *, clist(i)%loc, clist(i)%size
    end do

    print *, "allocating device memory"
#endif
    do i = 1, nc
        
        dlist(i)%size = clist(i)%size
#ifdef _memdiag        
        print *, dlist(i)%size
#endif
        allocate( dlist(i)%ar( dlist(i)%size + 2 ))
        dlist(i)%loc = c_devloc( dlist(i)%ar(1) )
    end do

   

end subroutine

subroutine dbcpr( p, pt )
    
    real, intent(in) :: p(*)
    integer*8 :: lc
    type(c_devptr), intent(out) :: pt


end subroutine

subroutine dbcpi1( p, pt )
    integer, intent(in) :: p(:)
    integer*8 :: lc
    type(c_devptr), intent(out) :: pt
    lc = loc(p(1))
    call dbcpg( lc, pt)
end subroutine

subroutine dbcpi2( p, pt )
    integer, intent(in) :: p(:,:)
    integer*8 :: lc
    type(c_devptr), intent(out) :: pt
    lc = loc(p(1,1))
    call dbcpg( lc, pt)
end subroutine

subroutine dbcpi3( p, pt )
    integer, intent(in) :: p(:,:,:)
    integer*8 :: lc
    type(c_devptr), intent(out) :: pt
    lc = loc(p(1,1,1))
    call dbcpg( lc, pt)
end subroutine

subroutine dbcpr1( p, pt )
    real, intent(in) :: p(:)
    integer*8 :: lc
    type(c_devptr), intent(out) :: pt
    lc = loc(p(1))
    call dbcpg( lc, pt)
end subroutine

subroutine dbcpr2( p, pt )
    real, intent(in) :: p(:,:)
    integer*8 :: lc
    type(c_devptr), intent(out) :: pt
    lc = loc(p(1,1))
    call dbcpg( lc, pt)
end subroutine

subroutine dbcpr3( p, pt )
    real, intent(in) :: p(:,:,:)
    integer*8 :: lc
    type(c_devptr), intent(out) :: pt
    lc = loc(p(1,1,1))
    call dbcpg( lc, pt)
end subroutine



subroutine dbcpg( lc, pt )
    integer*8, intent(in) :: lc
    type(c_devptr), intent(out) :: pt
    integer :: fl
    fl = 0
    do i = 1, np

        if (plist(i)%loc .eq. lc) then
#ifdef _memdiag
            print *, "pointer found at index ", i
#endif
            pt = c_devloc( dlist( plist(i)%cnum )%ar( plist(i)%cindex+1 )) 
            fl = 1
            plist(i)%agn = 1
        end if
    end do

    if (fl == 0) then
        print *, "ERROR! pointer not found!"
        stop
    end if

end subroutine


subroutine dbflushcp
    integer :: i
    integer :: err
#ifdef _memdiag   
    print  *, "checking that all pointers are assigned"
#endif
    do i = 1, np
        if (plist(i)%agn == 0) then
            print *, "ERROR! pointer not assigned at index ", plist(i)%oindex
            stop
        end if
    end do
#ifdef _memdiag
    print *, "pointers are OK"
#endif
    do i=1, nc
        err = cudaMemCpyAsync( dlist(i)%loc, plist(clist(i)%cindex)%locp , clist(i)%size+1)
        if (err <> 0) then
            print *, "ERROR! there was an error with a memory copy"
            stop
        end if
    end do
#ifdef _memdiag
    print *, "memory copied successfully"
#endif
end subroutine

subroutine dbclean
    integer :: i
   
    do i=1, nc
        dlist(i)%size=0
        clist(i)%size=0

        deallocate( dlist(i)%ar )
    end do
    nc = 0
    np = 0

end subroutine
#endif
end module

