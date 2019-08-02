#define I_AM_MAIN
#define VERIFY_(A) if(MAPL_VRFY(A,Iam,__LINE__,RC))call MAPL_Abort
!#include "MAPL_Generic.h"

    program ESMF_GenerateCSGridDescription

! ESMF Framework module
    use NetCDF
    use ESMF
    use MAPL_Mod
    use fv_grid_utils_mod, only: gnomonic_grids, cell_center2,get_area, direct_transform
    use fv_grid_tools_mod, only: mirror_grid
    use, intrinsic :: iso_fortran_env, only: REAL64,REAL32

    implicit none
    
! Local variables to be used in the Regrid method calls.
! The code creating and filling these variables is not included in the
! example documentation because those interfaces are not specific to
! Regrid.
    type(ESMF_Grid) ::  dstgrid
    integer :: rc

!EOC

    real(REAL64), parameter           :: PI = 3.14159265358979323846
    character(len=ESMF_MAXSTR)        :: IAm
    integer                           :: npets, localPet
    integer                           :: i, j, k
    type(ESMF_VM)                     :: vm
    real (ESMF_KIND_R8), dimension(2) :: mincoord
    real (ESMF_KIND_R8)               :: deltaX, deltaY
    integer                           :: IM_World, JM_World, scrip_size
    integer, parameter                :: grid_type = 0
    integer, parameter                :: KM_WORLD=1
    integer, parameter                :: NX=1
    integer, parameter                :: NY=6
    integer, parameter                :: ntiles=6
    integer, parameter                :: ndims=2
    integer                           :: N
    integer                           :: info
    integer                           :: start(2), cnt(2), hull_num, hull(4)
    integer                           :: UNIT
    real(ESMF_KIND_R8), allocatable   :: grid_global(:,:,:,:)
    real(ESMF_KIND_R8), allocatable   :: SCRIP_CenterLat(:), SCRIP_CenterLon(:)
    real(ESMF_KIND_R8), allocatable   :: SCRIP_CornerLat(:,:), SCRIP_CornerLon(:,:)
    real(ESMF_KIND_R8), allocatable   :: SCRIP_Area(:)
    real(ESMF_KIND_R8)                :: alocs(2), node_xy(2,4), node_xy_tmp(2,4), lon_e, lon_w
    integer                           :: cornerdim, gridsize, griddim, rankdim, mask
    integer                           :: cornerlon, cornerlat, centerlon, centerlat,cellarea
    integer, allocatable              :: IMS(:), JMS(:), sendData(:), GlobalCounts(:), recvCounts(:), recvOffsets(:)
    integer, allocatable              :: grid_imask(:)
    character(len=ESMF_MAXSTR)        :: gridname, FMT, FMTIM, FMTJM, FMTPET, title
    integer                           :: NPX, NPY
    integer                           :: isg, ieg
    integer                           :: jsg, jeg
    integer                           :: is, ie
    integer                           :: js, je
    integer                           :: myTile
    integer                           :: npts, tmp, mpiC, mpiC2
    integer                           :: IG, JG
    logical                           :: do_schmidt
    real(ESMF_KIND_R8)                :: target_lon, target_lat, stretch_fac
    type(ESMF_Config)                 :: CF
    integer                           :: status

#include "mpif.h"

    Iam = "ESMF_GenerateCSGridDescription"

    call ESMF_Initialize(logKindFlag=ESMF_LOGKIND_NONE,rc=status)
    VERIFY_(status)

    call ESMF_VMGetGlobal(vm, rc=status)
    VERIFY_(status)

! Get number of PETs we are running with
! --------------------------------------
    call ESMF_VMGet(vm, localPet=localPet, petCount=npets, mpiCommunicator=mpiC, rc=status)
    VERIFY_(status)
    call MPI_Comm_dup(mpic, mpic2, status)
    VERIFY_(status)
! Duplicate the MPI communicator not to interfere with ESMF communications.
! The duplicate MPI communicator can be used in any MPI call in the user
! code. Here the MPI_Barrier() routine is called.
    call MPI_Barrier(mpic2, status)
    VERIFY_(status)

    CF = ESMF_ConfigCreate(rc=status)
    VERIFY_(STATUS)
    call ESMF_ConfigLoadFile(CF,filename='GenScrip.rc',rc=status)
    VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute(CF, IM_World, Label='CUBE_DIM:',rc=status)
    VERIFY_(STATUS)
    JM_WORLD = 6 * IM_WORLD
    call ESMF_ConfigGetAttribute(CF, do_schmidt, Label='DO_SCHMIDT:',Default=.false.,rc=status)
    VERIFY_(STATUS)
    call ESMF_ConfigGetAttribute(CF, stretch_fac, Label='STRETCH_FAC:',rc=status)
    if (status /=0) then 
       if (do_schmidt) then
          write(*,*)"Asking for stretch grid without supplying stretch factor"
          call MPI_Abort(mpiC,status)
       end if
    end if
    call ESMF_ConfigGetAttribute(CF, target_lon, Label='TARGET_LON:',rc=status)
    if (status /=0) then 
       if (do_schmidt) then
          write(*,*)"Asking for stretch grid without supplying target lon"
          call MPI_Abort(mpiC,status)
       end if
    end if
    call ESMF_ConfigGetAttribute(CF, target_lat, Label='TARGET_LAT:',rc=status)
    if (status /=0) then 
       if (do_schmidt) then
          write(*,*)"Asking for stretch grid without supplying target lat"
          call MPI_Abort(mpiC,status)
       end if
    end if

! Create destination Cubed-Sphere grid
! ------------------------------------
    allocate( IMS(0:NX-1) )
    allocate( JMS(0:NY-1) )

    call DecomposeDim_ ( IM_WORLD  , IMS             , NX   )
    call DecomposeDim_ ( JM_WORLD/6, JMS(0:NY/6 -1)  , NY/6 )
    do n=2,6
       JMS((n-1)*NY/6 : n*NY/6 -1) = JMS(0:NY/6 -1)
    enddo

    deltaX = 2.0*PI/IM_WORLD
    deltaY = PI/(JM_WORLD  )
    minCoord(1) = 0.0
    minCoord(2) = -PI/2

    call GET_INT_FORMAT_(IM_WORLD, FMTIM)
    call GET_INT_FORMAT_(JM_WORLD, FMTJM)
    FMT = '(A,' // trim(FMTIM) //',A,' // trim(FMTJM) // ',A)'
    write(gridname,trim(FMT)) 'PE',IM_WORLD,'x',JM_WORLD,'-CF'
    dstgrid = ESMF_GridCreate(    &
          name=gridname,                 &
          countsPerDEDim1=ims,           &
          countsPerDEDim2=jms,           &
          indexFlag = ESMF_INDEX_GLOBAL,   &
          coordDep1 = (/1,2/),           &
          coordDep2 = (/1,2/),           &
          rc=status)
    VERIFY_(status)

! Allocate coords at default stagger location
! -------------------------------------------
    call ESMF_GridAddCoord(dstgrid, rc=status)
    VERIFY_(status)

    call ESMF_AttributeSet(dstgrid, 'GridType', 'Cubed-Sphere', rc=status)
    VERIFY_(status)

    call ESMF_GRID_INTERIOR(dstgrid,isg,ieg,jsg,jeg)
!    print*, 'ESMF_GRID_INTERIOR: ', isg, ieg, jsg, jeg

    npx = IM_WORLD
    npy = JM_WORLD
    myTile = jsg/(npy/ntiles)

    is = isg
    ie = ieg
    js = jsg - myTile*(npy/ntiles)
    je = jeg - myTile*(npy/ntiles)
    npts = (npy/ntiles)
    if (npts /= npx) then
       print*, 'Error npts /= npx', npts, npx
       status=1
    endif
    VERIFY_(status)

    print*, 'AppGridCreate: ', myTile, is, ie, js, je, npts

    allocate( grid_global(npts+1,npts+1,ndims,ntiles) )
    grid_global=Z'7FFC000000000000'

    call gnomonic_grids(grid_type, npts, grid_global(:,:,1,1), grid_global(:,:,2,1))

! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi]
! ---------------------------------------------------------------------------------------------
    call mirror_grid(grid_global, 0, npts+1, npts+1, 2, 6)

    do n=1,ntiles
       do j=1,npts+1
           do i=1,npts+1

! Shift the corner away from Japan close to east coast of China
! -------------------------------------------------------------
             if (.not. do_schmidt) grid_global(i,j,1,n) = grid_global(i,j,1,n) - pi/18.
             if ( grid_global(i,j,1,n) < 0. )              &
                 grid_global(i,j,1,n) = grid_global(i,j,1,n) + 2.*pi
             if (ABS(grid_global(i,j,1,n)) < 1.e-10) grid_global(i,j,1,n) = 0.0
             if (ABS(grid_global(i,j,2,n)) < 1.e-10) grid_global(i,j,2,n) = 0.0
          enddo
       enddo
    enddo

! Clean Up Corners
! ----------------
    grid_global(  1,1:npts+1,:,2)=grid_global(npts+1,1:npts+1,:,1)
    grid_global(  1,1:npts+1,:,3)=grid_global(npts+1:1:-1,npts+1,:,1)
    grid_global(1:npts+1,npts+1,:,5)=grid_global(1,npts+1:1:-1,:,1)
    grid_global(1:npts+1,npts+1,:,6)=grid_global(1:npts+1,1,:,1)
    grid_global(1:npts+1,  1,:,3)=grid_global(1:npts+1,npts+1,:,2)
    grid_global(1:npts+1,  1,:,4)=grid_global(npts+1,npts+1:1:-1,:,2)
    grid_global(npts+1,1:npts+1,:,6)=grid_global(npts+1:1:-1,1,:,2)
    grid_global(  1,1:npts+1,:,4)=grid_global(npts+1,1:npts+1,:,3)
    grid_global(  1,1:npts+1,:,5)=grid_global(npts+1:1:-1,npts+1,:,3)
    grid_global(npts+1,1:npts+1,:,3)=grid_global(1,1:npts+1,:,4)
    grid_global(1:npts+1,  1,:,5)=grid_global(1:npts+1,npts+1,:,4)
    grid_global(1:npts+1,  1,:,6)=grid_global(npts+1,npts+1:1:-1,:,4)
    grid_global(  1,1:npts+1,:,6)=grid_global(npts+1,1:npts+1,:,5)

   !------------------------
   ! Schmidt transformation:
   !------------------------
    if ( do_schmidt ) then
      target_lon = target_lon*PI/180._8
      target_lat = target_lat*PI/180._8
      do n=1,ntiles
         call direct_transform(stretch_fac, 1, npts+1, 1, npts+1, target_lon, target_lat, &
                               n, grid_global(1:npts+1,1:npts+1,1,n), grid_global(1:npts+1,1:npts+1,2,n))
      enddo
    endif

    tmp = (ieg-isg+1)*(jeg-jsg+1)
    allocate(SCRIP_CenterLat(tmp),stat=status)
    VERIFY_(status)
    SCRIP_CenterLat=Z'7FFC000000000000'
    allocate(SCRIP_CenterLon(tmp),stat=status)
    VERIFY_(status)
    SCRIP_CenterLon=Z'7FFC000000000000'
    allocate(SCRIP_CornerLat(4,tmp),stat=status)
    VERIFY_(status)
    SCRIP_CornerLat=Z'7FFC000000000000'
    allocate(SCRIP_CornerLon(4,tmp),stat=status)
    VERIFY_(status)
    SCRIP_CornerLon=Z'7FFC000000000000'
    allocate(SCRIP_Area(tmp),stat=status)
    VERIFY_(status)
    SCRIP_Area=Z'7FFC000000000000'

! write a separate GMT format multi-segment polygon file from each process
! can be used later to inspect the grid http://gmt.soest.hawaii.edu/
! ------------------------------------------------------------------------
    call GET_INT_FORMAT_(localPet, FMTPET)
    !FMT = '(A,'// trim(FMTIM) //',A,' // trim(FMTPET) //',A)'
    !write(filename,trim(FMT)) 'c', IM_WORLD,'.',localPet,'.gmt'
    !UNIT=28+localPet
    !open(UNIT,file=trim(filename),form='formatted',status='new')
    n=1
    do jg=jsg,jeg
       do ig=isg,ieg
          i=ig
          j=jg-myTile*npts
          call cell_center2(grid_global(i  ,j  ,1:2,myTile+1), grid_global(i+1,j  ,1:2,myTile+1),   &
                            grid_global(i+1,j+1,1:2,myTile+1), grid_global(i  ,j+1,1:2,myTile+1),   &
                            alocs)
          SCRIP_CenterLon(n) = alocs(1)*(180._8/PI)
          SCRIP_CenterLat(n) = alocs(2)*(180._8/PI)

          node_xy(1,1:4) = (/grid_global(i  ,j,1,myTile+1),grid_global(i+1,j,1,myTile+1),grid_global(i+1,j+1,1,myTile+1),grid_global(i,j+1,1,myTile+1)/)
          node_xy(2,1:4) = (/grid_global(i  ,j,2,myTile+1),grid_global(i+1,j,2,myTile+1),grid_global(i+1,j+1,2,myTile+1),grid_global(i,j+1,2,myTile+1)/)
          node_xy_tmp = node_xy

! Correct for the periodic boundary at 0/360
! ------------------------------------------
          lon_w = min( grid_global(i  ,j,1,myTile+1),grid_global(i+1,j,1,myTile+1),grid_global(i+1,j+1,1,myTile+1),grid_global(i,j+1,1,myTile+1) ) 
          lon_e = max( grid_global(i  ,j,1,myTile+1),grid_global(i+1,j,1,myTile+1),grid_global(i+1,j+1,1,myTile+1),grid_global(i,j+1,1,myTile+1) ) 
          if ( abs(lon_e - lon_w) > 1.5_8*pi .and. (SCRIP_CenterLon(n) < pi) ) then
             where(node_xy(1,:) > pi) node_xy_tmp(1,:) = node_xy(1,:) - 2._8*pi
          elseif ( abs(lon_e - lon_w) > 1.5_8*pi .and. (SCRIP_CenterLon(n) > pi) ) then
             where(node_xy(1,:) < pi) node_xy_tmp(1,:) = node_xy(1,:) + 2._8*pi
          endif
          call points_hull_2d(4, node_xy_tmp, hull_num, hull)
          if(ANY(hull==0)) then
             write(*,100)'Zero Hull ', grid_global(i  ,j,1,myTile+1),grid_global(i+1,j,1,myTile+1),grid_global(i+1,j+1,1,myTile+1),grid_global(i,j+1,1,myTile+1)
             write(*,100)'Zero Hull ', node_xy_tmp(1,:)
          endif
          !write(UNIT,103) '>'
          do k=1,4
            SCRIP_CornerLon(k,n) = node_xy(1,hull(k))*(180._8/PI)
            SCRIP_CornerLat(k,n) = node_xy(2,hull(k))*(180._8/PI) 
            !write(UNIT,102) SCRIP_CornerLon(k,n), SCRIP_CornerLat(k,n)
          enddo
          !SCRIP_Area(n) = get_area(node_xy(:,hull(1)),node_xy(:,hull(2)),node_xy(:,hull(3)),node_xy(:,hull(4)),1.0d0)
          SCRIP_Area(n) = get_area(grid_global(i  ,j  ,1:2,myTile+1), grid_global(i,j+1  ,1:2,myTile+1),   &
                            grid_global(i+1,j,1:2,myTile+1), grid_global(i+1,j+1,1:2,myTile+1),1.0d0)
          n=n+1
       enddo
    enddo

 100     format(a,4f20.15)
 101     format(a,f20.15)
 102     format(2f20.15)
 103     format(a)

    deallocate( grid_global )
    deallocate( IMS )
    deallocate( JMS )

    scrip_size = IM_World*JM_World
    call MPI_Info_create(info, status)
    VERIFY_(status)
    call MPI_Info_set(info, "cb_buffer_size", "1048576", status)
    VERIFY_(status)

    status = nf90_create(trim(gridname)//'.nc4', IOR(NF90_MPIIO,IOR(NF90_CLOBBER,NF90_NETCDF4)),unit,comm=mpic2,info=info)
    if(status /= nf90_noerr) then
       print*,'Error creating file ',status
       print*, NF90_STRERROR(status)
       stop
    endif

    FMT = '(A,' // ',A,' //',A)'
    write(title,trim(FMT)) 'GMAO ',trim(gridname),' Grid'
    status = nf90_put_att(unit,NF90_GLOBAL, 'title', trim(title))
    if(status /= nf90_noerr) then
       print*,'Error setting title',status
       print*, NF90_STRERROR(status)
       stop
    endif
    status = nf90_put_att(unit,NF90_GLOBAL, 'GridDescriptionFormat','SCRIP')
    if(status /= nf90_noerr) then
       print*,'Error setting GridDescriptionFormat',status
       print*, NF90_STRERROR(status)
       stop
    endif
    if (do_Schmidt) then
       status = nf90_put_att(unit,NF90_GLOBAL, 'TargetLon', target_lon*180._8/PI)
       if(status /= nf90_noerr) then
          print*,'Error setting TargetLon',status
          print*, NF90_STRERROR(status)
          stop
       endif
       status = nf90_put_att(unit,NF90_GLOBAL, 'TargetLat', target_lat*180._8/PI)
       if(status /= nf90_noerr) then
          print*,'Error setting TargetLat',status
          print*, NF90_STRERROR(status)
          stop
       endif
       status = nf90_put_att(unit,NF90_GLOBAL, 'StretchFactor', stretch_fac)
       if(status /= nf90_noerr) then
          print*,'Error setting StretchFactor',status
          print*, NF90_STRERROR(status)
          stop
       endif
    end if

    status = NF90_def_dim(UNIT, 'grid_size', SCRIP_size, gridsize) 
    if(status /= nf90_noerr) then
       print*,'Error defining grid_size',status
       print*, NF90_STRERROR(status)
       stop
    endif
    status = nf90_def_dim(unit,'grid_corners', 4, cornerdim)
    if(status /= nf90_noerr) then
       print*,'Error defining grid_corners',status
       print*, NF90_STRERROR(status)
       stop
    endif

! Peggy Li suggested setting grid_rank=1 and grid_dims=1
! so that ESMF will treat the grid as unstructured
! ------------------------------------------------------
    status = NF90_DEF_DIM(UNIT, 'grid_rank'  , 1, rankdim)
    if(status /= nf90_noerr) then
       print*,'Error defining grid_rank',status
       print*, NF90_STRERROR(status)
       stop 
    endif

! Grid dimensions
! ---------------
    status = nf90_def_var(UNIT, "grid_dims", NF90_INT, [rankdim], griddim)
    if(status /= nf90_noerr) then
       print*,'Error defining grid_dims',status
       print*, NF90_STRERROR(status)
       stop 
    endif

! Grid mask
! ---------
    status = nf90_def_var(UNIT, "grid_imask", NF90_INT, [gridsize], mask)
    if(status /= nf90_noerr) then
       print*,'Error defining grid_imask',status
       print*, NF90_STRERROR(status)
       stop
    endif
    status = nf90_put_att(unit,mask,"units","unitless")

! cell center Longitude variable
! ------------------------------
    status = nf90_def_var(UNIT, "grid_center_lon", NF90_DOUBLE, [gridsize], centerlon)
    if(status /= nf90_noerr) then
       print*,'Error defining cell center lons',status
       print*, NF90_STRERROR(status)
       stop 
    endif
    status = nf90_put_att(UNIT, centerlon, "units"    ,  "degrees")

! cell center Latitude variable
! -----------------------------
    status = nf90_def_var(UNIT, "grid_center_lat", NF90_DOUBLE, [gridsize], centerlat)
    if(status /= nf90_noerr) then
       print*,'Error defining cell center lats',status
       print*, NF90_STRERROR(status)
       stop
    endif
    status = nf90_put_att(UNIT, centerlat, "units"    , "degrees")

! cell corner Longitude variable
! ------------------------------
    status = nf90_def_var(UNIT, "grid_corner_lon", NF90_DOUBLE, [cornerdim,gridsize], cornerlon)
    if(status /= nf90_noerr) then
       print*,'Error defining cell corner lons',status
       print*, NF90_STRERROR(status)
       stop
    endif
    status = nf90_put_att(UNIT, cornerlon, "units"  ,   "degrees")

! cell corner Latitude variable
! -----------------------------
    status = nf90_def_var(UNIT, "grid_corner_lat", NF90_DOUBLE, [cornerdim,gridsize], cornerlat)
    if(status /= nf90_noerr) then
       print*,'Error defining cell corner lats',status
       print*, NF90_STRERROR(status)
       stop
    endif
    status = nf90_put_att(UNIT, cornerlat, "units"    ,"degrees")

    status = nf90_def_var(UNIT, "grid_area", NF90_DOUBLE, [gridsize], cellarea)
    if(status /= nf90_noerr) then
       print*,'Error defining cell area',status
       print*, NF90_STRERROR(status)
       stop
    endif
    status = nf90_put_att(UNIT, cellarea, "units"    , "radians^2")

    status = nf90_enddef(UNIT)
    if(status /= nf90_noerr) then
       print*,'Error exitting define mode',status
       print*, NF90_STRERROR(status)
       stop
    endif

    rc = NF90_PUT_VAR(UNIT, griddim, 1)
    if(rc /= nf90_noerr) then
       print*,'Error writing griddim',rc
       print*, NF90_STRERROR(rc)
       stop
    endif

    allocate (sendData(1),GlobalCounts(npets), recvCounts(npets), recvOffsets(npets), stat=status)
    VERIFY_(status)
    sendData = tmp
    recvCounts=1
    recvOffsets=0
    do i=2, npets
      recvOffsets(i) = recvOffsets(i-1) + recvCounts(i-1)
    enddo
    call ESMF_VMGatherV(vm,sendData=sendData,sendCount=1,recvData=GlobalCounts,recvCounts=recvCounts,recvOffsets=recvOffsets,rootPet=0,rc=status)
    VERIFY_(status)
    call ESMF_VMBroadcast(vm,bcstData=GlobalCounts,count=npets,rootPet=0, rc=status)
    VERIFY_(status)
    if (localPet == 0) print*,GlobalCounts

    start=1
    do i=1,localPet
       start(1) = start(1)+GlobalCounts(i)
    enddo

    cnt(1) = tmp; cnt(2)=1
    status = NF90_PUT_VAR(UNIT, centerlon, SCRIP_CenterLon, start=start, count=cnt)
    if(status /= nf90_noerr) then
       print*,'Error writing CenterLons ',status
       print*, NF90_STRERROR(status)
       stop
    endif

    status = NF90_PUT_VAR(UNIT, centerlat, SCRIP_CenterLat, start=start, count=cnt)
    if(status /= nf90_noerr) then
       print*,'Error writing CenterLats ',status
       print*, NF90_STRERROR(status)
       stop
    endif

    status = NF90_PUT_VAR(UNIT, cellarea, SCRIP_Area, start=start, count=cnt)
    if(status /= nf90_noerr) then
       print*,'Error writing CenterLats ',status
       print*, NF90_STRERROR(status)
       stop
    endif
! Grid mask
! ---------
    allocate(grid_imask(tmp), stat=status)
    VERIFY_(status)
    grid_imask = 1
    status = NF90_PUT_VAR(UNIT, mask, grid_imask, start=start, count=cnt)
    if(status /= nf90_noerr) then
       print*,'Error writing grid_imask',status
       print*, NF90_STRERROR(status)
       stop
    endif
    deallocate(grid_imask)

    start(2)=start(1)
    start(1)=1
    cnt(1)=4
    cnt(2)=tmp

    status = NF90_PUT_VAR(UNIT, cornerlat, SCRIP_CornerLat, start=start, count=cnt)
    if(status /= nf90_noerr) then
       print*,'Error writing CornerLats ',status
       print*, NF90_STRERROR(status)
       stop
    endif

    status = NF90_PUT_VAR(UNIT, cornerlon, SCRIP_CornerLon, start=start, count=cnt)
    if(status /= nf90_noerr) then
       print*,'Error writing CornerLons ',status
       print*, NF90_STRERROR(status)
       stop
    endif

    status = NF90_CLOSE(UNIT)
    if(status /= nf90_noerr) then
       print*,'Error closing output file ', status
       print*, NF90_STRERROR(status)
       stop
    endif
    call ESMF_VMBarrier(vm, rc=status)

    deallocate(SCRIP_CenterLat)
    deallocate(SCRIP_CenterLon)
    deallocate(SCRIP_CornerLat)
    deallocate(SCRIP_CornerLon)
    deallocate(sendData)
    deallocate(GlobalCounts)
    deallocate(recvCounts)
    deallocate(recvOffsets)

    call ESMF_Finalize(rc=status)

    contains

    subroutine GET_INT_FORMAT_(N, FMT)
       integer          :: N
       character(len=*) :: FMT

       IF(N < 10) THEN
          FMT = 'I1'
       ELSE IF (N< 100) THEN
          FMT = 'I2'
       ELSE IF (N< 1000) THEN
          FMT = 'I3'
       ELSE IF (N< 10000) THEN
          FMT = 'I4'
       else
          FMT = 'I5'
       end IF
    end subroutine GET_INT_FORMAT_

    subroutine DecomposeDim_( dim_world,dim,NDEs )
!
! From FMS/MPP
!
       implicit   none
       integer    dim_world, NDEs
       integer    dim(0:NDEs-1)

       integer :: is,ie,isg,ieg
       integer :: ndiv,ndivs,imax,ndmax,ndmirror,n
       integer :: ibegin(0:NDEs-1)
       integer :: iend(0:NDEs-1)

       logical :: symmetrize
       logical :: even, odd
       even(n) = (mod(n,2).EQ.0)
       odd (n) = (mod(n,2).EQ.1)

       isg = 1
       ieg = dim_world
       ndivs = NDEs
 
       is = isg
       n = 0
       do ndiv=0,ndivs-1
          !modified for mirror-symmetry
          !original line
          !                 ie = is + CEILING( float(ieg-is+1)/(ndivs-ndiv) ) - 1

          !problem of dividing nx points into n domains maintaining symmetry
          !i.e nx=18 n=4 4554 and 5445 are solutions but 4455 is not.
          !this will always work for nx even n even or odd
          !this will always work for nx odd, n odd
          !this will never  work for nx odd, n even: for this case we supersede the mirror calculation
          !                 symmetrize = .NOT. ( mod(ndivs,2).EQ.0 .AND. mod(ieg-isg+1,2).EQ.1 )
          !nx even n odd fails if n>nx/2
          symmetrize = ( even(ndivs) .AND. even(ieg-isg+1) ) .OR. &
                       (  odd(ndivs) .AND.  odd(ieg-isg+1) ) .OR. &
                       (  odd(ndivs) .AND. even(ieg-isg+1) .AND. ndivs.LT.(ieg-isg+1)/2 )

          !mirror domains are stored in the list and retrieved if required.
          if( ndiv.EQ.0 )then
             !initialize max points and max domains
             imax = ieg
             ndmax = ndivs
          end if
          !do bottom half of decomposition, going over the midpoint for odd ndivs
          if( ndiv.LT.(ndivs-1)/2+1 )then
             !domain is sized by dividing remaining points by remaining domains
             ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
             ndmirror = (ndivs-1) - ndiv !mirror domain
             if( ndmirror.GT.ndiv .AND. symmetrize )then !only for domains over the midpoint
                !mirror extents, the max(,) is to eliminate overlaps
                ibegin(ndmirror) = max( isg+ieg-ie, ie+1 )
                iend(ndmirror)   = max( isg+ieg-is, ie+1 )
                imax = ibegin(ndmirror) - 1
                ndmax = ndmax - 1
             end if
          else
             if( symmetrize )then
                !do top half of decomposition by retrieving saved values
                is = ibegin(ndiv)
                ie = iend(ndiv)
             else
                ie = is + CEILING( REAL(imax-is+1)/(ndmax-ndiv) ) - 1
             end if
          end if
          dim(ndiv) = ie-is+1
          is = ie + 1
       end do

  end subroutine DecomposeDim_

!logical function MAPL_VRFY(A,iam,line,rc)
   !integer,           intent(IN ) :: A
   !character*(*),     intent(IN ) :: iam
   !integer,           intent(IN ) :: line
   !integer, optional, intent(OUT) :: RC
     !MAPL_VRFY = A/=ESMF_SUCCESS
     !if(MAPL_VRFY)then
       !if(present(RC)) then
         !print'(A40,I10)',Iam,line
         !RC=A
       !endif
     !endif
!end function MAPL_VRFY

  !subroutine MAPL_Abort
    !call ESMF_Finalize(endflag = ESMF_END_ABORT)
  !end subroutine MAPL_Abort

  !subroutine ESMF_GRID_INTERIOR(GRID,I1,IN,J1,JN)
    !type (ESMF_Grid), intent(IN) :: grid
    !integer, intent(OUT)         :: I1, IN, J1, JN

!! local vars
    !integer                               :: status
    !character(len=ESMF_MAXSTR)            :: IAm='ESMF_GridInterior'

    !type (ESMF_DistGrid)                  :: distGrid
    !type(ESMF_DELayout)                   :: LAYOUT
    !integer,               allocatable    :: AL(:,:)
    !integer,               allocatable    :: AU(:,:)
    !integer                               :: nDEs
    !integer                               :: deId
    !integer                               :: gridRank
    !integer                               :: deList(1)

    !call ESMF_GridGet    (GRID, dimCount=gridRank, distGrid=distGrid, rc=STATUS)
    !call ESMF_DistGridGet(distGRID, delayout=layout, rc=STATUS)
    !call ESMF_DELayoutGet(layout, deCount =nDEs, localDeList=deList, rc=status)
    !deId = deList(1)

    !allocate (AL(gridRank,0:nDEs-1),  stat=status)
    !allocate (AU(gridRank,0:nDEs-1),  stat=status)

    !call ESMF_DistGridGet(distgrid, &
         !minIndexPDe=AL, maxIndexPDe=AU, rc=status)

    !I1 = AL(1, deId)
    !IN = AU(1, deId)
!!    ASSERT_(gridRank > 1) !ALT: tilegrid is 1d (without RC this only for info)
    !J1 = AL(2, deId)
    !JN = AU(2, deId)
    !deallocate(AU, AL)

  !end subroutine ESMF_GRID_INTERIOR

!function angle_rad_2d ( p1, p2, p3 )
subroutine angle_rad_2d ( p1, p2, p3, res )
!*****************************************************************************80
!
!! ANGLE_RAD_2D returns the angle swept out between two rays in 2D.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
!
!        P1
!        /
!       /
!      /
!     /
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( REAL64 ) P1(2), P2(2), P3(2), define the rays
!    P1 - P2 and P3 - P2 which define the angle.
!
!    Output, real ( REAL64 ) ANGLE_RAD_2D, the angle swept out by the rays,
!    in radians.  0 <= ANGLE_RAD_2D < 2 * PI.  If either ray has zero
!    length, then ANGLE_RAD_2D is set to 0.
  implicit none

  integer, parameter :: dim_num = 2

  real    (REAL64), parameter :: pi = 3.141592653589793D+00
  real    (REAL64) p(dim_num)
  real    (REAL64) p1(dim_num)
  real    (REAL64) p2(dim_num)
  real    (REAL64) p3(dim_num)
  real    (REAL64) res

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )


  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( p(1) == 0.0D+00 .and. p(2) == 0.0D+00 ) then
    res = 0.0D+00
    return
  end if

  res = atan2 ( p(2), p(1) )

  if ( res < 0.0D+00 ) then
    res = res + 2.0D+00 * pi
  end if

  return
end

subroutine points_hull_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_HULL_2D computes the convex hull of 2D points.
!
!  Discussion:
!
!    The work involved is N*log(H), where N is the number of points, and H is
!    the number of points that are on the hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NODE_NUM, the number of nodes.
!
!    Input, real ( REAL64 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Output, integer  HULL_NUM, the number of nodes that lie on
!    the convex hull.
!
!    Output, integer  HULL(NODE_NUM).  Entries 1 through HULL_NUM
!    contain the indices of the nodes that form the convex hull, in order.
!
  implicit none

  integer node_num
  real    (REAL64) angle
  real    (REAL64) angle_max
  real    (REAL64) di
  real    (REAL64) dr
  integer  first
  integer  hull(node_num)
  integer  hull_num
  integer  i
  real    (REAL64) node_xy(2,node_num)
  real    (REAL64) p_xy(2)
  integer  q
  real    (REAL64) q_xy(2)
  integer  r
  real    (REAL64) r_xy(2)

  if ( node_num < 1 ) then
    hull_num = 0
    return
  end if
!
!  If NODE_NUM = 1, the hull is the point.
!
  if ( node_num == 1 ) then
    hull_num = 1
    hull(1) = 1
    return
  end if
!
!  If NODE_NUM = 2, then the convex hull is either the two distinct points,
!  or possibly a single (repeated) point.
!
  if ( node_num == 2 ) then

    if ( node_xy(1,1) /= node_xy(1,2) .or. node_xy(2,1) /= node_xy(2,2) ) then
      hull_num = 2
      hull(1) = 1
      hull(2) = 2
    else
      hull_num = 1
      hull(1) = 1
    end if

    return

  end if
!
!  Find the leftmost point and call it "Q".
!  In case of ties, take the bottom-most.
!
  q = 1
  do i = 2, node_num
    if ( node_xy(1,i) < node_xy(1,q) .or. &
       ( node_xy(1,i) == node_xy(1,q) .and. node_xy(2,i) < node_xy(2,q) ) ) then
      q = i
    end if
  end do

  q_xy(1:2) = node_xy(1:2,q)
!
!  Remember the starting point, so we know when to stop!
!
  first = q
  hull_num = 1
  hull(1) = q
!
!  For the first point, make a dummy previous point, 1 unit south,
!  and call it "P".
!
  p_xy(1) = q_xy(1)
  p_xy(2) = q_xy(2) - 1.0D+00
!
!  Now, having old point P, and current point Q, find the new point R
!  so the angle PQR is maximal.
!
!  Watch out for the possibility that the two nodes are identical.
!
  do

    r = 0
    angle_max = 0.0D+00

    do i = 1, node_num

      if ( i /= q .and. &
           ( node_xy(1,i) /= q_xy(1) .or. node_xy(2,i) /= q_xy(2) ) ) then

        call angle_rad_2d(p_xy, q_xy, node_xy(1:2,i),angle)
        !angle = angle_rad_2d ( p_xy, q_xy, node_xy(1:2,i) )

        if ( r == 0 .or. angle_max < angle ) then

          r = i
          r_xy(1:2) = node_xy(1:2,r)
          angle_max = angle
!
!  In case of ties, choose the nearer point.
!
        else if ( r /= 0 .and. angle == angle_max ) then

          di = ( node_xy(1,i) - q_xy(1) )**2 + ( node_xy(2,i) - q_xy(2) )**2
          dr = ( r_xy(1)      - q_xy(1) )**2 + ( r_xy(2)      - q_xy(2) )**2

          if ( di < dr ) then
            r = i
            r_xy(1:2) = node_xy(1:2,r)
            angle_max = angle
          end if

        end if

      end if

    end do
!
!  We are done when we have returned to the first point on the convex hull.
!
    if ( r == first ) then
      exit
    end if

    hull_num = hull_num + 1
    if ( node_num < hull_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POINTS_HULL_2D - Fatal error!'
      write ( *, '(a)' ) '  The algorithm has failed.'
      stop
    end if
!
!  Add point R to convex hull.
!
    hull(hull_num) = r
!
!  Set P := Q, Q := R, and prepare to search for next point R.
!
    q = r

    p_xy(1:2) = q_xy(1:2)
    q_xy(1:2) = r_xy(1:2)

  end do

  return
end

    end program ESMF_GenerateCSGridDescription
