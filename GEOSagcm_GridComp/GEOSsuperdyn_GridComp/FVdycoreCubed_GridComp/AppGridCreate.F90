! routine to create the global edges and centers of the cubed-sphere grid
! but not the ESMF grid
subroutine AppCSEdgeCreateF(IM_WORLD, LonEdge,LatEdge, LonCenter, LatCenter, rc)
#include "MAPL_Generic.h"

   use ESMF
   use MAPL_BaseMod
   use MAPL_GenericMod
   use MAPL_ConstantsMod, only : pi=> MAPL_PI_R8
   use MAPL_ErrorHandlingMod
   use fv_arrays_mod,     only: REAL4, REAL8, R_GRID
   use fv_grid_utils_mod, only: gnomonic_grids, cell_center2, direct_transform
   use fv_grid_tools_mod, only: mirror_grid
   use FV_StateMod, only: FV_Atm
   implicit none

   ! !ARGUMENTS:
   integer,           intent(IN)     :: IM_WORLD
   integer, optional, intent(OUT)    :: rc
   real(ESMF_KIND_R8), intent(inout) :: LonEdge(IM_World+1,IM_World+1,6)
   real(ESMF_KIND_R8), intent(inout) :: LatEdge(IM_World+1,IM_World+1,6)
   real(ESMF_KIND_R8), optional, intent(inout) :: LonCenter(IM_World,IM_World)
   real(ESMF_KIND_R8), optional, intent(inout) :: LatCenter(IM_World,IM_World)

   ! ErrLog variables
   !-----------------

   integer                      :: STATUS
   character(len=ESMF_MAXSTR), parameter :: Iam="AppCSEdgeCreateF"

   ! Local variables
   !-----------------
#ifdef EIGHT_BYTE
   integer, parameter:: f_p = selected_real_kind(15)   ! same as 12 on Altix
#else
   ! Higher precisions for grid geometrical factors:
   integer, parameter:: f_p = selected_real_kind(20)
#endif
   integer, parameter            :: grid_type = 0
   integer                       :: npts
   integer                       :: ntiles=6
   integer                       :: ndims=2
   integer                       :: I, J, N
   integer                       :: IG, JG
   real(ESMF_KIND_R8)            :: alocs(2)

   real(R_GRID), allocatable :: grid_global(:,:,:,:)

   integer                         :: L

   npts = IM_World
   allocate( grid_global(npts+1,npts+1,ndims,ntiles) )
   call gnomonic_grids(grid_type, npts, grid_global(:,:,1,1), grid_global(:,:,2,1))
   ! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi]
   call mirror_grid(grid_global, 0, npts+1, npts+1, 2, 6)

   if (allocated(FV_Atm)) then
      if (.not. FV_Atm(1)%flagstruct%do_schmidt) &
            grid_global(:,:,1,:) = grid_global(:,:,1,:) - PI/FV_Atm(1)%flagstruct%shift_fac
   else
      grid_global(:,:,1,:) = grid_global(:,:,1,:) - PI/18.0
   end if

   where(grid_global(:,:,1,:) < 0.) &
         grid_global(:,:,1,:) = grid_global(:,:,1,:) + 2.* PI

   ! Keep Equator and Greenwich exact
   !---------------------------------

   where(abs(grid_global(:,:,:,1)) < 1.d-10) grid_global(:,:,:,1) = 0.0

   !---------------------------------
   ! Clean Up Corners
   !---------------------------------
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
   if (allocated(FV_Atm)) then
      if ( FV_Atm(1)%flagstruct%do_schmidt ) then
         do n=1,ntiles
            call direct_transform(FV_Atm(1)%flagstruct%stretch_fac, 1, npts+1, 1, npts+1, FV_Atm(1)%flagstruct%target_lon, FV_Atm(1)%flagstruct%target_lat, &
                  n, grid_global(1:npts+1,1:npts+1,1,n), grid_global(1:npts+1,1:npts+1,2,n))
         enddo
      endif
   end if

   if (present(LonCenter) .and. present(LatCenter)) then
      do n=1,ntiles
         do j=1,npts
            do i=1,npts
               call cell_center2(grid_global(i,j,  1:2,n), grid_global(i+1,j,  1:2,n),   &
                     grid_global(i,j+1,1:2,n), grid_global(i+1,j+1,1:2,n),   &
                     alocs)
               jg = (n-1)*npts + j
               LonCenter(i,jg) = alocs(1)
               LatCenter(i,jg) = alocs(2)
            enddo
         enddo
      enddo
   end if

   LonEdge = grid_global(:,:,1,:)
   LatEdge = grid_global(:,:,2,:)

   deallocate( grid_global )

   _RETURN(ESMF_SUCCESS)
end subroutine AppCSEdgeCreateF


!!!!!!!!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function AppGridCreateF(IM_WORLD, JM_WORLD, LM, NX, NY, rc) result(esmfgrid)
#include "MAPL_Generic.h"
#define DEALLOCGLOB_(A) if(associated(A))then;A=0;if(MAPL_ShmInitialized)then; call MAPL_DeAllocNodeArray(A,rc=STATUS);else; deallocate(A,stat=STATUS);endif;_VERIFY(STATUS);NULLIFY(A);endif

   use ESMF
   use MAPL_Mod
   use MAPL_BaseMod
   use MAPL_GenericMod
   use MAPL_ConstantsMod, only : pi=> MAPL_PI_R8
   use MAPL_ShmemMod

   use fv_arrays_mod,     only: REAL4, REAL8, R_GRID
   use fv_grid_utils_mod, only: gnomonic_grids, cell_center2, direct_transform
   use fv_grid_tools_mod, only: mirror_grid
   use FV_StateMod, only: FV_Atm

   implicit none

   ! !ARGUMENTS:
   integer,           intent(IN)    :: IM_WORLD, JM_WORLD, LM
   integer,           intent(IN)    :: NX, NY
   integer, optional, intent(OUT)   :: rc
   type (ESMF_Grid)                 :: esmfgrid

   ! ErrLog variables
   !-----------------

   integer                      :: STATUS
   character(len=ESMF_MAXSTR), parameter :: Iam="AppGridCreateF"

   ! Local variables
   !-----------------

   type(ESMF_DistGrid)           :: distGrid
   integer                       :: DECOUNT(2)
   integer                       :: GLOBAL_DE_START(2)
   integer                       :: NX2, NY2
   integer                       :: NPES, NPES_X, NPES_Y
   integer                       :: NPX, NPY
   integer                       :: isg, ieg
   integer                       :: jsg, jeg
   integer                       :: is, ie
   integer                       :: js, je
   integer                       :: myTile
   integer                       :: npts
   integer                       :: ntiles, grid_type
   integer                       :: ndims=2
   integer                       :: I, J, N
   integer                       :: IG, JG

   real(REAL8) :: deglon=0.0
   type (ESMF_Array),  pointer     :: coords(:)
   real(REAL8), pointer     :: lons(:,:)
   real(REAL8), pointer     :: lats(:,:)
   type (ESMF_Array),  target      :: tarray(2)
   real(REAL8)              :: alocs(2), dx, dy

   integer                         :: L
   integer, allocatable            :: IMS(:), JMS(:)
   character(len=ESMF_MAXSTR)      :: gridname, FMT, FMTIM, FMTJM
   type (ESMF_VM)                  :: VM
   type (ESMF_Config)              :: cf
   real(REAL8), allocatable :: gridCornerLats(:)
   real(REAL8), allocatable :: gridCornerLons(:)
   integer                         :: idx
   real(R_GRID), pointer     :: grid_global(:,:,:,:) => null()
   logical                  :: AppGridStandAlone

   ! We need the VM to create the grid ??
   !------------------------------------

   AppGridStandAlone = .false.
   call ESMF_VMGetCurrent(vm, rc=STATUS)
   _VERIFY(STATUS)            

   ! Check FV3 grid_type
   ! -------------------
   cf = ESMF_ConfigCreate(rc=rc)
   call ESMF_ConfigLoadFile( cf, 'fvcore_layout.rc', rc = rc )
   call ESMF_ConfigGetAttribute( cf, ntiles, label = 'ntiles:', default=6, rc = rc )
   call ESMF_ConfigGetAttribute( cf, grid_type, label = 'grid_type:', default=0, rc = rc )
   call ESMF_ConfigDestroy(cf, rc=rc)

   ! Give the IMS and JMS the MAPL default distribution
   ! --------------------------------------------------

   allocate( IMS(0:NX-1) )
   allocate( JMS(0:NY-1) )

   call GET_INT_FORMAT(IM_WORLD, FMTIM)
   call GET_INT_FORMAT(JM_WORLD, FMTJM)
   FMT = '(A,' // trim(FMTIM) //',A,' // trim(FMTJM) // ',A)'
   if (grid_type <= 3) then
      ! Cubed-Sphere
      call DecomposeDim ( IM_WORLD  , IMS             , NX   )
      call DecomposeDim ( JM_WORLD/6, JMS(0:NY/6 -1)  , NY/6 )
      do n=2,6
         JMS((n-1)*NY/6 : n*NY/6 -1) = JMS(0:NY/6 -1)
      enddo
      write(gridname,trim(FMT)) 'PE',IM_WORLD,'x',JM_WORLD,'-CF'
   else
      ! Doubly Periodic Cartesian
      call DecomposeDim ( IM_WORLD  , IMS             , NX   )
      call DecomposeDim ( JM_WORLD  , JMS             , NY   )
      write(gridname,trim(FMT)) 'PE',IM_WORLD,'x',JM_WORLD,'-DP'
   endif

   ! We should have a much simpler create with the next ESMF grid design
   !--------------------------------------------------------------------

   esmfgrid = ESMF_GridCreate(             &
         name=gridname,                 &
         countsPerDEDim1=IMS,           &
         countsPerDEDim2=JMS,           &
         indexFlag = ESMF_INDEX_DELOCAL,&
         gridEdgeLWidth = (/0,0/),      &
         gridEdgeUWidth = (/0,0/),      &
         coordDep1 = (/1,2/),           &
         coordDep2 = (/1,2/),           &
         rc=status)
   _VERIFY(STATUS)

   ! Allocate coords at default stagger location
   call ESMF_GridAddCoord(esmfgrid, rc=status)
   _VERIFY(STATUS)

   call ESMF_AttributeSet(esmfgrid, name='GRID_LM', value=LM, rc=status)
   _VERIFY(STATUS)

   if (grid_type <= 3) then
      call ESMF_AttributeSet(esmfgrid, 'GridType', 'Cubed-Sphere', rc=STATUS)
      _VERIFY(STATUS)
   else
      call ESMF_AttributeSet(esmfgrid, 'GridType', 'Doubly-Periodic', rc=STATUS)
      _VERIFY(STATUS)
   endif

   ! -------------

   NPES_X = size(IMS)
   NPES_Y = size(JMS)
   NPES = NPES_X+NPES_Y

 call MAPL_GRID_INTERIOR(esmfgrid,isg,ieg,jsg,jeg)

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
      STATUS=1
   endif
   _VERIFY(STATUS)

   if (.not.allocated(FV_Atm)) then 

      AppGridStandAlone = .true.
      if(MAPL_ShmInitialized) then
         call MAPL_AllocNodeArray(grid_global,Shp=(/npts+1,npts+1,ndims,ntiles/),rc=status)
         _VERIFY(STATUS)
      else
         allocate( grid_global(npts+1,npts+1,ndims,ntiles) )
      end if

      if (grid_type <= 3) then
         ! Cubed-Sphere
         call gnomonic_grids(grid_type, npts, grid_global(:,:,1,1), grid_global(:,:,2,1))
         ! mirror_grid assumes that the tile=1 is centered on equator and greenwich meridian Lon[-pi,pi]
         call mirror_grid(grid_global, 0, npts+1, npts+1, 2, 6)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! ! MAT BMA You cannot refer to FV_Atm here because line 288 says this                  !
         ! !         is only for FV_Atm not allocated                                            !
         !                                                                                       !
         ! if (.not. FV_Atm(1)%flagstruct%do_schmidt) &                                          !
         !       grid_global(:,:,1,:) = grid_global(:,:,1,:) - PI/FV_Atm(1)%flagstruct%shift_fac !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         grid_global(:,:,1,:) = grid_global(:,:,1,:) - PI/18.0

         where(grid_global(:,:,1,:) < 0.) &
               grid_global(:,:,1,:) = grid_global(:,:,1,:) + 2.* PI

         ! Keep Equator and Greenwich exact
         !---------------------------------

         where(abs(grid_global(:,:,:,1)) < 1.d-10) grid_global(:,:,:,1) = 0.0

         !---------------------------------
         ! Clean Up Corners
         !---------------------------------
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

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! ! MAT BMA You cannot refer to FV_Atm here because line 288 says this                                                                                    !
         ! !         is only for FV_Atm not allocated                                                                                                              !
         !                                                                                                                                                         !
         ! !------------------------                                                                                                                               !
         ! ! Schmidt transformation:                                                                                                                               !
         ! !------------------------                                                                                                                               !
         ! if ( FV_Atm(1)%flagstruct%do_schmidt ) then                                                                                                             !
         !    do n=1,ntiles                                                                                                                                        !
         !       call direct_transform(FV_Atm(1)%flagstruct%stretch_fac, 1, npts+1, 1, npts+1, FV_Atm(1)%flagstruct%target_lon, FV_Atm(1)%flagstruct%target_lat, & !
         !             n, grid_global(1:npts+1,1:npts+1,1,n), grid_global(1:npts+1,1:npts+1,2,n))                                                                  !
         !    enddo                                                                                                                                                !
         ! endif                                                                                                                                                   !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      else

         if (MAPL_AM_I_ROOT()) write(*,*)'AppGridCreate can not make DP grid by itself'
         _ASSERT(.false.,'needs informative message')

      end if

   end if

   ! Fill lat/lons at cell center locations
   ! Retrieve the coordinates so we can set them
   call ESMF_GridGetCoord(esmfgrid, coordDim=1, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         farrayPtr=lons, rc=status)
   _VERIFY(STATUS)

   call ESMF_GridGetCoord(esmfgrid, coordDim=2, localDE=0, &
         staggerloc=ESMF_STAGGERLOC_CENTER, &
         farrayPtr=lats, rc=status)
   _VERIFY(STATUS)

   ! Save corners into the esmfgrid
   !---------------------------
   allocate(gridCornerLons((size(lons,1)+1)*(size(lons,2)+1)), stat=status)
   _VERIFY(STATUS)
   allocate(gridCornerLats((size(lats,1)+1)*(size(lats,2)+1)), stat=status)
   _VERIFY(STATUS)

   idx = 0
   do jg=jsg,jeg+1
      do ig=isg,ieg+1
         i=ig
         j=jg-myTile*npts
         idx = idx + 1
         if (AppGridStandAlone) then
            gridCornerLons(idx) = grid_global(i, j, 1, myTile+1)
            gridCornerLats(idx) = grid_global(i, j, 2, myTile+1)
         else
            gridCornerLons(idx) = FV_Atm(1)%gridstruct%grid(i,j,1) 
            gridCornerLats(idx) = FV_Atm(1)%gridstruct%grid(i,j,2)
         end if
      end do
   end do

   call ESMF_AttributeSet(esmfgrid, name='GridCornerLons:', &
         itemCount = idx, valueList=gridCornerLons, rc=status)
   _VERIFY(STATUS)
   call ESMF_AttributeSet(esmfgrid, name='GridCornerLats:', &
         itemCount = idx, valueList=gridCornerLats, rc=status)
   _VERIFY(STATUS)

   deallocate(gridCornerLats, stat=status)
   _VERIFY(STATUS)
   deallocate(gridCornerLons, stat=status)
   _VERIFY(STATUS)

   do jg=jsg,jeg
      do ig=isg,ieg
         i=ig
         j=jg-myTile*npts
         if (AppGridStandAlone) then
            call cell_center2(grid_global(i,j,  1:2,myTile+1), grid_global(i+1,j,  1:2,myTile+1),   &
                  grid_global(i,j+1,1:2,myTile+1), grid_global(i+1,j+1,1:2,myTile+1),   &
                  alocs)
         else
            call cell_center2(dble(FV_Atm(1)%gridstruct%grid(i,j,  1:2)),     &
                  dble(FV_Atm(1)%gridstruct%grid(i+1,j,  1:2)),   &
                  dble(FV_Atm(1)%gridstruct%grid(i,j+1,1:2)),     &
                  dble(FV_Atm(1)%gridstruct%grid(i+1,j+1,1:2)),   &
                  alocs)
         end if
         i=ig-isg+1
         j=jg-jsg+1
         lons(i,j) = alocs(1)
         lats(i,j) = alocs(2)
      enddo
   enddo

   if (AppGridStandAlone) then
      DEALLOCGLOB_(grid_global)
   end if
   call MAPL_MemUtilsWrite(VM, trim(Iam), RC=STATUS )
   _VERIFY(STATUS)

   _RETURN(ESMF_SUCCESS)
end function AppGridCreateF

subroutine AppGridCreate (META, esmfgrid, RC)

#include "MAPL_Generic.h"

   use ESMF
   use MAPL_Mod
   use MAPL_BaseMod
   use MAPL_GenericMod
   use fv_arrays_mod,     only: REAL4, REAL8, R_GRID
   implicit none

   ! !ARGUMENTS:

   type(MAPL_MetaComp), intent(INOUT) :: META
   type (ESMF_Grid),    intent(  OUT) :: esmfgrid
   integer, optional,   intent(  OUT) :: rc

   ! ErrLog variables
   !-----------------

   integer                      :: STATUS
   character(len=ESMF_MAXSTR), parameter :: Iam="AppGridCreate"

   ! Local variables
   !-----------------

   integer                         :: IM_WORLD
   integer                         :: JM_WORLD
   integer                         :: LM
   integer                         :: NX, NY
   character(len=ESMF_MAXSTR)      :: tmpname, gridname
   character(len=ESMF_MAXSTR)      :: FMT, FMTIM, FMTJM
   real(REAL8) :: deglon=0.0

   interface
      function AppGridCreateF(IM_WORLD, JM_WORLD, LM, NX, NY, rc)
         use ESMF
         implicit none
         type(ESMF_Grid)                  :: AppGridCreateF
         integer,           intent(IN)    :: IM_WORLD, JM_WORLD, LM
         integer,           intent(IN)    :: NX, NY
         integer, optional, intent(OUT)   :: rc
      end function AppGridCreateF
   end interface

   ! Get Decomposition from CF
   !--------------------------

   ! !RESOURCE_ITEM: none :: Processing elements in 1st dimension
   call MAPL_GetResource( META, NX,       label ='NX:', default=1, rc = status )
   _VERIFY(STATUS)
   ! !RESOURCE_ITEM: none :: Processing elements in 2nd dimension
   call MAPL_GetResource( META, NY,       label ='NY:', default=1, rc = status )
   _VERIFY(STATUS)

   ! Get World problem size from CF
   !-------------------------------

   ! !RESOURCE_ITEM: none :: Grid size in 1st dimension
   call MAPL_GetResource( META, IM_WORLD, 'AGCM_IM:',            rc = status )
   _VERIFY(STATUS)
   ! !RESOURCE_ITEM: none :: Grid size in 2nd dimension
   call MAPL_GetResource( META, JM_WORLD, 'AGCM_JM:',            rc = status )
   _VERIFY(STATUS)
   ! !RESOURCE_ITEM: none :: Grid size in 3rd dimension
   call MAPL_GetResource( META, LM,       'AGCM_LM:', default=1, rc = status )
   _VERIFY(STATUS)

   ! The grid's name is optional
   !----------------------------

   call GET_INT_FORMAT(IM_WORLD, FMTIM)
   call GET_INT_FORMAT(JM_WORLD, FMTJM)
   FMT = '(A,' // trim(FMTIM) //',A,' // trim(FMTJM) // ',A)'
   write(tmpname,trim(FMT)) 'PE',IM_WORLD,'x',JM_WORLD,'-CF'
   ! !RESOURCE_ITEM: none :: Optional grid name
   call MAPL_GetResource( META, GRIDNAME, 'AGCM_GRIDNAME:', default=trim(tmpname), rc = status )
   _VERIFY(STATUS)

   esmfgrid = AppGridCreateF(IM_WORLD, JM_WORLD, LM, NX, NY, STATUS)
   _VERIFY(STATUS)

   _RETURN(ESMF_SUCCESS)

end subroutine AppGridCreate

subroutine GET_INT_FORMAT(N, FMT)
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
end subroutine GET_INT_FORMAT

subroutine DecomposeDim ( dim_world,dim,NDEs )

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

end subroutine DecomposeDim


