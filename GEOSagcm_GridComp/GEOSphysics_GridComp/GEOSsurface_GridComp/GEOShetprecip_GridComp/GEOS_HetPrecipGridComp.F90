
#include "MAPL_Generic.h"

!=============================================================================
module GEOS_HetPrecipGridCompMod
  use ESMF
  use MAPL
  use random
  implicit none

  public SetServices
  ! some module variables
  integer, parameter :: SECONDS_PER_HOUR = 3600

  real :: rho ! autocorrelation

  real :: sqrho ! = sqrt(1-rho**2) autocorrelation

  real, pointer :: norm_tile_area(:) => null()

  type MyTileIndex
     integer, pointer :: tileId(:) => null()
  end type MyTileIndex

  type (MyTileIndex), allocatable :: tileIndex(:,:)
  
contains
  subroutine SetServices ( GC, RC )

    type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
    integer,             intent(  OUT) :: RC  ! return code

    integer :: status

    character(len=ESMF_MAXSTR) :: COMP_NAME

    call ESMF_GridCompGet( GC, NAME=COMP_NAME, __RC__ )


     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'liquid_water_convective_precipitation', &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'PCU',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       __RC__ )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'liquid_water_large_scale_precipitation', &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'PLS',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       __RC__ )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'snowfall',                          &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'SNO',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       RC=STATUS  )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'icefall',                           &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'ICE',                               &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       __RC__ )

     call MAPL_AddImportSpec(GC,                             &
        LONG_NAME          = 'freezing_rain_fall',                &
        UNITS              = 'kg m-2 s-1',                        &
        SHORT_NAME         = 'FRZR',                              &
        DIMS               = MAPL_DimsHorzOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       __RC__ )

     call MAPL_AddInternalSpec(GC,                             &
        LONG_NAME          = 'prognostic_memory_state_var',                &
        UNITS              = 'N/A',                        &
        SHORT_NAME         = 'QVAR',                              &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       __RC__ )

     call MAPL_AddExportSpec(GC,                             &
        LONG_NAME          = 'weighting_factor_to_scale_precip',  &
        UNITS              = 'N/A',                        &
        SHORT_NAME         = 'WEIGHT',                              &
        DIMS               = MAPL_DimsTileOnly,                   &
        VLOCATION          = MAPL_VLocationNone,                  &
                                                       __RC__ )

! Set the Run entry point
! -----------------------

    call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_INITIALIZE, Initialize, __RC__)

    call MAPL_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN,  Run, __RC__)

    call MAPL_GenericSetServices    ( GC, __RC__ )
 
    _RETURN(ESMF_SUCCESS)

  end subroutine SetServices

  subroutine Initialize ( GC, IMPORT, EXPORT, CLOCK, RC )

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

    integer :: status
    real :: E ! this is Eulers constant, about 2.718281828

    type (MAPL_MetaComp), pointer :: MAPL
    type (ESMF_State) :: INTERNAL
    type (MAPL_LocStream) :: LocStream
    real :: DT, TAU_HETPR, SEED_OPT
    integer :: IM, JM, NT, I, J
    real, pointer :: qvar(:) => null()
    real :: totalArea
    integer, pointer :: tiletypes(:) => null()
    integer, allocatable :: the_seed(:)
    integer, allocatable :: tcount(:,:)
    integer, pointer :: local_i(:) => null()
    integer, pointer :: local_j(:) => null()
    integer :: it
    integer :: seed_len
    
    call MAPL_GenericInitialize(GC, IMPORT, EXPORT, CLOCK, __RC__)

    E = EXP(1.0)

    ! get MAPL
    call MAPL_GetObjectFromGC ( GC, MAPL, __RC__)

    ! get from MAPL: IM, JM, LocStream
    ! note that this locsteam in covers all surface types (land included)
    ! we could modify surface to pass a different LocStream, but that 
    ! would require more changes in Surface

    call MAPL_Get(MAPL, IM=IM, JM=JM, ExchangeGrid=LocStream, __RC__)
    call MAPL_Get(MAPL, &
         INTERNAL_ESMF_STATE=INTERNAL, &
         TILEAREA  = norm_tile_area,   &
         TILETYPES = TILETYPES,        &
                                __RC__ )

!ALT: For now restrict usage only for single column mode
!@@    _ASSERT(IM==1 .and. JM == 1, 'Only single column supported')

    call MAPL_LocStreamGet(LocStream, NT_LOCAL=NT, __RC__)

    ! get DT
    call MAPL_GetResource(MAPL, DT, Label="RUN_DT:", __RC__)

    ! get autocorrelation timescale
    call MAPL_GetResource(MAPL, TAU_HETPR, Label="HETPRECIP_TAU:", default=3600., __RC__)

    call MAPL_GetResource(MAPL, SEED_OPT, Label="HETPRECIP_SEED:", default=5., __RC__)

    RHO = (1.0/E) **(1.0*DT/TAU_HETPR)
    SQRHO = SQRT(1.0-RHO**2)


    !get local_i, local_j
    call MAPL_LocStreamGet(LocStream, NT_LOCAL=NT, __RC__)
    allocate(tcount(im,jm), _STAT) ! integer
    tcount = 0
    !loop over tiles
    do it=1,NT
       tcount(local_i(it), local_j(it)) = tcount(local_i(it), local_j(it)) + 1
    enddo

!@    allocate(mystate%tileIndex(im,jm), _STAT)
    allocate(tileIndex(im,jm), _STAT)
    do j=1,jm
       do i=1,im
!@          allocate(mystate%tileIndex(i,j)%tileId(tcount(i,j)), _STAT)
          allocate(tileIndex(i,j)%tileId(tcount(i,j)), _STAT)
       enddo
    enddo
    tcount = 0
    do it=1,NT
       i = local_i(it)
       j = local_j(it)
       tcount(i,j) = tcount(i,j) + 1
!@       mystate%tileIndex(i,j)%tileId(tcount(i,j)) = it
       tileIndex(i,j)%tileId(tcount(i,j)) = it
    enddo

    deallocate(tcount)
#ifdef __GFORTRAN__
    deallocate(local_i, local_j)
#endif


    ! get qvar

!@@@@ hetprecip-ovel-land restriction
    _ASSERT(all(tiletypes == MAPL_Land), 'Currently supporting cells exclusively covered by land')
    !@@@@ single column assumption (next 2 lines)
    do j=1,jm
       do i=1,im
!@          call normalizeTileArea(norm_tile_area(mystate%tileIndex(i,j)%tileId))
          call normalizeTileArea(norm_tile_area(tileIndex(i,j)%tileId))
       enddo
    enddo

    call MAPL_GetPointer(INTERNAL, QVAR, 'QVAR', __RC__)

    call random_seed(SIZE=seed_len)
    allocate(THE_SEED(seed_len))

    THE_SEED(1) = INT(SEED_OPT)
    THE_SEED(2) = 2*INT(SEED_OPT)

    ! Gfortran uses longer seeds, so fill the rest with zero
    if (seed_len > 2) THE_SEED(3:) = 0

    call random_seed(PUT=THE_SEED)

    if (all(qvar == 0.0)) then
       call RANDOM_NORMAL(qvar)
    end if

    deallocate(THE_SEED)

    _RETURN(ESMF_SUCCESS)

  contains
    subroutine normalizeTileArea(norm_tile_area)
      real  :: norm_tile_area(:)

      real :: totalArea
      totalArea = sum(norm_tile_area)
      norm_tile_area = norm_tile_area / totalArea

      return
    end subroutine normalizeTileArea

  end subroutine Initialize

  subroutine Run ( GC, IMPORT, EXPORT, CLOCK, RC )

    type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type(ESMF_State),    intent(inout) :: IMPORT ! Import state
    type(ESMF_State),    intent(inout) :: EXPORT ! Export state
    type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,   intent(  out) :: RC     ! Error code

    integer :: status

    type (MAPL_MetaComp), pointer :: MAPL
    real, pointer :: qvar(:) => null()
    real, pointer :: pcu(:,:) => null()
    real, pointer :: pls(:,:) => null()
    real, pointer :: sno(:,:) => null()
    real, pointer :: ice(:,:) => null()
    real, pointer :: frzr(:,:) => null()
    real, pointer :: psub(:) => null()
    type (ESMF_State) :: INTERNAL
    integer :: i, j, im, jm
    integer, pointer :: tileIdx(:)
    real :: totalPrecip
    

! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, __RC__)

    ! get qvar
    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, IM=IM, JM=JM, __RC__)
    call MAPL_GetPointer(INTERNAL, QVAR, 'QVAR', __RC__)


    call MAPL_GetPointer(IMPORT, PCU, 'PCU', __RC__)
    call MAPL_GetPointer(IMPORT, PLS, 'PLS', __RC__)
    call MAPL_GetPointer(IMPORT, SNO, 'SNO', __RC__)
    call MAPL_GetPointer(IMPORT, ICE, 'ICE', __RC__)
    call MAPL_GetPointer(IMPORT, FRZR, 'FRZR', __RC__)

    !ALT: the alloc-.true. is there from laziness
    call MAPL_GetPointer(EXPORT, psub, 'WEIGHT', alloc=.true., __RC__)

    DO J=1,JM
       DO I=1,IM

          totalPrecip = PCU(I,J)+PLS(I,J)+SNO(I,J)+ICE(I,J)+FRZR(I,J)
!@          tileIdx => mystate%tileIndex(i,j)%tileId
          tileIdx => tileIndex(i,j)%tileId
          
          call generateHetWeight(qvar(tileIdx), totalPrecip, psub(tileIdx), _RC)
       END DO
    END DO
    _RETURN(ESMF_SUCCESS)


  contains
     subroutine rankdata (nt, qtest, krank, rc)
       integer, intent(in) :: nt
       real, intent(in) :: qtest(:)
       integer, intent(out) :: krank(:)
       integer, optional, intent(out) :: rc
       
       integer :: status
       integer, allocatable :: ptest(:)
       integer, parameter :: large_int = 2**21
       
       allocate(ptest(nt), __STAT__)
       ptest=int(qtest*large_int)
       krank = [1:nt]

       call MAPL_Sort(ptest, krank)

       deallocate(ptest)

       _RETURN(ESMF_SUCCESS)
     end subroutine rankdata

    subroutine generateHetWeight(qvar, totalPrecip, psub, rc)
      real :: qvar(:)
      real, intent(in) :: totalPrecip
      real :: psub(:)
      integer, optional, intent(out) :: rc
      
    real, allocatable :: rn(:)
    real :: total
    real, parameter :: xlo=log10(1.5e-2), xhi=log10(10.) ! parameters for 1 deg domain
    real, parameter :: ylo=.973, yhi=0.0, fracdrymax=0.9722
    real, parameter :: a1opt=0.65, a2opt=0.96 ! for relationship between CDF and precip. scaling
    real, parameter :: piby2=3.14159/2.
    real :: psum0
    real :: alow, ahigh, aintegral1, aintegral2, aintegral
    real :: fracdry, yinterp
    real :: pperday
    integer :: NT
    integer :: i, iopt, n, itile
    integer :: ilargest
    integer, allocatable :: krank(:)
    real, allocatable :: psum(:)
    real, allocatable :: lowerlimit(:)
    real, allocatable :: upperlimit(:)
    integer :: status



    NT = size(qvar)
    allocate(rn(NT), __STAT__)

    call random_normal(rn) ! rn is array of size NT

!   generate new QVAR
    qvar=rho*qvar+sqrho*rn

    if (totalPrecip == 0.0) then
       print *,'WARNING: zero precip.'
       psub = 1.0
       _RETURN(ESMF_SUCCESS)
    end if


! determine fracdry

    pperday=log10(totalPrecip*SECONDS_PER_HOUR)

    yinterp=ylo+((pperday-xlo)/(xhi-xlo))*(yhi-ylo)
    if (pperday > xhi) yinterp = yhi
    fracdry=yinterp
    if(fracdry > fracdrymax) then 
       fracdry=fracdrymax
    end if
!    fracdry = 0.   !!! Temporary test

    allocate(psum(NT), krank(NT), __STAT__)
    allocate(lowerlimit(NT),upperlimit(NT), __STAT__ )

    ! sort QVAR
    call rankdata(nt, qvar, krank, __RC__)
!    print *,'krank=',krank

    psum0=0.
    ! compute a partial sum of tile areas in the ranked order;
    ! also set lower and upper limits of integral calculation
    do n = 1,NT
       itile = krank(n)
       lowerlimit(itile)=psum0
       psum0 = psum0 + norm_tile_area(itile)
       upperlimit(itile)=psum0      
    end do


    ! a1opt*(tan(a2opt x)) is a fit to the CDF - precip. scaling factor 
    ! relationship. Compute integral of a1opt*(tan(a2opt x)) analytically

    do n = 1,nt
       itile = krank(n)
       if(upperlimit(itile) < fracdry) then
          psub(itile)=0.
       else
          alow=lowerlimit(itile)
          if(lowerlimit(itile) < fracdry) alow=fracdry
          ahigh=upperlimit(itile)

          alow=(alow-fracdry)/(1.-fracdry)
          ahigh=(ahigh-fracdry)/(1.-fracdry)

          aintegral1=log(abs(1./cos(a2opt*ahigh*piby2))) ! natural logarithm
          aintegral2=log(abs(1./cos(a2opt*alow*piby2)))
          aintegral=(a1opt/(a2opt*piby2))*(aintegral1-aintegral2)
          psub(itile)=aintegral/(norm_tile_area(itile)*(1.-fracdry))

          iopt=-1
       end if
     end do

     ! add a protection for the case of non-zero grid precip but the scheme
     ! resulted in zero precip for all the tiles. Then we assign all precip to
     ! the tile with the largest qvar
     total = 0.0
     ilargest = -1
     do i=1,nt
        total = total + psub(i)
        ilargest = max(krank(i), ilargest)
     end do
     _ASSERT(ilargest > 0, 'Could not find largest QVAR')
     if(total == 0.0) psub(ilargest) = 1.0

     ! empirical rescaling to reduce time-mean precip bias
!     print *,'factor=',(4.2*EXP(-25.*norm_tile_area)+0.6)
!     psub = psub / (4.2*EXP(-25.*norm_tile_area)+0.6)

     ! scale the precip factor to make sure we preserve the grid box precip
!ALT: old scaling, assumes uniform area tiles     total = sum(psub)/nt
     total = 0.0
     do i=1,nt
        total = total + psub(i)*norm_tile_area(i)
     end do
     psub = psub/total


! all done
     deallocate(psum,krank, rn) 

    _RETURN(ESMF_SUCCESS)

   end subroutine generateHetWeight
  end subroutine Run

end module GEOS_HetPrecipGridCompMod

subroutine SetServices(gc, rc)
   use ESMF
   use GEOS_HetPrecipGridCompMod, only : mySetservices=>SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc,rc=rc)
end subroutine
