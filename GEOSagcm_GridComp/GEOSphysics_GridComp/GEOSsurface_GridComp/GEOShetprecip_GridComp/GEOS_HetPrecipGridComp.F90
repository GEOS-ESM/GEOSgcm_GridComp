
#include "MAPL_Generic.h"

!=============================================================================
module GEOS_HetPrecipGridCompMod
  use ESMF
  use MAPL
  implicit none

  public SetServices
  ! some module variables
  integer, parameter :: NPDF=24
  integer :: STEPS_PER_HOUR

  real :: rho ! autocorrelation
  real :: sqrho ! = sqrt(1-rho**2) autocorrelation
  real :: pdfw(0:NPDF)
  real :: pdfv(0:NPDF)
  real :: cdfw(0:NPDF)
  real :: cdfv(0:NPDF)
  real, pointer :: norm_tile_area(:) => null()

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
    integer :: RUN_DT
    real :: DT
    integer :: IM, JM, NT, I
    real, pointer :: qvar(:) => null()
    real :: totalArea
    integer, pointer :: tiletypes(:) => null()

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
    _ASSERT(IM==1 .and. JM == 1, 'Only single column supported')

    call MAPL_LocStreamGet(LocStream, NT_LOCAL=NT, __RC__)

    ! get DT
    call MAPL_GetResource(MAPL, DT, Label="RUN_DT:", __RC__)

    RUN_DT = nint(DT)
    ! compute number of steps in an hour
    STEPS_PER_HOUR = 3600/RUN_DT

    ! autocorrelation for 1 hour
    RHO = (1.0/E) **(1.0/STEPS_PER_HOUR)
    SQRHO = SQRT(1.0-RHO**2)

! pdf / cdf
!ALT: This is direct port from Randy's idl code
    ! piecewise fit to pdf of weights:
    pdfw=[ 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, &
           2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, &
           4.1, 4.3, 4.5, 4.7, 4.9] ! weights
    pdfv=[.121,.128,.137,.139,.119, .081,.059,.041,.030,.025, &
          .019,.014,.012,.009,.007, .005,.003,.002,.001,.001, &
          .0008,.0006,.0004,.0002,.0001] ! probability density
    pdfv=pdfv/sum(pdfv) ! scale to ensure it sums to unity
    cdfv(0)=pdfv(0)
    do i=1,NPDF
       cdfv(i)=cdfv(i-1)+pdfv(i) ! construct CDF
    end do
    cdfv=cdfv/cdfv(NPDF) ! final scaling, just to be sure

    ! get qvar

    _ASSERT(all(tiletypes == MAPL_Land), 'Currently supporting cells exclusively covered by land')
    totalArea = sum(norm_tile_area)
    norm_tile_area = norm_tile_area / totalArea

    call MAPL_GetPointer(INTERNAL, QVAR, 'QVAR', __RC__)

    ! optioanally deal with seeding the random number generator
    ! ...TBD...

    if (all(qvar == 0.0)) then
       call RANDOM_NUMBER(qvar)
    end if

    _RETURN(ESMF_SUCCESS)
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
    type (ESMF_State) :: INTERNAL

    real, allocatable :: rn(:)
    real :: totalPrecip
    real :: total
    real, parameter :: xlo=log10(0.07), xhi=log10(20.)
    real, parameter :: ylo=0.9, yhi=0.0, fracdrymax=0.95
    real :: fracdry, yinterp
    real :: pperhr
    integer :: NT
    integer :: i, iopt, n, itile
    integer, allocatable :: krank(:)
    real, pointer :: psub(:) => null()
    real, allocatable :: psum(:)
    real :: qfrac, wtx, wtpdf


! Get my internal MAPL_Generic state
!-----------------------------------

    call MAPL_GetObjectFromGC ( GC, MAPL, __RC__)

    ! get qvar
    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL, __RC__)
    call MAPL_GetPointer(INTERNAL, QVAR, 'QVAR', __RC__)

    NT = size(QVAR)
    allocate(rn(NT), __STAT__)

    call random_number(rn) ! rn is array of size NT

!   generate new QVAR
    qvar=rho*qvar+sqrho*rn

    call MAPL_GetPointer(IMPORT, PCU, 'PCU', __RC__)
    call MAPL_GetPointer(IMPORT, PLS, 'PLS', __RC__)
    call MAPL_GetPointer(IMPORT, SNO, 'SNO', __RC__)
    call MAPL_GetPointer(IMPORT, ICE, 'ICE', __RC__)
    call MAPL_GetPointer(IMPORT, FRZR, 'FRZR', __RC__)

    !ALT: the alloc-.true. is there from laziness
    call MAPL_GetPointer(EXPORT, psub, 'WEIGHT', alloc=.true., __RC__)

!ALT: single column assumption: 1 grid cell
    totalPrecip = PCU(1,1)+PLS(1,1)+SNO(1,1)+ICE(1,1)+FRZR(1,1)
    if (totalPrecip == 0.0) then
       print *,'WARNING: zero precip.'
       psub = 1.0
       _RETURN(ESMF_SUCCESS)
    end if


! determine fracdry

    pperhr=log10(totalPrecip*STEPS_PER_HOUR)

    yinterp=ylo+((pperhr-xlo)/(xhi-xlo))*(yhi-ylo)
    if (pperhr > xhi) yinterp = yhi
    fracdry=yinterp
    if(fracdry > fracdrymax) then 
       fracdry=fracdrymax
    end if

    allocate(psum(0:NT), krank(NT), __STAT__)

    ! sort QVAR
    call rankdata(nt, qvar, krank, __RC__)

    ! compute a partial sum of tile areas in the ranked order
    psum(0) = 0.0
    do n = 1,NT
       itile = krank(n)
       psum(n) =  psum(n-1) + norm_tile_area(itile)
    end do

    ! bring the partial sum to the "middle" of the last ranked tile
    do n = 1,NT
       itile = krank(n)
       psum(n) =  psum(n) - 0.5*norm_tile_area(itile)
    end do

    ! Note psum is the same variable as qranked in Randy's IDL code

    ! use that and the cdfw/v to calculate heterogeneneous factor

    do itile=1,nt
       if(psum(itile) < fracdry) then
          psub(itile)=0.
       else
          qfrac=(psum(itile)-fracdry)/(1.-fracdry)

          iopt=-1
          do i=0,NPDF 
             if(qfrac < cdfv(i)) cycle
             iopt=i ! find location on CDF
             if(iopt < NPDF) then
                wtx=(qfrac-cdfv(iopt))/(cdfv(iopt+1)-cdfv(iopt))
                wtpdf=pdfw(iopt)+wtx*(pdfw(iopt+1)-pdfw(iopt))
                ! actual calculation of the subgrid precipitation factor
                psub(itile)=wtpdf*(1.-fracdry)
             else
                !ALT: this logic need to be double checked
                wtpdf=pdfw(iopt)
                psub(itile)=wtpdf*(1.-fracdry)
             endif
          end do
       end if
     end do

     ! scale the precip factor to make sure we preserve the grid box precip
     total = sum(psub)/nt
     if (total /= 0.0) psub = psub/total

! all done
     deallocate(psum,krank, rn) 

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

  end subroutine Run

end module GEOS_HetPrecipGridCompMod

subroutine SetServices(gc, rc)
   use ESMF
   use GEOS_HetPrecipGridCompMod, only : mySetservices=>SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc,rc=rc)
end subroutine
