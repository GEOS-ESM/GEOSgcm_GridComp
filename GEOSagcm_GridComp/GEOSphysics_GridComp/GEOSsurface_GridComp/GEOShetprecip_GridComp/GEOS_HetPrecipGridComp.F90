
#include "MAPL_Generic.h"

!=============================================================================
module GEOS_HetPrecipGridCompMod
  use ESMF
  use MAPL
  use random
  implicit none

  public SetServices
  ! some module variables
  integer, parameter :: NPDF=28
  integer, parameter :: SECONDS_PER_DAY = 24*3600 ! should be 86400

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
    real :: DT, TAU_HETPR
    integer :: IM, JM, NT, I
    real, pointer :: qvar(:) => null()
    real :: totalArea
    integer, pointer :: tiletypes(:) => null()
    integer, allocatable :: the_seed(:)
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
    _ASSERT(IM==1 .and. JM == 1, 'Only single column supported')

    call MAPL_LocStreamGet(LocStream, NT_LOCAL=NT, __RC__)

    ! get DT
    call MAPL_GetResource(MAPL, DT, Label="RUN_DT:", __RC__)

    ! get autocorrelation timescale
    call MAPL_GetResource(MAPL, TAU_HETPR, Label="HETPRECIP_TAU:", default=3600., __RC__)

    RHO = (1.0/E) **(1.0*DT/TAU_HETPR)
    SQRHO = SQRT(1.0-RHO**2)

! pdf / cdf
!ALT: This is direct port from Randy's idl code
    ! piecewise fit to pdf of weights:
    pdfw=[ 0.00793, 0.01790, 0.02253, 0.02837, 0.03572, 0.04496, &
           0.05661, 0.07126, 0.08972, 0.11295, 0.14219, 0.17901, &
           0.22536, 0.28371, 0.35717, 0.44965, 0.56607, 0.71264, &
           0.89712, 1.12946, 1.42191, 1.79008, 2.25357, 2.83708, &
           3.57167, 4.49647, 5.66072, 7.12643, 8.97164 ]
    pdfv=[ 0.04172, 0.01285, 0.01469, 0.01668, 0.02590, 0.02158, &   ! 3-degree to 16km
           0.02514, 0.02773, 0.03141, 0.03369, 0.03649, 0.04119, &
           0.04404, 0.04917, 0.04951, 0.05514, 0.05683, 0.05903, &
           0.05865, 0.05544, 0.05261, 0.04915, 0.04156, 0.03439, &
           0.02713, 0.01960, 0.01271, 0.00767, 0.00349 ]
!    pdfv=[ 0.00956, 0.00351, 0.00499, 0.00618, 0.00776, 0.00969, &  ! 0.5-degree to 16km
!           0.01183, 0.01370, 0.01745, 0.02012, 0.02375, 0.02761, &
!           0.03253, 0.03826, 0.04460, 0.05401, 0.06462, 0.07981, &
!           0.09678, 0.13621, 0.09749, 0.08644, 0.05649, 0.03426, &
!           0.01610, 0.00517, 0.00105, 0.00004, 0.0 ]
!    pdfw=[ 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, &
!           2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9, &
!           4.1, 4.3, 4.5, 4.7, 4.9] ! weights
!    pdfv=[.121,.128,.137,.139,.119, .081,.059,.041,.030,.025, &   ! original
!          .019,.014,.012,.009,.007, .005,.003,.002,.001,.001, &
!          .0008,.0006,.0004,.0002,.0001] ! probability density
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

    call random_seed(SIZE=seed_len)
    allocate(THE_SEED(seed_len))

    THE_SEED(1)=1
    THE_SEED(2)=5

    ! Gfortran uses longer seeds, so fill the rest with zero
    if (seed_len > 2) THE_SEED(3:) = 0

    call random_seed(PUT=THE_SEED)


    if (all(qvar == 0.0)) then
       call RANDOM_NORMAL(qvar)
    end if

    deallocate(THE_SEED)

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
    real, parameter :: xlo=log10(5e-2), xhi=log10(10.)
    real, parameter :: ylo=0.9, yhi=0.0, fracdrymax=0.95
    real :: fracdry, yinterp
    real :: pperday
    integer :: NT
    integer :: i, iopt, n, itile
    integer :: ilargest
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

    call random_normal(rn) ! rn is array of size NT

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

    pperday=log10(totalPrecip*SECONDS_PER_DAY)

    yinterp=ylo+((pperday-xlo)/(xhi-xlo))*(yhi-ylo)
    if (pperday > xhi) yinterp = yhi
    fracdry=yinterp
    if(fracdry > fracdrymax) then 
       fracdry=fracdrymax
    end if

    allocate(psum(NT), krank(NT), __STAT__)

    ! sort QVAR
    call rankdata(nt, qvar, krank, __RC__)

    ! compute a partial sum of tile areas in the ranked order
    do n = 1,NT
       itile = krank(n)
       if (n == 1) then
          psum(itile) = norm_tile_area(itile)
       else
          psum(itile) = psum(krank(n-1)) + norm_tile_area(itile)
       end if
    end do

    ! bring the partial sum to the "middle" of the last ranked tile
!    do n = 1,NT
!       itile = krank(n)
!       psum(itile) =  psum(itile) - 0.5*norm_tile_area(itile)
!    end do

    ! Note psum is the same variable as qranked in Randy's IDL code

    ! use that and the cdfw/v to calculate heterogeneneous factor

    do n = 1,nt
       itile = krank(n)
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
