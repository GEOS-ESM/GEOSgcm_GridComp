#include "MAPL_Generic.h"

module Gigatraj_UtilsMod
   use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_null_ptr, c_associated, c_null_char
   use, intrinsic :: iso_c_binding, only : c_loc

   use ESMF
   use MAPL
   use mpi
   implicit none
   public :: parseCompsAndFieldsName
   public :: create_new_vars
   public :: get_levels
   public :: get_latlon_centers
   public :: get_cube_centers
   public :: horde
   public :: GigaTrajInternal
   public :: GigatrajInternalWrap

   type horde
     integer :: num_parcels
     integer, allocatable :: IDS(:)
     real, allocatable    :: lats(:), lons(:), zs(:)
   end type

   type GigaTrajInternal
     integer :: npes
     integer :: npz ! number of pressure levels
     type (ESMF_Grid) :: LatLonGrid
     type (ESMF_Grid) :: CubedGrid
     class (AbstractRegridder), pointer :: cube2latlon => null()
     integer, allocatable :: CellToRank(:,:)
     type(horde) :: parcels
     type(c_ptr) :: metSrc
     type(ESMF_Time) :: startTime
     type(ESMF_TimeInterval)  :: Integrate_DT
     character(len=ESMF_MAXSTR), allocatable :: ExtraFieldNames(:)
     character(len=ESMF_MAXSTR), allocatable :: ExtraCompNames(:)
     character(len=ESMF_MAXSTR), allocatable :: ExtrabundleNames(:)
     character(len=ESMF_MAXSTR), allocatable :: ExtraAliasNames(:)

     character(len=:), allocatable :: vCoord
     character(len=:), allocatable :: vAlias
     character(len=:), allocatable :: vTendency

     logical :: regrid_to_latlon
   end type

   type GigatrajInternalWrap
     type (GigaTrajInternal), pointer :: PTR
   end type

contains

   subroutine parseCompsAndFieldsName(fields_line, CompNames, BundleNames, FieldNames, AliasNames)
      character(*), intent(in) :: fields_line
      character(len=ESMF_MAXSTR), allocatable, intent(out) :: CompNames(:)
      character(len=ESMF_MAXSTR), allocatable, intent(out) :: BundleNames(:)
      character(len=ESMF_MAXSTR), allocatable, intent(out) :: FieldNames(:)
      character(len=ESMF_MAXSTR), allocatable, intent(out) :: AliasNames(:)
      integer :: num_field, i, j, k, l, endl, num_
      character(len=:), allocatable :: tmp, tmp_bnf, tmp_f, tmp_alias
      num_field = 1
      k = 1
      do
        i = index(fields_line(k:),';')
        if (i == 0) exit
        if (trim(fields_line(i+1:)) =='') exit ! take care of the last unnecessay ";"
        k = k+i
        num_field = num_field+1
      enddo

      allocate(Fieldnames(num_field))
      allocate(Compnames(num_field))
      allocate(BundleNames(num_field))
      allocate(AliasNames(num_field))

      k    = 1
      num_ = 1

      do
        i = index(fields_line(k:),';')
        if (i == 0) then
           endl = len(fields_line)
        else
           endl = (k-1)+i-1
        endif
        tmp = fields_line(k:endl)

        j = index(tmp, '%%')
        if (j == 0) print*, "Wrong format of the comp%%field"
        Compnames(num_) = trim(adjustl(tmp(1:j-1)))
        tmp_bnf  = trim(adjustl(tmp(j+2:)))

       l = index(tmp_bnf, '%')
        if (l /=0) then
           BundleNames(num_) = tmp_bnf(1:l-1)
           tmp_f  = tmp_bnf(l+1:)
        else
           BundleNames(num_) = 'NONE'
           tmp_f  = tmp_bnf
        endif

        ! Aliasing....., Hard coded here
        l = index(tmp_f, '|')
        if (l /=0) then
           FieldNames(num_) = tmp_f(1:l-1)
           tmp_alias  = tmp_f(l+1:)
        else
           FieldNames(num_) = tmp_f
           tmp_alias  = tmp_f
        endif

        AliasNames(num_) = tmp_alias

        num_ = num_ + 1
        k = endl + 2
        if (num_ > num_field) exit
      enddo
   end subroutine parseCompsAndFieldsName

   subroutine create_new_vars(meta, formatter, long_name, short_name, units)
     type(FileMetadata), intent(inout) :: meta
     type(Netcdf4_fileformatter), intent(inout) :: formatter
     character(*), intent(in) :: long_name
     character(*), intent(in) :: short_name
     character(*), intent(in) :: units
     type(Variable) :: var
     character(len=:), allocatable :: var_name
     if (MAPL_AM_I_Root()) then
       if( meta%has_variable(short_name)) return
       var_name = short_name
       var = variable(type=pFIO_REAL32, dimensions='id,time')
       call var%add_attribute('long_name', long_name)
       call var%add_attribute('units', units)
       call var%add_attribute('positive', "up")
       call var%add_attribute('_FillValue', -999.99)
       call var%add_attribute('missing_value', -999.99)
       call meta%add_variable(var_name, var)
       call formatter%add_variable(meta, short_name)
     endif
   end subroutine create_new_vars

   subroutine get_levels(P, func, levels, rc)
     real, dimension(:,:,:), intent(in) :: P
     character(*), intent(in) :: func
     real, dimension(:),    intent(out) :: levels
     integer, optional, intent(out) :: rc
     logical :: positive
     type (ESMF_VM) :: vm
     integer :: comm, lm, status, i, ll
     real :: local_min_val, local_max_val, lev01, levLm, delt
     real, allocatable :: temp(:,:)
     character(:), allocatable :: Iam

     Iam = "get_levels"

     lm = size(P,3)

     call ESMF_VMgetCurrent(vm)
     call ESMF_VMGet(vm, mpiCommunicator = comm, rc = status)
     positive = P(1,1,1) < P(1,1,2)      
     print*, "wjiang:positive", positive
     if (positive) then
        local_min_val = maxval(P(:,:,1))
        print*, "local_min_val:", local_max_val
        call MPI_Allreduce(lev01, local_min_val,1, MPI_FLOAT, MPI_MIN, comm, status)
        temp = P(:,:,lm)
        where(temp >= MAPL_UNDEF) temp = -MAPL_UNDEF
        local_max_val = maxval(temp)
        print*, "local_max_val:", local_max_val
        call MPI_Allreduce(levLm, local_max_val,1, MPI_FLOAT, MPI_MAX, comm, status)
     else
        local_min_val = minval(P(:,:,lm))
        print*, "local_min_val:", local_max_val
        call MPI_Allreduce(levLm, local_min_val,1, MPI_FLOAT, MPI_MIN, comm, status)
        temp = P(:,:,1)
        where(temp >= MAPL_UNDEF) temp = -MAPL_UNDEF
        local_max_val = maxval(temp)
        print*, "local_max_val:", local_max_val
        call MPI_Allreduce(lev01, local_max_val,1, MPI_FLOAT, MPI_MAX, comm, status)
     endif

     ll = size(levels)

     if (trim(func) == 'log') then
        delt = (log(levLm)-log(lev01))/(lm-1)
        levels =[ (exp(log(lev01)+i*delt), i=0, ll-1)]
     else
        delt = (levLm-lev01)/(lm-1)
        levels =[ (lev01 + i*delt, i=0, ll-1)]
     endif

   end subroutine get_levels

   subroutine get_cube_centers(GC, lon_center, lat_center, rc)
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      real, allocatable,     intent(out) :: lat_center(:,:), lon_center(:,:)
      integer, optional,   intent(  out) :: RC  
      integer :: i1, i2, j1, j2, imc, jmc, status 
      real(ESMF_KIND_R8), pointer     :: centerX(:,:)
      real(ESMF_KIND_R8), pointer     :: centerY(:,:)
      real(ESMF_KIND_R8), pointer     :: ptr(:,:)
      type(ESMF_Field)   :: field
      type(ESMF_RouteHandle) :: rh
      type (GigaTrajInternal), pointer :: GigaTrajInternalPtr
      type (GigatrajInternalWrap)   :: wrap
      type(ESMF_Grid) :: grid_
      character(:), allocatable :: Iam
      Iam="get_cube_centers,cube with halo"

      call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, status); _VERIFY(STATUS)
      GigaTrajInternalPtr => wrap%ptr
      grid_ = GigaTrajInternalPtr%CubedGrid
      call MAPL_Grid_interior(Grid_, i1,i2,j1,j2)
      imc = i2-i1 + 1
      jmc = j2-j1 + 1

      allocate(lon_center(imc+2, jmc+2))
      allocate(lat_center(imc+2, jmc+2))

      field = ESMF_FieldCreate(Grid_, ESMF_TYPEKIND_R8, name='halo', staggerLoc=ESMF_STAGGERLOC_CENTER,totalLWidth=[1,1],totalUWidth=[1,1],_RC)
      call ESMF_FieldGet(field,farrayPtr=ptr,_RC)
      ptr = 0.0d0
      call ESMF_FieldHaloStore(field,rh,_RC)

      call ESMF_GridGetCoord(grid_ , coordDim=1, localDE=0, &
              staggerloc=ESMF_STAGGERLOC_CENTER, &
              farrayPtr=centerX, _RC)

      ptr(1:imc,1:jmc)=centerX
      call ESMF_FieldHalo(field,rh, _RC)
      lon_center = ptr

      call ESMF_GridGetCoord(grid_ , coordDim=2, localDE=0, &
            staggerloc=ESMF_STAGGERLOC_CENTER, &
            farrayPtr=centerY, _RC)
      ptr = 0.0d0
      ptr(1:imc,1:jmc)=centerY
      call ESMF_FieldHalo(field,rh, _RC)
      lat_center = ptr
     
      lon_center = lon_center/MAPL_PI*180.0
      lat_center = lat_center/MAPL_PI*180.0
      where (lon_center < -180.) lon_center = lon_center + 360.
      where (lon_center >  180.) lon_center = lon_center - 360.
      call ESMF_FieldDestroy(field,_RC)
      call ESMF_FieldHaloRelease(rh,_RC)
      _RETURN(_SUCCESS)
   end subroutine get_cube_centers

   subroutine get_latlon_centers(GC, lon_center, lat_center, rc)
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      real, allocatable,     intent(out) :: lat_center(:), lon_center(:)
      integer, optional,   intent(  out) :: RC  
      integer :: i1, i2, j1, j2, imc, jmc, i, j, status, DIMS(3)
      real :: dlon, dlat
      type (GigaTrajInternal), pointer :: GigaTrajInternalPtr
      type (GigatrajInternalWrap)   :: wrap
      type(ESMF_Grid) :: grid_
      character(len=:), allocatable :: Iam
      Iam="get_latlon_centers, latlon with halo"


      call ESMF_UserCompGetInternalState(GC, 'GigaTrajInternal', wrap, status); _VERIFY(STATUS)
      GigaTrajInternalPtr => wrap%ptr
      grid_ = GigaTrajInternalPtr%LatLonGrid
      call MAPL_GridGet(Grid_, globalCellCountPerDim=DIMS, _RC)
      call MAPL_Grid_interior(Grid_, i1,i2,j1,j2)
      imc = i2-i1 + 1
      jmc = j2-j1 + 1

      allocate(lon_center(imc+2))
      allocate(lat_center(jmc+2))

      dlon = 360.0/dims(1)
        ! DE
        !lons_center = [(dlon*(i-1)+dlon/2., i= i1-1, i2+1)]
        ! DC   
      lon_center = [(dlon*(i-1) - 180.0 , i= i1-1, i2+1)]
        !PE
        !dlat = 180.0/dims(2) 
        !lats_center = [(-dlat/2. + dlat*j-90.0, j= j1-1, j2+1)] 
        !PC
      dlat = 180.0/(dims(2)-1)  ! PC
      lat_center = [(-90.0 + (j-1)*dlat, j= j1-1, j2+1)]
      where(lat_center <-90.) lat_center = -90.
      where(lat_center >90. ) lat_center =  90.
      _RETURN(_SUCCESS)
   end subroutine get_latlon_centers

end module
