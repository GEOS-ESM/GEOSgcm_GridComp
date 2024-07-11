module read_env

   use netcdf
   use catch_constants
   use pso_params, only: pso_vals

   implicit none

   private
   public :: read_env_data

contains

   SUBROUTINE check(istatus)
      USE netcdf
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: istatus
      IF (istatus /= nf90_noerr) THEN
        write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
      END IF
   END SUBROUTINE check 


   subroutine read_env_data(precip_fname,lai_fname,sand_fname,k0_fname,canopy_fname)
      
      ! intent in:
      character(*), intent(in) :: precip_fname ! path to the precip .nc4 file
      character(*), intent(in) :: lai_fname ! path to the lai .nc4 file
      character(*), intent(in) :: sand_fname ! path to the sand % .nc4 file
      character(*), intent(in) :: k0_fname ! path to the k0 .nc4 file
      character(*), intent(in) :: canopy_fname ! path to the canopy height .nc4 file


      !etc:
      integer :: ncid, ndim, nvar, natt, unlid, len1, varid, test
      real, allocatable, dimension(:) :: values, tile_nums
      
      ! lets get and assign the precip files
      ! open the precip .nc4
      call check(nf90_open(precip_fname, nf90_nowrite, ncid))
      call check(nf90_inquire(ncid, ndim, nvar, natt, unlid))
      call check(nf90_inquire_dimension(ncid,1,len=len1))
      ! allocate inside pso_params if needed
      if (.not. allocated(pso_vals%map_vals)) then
        allocate(pso_vals%map_vals(len1))
        allocate(pso_vals%all_tile_nums(len1))
      endif
      ! allocate the len and values for use here if needed
      if (.not. allocated(values)) then
         allocate(values(len1))
         allocate(tile_nums(len1))
      endif
      ! get the precipitation values and assign it
      call check(nf90_inq_varid(ncid,'vals',varid))
      call check(nf90_get_var(ncid, varid, values))
      pso_vals%map_vals = values
      ! get the tile information and assign it
      call check(nf90_inq_varid(ncid,'tile',varid))
      call check(nf90_get_var(ncid, varid, tile_nums))
      pso_vals%all_tile_nums = tile_nums
      ! close the precipitation data
      call check(nf90_close(ncid))

      ! let's do it again for lai
      ! open the lai .nc4
      call check(nf90_open(lai_fname, nf90_nowrite, ncid))
      call check(nf90_inquire(ncid, ndim, nvar, natt, unlid))
      call check(nf90_inquire_dimension(ncid,1,len=len1))
      ! allocate inside pso_vals if needed
      if (.not. allocated(pso_vals%lai_vals)) then
        allocate(pso_vals%lai_vals(len1))
      endif
      ! get the lai values and assign it
      call check(nf90_inq_varid(ncid,'lai',varid))
      call check(nf90_get_var(ncid, varid, values))
      pso_vals%lai_vals = values
      ! close the lai data
      call check(nf90_close(ncid))
      
      ! let's do it again for sand
      ! open the sand .nc4
      call check(nf90_open(sand_fname, nf90_nowrite, ncid))
      call check(nf90_inquire(ncid, ndim, nvar, natt, unlid))
      call check(nf90_inquire_dimension(ncid,1,len=len1))
      ! allocate inside pso_vals if needed
      if (.not. allocated(pso_vals%sand_vals)) then
        allocate(pso_vals%sand_vals(len1))
      endif
      ! get the lai values and assign it
      call check(nf90_inq_varid(ncid,'sand_perc',varid))
      call check(nf90_get_var(ncid, varid, values))
      pso_vals%sand_vals = values
      ! close the lai data
      call check(nf90_close(ncid))
      
      ! let's do it again for k0
      ! open the k0 .nc4
      call check(nf90_open(k0_fname, nf90_nowrite, ncid))
      call check(nf90_inquire(ncid, ndim, nvar, natt, unlid))
      call check(nf90_inquire_dimension(ncid,1,len=len1))
      ! allocate inside pso_vals if needed
      if (.not. allocated(pso_vals%k0_vals)) then
        allocate(pso_vals%k0_vals(len1))
      endif
      ! get the lai values and assign it
      call check(nf90_inq_varid(ncid,'k_sat',varid))
      call check(nf90_get_var(ncid, varid, values))
      pso_vals%k0_vals = values
      ! close the lai data
      call check(nf90_close(ncid))

      ! let's do it again for canopy height
      ! open the canopy height .nc4
      call check(nf90_open(canopy_fname, nf90_nowrite, ncid))
      call check(nf90_inquire(ncid, ndim, nvar, natt, unlid))
      call check(nf90_inquire_dimension(ncid,1,len=len1))
      ! allocate inside pso_vals if needed
      if (.not. allocated(pso_vals%canopy_vals)) then
        allocate(pso_vals%canopy_vals(len1))
      endif
      ! get the canopy values and assign it
      call check(nf90_inq_varid(ncid,'vals',varid))
      call check(nf90_get_var(ncid, varid, values))
      pso_vals%canopy_vals = values
      ! close the lai data
      call check(nf90_close(ncid))

  end subroutine read_env_data

end module read_env
