module read_env

   use netcdf
   use pso_params_mod_landshared

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


   subroutine read_env_data(met_path)
      
      ! intent in:
      character(*), intent(in) :: met_path

      !etc:
      integer :: ncid, ndim, nvar, natt, unlid, len1, varid, test
      real, allocatable, dimension(:) :: values
      
      write(*,*) 'inside function'

      call check(nf90_open(met_path, nf90_nowrite, ncid))
      write(*,*) 'after open'
      call check(nf90_inquire(ncid, ndim, nvar, natt, unlid))
      write(*,*) 'after inquire'
      call check(nf90_inquire_dimension(ncid,1,len=len1))
      
      write(*,*) 'file opened'

      if (.not. allocated(pso_params%map_vals)) then
        allocate(pso_params%map_vals(len1))
      endif
        
      write(*,*) 'map vals allocated'

      if (.not. allocated(values)) then
         allocate(values(len1))
      endif

      write(*,*) 'values allocated'

      call check(nf90_inq_varid(ncid,'MAP',varid))
      call check(nf90_get_var(ncid, varid, values))
      pso_params%map_vals = values

      write(*,*) 'values assigned to map_vals'
      write(*,*) 'pso_params%map_vals'
      write(*,*) pso_params%map_vals
      call check(nf90_close(ncid))
   end subroutine read_env_data

end module read_env
