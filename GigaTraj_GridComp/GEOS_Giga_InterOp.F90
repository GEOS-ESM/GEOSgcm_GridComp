module GEOS_Giga_InterOpMod
   use, intrinsic :: iso_c_binding, only : c_int, c_ptr, c_null_ptr, c_associated
   implicit none
   private

   interface

     function initGigaGridField3d(nlons, nlats, nzs, lons_ptr, lats_ptr, levs_ptr) result (field_ptr) bind(C, name="initGigaGridField3d")
       import :: c_int, c_ptr
       implicit none
       integer(c_int), intent(in), value :: nlons, nlats, nzs
       type(c_ptr) :: lons_ptr, lats_ptr, levs_ptr, field_ptr
     end function

   end interface

type(c_ptr), save :: field3d_ptr = c_null_ptr


contains

   subroutine init_gigatraj_obj(nlons, nlats, nzs, lons, lats, levs)
     integer, intent(in) :: nlons, nlats, nzs
     real, dimension(:), intent(in) :: lons, lats, levs

 

   end subroutine

end module GEOS_Giga_InterOpMod
