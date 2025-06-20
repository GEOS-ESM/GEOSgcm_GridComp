module findloc_mod

   implicit none

   private
   public :: findloc

   contains

   function findloc(array, value)

      integer, intent(in) :: array(:)
      integer, intent(in) :: value
      integer :: findloc(1)

      integer :: num_elements, i

      num_elements = size(array)

      findloc(1) = 0
      do i = 1, num_elements
         if (array(i) == value) then
               findloc(1) = i
               exit
         endif
      end do

   end function findloc

end module findloc_mod
