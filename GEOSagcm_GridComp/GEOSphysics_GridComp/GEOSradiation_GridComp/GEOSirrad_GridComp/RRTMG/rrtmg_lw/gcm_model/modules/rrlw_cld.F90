!------------------------------------------------------------------
! rrtmg_lw cloud property coefficients

! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

   module rrlw_cld

      implicit none
      save

      real, dimension(2)      :: absice0
      real, dimension(2,5)    :: absice1
      real, dimension(43,16)  :: absice2
      real, dimension(46,16)  :: absice3
      real, dimension(200,16) :: absice4
      real, dimension(58,16)  :: absliq1

      ! band to second index into absice1
      integer, parameter :: ice1b(16) = &
        [1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5]

   end module rrlw_cld

