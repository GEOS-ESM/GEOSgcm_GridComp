module ProcessLibraryParameters

  use MAPL_Constants, only: MAPL_ALHL, MAPL_GRAV, MAPL_CP

  implicit none

  private
  public :: alhlbcp, gravbcp

  real, parameter :: alhlbcp = MAPL_ALHL/MAPL_CP
  real, parameter :: gravbcp = MAPL_GRAV/MAPL_CP

end module ProcessLibraryParameters
