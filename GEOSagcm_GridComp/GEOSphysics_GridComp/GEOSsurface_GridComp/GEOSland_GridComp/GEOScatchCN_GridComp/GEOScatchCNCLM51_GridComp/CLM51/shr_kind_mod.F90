MODULE shr_kind_mod

  use MAPL_ConstantsMod, ONLY:  &
       MAPL_R8                  &
       MAPL_R4                  &
       MAPL_RN                  &
       MAPL_I8                  &
       MAPL_I4                  &
       MAPL_IN                
  
  !----------------------------------------------------------------------------
  ! precision/kind constants add data public
  !----------------------------------------------------------------------------
  public
  integer,parameter :: SHR_KIND_R8 = MAPL_R8                ! 8 byte real
  integer,parameter :: SHR_KIND_R4 = MAPL_R4                ! 4 byte real
  integer,parameter :: SHR_KIND_RN = MAPL_RN                ! native real
  integer,parameter :: SHR_KIND_I8 = MAPL_I8                ! 8 byte integer
  integer,parameter :: SHR_KIND_I4 = MAPL_I4                ! 4 byte integer
  integer,parameter :: SHR_KIND_IN = MAPL_IN                ! native integer
  integer,parameter :: SHR_KIND_CS = 80                     ! short char
  integer,parameter :: SHR_KIND_CM = 160                    ! mid-sized char
  integer,parameter :: SHR_KIND_CL = 256                    ! long char
  integer,parameter :: SHR_KIND_CX = 512                    ! extra-long char
  integer,parameter :: SHR_KIND_CXX= 4096                   ! extra-extra-long char

END MODULE shr_kind_mod
