module SHOCPARAMS

 implicit none

 type SHOCPARAMS_TYPE
    integer :: CLDLEN
    integer :: SUS12LEN
    integer :: BUOYOPT
    real    :: LAMBDA
    real    :: TSCALE
    real    :: VONK
    real    :: CKVAL
    real    :: CEFAC
    real    :: CESFAC
 endtype SHOCPARAMS_TYPE

end module SHOCPARAMS
