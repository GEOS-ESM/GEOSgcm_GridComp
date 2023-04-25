module SHOCPARAMS

 implicit none

 type SHOCPARAMS_TYPE
    integer :: LENOPT
    integer :: BUOYOPT
    real    :: LAMBDA
    real    :: TSCALE
    real    :: VONK
    real    :: CKVAL
    real    :: CEFAC
    real    :: CESFAC
    real    :: LENFAC1
    real    :: LENFAC2
    real    :: LENFAC3
    real    :: KRADFAC
    real    :: CLDLEN
    real    :: BRUNTMIN
 endtype SHOCPARAMS_TYPE

end module SHOCPARAMS
