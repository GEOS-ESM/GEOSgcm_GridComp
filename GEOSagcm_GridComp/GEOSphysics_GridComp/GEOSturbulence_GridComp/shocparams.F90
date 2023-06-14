module SHOCPARAMS

 implicit none

 type SHOCPARAMS_TYPE
    integer :: LENOPT
    integer :: BUOYOPT
    real    :: PRNUM
    real    :: LAMBDA
    real    :: TSCALE
    real    :: CKVAL
    real    :: CEFAC
    real    :: CESFAC
    real    :: LENFAC1
    real    :: LENFAC2
    real    :: LENFAC3
 endtype SHOCPARAMS_TYPE

end module SHOCPARAMS
