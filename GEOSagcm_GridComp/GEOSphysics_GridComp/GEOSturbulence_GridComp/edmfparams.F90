module EDMFPARAMS

 implicit none

 type EDMFPARAMS_TYPE
    integer :: DISCRETE
    integer :: IMPLICIT
    integer :: ENTRAIN
    integer :: DOCLASP
    integer :: THERMAL_PLUME
    integer :: TEST
    integer :: DEBUG
    real    :: STOCHFRAC
    real    :: ENTWFAC
    real    :: EDFAC
    real    :: ENT0
    real    :: ET
    real    :: ALPHATH
    real    :: ALPHAQT
    real    :: ALPHAW
    real    :: PWMAX
    real    :: PWMIN
    real    :: WA
    real    :: WB
    real    :: AU0
    real    :: CTH1
    real    :: CTH2
    real    :: RH0_QB
    real    :: C_KH_MF
    real    :: ICE_RAMP
 endtype EDMFPARAMS_TYPE

end module EDMFPARAMS
