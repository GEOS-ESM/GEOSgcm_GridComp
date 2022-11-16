module module_gate
  IMPLICIT NONE

  !-for BRAMS runs, set use_gate=.false.
  LOGICAL, PARAMETER :: use_gate = .false. 
  
  !- Here are the place for data related with the GATE soundings
  INTEGER, PARAMETER :: gate = 1   ! flag to turn on/off : 1/0
  INTEGER, PARAMETER :: klon = 161 ! number of soundings for gate
  INTEGER, PARAMETER :: klev = 41  ! number of vertical levels
  INTEGER, PARAMETER :: ktrac= 2   ! number of chemical tracers
  INTEGER, PARAMETER :: levs=klev  
  INTEGER, PARAMETER :: nvar_grads=200
  
  TYPE cupout_vars
    REAL             ,POINTER     :: varp(:,:)
    CHARACTER(LEN=80),ALLOCATABLE :: varn(:)
  END  TYPE cupout_vars
  
  TYPE (cupout_vars), ALLOCATABLE :: cupout(:)
  
  REAL, DIMENSION(klon,klev):: pgeo,ppres,ptemp,pq,pu,pv,pvervel, &
                               zrvten,ztten,zq1,zq2,zqr,zadvt,zadvq

  INTEGER :: JL
  character (len=128) :: runname, runlabel

end module module_gate
