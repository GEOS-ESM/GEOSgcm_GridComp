#include "MAPL_Generic.h"

module reservoirMod

  use MAPL

  implicit none
  private
  public :: res_state, Reservoir

  type RES_STATE ! reservoir-related variables
     logical              :: use_res
     integer, allocatable :: active_res(:)
     real,    allocatable :: cap_res(   :)    ! m3
     real,    allocatable :: alpha_res( :)
   contains
     procedure :: calc
  end type RES_STATE

  interface Reservoir
     module procedure new_Reservoir
  end interface Reservoir

  !-----------------------------------------

contains

  !---------------------------------------------------------------------------------
  ! Initialization subroutine for reservoirs
  
  function new_Reservoir(GC, use_res, rc) result (res)

    type(RES_STATE)                      :: res
    type(ESMF_GridComp),  intent(inout)  :: GC 
    logical,              intent(in)     :: use_res

    integer, optional,    intent(out)    :: rc

    ! ---------------
    type(MAPL_MetaComp), pointer   :: MAPL

    type (ESMF_State       )    :: INTERNAL
    real, dimension(:), pointer :: ACTIVE_RES_RS
    real, dimension(:), pointer :: CAP_RES_RS
    real, dimension(:), pointer :: RRM_ALPHA_RES_RS 

    integer :: STATUS

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL,  RC=STATUS)
    call MAPL_GetPointer(INTERNAL, ACTIVE_RES_RS,     'ACTIVE_RES',     RC=STATUS ) 
    call MAPL_GetPointer(INTERNAL, CAP_RES_RS,        'CAP_RES',        RC=STATUS )
    call MAPL_GetPointer(INTERNAL, RRM_ALPHA_RES_RS,  'RRM_ALPHA_RES',  RC=STATUS ) 

    allocate(res%active_res, source=int(ACTIVE_RES_RS))
    allocate(res%cap_res,    source=    CAP_RES_RS    )
    allocate(res%alpha_res,  source=RRM_ALPHA_RES_RS  )

    if(use_res.eqv..True.)then
      res%use_res   = .True.
    else
      res%use_res   = .False.
      res%active_res=0
    endif

    _RETURN(_SUCCESS)

  end function new_Reservoir
  
  !---------------------------------------------------------------------------------
  ! Reservoir calculation subroutine
  
  subroutine calc( this, Q_riv_out, Q_res, Wr_res, dt, rc)

    class(RES_STATE),  intent(in)    :: this
    real,              intent(in)    :: Q_riv_out(:)    ! outflow from river (going into reservoir)   [m3/s]
    real,              intent(out)   :: Q_res(:)        ! reservoir outflow                           [m3/s]
    real,              intent(inout) :: Wr_res(:)       ! reservoir storage                           [m3  ]
    real,              intent(in)    :: dt              ! routing model time step                     [s   ]
    integer, optional, intent(out)   :: rc

    ! ------------------------

    real    :: Qin_res, alp_res, tmpfac   ! Variables for inflow, coefficients, and factors
    integer :: i, status

    Q_res = 0.

    do i = 1, size(this%active_res)

       if (this%active_res(i) == 1) then                  ! reservoir is active

          ! Inflow to reservoir is outflow from river

          Qin_res = Q_riv_out(i)

          Q_res(i) = this%alpha_res(i) * Wr_res(i)

          ! Ensure outflow is within reasonable bounds
          Q_res(i) = max(0.0, Q_res(i))                      ! Ensure non-negative outflow
          Q_res(i) = min(Q_res(i), Wr_res(i) / dt + Qin_res) ! Limit outflow to prevent exceeding inflow and storage

          !if (fld_res == 1) Q_res = min(Q_res, Qfld_thres)  ! Limit outflow for flood control

          Wr_res(i) = Wr_res(i) + dt * (Qin_res - Q_res(i))  ! Update water storage in the reservoir
          Wr_res(i) = max(0.0, Wr_res(i))                    ! Ensure non-negative storage
          
          ! If the storage exceeds capacity, adjust outflow and storage
          if (this%cap_res(i) > 0. .and. Wr_res(i) > this%cap_res(i)) then
             Q_res(i)  = Q_res(i) + (Wr_res(i) - this%cap_res(i)) / dt      ! Adjust outflow for overflow
             Wr_res(i) = this%cap_res(i)                                    ! Limit storage to reservoir capacity
          endif
          
       endif
    enddo
    
    _RETURN(_SUCCESS)
    
  end subroutine calc
  
end module reservoirMod

! =============== EOF ===================================================================
