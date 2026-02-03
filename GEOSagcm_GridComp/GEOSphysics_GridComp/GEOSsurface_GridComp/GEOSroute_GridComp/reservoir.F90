#include "MAPL_Generic.h"

module reservoirMod

  use MAPL

  implicit none
  private
  public :: res_state, Reservoir

  !----Reservoir module constants----------

  real,    parameter :: fac_elec_a   = 0.30   ! Coefficient for hydropower calculation
  real,    parameter :: fac_elec_b   = 2.00   ! Exponent for hydropower calculation
  real,    parameter :: fac_irr_a    = 0.01   ! Coefficient for irrigation calculation (arid areas)
  real,    parameter :: fac_irr_b    = 3.00   ! Scaling factor for irrigation (arid areas)
  real,    parameter :: fac_sup_a    = 0.03   ! Coefficient for water supply calculation
  real,    parameter :: fac_sup_b    = 2.00   ! Exponent for water supply calculation
  real,    parameter :: fac_other_a  = 0.20   ! Coefficient for other reservoir types
  real,    parameter :: fac_other_b  = 2.00   ! Exponent for other reservoir types

  integer, parameter :: fac_fld      = 1      ! Flood control parameter

  real,    parameter :: fac_a_slake  = 0.003  ! Factor for small lakes
  real,    parameter :: fac_b_slake  = 0.40   ! Exponent for small lakes
  real,    parameter :: fac_a_llake  = 0.01   ! Factor for large lakes
  real,    parameter :: fac_b_llake  = 0.60   ! Exponent for large lakes
  
  real,    parameter :: thr_wid_lake = 1.e5   ! Threshold lake width (in m)


  type RES_STATE ! reservoir-related variables
     logical              :: use_res
     integer, allocatable :: active_res(:)
     integer, allocatable :: active_up( :,:)
     integer, allocatable :: type_res(  :)
     real,    allocatable :: cap_res(   :)    ! m3
     real,    allocatable :: wid_res(   :)    ! m
     integer, allocatable :: fld_res(   :)
     real,    allocatable :: Qfld_thres(:)    ! m3/s
     integer, allocatable :: cat2res(   :)
     real,    allocatable :: qres_acc(  :)
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
    real, dimension(:), pointer :: TYPE_RES_RS
    real, dimension(:), pointer :: CAP_RES_RS
    real, dimension(:), pointer :: WID_RES_RS
    real, dimension(:), pointer :: FLD_RES_RS

    integer :: STATUS

    call MAPL_GetObjectFromGC ( GC, MAPL, RC=STATUS)
    call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=INTERNAL,  RC=STATUS)
    call MAPL_GetPointer(INTERNAL, ACTIVE_RES_RS,  'ACTIVE_RES',  RC=STATUS ) 
    call MAPL_GetPointer(INTERNAL, TYPE_RES_RS,    'TYPE_RES',    RC=STATUS ) 
    call MAPL_GetPointer(INTERNAL, CAP_RES_RS,     'CAP_RES',     RC=STATUS ) 
    call MAPL_GetPointer(INTERNAL, WID_RES_RS,     'WID_RES',     RC=STATUS ) 
    call MAPL_GetPointer(INTERNAL, FLD_RES_RS,     'FLD_RES',     RC=STATUS ) 

    allocate(res%active_res, source=int(ACTIVE_RES_RS))
    allocate(res%type_res,   source=int(TYPE_RES_RS  ))
    allocate(res%cap_res,    source=    CAP_RES_RS   )
    allocate(res%wid_res,    source=    WID_RES_RS   )
    allocate(res%fld_res,    source=int(FLD_RES_RS   ))

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

          ! Determine alpha coefficient based on reservoir type  [1/s]

          tmpfac = 1.0 / (this%wid_res(i) / 1.e3) 
          
          select case (this%type_res(i))

          case  (1)               ! Irrigation reservoir 

             alp_res    = fac_irr_a   * ( tmpfac**fac_irr_b   ) / 3600.0
                                                                                         
          case  (2)               ! Hydropower reservoir                                 
                                                                                         
             alp_res    = fac_elec_a  * ( tmpfac**fac_elec_b  ) / 3600.0
                                                                                         
          case  (3)               ! Water supply reservoir                               
                                                                                         
             alp_res    = fac_sup_a   * ( tmpfac**fac_sup_b   ) / 3600.0
                                                                                         
          case  (0, 4, 5, 6)      ! Other reservoir types (generic alpha)                
                                                                                         
             alp_res    = fac_other_a * ( tmpfac**fac_other_b ) / 3600.0 
                                                                                         
          case  (-1)              ! Natural lake                                         
                                                                                         
             ! Determine lake type based on area and calculate alpha                     
             if (this%wid_res(i) >= thr_wid_lake) then                                   
                alp_res = fac_a_llake * ( tmpfac**fac_b_llake ) / 3600.0
             else                                                                        
                alp_res = fac_a_slake * ( tmpfac**fac_b_slake ) / 3600.0
             endif

          case default

             _ASSERT(.False., "ERROR: unknown reservoir type!")             

          end select

          Q_res(i) = alp_res * Wr_res(i)

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
