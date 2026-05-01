#include "MAPL_Generic.h"
module gc2gc
  use ESMF
  use MAPL
  use MAPL_GenericCplCompMod

  implicit none

  type T_CVARS
     character(len=ESMF_MAXSTR) :: vIn, vOut
  end type T_CVARS

  type T_GC2GC_STATE
     private
     character(len=ESMF_MAXSTR) :: name
     type (MAPL_LocStreamXFORM) :: XFORM
     integer :: nvars
     integer :: DT, COUPLE_DT
     logical :: average=.false.
     logical :: alreadyPrinted=.false.
     type(T_CVARS), pointer :: vars(:) => null()
     type(ESMF_State) :: StateIn, StateOut, StateTmp
     type(ESMF_CplComp) :: ccs
     type(ESMF_Alarm) :: alarm
   contains
     procedure :: g2ginitialize
     procedure :: accumulate
     procedure :: regrid
     procedure, private :: MAPL_ConnectSetName
     procedure, private :: MAPL_ConnectSetXform
     procedure, private :: MAPL_ConnectSetVars
     procedure, private :: MAPL_ConnectSetDt
     procedure, private :: MAPL_ConnectSetAverage
     procedure, private :: MAPL_ConnectSetAlarm
     generic :: set => MAPL_ConnectSetName
     generic :: set => MAPL_ConnectSetXform
     generic :: set => MAPL_ConnectSetVars
     generic :: set => MAPL_ConnectSetDt
     generic :: set => MAPL_ConnectSetAverage
     generic :: set => MAPL_ConnectSetAlarm
  end type T_GC2GC_STATE

  public T_GC2GC_STATE

contains
  subroutine g2ginitialize(this, clock, rc)
    class (T_GC2GC_STATE), intent(inout) :: this
    type(ESMF_Clock), intent(inout) :: clock
    integer, optional, intent(out) :: rc

    integer :: status
    integer :: n, k, rank, dims
    type (MAPL_VarSpec), pointer :: SSPEC(:)
    type (ESMF_Field) :: field, fieldCopy
    type (ESMF_TypeKind_flag) :: tk
    character (len=ESMF_MAXSTR), allocatable  :: itemNameList(:)

    _RETURN_UNLESS(this%average)

    ! create couplers as needed
    this%CCS = ESMF_CplCompCreate (                  &
         NAME = this%name, &
         contextFlag = ESMF_CONTEXT_PARENT_VM,              &
         _RC )
    this%stateTmp = ESMF_StateCreate(name='TmpInOut',_RC)

    ! fill the export state of the coupler
    call ESMF_StateGet(this%stateIn, ITEMCOUNT=N, _RC)
    allocate(itemNameList(N), _STAT)

    call ESMF_StateGet(this%stateIn, ITEMNAMELIST=itemNamelist, _RC)

    NULLIFY(SSPEC)
    k=0
    do n=1,size(this%vars)
       if (.not. any(itemNameList==this%vars(n)%vIn)) cycle
       k = k+1

       call ESMF_StateGet(this%stateIn, itemName=this%vars(n)%vIn, field=field, _RC)
       fieldCopy = MAPL_FieldCreate(field, name=this%vars(n)%vIn, doCopy=.true.,_RC)
       call MAPL_StateAdd(this%stateTmp, field=fieldCopy, _RC)
       ! get rank,typekind to use the appropriate overload
       call ESMF_FieldGet(field, rank=rank, typekind=tk, _RC)

       _ASSERT(tk==ESMF_TYPEKIND_R4, &
            "Currently the averaging coupler supports only R4. Var: "//&
            trim(this%vars(n)%vIn))

       select case (rank)
       case (1)
          DIMS = MAPL_DimsTileOnly
       case (2)
          DIMS = MAPL_DimsHorzOnly
          ! we are "lying". this is TileOnly + ungridded dims (but the coupler does not care)
       case (3)
          ! we are "lying". this is TileOnly + ungridded dims (but the coupler does not care)
          DIMS = MAPL_DimsHorzVert
       case default
          _FAIL('unsupported rank')
       end select

       call MAPL_VarSpecCreateInList(SSPEC,   &
            SHORT_NAME = this%vars(n)%vIn,                          &
            DIMS       = DIMS,                                &
            ACCMLT_INTERVAL= this%DT,                          &
            COUPLE_INTERVAL= this%COUPLE_DT,                         &
            _RC)

       if (MAPL_AM_I_ROOT()) then
          print *,'DEBUG:g2gini:',trim(this%vars(n)%vIn),&
               this%DT, this%COUPLE_DT, this%average, trim(this%name)
       end if

    end do
    deallocate(itemNameList, _STAT)

    if (MAPL_AM_I_Root()) then
       print *, "DEBUG: Averaging the following ",K," vars:"
!       call MAPL_VarSpecPrint(sspec, _RC)
    end if

!         CCSetServ
    call ESMF_CplCompSetServices (this%CCS, &
         GenericCplSetServices, _RC )

    call MAPL_CplCompSetVarSpecs(this%CCS, &
                                      SSPEC,&
                                      SSPEC,_RC)
    
    call MAPL_CplCompSetAlarm(this%CCS, this%alarm, _RC)
    
    call ESMF_CplCompInitialize (this%CCS, &
         importState=this%stateIn, & 
         exportState=this%stateTmp, &
         clock=CLOCK, userRc=status)
    _VERIFY(status)

    _RETURN(ESMF_SUCCESS)
  end subroutine g2ginitialize

  subroutine accumulate(this, clock, rc)
    class (T_GC2GC_STATE), intent(inout) :: this
    type(ESMF_Clock), intent(inout) :: clock
    integer, optional, intent(out) :: rc

    integer :: status

    _RETURN_UNLESS(this%average)

    call ESMF_CplCompRun(this%ccs, importState=this%stateIn, &
         exportState=this%stateTmp, clock=clock, userRc=status)
    _VERIFY(status)

    _RETURN(ESMF_SUCCESS)
  end subroutine accumulate

  subroutine regrid(this, rc)
    class (T_GC2GC_STATE), intent(inout) :: this
    integer, optional, intent(out) :: rc

    integer :: status
    integer :: n, k, m, rank
    real, pointer :: ptr1In(:)
    real, pointer :: ptr1Out(:)
    real, pointer :: ptr2In(:,:)
    real, pointer :: ptr2Out(:,:)
    real, pointer :: ptr3In(:,:,:)
    real, pointer :: ptr3Out(:,:,:)
    type (ESMF_Field) :: field
!    type (ESMF_TypeKind_flag) :: tk
    type (ESMF_State) :: StateIn
    character (len=ESMF_MAXSTR), allocatable  :: itemNameListIn(:)
    character (len=ESMF_MAXSTR), allocatable  :: itemNameListOut(:)
    logical :: skip
    
    stateIn = this%stateIn
    if (this%average) stateIn = this%stateTmp

    !??? VM barrier???

    call ESMF_StateGet(this%stateOut, ITEMCOUNT=N, _RC)
    allocate(itemNameListOut(N), _STAT)
    call ESMF_StateGet(this%stateOut, ITEMNAMELIST=itemNamelistOut, _RC)

    call ESMF_StateGet(this%stateIn, ITEMCOUNT=N, _RC)
    allocate(itemNameListIn(N), _STAT)
    call ESMF_StateGet(this%stateIn, ITEMNAMELIST=itemNamelistIn, _RC)

    do n=1,size(this%vars)
       ! we might need to check the field is in the state. If not, it might be Ok
       skip = .false.
       if (.not. any(itemNameListIn==this%vars(n)%vIn)) skip=.true.
       if (.not. any(itemNameListOut==this%vars(n)%vOut)) skip=.true.
       if (skip) cycle

       if (.not. this%alreadyPrinted) then
          if (MAPL_Am_I_Root()) then
             print *, "DEBUG: var "//trim(this%vars(n)%vIn),&
                  "  "//trim(this%name)
          end if
       end if

       call ESMF_StateGet(this%stateIn, this%vars(n)%vIn, field, _RC)
       ! get rank,typekind to use the appropriate overload
       call ESMF_FieldGet(field, rank=rank, _RC)
       select case (rank)
          case (1)
             call MAPL_GetPointer(this%stateOut, ptr1out, this%vars(n)%Vout, notFoundOK=.true., _RC)
             call MAPL_GetPointer(this%stateIn, ptr1in, this%vars(n)%Vin, notFoundOK=.true., _RC)       
             if (associated(ptr1Out) .and. associated(ptr1In)) then
                call MAPL_LocStreamTransform( ptr1Out, &
                     this%XFORM, ptr1In, _RC )
             end if
          case (2)
             call MAPL_GetPointer(this%STATEout, ptr2out, this%vars(n)%Vout, notFoundOK=.true., _RC)
             call MAPL_GetPointer(this%STATEin, ptr2in, this%vars(n)%Vin, notFoundOK=.true., _RC)

             if (associated(ptr2Out) .and. associated(ptr2In)) then
                do K=1,size(ptr2Out,2)
                   call MAPL_LocStreamTransform( ptr2Out(:,K), &
                        this%XFORM, ptr2In(:,K), _RC )
                enddo
             end if
          case (3)
             call MAPL_GetPointer(this%STATEout, ptr3out, this%vars(n)%Vout, notFoundOK=.true., _RC)
             call MAPL_GetPointer(this%STATEin, ptr3in, this%vars(n)%Vin, notFoundOK=.true., _RC)

             if (associated(ptr3Out) .and. associated(ptr3In)) then
                do M=1,size(ptr3Out,3)
                   do K=1,size(ptr3Out,2)
                      call MAPL_LocStreamTransform( ptr3Out(:,K,M), &
                           this%XFORM, ptr3In(:,K,M), _RC )
                   enddo
                enddo
             end if
       case default
          _FAIL("Unsupported rank")
       end select
    end do
    if (.not. this%alreadyPrinted) this%alreadyPrinted = .true.
    deallocate(itemNameListIn, itemNameListOut)
    _RETURN(ESMF_SUCCESS)

  end subroutine regrid

  subroutine MAPL_ConnectSetVars(this, varOut, varIn, rc)
    class (T_GC2GC_STATE), intent(INOUT) :: this
    character(len=*), intent(IN) :: varIn, varOut
    integer, optional, intent(OUT) :: rc

    integer :: status, I
    type (T_CVARS), pointer :: TMP(:)

    I = 0
    if (associated(this%vars)) I = size(this%vars)

    allocate(TMP(I+1),_STAT)

    TMP(1:I) = this%vars
    if(I>0) deallocate(this%vars)

    TMP(I+1)%Vin =  varIn
    TMP(I+1)%VOut =  varOut

    !      call move_alloc(from=tmp, to=this%vars)
    this%vars => tmp
    _RETURN(ESMF_SUCCESS)
  end subroutine MAPL_ConnectSetVars

  subroutine MAPL_ConnectSetXform(this, xform, rc)
    class (T_GC2GC_STATE), intent(INOUT) :: this
    type (MAPL_LocStreamXform), intent(IN) :: xform
    integer, optional, intent(OUT) :: rc

    this%xform = xform

    _RETURN(ESMF_SUCCESS)
  end subroutine MAPL_ConnectSetXform

  subroutine MAPL_ConnectSetName(this, stateIn, stateOut, name, rc)
    class (T_GC2GC_STATE), intent(INOUT) :: this
    type(ESMF_State), intent(IN) :: stateIn, stateOut
    character(len=*), optional, intent(IN) :: name
    integer, optional, intent(OUT) :: rc

    integer :: status

    this%name = ''
    this%stateIn = stateIn
    this%stateOut = stateOut
    if (present(name)) this%name = name

    _RETURN(ESMF_SUCCESS)
  end subroutine MAPL_ConnectSetName

  subroutine MAPL_ConnectSetDt(this, Dt, Couple_Dt, rc)
    class (T_GC2GC_STATE), intent(INOUT) :: this
    integer, intent(IN) :: dt, couple_dt
    integer, optional, intent(OUT) :: rc

    integer :: status
    
    this%dt = dt
    this%couple_dt = couple_dt

    this%average = dt/=couple_dt
    
    _RETURN(ESMF_SUCCESS)
  end subroutine MAPL_ConnectSetDt

  subroutine MAPL_ConnectSetAverage(this, average, rc)
    class (T_GC2GC_STATE), intent(INOUT) :: this
    logical, intent(IN) :: average
    integer, optional, intent(OUT) :: rc

    integer :: status
    
    this%average = average
    
    _RETURN(ESMF_SUCCESS)
  end subroutine MAPL_ConnectSetAverage

  subroutine MAPL_ConnectSetAlarm(this, alarm, rc)
    class (T_GC2GC_STATE), intent(INOUT) :: this
    type(ESMF_Alarm), intent(IN) :: alarm
    integer, optional, intent(OUT) :: rc

    integer :: status
    
    this%alarm = alarm
    
    _RETURN(ESMF_SUCCESS)
  end subroutine MAPL_ConnectSetAlarm

end module gc2gc

