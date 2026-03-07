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
     type(T_CVARS), pointer :: vars(:) => null()
     type(ESMF_State) :: StateIn, StateOut, StateTmp
     type(ESMF_CplComp) :: ccs
   contains
     procedure :: g2ginitialize
     procedure :: accumulate
     procedure :: regrid
     procedure, private :: MAPL_ConnectSetName
     procedure, private :: MAPL_ConnectSetXform
     procedure, private :: MAPL_ConnectSetVars
     procedure, private :: MAPL_ConnectSetDt
     generic :: set => MAPL_ConnectSetName
     generic :: set => MAPL_ConnectSetXform
     generic :: set => MAPL_ConnectSetVars
     generic :: set => MAPL_ConnectSetDt
  end type T_GC2GC_STATE

  public T_GC2GC_STATE

contains
  subroutine g2ginitialize(this, clock, rc)
    class (T_GC2GC_STATE), intent(inout) :: this
    type(ESMF_Clock), intent(inout) :: clock
    integer, optional, intent(out) :: rc

    integer :: status
    integer :: n, k, rank, dims
    type (MAPL_VarSpec), pointer :: SSPEC(:)=> NULL()
    type (ESMF_Field) :: field, fieldCopy
    type (ESMF_TypeKind_flag) :: tk
    character (len=ESMF_MAXSTR), allocatable  :: itemNameList(:)

    _RETURN_UNLESS(this%average)

    ! fill the export state of the coupler
    call ESMF_StateGet(this%stateIn, ITEMCOUNT=N, _RC)
    allocate(itemNameList(N), _STAT)

    call ESMF_StateGet(this%stateIn, ITEMNAMELIST=itemNamelist, _RC)

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



    end do
    deallocate(itemNameList, _STAT)
    
!         CCSetServ
    call ESMF_CplCompSetServices (this%CCS, &
         GenericCplSetServices, _RC )

    call MAPL_CplCompSetVarSpecs(this%CCS, &
                                      SSPEC,&
                                      SSPEC,_RC)
    
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

    stateIn = this%stateIn
    if (this%average) stateIn = this%stateTmp

    !??? VM barrier???

    do n=1,size(this%vars)
       ! we might need to check the field is in the state. If not, it might be Ok
       call ESMF_StateGet(stateIn, itemName=this%vars(n)%vIn, field=field, _RC)
       !!      call ESMF_StateGet(g2g%stateIn, g2g%vars(n)%vIn, field, _RC)
       ! get rank,typekind to use the appropriate overload
       call ESMF_FieldGet(field, rank=rank, _RC)
       select case (rank)
          case (1)
             call MAPL_GetPointer(this%STATEout, ptr1out, this%vars(n)%Vout, notFoundOK=.true., _RC)
             call MAPL_GetPointer(this%STATEin, ptr1in, this%vars(n)%Vin, notFoundOK=.true., _RC)

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
    deallocate(this%vars)

    allocate(TMP(I+1),stat=STATUS)

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
    if (present(name)) this%name = name
    this%stateIn = stateIn
    this%stateOut = stateOut

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
    
    if (this%average) then
       ! create couplers as needed
       this%CCS = ESMF_CplCompCreate (                  &
            NAME = this%name, &
            contextFlag = ESMF_CONTEXT_PARENT_VM,              &
            _RC )
       this%stateTmp = ESMF_StateCreate(name='TmpInOut',_RC)
    end if

    _RETURN(ESMF_SUCCESS)
  end subroutine MAPL_ConnectSetDt

end module gc2gc

