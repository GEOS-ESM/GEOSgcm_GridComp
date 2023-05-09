#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!                       NASA/GSFC, GMAO, Code 610.1                      !
!-------------------------------------------------------------------------
!BOP
!
!
! !MODULE:  Bundle_IncrementMod --- Computes tracer increments and puts them into a bundle
! 
! !INTERFACE:
!

MODULE  MBundle_IncrementMod

! !USES:



!
! !DESCRIPTION: This module can be used to compute model variable increments if the
!               following conditions are met:
!                  1) A non-increment, or "parent", bundle exists (e.g. TR, TRADV, MTR)
!                  2) The increment is computed as: (X_t2 - X_t1)/(t2-t1)
!                For new increment bundles the user must do at least the following:
!                  1) Define an increment bundle in Set Services of a gridded component.



!EOP
!-------------------------------------------------------------------------

!BOC

   USE ESMF
   USE MAPL

   IMPLICIT NONE
   PRIVATE
 
! !PUBLIC MEMBER FUNCTIONS: 
   PUBLIC Initialize_IncMBundle_init
   PUBLIC Initialize_IncMBundle_run
   PUBLIC Compute_IncMBundle

   type MAPPING
      integer :: n
      character(len=ESMF_MAXSTR), allocatable :: fieldName(:)
   end type MAPPING

   type (Mapping), allocatable :: map(:)
   
   character(len=ESMF_MAXSTR)                :: fieldname
   character(len=ESMF_MAXSTR)                :: org_bundle       ! Original bundle name
   character(len=ESMF_MAXSTR)                :: inc_bundle       ! Increment bundle name
   integer                                   :: i, ppos, NQ

 CONTAINS

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!BOP
!
!  !ROUTINE:  Initialize_IncMBundle_init - Initialize increment bundle with fields from a 
!             "parent" bundle. This must be called in the Initialize method of a gridded component.

   SUBROUTINE Initialize_IncMBundle_init(GC, state1, state2, RC)
   
     IMPLICIT NONE
  
     ! ARGUMENTS
     type(ESMF_GridComp),        intent(in)        :: GC          ! Gridded component 
     type(ESMF_State),           intent(in)        :: state1      ! Original bundle state
     type(ESMF_State),           intent(in)        :: state2      ! Increment bundle state
     integer, optional,          intent(  out)     :: RC          ! Error code


     ! TYPES and VARIABLES
     integer                                    :: STATUS
     type(ESMF_FieldBundle)                     :: BUNDLE, BUNDLEi
     type(ESMF_Config)                          :: cf               ! AGCM.rc
     character(len=ESMF_MAXSTR)                 :: IAm, valueOld
     character(len=ESMF_MAXSTR)                 :: longname         ! longname metadata description
     character(len=2)                           :: suffix           ! suffix appended to variable name
     integer                                    :: nCols
     integer :: n, k, kg, fieldRank
     integer :: j, j1, j2, n0
     integer, parameter :: nmax=9
     logical :: found
     type(ESMF_Field), pointer :: fields(:)
     character(len=ESMF_MAXSTR), allocatable    :: NAMES(:), RULES(:)
     character(len=ESMF_MAXSTR)    :: splitNames(nmax)
     type(ESMF_Field)                          :: field, TempField
     character(len=:), allocatable :: s
     
     Iam = "Initialize_IncMBundle_init"

! ============================================================================

! Begin...
     org_bundle = 'MTR'
     inc_bundle = 'MCHEMTRI'
     suffix = 'IM'
     longname = '_due_to_moist_processes'

     call ESMF_GridCompGet ( GC, config=cf, RC=STATUS )
     VERIFY_(STATUS)

     call ESMF_ConfigGetDim (cf, NQ, nCols, label=(trim(inc_bundle)//'_increments::'), rc=status)

     if (NQ > 0) then
        call ESMF_ConfigFindLabel (cf, (trim(inc_bundle)//'_increments::'), rc=STATUS)
        VERIFY_(STATUS)

        allocate (NAMES(NQ), RULES(NQ), _STAT)

        do i = 1, NQ
           call ESMF_ConfigNextLine(cf, _RC)
           call ESMF_ConfigGetAttribute(cf, NAMES(i), _RC)
           call ESMF_ConfigGetAttribute(cf, RULES(i), _RC)
        enddo

! Fill the increments bundle with fields from "parent" bundle 
!------------------------------------------------------------------
        call ESMF_StateGet(state1, org_bundle, BUNDLE,   rc=STATUS)
        VERIFY_(STATUS)

        call ESMF_StateGet(state2, inc_bundle, BUNDLEi, rc=STATUS)
        VERIFY_(STATUS)
  
        allocate(map(NQ), _STAT)
        
        do i = 1, NQ

           call ESMF_FieldBundleGet (BUNDLE, NAMES(i), field=field, rc=status)
           _ASSERT(status==0,trim(NAMES(i))//' is not valid. It likely does not exist in '//trim(org_bundle))

           call ESMF_FieldGet (field, name=fieldname, _RC)

        ! If the field rank is 3, create rank 2 array (reduction is due to vertical integration
           call ESMF_FieldGet(FIELD, dimCount=fieldRank, _RC)
           _ASSERT(fieldRank == 3, "Supporting only rank 3 fields")

           !KLUGDE: the only reason we call split to reduce rank 3 to rank 2
           ! we are using only 1 of the split fields,
           ! but at least we get the attributes correctly set
           call MAPL_FieldSplit(field, fields, _RC) 

           TempField = MAPL_FieldCreate (fields(1), name=(trim(fieldname)//suffix) ,DoCopy=.true., _RC)

           call ESMF_FieldGet(tempFIELD, dimCount=fieldRank, _RC)
           _ASSERT(fieldRank == 2, "Expecting rank 2 field")

           
           call ESMF_AttributeSet(tempField,name="DIMS",value=MAPL_DimsHorzOnly,_RC)
           call ESMF_AttributeSet(tempField,name="VLOCATION",value=MAPL_VLocationNone,_RC)
           call MAPL_FieldBundleAdd (BUNDLEi, TempField, rc=status)
           VERIFY_(STATUS)

           !ALT: for each of the entries in NAMES:
           ! generate "split" name (i.e. the entry in the original bundle),

           if (rules(i) == 'default') then
              call generateSplitName(NAMES(i), splitNames)

              ! the test if we even have splitting, and if YES, how many (n); then loop over n
              n = 1
              do k=2,nmax
                 call ESMF_FieldBundleGet(bundle,splitNames(k), ispresent=found, _RC)
                 if (.not.found) exit
                 n=n+1
              end do

              ! populate a mapping
              map(i)%n=n
              allocate(map(i)%fieldName(n),_STAT)
              do k=1,n
                 map(i)%fieldName(k) = splitNames(k)
              end do
              splitNames=""

           else
              ! process the list of aggregation names
              s = rules(i)
              ! 2 in the next expession refrects that we 1 fewer separator
              ! than items, and the original name is not the rules' list
              n=count([(s(k:k),k=1,len_trim(s))] == ',') + 2 

              ! populate a mapping
              map(i)%n=n
              allocate(map(i)%fieldName(n),_STAT)
              map(i)%fieldName(1) = fieldname
              j1=1
              n0 = len_trim(s)
              do k=2,n
                 j=index(s(j1:),',')
                 j2 = n0
                 if (j /= 0) j2= j-1
                 map(i)%fieldName(k) = s(j1:j2)
                 j1= j+1
              end do
           end if
        end do

        deallocate (NAMES, RULES, _STAT)



! Set Species Attributes
!----------------------------------
        do i = 1, NQ
           call ESMF_FieldBundleGet (BUNDLEi, fieldIndex=i, field=field, _RC )
           call ESMF_FieldGet (field, name=fieldname, _RC)
        
           if (fieldname==('AOADAYS'//suffix)) then
              call ESMF_AttributeSET (field, name='UNITS', value='days s-1', _RC)
           else
              call ESMF_AttributeGET (field, name='UNITS', value=valueOld, _RC)
              ! Remove kg-1
              kg = index(valueOld,'kg-1')
              if (kg == 0) then ! did not find kg-1, use oldvalue
                 call ESMF_AttributeSET (field, name='UNITS', value=trim(valueOld)//' s-1', _RC)
              else
                 call ESMF_AttributeSET (field, name='UNITS', value=valueOld(1:kg-1)//'m-2 s-1', _RC)
              end if
           end if

           ppos = len(trim(fieldname))
           call ESMF_AttributeSET (field, name='LONG_NAME', value=('tendency_of_'//fieldname(1:ppos-2)//trim(longname)), _RC)
        end do

     end if ! NQ > 0

     _RETURN(ESMF_SUCCESS)

   contains
     subroutine generateSplitName(name, splitNameArray)
      character(len=*) :: name
      character(len=*) :: splitNameArray(:)

      integer :: n
      integer :: i

      n = size(splitNameArray)

      splitNameArray(1) = name
      do i=2,n
         write(splitNameArray(i),'(A,I3.3)') trim(name), i
      end do
      
      _RETURN(ESMF_SUCCESS)
    end subroutine GenerateSplitName

   END SUBROUTINE Initialize_IncMBundle_init
!==============================================================================================================

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 610.1 GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
!  !ROUTINE:  Initialize_IncMBundle_run - Re-initialize increment bundle with data from the 
!  !          "parent" bundle within the Run method. The "parent" bundle is the non-increment 
!              bundle (e.g. TR, MTR, TRADV)

   SUBROUTINE Initialize_IncMBundle_run(state1, state2, dm, RC)

     IMPLICIT NONE

     ! ARGUMENTS
     type(ESMF_State),  intent(in)            :: state1           ! Original bundle state
     type(ESMF_State),  intent(in)            :: state2           ! Increment bundle state
     real,              intent(in)            :: dm(:,:,:) ! delp/mapl_grav
     integer, optional, intent(  out)         :: RC               ! Error code

     ! TYPES and VARIABLES
     integer                                            :: STATUS
     type(ESMF_FieldBundle)                             :: BUNDLE, BUNDLEi
     character(len=ESMF_MAXSTR)                         :: IAm
     character(len=ESMF_MAXSTR), allocatable            :: NAMES(:)
     real, dimension(:,:,:), pointer                    :: org_ptr
     real, dimension(:,:), pointer                      :: inc_ptr
     integer :: k
     
     Iam = "Initialize_IncMBundle_run"

! ============================================================================

! Begin...
     org_bundle = 'MTR'
     inc_bundle = 'MCHEMTRI'
     

!  !Initialize increment bundle in Run method before the child is called
!  !--------------------------------------------------------------------
     call ESMF_StateGet (state2, inc_bundle, BUNDLEi, _RC)
     call ESMF_StateGet (state1, org_bundle, BUNDLE, _RC)
     call ESMF_FieldBundleGet (BUNDLEi, fieldCount=NQ, _RC)

!  !Check if there is anything in the bundle.
     _RETURN_IF (NQ == 0)
     allocate (NAMES(NQ), _STAT)

     call ESMF_FieldBundleGet(BUNDLEi, itemOrderFlag=ESMF_ITEMORDER_ADDORDER, &
          fieldNameList=NAMES, _RC)

!    !Get increment data pointer and initialize value
     do i = 1, NQ

        call ESMFL_BundleGetPointerToData (BUNDLEi, trim(NAMES(i)), inc_ptr, _RC)

        inc_ptr = 0.0
        ! aggregage bins
        do k=1,map(i)%n
           call ESMFL_BundleGetPointerToData (BUNDLE, map(i)%fieldName(k),&
                org_ptr, _RC)
           inc_ptr = inc_ptr + sum(org_ptr*dm,3) ! vertical integration, mass weighted
        end do
     end do
     deallocate(NAMES)

    _RETURN(ESMF_SUCCESS)
  END SUBROUTINE Initialize_IncMBundle_run
!=======================================================================================
!
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 610.1 GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
!  !ROUTINE:  Compute_IncMBundle - Compute the increment after the child has run.

  SUBROUTINE Compute_IncMBundle(state1, state2, META, DM, RC)

    IMPLICIT NONE

    ! ARGUMENTS
    type(ESMF_State),               intent(in)         :: state1           ! Original bundle state
    type(ESMF_State),               intent(in)         :: state2           ! Increment bundle state
    type(MAPL_MetaComp), pointer,   intent(in)         :: META
    real,              intent(in)            :: dm(:,:,:) ! delp/mapl_grav
    integer, optional,              intent(  out)      :: RC               ! Error code


    ! TYPES and VARIABLES
    integer                                            :: STATUS
    type(ESMF_FieldBundle)                             :: BUNDLE, BUNDLEi
    character(len=ESMF_MAXSTR)                         :: IAm
    character(len=ESMF_MAXSTR), allocatable            :: NAMES(:)
    real                                               :: DT
    real, dimension(:,:,:), pointer                    :: org_ptr
    real, dimension(:,:),   pointer                    :: inc_ptr
    integer :: k
    real, allocatable :: tmp(:,:)

    Iam = "Compute_IncMBundle"
! ============================================================================
! Begin...
    org_bundle = 'MTR'
    inc_bundle = 'MCHEMTRI'

    call ESMF_StateGet (state2, inc_bundle, BUNDLEi, _RC)
    call ESMF_FieldBundleGet (BUNDLEi, fieldCount=NQ, _RC)

!  !Check if there is anything in the bundle.
    _RETURN_IF (NQ == 0)
    call ESMF_StateGet (state1, org_bundle, BUNDLE, _RC)

    call MAPL_GetResource(META, DT, label="RUN_DT:", _RC)
    allocate (NAMES(NQ), _STAT)

    call ESMF_FieldBundleGet(BUNDLEi, itemOrderFlag=ESMF_ITEMORDER_ADDORDER, &
         fieldNameList=NAMES, _RC)

!    !Get pointers to data
    do i = 1, NQ
       call ESMFL_BundleGetPointerToData (BUNDLEi, trim(NAMES(i)), inc_ptr, _RC)
       if (.not. allocated(tmp)) then
          allocate(tmp(size(inc_ptr,1), size(inc_ptr,2)), _STAT)
       end if
       tmp = 0.0
       ! aggregage bins
       do k=1,map(i)%n
          call ESMFL_BundleGetPointerToData (BUNDLE, map(i)%fieldName(k),&
               org_ptr, _RC)
          tmp = tmp + sum(org_ptr*dm,3) ! vertical integration, mass weighted
       end do
         
       !Compute increment and update pointer
       inc_ptr = (tmp-inc_ptr)/DT

    end do
    deallocate(tmp)
    deallocate(NAMES)
 
    _RETURN(ESMF_SUCCESS)
  END SUBROUTINE Compute_IncMBundle


END MODULE MBundle_IncrementMod
