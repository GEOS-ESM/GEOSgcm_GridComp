module MAPL_ForcingSpecificationMod
  implicit none
  private

  public :: ForcingSpecification

  character(*), parameter :: NOT_FOUND = 'NOT_FOUND'

  type :: ForcingSpecification
     type(ESMF_State) :: state
     type(ESMF_Time) :: time
     character(:), allocatable :: label_format
     character(:), allocatable :: var_name_format
   contains
     procedure :: format_name
     procedure :: get_label
     procedure :: get_var_name

     procedure :: read_base
     procedure :: read_multi
     generic :: read_forcing => read_base, read_multi
  end type ForcingSpecification

  interface ForcingSpecification
     module procedure new_reader
  end interface ForcingSpecification

contains


  function new_reader(state, time, unused, label_format, var_name_format) result(reader)
    type(ForcingSpec) :: reader

    type(ESMF_State), intent(in) :: state
    type(ESMF_Time), intent(in) :: time
    class(KeywordEnforcer), optional, intent(in) :: unused
    character(*), optional, intent(in) :: label_format
    character(*), optional, intent(in) :: var_name_format

    this%state = state
    this%time = time
    
    if (present(label_format)) then
       reader%label_format = label_format
    else
       reader%label_format = '(%a)'
    end if

    if (present(var_name_format)) then
       reader%var_name_format = var_name_format
    else
       reader%var_name_format = '(%a)')
    end if

  end function new_reader
  
  function format_name(this, base_name, fmt, bin) result(str)
    character(:), allocatable :: str
    class(ForcingSpecification), intent(in) :: this
    character(*), intent(in) :: base_name
    integer, optional, intent(in) :: bin

    character(ESMF_MAXSTR) :: buffer
    integer :: bin_

    if (present(bin)) then
       bin_ = bin
    else
       bin_ = 0
    end if

    read(buffer,fmt) base_name, bin_
    str = trim(buffer)

  end function format_name

  !
  ! Construct label for datafile in config
  !
  ! Note: Although in theory
  ! this function could encounter an exceptional condition, it would
  ! be from a typo in source rather than a missing/failing resource.
  !
  function get_label(this, base_name, bin) result(label)
    character(:), allocatable :: label
    class(ForcingSpecification), intent(in) :: this
    character(*), intent(in) :: base_name
    integer, optional, intent(in) :: bin

    label = this%format(base_name, this%label_format, bin)

  end function get_label

  ! Construct variable name in forcing file
  function get_var_name(this, base_name, bin) result(var_name)
    character(:), allocatable :: label
    class(ForcingSpecification), intent(in) :: this
    character(*), intent(in) :: base_name
    integer, optional, intent(in) :: bin

    label = this%format(base_name, this%label_format, bin)

  end function get_var_name


  subroutine get_datafile(this, base_name, datafile, bin, rc)
    class(ForcingSpecification), intent(in) :: this
    character(*), intent(in) :: base_name
    character(:), allocatable, intent(out) :: datafile
    integer, optional, intent(in) :: bin
    integer, optional, intent(out) :: rc

    integer :: status
    character(:), allocatable :: label
    character(ESMF_MAXSTR) :: buffer

    label = this%get_label(base_name, bin)
    call MAPL_GetResource(state, buffer, label=label, default=NOT_FOUND, rc=status)
    ASSERT_(status==0, 'label not found: "'//label//'"')

    datafile = trim(buffer)
    RETURN_(SUCCESS_)

  end subroutine get_datafile


  subroutine read_base(this, array, base_name, default, bin, rc)
    class(ForcingSpecification), intent(in) :: this
    real, intent(out) :: array(:)
    character(:), intent(in) :: base_name
    real, optional, intent(in) :: default
    integer, optional, intent(in) :: bin
    integer, intent(out), optional :: rc

    integer :: status
    character(:), allocatable :: datafile
    character(:), allocatable :: var_name

    call this%get_datafile(base_name, datafile, bin=bin, rc=status)
    _VERIFY(status)
    
    if (filename == NOT_FOUND) then
       if (present(default_)) then
          array = default
       else
          array = 0.
       end if
    else
       var_name = this%get_var_name(base_name, bin=bin)
       call MAPL_ReadForcing(this%state, var_name, filename, this%time, array, rc=status)
       VERIFY_(status)
    end if

    RETURN_(SUCCESS_)
  end subroutine read_base

  subroutine read_multi(this, array, base_name, default, rc)
    class(ForcingSpecification), intent(in) :: this
    real, intent(out) :: array(:,:)
    character(:), intent(in) :: base_name
    real, optional, intent(in) :: default
    integer, intent(out), optional :: rc

    integer :: status
    integer :: bin, num_bins
    character(:), allocatable :: datafile
    character(:), allocatable :: var_name

    num_bins = size(array,2)
    do bin = 1, num_bins
       call this%read(array(:,i), base_name, default, bin=bin, rc=staus)
       VERIFY_(status)
    end do
    
  end subroutine read_multi

  
end module MAPL_ForcingSpecificationMod
