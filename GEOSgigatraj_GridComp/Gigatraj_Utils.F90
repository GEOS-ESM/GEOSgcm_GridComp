module Gigatraj_UtilsMod
   use MAPL
   implicit none
   public :: parseCompsAndFieldsName
   public :: create_new_vars

contains

     subroutine parseCompsAndFieldsName(fields_line, CompNames, BundleNames, FieldNames, AliasNames)
        character(*), intent(in) :: fields_line
        character(len=ESMF_MAXSTR), allocatable, intent(out) :: CompNames(:)
        character(len=ESMF_MAXSTR), allocatable, intent(out) :: BundleNames(:)
        character(len=ESMF_MAXSTR), allocatable, intent(out) :: FieldNames(:)
        character(len=ESMF_MAXSTR), allocatable, intent(out) :: AliasNames(:)
        integer :: num_field, i, j, k, l, endl, num_
        character(len=:), allocatable :: tmp, tmp_bnf, tmp_f, tmp_alias
        num_field = 1
        k = 1
        do
          i = index(fields_line(k:),';')
          if (i == 0) exit
          if (trim(fields_line(i+1:)) =='') exit ! take care of the last unnecessay ";"
          k = k+i
          num_field = num_field+1
        enddo

        allocate(Fieldnames(num_field))
        allocate(Compnames(num_field))
        allocate(BundleNames(num_field))
        allocate(AliasNames(num_field))

        k    = 1
        num_ = 1

        do
          i = index(fields_line(k:),';')
          if (i == 0) then
             endl = len(fields_line)
          else
             endl = (k-1)+i-1
          endif
          tmp = fields_line(k:endl)

          j = index(tmp, '%%')
          if (j == 0) print*, "Wrong format of the comp%%field"
          Compnames(num_) = trim(adjustl(tmp(1:j-1)))
          tmp_bnf  = trim(adjustl(tmp(j+2:)))

         l = index(tmp_bnf, '%')
          if (l /=0) then
             BundleNames(num_) = tmp_bnf(1:l-1)
             tmp_f  = tmp_bnf(l+1:)
          else
             BundleNames(num_) = 'NONE'
             tmp_f  = tmp_bnf
          endif

          ! Aliasing....., Hard coded here
          l = index(tmp_f, '|')
          if (l /=0) then
             FieldNames(num_) = tmp_f(1:l-1)
             tmp_alias  = tmp_f(l+1:)
          else
             FieldNames(num_) = tmp_f
             tmp_alias  = tmp_f
          endif

          AliasNames(num_) = tmp_alias

          num_ = num_ + 1
          k = endl + 2
          if (num_ > num_field) exit
        enddo
     end subroutine parseCompsAndFieldsName

     subroutine create_new_vars(meta, formatter, long_name, short_name, units)
       type(FileMetadata), intent(inout) :: meta
       type(Netcdf4_fileformatter), intent(inout) :: formatter
       character(*), intent(in) :: long_name
       character(*), intent(in) :: short_name
       character(*), intent(in) :: units
       type(Variable) :: var
       character(len=:), allocatable :: var_name
       if (MAPL_AM_I_Root()) then
         if( meta%has_variable(short_name)) return
         var_name = short_name
         var = variable(type=pFIO_REAL32, dimensions='id,time')
         call var%add_attribute('long_name', long_name)
         call var%add_attribute('units', units)
         call var%add_attribute('positive', "up")
         call var%add_attribute('_FillValue', -999.99)
         call var%add_attribute('missing_value', -999.99)
         call meta%add_variable(var_name, var)
         call formatter%add_variable(meta, short_name)
       endif
     end subroutine create_new_vars

end module
