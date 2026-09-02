#include "MAPL_Generic.h"

module LisFieldsMod

   !DESCRIPTION:
   ! Import/Export/Internal field catalog for the LISF Plug, built directly
   ! from LisFieldCatalogMod's LIS_FieldList, whose entries carry
   ! role classification baked in at declaration (mirroring LIS_NUOPC_Gluecode's
   ! LIS_HookupInit), so this Plug (which never goes through LISF's own NUOPC
   ! coupling layer) doesn't need to link LISF's NUOPC cap or
   ! LIS_FORC_AttributesMod just to get this static data. Call
   ! LisFieldsInit() to get the import/export/internal field specs.
   ! Classification:
   !   - INTERNAL: directConn (LSM-internal state fields that bypass generic
   !     exchange; wiring these into Run is LSM-specific -- see
   !     LIS_CopyToNoah_3_3/NoahMP_* in LIS_NUOPC_DataCopy.F90).
   !   - IMPORT: adImport and not directConn.
   !   - EXPORT: adExport and not directConn.
   ! SHORT_NAME is the field's LIS_FieldList "stateName", uppercased.

   use MAPL_ErrorHandlingMod, only: MAPL_RTRN
   use LisFieldCatalogMod, only: LIS_Field, LIS_FieldList

   implicit none
   private

   public :: LIS_FieldSpec
   public :: LisFieldsInit

   type :: LIS_FieldSpec
      character(len=:), allocatable :: short_name
      character(len=:), allocatable :: lis_name
      character(len=:), allocatable :: long_name
      character(len=:), allocatable :: units
   end type LIS_FieldSpec

contains

   ! Sorts every LIS_FieldList entry by role into import/export/internal fields.
   subroutine LisFieldsInit(import_fields, export_fields, internal_fields, rc)
      type(LIS_FieldSpec), allocatable, intent(out) :: import_fields(:)
      type(LIS_FieldSpec), allocatable, intent(out) :: export_fields(:)
      type(LIS_FieldSpec), allocatable, intent(out) :: internal_fields(:)
      integer, optional, intent(out) :: rc

      character(len=*), parameter :: Iam = 'LisFieldsInit'
      type(LIS_FieldList) :: catalog
      integer :: iter
      integer :: num_imports, num_exports, num_internals, catalog_size

      catalog = LIS_FieldList()
      catalog_size = size(catalog%fields)

      num_imports = 0
      num_exports = 0
      num_internals = 0
      do iter = 1, catalog_size
         if (catalog%fields(iter)%directConn) then
            num_internals = num_internals + 1
         else
            if (catalog%fields(iter)%adImport) num_imports = num_imports + 1
            if (catalog%fields(iter)%adExport) num_exports = num_exports + 1
         end if
      end do

      allocate(import_fields(num_imports))
      allocate(export_fields(num_exports))
      allocate(internal_fields(num_internals))

      num_imports = 0
      num_exports = 0
      num_internals = 0
      do iter = 1, catalog_size
         if (catalog%fields(iter)%directConn) then
            num_internals = num_internals + 1
            call FillSpec(catalog%fields(iter), internal_fields(num_internals))
         else
            if (catalog%fields(iter)%adImport) then
               num_imports = num_imports + 1
               call FillSpec(catalog%fields(iter), import_fields(num_imports))
            end if
            if (catalog%fields(iter)%adExport) then
               num_exports = num_exports + 1
               call FillSpec(catalog%fields(iter), export_fields(num_exports))
            end if
         end if
      end do

      RETURN_(_SUCCESS)
   end subroutine LisFieldsInit

   subroutine FillSpec(field, spec)
      type(LIS_Field), intent(in) :: field
      type(LIS_FieldSpec), intent(out) :: spec

      spec%lis_name = trim(field%stateName)
      spec%short_name = ToUpper(trim(field%stateName))
      spec%long_name = trim(field%stdName)
      spec%units = trim(field%units)
   end subroutine FillSpec

   function ToUpper(str) result(upper)
      character(len=*), intent(in) :: str
      character(len=len(str)) :: upper

      integer :: i, c

      do i = 1, len(str)
         c = iachar(str(i:i))
         if (c >= iachar('a') .and. c <= iachar('z')) then
            upper(i:i) = achar(c - 32)
         else
            upper(i:i) = str(i:i)
         end if
      end do
   end function ToUpper

end module LisFieldsMod
