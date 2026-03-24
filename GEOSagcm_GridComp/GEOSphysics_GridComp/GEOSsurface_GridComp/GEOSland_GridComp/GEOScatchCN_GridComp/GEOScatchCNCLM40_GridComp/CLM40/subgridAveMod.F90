module subgridAveMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: subgridAveMod
!
! !DESCRIPTION:
! Utilities to perfrom subgrid averaging
!
! !USES:
  use clmtype
  use clm_varcon, only : spval
! use abortutils, only : endrun

! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: p2c   ! Perfrom an average from pfts to columns

  interface p2c
     module procedure p2c_1d
     module procedure p2c_2d
     module procedure p2c_1d_filter
     module procedure p2c_2d_filter
  end interface
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_1d
!
! !INTERFACE:
  subroutine p2c_1d (lbp, ubp, lbc, ubc, parr, carr, p2c_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from pfts to columns.
! Averaging is only done for points that are not equal to "spval".
!
! !USES:
    use clm_varpar, only : max_pft_per_col
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp              ! beginning and ending pft
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column
    real, intent(in)  :: parr(lbp:ubp)         ! pft array
    real, intent(out) :: carr(lbc:ubc)         ! column array
    character(len=*), intent(in) :: p2c_scale_type ! scale type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: pi,p,c,index           ! indices
    real :: scale_p2c(lbp:ubp)     ! scale factor for column->landunit mapping
    logical  :: found                  ! temporary for error check
    real :: sumwt(lbc:ubc)         ! sum of weights
    real, pointer :: wtcol(:)      ! weight of pft relative to column
    integer , pointer :: pcolumn(:)    ! column index of corresponding pft
    integer , pointer :: npfts(:)      ! number of pfts in column
    integer , pointer :: pfti(:)       ! initial pft index in column
!------------------------------------------------------------------------

    wtcol    => clm3%g%l%c%p%wtcol
    pcolumn  => clm3%g%l%c%p%column
    npfts    => clm3%g%l%c%npfts
    pfti     => clm3%g%l%c%pfti

    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0
       end do
    else
       write(6,*)'p2c_1d error: scale type ',p2c_scale_type,' not supported'
       stop 12 ! call endrun()
    end if

    carr(lbc:ubc) = spval
    sumwt(lbc:ubc) = 0.
    do p = lbp,ubp
       if (wtcol(p) /= 0.) then
          if (parr(p) /= spval) then
             c = pcolumn(p)
             if (sumwt(c) == 0.) carr(c) = 0.
             carr(c) = carr(c) + parr(p) * scale_p2c(p) * wtcol(p)
             sumwt(c) = sumwt(c) + wtcol(p)
          end if
       end if
    end do
    found = .false.
    do c = lbc,ubc
       if (sumwt(c) > 1.0 + 1.e-6) then
          found = .true.
          index = c
       else if (sumwt(c) /= 0.) then
          carr(c) = carr(c)/sumwt(c)
       end if
    end do
    if (found) then
       write(6,*)'p2c error: sumwt is greater than 1.0 at c= ',index
       stop 13 ! call endrun()
    end if

  end subroutine p2c_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_2d
!
! !INTERFACE:
  subroutine p2c_2d (lbp, ubp, lbc, ubc, num2d, parr, carr, p2c_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from landunits to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !USES:
    use clm_varpar, only : max_pft_per_col
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp              ! beginning and ending pft
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column
    integer , intent(in)  :: num2d                 ! size of second dimension
    real, intent(in)  :: parr(lbp:ubp,num2d)   ! pft array
    real, intent(out) :: carr(lbc:ubc,num2d)   ! column array
    character(len=*), intent(in) :: p2c_scale_type ! scale type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: j,pi,p,c,index         ! indices
    real :: scale_p2c(lbp:ubp)     ! scale factor for column->landunit mapping
    logical  :: found                  ! temporary for error check
    real :: sumwt(lbc:ubc)         ! sum of weights
    real, pointer :: wtcol(:)      ! weight of pft relative to column
    integer , pointer :: pcolumn(:)    ! column index of corresponding pft
    integer , pointer :: npfts(:)      ! number of pfts in column
    integer , pointer :: pfti(:)       ! initial pft index in column
!------------------------------------------------------------------------

    wtcol    => clm3%g%l%c%p%wtcol
    pcolumn  => clm3%g%l%c%p%column
    npfts    => clm3%g%l%c%npfts
    pfti     => clm3%g%l%c%pfti

    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0
       end do
    else
       write(6,*)'p2c_2d error: scale type ',p2c_scale_type,' not supported'
       stop 12 ! call endrun()
    end if

    carr(:,:) = spval
    do j = 1,num2d
       sumwt(:) = 0.
       do p = lbp,ubp
          if (wtcol(p) /= 0.) then
             if (parr(p,j) /= spval) then
                c = pcolumn(p)
                if (sumwt(c) == 0.) carr(c,j) = 0.
                carr(c,j) = carr(c,j) + parr(p,j) * scale_p2c(p) * wtcol(p)
                sumwt(c) = sumwt(c) + wtcol(p)
             end if
          end if
       end do
       found = .false.
       do c = lbc,ubc
          if (sumwt(c) > 1.0 + 1.e-6) then
             found = .true.
             index = c
          else if (sumwt(c) /= 0.) then
             carr(c,j) = carr(c,j)/sumwt(c)
          end if
       end do
       if (found) then
          write(6,*)'p2c_2d error: sumwt is greater than 1.0 at c= ',index,' lev= ',j
          stop 13 ! call endrun()
       end if
    end do 
  end subroutine p2c_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_1d_filter
!
! !INTERFACE:
  subroutine p2c_1d_filter (numfc, filterc, pftarr, colarr)
!
! !DESCRIPTION:
! perform pft to column averaging for single level pft arrays
!
! !USES:
    use clm_varpar, only : max_pft_per_col
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: numfc
    integer , intent(in)  :: filterc(numfc)
    real, pointer     :: pftarr(:)
    real, pointer     :: colarr(:)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: fc,c,pi,p           ! indices
    integer , pointer :: npfts(:)
    integer , pointer :: pfti(:)
    integer , pointer :: pftf(:)
    real, pointer :: wtcol(:)
    real, pointer :: wtgcell(:)
!-----------------------------------------------------------------------

    npfts   => clm3%g%l%c%npfts
    pfti    => clm3%g%l%c%pfti
    pftf    => clm3%g%l%c%pftf
    wtcol   => clm3%g%l%c%p%wtcol
    wtgcell => clm3%g%l%c%p%wtgcell

    do fc = 1,numfc
       c = filterc(fc)
       colarr(c) = 0.
       do p = pfti(c), pftf(c)
          if (wtgcell(p) > 0.) colarr(c) = colarr(c) + pftarr(p) * wtcol(p)
       end do
    end do

  end subroutine p2c_1d_filter

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_2d_filter
!
! !INTERFACE:
  subroutine p2c_2d_filter (lev, numfc, filterc, pftarr, colarr)
!
! !DESCRIPTION:
! perform pft to column averaging for multi level pft arrays
!
! !USES:
    use clm_varpar, only : max_pft_per_col

! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lev
    integer , intent(in)  :: numfc
    integer , intent(in)  :: filterc(numfc)
    real, pointer     :: pftarr(:,:)
    real, pointer     :: colarr(:,:)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: fc,c,pi,p,j    ! indices
    integer , pointer :: npfts(:)
    integer , pointer :: pfti(:)
    integer , pointer :: pftf(:)
    real, pointer :: wtcol(:)
!-----------------------------------------------------------------------

    npfts => clm3%g%l%c%npfts
    pfti  => clm3%g%l%c%pfti
    pftf  => clm3%g%l%c%pftf
    wtcol => clm3%g%l%c%p%wtcol

    do j = 1,lev
       do fc = 1,numfc
          c = filterc(fc)
          colarr(c,j) = 0.
          do p = pfti(c), pftf(c)
             colarr(c,j) = colarr(c,j) + pftarr(p,j) * wtcol(p)
          end do
       end do
    end do

  end subroutine p2c_2d_filter

end module subgridAveMod
