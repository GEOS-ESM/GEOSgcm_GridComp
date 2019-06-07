!-*- F90 -*-
module GHOST_CUBSPH_mod
  !--------------------------------------------------------------------!
  ! author:  Michael Herzog                                            !
  ! email:   Michael.Herzog@noaa.gov                                   !
  ! date:    May 2006                                                  !
  ! version: 0.1                                                       !
  !                                                                    !
  ! ghost cell update for cubed sphere, serial version                 !
  !--------------------------------------------------------------------!

  use fv_arrays_mod, only : REAL4, REAL8

  implicit none

  integer, parameter :: A_grid=0,                                       &
                        B_grid=1

  private
  public :: A_grid, B_grid, ghost_cubsph_update

  interface ghost_cubsph_update
     module procedure ghost_cubsph_update_3d_r8
     module procedure ghost_cubsph_update_2d_r8
     module procedure ghost_cubsph_update_3d_r4
     module procedure ghost_cubsph_update_2d_r4
  end interface

contains

!======================================================================!
  subroutine ghost_cubsph_update_3d_r8(var, n1x, nx, n1y, ny, nz, ng, ntiles, &
                                 ng_update, nlev, tile, grid_type)
    !------------------------------------------------------------------!
    ! serial ghost cell update for cubed sphere for global array var   !
    ! update ng (inner) ghost cells                                    !
    !------------------------------------------------------------------!
    integer, intent(in) :: n1x, nx, n1y, ny, nz, ng, ntiles,            &
                                 ng_update, nlev, tile, grid_type

    real(REAL8), dimension(n1x:nx,n1y:ny,nz,ntiles), intent(inout) :: var

    integer :: k

    do k=1,nlev
       call ghost_cubsph_update(var(:,:,k,:), n1x, nx, n1y, ny, ng, ntiles, &
                                 ng_update, tile, grid_type)
    enddo

  end subroutine ghost_cubsph_update_3d_r8
  !====================================================================!

!======================================================================!
  subroutine ghost_cubsph_update_2d_r8(var, n1x, nx, n1y, ny, ng, ntiles, &
                                 ng_update, tile, grid_type)
    !------------------------------------------------------------------!
    ! serial ghost cell update for cubed sphere for global array var   !
    ! update ng (inner) ghost cells                                    !
    !------------------------------------------------------------------!
    integer, intent(in) :: n1x, nx, n1y, ny, ng, ntiles,            &
                                 ng_update, tile, grid_type

    real(REAL8), dimension(n1x:nx,n1y:ny,ntiles), intent(inout) :: var
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: i, j, i_in, j_in, i_out, j_out, ig,                   &
               n1x_comp, nx_comp, n1y_comp, ny_comp

    n1x_comp=n1x+ng
    nx_comp=nx-ng
    n1y_comp=n1y+ng
    ny_comp=ny-ng
    if (grid_type==B_grid) then
       select case (tile)

       case (1)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig
                   var(i_in,j_in,tile)=var(i_out,j_out,6)
                
                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,3)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig
                   var(i_in,j_in,tile)=var(i_out,j_out,5)
                
                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,2)
                enddo
             enddo
          
       case (2)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,6)
                
                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig
                   var(i_in,j_in,tile)=var(i_out,j_out,3)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,1)
                
                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig
                   var(i_in,j_in,tile)=var(i_out,j_out,4)
                enddo
             enddo
          
       case (3)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig
                   var(i_in,j_in,tile)=var(i_out,j_out,2)
                
                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,5)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig
                   var(i_in,j_in,tile)=var(i_out,j_out,1)
                
                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,4)
                enddo
             enddo
          
       case (4)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,2)
                
                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig
                   var(i_in,j_in,tile)=var(i_out,j_out,5)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,3)
                   
                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig
                   var(i_in,j_in,tile)=var(i_out,j_out,6)
                enddo
             enddo
          
       case (5)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig
                   var(i_in,j_in,tile)=var(i_out,j_out,4)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,1)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig
                   var(i_in,j_in,tile)=var(i_out,j_out,3)
                
                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,6)
                enddo
             enddo
          
       case (6)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,4)
                
                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig
                   var(i_in,j_in,tile)=var(i_out,j_out,1)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,5)
                
                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig
                   var(i_in,j_in,tile)=var(i_out,j_out,2)
                enddo
             enddo
       end select

    elseif (grid_type==A_grid) then

       select case (tile)

       case (1)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig+1
                   var(i_in,j_in,tile)=var(i_out,j_out,6)
                
                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig-1
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,3)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig+1
                   var(i_in,j_in,tile)=var(i_out,j_out,5)
                
                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig-1
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,2)
                enddo
             enddo
          
       case (2)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig+1
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,6)
                
                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in,tile)=var(i_out,j_out,3)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig+1
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,1)
                
                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in,tile)=var(i_out,j_out,4)
                enddo
             enddo
          
       case (3)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig+1
                   var(i_in,j_in,tile)=var(i_out,j_out,2)
                
                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig-1
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,5)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig+1
                   var(i_in,j_in,tile)=var(i_out,j_out,1)
                
                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig-1
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,4)
                enddo
             enddo
          
       case (4)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig+1
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,2)
                
                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in,tile)=var(i_out,j_out,5)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig+1
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,3)
                
                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in,tile)=var(i_out,j_out,6)
                enddo
             enddo
          
       case (5)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig+1
                   var(i_in,j_in,tile)=var(i_out,j_out,4)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig-1
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,1)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig+1
                   var(i_in,j_in,tile)=var(i_out,j_out,3)
                
                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig-1
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,6)
                enddo
             enddo
          
       case (6)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig+1
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,4)
                
                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in,tile)=var(i_out,j_out,1)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig+1
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,5)
                   
                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in,tile)=var(i_out,j_out,2)
                enddo
             enddo
       end select
    else
       stop 'grid_type incorrect'
    end if

  end subroutine ghost_cubsph_update_2d_r8
  !====================================================================!

!======================================================================!
  subroutine ghost_cubsph_update_3d_r4(var, n1x, nx, n1y, ny, nz, ng, ntiles, &
                                 ng_update, nlev, tile, grid_type)
    !------------------------------------------------------------------!
    ! serial ghost cell update for cubed sphere for global array var   !
    ! update ng (inner) ghost cells                                    !
    !------------------------------------------------------------------!
    integer, intent(in) :: n1x, nx, n1y, ny, nz, ng, ntiles,            &
                                 ng_update, nlev, tile, grid_type

    real(REAL4), dimension(n1x:nx,n1y:ny,nz,ntiles), intent(inout) :: var

    integer :: k

    do k=1,nlev
       call ghost_cubsph_update(var(:,:,k,:), n1x, nx, n1y, ny, ng, ntiles, &
                                 ng_update, tile, grid_type)
    enddo

  end subroutine ghost_cubsph_update_3d_r4
  !====================================================================!

  subroutine ghost_cubsph_update_2d_r4(var, n1x, nx, n1y, ny, ng, ntiles, &
                                 ng_update, tile, grid_type)
    !------------------------------------------------------------------!
    ! serial ghost cell update for cubed sphere for global array var   !
    ! update ng (inner) ghost cells                                    !
    !------------------------------------------------------------------!
    integer, intent(in) :: n1x, nx, n1y, ny, ng, ntiles,            &
                                 ng_update, tile, grid_type

    real(REAL4), dimension(n1x:nx,n1y:ny,ntiles), intent(inout) :: var
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: i, j, i_in, j_in, i_out, j_out, ig,                   &
               n1x_comp, nx_comp, n1y_comp, ny_comp

    n1x_comp=n1x+ng
    nx_comp=nx-ng
    n1y_comp=n1y+ng
    ny_comp=ny-ng
    if (grid_type==B_grid) then
       select case (tile)

       case (1)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig
                   var(i_in,j_in,tile)=var(i_out,j_out,6)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,3)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig
                   var(i_in,j_in,tile)=var(i_out,j_out,5)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,2)
                enddo
             enddo

       case (2)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,6)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig
                   var(i_in,j_in,tile)=var(i_out,j_out,3)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,1)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig
                   var(i_in,j_in,tile)=var(i_out,j_out,4)
                enddo
             enddo

       case (3)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig
                   var(i_in,j_in,tile)=var(i_out,j_out,2)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,5)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig
                   var(i_in,j_in,tile)=var(i_out,j_out,1)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,4)
                enddo
             enddo

       case (4)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,2)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig
                   var(i_in,j_in,tile)=var(i_out,j_out,5)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,3)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig
                   var(i_in,j_in,tile)=var(i_out,j_out,6)
                enddo
             enddo

       case (5)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig
                   var(i_in,j_in,tile)=var(i_out,j_out,4)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,1)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig
                   var(i_in,j_in,tile)=var(i_out,j_out,3)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,6)
                enddo
             enddo

       case (6)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,4)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig
                   var(i_in,j_in,tile)=var(i_out,j_out,1)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,5)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig
                   var(i_in,j_in,tile)=var(i_out,j_out,2)
                enddo
             enddo
       end select

    elseif (grid_type==A_grid) then

       select case (tile)

       case (1)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig+1
                   var(i_in,j_in,tile)=var(i_out,j_out,6)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig-1
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,3)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig+1
                   var(i_in,j_in,tile)=var(i_out,j_out,5)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig-1
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,2)
                enddo
             enddo

       case (2)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig+1
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,6)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in,tile)=var(i_out,j_out,3)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig+1
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,1)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in,tile)=var(i_out,j_out,4)
                enddo
             enddo

       case (3)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig+1
                   var(i_in,j_in,tile)=var(i_out,j_out,2)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig-1
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,5)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig+1
                   var(i_in,j_in,tile)=var(i_out,j_out,1)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig-1
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,4)
                enddo
             enddo

       case (4)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig+1
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,2)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in,tile)=var(i_out,j_out,5)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig+1
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,3)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in,tile)=var(i_out,j_out,6)
                enddo
             enddo

       case (5)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=i
                   j_out=ny_comp-ig+1
                   var(i_in,j_in,tile)=var(i_out,j_out,4)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=n1x_comp+ig-1
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,1)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=ny_comp-ig+1
                   var(i_in,j_in,tile)=var(i_out,j_out,3)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=n1x_comp+ig-1
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,6)
                enddo
             enddo

       case (6)
             do i=n1x_comp,nx_comp
                do ig=1,ng_update
                   i_in =i
                   j_in =n1y_comp-ig
                   i_out=nx_comp-ig+1
                   j_out=ny_comp-i+1
                   var(i_in,j_in,tile)=var(i_out,j_out,4)

                   i_in =i
                   j_in =ny_comp+ig
                   i_out=i
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in,tile)=var(i_out,j_out,1)
                enddo
             enddo
             do j=n1y_comp,ny_comp
                do ig=1,ng_update
                   i_in =n1x_comp-ig
                   j_in =j
                   i_out=nx_comp-ig+1
                   j_out=j
                   var(i_in,j_in,tile)=var(i_out,j_out,5)

                   i_in =nx_comp+ig
                   j_in =j
                   i_out=nx_comp-j+1
                   j_out=n1y_comp+ig-1
                   var(i_in,j_in,tile)=var(i_out,j_out,2)
                enddo
             enddo
       end select
    else
       stop 'grid_type incorrect'
    end if

  end subroutine ghost_cubsph_update_2d_r4
  !====================================================================!

end module GHOST_CUBSPH_mod
