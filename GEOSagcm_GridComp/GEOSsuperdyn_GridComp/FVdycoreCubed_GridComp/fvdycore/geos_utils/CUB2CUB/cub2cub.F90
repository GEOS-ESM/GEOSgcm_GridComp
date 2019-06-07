!-*- F90 -*-
module CUB2CUB_mod
  !--------------------------------------------------------------------!
  ! author:  Michael Herzog                                            !
  ! email:   Michael.Herzog@noaa.gov                                   !
  ! date:    Feb 2007                                                  !
  ! version: 0.1                                                       !
  !                                                                    !
  ! routines for interpolation between two cubed sphere grids          !
  ! with different spatial resolution                                  !
  !--------------------------------------------------------------------!
  use fv_arrays_mod,  only: REAL4, REAL8, FVPRC
  implicit none

  private
  public :: read_c2c_namelist,                                          &
            get_c2c_weight,                                             &
            interpolate_c2c


  real(FVPRC), parameter :: pi = 3.141592653589793

  real(FVPRC), dimension(:,:,:,:), allocatable :: xyz_corner_in, xyz_corner_out

  logical :: west_edge  = .true., east_edge  = .true.,                  &
             south_edge = .true., north_edge = .true.
  logical :: sw_corner = .true., se_corner = .true.,                    &
             nw_corner = .true., ne_corner = .true.
  logical :: edge_interp = .false.

contains
!======================================================================!
  subroutine read_c2c_namelist(ntiles, grid_in, grid_out, dir_in, dir_out)
    !------------------------------------------------------------------!
    ! read namelist files                                              !
    !------------------------------------------------------------------!
    integer, intent(out) :: ntiles
    character(len=120), intent(out) :: grid_in , grid_out, dir_in, dir_out
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: ios, l
    logical :: exists
    character(len=80) :: filename

    namelist /c2c_nml/ ntiles, grid_in , grid_out, dir_in, dir_out
    !------------------------------------------------------------------!
    ! read c2c.nml namelist file                                       !
    !------------------------------------------------------------------!
    filename = "c2c.nml"
    inquire(file=trim(filename),exist=exists)

    if (.not. exists) then
       write(6,100) trim(filename)
100    format (/,"namelist file ",a," doesn't exist",/)
       stop
    else
       ntiles=6
       grid_in="grid_spec"
       grid_out="grid_spec"
       dir_in="IN"
       dir_out="OUT"

       open (10,file=filename)
       read (10,nml=c2c_nml,iostat=ios)
       close(10)
       if (ios > 0) then
          write(6,101) trim(filename), ios
101       format(/,"c2c_nml ERROR: reading ",a,", iostat=",i4,/)
          stop
       endif

       if (trim(dir_in)==trim(dir_out)) then
          write(6,110) trim(dir_in)
110       format (/,"FATAL: input/output directories are identical: ",a,/)
          stop
       endif

       grid_in =trim(dir_in )//"/"//trim(grid_in )
       filename=trim(grid_in)//".tile1.nc"
       inquire(file=trim(filename),exist=exists)
       if (.not. exists) then
          write(6,120) trim(grid_in)
120       format (/,"input grid_spec file ",a," doesn't exist",/)
          stop
       endif

       grid_out=trim(dir_out)//"/"//trim(grid_out)
       filename=trim(grid_out)//".tile1.nc"
       inquire(file=trim(filename),exist=exists)
       if (.not. exists) then
          write(6,120) trim(grid_out)
130       format (/,"output grid_spec file ",a," doesn't exist",/)
          stop
       endif

    endif

  end subroutine read_c2c_namelist
  !====================================================================!
  subroutine init_cubsph_sgrid(npx, npy, ntiles, grid_file, super_grid)
    !------------------------------------------------------------------!
    ! read/calculate supergrid                                         !
    !------------------------------------------------------------------!
    use CUB2LATLON_mod,  only: read_netcdf_grid
    use GRID_UTILS_mod, only: latlon2xyz

    integer, intent(in) :: npx, npy, ntiles
    character(len=120), intent(in) :: grid_file
    real(FVPRC), dimension(3,2*(npx-1)+1,2*(npy-1)+1,ntiles), intent(out) :: super_grid
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(FVPRC) :: abs_vector
    real(FVPRC), dimension(2,1:npx,1:npy,ntiles) :: sph_corner
    integer :: i, j, l, is, js
    !------------------------------------------------------------------!
    ! read sph_corner from file                                        !
    ! write to super_grid as xyz                                       !
    !------------------------------------------------------------------!
    call read_netcdf_grid(npx, npy, ntiles, grid_file, sph_corner)
    do l=1,ntiles
       do j=1,npy
          js=2*(j-1)+1
          do i=1,npx
             is=2*(i-1)+1
             call latlon2xyz(sph_corner(:,i,j,l), super_grid(:,is,js,l))
          enddo
       enddo
    enddo
    !------------------------------------------------------------------!
    ! calculate agrid                                                  !
    !------------------------------------------------------------------!
    do l=1,ntiles
       do j=1,2*(npy-1)-1,2
          do i=1,2*(npx-1)-1,2
             super_grid(:,i+1,j+1,l)=0.25*(super_grid(:,i  ,j  ,l)      &
                                          +super_grid(:,i+2,j  ,l)      &
                                          +super_grid(:,i+2,j+2,l)      &
                                          +super_grid(:,i  ,j+2,l))
             abs_vector=super_grid(1,i+1,j+1,l)*super_grid(1,i+1,j+1,l) &
                       +super_grid(2,i+1,j+1,l)*super_grid(2,i+1,j+1,l) &
                       +super_grid(3,i+1,j+1,l)*super_grid(3,i+1,j+1,l)
             if (abs_vector>0.) super_grid(:,i+1,j+1,l)=super_grid(:,i+1,j+1,l)/sqrt(abs_vector)
          enddo
       enddo
    enddo
    !------------------------------------------------------------------!
    ! calculate dgridu                                                 !
    !------------------------------------------------------------------!
    do l=1,ntiles
       do j=1,2*(npy-1)-1,2
          do i=1,2*(npx-1),2
             super_grid(:,i,j+1,l)=0.5*(super_grid(:,i  ,j  ,l)      &
                                       +super_grid(:,i  ,j+2,l))
             abs_vector=super_grid(1,i,j+1,l)*super_grid(1,i,j+1,l) &
                       +super_grid(2,i,j+1,l)*super_grid(2,i,j+1,l) &
                       +super_grid(3,i,j+1,l)*super_grid(3,i,j+1,l)
             if (abs_vector>0.) super_grid(:,i,j+1,l)=super_grid(:,i,j+1,l)/sqrt(abs_vector)
          enddo
       enddo
    enddo
    !------------------------------------------------------------------!
    ! calculate dgridv                                                 !
    !------------------------------------------------------------------!
    do l=1,ntiles
       do j=1,2*(npy-1),2
          do i=1,2*(npx-1)-1,2
             super_grid(:,i+1,j,l)=0.5*(super_grid(:,i  ,j,l)      &
                                       +super_grid(:,i+2,j,l))
             abs_vector=super_grid(1,i+1,j,l)*super_grid(1,i+1,j,l) &
                       +super_grid(2,i+1,j,l)*super_grid(2,i+1,j,l) &
                       +super_grid(3,i+1,j,l)*super_grid(3,i+1,j,l)
             if (abs_vector>0.) super_grid(:,i+1,j,l)=super_grid(:,i+1,j,l)/sqrt(abs_vector)
          enddo
       enddo
    enddo

  end subroutine init_cubsph_sgrid
  !====================================================================!
  subroutine get_c2c_weight(ntiles, npx_in, npy_in, sph_in,             &
                            npx_out, npy_out, sph_out,                  &
                            index_c2c,  weight_c2c)
    !------------------------------------------------------------------!
    ! calculate indices and weight for cub2cub interpolation           !
    !------------------------------------------------------------------!
    use GRID_UTILS_mod, only: latlon2xyz,      great_circle,            &
                              spherical_angle, dist2side

    integer, intent(in) :: ntiles, npx_in, npy_in, npx_out, npy_out
    real(FVPRC), intent(in) :: sph_in (2, 0:npx_in +1, 0:npy_in +1, ntiles),   &
                        sph_out(2, 0:npx_out+1, 0:npy_out+1, ntiles)
    integer, intent(out) :: index_c2c(3, npx_out-1, npy_out-1, ntiles)
    real(FVPRC),    intent(out) :: weight_c2c (4, npx_out-1, npy_out-1, ntiles)
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real(FVPRC) :: xyz_center_in (3, 0:npx_in  ,  0:npy_in  ,  ntiles),        &
            xyz_center_out(3, 0:npx_out  , 0:npy_out  , ntiles)
    real(FVPRC) :: abs_center, distance, shortest, sum, dist1, dist2, dist3, dist4
    integer :: i, j, l, ic, jc, i_guess, j_guess, lc, iter, index(3)
    logical :: found
    !------------------------------------------------------------------!
    ! calculate xyz cell corners and cell centers                      !
    !------------------------------------------------------------------!
    allocate(xyz_corner_in (3, 0:npx_in+1,  0:npy_in+1,  ntiles),        &
             xyz_corner_out(3, 0:npx_out+1, 0:npy_out+1, ntiles))
    do l=1,ntiles
       do j=0,npy_in+1
          do i=0,npx_in+1
            call latlon2xyz(sph_in(:,i,j,l), xyz_corner_in(:,i,j,l))
          enddo
       enddo
       do j=0,npx_in
          do i=0,npy_in
             xyz_center_in(:,i,j,l)=0.25*(xyz_corner_in(:,i  ,j  ,l)+xyz_corner_in(:,i+1,j  ,l)  &
                                         +xyz_corner_in(:,i  ,j+1,l)+xyz_corner_in(:,i+1,j+1,l))
             abs_center=xyz_center_in(1,i,j,l)*xyz_center_in(1,i,j,l)   &
                       +xyz_center_in(2,i,j,l)*xyz_center_in(2,i,j,l)   &
                       +xyz_center_in(3,i,j,l)*xyz_center_in(3,i,j,l)
             if (abs_center>0.) xyz_center_in(:,i,j,l)=xyz_center_in(:,i,j,l)/sqrt(abs_center)
          enddo
       enddo
    enddo

    do l=1,ntiles
       do j=0,npy_out+1
          do i=0,npx_out+1
            call latlon2xyz(sph_out(:,i,j,l), xyz_corner_out(:,i,j,l))
          enddo
       enddo
       do j=0,npx_out
          do i=0,npy_out
             xyz_center_out(:,i,j,l)=0.25*(xyz_corner_out(:,i  ,j  ,l)+xyz_corner_out(:,i+1,j  ,l)  &
                                          +xyz_corner_out(:,i  ,j+1,l)+xyz_corner_out(:,i+1,j+1,l))
             abs_center=xyz_center_out(1,i,j,l)*xyz_center_out(1,i,j,l) &
                       +xyz_center_out(2,i,j,l)*xyz_center_out(2,i,j,l) &
                       +xyz_center_out(3,i,j,l)*xyz_center_out(3,i,j,l)
             if (abs_center>0.) xyz_center_out(:,i,j,l)=xyz_center_out(:,i,j,l)/sqrt(abs_center)
          enddo
       enddo
    enddo
    !------------------------------------------------------------------!
    ! find lower left corner on grid_in for given grid_out location    !
    !                                                                  !
    ! NOTE: tile corner locations must match                           !
    !------------------------------------------------------------------!
    do l=1,ntiles
       do j=1,npy_out-1
          j_guess=int(real(j)*real(npy_in)/real(npy_out))
          do i=1,npx_out-1
             i_guess=int(real(i)*real(npx_in)/real(npx_out))

             found=.false.
             iter_loop: do iter=1,10
                if (.not. found) then
                   shortest=pi+pi
                   do jc=max(1, j_guess-iter+1), min(npy_in-1, max(1,j_guess+iter-1))
                      do ic=max(1, i_guess-iter+1), min(npx_in-1, max(1,i_guess+iter-1))
                         distance=great_circle(xyz_center_in(:,ic,jc,l),xyz_center_out(:,i,j,l))
                         if (distance<shortest) then
                            shortest=distance
                            index(1)=ic
                            index(2)=jc
                            index(3)=l
                         endif
                      enddo
                   enddo
                   call set_closest_index(index(1), index(2), index(3), found)
                   if (found) exit iter_loop
                endif
             enddo iter_loop
             if (.not. found) then
                write(6, 100) i, j, l
100             format ("FATAL: index search for ", 3i4," failed")
                stop
             endif
          enddo
       enddo
    enddo
    !------------------------------------------------------------------!
    ! calculate weights for interpolation:                             !
    ! calculate shortest distance to each side of rectangle            !
    ! formed by cubed sphere cell centers                              !
    ! special corner treatment                                         !
    !------------------------------------------------------------------!
    do l=1,ntiles
       do j=1,npy_out-1
          do i=1,npx_out-1
             ic=index_c2c(1,i,j,l)
             jc=index_c2c(2,i,j,l)
             lc=index_c2c(3,i,j,l)    ! l==lc!
             if (ic==npx_in-1 .and. jc==npy_in-1) then
                !---------------------------------------------------------!
                ! calculate weights for bilinear interpolation near corner!
                !---------------------------------------------------------!
                dist1=dist2side(xyz_center_in(:,ic+1,jc,l),xyz_center_in(:,ic,jc+1,l),xyz_center_out(:,i,j,l))
                dist2=dist2side(xyz_center_in(:,ic+1,jc,l),xyz_center_in(:,ic,jc  ,l),xyz_center_out(:,i,j,l))
                dist3=dist2side(xyz_center_in(:,ic  ,jc,l),xyz_center_in(:,ic,jc+1,l),xyz_center_out(:,i,j,l))
                
                weight_c2c(1,i,j,l)=dist1      ! ic,   jc    weight
                weight_c2c(2,i,j,l)=dist2      ! ic,   jc+1  weight
                weight_c2c(3,i,j,l)=0.         ! ic+1, jc+1  weight
                weight_c2c(4,i,j,l)=dist3      ! ic+1, jc    weight
                
                sum=weight_c2c(1,i,j,l)+weight_c2c(2,i,j,l)+weight_c2c(4,i,j,l)
                weight_c2c(1,i,j,l)=weight_c2c(1,i,j,l)/sum
                weight_c2c(2,i,j,l)=weight_c2c(2,i,j,l)/sum
                weight_c2c(4,i,j,l)=weight_c2c(4,i,j,l)/sum
                
             elseif (ic==0 .and. jc==npy_in-1) then
                !---------------------------------------------------------!
                ! calculate weights for bilinear interpolation near corner!
                !---------------------------------------------------------!
                dist1=dist2side(xyz_center_in(:,ic+1,jc+1,l),xyz_center_in(:,ic+1,jc,l),xyz_center_out(:,i,j,l))
                dist2=dist2side(xyz_center_in(:,ic+1,jc  ,l),xyz_center_in(:,ic  ,jc,l),xyz_center_out(:,i,j,l))
                dist3=dist2side(xyz_center_in(:,ic+1,jc+1,l),xyz_center_in(:,ic  ,jc,l),xyz_center_out(:,i,j,l))
                
                weight_c2c(1,i,j,l)=dist1      ! ic,   jc    weight
                weight_c2c(2,i,j,l)=0.         ! ic,   jc+1  weight
                weight_c2c(3,i,j,l)=dist2      ! ic+1, jc+1  weight
                weight_c2c(4,i,j,l)=dist3      ! ic+1, jc    weight
                
                sum=weight_c2c(1,i,j,l)+weight_c2c(3,i,j,l)+weight_c2c(4,i,j,l)
                weight_c2c(1,i,j,l)=weight_c2c(1,i,j,l)/sum
                weight_c2c(3,i,j,l)=weight_c2c(3,i,j,l)/sum
                weight_c2c(4,i,j,l)=weight_c2c(4,i,j,l)/sum
                
             elseif (jc==0 .and. ic==npx_in-1) then
                !---------------------------------------------------------!
                ! calculate weights for bilinear interpolation near corner!
                !---------------------------------------------------------!
                dist1=dist2side(xyz_center_in(:,ic,jc+1,l),xyz_center_in(:,ic+1,jc+1,l),xyz_center_out(:,i,j,l))
                dist2=dist2side(xyz_center_in(:,ic,jc  ,l),xyz_center_in(:,ic+1,jc+1,l),xyz_center_out(:,i,j,l))
                dist3=dist2side(xyz_center_in(:,ic,jc  ,l),xyz_center_in(:,ic  ,jc+1,l),xyz_center_out(:,i,j,l))
                
                weight_c2c(1,i,j,l)=dist1      ! ic,   jc    weight
                weight_c2c(2,i,j,l)=dist2      ! ic,   jc+1  weight
                weight_c2c(3,i,j,l)=dist3      ! ic+1, jc+1  weight
                weight_c2c(4,i,j,l)=0.         ! ic+1, jc    weight
                
                sum=weight_c2c(1,i,j,l)+weight_c2c(2,i,j,l)+weight_c2c(3,i,j,l)
                weight_c2c(1,i,j,l)=weight_c2c(1,i,j,l)/sum
                weight_c2c(2,i,j,l)=weight_c2c(2,i,j,l)/sum
                weight_c2c(3,i,j,l)=weight_c2c(3,i,j,l)/sum
                
             elseif (jc==0 .and. ic==0) then
                !---------------------------------------------------------!
                ! calculate weights for bilinear interpolation near corner!
                !---------------------------------------------------------!
                dist2=dist2side(xyz_center_in(:,ic+1,jc  ,l),xyz_center_in(:,ic+1,jc+1,l),xyz_center_out(:,i,j,l))
                dist3=dist2side(xyz_center_in(:,ic  ,jc+1,l),xyz_center_in(:,ic+1,jc  ,l),xyz_center_out(:,i,j,l))
                dist4=dist2side(xyz_center_in(:,ic+1,jc+1,l),xyz_center_in(:,ic  ,jc+1,l),xyz_center_out(:,i,j,l))
                
                weight_c2c(1,i,j,l)=0.         ! ic,   jc    weight
                weight_c2c(2,i,j,l)=dist2      ! ic,   jc+1  weight
                weight_c2c(3,i,j,l)=dist3      ! ic+1, jc+1  weight
                weight_c2c(4,i,j,l)=dist4      ! ic+1, jc    weight
                
                sum=weight_c2c(2,i,j,l)+weight_c2c(3,i,j,l)+weight_c2c(4,i,j,l)
                weight_c2c(2,i,j,l)=weight_c2c(2,i,j,l)/sum
                weight_c2c(3,i,j,l)=weight_c2c(3,i,j,l)/sum
                weight_c2c(4,i,j,l)=weight_c2c(4,i,j,l)/sum
                
             else
                !----------------------------------------------------------!
                ! calculate weights for bilinear interpolation if no corner!
                !----------------------------------------------------------!
                dist1=dist2side(xyz_center_in(:,ic  ,jc  ,l),xyz_center_in(:,ic  ,jc+1,l),xyz_center_out(:,i,j,l))
                dist2=dist2side(xyz_center_in(:,ic  ,jc+1,l),xyz_center_in(:,ic+1,jc+1,l),xyz_center_out(:,i,j,l))
                dist3=dist2side(xyz_center_in(:,ic+1,jc+1,l),xyz_center_in(:,ic+1,jc  ,l),xyz_center_out(:,i,j,l))
                dist4=dist2side(xyz_center_in(:,ic+1,jc  ,l),xyz_center_in(:,ic  ,jc  ,l),xyz_center_out(:,i,j,l))
                
                weight_c2c(1,i,j,l)=dist2*dist3      ! ic,   jc    weight
                weight_c2c(2,i,j,l)=dist3*dist4      ! ic,   jc+1  weight
                weight_c2c(3,i,j,l)=dist4*dist1      ! ic+1, jc+1  weight
                weight_c2c(4,i,j,l)=dist1*dist2      ! ic+1, jc    weight
                
                sum=weight_c2c(1,i,j,l)+weight_c2c(2,i,j,l)+weight_c2c(3,i,j,l)+weight_c2c(4,i,j,l)
                weight_c2c(:,i,j,l)=weight_c2c(:,i,j,l)/sum
             endif
          enddo
       enddo
    enddo

  contains
    !------------------------------------------------------------------!
    subroutine set_closest_index(ig, jg, lg, ok)
      !----------------------------------------------------------------!
      ! determine lower left corner                                    !
      !----------------------------------------------------------------!
      integer, intent(in)  :: ig, jg, lg
      logical, intent(out) :: ok

      real(FVPRC) :: angle_1, angle_1a, angle_1b,                              & 
              angle_2, angle_2a, angle_2b,                              & 
              angle_3, angle_3a, angle_3b,                              &
              angle_4, angle_4a, angle_4b
      real(FVPRC) :: angle_11, angle_11a, angle_11b,                           & 
              angle_22, angle_22a, angle_22b,                           & 
              angle_33, angle_33a, angle_33b,                           &
              angle_44, angle_44a, angle_44b


      ok=.false.
      angle_1 =spherical_angle(xyz_center_in(:,ig,jg,lg),xyz_center_in(:,ig+1,jg  ,lg),xyz_center_in(:,ig  ,jg+1,lg))
      angle_1a=spherical_angle(xyz_center_in(:,ig,jg,lg),xyz_center_in(:,ig+1,jg  ,lg),xyz_center_out(:,i,j,lg))
      angle_1b=spherical_angle(xyz_center_in(:,ig,jg,lg),xyz_center_in(:,ig  ,jg+1,lg),xyz_center_out(:,i,j,lg))
      if (max(angle_1a,angle_1b)<=angle_1) then
          if (ig+1==npx_in .and. jg+1==npy_in) then
            angle_11 =spherical_angle(xyz_center_in(:,ig+1,jg,lg),xyz_center_in(:,ig,jg+1,lg),xyz_center_in(:,ig,jg,lg))
            angle_11a=spherical_angle(xyz_center_in(:,ig+1,jg,lg),xyz_center_in(:,ig,jg  ,lg),xyz_center_out(:,i,j,lg))
            angle_11b=spherical_angle(xyz_center_in(:,ig+1,jg,lg),xyz_center_in(:,ig,jg+1,lg),xyz_center_out(:,i,j,lg))
         else
            angle_11 =spherical_angle(xyz_center_in(:,ig+1,jg+1,lg),xyz_center_in(:,ig  ,jg+1,lg),xyz_center_in(:,ig+1,jg,lg))
            angle_11a=spherical_angle(xyz_center_in(:,ig+1,jg+1,lg),xyz_center_in(:,ig+1,jg  ,lg),xyz_center_out(:,i,j,lg))
            angle_11b=spherical_angle(xyz_center_in(:,ig+1,jg+1,lg),xyz_center_in(:,ig  ,jg+1,lg),xyz_center_out(:,i,j,lg))
         endif
         if (max(angle_11a,angle_11b)<=angle_11) then
            ok=.true.
            index_c2c(1,i,j,l)=ig
            index_c2c(2,i,j,l)=jg
            index_c2c(3,i,j,l)=lg
         endif
      else
         angle_2 =spherical_angle(xyz_center_in(:,ig,jg,lg),xyz_center_in(:,ig,jg+1,lg),xyz_center_in(:,ig-1,jg,lg))
         angle_2a=angle_1b
         angle_2b=spherical_angle(xyz_center_in(:,ig,jg,lg),xyz_center_in(:,ig-1,jg,lg),xyz_center_out(:,i,j,lg))
         if (max(angle_2a,angle_2b)<=angle_2) then
            if (ig-1==0 .and. jg+1==npy_in) then
               angle_22 =spherical_angle(xyz_center_in(:,ig,jg+1,lg),xyz_center_in(:,ig  ,jg,lg),xyz_center_in(:,ig-1,jg,lg))
               angle_22a=spherical_angle(xyz_center_in(:,ig,jg+1,lg),xyz_center_in(:,ig-1,jg  ,lg),xyz_center_out(:,i,j,lg))
               angle_22b=spherical_angle(xyz_center_in(:,ig,jg+1,lg),xyz_center_in(:,ig  ,jg,lg),xyz_center_out(:,i,j,lg))
            else
               angle_22 =spherical_angle(xyz_center_in(:,ig-1,jg+1,lg),xyz_center_in(:,ig  ,jg+1,lg),xyz_center_in(:,ig-1,jg,lg))
               angle_22a=spherical_angle(xyz_center_in(:,ig-1,jg+1,lg),xyz_center_in(:,ig-1,jg  ,lg),xyz_center_out(:,i,j,lg))
               angle_22b=spherical_angle(xyz_center_in(:,ig-1,jg+1,lg),xyz_center_in(:,ig  ,jg+1,lg),xyz_center_out(:,i,j,lg))
            endif
            if (max(angle_22a,angle_22b)<=angle_22) then
               ok=.true.
               index_c2c(1,i,j,l)=ig-1
               index_c2c(2,i,j,l)=jg
               index_c2c(3,i,j,l)=lg
            endif
         else
            angle_3 =spherical_angle(xyz_center_in(:,ig,jg,lg),xyz_center_in(:,ig-1,jg,lg),xyz_center_in(:,ig,jg-1,lg))
            angle_3a=angle_2b
            angle_3b=spherical_angle(xyz_center_in(:,ig,jg,lg),xyz_center_in(:,ig,jg-1,lg),xyz_center_out(:,i,j,lg))
            if (max(angle_3a,angle_3b)<=angle_3) then
               if (ig==1 .and. jg==1) then
                  angle_33 =spherical_angle(xyz_center_in(:,ig,jg-1,lg),xyz_center_in(:,ig  ,jg,lg),xyz_center_in(:,ig-1,jg,lg))
                  angle_33a=spherical_angle(xyz_center_in(:,ig,jg-1,lg),xyz_center_in(:,ig-1,jg  ,lg),xyz_center_out(:,i,j,lg))
                  angle_33b=spherical_angle(xyz_center_in(:,ig,jg-1,lg),xyz_center_in(:,ig  ,jg,lg),xyz_center_out(:,i,j,lg))
               else
                  angle_33 =spherical_angle(xyz_center_in(:,ig-1,jg-1,lg),xyz_center_in(:,ig  ,jg-1,lg),xyz_center_in(:,ig-1,jg,lg))
                  angle_33a=spherical_angle(xyz_center_in(:,ig-1,jg-1,lg),xyz_center_in(:,ig-1,jg  ,lg),xyz_center_out(:,i,j,lg))
                  angle_33b=spherical_angle(xyz_center_in(:,ig-1,jg-1,lg),xyz_center_in(:,ig  ,jg-1,lg),xyz_center_out(:,i,j,lg))
               endif
               if (max(angle_33a,angle_33b)<=angle_33) then
                  ok=.true.
                  index_c2c(1,i,j,l)=ig-1
                  index_c2c(2,i,j,l)=jg-1
                  index_c2c(3,i,j,l)=lg
               endif
            else
               angle_4 =spherical_angle(xyz_center_in(:,ig,jg,lg),xyz_center_in(:,ig,jg-1,lg),xyz_center_in(:,ig+1,jg,lg))
               angle_4a=angle_3b
               angle_4b=spherical_angle(xyz_center_in(:,ig,jg,lg),xyz_center_in(:,ig+1,jg,lg),xyz_center_out(:,i,j,l))
               if (max(angle_4a,angle_4b)<=angle_4) then
                  if (ig+1==npx_in .and. jg-1==0) then
                     angle_44 =spherical_angle(xyz_center_in(:,ig+1,jg,lg),xyz_center_in(:,ig,jg  ,lg),xyz_center_in(:,ig,jg-1,lg))
                     angle_44a=spherical_angle(xyz_center_in(:,ig+1,jg,lg),xyz_center_in(:,ig,jg-1,lg),xyz_center_out(:,i,j,lg))
                     angle_44b=spherical_angle(xyz_center_in(:,ig+1,jg,lg),xyz_center_in(:,ig,jg  ,lg),xyz_center_out(:,i,j,lg))
                  else
                     angle_44 =spherical_angle(xyz_center_in(:,ig+1,jg-1,lg),xyz_center_in(:,ig+1,jg  ,lg),xyz_center_in(:,ig,jg-1,lg))
                     angle_44a=spherical_angle(xyz_center_in(:,ig+1,jg-1,lg),xyz_center_in(:,ig  ,jg-1,lg),xyz_center_out(:,i,j,lg))
                     angle_44b=spherical_angle(xyz_center_in(:,ig+1,jg-1,lg),xyz_center_in(:,ig+1,jg  ,lg),xyz_center_out(:,i,j,lg))
                  endif
                  if (max(angle_44a,angle_44b)<=angle_44) then
                     ok=.true.
                     index_c2c(1,i,j,l)=ig
                     index_c2c(2,i,j,l)=jg-1
                     index_c2c(3,i,j,l)=lg
                  endif
               endif
            endif
         endif
      endif

    end subroutine set_closest_index
    !------------------------------------------------------------------!
  end subroutine get_c2c_weight
  !====================================================================!
  subroutine interpolate_c2c(ntiles, npx_in, npy_in, npx_out, npy_out,  &
                             dir_in, dir_out, index_c2c, weight_c2c)
    !------------------------------------------------------------------!
    ! interpolate data from restart files from one cubed sphere grid   !
    ! to another                                                       !
    !------------------------------------------------------------------!
#include "netcdf.inc"
    integer, intent(in) :: ntiles, npx_in, npy_in, npx_out, npy_out
    integer, intent(in) :: index_c2c (3, npx_out-1, npy_out-1, ntiles)
    real(FVPRC), intent(in) ::   weight_c2c (4, npx_out-1, npy_out-1, ntiles)

    character(len=120), intent(in) :: dir_in, dir_out
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer, dimension(:), allocatable :: dimlen_in ,dimlen_out, vartype
    integer, dimension(ntiles) :: ncid_in, ncid_out
    integer :: status, ndims, nvars, ngatts, time_dim
    integer :: ires, itile, istart, istop, ivar
    integer :: nz
    logical :: exists
    character(len=120), dimension(:), allocatable :: varname
    character(len=120) :: filename, basename
    character(len=120), save :: restart_file(3)=(/"fv_core.res   ",     &
                                                  "fv_srf_wnd.res",     &
                                                  "fv_tracer.res "/)
    !------------------------------------------------------------------!
    ! loop over restart files                                          !
    !------------------------------------------------------------------!
    do ires=1,3
       !---------------------------------------------------------------!
       ! loop over tiles                                               !
       !---------------------------------------------------------------!
       do itile=1,ntiles
          !------------------------------------------------------------!
          ! open input restart file                                    !
          !------------------------------------------------------------!
          write(filename,100) trim(dir_in), trim(restart_file(ires)), itile
100       format(a,"/",a,".tile",i1,".nc")
          inquire(file=filename,exist=exists)
          if (.not. exists) then
             print 101, trim(filename)
101          format("FATAL: restart file ",a," doesn't exist" )
             stop
          endif

          status = nf_open(trim(filename), 0, ncid_in(itile))
          if (status /= nf_noerr) then
             print*,"FATAL: nf_open: could not open file ",trim(filename)
             stop
          endif
          !------------------------------------------------------------!
          ! open output restart file                                   !
          !------------------------------------------------------------!
          write(filename,100) trim(dir_out), trim(restart_file(ires)), itile
          status = nf_create(filename, nf_clobber, ncid_out(itile))
          if (status /= nf_noerr) then
             print*,"FATAL: nf_create: could not create file ",trim(filename)
             stop
          endif
       enddo
       !---------------------------------------------------------------!
       ! set dimensions and variables for restart file                 !
       !---------------------------------------------------------------!
       call set_dims_and_vars(ires)
       !---------------------------------------------------------------!
       ! define dimensions and variables                               !
       ! copy global attributes                                        !
       !---------------------------------------------------------------!
       call define_dims_and_vars(ires)
       call copy_global_att()
       !---------------------------------------------------------------!
       ! write variables corresponding to dimensions                   !
       !---------------------------------------------------------------!
       call write_dimvar()
       !---------------------------------------------------------------!
       ! interpolate and write data                                    !
       !---------------------------------------------------------------!
       call write_var()
       deallocate(varname, vartype, dimlen_in, dimlen_out)
       !---------------------------------------------------------------!
       ! close netcdf files                                            !
       !---------------------------------------------------------------!
       do itile=1,ntiles
          status = nf_close(ncid_in (itile))
          status = nf_close(ncid_out(itile))
       enddo
    enddo

  contains
    !------------------------------------------------------------------!
    subroutine set_dims_and_vars(ifile)
      !----------------------------------------------------------------!
      ! set names of variables and dimension of output file            !
      !----------------------------------------------------------------!
      integer, intent(in) :: ifile

      integer :: z_id

      status = nf_inq(ncid_in, ndims, nvars, ngatts, time_dim)
      status = nf_inq_dimid(ncid_in, "zaxis_1", z_id)
      status = nf_inq_dimlen(ncid_in, z_id, nz)

      allocate(varname(nvars), vartype(nvars), dimlen_in(ndims), dimlen_out(ndims))

      select case (ifile)
      case (1)
         !-------------------------------------------------------------!
         !  fv_core.res                                                !
         !-------------------------------------------------------------!
         if (nz==1) then
            if (nvars/=11 .or. ndims/=6) then
               write(6,200) nvars, ndims
200            format ("FATAL: nvars, ndims are not equal ( 11, 6 ): ",2i3)
               stop
            endif
            
            varname(1:11)=(/"xaxis_1", "xaxis_2", "yaxis_1", "yaxis_2", &
                            "zaxis_1", "Time   ", "u      ", "v      ", &
                            "T      ", "delp   ", "phis   "/)
            dimlen_in (1:6)=(/npx_in-1 , npx_in , npy_in , npy_in-1 , nz, 1/)
            dimlen_out(1:6)=(/npx_out-1, npx_out, npy_out, npy_out-1, nz, 1/)
         else
            if (nvars/=12 .or. ndims/=7) then
               write(6,201) nvars, ndims
201            format ("FATAL: nvars, ndims are not equal ( 12, 7 ): ",2i3)
               stop
            endif
            
            varname(1:12)=(/"xaxis_1", "xaxis_2", "yaxis_1", "yaxis_2", &
                            "zaxis_1", "zaxis_2", "Time   ", "u      ", &
                            "v      ", "T      ", "delp   ", "phis   "/)
            dimlen_in (1:7)=(/npx_in-1 , npx_in , npy_in , npy_in-1 , nz, 1, 1/)
            dimlen_out(1:7)=(/npx_out-1, npx_out, npy_out, npy_out-1, nz, 1, 1/)
         endif

      case (2)
         !-------------------------------------------------------------!
         !  fv_srf_wnd.res                                             !
         !-------------------------------------------------------------!
         if (nvars/=6 .or. ndims/=4) then
            write(6,200) nvars, ndims
210         format ("FATAL: nvars, ndims are not equal ( 6, 4 ): ",2i3)
            stop
         endif

         varname(1:6)=(/"xaxis_1", "yaxis_1", "zaxis_1", "Time   ",     &
                        "u_srf  ", "v_srf  "/)
         dimlen_in (1:4)=(/npx_in-1,  npy_in-1,  1, 1/)
         dimlen_out(1:4)=(/npx_out-1, npy_out-1, 1, 1/)

      case (3)
         !-------------------------------------------------------------!
         !  fv_tracer.res                                              !
         !-------------------------------------------------------------!
         if (ndims/=4) then
            write(6,200) ndims
220         format ("FATAL: ndims are not equal 4: ",i3)
            stop
         endif

         do ivar=1,nvars
            write(varname(ivar),230) ivar
         enddo
230      format ("atm_T",i1)   

         dimlen_in (1:4)=(/npx_in-1,  npy_in-1,  nz, 1/)
         dimlen_out(1:4)=(/npx_out-1, npy_out-1, nz, 1/)
         
      case default
         write(6,240)
240      format("FATAL: unknown restart file")
         stop
      end select

    end subroutine set_dims_and_vars
    !------------------------------------------------------------------!
    subroutine define_dims_and_vars(ifile)
      !----------------------------------------------------------------!
      ! define variables and dimensions in output file                 !
      !----------------------------------------------------------------!
      integer, intent(in) :: ifile

      integer :: idim, iatt, varndims, vardimids(4), varnatts, dimid, varid
      character(len=120) :: name, att_name
      !----------------------------------------------------------------!
      ! define dimensions                                              !
      !----------------------------------------------------------------!
      do idim=1,ndims-1
         status = nf_inq_dimname(ncid_in(1), idim, name)
         do itile=1,ntiles
            status = nf_def_dim(ncid_out(itile), trim(name), dimlen_out(idim), dimid)
         enddo
      enddo
      if (time_dim/=ndims) then
         write(6,300)
300      format("FATAL: last dimension is not unlimited dimension")
         stop
      endif
      do itile=1,ntiles
         status = nf_def_dim(ncid_out(itile), "Time", nf_unlimited, dimid) 
      enddo
      !----------------------------------------------------------------!
      ! define variables                                               !
      !----------------------------------------------------------------!
      do ivar=1,nvars
         if (ifile<3) then
            status = nf_inq_varname(ncid_in(1), ivar, name)
            if (trim(name)/=trim(varname(ivar))) then
               write(6,310) trim(varname(ivar)), trim(name)
310            format("FATAL: varname different from expected name( ",a," ): ",a)
               stop
            endif
         endif
         status = nf_inq_var(ncid_in(1), ivar, varname(ivar), vartype(ivar),    &
                             varndims, vardimids, varnatts)
         do itile=1,ntiles
            status = nf_def_var(ncid_out(itile), trim(varname(ivar)), vartype(ivar),   &
                                varndims, vardimids, varid)
         enddo
         do iatt=1,varnatts
             status = nf_inq_attname(ncid_in(1), ivar, iatt, att_name)
             do itile=1,ntiles
                status = nf_copy_att(ncid_in(1), ivar, att_name, ncid_out(itile), varid)
             enddo
          enddo
      enddo

    end subroutine define_dims_and_vars
    !------------------------------------------------------------------!
    subroutine copy_global_att()
      !----------------------------------------------------------------!
      ! copy global attributes from input to output file               !
      !----------------------------------------------------------------!
      integer :: iatt
      character(len=120) :: att_name

      do iatt=1,ngatts
         status = nf_inq_attname(ncid_in(1), nf_global, iatt, att_name)
         do itile=1,ntiles
            status = nf_copy_att(ncid_in(1), nf_global, att_name, ncid_out(itile), nf_global)
         enddo
      enddo
      do itile=1,ntiles
         status = nf_put_att_text(ncid_out(itile), nf_global, "description", 61,     &
              "data interpolated from cubed sphere with different resolution")
         status = nf_enddef(ncid_out(itile))
      enddo

    end subroutine copy_global_att
    !------------------------------------------------------------------!
    subroutine write_dimvar()
      !----------------------------------------------------------------!
      ! write variables that describe dimensions                       !
      !----------------------------------------------------------------!
      real*4, dimension(:), allocatable :: var_r4
      integer :: ivar, i, dimlen

      do ivar=1,ndims-1
         status = nf_inq_dimlen(ncid_in(1), ivar, dimlen)
         if (dimlen/=dimlen_in(ivar)) then
            write(6, 400) dimlen, dimlen_in(ivar)
400         format("FATAL: actual dimension length does not equal to expected length: ", 2i5)
            stop
         endif
         if (vartype(ivar)/=nf_float) then
            write(6,410) vartype(ivar)
410         format("FATAL: vartype not equal nf_float: ", i3)
            stop
         endif
         allocate(var_r4(dimlen_out(ivar)))
         do i=1,dimlen_out(ivar)
            var_r4(i)=i
         enddo
         do itile=1,ntiles
            status = nf_put_vara_real(ncid_out(itile), ivar, 1, dimlen_out(ivar), var_r4)
         enddo
         deallocate(var_r4)
      enddo
      do itile=1,ntiles
         status = nf_put_var1_double(ncid_out(itile), ndims, 1, 1.d0)
      enddo

    end subroutine write_dimvar
    !------------------------------------------------------------------!
    subroutine write_var()
      !----------------------------------------------------------------!
      ! interpolate and write data                                     !
      !----------------------------------------------------------------!
      use GHOST_CUBSPH_mod, only: ghost_cubsph_update, A_grid
      use GRID_UTILS_mod,   only: get_dx, get_dxa, get_dy, get_dya,     &
                                  get_center_vect, get_west_vect,       &
                                  get_south_vect, get_cosa_center
      use FLOW_PROJ_mod,    only: d2a_vect, a2d_vect
      
      real*8, dimension(:,:,:), allocatable :: var_r8
      real(FVPRC), dimension(:,:,:,:,:), allocatable :: vxyz_in, vxyz_out
      real(FVPRC), dimension(:,:,:,:), allocatable :: var_in, var_out
      real(FVPRC), dimension(:,:,:), allocatable :: u, v ,ec1, ec2, ew1, ew2, es1, es2
      real(FVPRC), dimension(:,:), allocatable :: dx, dy, dxa, dya, rdxa, rdya, cosa_s, sina_s
      real(FVPRC), dimension(:), allocatable :: edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n

      integer, dimension(4) :: start, count
      integer :: itile, k

      start(1:4)=(/1,1,1,1/)
      do ivar=ndims+1,nvars
         if (vartype(ivar)/=nf_double) then
            write(6,510) vartype(ivar)
510         format("FATAL: vartype not equal nf_double: ", i3)
            stop
         endif
         
         if (trim(varname(ivar))=="u") then
            !----------------------------------------------------------!
            ! horizontal flow                                          !
            !----------------------------------------------------------!
            allocate(u(0:npx_in,0:npy_in+1,nz), v(0:npx_in+1,0:npy_in,nz), &
                     vxyz_in(3,0:npx_in,0:npy_in,nz,ntiles))
            allocate(dx(0:npx_in,0:npy_in+1), dxa(0:npx_in,0:npy_in), rdxa(0:npx_in,0:npy_in))
            allocate(dy(0:npx_in+1,0:npy_in), dya(0:npx_in,0:npy_in), rdya(0:npx_in,0:npy_in))
            allocate(ec1(3,0:npx_in,0:npy_in), ec2(3,0:npx_in,0:npy_in))
            allocate(cosa_s(0:npx_in,0:npy_in), sina_s(0:npx_in,0:npy_in))
            !----------------------------------------------------------!
            ! loop over tiles                                          !
            !----------------------------------------------------------!
            do itile=1,ntiles
               !-------------------------------------------------------!
               ! read horizontal flow                                  !
               !-------------------------------------------------------!
               allocate(var_r8(npx_in-1,npy_in,nz))
               count(1:4)=(/npx_in-1, npy_in, nz, 1/)
               status = nf_get_vara_double(ncid_in(itile), ivar, start, count, var_r8)
               u(1:npx_in-1,1:npy_in,:)=var_r8(:,:,:)
               deallocate(var_r8)
               allocate(var_r8(npx_in,npy_in-1,nz))
               count(1:4)=(/npx_in, npy_in-1, nz, 1/)
               status = nf_get_vara_double(ncid_in(itile), ivar+1, start, count, var_r8)
               v(1:npx_in,1:npy_in-1,:)=var_r8(:,:,:)
               deallocate(var_r8)
               !-------------------------------------------------------!
               ! geometrical properties of input grid                  !
               !-------------------------------------------------------!
               call get_dx (xyz_corner_in(:,:,:,itile), 0, npx_in, 0, npy_in,             &
                                                        0, npx_in, 0, npy_in, dx)
               call get_dxa(xyz_corner_in(:,:,:,itile), 0, npx_in, 0, npy_in,             &
                                                        0, npx_in, 0, npy_in, dxa, rdxa=rdxa)
               call get_dy (xyz_corner_in(:,:,:,itile), 0, npx_in, 0, npy_in,             &
                                                        0, npx_in, 0, npy_in, dy)
               call get_dya(xyz_corner_in(:,:,:,itile), 0, npx_in, 0, npy_in,             &
                                                        0, npx_in, 0, npy_in, dya, rdya=rdya)
               call get_center_vect(xyz_corner_in(:,:,:,itile), 0, npx_in, 0, npy_in,     &
                                                                0, npx_in, 0, npy_in, ec1, ec2)
               call get_cosa_center(ec1, ec2, 0, npx_in, 0, npy_in,     &
                                              0, npx_in, 0, npy_in, cosa_s, sina_s)
               !-------------------------------------------------------!
               ! calculate flow vector for a-grid                      !
               !-------------------------------------------------------!
               call d2a_vect(u, v, dx, dy, rdxa, rdya, cosa_s, ec1, ec2, &
                             0, npx_in  , 0, npy_in  , 1, nz,            &
                             1, npx_in-1, 1, npy_in-1, 1, nz,            &
                             vxyz_in(:,:,:,:,itile))
            enddo
            deallocate(u, v, dx, dy, dxa, dya, rdxa, rdya, ec1, ec2, cosa_s, sina_s)
            allocate(vxyz_out(3,0:npx_out,0:npy_out,nz,ntiles))
            !----------------------------------------------------------!
            ! ghost cell update of vxyz_in                               !
            !----------------------------------------------------------!
            vxyz_in(:,0,     0,     :,:)=0.
            vxyz_in(:,npx_in,0     ,:,:)=0.
            vxyz_in(:,npx_in,npy_in,:,:)=0.
            vxyz_in(:,0,     npy_in,:,:)=0.
            do k=1,3
               do itile=1,ntiles
                  call ghost_cubsph_update(vxyz_in(k,:,:,:,:), 0, npx_in, 0, npy_in, nz, 1, ntiles, &
                                           1, nz, itile, A_grid)
               enddo
            enddo
            !----------------------------------------------------------!
            ! do interpolation of flow vector                          !
            !----------------------------------------------------------!
            do k=1,3
               call do_c2c_interpolation(vxyz_in(k,:,:,:,:), 0, npx_in, 0, npy_in, nz, ntiles, &
                                         index_c2c, weight_c2c, npx_out-1, npy_out-1,          &
                                         vxyz_out(k,1:npx_out-1,1:npy_out-1,:,:))
            enddo
            deallocate(vxyz_in)
            !----------------------------------------------------------!
            ! loop over tiles                                          !
            !----------------------------------------------------------!
            allocate(u(0:npx_out,0:npy_out+1,nz), v(0:npx_out+1,0:npy_out,nz))
            allocate(ew1(3,0:npx_out,0:npy_out+1), ew2(3,0:npx_out,0:npy_out+1))
            allocate(es1(3,0:npx_out+1,0:npy_out), es2(3,0:npx_out+1,0:npy_out))
            allocate(edge_vect_w(0:npy_out), edge_vect_e(0:npy_out),    &
                     edge_vect_s(0:npx_out), edge_vect_n(0:npx_out))
            do itile=1,ntiles
               !-------------------------------------------------------!
               ! ghost cell update of vxyz_out                           !
               !-------------------------------------------------------!
               do k=1,3
                  call ghost_cubsph_update(vxyz_out(k,:,:,:,:), 0, npx_out, 0, npy_out, nz, 1, ntiles, &
                                           1, nz, itile, A_grid)
               enddo
               !-------------------------------------------------------!
               ! geometrical properties of output grid                 !
               !-------------------------------------------------------!
               call get_west_vect(xyz_corner_out(:,:,:,itile),                       &
                    0, npx_out, 0, npy_out, 1, npx_out-1, 1, npy_out-1,              &
                    west_edge, east_edge, 1, npx_out-1, ew1, ew2)
               call get_south_vect(xyz_corner_out(:,:,:,itile),                      &
                    0, npx_out, 0, npy_out, 1, npx_out-1, 1, npy_out-1,              &
                    west_edge, east_edge, 1, npy_out-1, es1, es2)
               !-------------------------------------------------------!
               ! calculate co-variant flow components on d_grid        !
               !-------------------------------------------------------!
               call a2d_vect(vxyz_out(:,:,:,:,itile), ew1, ew2, es1, es2,            &
                    0, npx_out,   0, npy_out,   1, nz,                               &
                    1, npx_out-1, 1, npy_out-1, 1, nz,                               &
                    edge_interp, edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n, &
                    west_edge, east_edge, south_edge, north_edge,                    &
                    sw_corner, se_corner, nw_corner, ne_corner,                      &
                    1, npx_out-1, 1, npy_out-1, u, v)
               !-------------------------------------------------------!
               ! write flow components                                 !
               !-------------------------------------------------------!
               allocate(var_r8(npx_out-1,npy_out,nz))
               var_r8(1:npx_out-1, 1:npy_out, 1:nz)=u(1:npx_out-1, 1:npy_out, 1:nz)
               count(1:4)=(/npx_out-1, npy_out, nz, 1/)
               status = nf_put_vara_double(ncid_out(itile), ivar, start, count, var_r8)
               deallocate(var_r8)

               allocate(var_r8(npx_out,npy_out-1,nz))
               var_r8(1:npx_out, 1:npy_out-1, 1:nz)=v(1:npx_out, 1:npy_out-1, 1:nz)
               count(1:4)=(/npx_out, npy_out-1, nz, 1/)
               status = nf_put_vara_double(ncid_out(itile), ivar+1, start, count, var_r8)
               deallocate(var_r8)
            enddo
            deallocate(vxyz_out, u, v, ew1, ew2, es1, es2,                 &
                       edge_vect_w, edge_vect_e, edge_vect_s, edge_vect_n)
         elseif (trim(varname(ivar))=="phis") then
            !----------------------------------------------------------!
            ! read topography                                          !
            !----------------------------------------------------------!
            allocate(var_r8(npx_in-1,npy_in-1,1), var_in(0:npx_in,0:npy_in,1,ntiles))
            count(1:4)=(/npx_in-1, npy_in-1, 1, 1/)
            do itile=1,ntiles
               status = nf_get_vara_double(ncid_in(itile), ivar, start, count, var_r8)
               var_in(1:npx_in-1,1:npy_in-1,1,itile)=var_r8(:,:,1)
            enddo
            deallocate(var_r8)
            !----------------------------------------------------------!
            ! interpolate topography                                   !
            !----------------------------------------------------------!
            var_in(0,     0,     :,:)=0.
            var_in(npx_in,0     ,:,:)=0.
            var_in(npx_in,npy_in,:,:)=0.
            var_in(0,     npy_in,:,:)=0.
            do itile=1,ntiles
               call ghost_cubsph_update(var_in(:,:,:,:), 0, npx_in, 0, npy_in, 1, 1, ntiles, &
                                        1, 1, itile, A_grid)
            enddo
            allocate(var_out(0:npx_out,0:npy_out,1,ntiles))
            call do_c2c_interpolation(var_in(:,:,:,:), 0, npx_in, 0, npy_in, 1, ntiles, &
                                      index_c2c, weight_c2c, npx_out-1, npy_out-1,       &
                                      var_out(1:npx_out-1,1:npy_out-1,:,:))
            deallocate(var_in)
            !----------------------------------------------------------!
            ! write topography                                         !
            !----------------------------------------------------------!
            allocate(var_r8(npx_out-1,npy_out-1,1))
            count(1:4)=(/npx_out-1, npy_out-1, 1, 1/)
            do itile=1,ntiles
               var_r8(1:npx_out-1,1:npy_out-1,1)=var_out(1:npx_out-1,1:npy_out-1,1,itile)
               status = nf_put_vara_double(ncid_out(itile), ivar, start, count, var_r8)
            enddo
            deallocate(var_r8, var_out) 
         elseif (trim(varname(ivar))/="v") then
            !----------------------------------------------------------!
            ! read variable                                            !
            !----------------------------------------------------------!
            allocate(var_r8(npx_in-1,npy_in-1,nz), var_in(0:npx_in,0:npy_in,nz,ntiles))
            count(1:4)=(/npx_in-1, npy_in-1, nz, 1/)
            do itile=1,ntiles
               status = nf_get_vara_double(ncid_in(itile), ivar, start, count, var_r8)
               var_in(1:npx_in-1,1:npy_in-1,:,itile)=var_r8(:,:,:)
            enddo
            deallocate(var_r8)
            !----------------------------------------------------------!
            ! interpolate variable                                     !
            !----------------------------------------------------------!
            var_in(0,     0,     :,:)=0.
            var_in(npx_in,0     ,:,:)=0.
            var_in(npx_in,npy_in,:,:)=0.
            var_in(0,     npy_in,:,:)=0.
            do itile=1,ntiles
               call ghost_cubsph_update(var_in(:,:,:,:), 0, npx_in, 0, npy_in, nz, 1, ntiles, &
                                        1, nz, itile, A_grid)
            enddo
            allocate(var_out(0:npx_out,0:npy_out,nz,ntiles))
            call do_c2c_interpolation(var_in(:,:,:,:), 0, npx_in, 0, npy_in, nz, ntiles, &
                                      index_c2c, weight_c2c, npx_out-1, npy_out-1,       &
                                      var_out(1:npx_out-1,1:npy_out-1,:,:))
            deallocate(var_in)
            !----------------------------------------------------------!
            ! write variable                                           !
            !----------------------------------------------------------!
            allocate(var_r8(npx_out-1,npy_out-1,nz))
            count(1:4)=(/npx_out-1, npy_out-1, nz, 1/)
            do itile=1,ntiles
               var_r8(1:npx_out-1,1:npy_out-1,1:nz)=var_out(1:npx_out-1,1:npy_out-1,1:nz,itile)
               status = nf_put_vara_double(ncid_out(itile), ivar, start, count, var_r8)
            enddo
            deallocate(var_r8, var_out) 
         endif
      enddo

    end subroutine write_var
    !------------------------------------------------------------------!
  end subroutine interpolate_c2c
  !====================================================================!
  subroutine do_c2c_interpolation(var_in, isd, ied, jsd, jed, nz, ntiles, &
                                  index_c2c, weight_c2c, nx, ny, var_out)
    !------------------------------------------------------------------!
    ! do interpolation from var_in to var_out                          !
    !------------------------------------------------------------------!
    integer, intent(in) :: isd, ied, jsd, jed, nz, ntiles, nx, ny
    real(FVPRC), dimension(isd:ied,jsd:jed,nz,ntiles), intent(in)  :: var_in
    real(FVPRC), dimension(4,nx,ny,ntiles),            intent(in)  :: weight_c2c
    integer, dimension(3,nx,ny,ntiles),         intent(in)  :: index_c2c
    real(FVPRC), dimension(nx,ny,nz,ntiles),           intent(out) :: var_out
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    integer :: i, j, k, l, ic, jc, lc

    do l=1,ntiles
       do k=1,nz
          do j=1,ny
             do i=1,nx
                ic=index_c2c(1,i,j,l)
                jc=index_c2c(2,i,j,l)
                lc=index_c2c(3,i,j,l)

                var_out(i,j,k,l)=weight_c2c(1,i,j,l)*var_in(ic  ,jc  ,k, lc)  &
                                +weight_c2c(2,i,j,l)*var_in(ic  ,jc+1,k, lc)  &
                                +weight_c2c(3,i,j,l)*var_in(ic+1,jc+1,k, lc)  &
                                +weight_c2c(4,i,j,l)*var_in(ic+1,jc  ,k, lc)
             enddo
          enddo
       enddo
    enddo

  end subroutine do_c2c_interpolation
  !====================================================================!
end module CUB2CUB_mod
