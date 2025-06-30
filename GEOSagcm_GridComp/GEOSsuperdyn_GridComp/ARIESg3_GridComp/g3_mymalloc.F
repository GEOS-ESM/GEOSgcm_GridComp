      subroutine malloc_1d_r (a,im)
      implicit none
      real(kind=8), dimension(:), pointer :: a
      integer i,im,m
      if(.not.associated(a)) then
      allocate(a(im))
      do i=1,im
      a(i) = 0.0
      enddo
      else
             m=size(a)
          if(m.ne.im) then
          print *, 'Allocated Array Size (',m,') does not match request (',im,')!'  
          call my_finalize
          call my_exit (101)
          endif
      endif
      return
      end


      subroutine malloc_2d_r (a,im,jm)
      implicit none
      real(kind=8), dimension(:,:), pointer :: a
      integer i,j,m,im,jm
      if(.not.associated(a)) then
      allocate(a(im,jm))
      do j=1,jm
      do i=1,im
      a(i,j) = 0.0
      enddo
      enddo
      else
             m=size(a)
          if(m.ne.im*jm) then
          print *, 'Allocated Array Size (',m,') does not match request (',im*jm,')!'
          call my_finalize
          call my_exit (101)
          endif
      endif
      return
      end


      subroutine malloc_3d_r (a,im,jm,lm)
      implicit none
      real(kind=8), dimension(:,:,:), pointer :: a
      integer i,j,l,im,jm,lm,m
      if(.not.associated(a)) then
      allocate(a(im,jm,lm))
      do l=1,lm
      do j=1,jm
      do i=1,im
      a(i,j,l) = 0.0
      enddo
      enddo
      enddo
      else
             m=size(a)
          if(m.ne.im*jm*lm) then
          print *, 'Allocated Array Size (',m,') does not match request (',im*jm*lm,')!'
          call my_finalize
          call my_exit (101)
          endif
      endif
      return
      end


      subroutine malloc_4d_r (a,im,jm,lm,nm)
      implicit none
      real(kind=8), dimension(:,:,:,:), pointer :: a
      integer i,j,l,n,im,jm,lm,nm,m
      if(.not.associated(a)) then
      allocate(a(im,jm,lm,nm))
      do n=1,nm
      do l=1,lm
      do j=1,jm
      do i=1,im
      a(i,j,l,n) = 0.0
      enddo
      enddo
      enddo
      enddo
      else
             m=size(a)
          if(m.ne.im*jm*lm*nm) then
          print *, 'Allocated Array Size (',m,') does not match request (',im*jm*lm*nm,')!'
          call my_finalize
          call my_exit (101)
          endif
      endif
      return
      end


      subroutine malloc_1d_i (a,im)
      implicit none
      integer, dimension(:), pointer :: a
      integer i,im,m
      if(.not.associated(a)) then
      allocate(a(im))
      do i=1,im
      a(i) = 0
      enddo
      else
             m=size(a)
          if(m.ne.im) then
          print *, 'Allocated Array Size (',m,') does not match request (',im,')!'
          call my_finalize
          call my_exit (101)
          endif
      endif
      return
      end


      subroutine malloc_2d_i (a,im,jm)
      implicit none
      integer, dimension(:,:), pointer :: a
      integer i,j,m,im,jm
      if(.not.associated(a)) then
      allocate(a(im,jm))
      do j=1,jm
      do i=1,im
      a(i,j) = 0
      enddo
      enddo
      else
             m=size(a)
          if(m.ne.im*jm) then
          print *, 'Allocated Array Size (',m,') does not match request (',im*jm,')!'
          call my_finalize
          call my_exit (101)
          endif
      endif
      return
      end


      subroutine malloc_3d_i (a,im,jm,lm)
      implicit none
      integer, dimension(:,:,:), pointer :: a
      integer i,j,l,im,jm,lm,m
      if(.not.associated(a)) then
      allocate(a(im,jm,lm))
      do l=1,lm
      do j=1,jm
      do i=1,im
      a(i,j,l) = 0
      enddo
      enddo
      enddo
      else
             m=size(a)
          if(m.ne.im*jm*lm) then
          print *, 'Allocated Array Size (',m,') does not match request (',im*jm*lm,')!'
          call my_finalize
          call my_exit (101)
          endif
      endif
      return
      end


      subroutine malloc_4d_i (a,im,jm,lm,nm)
      implicit none
      integer, dimension(:,:,:,:), pointer :: a
      integer i,j,l,n,im,jm,lm,nm,m
      if(.not.associated(a)) then
      allocate(a(im,jm,lm,nm))
      do n=1,nm
      do l=1,lm
      do j=1,jm
      do i=1,im
      a(i,j,l,n) = 0
      enddo
      enddo
      enddo
      enddo
      else
             m=size(a)
          if(m.ne.im*jm*lm*nm) then
          print *, 'Allocated Array Size (',m,') does not match request (',im*jm*lm*nm,')!'
          call my_finalize
          call my_exit (101)
          endif
      endif
      return
      end


      subroutine malloc_1d_c8 (a,i)
      implicit none
      character*8, dimension(:), pointer :: a
      integer i,m
      if(.not.associated(a)) then
      allocate(a(i))
      else
             m=size(a)
          if(m.ne.i) then
          print *, 'Allocated Array Size (',m,') does not match request (',i,')!'  
          call my_finalize
          call my_exit (101)
          endif
      endif
      return
      end


      subroutine malloc_2d_c8 (a,i,j)
      implicit none
      character*8, dimension(:,:), pointer :: a
      integer i,j,m
      if(.not.associated(a)) then
      allocate(a(i,j))
      else
             m=size(a)
          if(m.ne.i*j) then
          print *, 'Allocated Array Size (',m,') does not match request (',i*j,')!'
          call my_finalize
          call my_exit (101)
          endif
      endif
      return
      end


      subroutine malloc_3d_c8 (a,i,j,k)
      implicit none
      character*8, dimension(:,:,:), pointer :: a
      integer i,j,k,m
      if(.not.associated(a)) then
      allocate(a(i,j,k))
      else
             m=size(a)
          if(m.ne.i*j*k) then
          print *, 'Allocated Array Size (',m,') does not match request (',i*j*k,')!'
          call my_finalize
          call my_exit (101)
          endif
      endif
      return
      end


      subroutine malloc_4d_c8 (a,i,j,k,l)
      implicit none
      character*8, dimension(:,:,:,:), pointer :: a
      integer i,j,k,l,m
      if(.not.associated(a)) then
      allocate(a(i,j,k,l))
      else
             m=size(a)
          if(m.ne.i*j*k*l) then
          print *, 'Allocated Array Size (',m,') does not match request (',i*j*k*l,')!'
          call my_finalize
          call my_exit (101)
          endif
      endif
      return
      end
