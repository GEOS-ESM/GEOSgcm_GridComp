      subroutine G3_AVRX ( U,IM,JM,LM,S,lattice )
      use G3_MPI_Util_Mod
      implicit none
      type ( dynamics_lattice_type ) lattice
      INTEGER         IM,JM,LM
      REAL(kind=8)  U(IM,JM,LM)
      REAL(kind=8)  S(lattice%imglobal+2,JM)

      integer  status(mpi_status_size)
      integer   stats(mpi_status_size,0:lattice%nx-1)
      integer   statr(mpi_status_size,0:lattice%nx-1)

      integer  sendquest(0:lattice%nx-1)
      integer  recvquest(0:lattice%nx-1)
      integer  ierror,peid

      logical first,  flag
      data    first /.true./

      INTEGER    FORWARD,    BACKWARD
      PARAMETER (FORWARD=-1, BACKWARD=1)

      REAL(kind=8) ONE
      REAL(kind=8) ZERO
      REAL(kind=8) SMIN
      PARAMETER ( ONE=1.0)
      PARAMETER (ZERO=0.0)

      real(kind=8),    allocatable :: sendbuf(:,:)
      real(kind=8),    allocatable :: recvbuf(:,:)
      real(kind=8),    allocatable ::       z(:,:)
      real(kind=8),    allocatable ::       b(:,:)
      integer, allocatable ::      j2(:)

      INTEGER I,J,L, J1(JM*LM), L1(JM*LM)
      INTEGER NFFTS
      INTEGER IM0
      DATA    IM0/0/

      integer n,num,rem,isum,lsum,len,len0
      integer,                   save :: ix(19)
      real(kind=8), allocatable, save :: tr(:)

      if (lattice%imglobal.ne.im0) then
        if(im0.ne.0) deallocate ( tr )
                       allocate ( tr(lattice%imglobal*2) )
        call fftfax (lattice%imglobal,ix,tr)
        im0=lattice%imglobal
      endif

! Compute Number of FFTs to Perform and Load Data into Buffer
! -----------------------------------------------------------
      allocate ( sendbuf(im,jm*lm) )

      nffts = 0
      do j=1,jm
            smin = minval( s(:,j) )
        if( smin.lt.0.9999 ) then
          do L=1,lm
             nffts  = nffts + 1
          j1(nffts) = j
          L1(nffts) = L
             do i=1,im
             sendbuf(i,nffts) = u(i,j,L)
             enddo
          enddo
        endif
      enddo

      num = nffts/lattice%nx
      rem = nffts-lattice%nx*num

                                 len0 = num    ! Define number of FFTs on myid
      if( lattice%pei.le.rem-1 ) len0 = num+1  ! Define number of FFTs on myid

#if 0
      if( first ) then
          if( lattice%myid.eq.0 ) then
          print *, 'FFT Load Balance for Upper-Air Field:'
          print *, '-------------------------------------'
          endif
            do n=0,lattice%nx*lattice%ny-1
            if( n.eq.lattice%myid ) then
            write(6,1000) lattice%myid,lattice%pei,lattice%pej,nffts,num,rem,len0
            if( mod(n+1,lattice%nx).eq.0 ) print *
            endif
            call my_barrier (lattice%comm)
            enddo
            if( lattice%myid.eq.lattice%nx*lattice%ny-1 ) print *
      first = .false.
      endif
 1000 format(1x,'absolute PE id: ',i3,'  relative (pei,pej): ',i2,',',i2,'  nffts: ',i6,2x, &
                'num: ',i6,2x,'rem: ',i6,2x,'len0: ',i6)
#endif

! Distribute Data Across PEs in X-direction
! -----------------------------------------
      if( nffts.ne.0 ) then

      allocate ( j2(len0) )
      allocate ( z(lattice%imglobal+2,len0  ) )
      allocate ( b(lattice%imglobal+2,len0*2) )

      isum = 0
      lsum = 0
      do n = 0,lattice%nx-1
      peid = n + lattice%pej*lattice%nx
                       len = num
      if( n.le.rem-1 ) len = num+1
              sendquest(n) = mpi_request_null
      if( len.ne.0   ) then
          if( peid.ne.lattice%myid ) then
              call mpi_isend ( sendbuf(1,1+lsum),im*len,mpi_double_precision,peid,peid,lattice%comm,sendquest(n),ierror )
          else
              do L=1,len0
                 j2(L) = j1(L+lsum)
                 do i=1,im
                 z(i+isum,L) = sendbuf(i,L+lsum)
                 enddo
              enddo
          endif
      endif
      isum = isum + lattice%im(n)
      lsum = lsum + len
      enddo

! Receive Data and Perform FFT
! ----------------------------
      if( len0.ne.0 ) then

      allocate ( recvbuf( lattice%imglobal*len0,0:lattice%nx-1 ) )
      do n = 0,lattice%nx-1
      peid = n + lattice%pej*lattice%nx
      if( peid.ne.lattice%myid ) then
          call mpi_irecv ( recvbuf(1,n),lattice%im(n)*len0,mpi_double_precision,peid,lattice%myid,lattice%comm,recvquest(n),ierror )
      else
          recvquest(n) = mpi_request_null
      endif
      enddo

      call mpi_waitall ( lattice%nx,sendquest(0:lattice%nx-1),stats(1,0),ierror )
      call mpi_waitall ( lattice%nx,recvquest(0:lattice%nx-1),statr(1,0),ierror )

      isum = 0
      do n = 0,lattice%nx-1
      peid = n + lattice%pej*lattice%nx
      if( peid.ne.lattice%myid ) then
           lsum = 0
           do L=1,len0
           do i=1,lattice%im(n)
           lsum = lsum+1
            z(i+isum,L) = recvbuf(lsum,n)
           enddo
           enddo
      endif
      isum = isum + lattice%im(n)
      enddo

      do L=1,len0
      z(lattice%imglobal+1,L) = 0.0
      z(lattice%imglobal+2,L) = 0.0
      enddo

! Perform FFT
! -----------
      call rfftmlt( Z,B,TR,IX,1,lattice%imglobal+2,lattice%imglobal,len0,FORWARD )
      do L=1,len0
      do i=1,lattice%imglobal+2
       z(i,L) = z(i,L)*s(i,j2(L))
      enddo
      enddo
      call rfftmlt( Z,B,TR,IX,1,lattice%imglobal+2,lattice%imglobal,len0,BACKWARD )

! Distribute Filtered Data Back to source PEs
! -------------------------------------------
      deallocate ( sendbuf )
      deallocate ( recvbuf )
        allocate ( sendbuf( lattice%imglobal*len0,0:lattice%nx-1 ) )
      isum = 0
      do n = 0,lattice%nx-1
      peid = n + lattice%pej*lattice%nx
      if( peid.ne.lattice%myid ) then

          lsum = 0
          do L=1,len0
          do i=1,lattice%im(n)
          lsum = lsum+1
          sendbuf(lsum,n) = z(i+isum,L)
          enddo
          enddo
          call mpi_isend ( sendbuf(1,n),lsum,mpi_double_precision,peid,peid,lattice%comm,sendquest(n),ierror )

      else
          sendquest(n) = mpi_request_null
      endif
      isum = isum + lattice%im(n)
      enddo

      endif  ! End len0.ne.0 check

! Receive Filtered Data
! ---------------------
      allocate ( recvbuf(im,nffts) )
      isum = 0
      lsum = 0
      do n = 0,lattice%nx-1
      peid = n + lattice%pej*lattice%nx
                       len = num
      if( n.le.rem-1 ) len = num+1
              recvquest(n) = mpi_request_null
      if( len.ne.0   ) then
          if( peid.ne.lattice%myid ) then
              call mpi_irecv ( recvbuf(1,1+lsum),im*len,mpi_double_precision,peid,lattice%myid,lattice%comm,recvquest(n),ierror )
          else
              do L=1,len0
              do i=1,im
                recvbuf(i,L+lsum) = z(i+isum,L)
              enddo
              enddo
          endif
      endif
      isum = isum + lattice%im(n)
      lsum = lsum + len
      enddo

      call mpi_waitall ( lattice%nx,sendquest(0:lattice%nx-1),stats(1,0),ierror )
      call mpi_waitall ( lattice%nx,recvquest(0:lattice%nx-1),statr(1,0),ierror )

! Reconstruct Filtered Field
! --------------------------
      do n=1,nffts
      do i=1,im
      u(i,j1(n),L1(n)) = recvbuf(i,n)
      enddo
      enddo

      deallocate ( z,b,j2  )
      deallocate ( sendbuf )
      deallocate ( recvbuf )
      else
      deallocate ( sendbuf )
      endif

      return
      end
