      module g3_dynamics_lattice_module

c Define Lattice
c --------------
      type dynamics_lattice_type
           integer            :: imglobal      ! Global  Size in X
           integer            :: jmglobal      ! Global  Size in Y
           integer            :: lm            ! Global  Size in Z (Note: Same as Local Size)
           integer            :: imax          ! Maximum local im
           integer            :: jmax          ! Maximum local jm
           integer            :: comm          ! MPI Communicator
           integer            :: nx,ny         ! Size of PE lattice
           integer            :: pei,pej       ! Relative address of local PE on lattice
           integer            :: myid          ! Absolute address of local PE
           integer            :: mpi_rkind     ! mpi_real or mpi_double_precision based on kind
           integer, pointer   :: npeg(:)       ! Number of pole PE ghosts per processor
           integer, pointer   ::  img(:,:)     ! Number of grid-points associated with each pole PE ghost
           integer, pointer   ::  im0(:,:)     ! Beginning address     associated with each pole PE ghost
           integer, pointer   :: ppeg(:)       ! Relative address of pole PE ghost for each iglobal
           integer, pointer   ::   im(:)       ! Array  of local zonal      dimension for each PE in x
           integer, pointer   ::   jm(:)       ! Array  of local meridional dimension for each PE in y
           integer, pointer   ::    ilocal(:)  ! Array  of  local i-index for global i-location
           integer, pointer   ::    jlocal(:)  ! Array  of  local j-index for global j-location
           integer, pointer   ::   iglobal(:)  ! Array  of global i-index for  local i-location
           integer, pointer   ::   jglobal(:)  ! Array  of global j-index for  local j-location
           integer, pointer   :: peiglobal(:)  ! Relative PE address associated with iglobal
           integer, pointer   :: pejglobal(:)  ! Relative PE address associated with jglobal
      endtype dynamics_lattice_type

      contains

C **********************************************************************
      subroutine create_dynamics_lattice (lattice,nx,ny)
C **********************************************************************
      implicit none
      type ( dynamics_lattice_type) lattice
      integer n,m,nx,ny

c Lattice%im
c ----------
      if(.not.associated(lattice%im)) then
      allocate ( lattice%im(0:nx-1) )
      do n=0,nx-1
      lattice%im(n) = 0
      enddo
      else
             m=size(lattice%im)
          if(m.ne.nx) then
          print *, 'Allocated Lattice Size (',m,') does not match request (',nx,') for lattice%im!'
          stop
          endif
      endif

c Lattice%jm
c ----------
      if(.not.associated(lattice%jm)) then
      allocate ( lattice%jm(0:ny-1) )
      do n=0,ny-1
      lattice%jm(n) = 0
      enddo
      else
             m=size(lattice%jm)
          if(m.ne.ny) then
          print *, 'Allocated Lattice Size (',m,') does not match request (',ny,') for lattice%jm!'
          stop
          endif
      endif

c Lattice%npeg
c ------------
      if(.not.associated(lattice%npeg)) then
      allocate ( lattice%npeg(0:nx-1) )
      do n=0,nx-1
      lattice%npeg(n) = 0
      enddo
      else
             m=size(lattice%npeg)
          if(m.ne.nx) then
          print *, 'Allocated Lattice Size (',m,') does not match request (',nx,') for lattice%npeg!'
          stop
          endif
      endif

      return
      end subroutine create_dynamics_lattice

C **********************************************************************
      subroutine destroy_dynamics_lattice (lattice)
C **********************************************************************
      implicit none
      type ( dynamics_lattice_type) lattice
      if(associated( lattice%im        )) deallocate ( lattice%im        )
      if(associated( lattice%jm        )) deallocate ( lattice%jm        )
      if(associated( lattice%npeg      )) deallocate ( lattice%npeg      )
      if(associated( lattice%ppeg      )) deallocate ( lattice%ppeg      )
      if(associated( lattice%img       )) deallocate ( lattice%img       )
      if(associated( lattice%im0       )) deallocate ( lattice%im0       )
      if(associated( lattice%ilocal    )) deallocate ( lattice%ilocal    )
      if(associated( lattice%jlocal    )) deallocate ( lattice%jlocal    )
      if(associated( lattice%iglobal   )) deallocate ( lattice%iglobal   )
      if(associated( lattice%jglobal   )) deallocate ( lattice%jglobal   )
      if(associated( lattice%peiglobal )) deallocate ( lattice%peiglobal )
      if(associated( lattice%pejglobal )) deallocate ( lattice%pejglobal )
      return
      end subroutine destroy_dynamics_lattice

      endmodule g3_dynamics_lattice_module

