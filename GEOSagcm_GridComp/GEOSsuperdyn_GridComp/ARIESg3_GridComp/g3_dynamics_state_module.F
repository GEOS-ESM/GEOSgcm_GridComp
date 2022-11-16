      module  g3_dynamics_state_module

c Declare Used Modules
c --------------------
      use g3_dynamics_lattice_module

      implicit none

c Dynamics State Variables
c ------------------------
      type    dynamics_vars_type
              real(kind=8) , pointer ::    p(:,:)
              real(kind=8) , pointer ::    u(:,:,:)
              real(kind=8) , pointer ::    v(:,:,:)
              real(kind=8) , pointer ::    t(:,:,:)
              real(kind=8) , pointer ::    q(:,:,:,:)
      endtype dynamics_vars_type


c Dynamics Grid
c -------------
      type    dynamics_grid_type
              type ( dynamics_lattice_type ) lattice
              real(kind=8)    ptop, pimean
              real(kind=8)    lam_np, phi_np, lam_0
              real(kind=8)    dl, dp
              integer im,jm,lm,n
              integer ntracer, ptracer
              real(kind=8)  , pointer :: dlam(:)
              real(kind=8)  , pointer :: dphi(:)
              real(kind=8)  , pointer :: sige(:)
              real(kind=8)  , pointer ::  sig(:)
              real(kind=8)  , pointer :: dsig(:)
              real(kind=8)  , pointer ::  alf(:)
              real(kind=8)  , pointer ::  bet(:)
      endtype dynamics_grid_type


c Define Dynamics State
c ---------------------
      type    dynamics_state_type
              type ( dynamics_grid_type ) grid
              type ( dynamics_vars_type ) vars(2)
      endtype dynamics_state_type

      contains

C **********************************************************************
      subroutine create_dynamics (state,im,jm,lm,nt)
C **********************************************************************
      implicit none
      type ( dynamics_state_type) state
      integer im,jm,lm,nt
      call create_dynamics_vars (state%vars(1),im,jm,lm,nt)
      call create_dynamics_vars (state%vars(2),im,jm,lm,nt)
      call create_dynamics_grid (state%grid,state%grid%lattice%imglobal,
     .                                      state%grid%lattice%jmglobal,lm)
      return
      end subroutine create_dynamics

C **********************************************************************
      subroutine create_dynamics_vars (vars,im,jm,lm,nt)
C **********************************************************************
      implicit none
      include 'g3_mymalloc_interface'
      type ( dynamics_vars_type) vars
      integer im,jm,lm,nt
      call mymalloc ( vars%p,im,jm )
      call mymalloc ( vars%u,im,jm,lm )
      call mymalloc ( vars%v,im,jm,lm )
      call mymalloc ( vars%t,im,jm,lm )
      call mymalloc ( vars%q,im,jm,lm,nt )
      return
      end subroutine create_dynamics_vars

C **********************************************************************
      subroutine create_dynamics_grid (grid,im,jm,lm)
C **********************************************************************
      implicit none
      include 'g3_mymalloc_interface'
      type ( dynamics_grid_type) grid
      integer im,jm,lm
      call mymalloc ( grid%dlam,im )
      call mymalloc ( grid%dphi,jm )
      call mymalloc ( grid%alf ,lm+1 )
      call mymalloc ( grid%bet ,lm+1 )
      call mymalloc ( grid%sige,lm+1 )
      call mymalloc ( grid% sig,lm )
      call mymalloc ( grid%dsig,lm )
      return
      end subroutine create_dynamics_grid

C **********************************************************************
      subroutine destroy_dynamics (state)
C **********************************************************************
      implicit none
      type ( dynamics_state_type) state
      call destroy_dynamics_vars (state%vars(1))
      call destroy_dynamics_vars (state%vars(2))
      call destroy_dynamics_grid (state%grid)
      return
      end subroutine destroy_dynamics

C **********************************************************************
      subroutine destroy_dynamics_vars (vars)
C **********************************************************************
      implicit none
      type ( dynamics_vars_type) vars
      if(associated( vars%p )) deallocate ( vars%p )
      if(associated( vars%u )) deallocate ( vars%u )
      if(associated( vars%v )) deallocate ( vars%v )
      if(associated( vars%t )) deallocate ( vars%t )
      if(associated( vars%q )) deallocate ( vars%q )
      return
      end subroutine destroy_dynamics_vars

C **********************************************************************
      subroutine destroy_dynamics_grid (grid)
C **********************************************************************
      implicit none
      type ( dynamics_grid_type) grid
      if(associated( grid%dlam )) deallocate ( grid%dlam )
      if(associated( grid%dphi )) deallocate ( grid%dphi )
      if(associated( grid%alf  )) deallocate ( grid%alf  )
      if(associated( grid%bet  )) deallocate ( grid%bet  )
      if(associated( grid%sige )) deallocate ( grid%sige )
      if(associated( grid%sig  )) deallocate ( grid%sig  )
      if(associated( grid%dsig )) deallocate ( grid%dsig )
      return
      end subroutine destroy_dynamics_grid


      endmodule g3_dynamics_state_module
